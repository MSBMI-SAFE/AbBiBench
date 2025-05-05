from collections.abc import Sequence
from typing import Optional, Tuple, Union
import torch
import torch.utils.checkpoint
from torch import nn
from torch.nn import BCEWithLogitsLoss, CrossEntropyLoss, MSELoss
from transformers.activations import ACT2FN
from transformers.modeling_outputs import (
    BaseModelOutput,
    MaskedLMOutput,
    SequenceClassifierOutput,
    TokenClassifierOutput,
)
from transformers.modeling_utils import PreTrainedModel
from .configuration_prosst import ProSSTConfig
import torch.nn.functional as F


def build_relative_position(query_size, key_size, device):
    """
    Build relative position according to the query and key

    We assume the absolute position of query \\(P_q\\) is range from (0, query_size) and the absolute position of key
    \\(P_k\\) is range from (0, key_size), The relative positions from query to key is \\(R_{q \\rightarrow k} = P_q -
    P_k\\)

    Args:
        query_size (int): the length of query
        key_size (int): the length of key

    Return:
        `torch.LongTensor`: A tensor with shape [1, query_size, key_size]

    """

    q_ids = torch.arange(query_size, dtype=torch.long, device=device)
    k_ids = torch.arange(key_size, dtype=torch.long, device=device)
    rel_pos_ids = q_ids[:, None] - k_ids.view(1, -1).repeat(query_size, 1)
    rel_pos_ids = rel_pos_ids[:query_size, :]
    rel_pos_ids = rel_pos_ids.unsqueeze(0)
    return rel_pos_ids


@torch.jit.script
def c2p_dynamic_expand(c2p_pos, query_layer, relative_pos):
    return c2p_pos.expand(
        [
            query_layer.size(0),
            query_layer.size(1),
            query_layer.size(2),
            relative_pos.size(-1),
        ]
    )


@torch.jit.script
def p2c_dynamic_expand(c2p_pos, query_layer, key_layer):
    return c2p_pos.expand(
        [
            query_layer.size(0),
            query_layer.size(1),
            key_layer.size(-2),
            key_layer.size(-2),
        ]
    )


@torch.jit.script
def pos_dynamic_expand(pos_index, p2c_att, key_layer):
    return pos_index.expand(
        p2c_att.size()[:2] + (pos_index.size(-2), key_layer.size(-2))
    )


def rotate_half(x):
    x1, x2 = x.chunk(2, dim=-1)
    return torch.cat((-x2, x1), dim=-1)


def apply_rotary_pos_emb(x, cos, sin):
    cos = cos[:, :, : x.shape[-2], :]
    sin = sin[:, :, : x.shape[-2], :]

    return (x * cos) + (rotate_half(x) * sin)


class RotaryEmbedding(torch.nn.Module):
    """
    Rotary position embeddings based on those in
    [RoFormer](https://huggingface.co/docs/transformers/model_doc/roformer). Query and keys are transformed by rotation
    matrices which depend on their relative positions.
    """

    def __init__(self, dim: int):
        super().__init__()
        # Generate and save the inverse frequency buffer (non trainable)
        inv_freq = 1.0 / (
            10000 ** (torch.arange(0, dim, 2, dtype=torch.int64).float() / dim)
        )
        inv_freq = inv_freq
        self.register_buffer("inv_freq", inv_freq)

        self._seq_len_cached = None
        self._cos_cached = None
        self._sin_cached = None

    def _update_cos_sin_tables(self, x, seq_dimension=2):
        seq_len = x.shape[seq_dimension]

        # Reset the tables if the sequence length has changed,
        # or if we're on a new device (possibly due to tracing for instance)
        if seq_len != self._seq_len_cached or self._cos_cached.device != x.device:
            self._seq_len_cached = seq_len
            t = torch.arange(x.shape[seq_dimension], device=x.device).type_as(
                self.inv_freq
            )
            freqs = torch.outer(t, self.inv_freq)
            emb = torch.cat((freqs, freqs), dim=-1).to(x.device)

            self._cos_cached = emb.cos()[None, None, :, :]
            self._sin_cached = emb.sin()[None, None, :, :]

        return self._cos_cached, self._sin_cached

    def forward(
        self, q: torch.Tensor, k: torch.Tensor
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        self._cos_cached, self._sin_cached = self._update_cos_sin_tables(
            k, seq_dimension=-2
        )

        return (
            apply_rotary_pos_emb(q, self._cos_cached, self._sin_cached),
            apply_rotary_pos_emb(k, self._cos_cached, self._sin_cached),
        )


class MaskedConv1d(nn.Conv1d):
    """A masked 1-dimensional convolution layer.

    Takes the same arguments as torch.nn.Conv1D, except that the padding is set automatically.

         Shape:
            Input: (N, L, in_channels)
            input_mask: (N, L, 1), optional
            Output: (N, L, out_channels)
    """

    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        kernel_size: int,
        stride: int = 1,
        dilation: int = 1,
        groups: int = 1,
        bias: bool = True,
    ):
        """
        :param in_channels: input channels
        :param out_channels: output channels
        :param kernel_size: the kernel width
        :param stride: filter shift
        :param dilation: dilation factor
        :param groups: perform depth-wise convolutions
        :param bias: adds learnable bias to output
        """
        padding = dilation * (kernel_size - 1) // 2
        super().__init__(
            in_channels,
            out_channels,
            kernel_size,
            stride=stride,
            dilation=dilation,
            groups=groups,
            bias=bias,
            padding=padding,
        )

    def forward(self, x, input_mask=None):
        if input_mask is not None:
            x = x * input_mask
        return super().forward(x.transpose(1, 2)).transpose(1, 2)


class Attention1dPooling(nn.Module):
    def __init__(self, config):
        super().__init__()
        self.layer = MaskedConv1d(config.hidden_size, 1, 1)

    def forward(self, x, input_mask=None):
        batch_szie = x.shape[0]
        attn = self.layer(x)
        attn = attn.view(batch_szie, -1)
        if input_mask is not None:
            attn = attn.masked_fill_(
                ~input_mask.view(batch_szie, -1).bool(), float("-inf")
            )
        attn = F.softmax(attn, dim=-1).view(batch_szie, -1, 1)
        out = (attn * x).sum(dim=1)
        return out


class MeanPooling(nn.Module):
    """Mean Pooling for sentence-level classification tasks."""

    def __init__(self):
        super().__init__()

    def forward(self, features, input_mask=None):
        if input_mask is not None:
            # Applying input_mask to zero out masked values
            masked_features = features * input_mask.unsqueeze(2)
            sum_features = torch.sum(masked_features, dim=1)
            mean_pooled_features = sum_features / input_mask.sum(dim=1, keepdim=True)
        else:
            mean_pooled_features = torch.mean(features, dim=1)
        return mean_pooled_features


class ContextPooler(nn.Module):
    def __init__(self, config):
        super().__init__()
        scale_hidden = getattr(config, "scale_hidden", 1)
        if config.pooling_head == "mean":
            self.mean_pooling = MeanPooling()
        elif config.pooling_head == "attention":
            self.mean_pooling = Attention1dPooling(config)
        self.dense = nn.Linear(
            config.pooler_hidden_size, scale_hidden * config.pooler_hidden_size
        )
        self.dropout = nn.Dropout(config.pooler_dropout)
        self.config = config

    def forward(self, hidden_states, input_mask=None):
        # We "pool" the model by simply taking the hidden state corresponding
        # to the first token.

        context_token = self.mean_pooling(hidden_states, input_mask)
        context_token = self.dropout(context_token)
        pooled_output = self.dense(context_token)
        pooled_output = torch.tanh(pooled_output)
        return pooled_output

    @property
    def output_dim(self):
        return self.config.hidden_size


class ProSSTLayerNorm(nn.Module):
    """LayerNorm module in the TF style (epsilon inside the square root)."""

    def __init__(self, size, eps=1e-12):
        super().__init__()
        self.weight = nn.Parameter(torch.ones(size))
        self.bias = nn.Parameter(torch.zeros(size))
        self.variance_epsilon = eps

    def forward(self, hidden_states):
        input_type = hidden_states.dtype
        hidden_states = hidden_states.float()
        mean = hidden_states.mean(-1, keepdim=True)
        variance = (hidden_states - mean).pow(2).mean(-1, keepdim=True)
        hidden_states = (hidden_states - mean) / torch.sqrt(
            variance + self.variance_epsilon
        )
        hidden_states = hidden_states.to(input_type)
        y = self.weight * hidden_states + self.bias
        return y


class DisentangledSelfAttention(nn.Module):

    def __init__(self, config: ProSSTConfig):
        super().__init__()
        if config.hidden_size % config.num_attention_heads != 0:
            raise ValueError(
                f"The hidden size ({config.hidden_size}) is not a multiple of the number of attention "
                f"heads ({config.num_attention_heads})"
            )
        self.num_attention_heads = config.num_attention_heads
        self.attention_head_size = int(config.hidden_size / config.num_attention_heads)
        self.all_head_size = self.num_attention_heads * self.attention_head_size

        # Q, K, V projection layers
        self.query = nn.Linear(config.hidden_size, self.all_head_size)
        self.key = nn.Linear(config.hidden_size, self.all_head_size)
        self.value = nn.Linear(config.hidden_size, self.all_head_size)

        # AA->SS, AA->POS, SS->AA, POS->AA and AA->AA attention layers
        self.pos_att_type = (
            config.pos_att_type if config.pos_att_type is not None else []
        )

        self.relative_attention = getattr(config, "relative_attention", False)
        self.position_embedding_type = getattr(
            config, "position_embedding_type", "relative"
        )
        if self.position_embedding_type == "rotary":
            self.rotary_embeddings = RotaryEmbedding(dim=self.attention_head_size)
            if self.relative_attention:

                if "aa2ss" in self.pos_att_type:
                    self.ss_proj = nn.Linear(
                        config.hidden_size, self.all_head_size, bias=False
                    )

                if "ss2aa" in self.pos_att_type:
                    self.ss_q_proj = nn.Linear(config.hidden_size, self.all_head_size)

        elif self.position_embedding_type == "relative":
            if self.relative_attention:
                self.max_relative_positions = getattr(
                    config, "max_relative_positions", -1
                )
                if self.max_relative_positions < 1:
                    self.max_relative_positions = config.max_position_embeddings
                self.pos_dropout = nn.Dropout(config.hidden_dropout_prob)

                # amino acid to position
                if "aa2pos" in self.pos_att_type:
                    self.pos_proj = nn.Linear(
                        config.hidden_size, self.all_head_size, bias=False
                    )  # Key

                if "pos2aa" in self.pos_att_type:
                    self.pos_q_proj = nn.Linear(
                        config.hidden_size, self.all_head_size
                    )  # Query

                if "aa2ss" in self.pos_att_type:
                    self.ss_proj = nn.Linear(
                        config.hidden_size, self.all_head_size, bias=False
                    )

                if "ss2aa" in self.pos_att_type:
                    self.ss_q_proj = nn.Linear(config.hidden_size, self.all_head_size)

        self.dropout = nn.Dropout(config.attention_probs_dropout_prob)

    def transpose_for_scores(self, x):
        # x [batch_size, seq_len, all_head_size]
        new_x_shape = x.size()[:-1] + (self.num_attention_heads, -1)
        # x [batch_size, seq_len, num_attention_heads, attention_head_size]
        x = x.view(new_x_shape)
        # x [batch_size, num_attention_heads, seq_len, attention_head_size]
        return x.permute(0, 2, 1, 3)

    def forward(
        self,
        hidden_states,
        attention_mask,
        output_attentions=False,
        query_states=None,
        relative_pos=None,
        rel_embeddings=None,
        ss_hidden_states=None,
    ):
        query_layer = self.transpose_for_scores(self.query(hidden_states))
        key_layer = self.transpose_for_scores(self.key(hidden_states))
        value_layer = self.transpose_for_scores(self.value(hidden_states))

        if self.position_embedding_type == "rotary":
            query_layer, key_layer = self.rotary_embeddings(query_layer, key_layer)

        rel_att = None
        scale_factor = 1 + len(self.pos_att_type)
        scale = torch.sqrt(
            torch.tensor(query_layer.size(-1), dtype=torch.float) * scale_factor
        )
        query_layer = query_layer / scale.to(dtype=query_layer.dtype)

        # [batch_size, num_attention_heads, seq_len, seq_len]
        attention_scores = torch.matmul(query_layer, key_layer.transpose(-1, -2))

        if self.relative_attention:
            if self.position_embedding_type == "relative":
                rel_embeddings = self.pos_dropout(rel_embeddings)
            rel_att = self.disentangled_att_bias(
                query_layer,
                key_layer,
                relative_pos,
                rel_embeddings,
                scale_factor,
                ss_hidden_states,
            )

        if rel_att is not None:
            attention_scores = attention_scores + rel_att

        rmask = ~(attention_mask.to(torch.bool))
        attention_probs = attention_scores.masked_fill(rmask, float("-inf"))
        attention_probs = torch.softmax(attention_probs, -1)
        attention_probs = attention_probs.masked_fill(rmask, 0.0)
        # attention_probs = XSoftmax.apply(attention_scores, attention_mask, -1)
        attention_probs = self.dropout(attention_probs)

        context_layer = torch.matmul(attention_probs, value_layer)
        context_layer = context_layer.permute(0, 2, 1, 3).contiguous()
        new_context_layer_shape = context_layer.size()[:-2] + (-1,)
        context_layer = context_layer.view(new_context_layer_shape)
        if output_attentions:
            return (context_layer, attention_probs)
        else:
            return context_layer

    def disentangled_att_bias(
        self,
        query_layer,
        key_layer,
        relative_pos,
        rel_embeddings,
        scale_factor,
        ss_hidden_states,
    ):
        if self.position_embedding_type == "relative":
            if relative_pos is None:
                q = query_layer.size(-2)
                relative_pos = build_relative_position(
                    q, key_layer.size(-2), query_layer.device
                )
            if relative_pos.dim() == 2:
                relative_pos = relative_pos.unsqueeze(0).unsqueeze(0)
            elif relative_pos.dim() == 3:
                relative_pos = relative_pos.unsqueeze(1)
            # bxhxqxk
            elif relative_pos.dim() != 4:
                raise ValueError(
                    f"Relative position ids must be of dim 2 or 3 or 4. {relative_pos.dim()}"
                )

            att_span = min(
                max(query_layer.size(-2), key_layer.size(-2)),
                self.max_relative_positions,
            )
            relative_pos = relative_pos.long().to(query_layer.device)
            rel_embeddings = rel_embeddings[
                self.max_relative_positions
                - att_span : self.max_relative_positions
                + att_span,
                :,
            ].unsqueeze(0)

            score = 0

            if "aa2pos" in self.pos_att_type:
                pos_key_layer = self.pos_proj(rel_embeddings)
                pos_key_layer = self.transpose_for_scores(pos_key_layer)
                aa2p_att = torch.matmul(query_layer, pos_key_layer.transpose(-1, -2))
                aa2p_pos = torch.clamp(relative_pos + att_span, 0, att_span * 2 - 1)
                aa2p_att = torch.gather(
                    aa2p_att,
                    dim=-1,
                    index=c2p_dynamic_expand(aa2p_pos, query_layer, relative_pos),
                )
                score += aa2p_att

            if "pos2aa" in self.pos_att_type:
                pos_query_layer = self.pos_q_proj(rel_embeddings)
                pos_query_layer = self.transpose_for_scores(pos_query_layer)
                pos_query_layer /= torch.sqrt(
                    torch.tensor(pos_query_layer.size(-1), dtype=torch.float)
                    * scale_factor
                )
                if query_layer.size(-2) != key_layer.size(-2):
                    r_pos = build_relative_position(
                        key_layer.size(-2), key_layer.size(-2), query_layer.device
                    )
                else:
                    r_pos = relative_pos
                p2aa_pos = torch.clamp(-r_pos + att_span, 0, att_span * 2 - 1)
                p2aa_att = torch.matmul(
                    key_layer,
                    pos_query_layer.transpose(-1, -2).to(dtype=key_layer.dtype),
                )
                p2aa_att = torch.gather(
                    p2aa_att,
                    dim=-1,
                    index=p2c_dynamic_expand(p2aa_pos, query_layer, key_layer),
                ).transpose(-1, -2)

                if query_layer.size(-2) != key_layer.size(-2):
                    pos_index = relative_pos[:, :, :, 0].unsqueeze(-1)
                    p2aa_att = torch.gather(
                        p2aa_att,
                        dim=-2,
                        index=pos_dynamic_expand(pos_index, p2aa_att, key_layer),
                    )
                score += p2aa_att

            # content -> structure
            if "aa2ss" in self.pos_att_type:
                assert ss_hidden_states is not None
                ss_key_layer = self.ss_proj(ss_hidden_states)
                ss_key_layer = self.transpose_for_scores(ss_key_layer)
                # [batch_size, num_attention_heads, seq_len, seq_len]
                aa2ss_att = torch.matmul(query_layer, ss_key_layer.transpose(-1, -2))
                score += aa2ss_att

            if "ss2aa" in self.pos_att_type:
                assert ss_hidden_states is not None
                ss_query_layer = self.ss_q_proj(ss_hidden_states)
                ss_query_layer = self.transpose_for_scores(ss_query_layer)
                ss_query_layer /= torch.sqrt(
                    torch.tensor(ss_query_layer.size(-1), dtype=torch.float)
                    * scale_factor
                )
                ss2aa_att = torch.matmul(
                    key_layer, query_layer.transpose(-1, -2).to(dtype=key_layer.dtype)
                )
                score += ss2aa_att
            return score
        elif self.position_embedding_type == "rotary":
            score = 0
            if "aa2ss" in self.pos_att_type:
                assert ss_hidden_states is not None
                ss_key_layer = self.ss_proj(ss_hidden_states)
                ss_key_layer = self.transpose_for_scores(ss_key_layer)
                aa2ss_att = torch.matmul(query_layer, ss_key_layer.transpose(-1, -2))
                score += aa2ss_att

            if "ss2aa" in self.pos_att_type:
                assert ss_hidden_states is not None
                ss_query_layer = self.ss_q_proj(ss_hidden_states)
                ss_query_layer = self.transpose_for_scores(ss_query_layer)
                ss_query_layer /= torch.sqrt(
                    torch.tensor(ss_query_layer.size(-1), dtype=torch.float)
                    * scale_factor
                )
                ss2aa_att = torch.matmul(
                    key_layer, query_layer.transpose(-1, -2).to(dtype=key_layer.dtype)
                )
                score += ss2aa_att
            return score


class ProSSTSelfOutput(nn.Module):
    def __init__(self, config):
        super().__init__()
        self.dense = nn.Linear(config.hidden_size, config.hidden_size)
        self.LayerNorm = ProSSTLayerNorm(config.hidden_size, config.layer_norm_eps)
        self.dropout = nn.Dropout(config.hidden_dropout_prob)

    def forward(self, hidden_states, input_tensor):
        hidden_states = self.dense(hidden_states)
        hidden_states = self.dropout(hidden_states)
        hidden_states = self.LayerNorm(hidden_states + input_tensor)
        return hidden_states


class ProSSTAttention(nn.Module):
    def __init__(self, config):
        super().__init__()
        self.self = DisentangledSelfAttention(config)
        self.output = ProSSTSelfOutput(config)
        self.config = config

    def forward(
        self,
        hidden_states,
        attention_mask,
        output_attentions=False,
        query_states=None,
        relative_pos=None,
        rel_embeddings=None,
        ss_hidden_states=None,
    ):
        self_output = self.self(
            hidden_states,
            attention_mask,
            output_attentions,
            query_states=query_states,
            relative_pos=relative_pos,
            rel_embeddings=rel_embeddings,
            ss_hidden_states=ss_hidden_states,
        )
        if output_attentions:
            self_output, att_matrix = self_output
        if query_states is None:
            query_states = hidden_states
        attention_output = self.output(self_output, query_states)

        if output_attentions:
            return (attention_output, att_matrix)
        else:
            return attention_output


class ProSSTIntermediate(nn.Module):
    def __init__(self, config):
        super().__init__()
        self.dense = nn.Linear(config.hidden_size, config.intermediate_size)
        if isinstance(config.hidden_act, str):
            self.intermediate_act_fn = ACT2FN[config.hidden_act]
        else:
            self.intermediate_act_fn = config.hidden_act

    def forward(self, hidden_states: torch.Tensor) -> torch.Tensor:
        hidden_states = self.dense(hidden_states)
        hidden_states = self.intermediate_act_fn(hidden_states)
        return hidden_states


class ProSSTOutput(nn.Module):
    def __init__(self, config):
        super().__init__()
        self.dense = nn.Linear(config.intermediate_size, config.hidden_size)
        self.LayerNorm = ProSSTLayerNorm(config.hidden_size, config.layer_norm_eps)
        self.dropout = nn.Dropout(config.hidden_dropout_prob)
        self.config = config

    def forward(self, hidden_states, input_tensor):
        hidden_states = self.dense(hidden_states)
        hidden_states = self.dropout(hidden_states)
        hidden_states = self.LayerNorm(hidden_states + input_tensor)
        return hidden_states


class ProSSTLayer(nn.Module):
    def __init__(self, config):
        super().__init__()
        self.attention = ProSSTAttention(config)
        self.intermediate = ProSSTIntermediate(config)
        self.output = ProSSTOutput(config)

    def forward(
        self,
        hidden_states,
        attention_mask,
        query_states=None,
        relative_pos=None,
        rel_embeddings=None,
        output_attentions=False,
        ss_hidden_states=None,
    ):
        attention_output = self.attention(
            hidden_states,
            attention_mask,
            output_attentions=output_attentions,
            query_states=query_states,
            relative_pos=relative_pos,
            rel_embeddings=rel_embeddings,
            ss_hidden_states=ss_hidden_states,
        )
        if output_attentions:
            attention_output, att_matrix = attention_output
        intermediate_output = self.intermediate(attention_output)
        layer_output = self.output(intermediate_output, attention_output)
        if output_attentions:
            return (layer_output, att_matrix)
        else:
            return layer_output


class ProSSTEncoder(nn.Module):
    """Modified BertEncoder with relative position bias support"""

    def __init__(self, config):
        super().__init__()
        self.layer = nn.ModuleList(
            [ProSSTLayer(config) for _ in range(config.num_hidden_layers)]
        )
        self.relative_attention = getattr(config, "relative_attention", False)
        if self.relative_attention:
            self.max_relative_positions = getattr(config, "max_relative_positions", -1)
            if self.max_relative_positions < 1:
                self.max_relative_positions = config.max_position_embeddings
            self.rel_embeddings = nn.Embedding(
                self.max_relative_positions * 2, config.hidden_size
            )
        self.gradient_checkpointing = False

    def get_rel_embedding(self):
        rel_embeddings = self.rel_embeddings.weight if self.relative_attention else None
        return rel_embeddings

    def get_attention_mask(self, attention_mask):
        if attention_mask.dim() <= 2:
            extended_attention_mask = attention_mask.unsqueeze(1).unsqueeze(2)
            attention_mask = extended_attention_mask * extended_attention_mask.squeeze(
                -2
            ).unsqueeze(-1)
        elif attention_mask.dim() == 3:
            attention_mask = attention_mask.unsqueeze(1)

        return attention_mask

    def get_rel_pos(self, hidden_states, query_states=None, relative_pos=None):
        if self.relative_attention and relative_pos is None:
            q = (
                query_states.size(-2)
                if query_states is not None
                else hidden_states.size(-2)
            )
            relative_pos = build_relative_position(
                q, hidden_states.size(-2), hidden_states.device
            )
        return relative_pos

    def forward(
        self,
        hidden_states,
        attention_mask,
        output_hidden_states=True,
        output_attentions=False,
        query_states=None,
        relative_pos=None,
        ss_hidden_states=None,
        return_dict=True,
    ):
        attention_mask = self.get_attention_mask(attention_mask)
        relative_pos = self.get_rel_pos(hidden_states, query_states, relative_pos)

        all_hidden_states = () if output_hidden_states else None
        all_attentions = () if output_attentions else None

        if isinstance(hidden_states, Sequence):
            next_kv = hidden_states[0]
        else:
            next_kv = hidden_states
        rel_embeddings = self.get_rel_embedding()
        for i, layer_module in enumerate(self.layer):
            if output_hidden_states:
                all_hidden_states = all_hidden_states + (hidden_states,)

            if self.gradient_checkpointing and self.training:

                def create_custom_forward(module):
                    def custom_forward(*inputs):
                        return module(*inputs, output_attentions)

                    return custom_forward

                hidden_states = torch.utils.checkpoint.checkpoint(
                    create_custom_forward(layer_module),
                    next_kv,
                    attention_mask,
                    query_states,
                    relative_pos,
                    rel_embeddings,
                    ss_hidden_states,
                )
            else:
                hidden_states = layer_module(
                    next_kv,
                    attention_mask,
                    query_states=query_states,
                    relative_pos=relative_pos,
                    rel_embeddings=rel_embeddings,
                    output_attentions=output_attentions,
                    ss_hidden_states=ss_hidden_states,
                )

            if output_attentions:
                hidden_states, att_m = hidden_states

            if query_states is not None:
                query_states = hidden_states
                if isinstance(hidden_states, Sequence):
                    next_kv = hidden_states[i + 1] if i + 1 < len(self.layer) else None
            else:
                next_kv = hidden_states

            if output_attentions:
                all_attentions = all_attentions + (att_m,)

        if output_hidden_states:
            all_hidden_states = all_hidden_states + (hidden_states,)

        if not return_dict:
            return tuple(
                v
                for v in [hidden_states, all_hidden_states, all_attentions]
                if v is not None
            )
        return BaseModelOutput(
            last_hidden_state=hidden_states,
            hidden_states=all_hidden_states,
            attentions=all_attentions,
        )


class ProSSTEmbeddings(nn.Module):
    """Construct the embeddings from word, position and token_type embeddings."""

    def __init__(self, config):
        super().__init__()
        pad_token_id = getattr(config, "pad_token_id", 0)
        self.embedding_size = getattr(config, "embedding_size", config.hidden_size)
        self.word_embeddings = nn.Embedding(
            config.vocab_size, self.embedding_size, padding_idx=pad_token_id
        )

        self.position_biased_input = getattr(config, "position_biased_input", False)
        if not self.position_biased_input:
            self.position_embeddings = None
        else:
            # assert getattr(config, "position_embedding_type", "relative") == "absolute"
            self.position_embeddings = nn.Embedding(
                config.max_position_embeddings, self.embedding_size
            )

        if config.type_vocab_size > 0:
            self.token_type_embeddings = nn.Embedding(
                config.type_vocab_size, self.embedding_size
            )

        if config.ss_vocab_size > 0:
            self.ss_embeddings = nn.Embedding(config.ss_vocab_size, self.embedding_size)
            self.ss_layer_norm = ProSSTLayerNorm(
                config.hidden_size, config.layer_norm_eps
            )

        if self.embedding_size != config.hidden_size:
            self.embed_proj = nn.Linear(
                self.embedding_size, config.hidden_size, bias=False
            )
        self.LayerNorm = ProSSTLayerNorm(config.hidden_size, config.layer_norm_eps)

        self.dropout = nn.Dropout(config.hidden_dropout_prob)
        self.config = config

        # position_ids (1, len position emb) is contiguous in memory and exported when serialized
        if self.position_biased_input:
            self.register_buffer(
                "position_ids",
                torch.arange(config.max_position_embeddings).expand((1, -1)),
                persistent=False,
            )

    def forward(
        self,
        input_ids=None,
        ss_input_ids=None,
        token_type_ids=None,
        position_ids=None,
        mask=None,
        inputs_embeds=None,
    ):
        if input_ids is not None:
            input_shape = input_ids.size()
        else:
            input_shape = inputs_embeds.size()[:-1]

        seq_length = input_shape[1]

        if position_ids is None and self.position_biased_input:
            position_ids = self.position_ids[:, :seq_length]
            if seq_length > position_ids.size(1):
                zero_padding = (
                    torch.zeros(
                        (input_shape[0], seq_length - position_ids.size(1)),
                        dtype=torch.long,
                        device=position_ids.device,
                    )
                    + 2047
                )
                position_ids = torch.cat([position_ids, zero_padding], dim=1)

        if token_type_ids is None:
            token_type_ids = torch.zeros(
                input_shape, dtype=torch.long, device=self.position_ids.device
            )

        if inputs_embeds is None:
            if self.config.token_dropout:
                inputs_embeds = self.word_embeddings(input_ids)
                inputs_embeds.masked_fill_(
                    (input_ids == self.config.mask_token_id).unsqueeze(-1), 0.0
                )
                mask_ratio_train = self.config.mlm_probability * 0.8
                src_lengths = mask.sum(dim=-1)
                mask_ratio_observed = (input_ids == self.config.mask_token_id).sum(
                    -1
                ).to(inputs_embeds.dtype) / src_lengths
                inputs_embeds = (
                    inputs_embeds
                    * (1 - mask_ratio_train)
                    / (1 - mask_ratio_observed)[:, None, None]
                )
            else:
                inputs_embeds = self.word_embeddings(input_ids)

        if self.position_embeddings is not None and self.position_biased_input:
            position_embeddings = self.position_embeddings(position_ids.long())
        else:
            position_embeddings = torch.zeros_like(inputs_embeds)

        embeddings = inputs_embeds
        if self.position_biased_input:
            embeddings += position_embeddings
        if self.config.type_vocab_size > 0:
            token_type_embeddings = self.token_type_embeddings(token_type_ids)
            embeddings += token_type_embeddings

        if self.embedding_size != self.config.hidden_size:
            embeddings = self.embed_proj(embeddings)

        embeddings = self.LayerNorm(embeddings)

        if mask is not None:
            if mask.dim() != embeddings.dim():
                if mask.dim() == 4:
                    mask = mask.squeeze(1).squeeze(1)
                mask = mask.unsqueeze(2)
            mask = mask.to(embeddings.dtype)
            embeddings = embeddings * mask

        embeddings = self.dropout(embeddings)

        if self.config.ss_vocab_size > 0:
            ss_embeddings = self.ss_embeddings(ss_input_ids)
            ss_embeddings = self.ss_layer_norm(ss_embeddings)
            if mask is not None:
                if mask.dim() != ss_embeddings.dim():
                    if mask.dim() == 4:
                        mask = mask.squeeze(1).squeeze(1)
                    mask = mask.unsqueeze(2)
                mask = mask.to(ss_embeddings.dtype)
                ss_embeddings = ss_embeddings * mask
                ss_embeddings = self.dropout(ss_embeddings)
            return embeddings, ss_embeddings

        return embeddings, None


class ProSSTPreTrainedModel(PreTrainedModel):
    """
    An abstract class to handle weights initialization and a simple interface for downloading and loading pretrained
    models.
    """

    config_class = ProSSTConfig
    base_model_prefix = "ProSST"
    _keys_to_ignore_on_load_unexpected = ["position_embeddings"]
    supports_gradient_checkpointing = True

    def _init_weights(self, module):
        """Initialize the weights."""
        if isinstance(module, nn.Linear):
            # Slightly different from the TF version which uses truncated_normal for initialization
            # cf https://github.com/pytorch/pytorch/pull/5617
            module.weight.data.normal_(mean=0.0, std=self.config.initializer_range)
            if module.bias is not None:
                module.bias.data.zero_()
        elif isinstance(module, nn.Embedding):
            module.weight.data.normal_(mean=0.0, std=self.config.initializer_range)
            if module.padding_idx is not None:
                module.weight.data[module.padding_idx].zero_()

    def _set_gradient_checkpointing(self, module, value=False):
        if isinstance(module, ProSSTEncoder):
            module.gradient_checkpointing = value


class ProSSTModel(ProSSTPreTrainedModel):
    def __init__(self, config):
        super().__init__(config)

        self.embeddings = ProSSTEmbeddings(config)
        self.encoder = ProSSTEncoder(config)
        self.config = config
        # Initialize weights and apply final processing
        self.post_init()

    def get_input_embeddings(self):
        return self.embeddings.word_embeddings

    def set_input_embeddings(self, new_embeddings):
        self.embeddings.word_embeddings = new_embeddings

    def _prune_heads(self, heads_to_prune):
        """
        Prunes heads of the model. heads_to_prune: dict of {layer_num: list of heads to prune in this layer} See base
        class PreTrainedModel
        """
        raise NotImplementedError(
            "The prune function is not implemented in DeBERTa model."
        )

    def forward(
        self,
        input_ids: Optional[torch.Tensor] = None,
        ss_input_ids: Optional[torch.Tensor] = None,
        attention_mask: Optional[torch.Tensor] = None,
        token_type_ids: Optional[torch.Tensor] = None,
        position_ids: Optional[torch.Tensor] = None,
        inputs_embeds: Optional[torch.Tensor] = None,
        output_attentions: Optional[bool] = None,
        output_hidden_states: Optional[bool] = None,
        return_dict: Optional[bool] = None,
    ) -> Union[Tuple, BaseModelOutput]:
        output_attentions = (
            output_attentions
            if output_attentions is not None
            else self.config.output_attentions
        )
        output_hidden_states = (
            output_hidden_states
            if output_hidden_states is not None
            else self.config.output_hidden_states
        )
        return_dict = (
            return_dict if return_dict is not None else self.config.use_return_dict
        )

        if input_ids is not None and inputs_embeds is not None:
            raise ValueError(
                "You cannot specify both input_ids and inputs_embeds at the same time"
            )
        elif input_ids is not None:
            self.warn_if_padding_and_no_attention_mask(input_ids, attention_mask)
            input_shape = input_ids.size()
        elif inputs_embeds is not None:
            input_shape = inputs_embeds.size()[:-1]
        else:
            raise ValueError("You have to specify either input_ids or inputs_embeds")

        device = input_ids.device if input_ids is not None else inputs_embeds.device

        if attention_mask is None:
            attention_mask = torch.ones(input_shape, device=device)
        if token_type_ids is None:
            token_type_ids = torch.zeros(input_shape, dtype=torch.long, device=device)

        embedding_output, ss_embeddings = self.embeddings(
            input_ids=input_ids,
            ss_input_ids=ss_input_ids,
            token_type_ids=token_type_ids,
            position_ids=position_ids,
            mask=attention_mask,
            inputs_embeds=inputs_embeds,
        )

        encoder_outputs = self.encoder(
            embedding_output,
            attention_mask,
            output_hidden_states=True,
            output_attentions=output_attentions,
            return_dict=return_dict,
            ss_hidden_states=ss_embeddings,
        )
        encoded_layers = encoder_outputs[1]

        sequence_output = encoded_layers[-1]

        if not return_dict:
            return (sequence_output,) + encoder_outputs[
                (1 if output_hidden_states else 2) :
            ]

        return BaseModelOutput(
            last_hidden_state=sequence_output,
            hidden_states=(
                encoder_outputs.hidden_states if output_hidden_states else None
            ),
            attentions=encoder_outputs.attentions,
        )


class ProSSTPredictionHeadTransform(nn.Module):
    def __init__(self, config):
        super().__init__()
        self.embedding_size = getattr(config, "embedding_size", config.hidden_size)

        self.dense = nn.Linear(config.hidden_size, self.embedding_size)
        if isinstance(config.hidden_act, str):
            self.transform_act_fn = ACT2FN[config.hidden_act]
        else:
            self.transform_act_fn = config.hidden_act
        self.LayerNorm = nn.LayerNorm(self.embedding_size, eps=config.layer_norm_eps)

    def forward(self, hidden_states):
        hidden_states = self.dense(hidden_states)
        hidden_states = self.transform_act_fn(hidden_states)
        hidden_states = self.LayerNorm(hidden_states)
        return hidden_states


class ProSSTLMPredictionHead(nn.Module):
    def __init__(self, config):
        super().__init__()
        self.transform = ProSSTPredictionHeadTransform(config)

        self.embedding_size = getattr(config, "embedding_size", config.hidden_size)
        # The output weights are the same as the input embeddings, but there is
        # an output-only bias for each token.
        self.decoder = nn.Linear(self.embedding_size, config.vocab_size, bias=False)

        # self.bias = nn.Parameter(torch.zeros(config.vocab_size))
        # Need a link between the two variables so that the bias is correctly resized with `resize_token_embeddings`
        # self.decoder.bias = self.bias

    def forward(self, hidden_states):
        hidden_states = self.transform(hidden_states)
        hidden_states = self.decoder(hidden_states)
        return hidden_states


class ProSSTOnlyMLMHead(nn.Module):
    def __init__(self, config):
        super().__init__()
        self.predictions = ProSSTLMPredictionHead(config)

    def forward(self, sequence_output):
        prediction_scores = self.predictions(sequence_output)
        return prediction_scores


class ProSSTPreTrainedModel(PreTrainedModel):
    """
    An abstract class to handle weights initialization and a simple interface for downloading and loading pretrained
    models.
    """

    config_class = ProSSTConfig
    base_model_prefix = "ProSST"
    _keys_to_ignore_on_load_unexpected = ["position_embeddings"]
    supports_gradient_checkpointing = True

    def _init_weights(self, module):
        """Initialize the weights."""
        if isinstance(module, nn.Linear):
            # Slightly different from the TF version which uses truncated_normal for initialization
            # cf https://github.com/pytorch/pytorch/pull/5617
            module.weight.data.normal_(mean=0.0, std=self.config.initializer_range)
            if module.bias is not None:
                module.bias.data.zero_()
        elif isinstance(module, nn.Embedding):
            module.weight.data.normal_(mean=0.0, std=self.config.initializer_range)
            if module.padding_idx is not None:
                module.weight.data[module.padding_idx].zero_()

    def _set_gradient_checkpointing(self, module, value=False):
        if isinstance(module, ProSSTEncoder):
            module.gradient_checkpointing = value


class ProSSTForMaskedLM(ProSSTPreTrainedModel):
    _tied_weights_keys = [
        "cls.predictions.decoder.weight",
        "cls.predictions.decoder.bias",
    ]

    def __init__(self, config):
        super().__init__(config)

        self.prosst = ProSSTModel(config)
        self.cls = ProSSTOnlyMLMHead(config)

        # Initialize weights and apply final processing
        self.post_init()
        
    def get_input_embeddings(self):
        return self.prosst.embeddings.word_embeddings

    def get_output_embeddings(self):
        return self.cls.predictions.decoder

    def set_output_embeddings(self, new_embeddings):
        self.cls.predictions.decoder = new_embeddings

    def forward(
        self,
        input_ids: Optional[torch.Tensor] = None,
        ss_input_ids: Optional[torch.Tensor] = None,
        attention_mask: Optional[torch.Tensor] = None,
        token_type_ids: Optional[torch.Tensor] = None,
        position_ids: Optional[torch.Tensor] = None,
        inputs_embeds: Optional[torch.Tensor] = None,
        labels: Optional[torch.Tensor] = None,
        output_attentions: Optional[bool] = None,
        output_hidden_states: Optional[bool] = None,
        return_dict: Optional[bool] = None,
    ) -> Union[Tuple, MaskedLMOutput]:
        r"""
        labels (`torch.LongTensor` of shape `(batch_size, sequence_length)`, *optional*):
            Labels for computing the masked language modeling loss. Indices should be in `[-100, 0, ...,
            config.vocab_size]` (see `input_ids` docstring) Tokens with indices set to `-100` are ignored (masked), the
            loss is only computed for the tokens with labels in `[0, ..., config.vocab_size]`
        """

        return_dict = (
            return_dict if return_dict is not None else self.config.use_return_dict
        )

        outputs = self.prosst(
            input_ids,
            ss_input_ids=ss_input_ids,
            attention_mask=attention_mask,
            token_type_ids=token_type_ids,
            position_ids=position_ids,
            inputs_embeds=inputs_embeds,
            output_attentions=output_attentions,
            output_hidden_states=output_hidden_states,
            return_dict=return_dict,
        )

        sequence_output = outputs[0]
        prediction_scores = self.cls(sequence_output)

        masked_lm_loss = None
        if labels is not None:
            loss_fct = CrossEntropyLoss()  # -100 index = padding token
            masked_lm_loss = loss_fct(
                prediction_scores.view(-1, self.config.vocab_size), labels.view(-1)
            )

        if not return_dict:
            output = (prediction_scores,) + outputs[1:]
            return (
                ((masked_lm_loss,) + output) if masked_lm_loss is not None else output
            )

        return MaskedLMOutput(
            loss=masked_lm_loss,
            logits=prediction_scores,
            hidden_states=outputs.hidden_states,
            attentions=outputs.attentions,
        )


class ProSSTForSequenceClassification(ProSSTPreTrainedModel):
    def __init__(self, config):
        super().__init__(config)

        num_labels = getattr(config, "num_labels", 2)
        self.num_labels = num_labels
        self.scale_hidden = getattr(config, "scale_hidden", 1)
        self.prosst = ProSSTModel(config)
        self.pooler = ContextPooler(config)
        output_dim = self.pooler.output_dim * self.scale_hidden

        self.classifier = nn.Linear(output_dim, num_labels)
        drop_out = getattr(config, "cls_dropout", None)
        drop_out = self.config.hidden_dropout_prob if drop_out is None else drop_out
        self.dropout = nn.Dropout(drop_out)

        # Initialize weights and apply final processing
        self.post_init()

    def get_input_embeddings(self):
        return self.prosst.get_input_embeddings()

    def set_input_embeddings(self, new_embeddings):
        self.prosst.set_input_embeddings(new_embeddings)

    def forward(
        self,
        input_ids: Optional[torch.Tensor] = None,
        ss_input_ids: Optional[torch.Tensor] = None,
        attention_mask: Optional[torch.Tensor] = None,
        token_type_ids: Optional[torch.Tensor] = None,
        position_ids: Optional[torch.Tensor] = None,
        inputs_embeds: Optional[torch.Tensor] = None,
        labels: Optional[torch.Tensor] = None,
        output_attentions: Optional[bool] = None,
        output_hidden_states: Optional[bool] = None,
        return_dict: Optional[bool] = None,
    ) -> Union[Tuple, SequenceClassifierOutput]:
        r"""
        labels (`torch.LongTensor` of shape `(batch_size,)`, *optional*):
            Labels for computing the sequence classification/regression loss. Indices should be in `[0, ...,
            config.num_labels - 1]`. If `config.num_labels == 1` a regression loss is computed (Mean-Square loss), If
            `config.num_labels > 1` a classification loss is computed (Cross-Entropy).
        """
        return_dict = (
            return_dict if return_dict is not None else self.config.use_return_dict
        )

        outputs = self.prosst(
            input_ids,
            ss_input_ids=ss_input_ids,
            token_type_ids=token_type_ids,
            attention_mask=attention_mask,
            position_ids=position_ids,
            inputs_embeds=inputs_embeds,
            output_attentions=output_attentions,
            output_hidden_states=output_hidden_states,
            return_dict=return_dict,
        )

        encoder_layer = outputs[0]
        pooled_output = self.pooler(encoder_layer, attention_mask)
        pooled_output = self.dropout(pooled_output)
        logits = self.classifier(pooled_output)

        loss = None
        if labels is not None:
            if self.config.problem_type is None:
                if self.num_labels == 1:
                    # regression task
                    loss_fn = nn.MSELoss()
                    logits = logits.view(-1).to(labels.dtype)
                    loss = loss_fn(logits, labels.view(-1))
                elif labels.dim() == 1 or labels.size(-1) == 1:
                    label_index = (labels >= 0).nonzero()
                    labels = labels.long()
                    if label_index.size(0) > 0:
                        labeled_logits = torch.gather(
                            logits,
                            0,
                            label_index.expand(label_index.size(0), logits.size(1)),
                        )
                        labels = torch.gather(labels, 0, label_index.view(-1))
                        loss_fct = CrossEntropyLoss()
                        loss = loss_fct(
                            labeled_logits.view(-1, self.num_labels).float(),
                            labels.view(-1),
                        )
                    else:
                        loss = torch.tensor(0).to(logits)
                else:
                    log_softmax = nn.LogSoftmax(-1)
                    loss = -((log_softmax(logits) * labels).sum(-1)).mean()
            elif self.config.problem_type == "regression":
                loss_fct = MSELoss()
                if self.num_labels == 1:
                    loss = loss_fct(logits.squeeze(), labels.squeeze())
                else:
                    loss = loss_fct(logits, labels)
            elif self.config.problem_type == "binary_classification":
                loss_fct = BCEWithLogitsLoss()
                loss = loss_fct(logits.squeeze(), labels.squeeze().to(logits.dtype))
            elif self.config.problem_type == "single_label_classification":
                loss_fct = CrossEntropyLoss()
                loss = loss_fct(logits.view(-1, self.num_labels), labels.view(-1))
            elif self.config.problem_type == "multi_label_classification":
                loss_fct = BCEWithLogitsLoss()
                loss = loss_fct(logits, labels.to(logits.dtype))
        if not return_dict:
            output = (logits,) + outputs[1:]
            return ((loss,) + output) if loss is not None else output

        return SequenceClassifierOutput(
            loss=loss,
            logits=logits,
            hidden_states=outputs.hidden_states,
            attentions=outputs.attentions,
        )


class ProSSTForTokenClassification(ProSSTPreTrainedModel):
    def __init__(self, config):
        super().__init__(config)
        self.num_labels = config.num_labels

        self.prosst = ProSSTModel(config)
        self.dropout = nn.Dropout(config.hidden_dropout_prob)
        self.classifier = nn.Linear(config.hidden_size, config.num_labels)

        # Initialize weights and apply final processing
        self.post_init()

    def forward(
        self,
        input_ids: Optional[torch.Tensor] = None,
        attention_mask: Optional[torch.Tensor] = None,
        token_type_ids: Optional[torch.Tensor] = None,
        position_ids: Optional[torch.Tensor] = None,
        inputs_embeds: Optional[torch.Tensor] = None,
        labels: Optional[torch.Tensor] = None,
        output_attentions: Optional[bool] = None,
        output_hidden_states: Optional[bool] = None,
        return_dict: Optional[bool] = None,
    ) -> Union[Tuple, TokenClassifierOutput]:
        r"""
        labels (`torch.LongTensor` of shape `(batch_size, sequence_length)`, *optional*):
            Labels for computing the token classification loss. Indices should be in `[0, ..., config.num_labels - 1]`.
        """
        return_dict = (
            return_dict if return_dict is not None else self.config.use_return_dict
        )

        outputs = self.prosst(
            input_ids,
            attention_mask=attention_mask,
            token_type_ids=token_type_ids,
            position_ids=position_ids,
            inputs_embeds=inputs_embeds,
            output_attentions=output_attentions,
            output_hidden_states=output_hidden_states,
            return_dict=return_dict,
        )

        sequence_output = outputs[0]

        sequence_output = self.dropout(sequence_output)
        logits = self.classifier(sequence_output)

        loss = None
        if labels is not None:
            loss_fct = CrossEntropyLoss()
            loss = loss_fct(logits.view(-1, self.num_labels), labels.view(-1))

        if not return_dict:
            output = (logits,) + outputs[1:]
            return ((loss,) + output) if loss is not None else output

        return TokenClassifierOutput(
            loss=loss,
            logits=logits,
            hidden_states=outputs.hidden_states,
            attentions=outputs.attentions,
        )


ProSSTModel.register_for_auto_class("AutoModel")
ProSSTForMaskedLM.register_for_auto_class("AutoModelForMaskedLM")
ProSSTForSequenceClassification.register_for_auto_class(
    "AutoModelForSequenceClassification"
)
ProSSTForTokenClassification.register_for_auto_class("AutoModelForTokenClassification")
