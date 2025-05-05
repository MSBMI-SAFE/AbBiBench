from transformers import PretrainedConfig

class ProSSTConfig(PretrainedConfig):
    model_type = "ProSST"

    def __init__(
        self,
        token_dropout=True,
        mlm_probability=0.15,
        vocab_size=1024,
        type_vocab_size=0,
        ss_vocab_size=0,
        hidden_size=768,
        num_hidden_layers=12,
        num_attention_heads=12,
        intermediate_size=3072,
        hidden_act="gelu",
        hidden_dropout_prob=0.1,
        attention_probs_dropout_prob=0.1,
        mask_token_id=24,
        initializer_range=0.02,
        layer_norm_eps=1e-7,
        pad_token_id=0,
        position_biased_input=False,
        pooler_dropout=0,
        pooler_hidden_act="gelu",
        pos_att_type=None,
        position_embedding_type="relative",
        max_position_embeddings=1024,
        max_relative_positions=-1,
        relative_attention=False,
        pooling_head="mean",
        scale_hidden=1,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.token_dropout = token_dropout
        self.mlm_probability = mlm_probability
        self.hidden_size = hidden_size
        self.num_hidden_layers = num_hidden_layers
        self.num_attention_heads = num_attention_heads
        self.intermediate_size = intermediate_size
        self.hidden_act = hidden_act
        self.hidden_dropout_prob = hidden_dropout_prob
        self.attention_probs_dropout_prob = attention_probs_dropout_prob
        self.max_position_embeddings = max_position_embeddings
        self.type_vocab_size = type_vocab_size
        self.ss_vocab_size = ss_vocab_size
        self.initializer_range = initializer_range
        self.relative_attention = relative_attention
        self.max_relative_positions = max_relative_positions
        self.pad_token_id = pad_token_id
        self.position_biased_input = position_biased_input
        self.mask_token_id = mask_token_id
        self.position_embedding_type = position_embedding_type
        self.pooling_head = pooling_head
        self.scale_hidden = scale_hidden

        # Backwards compatibility
        if type(pos_att_type) == str:
            pos_att_type = [x.strip() for x in pos_att_type.lower().split("|")]

        self.pos_att_type = pos_att_type
        self.vocab_size = vocab_size
        self.layer_norm_eps = layer_norm_eps

        self.pooler_hidden_size = kwargs.get("pooler_hidden_size", hidden_size)
        self.pooler_dropout = pooler_dropout
        self.pooler_hidden_act = pooler_hidden_act
        
ProSSTConfig.register_for_auto_class()