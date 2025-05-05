import argparse
import numpy as np
import pandas as pd
from utils.random_seed import setup_seed
from data.pdb_utils import VOCAB, AgAbComplex, Residue, Protein, Peptide
from data.dataset import _generate_chain_data 
from models.dyMEAN.dyMEANOpt_model import dyMEANOptModel

import torch
from copy import deepcopy
from torch.utils.data import DataLoader

from pathlib import Path
from prepare_data import get_metadata, get_structure_data

def to_cplx(ori_cplx, ab_x, ab_s) -> AgAbComplex: # copied from dyMEAN.api.generate
    heavy_chain, light_chain = [], []
    chain = None
    for residue, residue_x in zip(ab_s, ab_x):
        residue = VOCAB.idx_to_symbol(residue)
        if residue == VOCAB.BOA:
            continue
        elif residue == VOCAB.BOH:
            chain = heavy_chain
            continue
        elif residue == VOCAB.BOL:
            chain = light_chain
            continue
        if chain is None:  # still in antigen region
            continue
        coord, atoms = {}, VOCAB.backbone_atoms + VOCAB.get_sidechain_info(residue)

        for atom, x in zip(atoms, residue_x):
            coord[atom] = x
        chain.append(Residue(
            residue, coord, _id=(len(chain), ' ')
        ))
    heavy_chain = Peptide(ori_cplx.heavy_chain, heavy_chain)
    light_chain = Peptide(ori_cplx.light_chain, light_chain)
    for res, ori_res in zip(heavy_chain, ori_cplx.get_heavy_chain()):
        res.id = ori_res.id
    for res, ori_res in zip(light_chain, ori_cplx.get_light_chain()):
        res.id = ori_res.id

    peptides = {
        ori_cplx.heavy_chain: heavy_chain,
        ori_cplx.light_chain: light_chain
    }
    antibody = Protein(ori_cplx.pdb_id, peptides)
    cplx = AgAbComplex(
        ori_cplx.antigen, antibody, ori_cplx.heavy_chain,
        ori_cplx.light_chain, skip_epitope_cal=True,
        skip_validity_check=True
    )
    cplx.cdr_pos = ori_cplx.cdr_pos
    return cplx

class ComplexSummary: # copied from dyMEAN.api.optimize
    def __init__(self, pdb, heavy_chain, light_chain, antigen_chains) -> None:
        self.pdb = pdb
        self.heavy_chain = heavy_chain
        self.light_chain = light_chain
        self.antigen_chains = antigen_chains

class Dataset(torch.utils.data.Dataset): # copied from dyMEAN.api.optimize
    def __init__(self, pdb, heavy_chain, light_chain, antigen_chains, num_residue_changes, cdr='H3', pos_prob=None) -> None:
        super().__init__()
        self.pdb = pdb
        self.num_residue_changes = num_residue_changes
        self.cdr = cdr
        cplx = AgAbComplex.from_pdb(pdb, heavy_chain, light_chain, antigen_chains)
        
        # generate antigen data
        ag_residues = []
        for residue, chain, i in cplx.get_epitope():
            ag_residues.append(residue)
        ag_data = _generate_chain_data(ag_residues, VOCAB.BOA)
        

        hc, lc = cplx.get_heavy_chain(), cplx.get_light_chain()
        hc_residues, lc_residues = [], []

        # generate heavy chain data
        for i in range(len(hc)):
            hc_residues.append(hc.get_residue(i))
        hc_data = _generate_chain_data(hc_residues, VOCAB.BOH)
        

        # generate light chain data
        for i in range(len(lc)):
            lc_residues.append(lc.get_residue(i))
        lc_data = _generate_chain_data(lc_residues, VOCAB.BOL)
        

        data = { key: np.concatenate([ag_data[key], hc_data[key], lc_data[key]], axis=0) \
                 for key in hc_data}  # <X, S, residue_pos>
        
        cmask = [0 for _ in ag_data['S']] + [0] + [1 for _ in hc_data['S'][1:]] + [0] + [1 for _ in lc_data['S'][1:]]
        smask = [0 for _ in range(len(ag_data['S']) + len(hc_data['S']) + len(lc_data['S']))]
        data['cmask'], data['smask'] = cmask, smask

        self.cdr_idxs = []
        cdrs = [self.cdr] if type(self.cdr) == str else self.cdr
        for cdr in cdrs:
            cdr_range = cplx.get_cdr_pos(cdr)
            offset = len(ag_data['S']) + 1 + (0 if cdr[0] == 'H' else len(hc_data['S']))
            for idx in range(offset + cdr_range[0], offset + cdr_range[1] + 1):
                self.cdr_idxs.append(idx)
        self.pos_prob = []
        if pos_prob is not None: # add positional probabilities
            for cdr in cdrs:
                cdr_range = cplx.get_cdr_pos(cdr)
                for idx in range(cdr_range[0], cdr_range[1]+1):
                    self.pos_prob.append(pos_prob[idx])

        self.cplx, self.data = cplx, data

    def __getitem__(self, idx):
        data = deepcopy(self.data)
        num_residue_chain = min(self.num_residue_changes[idx], len(self.cdr_idxs))
        if num_residue_chain <= 0:
            num_residue_chain = np.random.randint(1, len(self.cdr_idxs))
        if len(self.pos_prob) > 0:
            prob_sum = np.sum(self.pos_prob)
            normalized_prob = self.pos_prob / prob_sum #probabilities needs to be sumed to 1
            mask_idxs = np.random.choice(self.cdr_idxs, size=num_residue_chain, replace=False, p=normalized_prob)
        else: 
            mask_idxs = np.random.choice(self.cdr_idxs, size=num_residue_chain, replace=False)
        for i in mask_idxs:
            data['smask'][i] = 1
        return data

    def __len__(self):
        return len(self.num_residue_changes)

    @classmethod
    def collate_fn(cls, batch):
        keys = ['X', 'S', 'smask', 'cmask', 'residue_pos']
        types = [torch.float, torch.long, torch.bool, torch.bool, torch.long]
        res = {}
        for key, _type in zip(keys, types):
            val = []
            for item in batch:
                val.append(torch.tensor(item[key], dtype=_type))
            res[key] = torch.cat(val, dim=0)
        lengths = [len(item['S']) for item in batch]
        res['lengths'] = torch.tensor(lengths, dtype=torch.long)
        return res

def main(args):
    # load model 
    device = torch.device('cpu' if args.gpu == -1 else f'cuda:{args.gpu}')
    loaded_model = torch.load(args.ckpt, map_location='cpu')
    state_dict = loaded_model.state_dict()
    model = dyMEANOptModel(embed_size=64, hidden_size=128, n_channel=VOCAB.MAX_ATOM_NUMBER,
                           num_classes=VOCAB.get_num_amino_acid_type(), mask_id=VOCAB.get_mask_idx())
    model.load_state_dict(state_dict)    
    model.to(device)
    model.eval()

    # load affinity data
    if args.name == "aayl49_ml":
        name = "aayl49_ML"
    else:
        name = args.name
    pdb_info = get_metadata()[name] 
    df = pd.read_csv(pdb_info["affinity_data"][0])

    # setup outputs
    fout_path = f'./notebooks/scoring_outputs/{args.name}_benchmarking_data_dyMEAN_scores.csv'
    
    # calculate sequence negative log-likelihood for mutant sequence
    fnll_list, nll_list = [], []
    for idx, row in df.iterrows():
        pdb_data = f"{name}_{idx}.pdb"
        get_structure_data(pdb_info["pdb_path"], pdb_info["heavy_chain"], row["mut_heavy_chain_seq"], pdb_data, renumber=True)
        
        
        cplx_summary = ComplexSummary(pdb=pdb_data,
                                      heavy_chain=pdb_info["heavy_chain"],
                                      light_chain=pdb_info["light_chain"],
                                      antigen_chains=pdb_info["antigen_chains"])
        
        num_residue_changes = [1] # one datapoint
        cdr_type = "H3"
        pos_prob = None
        batch_size = 1
        num_workers = 1
        dataset = Dataset(**cplx_summary.__dict__, num_residue_changes=num_residue_changes, cdr=cdr_type, pos_prob=pos_prob)
        dataloader = DataLoader(dataset, batch_size=batch_size,
                                num_workers=num_workers,
                                collate_fn=Dataset.collate_fn)
        batch = [data for data in dataloader][0] # only one datapoint, no need to iterate 
        filtered_batch = {k: batch[k].to(device) if hasattr(batch[k], 'to') else batch[k] 
                     for k in batch 
                     if k in ['X', 'S', 'cmask', 'smask', 'residue_pos', 'lengths']}
        
        #for k, v in filtered_batch.items():
        #    print(f"{k}: {v.shape}")
        fgen_X, fgen_S, fsnll = model.get_seq_nlls(**filtered_batch, fixbb=True)
        gen_X, gen_S, snll = model.get_seq_nlls(**filtered_batch, fixbb=False)
        fnll_list.append(fsnll.detach().cpu().item())
        nll_list.append(snll.detach().cpu().item())
        #print(f"idx={idx}: {imgt_pdb} | fixed snll={fsnll} | snll={snll}")

        # Clear GPU memory
        del filtered_batch  # Delete the filtered batch to free memory
        torch.cuda.empty_cache()  # Clear unused memory from the cache

        # Remove tmp files
        Path(pdb_data).unlink()
        
    
    # dyMEAN reports sequence negative log-likelihood, but we want log-likelihood
    df.loc[:, 'log-likelihood (fixed backbone)'] = [-fnll for fnll in fnll_list]
    df.loc[:, 'log-likelihood'] = [-nll for nll in nll_list]
    df.to_csv(fout_path, header=True, index=False)
    return df 

def parse():
    parser = argparse.ArgumentParser(description='Calculate sequence log-likelihood from dyMEAN')
    parser.add_argument('--name', type=str, default='3gbn')
    parser.add_argument('--gpu', type=int, default=4, help='-1 for cpu')
    parser.add_argument('--seed', type=int, default=42, help='seed number')
    parser.add_argument('--ckpt', type=str, default='./models/dyMEAN/checkpoints/cdrh3_opt.ckpt')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse()
    setup_seed(args.seed)
    main(args)