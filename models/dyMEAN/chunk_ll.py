import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from utils.random_seed import setup_seed
from data.pdb_utils import VOCAB, AgAbComplex, Residue, Protein, Peptide
from data.dataset import _generate_chain_data 
from models.dyMEAN.dyMEANOpt_model import dyMEANOptModel

import torch
from copy import deepcopy
from torch.utils.data import DataLoader
from prepare_data import get_metadata, renumber_pdb
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

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

def process_chunk_parallel(chunk_data):
    """Process a chunk on a specific GPU"""
    chunk_df, pdb_info, gpu_id, args = chunk_data
    
    # Set GPU for this process
    device = torch.device('cpu' if gpu_id == -1 else f'cuda:{gpu_id}')
    
    # Load model (each process loads its own copy)
    loaded_model = torch.load(args.ckpt, map_location='cpu')
    state_dict = loaded_model.state_dict()
    model = dyMEANOptModel(embed_size=64, hidden_size=128, n_channel=VOCAB.MAX_ATOM_NUMBER,
                           num_classes=VOCAB.get_num_amino_acid_type(), mask_id=VOCAB.get_mask_idx())
    model.load_state_dict(state_dict)    
    model.to(device)
    model.eval()

    # Process chunk
    fnll_list, nll_list = [], []
    for idx, row in chunk_df.iterrows():
        ori_pdb = f"data/mutant_structure/mutant_structure_{pdb_info['pdb']}/variant_{idx}.pdb"
        pdb_data = f"data/mutant_structure/mutant_structure_{pdb_info['pdb']}/variant_{idx}.imgt.pdb"
        renumber_pdb(ori_pdb, pdb_data, scheme='imgt', mute=True)
        
        cplx_summary = ComplexSummary(pdb=pdb_data,
                                      heavy_chain=pdb_info["heavy_chain"],
                                      light_chain=pdb_info["light_chain"],
                                      antigen_chains=pdb_info["antigen_chains"])
        
        # Create dataset and dataloader
        dataset = Dataset(**cplx_summary.__dict__, num_residue_changes=[1], cdr="H3", pos_prob=None)
        dataloader = DataLoader(dataset, batch_size=1, num_workers=1, collate_fn=Dataset.collate_fn)
        batch = next(iter(dataloader))
        
        filtered_batch = {k: batch[k].to(device) if hasattr(batch[k], 'to') else batch[k] 
                         for k in batch 
                         if k in ['X', 'S', 'cmask', 'smask', 'residue_pos', 'lengths']}
        
        with torch.inference_mode(): 
            fgen_X, fgen_S, fsnll = model.get_seq_nlls(**filtered_batch, fixbb=True)
            gen_X, gen_S, snll = model.get_seq_nlls(**filtered_batch, fixbb=False)
            
        fnll_list.append(fsnll.detach().cpu().item())
        nll_list.append(snll.detach().cpu().item())
        
        # Clear GPU memory
        del filtered_batch
        torch.cuda.empty_cache()

        # Remove tmp files
        Path(pdb_data).unlink()
        
    
    # Add results to chunk
    chunk_df = chunk_df.copy()
    chunk_df.loc[:, 'log-likelihood (fixed backbone)'] = [-fnll for fnll in fnll_list]
    chunk_df.loc[:, 'log-likelihood'] = [-nll for nll in nll_list]
    
    return chunk_df

def main(args):
    """Improved main function with true parallel processing"""
    
    # Load data
    if args.name == "aayl49_ml":
        name = "aayl49_ML"
    else:
        name = args.name
    pdb_info = get_metadata()[name] 
    if args.affinity_data:
        affinity_data = args.affinity_data
    else:
        affinity_data = pdb_info["affinity_data"][0]
    df = pd.read_csv(affinity_data)
    
    # Split dataframe into chunks for each GPU
    prefix = args.prefix
    num_gpus = len(args.gpu)
    chunk_size = len(df) // num_gpus
    chunks = []
    
    for i, gpu_id in enumerate(args.gpu):
        start_idx = i * chunk_size
        if i == num_gpus - 1:  # Last chunk gets remaining rows
            end_idx = len(df)
        else:
            end_idx = (i + 1) * chunk_size
            
        chunk_df = df.iloc[start_idx:end_idx].copy()
        chunks.append((chunk_df, pdb_info, gpu_id, args))
    
    # Process chunks in parallel
    print(f"Processing {len(df)} samples across {num_gpus} GPUs...")
    
    with ProcessPoolExecutor(max_workers=num_gpus) as executor:
        # Submit all chunks
        futures = {executor.submit(process_chunk_parallel, chunk): i 
                  for i, chunk in enumerate(chunks)}
        
        results = []
        for future in as_completed(futures):
            chunk_idx = futures[future]
            try:
                result_df = future.result()
                results.append((chunk_idx, result_df))
                print(f"Completed chunk {chunk_idx + 1}/{num_gpus}")
            except Exception as exc:
                print(f'Chunk {chunk_idx} generated an exception: {exc}')
    
    # Merge results in original order
    results.sort(key=lambda x: x[0])  # Sort by chunk index
    merged_df = pd.concat([result[1] for result in results], axis=0, ignore_index=True)
    
    # Save results
    fout_path = f'./notebooks/scoring_outputs/{prefix}_benchmarking_data_dyMEAN_scores.csv'
    merged_df.to_csv(fout_path, header=True, index=False)
    
    return merged_df

def parse():
    parser = argparse.ArgumentParser(description='Calculate sequence log-likelihood from dyMEAN using parallel processing')
    parser.add_argument('--name', type=str, default='3gbn', help='dataset name, default: 3gbn')
    parser.add_argument('--affinity_data', type=str, help='Path to binding_score csv file')
    parser.add_argument('--gpu', type=lambda s: [int(x) for x in s.split(',')], default=[0, 1], 
                       help='comma-separated list of GPU ids, e.g., "0,1,2,3"')
    parser.add_argument('--seed', type=int, default=42, help='seed number')
    parser.add_argument('--prefix', type=str, default='3gbn', help='prefix of the output scores.csv, default: 3gbn')
    parser.add_argument('--ckpt', type=str, default='./models/dyMEAN/checkpoints/cdrh3_opt.ckpt')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse()
    setup_seed(args.seed)
    
    print("="*50)
    print("PARALLEL DyMEAN SCORING")
    print("="*50)
    print(f"GPUs: {args.gpu}")
    print(f"Dataset: {args.name}")
    print(f"Seed: {args.seed}")
    print(f"Output prefix: {args.prefix}")
    print("="*50)

    start_time = datetime.now()
    main(args)
    end_time = datetime.now()
    print(f"Finished in {end_time-start_time}")