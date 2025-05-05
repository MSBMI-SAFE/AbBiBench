import os
import sys
import json
import copy
import argparse
import torch
import numpy as np
import pandas as pd
import random

from MEAN_.data.dataset import CustomEquiAACDataset
from MEAN_.models.MCAttGNN.mc_att_model import MyEfficientMCAttModel
from MEAN_.data import VOCAB
from utils import load_structure, extract_coords_from_complex, get_metadata

def setup_seed(seed):
     torch.manual_seed(seed)
     torch.cuda.manual_seed_all(seed)
     np.random.seed(seed)
     random.seed(seed)
     torch.backends.cudnn.deterministic = True
     torch.backends.cudnn.benchmark = False

def main(args):

    # load metadata
    if args.name == "aayl49_ml":
        name = "aayl49_ML"
    else:
        name = args.name
    info = get_metadata()[name] 
    pdb_id = info["pdb"].split('_')[0]
    excel_file = info["affinity_data"][0]
    pdb_file = info["pdb_path"]
    heavy_chain_id = info["heavy_chain"]
    light_chain_id = info["light_chain"]
    antigen_chains = info["antigen_chains"]
    chain_order = info["chain_order"]

    # setup outputs
    json_path = f'./models/MEAN/data/{pdb_id}/pdb_info.json'
    os.makedirs(os.path.dirname(json_path), exist_ok=True)
    info['pdb_data_path'] = info['pdb_path']
    with open(json_path, "w") as json_file:
        json.dump(info, json_file) #, indent=4)

    save_dir = f'./models/MEAN/data/{pdb_id}/preprocessed'
    outpath = f'./notebooks/scoring_outputs/{name}_benchmarking_data_MEAN_scores.csv'
    
    device = torch.device('cpu' if args.gpu==-1 else f'cuda:{args.gpu}')
    test_set = CustomEquiAACDataset(json_path, save_dir=save_dir)

    complex = test_set.data[0]
    n_hs = len(complex.get_heavy_chain())
    n_ls = len(complex.get_light_chain())
    offset = {f'{heavy_chain_id}':1, f'{light_chain_id}':n_hs+2}
    n_a = 0
    for i, c in enumerate(antigen_chains):
        offset[f'{c}'] = n_hs + n_ls + n_a + 3
        n_a += len(test_set.data[0].get_antigen_chains()[i]) 

    X = test_set[0]['X'].to(device)
    S = test_set[0]['S']
    n_tokens = test_set[0]['S'].shape[0]
    offsets = torch.tensor([0,n_tokens]).to(device)

    # load affinity data
    df = pd.read_csv(excel_file)

    sys.path.append('./models/MEAN/MEAN_/')
    loaded_model = torch.load(args.ckpt, map_location='cpu')
    state_dict = loaded_model.state_dict()    

    n_channel = test_set[0]['X'].shape[1]
    model = MyEfficientMCAttModel(
        embed_size=64, hidden_size=128, n_channel=n_channel, n_edge_feats=1,
        n_layers=3, alpha=0.8,
        n_iter=3
    )
    model.load_state_dict(state_dict)    
    model.to(device)
    model.eval()

    mut_S = S.clone()
    total_run = len(df)
    current_run = 0
    for idx, row in df.iterrows():
        current_run += 1
        progress = (current_run / total_run) * 100
        print(f"{progress:.2f}% ({current_run}/{total_run})")

        structure = load_structure(pdb_file)
        coords, native_seqs = extract_coords_from_complex(structure)

        mutated_seqs = {}
        mutated_seqs[heavy_chain_id] = row['mut_heavy_chain_seq']
        mutated_seqs[light_chain_id] = native_seqs[light_chain_id]
        for c in antigen_chains:
            mutated_seqs[c] = native_seqs[c]

        mutated_index = []
        for chain in chain_order:
            native_seq = native_seqs[chain]
            mut_seq = mutated_seqs[chain]
            mutations = [(i+1, native_seq[i], mut_seq[i]) for i in range(len(native_seq)) if native_seq[i] != mut_seq[i]]
            for single_mutation in mutations:
                pos, wt, mt = single_mutation
                ind = pos + offset[chain] - 1
                assert wt == VOCAB.idx_to_symbol(S[ind])
                mut_S[ind] = VOCAB.symbol_to_idx(mt)
                mutated_index.append(ind)

        mask = torch.zeros(n_tokens, dtype=torch.bool)
        mask[mutated_index] = True

        mut_S = mut_S.to(device)
        with torch.no_grad():
            _, ll_complex_fixed = model.generate_score(copy.deepcopy(X), copy.deepcopy(mut_S), copy.deepcopy(offsets), copy.deepcopy(mask), update_structure=False)            
            _, ll_complex = model.generate_score(copy.deepcopy(X), copy.deepcopy(mut_S), copy.deepcopy(offsets), copy.deepcopy(mask), update_structure=True)            

        df.at[idx, 'log-likelihood (fixed backbone)'] = ll_complex_fixed
        df.at[idx, 'log-likelihood'] = ll_complex

    df.to_csv(outpath, header=True, index=False)

def parse():
    parser = argparse.ArgumentParser(description='generation by MEAN')
    parser.add_argument('--name', type=str, default='3gbn', help='dataset name, default: 3gbn')
    parser.add_argument('--gpu', type=int, default=4, help='-1 for cpu')
    parser.add_argument('--seed', type=int, default=1, help='Batch size')
    parser.add_argument('--ckpt', type=str, default='./models/MEAN/MEAN_/checkpoints/ckpt/opt_cdrh3_mean.ckpt')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse()
    setup_seed(args.seed)
    main(args)