import os
import json
import numpy as np
import torch
import random
import pandas as pd
from argparse import ArgumentParser
from utils import load_structure, extract_coords_from_complex, get_metadata
from interface_extractor import InterfaceExtractor
from freesasa_wrapper import SASA


def setup_seed(seed):
     torch.manual_seed(seed)
     torch.cuda.manual_seed_all(seed)
     np.random.seed(seed)
     random.seed(seed)
     torch.backends.cudnn.deterministic = True
     torch.backends.cudnn.benchmark = False

def eval_epitope_sasa(pdb_path, antigen_chain, antibody_chain, contact_dist):
    ie = InterfaceExtractor(pdb_path, antigen_chain, antibody_chain, num_epitope_residues=96)
    interface_df = ie.extract_binding_epitope()
    filtered_interface_df = interface_df[interface_df['distance'] <= contact_dist]
    epitope_pos = filtered_interface_df['epitope_position'].values.tolist()

    fsa = SASA(pdb_path)
    sasa = fsa.get_sasa(chain_id=antigen_chain)
    ttl_sasa = sasa[sasa['resd_number'].astype('int').isin(epitope_pos)]['relative_mainchain'].sum()
    return ttl_sasa

def main(args):
    # load metadata
    if args.name == "aayl49_ml":
        name = "aayl49_ML"
    else:
        name = args.name
    info = get_metadata()[name]
    excel_file = info["affinity_data"][0]
    pdb_file = info["pdb_path"]
    heavy_chain_id = info["heavy_chain"]
    light_chain_id = info["light_chain"]
    antigen_chains = info["antigen_chains"]
    chain_order = info["chain_order"]
    pdb_name_offset = info["pdb"]
    epitope_chain = info["epitope_chain"]
    paratope_chain = info["paratope_chain"]

    # setup outputs
    mutant_pdb_dir = f"./metrics/mutant_structure/{name}"
    outpath = f'./notebooks/scoring_outputs/{name}_benchmarking_data_epitopeSA_scores.csv'

    # load affinity data
    df = pd.read_csv(excel_file)

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

        mut_info = []
        for chain in chain_order:
            native_seq = native_seqs[chain]
            mut_seq = mutated_seqs[chain]
            mutations = [(i+1, native_seq[i], mut_seq[i]) for i in range(len(native_seq)) if native_seq[i] != mut_seq[i]]

            for single_mutation in mutations:
                pos, wt, mt = single_mutation
                mut_info.append(f'{wt}{chain}{pos}{mt}')                
        mut_info = ','.join(mut_info)
        mut_info += ';'

        mutations = mut_info[:-1]

        pdb_name = pdb_name_offset + f'_{mutations}'
        mutated_path = os.path.join(mutant_pdb_dir, f'{pdb_name}.pdb') 
        sasa = eval_epitope_sasa(mutated_path, epitope_chain, paratope_chain, args.contact_dist)
        df.at[idx, 'EpitopeSASA (mut)'] = sasa
    
    df.to_csv(outpath, header=True, index=False)

def parse():
    parser = ArgumentParser(description='Calculate Epitope Surface Area')
    parser.add_argument('--name', type=str, default='1mlc')
    parser.add_argument('--seed', type=int, default=1, help='Seed to use')
    parser.add_argument('--contact_dist', type=float, default=5.0, help='Float defining distace between contacting epitope and paratope')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse()
    setup_seed(args.seed)
    main(args)