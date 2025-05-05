import os
import logging
import numpy as np
import torch
import random
import pandas as pd
from argparse import ArgumentParser
from utils import load_structure, extract_coords_from_complex, get_metadata
import shutil


def setup_seed(seed):
     torch.manual_seed(seed)
     torch.cuda.manual_seed_all(seed)
     np.random.seed(seed)
     random.seed(seed)
     torch.backends.cudnn.deterministic = True
     torch.backends.cudnn.benchmark = False


def binding_dg(pdb_path, out_dir, receptor_chain, ligand_chain):
    def exec_foldx(cmd):
        r = os.popen(cmd)
        text = r.read()
        r.close()
        return text

    pdb_dir = os.path.dirname(pdb_path)
    pdb = os.path.basename(pdb_path)
    pdb_id = pdb.split(".pdb")[0]
    rec_str, lig_str = ''.join(receptor_chain), ''.join(ligand_chain)
    cmd = f'./metrics/FoldX/FoldX --command AnalyseComplex --pdb-dir {pdb_dir} --pdb {pdb} --output-dir {out_dir} --output-file {pdb_id} --analyseComplexChains {rec_str},{lig_str}'
    exec_foldx(cmd)

    foldx_out = f'Interaction_{pdb_id}_AC.fxout'
    foldx_out = os.path.join(out_dir, foldx_out)
    with open(foldx_out, 'r') as f:
        data = f.readlines()

    header = data[8].strip().split('\t')
    data_lines = [line.strip().split('\t') for line in data[9:]]
    df = pd.DataFrame(data_lines, columns=header)
    dg = float(df.iloc[0]['Interaction Energy'])

    os.remove(foldx_out)
    os.remove(os.path.join(out_dir,f'Indiv_energies_{pdb_id}_AC.fxout'))
    os.remove(os.path.join(out_dir,f'Interface_Residues_{pdb_id}_AC.fxout'))
    os.remove(os.path.join(out_dir,f'Summary_{pdb_id}_AC.fxout'))

    return dg

def main(args):

    tmpdir = f'./tmp_foldx_{args.name}'
    os.makedirs(tmpdir, exist_ok=True)

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
    receptor_chains = info["receptor_chains"]
    ligand_chains = info["ligand_chains"]

    # setup outputs 
    mutant_pdb_dir = f"./metrics/mutant_structure/{name}"
    outpath = f'./notebooks/scoring_outputs/{name}_benchmarking_data_FoldX_scores.csv'
    
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
        dg = binding_dg(mutated_path, tmpdir, receptor_chains, ligand_chains)
        df.at[idx, 'dg'] = dg
        
    df.to_csv(outpath)
    shutil.rmtree(tmpdir)

def parse():
    parser = ArgumentParser(description='Generate antibody')
    parser.add_argument('--name', type=str, default='1mlc')
    parser.add_argument('--seed', type=int, default=1, help='Seed to use')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse()
    setup_seed(args.seed)
    main(args)