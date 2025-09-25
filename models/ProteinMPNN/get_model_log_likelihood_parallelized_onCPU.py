import os
import torch
import random
import subprocess
import numpy as np
import pandas as pd
from argparse import ArgumentParser
from utils import load_structure, extract_coords_from_complex, get_metadata
from copy import deepcopy
from Bio.PDB import PDBParser, PDBIO
from Bio.Data import IUPACData
import shutil
from pathlib import Path
import multiprocessing
from tqdm import tqdm # For a progress bar
import hashlib

one_to_three = IUPACData.protein_letters_1to3_extended
three_to_one = IUPACData.protein_letters_3to1_extended

def setup_seed(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    # The following two lines are not needed for CPU-only runs but are kept for consistency
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def numbered_to_sequential(input_pdb, output_pdb):
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("input_structure", input_pdb)
    
    mapping_dict = {}

    new_structure = deepcopy(structure)    
    for model in new_structure:
        for chain in model:
            new_residue_id = 1 # Start renumbering residues from 1
            for residue in chain:
                original_id = (chain.id, residue.id[0], residue.id[1], residue.id[2])
                residue.id = (residue.id[0], residue.id[1]+10000, residue.id[2])
                new_id = (residue.id[0], new_residue_id, ' ')
                mapping_dict[original_id] = new_id
                new_residue_id += 1

    for original_model, new_model in zip(structure, new_structure):
        for original_chain, new_chain in zip(original_model, new_model):
            for original_residue, new_residue in zip(original_chain, new_chain):
                original_id = (original_chain.id, original_residue.id[0], original_residue.id[1], original_residue.id[2])
                new_id = mapping_dict[original_id]
                new_residue.id = new_id

    io = PDBIO()
    io.set_structure(new_structure)
    io.save(output_pdb)
    return mapping_dict

def mutate_pdb_sequence(pdb_path, mut_info, output_pdb_path):
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_path)

    for mutation in mut_info:
        original_residue = mutation[0]
        chain_id = mutation[1]
        residue_index = int(mutation[2:-1])
        mutated_residue = mutation[-1]

        chain = structure[0][chain_id]

        found = False
        for residue in chain:
            if residue.get_id()[1] == residue_index:
                resname = residue.get_resname().upper()[0] + residue.get_resname().lower()[1:]
                resname = three_to_one[resname]
                assert resname == original_residue, \
                    f"Original residue {original_residue} expected at chain {chain_id} residue {residue_index}, " \
                    f"but found {resname}"
                
                residue.resname = one_to_three[mutated_residue].upper()
                found = True
                break
        
        if not found:
            raise ValueError(f"Residue {residue_index} in chain {chain_id} not found or mismatch in the PDB file.")

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb_path)
    return output_pdb_path


def eval_mpnn_score(pdb_file, output_dir, MPNN_PATH):

    result = subprocess.run([
        "python", f"{MPNN_PATH}/protein_mpnn_run.py",
        "--pdb_path", pdb_file,
        "--out_folder", output_dir,
        "--score_only", "1",
        "--seed", "37",
        "--gpu", "-1"  # Force CPU usage
    ], capture_output=True, text=True)

    if result.returncode == 0:
        filename = os.path.splitext(os.path.basename(pdb_file))[0]
        data = f'{output_dir}/score_only/{filename}_pdb.npz'
        if os.path.exists(data):
            loaded_data = np.load(data)
            score = loaded_data['score'][0]
            return -score
        else:
            print(f"Error: Output file does not exist: {data}")
            return None
    else:
        print(f"Error running ProteinMPNN: {result.stderr}")
        return None

# --- NEW WORKER FUNCTION ---
def process_row(args_tuple):
    """
    This function processes a single row of the DataFrame.
    It's designed to be called by a separate process.
    """
    idx, row, pdb_file, native_seqs, chain_order, heavy_chain_id, light_chain_id, antigen_chains, mpnn_path, base_dir_name, name_prefix, mut_location, col_names = args_tuple
    
    # Each process gets its own unique subdirectory to prevent file collisions
    pid = os.getpid()
    process_dir = os.path.join(base_dir_name, str(pid))
    os.makedirs(process_dir, exist_ok=True)

    if mut_location == 'heavy':
        mutated_seqs = {
        heavy_chain_id: row['heavy_chain_seq'],
        light_chain_id: native_seqs[light_chain_id]
        }
    elif mut_location == 'light':
        mutated_seqs = {
            heavy_chain_id: native_seqs[heavy_chain_id],
            light_chain_id: row['light_chain_seq']
        }
    else:
        mutated_seqs = {
            heavy_chain_id: row['heavy_chain_seq'],
            light_chain_id: row['light_chain_seq']
        }
    
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

    mut_info_ = ','.join(mut_info)
    mut_info_ += ';'
    mutations = mut_info_[:-1]
    
    # pdb_name = args.name + f'_{mutations}' # Getting too long file name issue
    mut_hash = hashlib.md5(mutations.encode()).hexdigest()
    pdb_name = f"{args.name}_{mut_hash}"
    temp_pdb_path = os.path.join(process_dir, f"{pdb_name}.pdb")
    
    _ = numbered_to_sequential(pdb_file, temp_pdb_path)
    mutated_path = mutate_pdb_sequence(temp_pdb_path, mut_info, temp_pdb_path)
    
    if not Path(mutated_path).exists():
        print(f'{mutated_path} not found')
        score = None
    else:
        score = eval_mpnn_score(mutated_path, process_dir, mpnn_path)

    shutil.rmtree(process_dir)
    
    return idx, score


def main(args):
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

    dir_name = f'./tmp_proteinmpnn_{name}'
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    outpath = f'./notebooks/scoring_outputs/{name}_benchmarking_data_ProteinMPNN_scores.csv'

    df = pd.read_csv(excel_file)
    structure = load_structure(pdb_file)
    _, native_seqs = extract_coords_from_complex(structure)
    
    tasks = []
    for idx, row in df.iterrows():
        task_args = (
            idx, row, pdb_file, native_seqs, chain_order,
            heavy_chain_id, light_chain_id, antigen_chains,
            args.mpnn_path, dir_name, args.name, args.mut_location, df.columns
        )
        tasks.append(task_args)

    print(f"Starting processing of {len(df)} mutations on {args.num_cores} CPU cores...")
    
    results = []
    with multiprocessing.Pool(processes=args.num_cores) as pool:
        results = list(tqdm(pool.imap_unordered(process_row, tasks), total=len(tasks)))

    print("Processing complete. Compiling results...")
    for idx, score in results:
        if score is not None:
            df.at[idx, 'log-likelihood'] = score
        else:
            df.at[idx, 'log-likelihood'] = np.nan

    df.to_csv(outpath, header=True, index=False)
    print(f"Results saved to: {outpath}")
    shutil.rmtree(dir_name)

def parse():
    parser = ArgumentParser(description='Generate antibody scores with ProteinMPNN')
    parser.add_argument('--name', type=str, default='3gbn', help='dataset name, default: 3gbn')
    parser.add_argument('--num_cores', type=int, default=1, help='Number of CPU cores to use for parallel processing')
    parser.add_argument('--seed', type=int, default=42, help='Seed to use')
    parser.add_argument('--mpnn_path', type=str, default='./models/ProteinMPNN', help='Path to ProteinMPNN directory')
    parser.add_argument('--mut_location', type=str, default='both', choices=['heavy', 'light', 'both'], 
                        help="Where mutations are applied: 'heavy', 'light', or 'both'")
    return parser.parse_args()

if __name__ == '__main__':
    args = parse()
    setup_seed(args.seed)
    main(args)