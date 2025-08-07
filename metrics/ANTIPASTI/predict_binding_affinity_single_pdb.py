""" 
cp modified preprocessing.py
cp env_w_r.yml

renumber pdb outside of ANTIPASTI
"""

import os
import sys
import argparse
import numpy as np
import torch
import tempfile
import shutil
from Bio import PDB
from Bio.PDB import PDBParser

# ANTIPASTI
from antipasti.preprocessing.preprocessing import Preprocessing
from antipasti.utils.torch_utils import load_checkpoint

def relabel_pdb_chains(input_pdb, output_pdb, chain_mapping):
    parser = PDB.PDBParser()
    structure = parser.get_structure("structure", input_pdb)
    #chain_mapping = {'A': 'A', 'B': 'B', 'C': 'H', 'D': 'L'}
    for model in structure:
        for chain in model:
            if chain.id in chain_mapping:
                chain.id = chain_mapping[chain.id]
    
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

def generate_residue_label_array(pdb_file, heavy_chain_id, light_chain_id, output_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    labels = ['START']
    for chain_id in [heavy_chain_id, light_chain_id]:
        if chain_id not in structure[0]:
            continue
        chain = structure[0][chain_id]
        for residue in chain:
            hetfield, resseq, icode = residue.id
            if hetfield != ' ' or 'CA' not in residue:
                continue
            icode = icode.strip()
            res_id_str = f"{chain_id}{resseq:>3}{icode}".ljust(5)
            labels.append(res_id_str)
    labels.append('END')
    np.save(output_file, np.array(labels, dtype='<U5'))
    print(f'saved list_of_residues to {output_file}')

def eval_binding_affinity(pdb_path, heavy_chain, light_chain, antigen_chains, alphafold):
    # Setup temporary environment
    temp_root = tempfile.mkdtemp()
    structure_dir = os.path.join(temp_root, "structure")
    dccm_map_dir = os.path.join(temp_root, "dccm_map")
    residue_dir = os.path.join(temp_root, "list_of_residues")
    os.makedirs(structure_dir)
    os.makedirs(dccm_map_dir)
    os.makedirs(residue_dir)

    test_pdb_id = 'test_af' if alphafold else 'test'
    tmp_pdb = os.path.join(structure_dir, f"{test_pdb_id}.pdb")
    if alphafold:
        antigen_remap = {chain: chr(67 + i) for i, chain in enumerate(antigen_chains)}
        chain_mapping = {
            heavy_chain: "A",
            light_chain: "B",
            **antigen_remap
        }
        relabel_pdb_chains(pdb_path, tmp_pdb, chain_mapping)
        print(f'relabeling chain ids...:{chain_mapping}')
        heavy_chain = "A"
        light_chain = "B"
    else:
        shutil.copy(pdb_path, tmp_pdb)
    print(f'copy input pdb to {tmp_pdb}')

    # Generate residue list
    residue_file = os.path.join(residue_dir, f"{test_pdb_id}.npy")
    generate_residue_label_array(tmp_pdb, heavy_chain, light_chain, residue_file)

    # Model settings
    modes = 'all'
    n_filters = 4
    filter_size = 4
    pooling_size = 1
    n_max_epochs = 1044
    pathological = ['5omm', '5i5k', '1uwx', '1mj7', '1qfw', '1qyg', '4ffz', '3ifl', '3lrh', '3pp4', '3ru8', '3t0w', '3t0x', '4fqr', '4gxu', '4jfx', '4k3h', '4jfz', '4jg0', '4jg1', '4jn2', '4o4y', '4qxt', '4r3s', '4w6y', '4w6y', '5ies', '5ivn', '5j57', '5kvd', '5kzp', '5mes', '5nmv', '5sy8', '5t29', '5t5b', '5vag', '3etb', '3gkz', '3uze', '3uzq', '4f9l', '4gqp', '4r2g', '5c6t', '3fku', '1oau', '1oay']
    scfv = ['4gqp', '3etb', '3gkz', '3uze', '3uzq', '3gm0', '4f9l', '6ejg', '6ejm', '1h8s', '5dfw', '6cbp', '4f9p', '5kov', '1dzb', '5j74', '5aaw', '3uzv', '5aam', '3ux9', '5a2j', '5a2k', '5a2i', '3fku', '5yy4', '3uyp', '5jyl', '1y0l', '1p4b', '3kdm', '4lar', '4ffy', '2ybr', '1mfa', '5xj3', '5xj4', '4kv5', '5vyf'] 
    pathological += scfv
    model_path = f'metrics/ANTIPASTI/checkpoints/full_ags_all_modes/model_epochs_{n_max_epochs}_modes_{modes}_pool_{pooling_size}_filters_{n_filters}_size_{filter_size}.pt' 
    
    # Run preprocessing
    preprocessed_data = Preprocessing(
        data_path='metrics/ANTIPASTI/data/',
        scripts_path='metrics/ANTIPASTI/scripts',
        dccm_map_path='dccm_maps_full_ags_all/',
        modes=modes,
        pathological=pathological,
        renew_maps=False,
        renew_residues=True,
        stage='predicting',
        test_data_path=f'{temp_root}/',
        test_dccm_map_path='dccm_map/',
        test_residues_path='list_of_residues/',
        test_structure_path='structure/',
        test_pdb_id=test_pdb_id,
        alphafold=alphafold
    )
    input_shape = preprocessed_data.test_x.shape[-1]

    #print(f"loading model:{model_path}")
    model, _, _, _, _ = load_checkpoint(model_path, input_shape)
    model.eval()
 
    # Predict log Kd
    test_sample = torch.from_numpy(preprocessed_data.test_x.reshape(1, 1, input_shape, input_shape).astype(np.float32))
    pred_logKd = model(test_sample)[0].detach().numpy()[0, 0]
    #print(f'\nPredicted log(Kd): {pred_logKd:.4f}')

    # Clean up temporary files
    shutil.rmtree(temp_root)
    return pred_logKd


def main():
    args = parse_args()
    pdb_path = args.pdb_path
    heavy_chain = args.heavy_chain
    light_chain = args.light_chain
    antigen_chains = args.antigen_chains
    alphafold = args.alphafold
    pred_logKd = eval_binding_affinity(pdb_path, heavy_chain, light_chain, antigen_chains, alphafold)
    print(f'{pdb_path}: predicted binding affinity:{pred_logKd}')
    return pred_logKd

def parse_args():
    parser = argparse.ArgumentParser(description='Predict log(Kd) from antibody-antigen complex PDB')
    parser.add_argument('--pdb_path', type=str, required=True, help='Path to Chothia-numbered PDB file')
    parser.add_argument('--heavy_chain', type=str, required=True, help='Heavy chain ID')
    parser.add_argument('--light_chain', type=str, required=True, help='Light chain ID')
    parser.add_argument('--antigen_chains', type=str, required=True, help='Antigen chains, eg AB or A,B')
    parser.add_argument('--no-alphafold', dest='alphafold', action='store_false',
                        help='Disable if structure is NOT from AlphaFold (default: True)')
    parser.set_defaults(alphafold=True)
    return parser.parse_args()

if __name__ == "__main__":
    main()

