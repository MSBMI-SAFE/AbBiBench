import os, sys
# Add the parent directory of the current script to the Python path
#current_dir = os.path.dirname(os.path.realpath(__file__))
#parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
#sys.path.append(parent_dir)

import json
import glob
from time import time
from Bio import PDB
import warnings
import configparser
import pandas as pd
import random

AA_CODE = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP',
           'C': 'CYS', 'E': 'GLU', 'Q': 'GLN', 'G': 'GLY',
           'H':'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
           'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER',
           'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}
def read_config():
    config = configparser.ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__), '..', 'config', 'config.ini'))
    return config

def built_in_pkgs():
    import os
    import json
    import argparse
    import numpy as np
    import pandas as pd
    from tqdm import tqdm
    from time import time
    from tqdm.contrib.concurrent import process_map
    from datetime import datetime		
    return os, json, argparse, np, pd, tqdm, time, process_map, datetime

def get_aa_code(code='one_letter'):
    """
    :param use: string representing which type of amino acid coding you want to return, option=['one_letter', 'three_letter']
    :return aa_list: list of amino acid codes
    """
    if code == 'one_letter':
        return list(AA_CODE.keys())
    else:
        return list(AA_CODE.values())

def load_pdb(pdb_path):
    """
    :param pdb_path: string representing path to the pdb file
    :return structure: BioPythnon's structure object
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parser = PDB.PDBParser(QUIET=True)  # Use QUIET=True to suppress warnings
        structure = parser.get_structure('structure', pdb_path)
    return structure
def get_residue_number(pdb_path, chain_id):
    """
    :param pdb_path: string representing path to the pdb file
    :param chain_id: string representing chain id where to look for residue
    :return resd_dict: dictionary storing index: (residue_name, residue_number)
    :return df: dataframe with headers=['chain_id', 'amino acid', 'amino acid number', 'residue', 'residue number']
    """
    # skip non-standard residues, e.g., HOH
    valid_codes = get_aa_code(code='three_letter')
        
    # get one-letter code
    aa_dict = {val: key for key, val in AA_CODE.items()}
    resd_list = []
    structure = load_pdb(pdb_path)
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                residues = list(chain.get_residues())
                for idx, residue in enumerate(residues):
                    # residue.id[1] is residue number
                    # residue.resname is residue in three-letter code, e.g., ALA, GLU
                    if residue.resname in valid_codes:
                        # if insertion code (icode) is not empty
                        # it will lead to multiple residues in the same positions!!!
                        #resd_list.append( (chain_id, residue.resname, (residue.id[1], residue.id[2])) )
                        resd_list.append( (chain_id, aa_dict[residue.resname], idx+1, residue.resname, str(residue.id[1])+residue.id[2]) )
    # list 2 df
    df = pd.DataFrame(resd_list, columns=['chain_id', 'amino acid', 'amino acid number', 'residue', 'residue number'])
    df['residue number'] = df['residue number'].str.strip()
    return df

def sub_pdb(pdb_path, chain_id, pos, mut_aa):
    """
    :param pdb_path: string representing path to the pdb file
    :param chain_id: string representing chain id in thd pdb fil
    :param pos: integer representing position to be substituted with mut_aa
    :param mut_aa: string representing the three-letter code of a amino acid
    :return mut_path: string representing path to mutant pdb file
    """
    # setup output filename
    basename = os.path.basename(pdb_path).split('.pdb')[0]
    mut_path = f'{os.path.dirname(pdb_path)}/{basename}_sub{mut_aa}_{pos}.pdb'

    # load pdb
    structure = load_pdb(pdb_path)

    # collect chain_ids
    chain_list = [chain.id for model in structure for chain in model]
    if chain_id not in chain_list:
        raise ValueError(f'{chain_id} not found!!! Found chains={chain_list}')

    # substitue the residue on given position with ALA
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                residues = list(chain.get_residues())
                for idx, residue in enumerate(residues):
                    # residue.id[1] is residue number
                    # residue.resname is residue in three-letter code, e.g., ALA, GLU
                    ori_aa = residue.resname
                    mut_aa = mut_aa
                    if residue.id[1] == int(pos):
                        residue.resname = mut_aa
                        #print(f'position={pos}: original residue={ori_aa} | mutant={residue.resname}')
    # save to pdb file
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(mut_path)
    return mut_path

def clean_up(path, substr=None):
    if substr is not None:
        pattern = os.path.join(path, f'*{substr}*')
        file_list = glob.glob(pattern)
    else:
        file_list = glob.glob(os.path.join(path, '*'))

    for f in file_list:
        try:
            if os.path.isfile(f):
                os.remove(f)
            elif os.path.isdir(f):
                os.rmdir(f)
        except OSError as e:
            print(f"Error deleting {f}: {e}")

def json2dict(json_path):
    with open(json_path, 'r') as fin:
        items = [json.loads(line) for line in fin.read().strip().split('\n')]
    info = { item['pdb']: item for item in items }
    return info

def get_time_sign() -> str:
    time_note = time()
    unique_id = round((time_note - int(time_note)) * 1_000_000)  # Microseconds
    pid = os.getpid()  # Get the process ID
    rand_num = random.randint(0, 9999)  # Small random number
    return f"{pid}_{unique_id}_{rand_num}"

def create_dir(new_dir):
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
        print(f'creating folder: {new_dir}')
    print(f'{new_dir} already exists, will overwirte')

def remove_chain(input_pdb, output_pdb, remove_chains):
    remove_chains = set(remove_chains)
    kept_chains = set()
    
    # First pass: identify kept chains
    with open(input_pdb, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain = line[21]
                if chain not in remove_chains:
                    kept_chains.add(chain)
    
    # Second pass: write kept chains to output file
    with open(input_pdb, 'r') as f_in, open(output_pdb, 'w') as f_out:
        for line in f_in:
            if line.startswith(('ATOM', 'HETATM', 'TER')):
                if line[21] not in remove_chains:
                    f_out.write(line)
            elif line.startswith(('HEADER', 'TITLE', 'REMARK')):
                f_out.write(line)
        f_out.write('END\n')

    print(f"Chains removed: {', '.join(sorted(remove_chains))}")
    print(f"Chains kept: {', '.join(sorted(kept_chains))}")
    print(f"Modified structure written to {output_pdb}")