import os
import pandas as pd
from argparse import ArgumentParser
from utils import load_structure, extract_coords_from_complex, get_metadata
from Bio.PDB import PDBParser, PDBIO
from copy import deepcopy
import shutil


def numbered_to_sequential(input_pdb, output_pdb):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("input_structure", input_pdb)
    
    mapping_dict = {}  # To store the mapping: (original_chain, original_resseq, original_icode) -> new_resseq

    new_structure = deepcopy(structure)    
    for model in new_structure:
        for chain in model:
            new_residue_id = 1 # Start renumbering residues from 1
            for residue in chain:
                # Get original residue info
                original_id = (chain.id, residue.id[0], residue.id[1], residue.id[2])
                residue.id = (residue.id[0], residue.id[1]+10000, residue.id[2])

                # Renumber residue
                new_id = (residue.id[0], new_residue_id, ' ')
                
                # residue.id = new_id  # Update only resseq

                # Store mapping
                mapping_dict[original_id] = new_id #new_residue_id
                new_residue_id += 1  # Increment residue ID sequentially

    # mapping만 먼저 만들어두고 그 뒤에 replace를 해야지 새로운 new_residue_id가 original_id의 index보다 빨라지면 index가 겹쳐서 발생하는 문제를 해결할 수 있음
    # new_structure = deepcopy(structure)
    for original_model, new_model in zip(structure, new_structure):
        for original_chain, new_chain in zip(original_model, new_model):
            for original_residue, new_residue in zip(original_chain, new_chain):
                original_id = (original_chain.id, original_residue.id[0], original_residue.id[1], original_residue.id[2])
                new_id = mapping_dict[original_id]
                new_residue.id = new_id

    # residue.id = new_id  # Update only resseq
    # Write renumbered structure to a new file
    io = PDBIO()
    io.set_structure(new_structure)
    io.save(output_pdb)

    return mapping_dict


def build_model(pdb_path, mutations, outname, outdir, tmpdir):
    def exec_foldx(cmd):
        r = os.popen(cmd)
        text = r.read()
        r.close()
        return text

    # Specify the output file path
    mut_path = os.path.join(tmpdir, 'individual_list.txt')

    # Write the content to the file
    with open(mut_path, 'w') as file:
        file.write(mutations.strip())

    pdb_dir = os.path.dirname(pdb_path)
    pdb = os.path.basename(pdb_path)
    pdb_id = pdb[:-4]

    cmd = f'./metrics/FoldX/FoldX --command BuildModel --pdb-dir {pdb_dir} --pdb {pdb} --mutant-file {mut_path} --output-dir {tmpdir}'# --output-file {outname}'
    exec_foldx(cmd)

    foldx_out = os.path.join(tmpdir, f'{pdb_id}_1.pdb')
    move_path = os.path.join(outdir, f'{outname}.pdb')
    shutil.move(foldx_out, move_path)

    # remove
    # os.remove(os.path.join(tmpdir, f'{pdb_id}_1.pdb'))
    os.remove(os.path.join(tmpdir, f'Average_{pdb_id}.fxout'))        
    os.remove(os.path.join(tmpdir, f'Dif_{pdb_id}.fxout'))        
    os.remove(mut_path)        
    os.remove(os.path.join(tmpdir, f'PdbList_{pdb_id}.fxout'))        
    os.remove(os.path.join(tmpdir, f'Raw_{pdb_id}.fxout'))        
    os.remove(os.path.join(tmpdir, f'WT_{pdb_id}_1.pdb'))        

    return move_path

def clean_folder(tmpdir, ext=".pdb"):
    # First remove all .pdb files in the directory
    for filename in os.listdir(tmpdir):
        if filename.endswith(ext):
            file_path = os.path.join(tmpdir, filename)
            try:
                os.remove(file_path)
                #print(f"Removed file: {file_path}")
            except Exception as e:
                print(f"Error removing {file_path}: {e}")

    # Then remove the directory itself
    shutil.rmtree(tmpdir, ignore_errors=True)
    print(f"Removed directory: {tmpdir}")

def build_pdb(args):

    tmpdir = f'./tmp_build_{args.name}'
    os.makedirs(tmpdir, exist_ok=True)

    tmpdir2 = f'./tmp_build2_{args.name}'
    os.makedirs(tmpdir2, exist_ok=True)

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
    
    # setup outputs 
    build_pdb_dir = f'./metrics/mutant_structure/{name}'
    os.makedirs(tmpdir, exist_ok=True)

    df = pd.read_csv(excel_file)     
    if args.start_idx is not None:
        df = df[args.start_idx:]

    total_run = len(df)
    current_run = 0
    for idx, row in df.iterrows():
        current_run += 1
        progress = (current_run / total_run) * 100
        print(f"{progress:.2f}% ({current_run}/{total_run})")

        structure = load_structure(pdb_file)
        _, native_seqs = extract_coords_from_complex(structure)

        mutated_seqs = {}
        mutated_seqs[heavy_chain_id] = row['mut_heavy_chain_seq']
        mutated_seqs[light_chain_id] = native_seqs[light_chain_id]
        for c in antigen_chains:
            mutated_seqs[c] = native_seqs[c]

        mut_info = []
        for chain in chain_order:
            native_seq = native_seqs[chain]
            mut_seq = mutated_seqs[chain]
            if len(native_seq) != len(mut_seq):
                raise ValueError(f"Error, length is not matched!!! chain={chain}, wt={len(native_seq)}, mut={len(mut_seq)}")
            mutations = [(i+1, native_seq[i], mut_seq[i]) for i in range(len(native_seq)) if native_seq[i] != mut_seq[i]]

            for single_mutation in mutations:
                pos, wt, mt = single_mutation
                mut_info.append(f'{wt}{chain}{pos}{mt}')                
        mut_info = ','.join(mut_info)
        mut_info += ';'

        mutations = mut_info[:-1]

        pdb_name = pdb_name_offset + f'_{mutations}'
        mutated_path = os.path.join(build_pdb_dir, f'{pdb_name}.pdb')
        if os.path.exists(mutated_path) == False:
            sequential_pdb_file = os.path.join(tmpdir2, pdb_name + '.pdb')
            # sequential_pdb_file = os.path.join(tmpdir, pdb_name + '.pdb')
            _ = numbered_to_sequential(pdb_file, sequential_pdb_file)
            mutated_path = build_model(sequential_pdb_file, mut_info, pdb_name, outdir=build_pdb_dir, tmpdir=tmpdir)

        df.at[idx, 'mutated_pdb_path'] = mutated_path

    # remove tmp folders
    clean_folder(tmpdir, ext=".pdb")
    clean_folder(tmpdir2, ext=".pdb")
    return mutated_path

def parse():
    parser = ArgumentParser(description='Generate antibody')
    parser.add_argument('--tmp_num', type=int, default=20)
    parser.add_argument('--start_idx', type=int, default=0)

    # data
    parser.add_argument('--name', type=str, default='1mlc')

    # model
    parser.add_argument('--gpu', type=int, default=7)

    return parser.parse_args()


if __name__ == '__main__':
    args = parse()
    # mapping_dict = imgt_to_sequential(args)
    mutated_path = build_pdb(args)
    