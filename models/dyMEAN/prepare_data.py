import json
from Bio import PDB
from anarci import run_anarci
from Bio.PDB import PDBParser, PDBIO
from Bio.SeqUtils import seq1


def get_metadata(json_path="./data/metadata.json"):
    with open(json_path, "r") as fin:
        data_dict = json.load(fin)
    return data_dict
    
def get_structure_data(pdb_path, chain_id, new_sequence, output_path, renumber=False):
    """Replace residues in a PDB chain with a new sequence."""
    
    # Load structure
    parser = PDB.PDBParser()
    structure = parser.get_structure('pdb', pdb_path)
    
    # Get target chain
    chain = structure[0][chain_id]
    
    if len(chain) != len(new_sequence):
        raise ValueError(f"New sequence length ({len(new_sequence)}) doesn't match chain length ({len(chain)})")
    
    # Create residue mapping
    res_mapping = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
        'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
        'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
        'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
    }
    
    # Replace residues
    for residue, new_aa in zip(chain, new_sequence):
        residue.resname = res_mapping[new_aa]
    
    # Save modified structure
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path)

    if renumber:
        # Renumber PDB using imgt scheme
        renumber_pdb(output_path, output_path, scheme='imgt', mute=True)


def renumber_seq(seq, scheme='imgt'): # copied from dyMEAN.renumbering.py
    _, numbering, details, _ = run_anarci([('A', seq)], scheme=scheme, allowed_species=['mouse', 'human'])
    numbering = numbering[0]
    fv, position = [], []
    if not numbering:  # not antibody
        return None
    chain_type = details[0][0]['chain_type']
    numbering = numbering[0][0]
    for pos, res in numbering:
        if res == '-':
            continue
        fv.append(res)
        position.append(pos)
    return ''.join(fv), position, chain_type

def renumber_pdb(pdb, out_pdb, scheme='imgt', mute=False): # copied from dyMEAN.renumbering.py
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('anonym', pdb)
    for chain in structure.get_chains():
        seq = []
        for residue in chain:
            hetero_flag, _, _ = residue.get_id()
            if hetero_flag != ' ':
                continue
            seq.append(seq1(residue.get_resname()))
        seq = ''.join(seq)
        res = renumber_seq(seq, scheme)
        if res is None:
            continue
        fv, position, chain_type = res
        if not mute:
            print(f'chain {chain.id} type: {chain_type}')
        start = seq.index(fv)
        end = start + len(fv)
        assert start != -1, 'fv not found'
        seq_index, pos_index = -1, 0
        for r in list(chain.get_residues()):
            hetero_flag, _, _ = r.get_id()
            if hetero_flag != ' ':
                continue
            seq_index += 1
            if seq_index < start or seq_index >= end:
                chain.__delitem__(r.get_id())
                continue
            assert fv[pos_index] == seq1(r.get_resname()), f'Inconsistent residue in Fv {fv[pos_index]} at {r._id}'
            r._id = (' ', *position[pos_index])
            pos_index += 1
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_pdb)