import json
import numpy as np
from typing import List
import biotite
from biotite.structure.io import pdbx, pdb
from biotite.structure import get_chains
from biotite.structure import filter_backbone
from biotite.structure.residues import get_residues
from biotite.sequence import ProteinSequence

def get_metadata(json_path="./data/metadata.json"):
    with open(json_path, "r") as fin:
        data_dict = json.load(fin)
    return data_dict
    
def load_structure(fpath, chain=None, include_O=False):
    """
    Args:
        fpath: filepath to either pdb or cif file
        chain: the chain id or list of chain ids to load
    Returns:
        biotite.structure.AtomArray
    """
    if fpath.endswith('cif'):
        with open(fpath) as fin:
            pdbxf = pdbx.PDBxFile.read(fin)
        structure = pdbx.get_structure(pdbxf, model=1)
    elif fpath.endswith('pdb'):
        with open(fpath) as fin:
            pdbf = pdb.PDBFile.read(fin)
        structure = pdb.get_structure(pdbf, model=1)

    if include_O:
        bbmask = np.isin(structure.atom_name, ["N", "CA", "C", "O"])
    else:
        bbmask = filter_backbone(structure)
    # bbmask = filter_backbone(structure)
    structure = structure[bbmask]
    all_chains = get_chains(structure)
    if len(all_chains) == 0:
        raise ValueError('No chains found in the input file.')
    if chain is None:
        chain_ids = all_chains
    elif isinstance(chain, list):
        chain_ids = chain
    else:
        chain_ids = [chain] 
    for chain in chain_ids:
        if chain not in all_chains:
            raise ValueError(f'Chain {chain} not found in input file')
    chain_filter = [a.chain_id in chain_ids for a in structure]
    structure = structure[chain_filter]
    return structure


def extract_coords_from_complex(structure: biotite.structure.AtomArray, include_O=False):
    """
    Args:
        structure: biotite AtomArray
    Returns:
        Tuple (coords_list, seq_list)
        - coords: Dictionary mapping chain ids to L x 3 x 3 array for N, CA, C
          coordinates representing the backbone of each chain
        - seqs: Dictionary mapping chain ids to native sequences of each chain
    """
    coords = {}
    seqs = {}
    all_chains = get_chains(structure)
    for chain_id in all_chains:
        chain = structure[structure.chain_id == chain_id]
        coords[chain_id], seqs[chain_id] = extract_coords_from_structure(chain, include_O)
        # coords[chain_id], seqs[chain_id] = extract_coords_from_structure(chain)
    return coords, seqs

def extract_coords_from_structure(structure: biotite.structure.AtomArray, include_O=False):
    """
    Args:
        structure: An instance of biotite AtomArray
    Returns:
        Tuple (coords, seq)
            - coords is an L x 3 x 3 array for N, CA, C coordinates
            - seq is the extracted sequence
    """
    # Filter out HETATM by checking the 'atom_id' field and standard amino acids - CGU같은게 나와서 에러 나서 내가 추가함
    atom_is_standard = structure.hetero == False  # Exclude HETATM
    structure = structure[atom_is_standard]

    if include_O:
        coords = get_atom_coords_residuewise(["N", "CA", "C", "O"], structure)
    else:
        coords = get_atom_coords_residuewise(["N", "CA", "C"], structure)    
    # coords = get_atom_coords_residuewise(["N", "CA", "C"], structure)
    residue_identities = get_residues(structure)[1]
    seq = ''.join([ProteinSequence.convert_letter_3to1(r) for r in residue_identities])
    return coords, seq

def get_atom_coords_residuewise(atoms: List[str], struct: biotite.structure.AtomArray):
    """
    Example for atoms argument: ["N", "CA", "C"]
    """
    def filterfn(s, axis=None):
        filters = np.stack([s.atom_name == name for name in atoms], axis=1)
        sum = filters.sum(0)
        if not np.all(sum <= np.ones(filters.shape[1])):
            raise RuntimeError("structure has multiple atoms with same name")
        index = filters.argmax(0)
        coords = s[index].coord
        coords[sum == 0] = float("nan")
        return coords

    return biotite.structure.apply_residue_wise(struct, struct, filterfn)