#!/usr/bin/python
# -*- coding:utf-8 -*-
from anarci import run_anarci
from Bio.PDB import PDBParser, PDBIO
from Bio.SeqUtils import seq1
import os
import json

# 1. basic variables
PROJ_DIR = os.path.split(__file__)[0]
RENUMBER = os.path.join(PROJ_DIR, 'utils', 'renumber.py')

# 3. antibody numbering, [start, end] of residue id, both start & end are included
# 3.1 IMGT numbering definition
class IMGT:
    # heavy chain
    HFR1 = (1, 26)
    HFR2 = (39, 55)
    HFR3 = (66, 104)
    HFR4 = (118, 129)

    H1 = (27, 38)
    H2 = (56, 65)
    H3 = (105, 117)

    # light chain
    LFR1 = (1, 26)
    LFR2 = (39, 55)
    LFR3 = (66, 104)
    LFR4 = (118, 129)

    L1 = (27, 38)
    L2 = (56, 65)
    L3 = (105, 117)

    Hconserve = {
        23: ['CYS'],
        41: ['TRP'],
        104: ['CYS']
    }

    Lconserve = {
        23: ['CYS'],
        41: ['TRP'],
        104: ['CYS']
    }

    @classmethod
    def renumber(cls, pdb, out_pdb):
        code = os.system(f'python {RENUMBER} {pdb} {out_pdb} imgt 0')
        return code

def extract_chain_sequence(chain, verbose=False):
    """Extract amino acid sequence from a PDB chain with detailed handling"""
    seq = []
    skipped_residues = []
    
    for residue in chain:
        hetero_flag, resseq, icode = residue.get_id()
        resname = residue.get_resname()
        
        # Skip water molecules
        if resname == 'HOH':
            continue
            
        try:
            aa = seq1(resname)
            # Keep residue if it's a standard amino acid, regardless of hetero flag
            if aa != "X":
                if verbose and hetero_flag.strip():
                    print(f"Including HETATM residue {resname} ({aa}) at position {resseq}")
                seq.append(aa)
            else:
                if verbose:
                    print(f"Skipping non-standard residue: {resname} at position {resseq}")
                skipped_residues.append((resname, resseq))
        except Exception as e:
            if verbose:
                print(f"Error converting residue {resname} at position {resseq}: {str(e)}")
            skipped_residues.append((resname, resseq))
            
    if verbose and skipped_residues:
        print(f"Skipped residues (resname, position): {skipped_residues}")
        
    return ''.join(seq)

def get_position_num(pos):
    """Extract numerical position from position tuple/number"""
    if isinstance(pos, tuple):
        return pos[0]  # Get the numerical part
    return pos

def get_cdr_info(numbering, seq, cdr_range):
    """Extract CDR position and sequence information
    Returns:
        tuple: (residue positions, index positions, sequence)
        - residue positions: PDB residue numbers
        - index positions: 1-indexed position in sequence
        - sequence: CDR sequence string
    """
    start, end = cdr_range
    cdr_positions = []
    cdr_seq = []
    idx_positions = []  # Track 1-indexed positions
    
    for i, (pos, aa) in enumerate(zip(numbering, seq), 1):  # 1-indexed
        pos_num = get_position_num(pos)
        if start <= pos_num <= end:
            cdr_positions.append(pos)
            cdr_seq.append(aa)
            idx_positions.append(i)
    
    if not cdr_positions:
        return [], [], ""
            
    return ([get_position_num(min(cdr_positions)), get_position_num(max(cdr_positions))], 
            [min(idx_positions), max(idx_positions)],
            ''.join(cdr_seq))

def renumber_seq(seq, scheme='imgt'):
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

def create_summary_json(pdb_file, structure, heavy_chain, light_chain, antigen_chains, scheme='imgt', pre_numbered=False):
    """Create summary JSON with antibody information"""
    summary = {
        "pdb": os.path.splitext(os.path.basename(pdb_file))[0],
        "heavy_chain": heavy_chain,
        "light_chain": light_chain,
        "antigen_chains": antigen_chains,
        "pdb_data_path": pdb_file,
        "numbering": scheme,
        "pre_numbered": pre_numbered
    }
    
    # Extract sequences
    h_chain = structure[0][heavy_chain]
    l_chain = structure[0][light_chain]
    
    # Extract sequences with verbose output for antigen chains
    h_seq = extract_chain_sequence(h_chain)
    l_seq = extract_chain_sequence(l_chain)
    antigen_seqs = []
    for chain_id in antigen_chains:
        if chain_id in structure[0]:
            chain = structure[0][chain_id]
            #print(f"\nExtracting sequence for antigen chain {chain_id}:")
            seq = extract_chain_sequence(chain, verbose=True)
            antigen_seqs.append(seq)
            #print(f"Chain {chain_id} sequence length: {len(seq)}")
        else:
            print(f"Warning: Chain {chain_id} not found in structure")
    
    # Get numbering information
    h_res = renumber_seq(h_seq, scheme)
    l_res = renumber_seq(l_seq, scheme)
    
    if h_res and l_res:
        h_fv, h_pos, _ = h_res
        l_fv, l_pos, _ = l_res
        
        summary.update({
            "heavy_chain_seq": h_seq,
            "light_chain_seq": l_seq,
            "antigen_seqs": antigen_seqs
        })
        
        # Extract CDR information
        if scheme.lower() == 'imgt':
            numbering_class = IMGT
        else:
            numbering_class = 'chothia'
            
        # Heavy chain CDRs
        for i in range(1, 4):
            cdr_range = getattr(numbering_class, f'H{i}')
            resd_pos, idx_pos, seq = get_cdr_info(h_pos, h_fv, cdr_range)
            summary[f'cdrh{i}_resd_pos'] = resd_pos
            summary[f'cdrh{i}_idx_pos'] = idx_pos
            summary[f'cdrh{i}_seq'] = seq
            
        # Light chain CDRs
        for i in range(1, 4):
            cdr_range = getattr(numbering_class, f'L{i}')
            resd_pos, idx_pos, seq = get_cdr_info(l_pos, l_fv, cdr_range)
            summary[f'cdrl{i}_resd_pos'] = resd_pos
            summary[f'cdrl{i}_idx_pos'] = idx_pos
            summary[f'cdrl{i}_seq'] = seq
    
    return summary

def renumber_pdb(pdb, out_pdb, scheme='imgt', summary_json=None, mute=False, heavy_chain='H', light_chain='L', antigen_chains=None):
    """Renumber PDB and optionally create summary JSON"""
    if antigen_chains is None:
        antigen_chains = []
        
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('anonym', pdb)
    
    # First pass: identify chains and renumber
    for chain in structure.get_chains():
        seq = extract_chain_sequence(chain)
        res = renumber_seq(seq, scheme)
        if res is None:
            continue
        fv, position, chain_type = res
        if not mute:
            print(f'chain {chain.id} type: {chain_type}')
        
        # Renumber the chain
        start = seq.index(fv)
        end = start + len(fv)
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
    
    # Save renumbered PDB
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_pdb)
    
    # Create and save summary JSON if requested
    if summary_json:
        summary = create_summary_json(out_pdb, structure, heavy_chain, light_chain, antigen_chains, scheme)
        with open(summary_json, 'w') as f:
            json.dump(summary, f, indent=2)

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 4:
        print("Usage: python renumber.py input.pdb output.pdb scheme [summary.json] [mute] [heavy_chain] [light_chain] [antigen_chains]")
        sys.exit(1)
        
    infile, outfile, scheme = sys.argv[1:4]
    summary_json = sys.argv[4] if len(sys.argv) > 4 else None
    mute = bool(sys.argv[5]) if len(sys.argv) > 5 else False
    heavy_chain = sys.argv[6] if len(sys.argv) > 6 else 'H'
    light_chain = sys.argv[7] if len(sys.argv) > 7 else 'L'
    def parse_chain_input(chain_str):
        # Handle empty input
        if not chain_str:
            return []
        
        # First try comma-separated format
        if ',' in chain_str:
            chains = [c.strip() for c in chain_str.split(',')]
        else:
            # If no commas, split into individual characters
            chains = list(chain_str.replace(' ', ''))
        
        # Remove any empty strings and return unique chains
        return list(filter(None, chains))

    antigen_chains = parse_chain_input(sys.argv[8]) if len(sys.argv) > 8 else []
    
    renumber_pdb(infile, outfile, scheme, summary_json, mute, heavy_chain, light_chain, antigen_chains)