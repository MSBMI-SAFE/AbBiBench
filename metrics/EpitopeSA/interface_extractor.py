# Add the dyMEAN package to the Python path
import os, sys


#import dyMEAN
from pdb_utils import Protein, VOCAB

# other imports
import json
import argparse
import numpy as np
import pandas as pd

import project_utils as puts


class InterfaceExtractor:
    """
    modify dyMEAN's binding_interface.py to get outputs in dataframe/json format
    """
    def __init__(self, pdb, receptor_chains, ligand_chains, num_epitope_residues=48):
        self.pdb = pdb
        self.receptor_chains = receptor_chains
        self.ligand_chains = ligand_chains
        self.num_epitope_residues = num_epitope_residues

    def extract_binding_epitope(self):
        epitope, dists = self.get_interface(self.pdb, self.receptor_chains, self.ligand_chains, self.num_epitope_residues)
        record_list = []
        for i, (e, p) in enumerate(zip(epitope, dists)):
            e_res, e_chain_name, _ = e
            p_res, p_chain_name, _, d = p
            record_list.append( (e_chain_name, e_res.get_id()[0], e_res.get_symbol(), p_chain_name, p_res.get_id()[0], p_res.get_symbol(), round(d, 3)) )
        col_list =['epitope_chain', 'epitope_position', 'epitope_residue', 'paratope_chain', 'paratope_position', 'paratope_residue', 'distance']
        df = pd.DataFrame.from_records(record_list, columns=col_list)
        return df

    def get_interface(self, pdb, receptor_chains, ligand_chains, num_epitope_residues):
        prot = Protein.from_pdb(self.pdb)
        #print(f'pdb={prot.get_id()}: {prot.peptides} | receptor_chains={receptor_chains} | ligand_chains={ligand_chains}')
        for c in self.receptor_chains:
            assert c in prot.peptides, f'Chain {c} not found for receptor'
        for c in self.ligand_chains:
            assert c in prot.peptides, f'Chain {c} not found for ligand'
        receptor = Protein(prot.get_id(), {c: prot.get_chain(c) for c in receptor_chains})
        ligand = Protein(prot.get_id(), {c: prot.get_chain(c) for c in ligand_chains})
        #print(receptor)
        rec_rids, rec_xs, lig_rids, lig_xs = [], [], [], []
        rec_mask, lig_mask = [], []
        for _type, protein in zip(['rec', 'lig'], [receptor, ligand]):
            is_rec = _type == 'rec'
            rids = []
            if is_rec:
                rids, xs, masks = rec_rids, rec_xs, rec_mask
            else:
                rids, xs, masks = lig_rids, lig_xs, lig_mask
            for chain_name, chain in protein:
                for i, residue in enumerate(chain):
                    bb_coord = residue.get_backbone_coord_map()
                    sc_coord = residue.get_sidechain_coord_map()
                    coord = {}
                    coord.update(bb_coord)
                    coord.update(sc_coord)
                    num_pad = VOCAB.MAX_ATOM_NUMBER - len(coord)
                    x = [coord[key] for key in coord] + [[0, 0, 0] for _ in range(num_pad)]
                    mask = [1 for _ in coord] + [0 for _ in range(num_pad)]
                    rids.append((chain_name, i))
                    xs.append(x)
                    masks.append(mask)
        # calculate distance
        rec_xs, lig_xs = np.array(rec_xs), np.array(lig_xs) # [Nrec/lig, M, 3], M == MAX_ATOM_NUM
        rec_mask, lig_mask = np.array(rec_mask).astype('bool'), np.array(lig_mask).astype('bool')  # [Nrec/lig, M]
        dist = np.linalg.norm(rec_xs[:, None] - lig_xs[None, :], axis=-1)  # [Nrec, Nlig, M]
        dist = dist + np.logical_not(rec_mask[:, None] * lig_mask[None, :]) * 1e6  # [Nrec, Nlig, M]
        dist_mat = np.min(dist, axis=-1)  # [Nrec, Nlig]
        min_dists = np.min(dist_mat, axis=-1)  # [rec_len]
        topk = min(len(min_dists), self.num_epitope_residues)
        ind = np.argpartition(-min_dists, -topk)[-topk:]
        lig_idxs = np.argmin(dist_mat, axis=-1)  # [Nrec]
        epitope, dists = [], []
        for idx in ind:
            # epitope
            chain_name, i = rec_rids[idx]
            residue = receptor.peptides[chain_name].get_residue(i)
            epitope.append((residue, chain_name, i))
            # nearest ligand residue
            chain_name, i = lig_rids[lig_idxs[idx]]
            residue = ligand.peptides[chain_name].get_residue(i)
            dists.append((residue, chain_name, i ,min_dists[idx]))
        return epitope, dists