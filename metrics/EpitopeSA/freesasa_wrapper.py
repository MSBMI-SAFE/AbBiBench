from Bio.PDB import PDBParser
import freesasa
import pandas as pd

class SASA:
    def __init__(self, pdb_path):
        # load structure
        structure = freesasa.Structure(pdb_path)
        # calculate residue-level solvent accessible surface area
        result = freesasa.calc(structure)
        self.sasa = result.residueAreas()

    def get_sasa(self, chain_id):
        chain_sasa_dict = self._collect_sasa_from_all_chains()
        if chain_id not in chain_sasa_dict.keys():
            raise ValueError(f'Error, {chain_id} not found!!!')
        else:
            return chain_sasa_dict[chain_id]

    def _collect_sasa_from_all_chains(self):
        chain_list = self.sasa.keys()
        df_dict = {} # chain_id: dataframe
        for chain in chain_list:
            data_list = self.sasa[chain].items()
            record_list = []
            for data in data_list:
                residueNumber = data[0]
                residueType = data[1].residueType
                totalSASA = data[1].total
                sidechainSASA = data[1].sideChain
                relativeMainChain = data[1].relativeMainChain
                relativeSideChain = data[1].relativeSideChain
                # append record
                record = (chain, residueNumber, residueType, totalSASA, sidechainSASA, relativeMainChain, relativeSideChain)
                record_list.append( record )
            # record to dataframe
            col_list = ['chain_id', 'resd_number', 'resd_type', 'total_sasa',
                        'sideChain_sasa', 'relative_mainchain', 'relative_sidechain']
            df = pd.DataFrame.from_records(record_list, columns=col_list)
            df_dict[chain] = df   
        return df_dict
