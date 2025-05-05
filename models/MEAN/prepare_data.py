import os
import json
import argparse

def get_antibody_metadata(pdb_id):
    if pdb_id == '3gbn':
        pdb = '3gbn_hlab'
        pdb_path = f"./data/complex_structure/{pdb_id}_hlab.pdb"
        json_path = f'./models/MEAN/data/{pdb_id}/pdb_info.json'
        heavy_chain_id = 'H' ; light_chain_id = 'L' ; antigen_chains = ['A','B']

    elif pdb_id == '4fqi':
        pdb = '4fqi_hlab'
        pdb_path = f"./data/complex_structure/{pdb_id}_hlab.pdb"
        json_path = f'./models/MEAN/data/{pdb_id}/pdb_info.json'
        heavy_chain_id = 'H' ; light_chain_id = 'L' ; antigen_chains = ['A','B']

    elif pdb_id == '2fjg':
        pdb = '2fjg_hlv'
        pdb_path = f"./data/complex_structure/{pdb_id}_hlv.pdb"
        json_path = f'./models/MEAN/data/{pdb_id}/pdb_info.json'
        heavy_chain_id = 'H' ; light_chain_id = 'L' ; antigen_chains = ['V']

    elif pdb_id == 'aayl49':
        pdb = 'AAYL49_bca'
        pdb_path = f"./data/complex_structure/AAYL49_bca.pdb"
        json_path = f'./models/MEAN/data/{pdb_id}/pdb_info.json'
        heavy_chain_id = 'B' ; light_chain_id = 'C' ; antigen_chains = ['A']

    elif pdb_id == 'aayl51':
        pdb = 'AAYL51_bca'
        pdb_path = f"./data/complex_structure/AAYL51_bca.pdb"
        json_path = f'./models/MEAN/data/{pdb_id}/pdb_info.json'
        heavy_chain_id = 'B' ; light_chain_id = 'C' ; antigen_chains = ['A']

    elif pdb_id == '1mlc':
        pdb = '1mlc_bae'
        pdb_path = f"./data/complex_structure/1mlc_bae.pdb"
        json_path = f'./models/MEAN/data/{pdb_id}/pdb_info.json'
        heavy_chain_id = 'B' ; light_chain_id = 'A' ; antigen_chains = ['E']

    elif pdb_id == '1n8z':
        pdb = f'{pdb_id}_bac'
        pdb_path = f"./data/complex_structure/{pdb_id}_bac.pdb" 
        json_path = f'./models/MEAN/data/{pdb_id}/pdb_info.json'
        heavy_chain_id = 'B' ; light_chain_id = 'A' ; antigen_chains = ['C']

    return pdb, pdb_path, json_path, heavy_chain_id, light_chain_id, antigen_chains

def main(args):    
    pdb, pdb_path, json_path, heavy_chain_id, light_chain_id, antigen_chains = get_antibody_metadata(args.pdb_id)
    
    data = {"pdb": pdb, "pdb_data_path": pdb_path, "heavy_chain": heavy_chain_id, "light_chain": light_chain_id, "antigen_chains": antigen_chains}
    os.makedirs(os.path.dirname(json_path), exist_ok=True)
    with open(json_path, "w") as json_file:
        json.dump(data, json_file) #, indent=4)

def parse():
    parser = argparse.ArgumentParser(description='generation by MEAN')
    parser.add_argument('--pdb_id', type=str, default='aayl49', choices=['3gbn','4fqi','2fjg','aayl49','aayl51', '1mlc', '1n8z'])
    return parser.parse_args()

if __name__ == '__main__':
    args = parse()
    main(args)