import json
import copy
import torch
import argparse
import functools
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from torch.utils.data import DataLoader
from prepare_data import get_metadata, get_structure_data

from diffab.utils.misc import seed_all
from diffab.models import get_model
from diffab.datasets.custom import preprocess_antibody_structure
from diffab.tools.renumber import renumber as renumber_antibody
from diffab.utils.train import recursive_to
from diffab.utils.misc import *
from diffab.utils.data import *
from diffab.utils.transforms import *
from diffab.utils.inference import *

def create_data_variants(config, structure_factory):
    structure = structure_factory()
    structure_id = structure['id']

    data_variants = []
    if config.mode == 'single_cdr':
        cdrs = sorted(list(set(find_cdrs(structure)).intersection(config.sampling.cdrs)))
        for cdr_name in cdrs:
            transform = Compose([
                MaskSingleCDR(cdr_name, augmentation=False),
                MergeChains(),
            ])
            data_var = transform(structure_factory())
            residue_first, residue_last = get_residue_first_last(data_var)
            data_variants.append({
                'data': data_var,
                'name': f'{structure_id}-{cdr_name}',
                'tag': f'{cdr_name}',
                'cdr': cdr_name,
                'residue_first': residue_first,
                'residue_last': residue_last,
            })
    elif config.mode == 'multiple_cdrs':
        cdrs = sorted(list(set(find_cdrs(structure)).intersection(config.sampling.cdrs)))
        transform = Compose([
            MaskMultipleCDRs(selection=cdrs, augmentation=False),
            MergeChains(),
        ])
        data_var = transform(structure_factory())
        data_variants.append({
            'data': data_var,
            'name': f'{structure_id}-MultipleCDRs',
            'tag': 'MultipleCDRs',
            'cdrs': cdrs,
            'residue_first': None,
            'residue_last': None,
        })
    elif config.mode == 'full':
        transform = Compose([
            MaskAntibody(),
            MergeChains(),
        ])
        data_var = transform(structure_factory())
        data_variants.append({
            'data': data_var,
            'name': f'{structure_id}-Full',
            'tag': 'Full',
            'residue_first': None,
            'residue_last': None,
        })
    elif config.mode == 'abopt':
        cdrs = sorted(list(set(find_cdrs(structure)).intersection(config.sampling.cdrs)))
        for cdr_name in cdrs:
            transform = Compose([
                MaskSingleCDR(cdr_name, augmentation=False),
                MergeChains(),
            ])
            data_var = transform(structure_factory())
            residue_first, residue_last = get_residue_first_last(data_var)
            for opt_step in config.sampling.optimize_steps:
                data_variants.append({
                    'data': data_var,
                    'name': f'{structure_id}-{cdr_name}-O{opt_step}',
                    'tag': f'{cdr_name}-O{opt_step}',
                    'cdr': cdr_name,
                    'opt_step': opt_step,
                    'residue_first': residue_first,
                    'residue_last': residue_last,
                })
    else:
        raise ValueError(f'Unknown mode: {config.mode}.')
    return data_variants

@functools.lru_cache(maxsize=1)
def load_model(checkpoint_path, device):
    ckpt = torch.load(checkpoint_path, map_location=device)
    model = get_model(ckpt['config'].model).to(device)
    model.load_state_dict(ckpt['model'])
    return model

def score_pdb(pdb, heavy_chain, light_chain, loaded_model, config_yml, device, no_renumber=False):
    # Load configs
    config, config_name = load_config(config_yml)
    
    # Structure loading
    data_id = os.path.basename(pdb)
    if no_renumber:
        pdb_path = pdb
    else:
        in_pdb_path = pdb
        out_pdb_path = os.path.splitext(in_pdb_path)[0] + '_chothia.pdb'
        heavy_chains, light_chains = renumber_antibody(in_pdb_path, out_pdb_path)
        pdb_path = out_pdb_path

        if heavy_chain is None and len(heavy_chains) > 0:
            heavy_chain = heavy_chains[0]
        if light_chain is None and len(light_chains) > 0:
            light_chain = light_chains[0]
    if heavy_chain is None and light_chain is None:
        raise ValueError("Neither heavy chain id (--heavy) or light chain id (--light) is specified.")

    get_structure = lambda: preprocess_antibody_structure({
        'id': data_id,
        'pdb_path': pdb_path,
        'heavy_id': heavy_chain,
        # If the input is a nanobody, the light chain will be ignores
        'light_id': light_chain,
    })

    # Make data variants
    data_variants = create_data_variants(
        config = config,
        structure_factory = get_structure,
    )

    # Start modeling
    collate_fn = PaddingCollate(eight=False)
    inference_tfm = [ PatchAroundAnchor(), ]
    if 'abopt' not in config.mode:  # Don't remove native CDR in optimization mode
        inference_tfm.append(RemoveNative(
            remove_structure = config.sampling.sample_structure,
            remove_sequence = config.sampling.sample_sequence,
        ))
    inference_tfm = Compose(inference_tfm)

    model = loaded_model
    model.eval()
    batch_size = 1
    with torch.no_grad():
        for variant in data_variants:
            data_cropped = inference_tfm(
                copy.deepcopy(variant['data'])
            )
            data_list_repeat = [ data_cropped ] * config.sampling.num_samples
            loader = DataLoader(data_list_repeat, batch_size=batch_size, shuffle=False, collate_fn=collate_fn,
                                num_workers=4, pin_memory=True)
            
            count = 0
            for batch in tqdm(loader, desc=variant['name'], dynamic_ncols=True):
                batch = recursive_to(batch, device)
                # get sequence log-likelihood
                loss_dict = model.forward(batch)  
                # Convert the tensors to Python scalars and create a DataFrame
                loss_dict = {key: value.item() for key, value in loss_dict.items()}  # Convert tensors to scalars

    return loss_dict['loglikelihood']

def main(args):
    # load model 
    device = torch.device('cpu' if args.gpu == -1 else f'cuda:{args.gpu}')
    opt_model = load_model(args.opt_ckpt, device)
    fixbb_model = load_model(args.fixbb_ckpt, device)

    # load affinity data
    if args.name == "aayl49_ml":
        name = "aayl49_ML"
    else:
        name = args.name
    pdb_info = get_metadata()[name]
    df = pd.read_csv(pdb_info["affinity_data"][0])
    fout_path = f'./notebooks/scoring_outputs/{name}_benchmarking_data_diffab_scores.csv'
    opt_yml = f"./models/diffab/configs/opt.yml"
    fixbb_yml = f"./models/diffab/configs/fixbb.yml"

    # calculate sequence negative log-likelihood for mutant sequence
    fll_list, ll_list = [], []
    for idx, row in df.iterrows():
        pdb_data = f"{name}_{idx}.pdb"
        get_structure_data(pdb_info["pdb_path"], pdb_info["heavy_chain"], row["mut_heavy_chain_seq"], pdb_data, renumber=True)
        ll = score_pdb(pdb_data, pdb_info["heavy_chain"], pdb_info["light_chain"], opt_model, opt_yml, device)
        fll = score_pdb(pdb_data, pdb_info["heavy_chain"], pdb_info["light_chain"], fixbb_model, fixbb_yml, device)
        fll_list.append(fll)
        ll_list.append(ll)

        # Clear GPU memory
        torch.cuda.empty_cache()  # Clear unused memory from the cache

        # remove tmp file    
        Path(pdb_data).unlink()            
        Path(f'{name}_{idx}_chothia.pdb').unlink()

    # get scores
    df.loc[:, 'log-likelihood (fixed backbone)'] = fll_list
    df.loc[:, 'log-likelihood'] = ll_list
    df.to_csv(fout_path, header=True, index=False)
    return df 

def parse():
    parser = argparse.ArgumentParser(description='Calculate sequence log-likelihood from diffab')
    parser.add_argument('--name', type=str, default='3gbn', help='dataset name, default: 3gbn')
    parser.add_argument('--gpu', type=int, default=4, help='-1 for cpu')
    parser.add_argument('--seed', type=int, default=42, help='seed number')
    parser.add_argument('--opt_ckpt', type=str, default='./models/diffab/trained_models/codesign_single.pt')
    parser.add_argument('--fixbb_ckpt', type=str, default='./models/diffab/trained_models/fixbb.pt')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse()
    seed_all(args.seed)
    main(args)