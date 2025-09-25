import os
import json
import copy
import torch
import argparse
import functools
import pandas as pd
import numpy as np
import multiprocessing as mp
from tqdm import tqdm
from pathlib import Path
from datetime import datetime
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


def load_model(checkpoint_path, device):
    """Load model without caching since each process needs its own model"""
    ckpt = torch.load(checkpoint_path, map_location=device)
    model = get_model(ckpt['config'].model).to(device)
    model.load_state_dict(ckpt['model'])
    return model


def score_pdb(pdb, heavy_chain, light_chain, loaded_model, config_yml, device, no_renumber=False):
    config, config_name = load_config(config_yml)
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
        'light_id': light_chain,
    })

    data_variants = create_data_variants(config=config, structure_factory=get_structure)

    # FIX 1: Set num_workers=0 to avoid subprocess creation in daemon processes
    collate_fn = PaddingCollate(eight=False)
    inference_tfm = [PatchAroundAnchor()]
    if 'abopt' not in config.mode:
        inference_tfm.append(RemoveNative(
            remove_structure=config.sampling.sample_structure,
            remove_sequence=config.sampling.sample_sequence,
        ))
    inference_tfm = Compose(inference_tfm)

    model = loaded_model
    model.eval()
    with torch.no_grad():
        for variant in data_variants:
            data_cropped = inference_tfm(copy.deepcopy(variant['data']))
            data_list_repeat = [data_cropped] * config.sampling.num_samples
            
            # CRITICAL FIX: Set num_workers=0 to avoid daemon process spawning children
            loader = DataLoader(data_list_repeat, batch_size=1, shuffle=False,
                               collate_fn=collate_fn, num_workers=0, pin_memory=False)

            for batch in tqdm(loader, desc=variant['name'], dynamic_ncols=True):
                batch = recursive_to(batch, device)
                loss_dict = model.forward(batch)
                loss_dict = {key: value.item() for key, value in loss_dict.items()}

    return loss_dict['loglikelihood']


def score_df_chunk(chunk_df, pdb_info, opt_model, opt_yml, fixbb_model, fixbb_yml, device, name):
    """Score a dataframe chunk on a specific device"""
    fll_list, ll_list = [], []
    
    for idx, row in chunk_df.iterrows():
        try:
            pdb_data = f"data/mutant_structure/mutant_structure_{pdb_info['pdb']}/variant_{idx}.pdb"
            print(f"Processing {pdb_data} on {device}")
            
            ll = score_pdb(pdb_data, pdb_info["heavy_chain"], pdb_info["light_chain"], 
                          opt_model, opt_yml, device)
            fll = score_pdb(pdb_data, pdb_info["heavy_chain"], pdb_info["light_chain"], 
                           fixbb_model, fixbb_yml, device)
            
            fll_list.append(fll)
            ll_list.append(ll)
            
            # Clean up GPU memory and temporary files
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
            
            # Clean up renumbered PDB files more carefully
            chothia_file = Path(f'{pdb_data[:-4]}_chothia.pdb')
            if chothia_file.exists():
                try:
                    chothia_file.unlink()
                except OSError:
                    pass  # File might be in use, ignore cleanup error
            
        except Exception as e:
            print(f"Error processing variant {idx}: {str(e)}")
            fll_list.append(None)
            ll_list.append(None)

    # Create a copy of the chunk to avoid modifying the original
    result_df = chunk_df.copy()
    result_df.loc[:, 'log-likelihood (fixed backbone)'] = fll_list
    result_df.loc[:, 'log-likelihood'] = ll_list
    
    return result_df


def process_chunk_worker(chunk_info):
    """Worker function for multiprocessing"""
    try:
        chunk_idx, chunk_df, pdb_info, args, name, opt_yml, fixbb_yml = chunk_info
        
        # Set up device for this worker
        gpu_id = args.gpu[chunk_idx]
        device = torch.device('cpu' if gpu_id == -1 else f'cuda:{gpu_id}')
        
        print(f"Worker {chunk_idx}: Starting on device {device} with {len(chunk_df)} samples")
        start_time = datetime.now()
        
        # Set torch thread settings for better multiprocessing
        torch.set_num_threads(1)
        
        # Load models on this worker's device
        opt_model = load_model(args.opt_ckpt, device=device)
        fixbb_model = load_model(args.fixbb_ckpt, device=device)
        
        # Process the chunk
        result_df = score_df_chunk(chunk_df, pdb_info, opt_model, opt_yml, 
                                  fixbb_model, fixbb_yml, device, name)
        
        end_time = datetime.now()
        print(f"Worker {chunk_idx}: Completed in {end_time - start_time}")
        
        # Clean up models
        del opt_model, fixbb_model
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        
        return result_df
        
    except Exception as e:
        print(f"Error in worker {chunk_idx}: {str(e)}")
        import traceback
        traceback.print_exc()
        # Return empty dataframe with same structure on error
        error_df = chunk_df.copy()
        error_df.loc[:, 'log-likelihood (fixed backbone)'] = [None] * len(chunk_df)
        error_df.loc[:, 'log-likelihood'] = [None] * len(chunk_df)
        return error_df


def main(args):
    print(f"Starting parallel processing with {len(args.gpu)} GPUs: {args.gpu}")
    
    # Setup data
    name = "aayl49_ML" if args.name == "aayl49_ml" else args.name
    pdb_info = get_metadata()[name]
    pdb = pdb_info["pdb"]
    if args.affinity_data:
        affinity_data = args.affinity_data
    else:
        affinity_data = pdb_info["affinity_data"][0]
    df = pd.read_csv(affinity_data)
    
    print(f"Total samples to process: {len(df)}")
    
    # Setup output and config paths
    prefix = args.prefix
    fout_path = f'./notebooks/scoring_outputs/{prefix}_benchmarking_data_diffab_scores.csv'
    opt_yml = "./models/diffab/configs/opt.yml"
    fixbb_yml = "./models/diffab/configs/fixbb.yml"
    
    # Split dataframe into chunks
    chunks = np.array_split(df, len(args.gpu))
    print(f"Split into {len(chunks)} chunks with sizes: {[len(chunk) for chunk in chunks]}")
    
    # Prepare chunk information for workers
    chunk_infos = []
    for chunk_idx, chunk_df in enumerate(chunks):
        chunk_info = (chunk_idx, chunk_df, pdb_info, args, name, opt_yml, fixbb_yml)
        chunk_infos.append(chunk_info)
    
    # Process chunks in parallel
    start_time = datetime.now()
    print(f"Starting parallel processing at {start_time}")
    
    try:
        # FIX 2: Use 'forkserver' or 'spawn' method to avoid CUDA context issues
        # 'forkserver' is generally safer for CUDA applications than 'spawn'
        if hasattr(mp, 'get_start_method'):
            current_method = mp.get_start_method()
            print(f"Current multiprocessing method: {current_method}")
            
        try:
            mp.set_start_method('forkserver', force=True)
            print("Using 'forkserver' multiprocessing method")
        except:
            mp.set_start_method('spawn', force=True)
            print("Using 'spawn' multiprocessing method")
        
        # FIX 3: Use context manager properly and handle exceptions
        with mp.Pool(processes=len(args.gpu)) as pool:
            try:
                df_list = pool.map(process_chunk_worker, chunk_infos)
            except Exception as e:
                print(f"Error during parallel processing: {str(e)}")
                pool.terminate()
                pool.join()
                raise
        
        # Merge results
        merged_df = pd.concat(df_list, ignore_index=True)
        
        # Sort by original index to maintain order
        merged_df = merged_df.sort_index()
        
        # Save results
        os.makedirs(os.path.dirname(fout_path), exist_ok=True)
        merged_df.to_csv(fout_path, header=True, index=False)
        
        end_time = datetime.now()
        total_time = end_time - start_time
        print(f"Parallel processing completed in {total_time}")
        print(f"Results saved to: {fout_path}")
        print(f"Successfully processed: {len(merged_df)} samples")
        
        # FIX 4: Better error reporting - check for actual successful processing
        successful_samples = (~merged_df['log-likelihood'].isna()).sum()
        failed_samples = merged_df['log-likelihood'].isna().sum()
        
        print(f"Actually successful samples: {successful_samples}")
        if failed_samples > 0:
            print(f"Warning: {failed_samples} samples failed to process")
        
        return merged_df
        
    except Exception as e:
        print(f"Error in main processing: {str(e)}")
        import traceback
        traceback.print_exc()
        raise


def parse():
    parser = argparse.ArgumentParser(description='Calculate sequence log-likelihood from diffab using parallel processing')
    parser.add_argument('--name', type=str, default='3gbn', help='dataset name, default: 3gbn')
    parser.add_argument('--affinity_data', type=str, help='Path to binding_score csv file')
    parser.add_argument('--gpu', type=lambda s: [int(x) for x in s.split(',')], default=[0, 1], 
                       help='comma-separated list of GPU ids, e.g., "0,1,2,3"')
    parser.add_argument('--seed', type=int, default=42, help='seed number')
    parser.add_argument('--prefix', type=str, default='3gbn', help='prefix of the output scores.csv, default: 3gbn')
    parser.add_argument('--opt_ckpt', type=str, default='./models/diffab/trained_models/codesign_single.pt',
                       help='path to optimization checkpoint')
    parser.add_argument('--fixbb_ckpt', type=str, default='./models/diffab/trained_models/fixbb.pt',
                       help='path to fixed backbone checkpoint')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse()
    seed_all(args.seed)
    
    print("="*50)
    print("PARALLEL DIFFAB SCORING")
    print("="*50)
    print(f"GPUs: {args.gpu}")
    print(f"Dataset: {args.name}")
    print(f"Seed: {args.seed}")
    print(f"Output prefix: {args.prefix}")
    print("="*50)
    
    main(args)