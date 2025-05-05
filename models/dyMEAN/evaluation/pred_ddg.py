#!/usr/bin/python
# -*- coding:utf-8 -*-
import os
import shutil

import torch

from utils.time_sign import get_time_sign

FILE_DIR = os.path.abspath(os.path.split(__file__)[0])
MODULE_DIR = os.path.join(FILE_DIR, 'ddg')
from .ddg.models.predictor import DDGPredictor
from .ddg.utils.misc import *
from .ddg.utils.data import *
from .ddg.utils.protein import *

CKPT = torch.load(os.path.join(MODULE_DIR, 'data', 'model.pt'))
MODEL = DDGPredictor(CKPT['config'].model)
MODEL.load_state_dict(CKPT['model'])
DEVICE = torch.device('cuda:0')
MODEL.to(DEVICE)
MODEL.eval()


def pred_ddg(mut_pdb, wt_pdb):
    batch = load_wt_mut_pdb_pair(wt_pdb, mut_pdb)
    batch = recursive_to(batch, DEVICE)
    wt_size = len(batch['wt']['aa'][0])
    mut_size = len(batch['mut']['aa'][0])
    if wt_size != mut_size:
        print(f'WARRNING, dimension not equal: wt={wt_size}, mut={mut_size}, setting ddg=0')
        ddg = 0
    else:
        #print(f'loading model to predict ddg.....')    
        with torch.no_grad():
            pred = MODEL(batch['wt'], batch['mut']).item()

    return pred

def _pred_ddg(mut_pdb, wt_pdb, device):
    """
    return predicted delta delta g
    """
    model_ckpt = os.path.join(MODULE_DIR, 'data', 'model.pt')
    batch = load_wt_mut_pdb_pair(wt_pdb, mut_pdb)
    batch = recursive_to(batch, device)
    wt_size = len(batch['wt']['aa'][0])
    mut_size = len(batch['mut']['aa'][0])
    if wt_size != mut_size:
        print(f'WARRNING, dimension not equal: wt={wt_size}, mut={mut_size}, setting ddg=0')
        ddg = 0
    else:
        print(f'loading model to predict ddg.....')

        ckpt = torch.load(model_ckpt, map_location=device)
        config = ckpt['config']
        weight = ckpt['model']
        model = DDGPredictor(config.model).to(device)
        model.load_state_dict(weight)

        with torch.no_grad():
            model.eval()
            pred = model(batch['wt'], batch['mut'])
            ddg = pred.item()
            #print(ddg)
    return ddg




from configs import FOLDX_BIN, CACHE_DIR


def foldx_minimize_energy(pdb_path, out_path=None):
    filename = get_time_sign() + os.path.basename(os.path.splitext(pdb_path)[0]) + '_foldx.pdb'
    tmpfile = os.path.join(CACHE_DIR, filename)
    shutil.copyfile(pdb_path, tmpfile)
    p = os.popen(f'cd {CACHE_DIR}; {FOLDX_BIN} --command=Optimize --pdb={filename}')
    p.read()
    p.close()
    os.remove(tmpfile)
    filename = 'Optimized_' + filename
    tmpfile = os.path.join(CACHE_DIR, filename)
    if out_path is None:
        out_path = os.path.join(os.path.split(pdb_path)[0], filename)
    shutil.copyfile(tmpfile, out_path)
    os.remove(tmpfile)
    return out_path


def foldx_dg(pdb_path, rec_chains, lig_chains):
    filename = get_time_sign() + os.path.basename(os.path.splitext(pdb_path)[0]) + '_foldx.pdb'
    tmpfile = os.path.join(CACHE_DIR, filename)
    shutil.copyfile(pdb_path, tmpfile)
    rec, lig = ''.join(rec_chains), ''.join(lig_chains)
    # p = os.popen(f'cd {CACHE_DIR}; {FOLDX_BIN} --command=Optimize --pdb={filename}')
    # p.read()
    # p.close()
    # os.remove(tmpfile)

    # filename = 'Optimized_' + filename
    # tmpfile = os.path.join(CACHE_DIR, filename)
    p = os.popen(f'cd {CACHE_DIR}; {FOLDX_BIN} --command=AnalyseComplex --pdb={filename} --analyseComplexChains={rec},{lig}')
    aff = float(p.read().split('\n')[-8].split(' ')[-1])
    p.close()
    os.remove(tmpfile)
    return aff


def foldx_ddg(wt_pdb, mut_pdb, rec_chains, lig_chains):
    wt_dg = foldx_dg(wt_pdb, rec_chains, lig_chains)
    mut_dg = foldx_dg(mut_pdb, rec_chains, lig_chains)
    return mut_dg - wt_dg
