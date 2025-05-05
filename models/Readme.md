Below, we show calculation of log-likelihood from computational models.
Here we use the **AAYL49.pdb** ab-ag complex and its mutant sequences **filtered_ML_data_no_duplicates.csv** as an example.
Please refer to data/binding_affinity/Readme.md for data information.

# ESM-2

```{bash}
$ python scoring/ESM-2/get_model_likelihood.py 

```

# ESM-IF

```{bash}
$ python models/ESM-IF/structural-evolution/bin/get_ESMIF1_likelihood.py \
--pdb_file /isilon/ytang4/FLU_Project/benchmark/data/complex_structure/4fqi_hlab.pdb \
--csv_file /isilon/ytang4/FLU_Project/benchmark/data/binding_affinity/4fqi_h1_benchmarking_data.csv

```

# MEAN

```{bash}
# activate environment
$ conda activate MEAN

# get model predicted likelihood
$ python scoring/MEAN/get_model_likelihood.py

```

# dyMEAN

```{bash}
# activate environment
$ conda activate dyMEAN

# get model predicted likelihhod
$ python scoring/dyMEAN/scripts/get_model_likelihood.py \
    --gpu 6 \
    --pdb_path data/complex_structure/AAYL49.pdb \
    --heavy_chain B \
    --light_chain C \
    --antigen_chain A \
    --dms_path data/binding_affinity/filtered_ML_data_no_duplicates.csv \
    --save_dir ./
```

# diffab

```{bash}
# activate environment
$ conda activate diffab

# get model predicted likelihood
$ python scoring/diffab/get_model_likelihood.py \
    --opt_ckpt scoring/diffab/trained_models/codesign_single.pt \
    --fixbb_ckpt scoring/diffab/trained_models/fixbb.pt
    --pdb_path data/complex_structure/AAYL49.pdb \
    --heavy_chain B \ 
    --light_chain C \
    --antigen_chain A \ 
    --dms_path data/binding_affinity/filtered_ML_data_no_duplicates.csv \
    --save_dir ./
```

# SaProt

```{bash}
# get SaProt predicted likelihood
$ python models/SaProt/get_SaProt_likelihood.py \
  --csv_file ./data/binding_affinity/aayl51_benchmarking_data.csv \
  --pdb_file ./data/complex_structure/AAYL51_bca.pdb 

```

# EpitopeSA

```{bash}
# activate environment
$ conda activate dyMEAN

# get total SASA
$ python scoring/dyMEAN/scripts/get_epitope_sasa.py \
    --native_pdb data/complex_structure/AAYL49.pdb \
    --antibody_chain B \
    --antigen_chain A \
    --dms_path data/binding_affinity/filtered_ML_data_no_duplicates.csv \
    --save_dir ./
```
