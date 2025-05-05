# 🧪 AbBiBench: Antibody Binding Benchmarking

This is the code for **AbBiBench** (*Anti*body *Bi*nding *Bench*marking), a benchmarking framework for optimizing antibody binding affinity. We use experimental antibody–antigen binding affinity measurements to evaluate the performance of widely used computational models for antibody sequence engineering, including **ESM-2**, **AntiBERTy**, **CurrAb**, **SaProt**, **ProSST**, **ESM-3**, **ProGen2**, **ProtGPT2**, **ProteinMPNN**, **ESM-IF**, **Antifold**, **DiffAb**, **MEAN**, and **dyMEAN**. We also compare several commonly used physics-based metrics, such as **−ΔG** and **−SASA**.


# Leaderboard
| Rank | Model Type        | Model          | 1mhp  | 1mlc  | 1n8z  | 2fjg  | 3gbn_h1 | 3gbn_h9 | 4fqi_h1 | 4fqi_h3 | aayl49 | aayl49_ML | aayl51 | Avg. Spearman ↑ |
|------|-------------------|----------------|-------|-------|-------|-------|---------|---------|----------|----------|--------|------------|--------|----------------|
| 🥇 1 | Inverse Folding   | ProteinMPNN    | -0.28 | -0.33 | -0.17 | 0.49  | 0.59    | 0.64    | 0.61     | 0.42     | 0.4    | 0.34       | 0.32   | 0.275454545    |
| 🥈 2 | Inverse Folding   | ESMIF1         | -0.36 | -0.36 | -0.11 | 0.55  | 0.59    | 0.54    | 0.65     | 0.49     | 0.39   | 0.27       | 0.34   | 0.271818182    |
| 🥉 3 | Biophysics        | -ΔG             | 0.24  | -0.05 | -0.28 | -0.01 | 0.59    | 0.32    | 0.64     | 0.29     | -0.01  | 0.24       | 0.07   | 0.185454545    |
| 4    | Inverse Folding   | Antifold       | -0.27 | -0.44 | 0.16  | 0.44  | 0.12    | 0.27    | 0.42     | 0.37     | 0.39   | 0.14       | 0.32   | 0.174545455    |
| 5    | Diffusion         | diffab_fixbb   | -0.14 | 0.01  | 0.04  | -0.05 | 0.54    | 0.76    | 0        | 0        | 0.18   | -0.01      | 0.19   | 0.138181818    |
| 6    | Masked LM         | SaProt         | -0.23 | 0.23  | 0.11  | -0.16 | 0.53    | 0.6     | 0.48     | 0.28     | -0.27  | 0.11       | -0.17  | 0.137272727    |
| 7    | Diffusion         | diffab         | -0.16 | 0.03  | -0.02 | -0.01 | 0.67    | 0.61    | 0        | -0.01    | 0.15   | 0          | 0.03   | 0.117272727    |
| 8    | Masked LM         | CurrAb         | 0.28  | 0.16  | -0.39 | 0.17  | 0.16    | 0.23    | 0.19     | 0.13     | 0.03   | 0.2        | 0.04   | 0.109090909    |
| 9    | Masked LM         | ESM2           | 0.13  | 0.02  | -0.13 | 0.08  | 0.23    | 0.38    | -0.02    | -0.02    | -0.04  | -0.14      | -0.11  | 0.034545455    |
| 10   | Graph Model       | dyMEAN_fixbb   | 0.28  | 0.03  | 0     | -0.02 | -0.02   | 0       | 0.04     | 0.02     | 0.02   | 0.02       | -0.02  | 0.031818182    |
| 11   | Graph Model       | dyMEAN         | -0.07 | 0.02  | 0     | 0.03  | -0.02   | -0.02   | 0.03     | 0.02     | -0.03  | -0.01      | -0.03  | -0.007272727   |
| 12   | Masked LM         | ESM3           | -0.26 | -0.3  | -0.12 | 0.13  | -0.24   | -0.22   | -0.2     | 0.03     | 0.39   | 0.12       | 0.26   | -0.037272727   |
| 13   | Masked LM         | ProSST         | -0.21 | 0.01  | -0.26 | 0.09  | -0.3    | -0.07   | -0.07    | 0.1      | 0.13   | -0.01      | 0.11   | -0.043636364   |
| 14   | Autoregressive LM | ProtGPT2       | -0.1  | -0.21 | 0.15  | 0.04  | -0.39   | -0.18   | -0.2     | 0        | 0.06   | 0.05       | 0.1    | -0.061818182   |
| 15   | Biophysics        | -SASA          | -0.17 | 0.1   | 0.17  | -0.02 | -0.26   | -0.2    | -0.14    | -0.15    | 0.05   | -0.18      | 0.02   | -0.070909091   |
| 16   | Graph Model       | MEAN_fixbb     | -0.21 | 0.02  | 0.16  | -0.18 | -0.2    | -0.04   | -0.36    | -0.21    | 0.06   | 0.02       | -0.05  | -0.09          |
| 17   | Graph Model       | MEAN           | -0.21 | 0.02  | 0.15  | -0.18 | -0.24   | 0       | -0.6     | -0.28    | 0.07   | 0.02       | -0.05  | -0.118181818   |
| 18   | Masked LM         | AntiBERTy      | -0.18 | -0.24 | -0.24 | 0.11  | -0.72   | -0.75   | -0.38    | -0.2     | 0.21   | -0.14      | 0.22   | -0.21          |
| 19   | Autoregressive LM | progen2-large  | -0.38 | -0.29 | -0.21 | 0.27  | -0.76   | -0.62   | -0.45    | -0.32    | 0.26   | -0.11      | 0.2    | -0.219090909   |

Each value in this table indicates the Spearman correlation between the model's predicted log-likelihood scores and the corresponding experimental measurement from a specific antibody–antigen dataset. They are ranked according to the average Spearman correlation coefficient across multiple datasets.

# Installation

We recommend create a conda environment for each tool:

```{bash}

$ conda env create --name ENV_NAME --file envs/ENV_FILE.yml

```
We have provided requirement files for each tools in __envs__ directory, including `diffab.yml`, `dyMEAN.yml`,
`esmif.yml`, `MEAN_ProteinMPNN.yml`, `prosst.yml`, `SaProt.yml`

# Data Resource
📂 The dataset used in this project is publicly available on [Hugging Face Datasets](https://huggingface.co/datasets/AbBibench/Antibody_Binding_Benchmark_Dataset).

# Model log-likelihood scoring

## Run the Script

```{bash}

cd ./scripts
python eval_seq.py --model [MODEL] --data [DATA]

```
Where MODEL ∈ { diffab, ESM-IF, AntiFold, ESM-2, ESM3-Open, AntiBERTy, CurrAb, dyMEAN, MEAN, ProteinMPNN, ProSST, ProGen2, ProSST, foldx, sasa }, and DATA ∈ { 3gbn, 4fqi, 2fjg, aayl49, aayl49_ml, aayl51, 1mlc, 1n8z , 1mph}

## Example

```{bash}

cd ./scripts
python eval_seq.py --model diffab --data 3gbn

```
This will:
1. Activate the Conda environment diffab.
2. Run models/diffab/get_model_log_likelihood.py --name 3gbn.
3. Save the output to: benchmark/notebooks/scoring_outputs/3gbn_benchmarking_data_diffab_scores.csv

# Correlation to antibody-antigen binding affinity
 
We provide a Jupyter Notebook in __notebooks/figure.ipynb__ to reproduce our correlation results shown in our paper.

