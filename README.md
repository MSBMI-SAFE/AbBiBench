# 🧪 AbBiBench: Antibody Binding Benchmarking

This is the code for **AbBiBench** (*Anti*body *Bi*nding *Bench*marking), a benchmarking framework for optimizing antibody binding affinity. We use experimental antibody–antigen binding affinity measurements to evaluate the performance of widely used computational models for antibody sequence engineering, including **ESM-2**, **AntiBERTy**, **CurrAb**, **SaProt**, **ProSST**, **ESM-3**, **ProGen2**, **ProtGPT2**, **ProteinMPNN**, **ESM-IF**, **Antifold**, **DiffAb**, **MEAN**, **dyMEAN**, **AF3**, and **Boltz-2**. We also compare several commonly used physics-based metrics, such as **−ΔG** and **−SASA**.


# Leaderboard
| Rank | Model Type           | Model               | 1mhp  | 1mlc  | 1n8z  | 2fjg  | 3gbn_h1 | 3gbn_h9 | 4fqi_h1 | 4fqi_h3 | 4d5_her2 | 5a12_ang2 | 5a12_vegf | aayl50 | aayl49 | aayl49_ML | aayl51 | aayl52 | Avg. Spearman ↑ |
|------|----------------------|---------------------|-------|-------|-------|-------|---------|---------|---------|---------|----------|-----------|-----------|--------|--------|-----------|--------|--------|----------------|
| 🥇 1 | Inverse Folding      | ProteinMPNN         | -0.02 | -0.21 | -0.17 | 0.46  | 0.59    | 0.64    | 0.61    | 0.42    | 0.32     | 0.13      | 0.54      | 0.19   | 0.40   | 0.34      | 0.32   | 0.23   | 0.30           |
| 🥈 2 | Inverse Folding      | ESMIF1              | 0.01  | -0.31 | -0.11 | 0.49  | 0.59    | 0.54    | 0.65    | 0.49    | 0.40     | 0.13      | 0.24      | 0.14   | 0.39   | 0.27      | 0.34   | 0.25   | 0.28           |
| 🥉 3 | Inverse Folding      | Antifold            | -0.02 | -0.31 | 0.16  | 0.41  | 0.12    | 0.27    | 0.42    | 0.37    | 0.34     | 0.18      | 0.21      | 0.07   | 0.39   | 0.14      | 0.32   | 0.24   | 0.21           |
| 4    | Structure Prediction | Boltz-2             | -0.22 | 0.02  | -0.06 | 0.08  | 0.71    | 0.56    | 0.40    | 0.27    | 0.03     | -0.02     | 0.31      | -0.02  | 0.02   | 0.02      | 0.05   | -0.04  | 0.13           |
| 5    | Biophysics           | FoldX               | -0.02 | 0.02  | -0.28 | -0.02 | 0.59    | 0.32    | 0.64    | 0.29    | 0.01     | -0.04     | 0.11      | -0.02  | -0.01  | 0.24      | 0.07   | 0.09   | 0.12           |
| 6    | Diffusion            | diffab              | 0.24  | 0.01  | -0.02 | -0.01 | 0.67    | 0.61    | 0.00    | -0.01   | -0.01    | 0.03      | 0.00      | 0.02   | 0.15   | 0.00      | 0.03   | 0.01   | 0.11           |
| 7    | Diffusion            | diffab_fixbb        | -0.09 | 0.01  | 0.04  | -0.02 | 0.54    | 0.76    | 0.00    | 0.00    | -0.01    | 0.02      | 0.01      | 0.00   | 0.18   | -0.01     | 0.19   | 0.00   | 0.10           |
| 8    | Masked LM            | CurrAb              | 0.12  | 0.11  | -0.39 | 0.14  | 0.16    | 0.23    | 0.19    | 0.13    | 0.14     | 0.05      | 0.03      | 0.01   | 0.03   | 0.20      | 0.04   | 0.01   | 0.07           |
| 9    | Masked LM            | ESM3-Open-structure | 0.06  | -0.28 | -0.12 | 0.17  | -0.24   | -0.22   | -0.20   | 0.03    | 0.25     | 0.16      | 0.13      | 0.09   | 0.39   | 0.12      | 0.26   | 0.22   | 0.05           |
| 10   | Masked LM            | SaProt              | -0.15 | 0.25  | 0.11  | -0.21 | 0.53    | 0.60    | 0.48    | 0.28    | -0.34    | -0.08     | -0.14     | -0.11  | -0.27  | 0.11      | -0.17  | -0.27  | 0.04           |
| 11   | Masked LM            | ESM2                | 0.20  | 0.06  | -0.13 | 0.14  | 0.23    | 0.38    | -0.02   | -0.02   | -0.20    | -0.07     | 0.20      | -0.03  | -0.04  | -0.14     | -0.11  | -0.16  | 0.02           |
| 12   | Graph Model          | dyMEAN_fixbb        | 0.15  | 0.01  | 0.00  | -0.02 | -0.02   | 0.00    | 0.04    | 0.02    | -0.02    | -0.01     | 0.02      | -0.01  | 0.02   | 0.02      | -0.02  | 0.02   | 0.01           |
| 13   | Autoregressive LM    | ProtGPT2            | 0.14  | -0.21 | 0.15  | 0.17  | -0.39   | -0.18   | -0.20   | 0.00    | 0.05     | 0.08      | 0.18      | 0.02   | 0.06   | 0.05      | 0.10   | -0.06  | 0.00           |
| 14   | Masked LM            | ProSST              | 0.09  | 0.02  | -0.26 | 0.08  | -0.30   | -0.07   | -0.07   | 0.10    | 0.07     | 0.16      | -0.06     | -0.02  | 0.13   | -0.01     | 0.11   | -0.03  | 0.00           |
| 15   | Graph Model          | dyMEAN              | -0.08 | 0.02  | 0.00  | 0.01  | -0.02   | -0.02   | 0.03    | 0.02    | 0.03     | 0.01      | 0.01      | -0.02  | -0.03  | -0.01     | -0.03  | 0.00   | -0.01          |
| 16   | Graph Model          | MEAN_fixbb          | -0.07 | 0.01  | 0.16  | -0.12 | -0.20   | -0.04   | -0.36   | -0.21   | 0.06     | 0.02      | 0.34      | 0.01   | 0.06   | 0.02      | -0.05  | 0.02   | -0.02          |
| 17   | Structure Prediction | AF3                 | -0.54 | -0.17 | -0.16 | 0.09  | -0.05   | 0.13    | 0.05    | 0.05    | 0.02     | -0.11     | -0.01     | -0.01  | 0.23   | 0.00      | 0.07   | 0.01   | -0.02          |
| 18   | Biophysics           | epitopeSA           | 0.08  | 0.09  | 0.17  | 0.02  | -0.26   | -0.20   | -0.14   | -0.15   | 0.00     | -0.03     | 0.12      | 0.01   | 0.05   | -0.18     | 0.02   | -0.18  | -0.04          |
| 19   | Graph Model          | MEAN                | 0.01  | 0.01  | 0.15  | -0.12 | -0.24   | 0.00    | -0.60   | -0.28   | 0.02     | 0.01      | 0.16      | 0.02   | 0.07   | 0.02      | -0.05  | 0.02   | -0.05          |
| 20   | Autoregressive LM    | progen2-large       | -0.01 | -0.28 | -0.21 | 0.26  | -0.76   | -0.62   | -0.45   | -0.32   | 0.07     | 0.15      | 0.09      | 0.11   | 0.26   | -0.11     | 0.20   | 0.20   | -0.09          |
| 21   | Masked LM            | AntiBERTy           | 0.04  | -0.13 | -0.24 | 0.01  | -0.72   | -0.75   | -0.38   | -0.20   | 0.13     | 0.17      | 0.02      | 0.04   | 0.21   | -0.14     | 0.22   | 0.01   | -0.11          |
<!-- | -    | Biophysics           | ANTIPASTI           | -0.18 | -0.05 | -0.43 | -     | -       | -       | -       | -       | -        | -         | -         | -      | 0.07   | -         | -0.01  | -      | -              |
| -    | Diffusion            | AbX                 | -0.16 | -0.24 | 0.14  | 0.11  | -0.63   | -0.47   | -0.53   | -0.34   | -        | -         | -         | -      | 0.16   | 0.01      | 0.11   | -      | -              | -->



Each value in this table indicates the Spearman correlation between the model's predicted log-likelihood scores and the corresponding experimental measurement from a specific antibody–antigen dataset. They are ranked according to the average Spearman correlation coefficient across multiple datasets.

# Installation

We recommend create a conda environment for each tool:

```{bash}

$ conda env create --name ENV_NAME --file envs/ENV_FILE.yml

```
We have provided requirement files for each tools in __envs__ directory, including `diffab.yml`, `dyMEAN.yml`,
`esmif.yml`, `MEAN_ProteinMPNN.yml`, `prosst.yml`, `SaProt.yml`

# Data Resource
📂 The dataset used in this project is publicly available on [Hugging Face Datasets](https://huggingface.co/datasets/AbBibench/Antibody_Binding_Benchmark_Dataset). Please place the downloaded data in the data folder under the project root directory to ensure the program runs correctly.

The latest AbBiBench dataset can be easily loaded via Hugging face. Below is an example that demonstrates the entire workflow—from listing and loading data, to filtering by antigen and downloading/parsing PDB files. We will also update this example in our GitHub repository and provide a PyTorch dataset version:

```python
from huggingface_hub import list_repo_files, hf_hub_download
from datasets import load_dataset, concatenate_datasets
from tqdm import tqdm
import biotite.structure.io as bsio

REPO = "AbBibench/Antibody_Binding_Benchmark_Dataset"

# 1. List all CSV files in the binding_affinity directory
csv_files = [
    f for f in list_repo_files(REPO, repo_type="dataset")
    if f.startswith("binding_affinity/") and f.endswith("_benchmarking_data.csv")
]

# 2. Load and concatenate all subsets
all_splits = []
for csv in tqdm(csv_files, desc="Loading CSVs"):
    ds = load_dataset(REPO, data_files={ "data": csv }, split="train")
    all_splits.append(ds)
full_ds = concatenate_datasets(all_splits)
print(full_ds)    # overview of the full dataset

# 3. Filter for samples belonging to influenza H1 (3gbn_h1)
h1_ds = full_ds.filter(lambda x: x["antigen_id"].endswith("3gbn_h1"))

# 4. List PDB structure files corresponding to this antigen
antigen_id     = "3gbn_h1"
base_id        = antigen_id.split("_")[0]
structure_files = [
    f for f in list_repo_files(REPO, repo_type="dataset")
    if f.startswith(f"structures/{base_id}") and f.endswith(".pdb")
]

# 5. Download and parse each PDB using Biotite
for pdb_file in structure_files:
    local_pdb = hf_hub_download(
        repo_id=REPO, filename=pdb_file, repo_type="dataset"
    )
    print("Downloaded to:", local_pdb)
    atom_array = bsio.load_structure(local_pdb)
    print("Chains:", atom_array.chain_id)
```

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

# 🏆 Contribute to the AbBiBench Leaderboard — We Welcome Your Model and Data!

We maintain a public **AbBiBench** leaderboard and **actively invite external submissions** that benchmark new models or datasets for antibody–antigen binding affinity.

---


## 🚀 Step‑by‑step guide for submitting model results

1. **Fork** this repository and create a new branch:

   ```bash
   git clone https://github.com/<your_username>/AbBiBench.git
   cd AbBiBench
   git checkout -b leaderboard-<your_model>
   ```

2. **Add your code and results**

    | Requirement | Details |
    |-------------|---------|
    | **Project layout** | Place all evaluation code inside **`models/<your_model>/`**. |
    | **CLI interface** | Your main script must accept **`--name $name`** (dataset name). |
    | **Output format** | For each mutant, write a CSV of scores to **`notebooks/scoring_outputs/`**. |
    | **Environment** | Put any `environment.yml` or `requirements.txt` in **`envs/`**. |
    | **Leaderboard row** | Append **one line** to `leaderboard/leaderboard.csv` (preserve column order). |

3. **Commit and push**

   ```bash
   git add models/<your_model> envs/ notebooks/scoring_outputs/<file>.csv README.md
   git commit -m "Leaderboard submission: <your_model>"
   git push -u origin leaderboard-<your_model>
   ```

4. **Open a Pull Request** to `master`  

   Title your PR:

   ```
   Leaderboard submission: <Your Model Name>
   ```

   and include the following template in the PR description:

   ```markdown
   ### Method name
   <Your model>

   ### Short description (≤ 100 words)
   …

   ### Reference
   arXiv / DOI / blog link (optional)

   ### Reproduction command
   python models/<your_model>/run.py --name 1mhp
   ```

5. **Review and merge**  
   We will verify your scores and code within ~7 days. Once merged, your model will appear automatically on the leaderboard.


## 📦 Contribute Data to `AbBibench/Antibody_Binding_Benchmark_Dataset`

We warmly welcome community contributions of new **antibody–antigen binding affinity datasets** to the AbBiBench benchmark on the Hugging Face Hub.  
Data **must be shared under an open license** (CC‑BY‑4.0 or a compatible license).

---


1. **Install Git LFS and sign in to Hugging Face**

   ```bash
   conda install -c conda-forge git-lfs
   git lfs install        # one‑time setup
   pip install -U huggingface_hub
   huggingface-cli login  # paste your HF access token
   ```

2. **Fork and clone the dataset repo**

   ```bash
   # Replace <username> with your HF account
   git clone https://huggingface.co/datasets/<username>/Antibody_Binding_Benchmark_Dataset
   cd Antibody_Binding_Benchmark_Dataset
   git remote add upstream https://huggingface.co/datasets/AbBibench/Antibody_Binding_Benchmark_Dataset
   git pull upstream main   # stay up to date
   ```

3. **Add your data**
   
   - Each **CSV inside `binding_affinity/`** must include at least:
   
      | column      | description                                                          |
      |-------------|----------------------------------------------------------------------|
      | `mut_heavy_chain_seq`  | Amino‑acid sequence for each mutant of heavy chain                                               |
      | `binding_score`  | Experimental affinity value |
   
   - Place every PDB/mmCIF file inside `complex_structure/`. 
   
   - Each study **must** provide a `metadata.json` at the root of its folder. The file should be a **dictionary keyed by complex ID** (typically the PDB code). For each complex include the fields below:
   
      | key            | type / example | description |
      |----------------|----------------|-------------|
      | `pdb`          | `"1mhp_hla"`   | PDB identifier (or custom) |
      | `pdb_path`     | `"./data/complex_structure/1mhp_hla.pdb"` | Relative path to the structure file |
      | `heavy_chain`  | `"H"`          | Heavy chain ID of the antibody |
      | `light_chain`  | `"L"`          | Light chain ID of the antibody |
      | `antigen_chains` | `["A"]`      | Antigen chain IDs |
      | `affinity_data`  | `["./data/binding_affinity/1mhp_benchmarking_data.csv"]` | Paths to corresponding affinity CSV files |
      | `receptor_chains` | `["A"]`     | Chains treated as receptor in docking (if applicable) |
      | `ligand_chains`   | `["H","L"]` | Chains treated as ligand in docking |
      | `chain_order`     | `["H","L","A"]` | Ordering of chains in the complex file |
      | `epitope_chain`   | `"A"`       | Chain containing the epitope residues |
      | `paratope_chain`  | `"H"`       | Chain containing the paratope residues |


4. **Commit and push**

   ```bash
   git checkout -b add-<your_study_name>
   git add data/<your_dataset>.csv metadata.json
   git commit -m "Add <your_study_name> dataset (n=1234 mutants)"
   git push -u origin add-<your_study_name>
   ```

5. **Open a Pull Request on the HF Hub**

   Use **Contribute → Pull request** on the repo page and fill out:

   ```markdown
   ### Study name
   <your_study_name>

   ### Description (≤ 100 words)
   Short summary of the experiment, antigen, number of mutants, and assay.

   ### Files added
   - data/<your_study_name>/binding_affinity/*.csv
   - data/<your_study_name>/complex_structure/*.pdb
   - …

   ### License
   CC-BY-4.0
   ```

We will review your PR—checking format, license, and basic biological plausibility—within **about 7 days**. Once merged, your data will appear in the next dataset snapshot and can be used immediately by AbBiBench.


🙏 **Thanks for contributing and helping improve antibody‑design benchmarks!**
