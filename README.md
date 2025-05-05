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
📂 The dataset used in this project is publicly available on [Hugging Face Datasets](https://huggingface.co/datasets/AbBibench/Antibody_Binding_Benchmark_Dataset). Please place the downloaded data in the data folder under the project root directory to ensure the program runs correctly.

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
