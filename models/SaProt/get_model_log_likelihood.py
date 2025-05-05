import argparse
import os
import sys
import json
import pandas as pd
import numpy as np
from transformers import AutoModelForMaskedLM, AutoTokenizer
import torch
from tqdm import tqdm
from Bio.SeqUtils import seq1

sys.path.append('./models/SaProt/utils')
from foldseek_util import get_struc_seq

foldseek_struc_vocab = "pynwrqhgdlvtmfsaeikc#"


def parse_args():
    parser = argparse.ArgumentParser(description="SaProt Benchmarking with JSON Metadata")
    parser.add_argument("--name", required=True,
                        help="Key name from the JSON file (e.g., 3gbn, 4fqi, etc.)")
    parser.add_argument("--json_file", default="./data/metadata.json",
                        help="Path to the JSON file containing dataset metadata")
    parser.add_argument("--gpu", type=int, default=0,
                        help="Which GPU to use; set -1 for CPU only.")
    parser.add_argument("--output_dir", default="./notebooks/scoring_outputs",
                        help="Directory for output CSV files")
    # foldseek_bin can come from JSON or override here:
    parser.add_argument("--foldseek_bin", default=None,
                        help="Path to foldseek binary (if not provided in JSON).")
    return parser.parse_args()

def extract_sequences_from_pdb(pdb_file, chains):
    """
    Extract heavy, light, and antigen chain sequences from a PDB file.
    chains[0] = heavy, chains[1] = light, chains[2..] = antigen(s)
    """
    from Bio import PDB
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("PDB", pdb_file)

    heavy_chain = chains[0] if len(chains) > 0 else None
    light_chain = chains[1] if len(chains) > 1 else None
    antigen_chains = chains[2:] if len(chains) > 2 else []

    heavy_seq, light_seq, antigen_seq = "", "", ""

    for model in structure:
        for chain in model:
            seq = ""
            for residue in chain:
                if residue.get_resname() in PDB.Polypeptide.standard_aa_names:
                    try:
                        aa = seq1(residue.get_resname())
                    except Exception as e:
                        print(f"Error converting residue {residue.get_resname()}: {e}")
                        aa = "X"
                    seq += aa

            # Assign to the correct chain
            if heavy_chain and chain.id == heavy_chain:
                heavy_seq = seq
            elif light_chain and chain.id == light_chain:
                light_seq = seq
            elif chain.id in antigen_chains:
                antigen_seq += seq  # concatenate multiple antigens
    return heavy_seq, light_seq, antigen_seq

def get_ll_average_complex(model, tokenizer, seq: str, mut_info: str) -> float:
    """
    Calculate the average log-likelihood of the entire mutated protein complex, including non-mutated positions.
    
    Args:
        model: The SaProt model.
        tokenizer: The SaProt tokenizer.
        seq (str): The wild-type sequence concatenated with structural sequence.
        mut_info (str): Mutation information in the format "A123B" or "A123B:C124D".
    
    Returns:
        float: The average log-likelihood score of the entire complex.
    """

    tokens = tokenizer.tokenize(seq)

    for single in mut_info.split(":"):
        wt = single[0]
        pos = int(single[1:-1])
        mut_aa = single[-1]

        # Check position
        if pos - 1 >= len(tokens):
            print(f"Mutation position {pos} out of bounds for sequence len {len(tokens)}.")
            return None

        token_aa = tokens[pos - 1][0]
        if token_aa != wt:
            print(f"Mismatch at pos {pos}: expected {wt}, found {token_aa}.")
            return None

        # Replace token with mutated token
        tokens[pos - 1] = mut_aa + tokens[pos - 1][-1]

    mutated_seq = " ".join(tokens)
    inputs = tokenizer(mutated_seq, return_tensors="pt")
    inputs = {k: v.to(model.device) for k, v in inputs.items()}

    with torch.no_grad():
        outputs = model(**inputs)
        logits = outputs.logits
        probs = logits.softmax(dim=-1)

    vocab = tokenizer.get_vocab()
    ll_scores = []

    for pos, token in enumerate(tokens):
        if token not in vocab:
            print(f"Token '{token}' not in vocab at pos {pos+1}.")
            return None
        idx = vocab[token]
        prob = probs[0, pos, idx]
        ll_scores.append(float(torch.log(prob)) if prob > 0 else float("-inf"))

    return np.mean(ll_scores)

def calc_fitness(foldseek_bin, model, data, tokenizer,
                 mutation_col='mutations', target_seq=None, pdb_file=None, chains=None):
    """
    Calculate fitness scores for a set of mutations by averaging log-likelihoods
    over the entire complex. 
    """
    # 1) Use foldseek to get structural sequence
    struc_seq_dict = get_struc_seq(foldseek_bin, pdb_file, chains,
                                   plddt_mask=False, plddt_threshold=70)

    # Concatenate the structure sequences
    struc_seq = ""
    for chain in chains:
        if chain in struc_seq_dict:
            seq_for_chain = struc_seq_dict[chain][1]  # (pdb_id, seq)
            struc_seq += seq_for_chain.lower()
        else:
            print(f"Warning: chain {chain} not found in struct seq dict. Using ''.")
            struc_seq += ""

    # Merge target_seq with structural sequence (character-by-character)
    # => final seq length = length of target_seq
    # e.g. if target_seq is "ABCDE" and struc_seq is "abcde", you form "Aa Bb Cc Dd Ee"
    seq = "".join([a + b for a, b in zip(target_seq, struc_seq)])

    log_proba_list = []

    for mut_info in tqdm(data[mutation_col], desc="Calculating fitness scores"):
        score = get_ll_average_complex(model, tokenizer, seq, mut_info)
        log_proba_list.append(score)

    return np.array(log_proba_list)

def identify_mutations(wt_hc: str, mutated_hc: str) -> str:
    """
    Identify all differences between WT heavy chain and mutated heavy chain,
    returning something like "A123B:C200E".
    """
    mutations = []
    for i, (wt_aa, mut_aa) in enumerate(zip(wt_hc, mutated_hc), start=1):
        if wt_aa != mut_aa:
            mutations.append(f"{wt_aa}{i}{mut_aa}")
    return ":".join(mutations) if mutations else ""

def main():
    args = parse_args()

    # 1) Load JSON metadata
    with open(args.json_file, "r") as f:
        dataset_info = json.load(f)

    # 2) Check name
    if args.name not in dataset_info:
        raise ValueError(f"'{args.name}' not found in {args.json_file}.")

    metadata = dataset_info[args.name]
    pdb_path = metadata["pdb_path"]
    heavy_chain = metadata["heavy_chain"]
    light_chain = metadata["light_chain"]
    antigen_chains = metadata["antigen_chains"]
    # chain_order = metadata["chain_order"]  # (If needed, but not used here)
    affinity_data_files = metadata["affinity_data"]
    foldseek_bin = args.foldseek_bin if args.foldseek_bin else metadata.get("foldseek_bin", "./models/SaProt/bin/foldseek")

    # Make sure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # 3) Load SaProt model
    print("Loading SaProt model...")
    tokenizer = AutoTokenizer.from_pretrained("westlake-repl/SaProt_650M_AF2")
    model = AutoModelForMaskedLM.from_pretrained("westlake-repl/SaProt_650M_AF2")
    model.eval()

    # 4) GPU?
    if torch.cuda.is_available() and args.gpu >= 0:
        torch.cuda.set_device(args.gpu)
        model = model.cuda()
        print(f"Model loaded onto GPU {args.gpu}")
    else:
        print("Running on CPU.")

    # 5) Loop over each CSV
    for csv_file in affinity_data_files:
        if not os.path.exists(csv_file):
            print(f"CSV {csv_file} not found, skipping.")
            continue

        df = pd.read_csv(csv_file)
        print(f"Loaded {len(df)} rows from {csv_file}")

        # 5a) Extract wild-type sequences from PDB
        if not os.path.exists(pdb_path):
            print(f"PDB file {pdb_path} not found, skipping.")
            continue

        # (Use the function we defined)
        wt_hc, wt_lc, wt_ag = extract_sequences_from_pdb(
            pdb_file=pdb_path,
            chains=[heavy_chain, light_chain] + antigen_chains
        )
        print(f"WT Heavy chain length: {len(wt_hc)}, Light: {len(wt_lc)}, Antigen(s): {len(wt_ag)}")

        # 5b) Ensure 'mut_heavy_chain_seq' is present
        if "mut_heavy_chain_seq" not in df.columns:
            print(f"'mut_heavy_chain_seq' col missing in {csv_file}, skipping.")
            continue

        # 5c) Build 'mutations' column
        df["mutations"] = df.apply(
            lambda row: identify_mutations(wt_hc, row["mut_heavy_chain_seq"])
                        if isinstance(row["mut_heavy_chain_seq"], str) else None,
            axis=1
        )

        # Build combined WT sequence for foldseek
        target_seq = (wt_hc.upper() + wt_lc.upper() + wt_ag.upper())

        # 5d) Score each row
        fitness_scores = calc_fitness(
            foldseek_bin=foldseek_bin,
            model=model,
            data=df,
            tokenizer=tokenizer,
            mutation_col="mutations",
            target_seq=target_seq,
            pdb_file=pdb_path,
            chains=[heavy_chain, light_chain] + antigen_chains
        )

        df["log-likelihood"] = fitness_scores

        # 5e) Save output
        basename = os.path.splitext(os.path.basename(csv_file))[0]
        output_file = os.path.join(
            args.output_dir, f"{basename}_SaProt_scores.csv"
        )
        df.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()
