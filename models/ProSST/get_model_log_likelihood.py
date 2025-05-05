import os
import sys
import argparse
import pandas as pd
import numpy as np
import logging
from tqdm import tqdm
from scipy.stats import spearmanr

import torch
from transformers import AutoModelForMaskedLM, AutoTokenizer
from prosst.structure.quantizer import PdbQuantizer
from Bio import SeqIO
from Bio.SeqUtils import seq1
import json

def get_metadata(json_path="../data/metadata.json"):
    with open(json_path, "r") as fin:
        data_dict = json.load(fin)
    return data_dict


def extract_sequences_from_pdb(pdb_file, chains):
    """
    Extracts heavy, light, and antigen chain sequences from a PDB file using provided chain IDs.
    
    Args:
        pdb_file (str): Path to the PDB file.
        chains (list of str): List of chain IDs, where:
            - chains[0] is the heavy chain,
            - chains[1] is the light chain,
            - chains[2:] (if any) are antigen chains, concatenated together.
    
    Returns:
        tuple: (heavy_chain_seq, light_chain_seq, antigen_chain_seq)
    """
    from Bio import PDB
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("PDB", pdb_file)
    
    heavy_chain = chains[0] if len(chains) > 0 else None
    light_chain = chains[1] if len(chains) > 1 else None
    antigen_chains = chains[2:] if len(chains) > 2 else []
    
    heavy_seq = ""
    light_seq = ""
    antigen_seq = ""
    
    # Iterate over all models and chains in the structure
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

            # Check which chain this is and assign appropriately
            if heavy_chain and chain.id == heavy_chain:
                heavy_seq = seq
            elif light_chain and chain.id == light_chain:
                light_seq = seq
            elif chain.id in antigen_chains:
                antigen_seq += seq  # Concatenate if multiple antigen chains are provided
    return heavy_seq, light_seq, antigen_seq


def infer_chain_ids(pdb_file_path):
    """
    Infers the heavy, light, and antigen chain IDs from the PDB filename.
    Example: '3gbn_hlab.pdb' -> Heavy: H, Light: L, Antigen: ['A','B']
    """
    filename = os.path.basename(pdb_file_path)
    parts = filename.split('_')
    
    # If the file is "xxxx_hlab.pdb", parts[0] = "xxxx", parts[1] = "hlab.pdb"
    chain_part = parts[1].split('.')[0]  # e.g. 'hlab'
    
    heavy_chain_id = chain_part[0].upper()         # 'H'
    light_chain_id = chain_part[1].upper()         # 'L'
    antigen_chain_ids = list(chain_part[2:].upper())  # ['A','B'] if "hlab"
    
    print(f"Inferred chain IDs -> Heavy: {heavy_chain_id}, Light: {light_chain_id}, Antigen: {antigen_chain_ids}")
    return heavy_chain_id, light_chain_id, antigen_chain_ids


def setup_logging(log_file):
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )


def identify_mutations(wildtype_hc: str, mutated_hc: str) -> str:
    """
    Compares the wildtype heavy chain (HC) with the mutated heavy chain (HC) to identify mutation sites.
    Returns a string in the format: "A123B:C124D", which represents multiple point mutations.
    """
    mutations = []
    for i, (wt_aa, mut_aa) in enumerate(zip(wildtype_hc, mutated_hc), start=1):
        if wt_aa != mut_aa:
            mutations.append(f"{wt_aa}{i}{mut_aa}")
    return ":".join(mutations) if mutations else ""


def get_ll_average_complex_prosst(model, tokenizer, structure_seq, residue_seq, pdb_file):
    """
    Uses ProSST to perform a forward pass on the mutated residue_seq and computes the average log-likelihood
    of the entire sequence.

    Parameters:
    - model: ProSST model (AutoModelForMaskedLM)
    - tokenizer: The corresponding ProSST tokenizer
    - residue_seq: Protein sequence to be evaluated (str), which should be the mutated sequence
    - pdb_file: The corresponding PDB file path; if it only contains the wildtype structure, ensure it matches the mutated sequence as closely as possible

    Returns:
    - avg_ll: float, the average log-likelihood of the entire sequence
    """
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)

    # ProSST convention: add +3 to each token in the structural quantization result (due to [CLS], [SEP], [PAD] in the tokenizer)
    structure_seq_offset = [x + 3 for x in structure_seq]

    # 2. Tokenize the amino acid sequence
    tokenized_res = tokenizer([residue_seq], return_tensors='pt')
    input_ids = tokenized_res['input_ids'].to(device)
    attention_mask = tokenized_res['attention_mask'].to(device)

    # 3. Construct the input for the structure sequence: wrapped with [CLS]=1 and [SEP]=2
    structure_input_ids = torch.tensor(
        [1, *structure_seq_offset, 2], dtype=torch.long
    ).unsqueeze(0).to(device)

    # 4. Perform the forward pass
    with torch.no_grad():
        outputs = model(
            input_ids=input_ids,
            attention_mask=attention_mask,
            ss_input_ids=structure_input_ids
        )
    # outputs.logits: [batch_size, seq_len, vocab_size]

    # 5. Usually ProSST input includes additional special tokens at the beginning and end, so remove the first and last positions.
    #    logits.shape -> [1, seq_len, vocab_size] -> squeeze -> [seq_len, vocab_size]
    logits = outputs.logits[:, 1:-1, :].log_softmax(dim=-1).squeeze(0)
    # At this point, logits.shape[0] typically matches the length of residue_seq (assuming no extra special tokens added by the tokenizer)

    # If the tokenizer adds tokens other than [SEP] at the end, you need to verify that the sequence lengths align.
    # Here we assume it exactly corresponds to residue_seq.

    # 6. Compute the logP for the true amino acid at each position and accumulate the values.
    vocab = tokenizer.get_vocab()
    seq_len = len(residue_seq)
    log_probs = []
    for i in range(seq_len):
        aa = residue_seq[i]
        if aa not in vocab:
            # If a rare amino acid or an unrecognized character appears, special handling is required.
            logging.warning(f"Amino acid '{aa}' not in tokenizer vocab. Assign log_prob=-inf.")
            log_probs.append(float('-inf'))
            continue

        aa_idx = vocab[aa]  # The index of the amino acid in the vocab
        ll_i = logits[i, aa_idx].item()  # logP(aa | context, structure)
        log_probs.append(ll_i)

    if all(np.isinf(x) and x < 0 for x in log_probs):
        # If all values are -inf, it indicates that the sequence does not match the tokenizer.
        return None

    # 7. Compute the average log-likelihood of the entire sequence.
    avg_ll = float(np.mean(log_probs))
    return avg_ll


def main():
    parser = argparse.ArgumentParser(description="ProSST Average Log-Likelihood Benchmarking Script")
    parser.add_argument('--name', type=str, default='3gbn', choices=['3gbn','4fqi','2fjg','aayl49','aayl49_ml','aayl51','1mlc', '1n8z', '1mhp'])
    args = parser.parse_args()

    # Log file
    base_name = os.path.splitext(os.path.basename(args.name))[0]
    log_file = f"./{base_name}_prosst_avgll_log.log"
    setup_logging(log_file)

    # load affinity data
    if args.name == "aayl49_ml":
        name = "aayl49_ML"
    else:
        name = args.name
    pdb_info = get_metadata()[name]
    print('pdb_info:', pdb_info)
    df = pd.read_csv('.'+pdb_info["affinity_data"][0])
    pdb_info = get_metadata()[name]
    pdb_file = '.'+pdb_info["pdb_path"]
    output_file = f'../notebooks/scoring_outputs/{name}_benchmarking_data_prosst_scores.csv'


    # Load ProSST
    logging.info("Loading ProSST model and tokenizer...")
    prosst_model = AutoModelForMaskedLM.from_pretrained("AI4Protein/ProSST-2048", trust_remote_code=True)
    prosst_tokenizer = AutoTokenizer.from_pretrained("AI4Protein/ProSST-2048", trust_remote_code=True)
    processor = PdbQuantizer()  # Used for structural quantization
    # 1. Quantize structure (assuming pdb_file corresponds to residue_seq)
    structure_seq = processor(pdb_file)  # Returns a list of integers, each representing a quantized structural token

    # Infer chain IDs from the PDB filename
    heavy_chain_id, light_chain_id, antigen_chain_ids = infer_chain_ids(pdb_file)
    chains = [heavy_chain_id, light_chain_id] + antigen_chain_ids

    print("Extracting wild-type sequences from the provided PDB file.")
    wt_hc, wt_lc, wt_ag = extract_sequences_from_pdb(pdb_file, chains)

    # Add these sequences to the DataFrame
    # df['wt_heavy_chain_seq'] = wt_hc
    df['LC'] = wt_lc
    df['Target'] = wt_ag

    # Iterate over each record in the dataset and compute the average log-likelihood of the entire mutated sequence
    model_scores = []
    for idx, row in tqdm(df.iterrows(), total=df.shape[0], desc="Processing rows"):
        mutated_hc = row['mut_heavy_chain_seq'].replace(';', '').upper()
        mutated_lc = row['LC'].replace(';', '').upper()  # How to retrieve this; the current data format no longer provides it
        mutated_ag = row['Target'].replace(';', '').upper()  # How to retrieve this; the current data format no longer provides it
        # Concatenate HC + LC + antigen to form the complete mutated sequence
        mutated_full_seq = mutated_hc + mutated_lc + mutated_ag

        try:
            avg_ll = get_ll_average_complex_prosst(
                model=prosst_model,
                tokenizer=prosst_tokenizer,
                structure_seq=structure_seq,
                residue_seq=mutated_full_seq,
                pdb_file=pdb_file
            )
            model_scores.append(avg_ll)
            logging.info(f"Row {idx} - Average LL: {avg_ll}")
        except Exception as e:
            logging.error(f"Error calculating fitness for row {idx}: {e}")
            model_scores.append(None)

    # Write the results into a new column, for example named "ProSST_joint_ll_average"
    df['ProSST_joint_ll_average'] = model_scores

    # Output to Excel
    try:
        df.to_csv(output_file, index=False)
        logging.info(f"Updated dataset saved to {output_file}")
    except Exception as e:
        logging.error(f"Failed to save updated dataset to {output_file}: {e}")

    # # Compute correlation (Spearman); assuming the table contains a column 'neg_log_Kd'
    # if 'binding_score' in df.columns:
    #     df_valid = df.dropna(subset=['ProSST_joint_ll_average', 'binding_score'])
    #     logging.info(f"Number of valid rows for correlation: {len(df_valid)}")
    
    #     if len(df_valid) > 1:
    #         spearman_corr, p_value = spearmanr(df_valid['ProSST_joint_ll_average'], df_valid['binding_score'])
    #         logging.info(f"Spearman correlation: {spearman_corr:.2f}, p-value: {p_value:.2e}")
    #         print(f"Spearman correlation between ProSST_joint_ll_average and neg_log_Kd: {spearman_corr:.2f}, p-value: {p_value:.2e}")
    
    #         # Calculate Top-K Recall
    #         k = 10
    #         try:
    #             df_sorted = df_valid.sort_values(by='ProSST_joint_ll_average', ascending=False)
    #             top_k_predicted = df_sorted.head(k).index
    #             df_actual_sorted = df_valid.sort_values(by='binding_score', ascending=True)
    #             top_k_actual = df_actual_sorted.head(k).index
    #             top_k_recall = len(set(top_k_predicted).intersection(set(top_k_actual))) / k
    #             logging.info(f"Top-{k} Recall: {top_k_recall:.2f}")
    #             print(f"Top-{k} Recall: {top_k_recall:.2f}")
    #         except Exception as e:
    #             logging.error(f"Error calculating Top-{k} Recall: {e}")
    #             print(f"Error calculating Top-{k} Recall: {e}")
    #     else:
    #         print("Not enough valid data points to compute correlation.")
    # else:
    #     print("Column 'neg_log_Kd' not found in the dataset. Skipping correlation.")


if __name__ == '__main__':
    main()

