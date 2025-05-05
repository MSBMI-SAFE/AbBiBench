import os
from Bio import PDB
from Bio.SeqUtils import seq1
import json

def extract_sequences_from_pdb(pdb_file, chains):
    """
    Extract chain sequences from a PDB file.

    Parameters:
      pdb_file (str): Path to the PDB file.
      chains (list): A list containing chain identifiers in the order:
                     [heavy_chain, light_chain, *antigen_chain(s)].
                     
    Returns:
      A dictionary mapping chain ids to sequence strings.
      For example:
          {
              "heavy": "EVQLVQSG...",
              "light": "DVVMTQSS...",
              "antigen": {"A": "MISSS...", "B": "LLKDS..."}
          }
      If a chain is not found, its entry is set to an empty string.
    """
    parser = PDB.PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("PDB_structure", pdb_file)
    except Exception as e:
        print(f"Error parsing {pdb_file}: {e}")
        return None

    # Prepare the dict for storing sequences.
    sequences = {}
    heavy_chain = chains[0] if len(chains) > 0 else None
    light_chain = chains[1] if len(chains) > 1 else None
    antigen_chains = chains[2:] if len(chains) > 2 else []

    sequences["heavy"] = ""
    sequences["light"] = ""
    sequences["antigen"] = {chain_id: "" for chain_id in antigen_chains}

    # Iterate over all models and chains in the structure.
    # (Typically there is just one model.)
    for model in structure:
        for chain in model:
            # Only consider chains we are interested in.
            if chain.id not in [heavy_chain, light_chain] and chain.id not in antigen_chains:
                continue

            seq_chars = []
            for residue in chain:
                # Only include standard amino acids.
                if residue.get_resname() in PDB.Polypeptide.standard_aa_names:
                    seq_chars.append(seq1(residue.get_resname()))
            chain_seq = "".join(seq_chars)

            if chain.id == heavy_chain:
                sequences["heavy"] = chain_seq
            elif chain.id == light_chain:
                sequences["light"] = chain_seq
            elif chain.id in antigen_chains:
                # If there are multiple occurrences of an antigen chain,
                # you could either overwrite or concatenate (here we concatenate).
                sequences["antigen"][chain.id] += chain_seq

    return sequences

def main():
    metadata_path = "./data/metadata.json"
    if not os.path.exists(metadata_path):
        print(f"Metadata file {metadata_path} not found!")
        return

    with open(metadata_path, "r") as f:
        metadata = json.load(f)

    data_seq_lengths = {}  # To store all chain lengths for each complex.

    # Iterate through each complex in the metadata.
    for complex_id, info in metadata.items():
        pdb_path = info.get("pdb_path")
        if not pdb_path or not os.path.exists(pdb_path):
            print(f"Complex {complex_id}: PDB file '{pdb_path}' not found. Skipping.")
            continue

        heavy_chain = info.get("heavy_chain")
        light_chain = info.get("light_chain")
        antigen_chains = info.get("antigen_chains", [])
        # Build the list of chains in the order expected.
        chains = [heavy_chain, light_chain] + antigen_chains

        sequences = extract_sequences_from_pdb(pdb_path, chains)
        if sequences is None:
            print(f"Complex {complex_id}: Could not extract sequences.")
            continue

        # Get lengths.
        heavy_len = len(sequences.get("heavy", ""))
        light_len = len(sequences.get("light", ""))
        antigen_lengths = {chain: len(seq)
                           for chain, seq in sequences.get("antigen", {}).items()}

        data_seq_lengths[complex_id] = {
            "heavy_length": heavy_len,
            "light_length": light_len,
            "antigen_lengths": antigen_lengths
        }

        print(f"Complex {complex_id}:")
        print(f"  Heavy chain ({heavy_chain}) length: {heavy_len} aa")
        print(f"  Light chain ({light_chain}) length: {light_len} aa")
        for ant_id, length in antigen_lengths.items():
            print(f"  Antigen chain ({ant_id}) length: {length} aa")
        print()

    # Optionally, save the lengths dictionary to a JSON file.
    output_file = "./models/AntiBERTy/sequence_lengths.json"
    with open(output_file, "w") as out_f:
        json.dump(data_seq_lengths, out_f, indent=2)
    print(f"Sequence lengths have been saved to '{output_file}'.")

if __name__ == "__main__":
    main()