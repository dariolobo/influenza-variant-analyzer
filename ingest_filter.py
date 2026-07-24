import os
from Bio import SeqIO

def filter_neuraminidase_sequences(input_fasta, output_fasta):
    """
    Reads a FASTA file. filters for Neuraminidade (NA) sequences,
    and writes the valid ones on a new FASTA file.
    """
    print(f"--- Starting new digestion of: {input_fasta} ---")
    
    total_records = 0
    valid_sequences = []

    if not os.path.exists(input_fasta):
        print(f"[❌ Error] Input file '{input_fasta}' not found.")
        return False
    
    # 1. Read and parse the input FASTA file
    with open(input_fasta, "r") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            total_records += 1
            description = record.description.lower()

            # Validation: Match 'neuraminidase' or 'na' as standalone words, preventing false positives from longer words containing 'na'
            if "neuraminidase" in description or "na" in description.split():
                valid_sequences.append(record)
            else:
                print(f"[⚠️ Skipping] Non-NA or low quality record: {record.id}")

    # 2. Write the cleaned/filtered sequences to the output file
    if valid_sequences:
        with open(output_fasta, "w") as outfile:
            SeqIO.write(valid_sequences, outfile, "fasta")
        print(f"\n--- Print Summary ---")
        print(f"Total sequences scanned: {total_records}")
        print(f"Valid NA sequences saved: {len(valid_sequences)} -> Written to '{output_fasta}'")
        return True
    else:
        print("\n[❌ Error] No valid Neuraminidase sequence were found after filtering.")
        return False
    
# This block allows you to run this module directly from the terminal for testing
if __name__ == "__main__":
    # We define standard test filenames
    input_file = "raw_influenza_data.fasta"
    output_file = "filtered_na_sequences.fasta"

    print("Running Module 1 locally...")
    filter_neuraminidase_sequences(input_file, output_file)