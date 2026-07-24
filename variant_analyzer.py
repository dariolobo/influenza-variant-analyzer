import os
import csv
from Bio import Entrez
from Bio import SeqIO
from Bio import Align

# --- NCBI ACCESS CONFIGURATION ---
# The email is loaded from enviroment variables for security and privacy
Entrez.email = "302364874+dariolobo@users.noreply.github.com"

def fetch_reference_protein(accession_id="ACP44189"):
    """
    Downloads the reference protein sequence from the NCBI Protein database
    using Entrez and caches it locally to prevent redundant network requests.
    """
    filename = f"reference_{accession_id}.fasta"

    # Check if the file already exists locally to save bandwith
    if os.path.exists(filename):
        print(f"[*] Reference sequence '{filename}' already exist locally. Loading...")
        return SeqIO.read(filename, "fasta")
    
    print(f"[*] Downloading reference protein '{accession_id}' from NCBI...")
    try:
        # Fetch the sequence from the NCBI Protein database in FASTA format
        with Entrez.efetch(db="protein", id=accession_id, rettype="fasta", retmode="text") as handle:
            record = SeqIO.read(handle, "fasta")

        # Cache the record locally for future executions
        with open(filename, "w") as outfile:
            SeqIO.write(record, outfile, "fasta")

        print(f"[✅ Success] Reference protein saved as '{filename}'")
        return record
    except Exception as e:
        print(f"[❌ Error] Failed to connect or download from NCBI: {e}")
        return None

def analyze_patient_sequences(reference_record, patient_fasta_path="filtered_na_sequences.fasta"):
    """
    Performs local pairwise alignment between the reference protein and patient sequences
    to accuraately detect the H275Y resistance mutation in the neuraminidase protein.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode= 'local' #Smith-Waterman algorithm for local alignment

    # Standard scoring parameters for protein alignments
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1

    print(f"\n[*] Starting variant analyzer against reference: {reference_record.id}")
    print(f"[*] Reading patient sequences from '{patient_fasta_path}'...\n")

    # Initialize list to store results for CSV export
    results = []

    try:
        for patient_record in SeqIO.parse(patient_fasta_path, "fasta"):
            try:
                alignments = aligner.align(reference_record.seq, patient_record.seq)
                best_alignment = alignments[0]
            except (IndexError, OverflowError, Exception):
                print(f"[⚠️ Warning] No alignment could be generated for sample: {patient_record.id}")
                continue

            # Clinical position 275 maps to 0-based index 274 index in Python
            target_index = 274

            patient_amino_acid = "?"

            # Alternative mapping logic using .aligned property
            try:
                ref_blocks, pat_blocks = best_alignment.aligned
                for ref_block, pat_block in zip(ref_blocks, pat_blocks):
                    ref_start, ref_end = ref_block
                    if ref_start <= target_index < ref_end:
                        offset = target_index - ref_start
                        pat_start, pat_end = pat_block
                        patient_coordinate = pat_start + offset
                        patient_amino_acid = patient_record.seq[patient_coordinate]
                        break
                else:
                    patient_amino_acid = "-"
            except Exception:
                patient_amino_acid = "?"

            ref_amino_acid = reference_record.seq[target_index]

            # Pharmacological classification of the variant
            if patient_amino_acid == "Y":
                status = "RESISTANT (H275Y Mutation Detected)"
                console_icon = "🚨"
            elif patient_amino_acid == "H":
                status = "Sensitive (Wild-Type Histidine)"
                console_icon = "✅"
            elif patient_amino_acid == "-":
                status = "Indeterminate (Deletion/Gap at position 275)"
                console_icon = "⚠️"
            else:
                status = f"Unknown Variant ({ref_amino_acid}275{patient_amino_acid})"
                console_icon = "❓"
            # Print formatted log to terminal with icon
            print(f"[Result] {patient_record.id[:25]}... | Pos 275: {ref_amino_acid} -> {patient_amino_acid} | {console_icon} {status}")

            # Collect patient metrics for CSV report
            results.append({
                "sample_id": patient_record.id,
                "ref_position": 275,
                "reference_aa": ref_amino_acid,
                "patient_aa": patient_amino_acid,
                "status": status
            })

    # Write extracted metrics to CSV file
        output_csv = "resistance_report.csv"
        if results:
            with open(output_csv, "w", newline="", encoding="utf-8") as csvfile:
                fieldnames = ["sample_id", "ref_position", "reference_aa", "patient_aa", "status"]
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(results)
            print(f"\n[✔] Report successfully saved to '{output_csv}'")

    except FileNotFoundError:
        print(f"[❌ Error] Target file '{patient_fasta_path} was not found. Please verify the path.")
    except Exception as e:
        print(f"[❌ Error] An unexpected error ocurred: {e}")

if __name__ == "__main__":
    # Execute an initial download test
    ref_seq = fetch_reference_protein()
    if ref_seq:
        print(f"Loaded reference: {ref_seq.id} | Length: {len(ref_seq.seq)} amino acids")

        # Run the secondary analysis pipeline against your filtered samples
        analyze_patient_sequences(ref_seq, "filtered_na_sequences.fasta")