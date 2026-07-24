# Influenza A Neuraminidase Variant Analyzer (H275Y Resistance Detection)

An automated Python bioinformatics pipeline designed to ingest, filter, and analyze Influenza A sequence data. The primary objective is to align patient sequences against a reference protein to accurately detect antiviral resistance mutations, specifically **H275Y** in the Neuraminidase (NA) protein.

---

## 🧬 Project Overview & Biological Context

Antiviral resistance in Influenza A (such as resistance to oseltamivir / Tamiflu) is frequently conferred by single nucleotide variations leading to amino acid substitutions in key viral proteins. 

This project automates the clinical variant screening workflow by:
1. Fetching or utilizing standard reference protein sequences (NCBI accession: `ACP44189.1`).
2. Filtering raw FASTA sequence files based on quality and completeness.
3. Performing pairwise protein sequence alignments using Biopython's `PairwiseAligner`.
4. Mapping target clinical positions (position 275) to output structured resistance reports in CSV format.

---

## 🏗️ Pipeline Architecture

The workflow is divided into two main modular scripts:

```text
raw_influenza_data.fasta
          │
          ▼
   ingest_filter.py  ──► Generates quality-filtered FASTA ('filtered_na_sequences.fasta')
          │
          ▼
 variant_analyzer.py ──► Fetches Reference + Performs Alignment
          │
          ▼
resistance_report.csv
```

### Module Breakdown:
* **`ingest_filter.py`**: Reads raw genomic/protein data, filters out non-NA or low-quality/incomplete sequences, and outputs a clean input dataset.
* **`variant_analyzer.py`**: Loads reference sequences via NCBI Entrez or local cache, executes pairwise alignment against filtered patient records, inspects clinical amino acid positions, and exports results.

---

## 🚀 Getting Started

### Prerequisites
* Python 3.8+
* `biopython` library

### Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/your-username/influenza-variant-analyzer.git
   cd influenza-variant-analyzer
   ```

2. **Set up a virtual environment (Recommended):**
   ```bash
   python3 -m venv bio_env
   source bio_env/bin/activate
   ```

3. **Install dependencies:**
   ```bash
   pip install biopython
   ```

---

## 💻 Usage

Run the pipeline sequentially in your terminal:

1. **Step 1: Filter raw sequences**
   ```bash
   python3 ingest_filter.py
   ```

2. **Step 2: Run variant analysis**
   ```bash
   python3 variant_analyzer.py
   ```

---

## 📊 Sample Output

Upon completion, `variant_analyzer.py` generates `resistance_report.csv` formatted as follows:

| sample_id | ref_position | reference_aa | patient_aa | status |
| :--- | :---: | :---: | :---: | :--- |
| `EPI_ISL_123456` | 275 | E | - | Indeterminate (Gap/Deletion) |
| `EPI_ISL_789101` | 275 | E | E | Wild Type / Standard Variant |

---

## 🛠️ Built With

* **[Python 3](https://www.python.org/)** - Core programming language.
* **[Biopython](https://biopython.org/)** - For FASTA parsing, Entrez fetching, and pairwise sequence alignment (`Bio.Align.PairwiseAligner`).