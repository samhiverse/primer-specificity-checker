# Primer Specificity and Cross-Reactivity Checker

## Project Overview
This project is a Python-based tool designed to check primer specificity and identify potential PCR cross-reactivity. It uses NCBI databases and BLAST searches to examine whether forward and reverse primers can bind to other DNA sequences.

## Aim
To identify potential unintended primer-binding sites and predict possible off-target PCR products.

## Technologies Used
- Python
- NCBI E-utilities
- NCBI BLAST
- DNA sequence analysis

## Project Workflow
1. Input forward and reverse primer sequences.
2. Check primer compatibility with the intended target sequence.
3. Retrieve nucleotide sequences from NCBI.
4. Perform BLAST searches for both primers.
5. Analyse common sequence hits and primer orientations.
6. Verify candidate binding sites and estimate potential PCR product sizes.
7. Generate a cross-reactivity analysis report.

## Test Case
- Target gene: GAPDH
- Organism: Homo sapiens
- Accession: NM_002046.7
- Expected amplicon: Approximately 63 bp

## Preliminary Results

The analysis identified three candidate sequence records:

| Accession | Potential product |
|---|---:|
| PX660077.1 | 63 bp |
| NG_125030.1 | 303 bp |
| OQ420599.1 | 63 bp |

These are computational predictions and do not confirm experimental amplification.

## Files
- `CrossReactivityCheck.py` – Basic primer checking.
- `target_sequence.txt` – Target DNA sequence.
- `ncbi_connection_test.py` – NCBI sequence retrieval.
- `ncbi_blast_test.py` – BLAST search for primer sequences.
- `ncbi_blast_results.txt` – BLAST output.
- `cross_reactivity_analysis.py` – Analysis of common BLAST hits.
- `cross_reactivity_report.txt` – Preliminary analysis report.
- `verify_candidate_sequences.py` – Verification of candidate binding sites.

## How to Run
Install Python 3, then run the required scripts from the project directory:

```bash
python CrossReactivityCheck.py
python ncbi_connection_test.py
python ncbi_blast_test.py
python cross_reactivity_analysis.py
python verify_candidate_sequences.py
```

An internet connection is required for the NCBI-based scripts.

## Limitations
- The analysis is limited to the BLAST hits returned by the search.
- Potential PCR products are computational predictions.
- Experimental PCR validation is not included.

## Conclusion
The project demonstrates a basic computational approach to checking primer specificity and identifying potential cross-reactivity using Python and NCBI BLAST.
