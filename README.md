# GeneWeaver: GPU-Accelerated CRISPR Alignment Engine

**GeneWeaver** is a Python-based bioinformatics prototype designed to search DNA sequences for potential CRISPR-Cas9 target matches. 

## Project Overview
The core idea behind GeneWeaver is to computationally screen and rank potential DNA target sites. It takes a user-provided guide RNA sequence and a reference DNA sequence (or FASTA file), locates suitable PAM (Protospacer Adjacent Motif) sites, extracts candidate target sequences, and compares them against the guide to evaluate their similarity.

This tool is built to accelerate the candidate-search stage of CRISPR alignment by utilizing GPU parallel processing. This allows for the rapid screening of large genetic datasets before executing more expensive detailed alignments.

## Core Objectives
- **Target Identification:** Scan both forward and reverse-complement DNA strands to identify NGG PAM-compatible regions.
- **Sequence Alignment:** Compare extracted 20-nucleotide candidates with the guide sequence to calculate mismatch counts and alignment scores.
- **Performance Benchmarking:** Provide a direct, measurable comparison between standard CPU execution and GPU-accelerated processing for genome alignment tasks.

## Project Scope
This initial prototype is strictly focused on the **SpCas9 system** (using the NGG PAM pattern) and 20-nucleotide guide comparisons. 

**Note:** GeneWeaver is designed to act as a computational search metric and candidate evaluator, rather than a biological cleavage-risk predictor. Advanced features like bulge-aware alignment or support for alternative Cas enzymes are outside the current baseline scope and may be added in future iterations.
