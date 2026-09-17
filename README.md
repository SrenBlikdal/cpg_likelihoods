# Rust implementation of CpG-likelihoods

Single-molecule sequencing is readily adapted for the analysis of epigenetic DNA modifications such as 5-methylcytosine (5mC) and 5-hydroxymethylcytosine (5hmC). One major advantage of Oxford Nanopore (Nanopore) and Pacific Biosciences (PacBio) sequencing platforms is that they can capture both the primary DNA sequence and epigenetic modification state in modified BAM (modBAM) format.

This enables analysis at sample-specific CpG loci that are not represented in the reference genome. These non-reference CpG sites often account for approximately 50% of heterozygous CpG loci in a sample, and also include a number of homozygous CpG loci depending on the evolutionary distance between the sample and the reference genome.

Because these sample-specific loci can be important for accurate modification analysis, this project provides a workflow for identifying and analyzing them.

## Installation
{add installation guide}

## Example usage:
{add example}

```bash
modkit pileup --cpg-likelihood --interactive
```

## Workflow overview
The workflow consists of three modules:

1. Estimate error rate
This step estimates the sequencing/mapping error rate at putative CpG sites in the bedmethyl file. This error rate is used to model which sites represent homozygous, heterozygous, and non-CpG sites.

{insert formula}

{insert example}

2. CpG-likelihoods
This module computes CpG likelihoods for all sites represented in the bedmethyl file, including CpG loci that are not present in the reference genome.

{insert formula}

{insert example}

3. Summarize likelihoods
This step summarizes the computed likelihoods for downstream interpretation and comparison across loci or samples.

{insert formula}

{insert example}
