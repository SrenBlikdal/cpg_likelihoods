# Rust implementation of CpG-likelihoods

Single molecule sequencing is redily addapted for the analysis of epigenetic modification of DNA methylation (5mC) and hydroxymethylation. 
One advantage of single-molecule sequencing platforms from Oxford Nanopore (Nanopore) and Pacific bio.. (PacBio) is the readout of both the primary sequence
and epigenetic modificatation in modbam format. This enables extending analysis to sample specific CpG loci, not represented in the reference genome. 
These non-refererence CpG loci typically represent ~50 % of the heterozygous CpG loci in a sample and a number of homozygous CpG loci depending on the 
evolutionary distance between a sample and the reference genome. As these sample-specific loci are important for accurate modification analysis we here 
provide a plugin for identifying sample-specific CpG loci in a sample. 

Example usage modkit pileup --cpg-likelihood --interactive

Use this error rate or manually type errorrate. 

Calculating likelihoods for file 
Save output as [inputfile].bed 

# How it works

## Estimate error rate
First part of the program is to estimate the sample-specific CpG error-rate directly from the bedMethyl file after modkit pileup.  
While the errorrate depend on sequencing context, mapping, and we use a simple approach of... 

{Insert Formula}

Where X,Y and Z represent. 

This module can be skipped by manually assigning a error rate using --error_c2n

## CpG-likelihoods 
The likelihood that a given site is...

{insert formula here}

{Work through an example}

IGV screenshot - 

Loci xx has

The PL scores are added to col19 col20 col21  

## Summarize likelihoods
A bedmethyl file with PL fields in columns 19-21 can be summarize to check the number of homozygous given the type of CpG state and a GQ threshold..








