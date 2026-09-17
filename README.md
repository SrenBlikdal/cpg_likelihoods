# Rust implementation of CpG-likelihoods

Single molecule sequencing is redily addapted for the analysis of epigenetic modification of DNA methylation (5mC) and hydroxymethylation. 
One advantage of single-molecule sequencing platforms from Oxford Nanopore (Nanopore) and Pacific bio.. (PacBio) is the readout of both the primary sequence
and epigenetic modificatation in modbam format. This enables extending analysis to sample specific CpG loci, not represented in the reference genome. 
These non-refererence CpG loci typically represent ~50 % of the heterozygous CpG loci in a sample and a number of homozygous CpG loci depending on the 
evolutionary distance between a sample and the reference genome. As these sample-specific loci are important for accurate modification analysis we here 
provide a plugin for identifying sample-specific CpG loci in a sample. 

Example usage modkit pileup --cpg-likelihood --interactive




The workflow has three modules:

Estimate error rate -

CpG-likelihoods -

Summarize likelihoods -



