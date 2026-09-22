# Rust implementation of CpG-likelihoods
Single-molecule sequencing is readily adapted for the analysis of epigenetic DNA modifications such as 5-methylcytosine (5mC) and 5-hydroxymethylcytosine (5hmC). One major advantage of Oxford Nanopore (Nanopore) and Pacific Biosciences (PacBio) sequencing platforms is that they capture both the primary DNA sequence and epigenetic modification states. This enables epigenetic analysis of sample-specific CpG loci that are not represented in the reference genome which typically account for  50% of heterozygous CpG loci in a sample, and a number of homozygous CpG loci depending on the evolutionary distance between the sample and the reference genome. Because these sample-specific loci are important for accurate modification analysis, this project provides a workflow for identifying them.

## Installation
{add installation guide}

## Example usage:
{add example}

```bash
modkit pileup --cpg-likelihood --interactive
```

# Workflow overview
The workflow consists of three modules:

## Estimate error rate
This step estimates the combined error rate for sequencing and mapping errors at CpG loci. The input is a `bedMethyl` file containing read-based pileup information for all candidate CpG sites. 

For each locus, the pileup has recorded genotypic information as:

- The number of mapped CpG-sites (C<sub>i</sub>) : N<sub>valid_cov</sub>
- The number of mapped non-CpGs (D<sub>i</sub>) :
  - Number of deletions : N<sub>delete</sub> 
  - Number of DN : N<sub>diff</sub>
  - Number of CH : N<sub>nocall</sub> 

The function first identifies a subset of loci that are likely to represent homozygous CpG sites. By default the loci are retained if: 

$$
10 \leq C_i \leq 100 \qquad \text{and} \qquad \frac{C_i}{C_i + D_i} > 0.8,
$$

The corresponding thresholds can be adjusted using the command-line parameters `--min-cov`, `--max-cov`, and `--min-ratio`.

Assuming these loci are true homozygous CpG sites, the global CpG-to-non-CpG miscall rate is estimated as:

$$
\hat{\varepsilon} = \frac{\sum_{i} D_i} {\sum_{i} (C_i + D_i)}.
$$

This estimate is treated as the sample-wide CpG error rate and is used in subsequent likelihood calculations. The approach assumes that sequencing and mapping errors are approximately symmetric across CpG and non-CpG states throughout the genome.

**Example**

| chrom | start | end | ... | N<sub>valid_cov | N<sub>delete | N<sub>diff | N<sub>nocall | Retained |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr1 | 1000 | 1001 | ... | 30 | 0 | 0 | 2 | ✓ |
| chr1 | 1001 | 1002 | ... | 28 | 0 | 0 | 0 | ✓ |
| chr1 | 1011 | 1012 | ... | 15 | 0 | 15 | 1 | ✗ (`min-ratio`) |
| chr1 | 1012 | 1013 | ... | 13 | 0 | 15 | 1 | ✗ (`min-ratio`) |
| chr1 | 1021 | 1022 | ... | 1 | 0 | 28 | 1 | ✗ (`min-cov`,`min-ratio`) |
| chr1 | 1033 | 1034 | ... | 30 | 6 | 1 | 1 | ✗ (`min-ratio`) |

$$
\hat{\varepsilon} = \frac{2} {(58+2)} = 0.03333
$$

## Likelihoods
This step models the likelihood of observing the read data given an underlying CpG state. Under the assumption that all loci are biallelic and diploid the three possible CpG states (z) are:

- Homozygues CpG (CpG<sub>hom</sub>): Both alleles at the locus are CpG
- Heterozygues CpG (CpG<sub>het</sub>): One allele at the locus is CpG, the other in non-CpG - non-CpG (nonCpG<sub>het</sub>): Neither alleles at the locus are CpGs

The likelihood of the data given the genotype and error is modeled assuming reads and errors are independent. Specifically the likelihoods for the different CpG states at each locus are computed as: 

$$
L({\mathrm{CpG_hom}})=
(1-\hat{\varepsilon})^{\ C_i} \ 
\hat{\varepsilon}\ ^{D_i}
$$

$$
L({\mathrm{CpG_het}})=
0.5^{\ C_i+D_i}
$$

$$
L({\mathrm{non-CpG}})=
\hat{\varepsilon}^{\ C_i} \ 
(1-\hat{\varepsilon})^{\ D_i}
$$

These values are transformed to raw Phred-scaled likelihoods for the CpG state L(z|Data) by -10*log(L(z|Data)) and normalized such that the lowest PL for each locus is 0 and the rest are scaled relative to that. 

**Example**

The first locus in the example above (`chr1:1000-1001`) had 30 read supporting CpG and 2 reads supporting non-CpG and a estimated error rate of 0.03333:

$$
C_i = 30
\qquad
D_i = 2
\qquad
\hat{\varepsilon} = 0.03333
$$

The raw likelihoods are:

$$
L(\mathrm{CpG_{hom}})=
(1-0.03333)^{30}
\times
0.03333^{2} =
4.00 \times 10^{-4}
$$

$$
L(\mathrm{CpG_{het}})=
0.5^{30+2} =
2.33 \times 10^{-10}
$$

$$
L(\mathrm{nonCpG})=
0.03333^{30}
\times
(1-0.03333)^{2} = 
6.46 \times 10^{-45}
$$

The raw Phred-scaled likelihoods are:

$$
PL(\mathrm{CpG_{hom}})=
-10\log_{10}(4.00 \times 10^{-4})=
33.98
$$

$$
PL(\mathrm{CpG_{het}})=
-10\log_{10}(2.33 \times 10^{-10})=
96.33
$$

$$
PL(\mathrm{nonCpG})=
-10\log_{10}(6.46 \times 10^{-45})=
441.90
$$

The normalized PL where the minimum PL is then subtracted from all values are:

$$
PL(\mathrm{CpG_{hom}})=0
$$
$$
PL(\mathrm{CpG_{het}})=62
$$
$$
PL(\mathrm{CpG_{het}})=408
$$

Indicating strong support for a homozygous CpG at (`chr1:1000-1001`). 

The example table with normalized PL added:

| chrom | start | end | ... | N_valid_cov | N_delete | N_diff | N_nocall | PL(CpG_hom) | PL(CpG_het) | PL(non-CpG) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr1 | 1000 | 1001 | ... | 30 | 0 | 0 | 2 | 0 | 62 | 408 |
| chr1 | 1001 | 1002 | ... | 28 | 0 | 0 | 0 | 0 | 77 | 412 |
| chr1 | 1011 | 1012 | ... | 15 | 0 | 15 | 1 | 193 | 0 | 47 |
| chr1 | 1012 | 1013 | ... | 13 | 0 | 15 | 1 | 181 | 0 | 34 |
| chr1 | 1021 | 1022 | ... | 1 | 0 | 28 | 1 | 409 | 71 | 0 |
| chr1 | 1033 | 1034 | ... | 30 | 6 | 1 | 1 | 8 | 0 | 330 |

### Summarize/filter likelihoods
This step summarizes and filter across the inferred CpG states in a `bedMethyl` file after likelihood calculations.

For each locus, the function: 
1. Identifies the most likely CpG state based on the normalized PL values.
2. Identifies the second most likely CpG state.
3. Calculates a quality score (GQ) for the CpG state as the difference between the lowest (0) and second-lowest PL value.

It is often hard to distinguish homozygous and heterozygous CpGs leading to a high number of unclassified loci with CpG<sub>hom</sub> and CpG<sub>het</sub> as the two most likely CpG states, but low QC. To keep these loci for downstream analysis CpG<sub>hom</sub> and CpG<sub>het</sub> can be collapsed to combined CpG state (CpG<sub>hom/het</sub>) with the corresponding GQ describing the certainty that the locus is not a non-CpG. This setting classifies more loci than CpG<sub>hom</sub> and CpG<sub>het</sub>.  

**Example**

| chrom | start | end | ... | N_valid_cov | N_delete | N_diff | N_nocall | PL(CpG_hom) | PL(CpG_het) | PL(non-CpG) | state | GQ
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr1 | 1000 | 1001 | ... | 30 | 0 | 0 | 2 | 0 | 62 | 408 | CpG_hom | 62 |
| chr1 | 1001 | 1002 | ... | 28 | 0 | 0 | 0 | 0 | 77 | 412 | CpG_hom | 77 |
| chr1 | 1011 | 1012 | ... | 15 | 0 | 15 | 1 | 193 | 0 | 47 | CpG_het | 47 |
| chr1 | 1012 | 1013 | ... | 13 | 0 | 15 | 1 | 181 | 0 | 34 | CpG_het | 34 |
| chr1 | 1021 | 1022 | ... | 1 | 0 | 28 | 1 | 409 | 71 | 0 | non-CpG | 71 |
| chr1 | 1033 | 1034 | ... | 30 | 6 | 1 | 1 | 8 | 0 | 330 | CpG_het | 8 |

Given a minimum GQ threshold the summary reports: 
- Total number of loci in the bedMethyl file.
Classification with separate homozygous and heterozygous CpG states:
- Number and percent of loci classified as CpG_hom.
- Number and percent of loci classified as CpG_het.
- Number and percent of loci classified as non-CpG.
- Number and percent of unclassified loci.   
Classification with collapsed homozygous and heterozygous CpG states:
- Number and percent of loci classified as CpG_hom/het.
- Number and percent of loci classified as non-CpG.
- Number and percent of unclassified loci. 

To filter a bedMethyl file after likelihood calculations we can use the filter likelihood function which follow the same logic as the summary function

**Example**

| chrom | start | end | ... | N_valid_cov | N_delete | N_diff | N_nocall | PL(CpG_hom) | PL(CpG_het) | PL(non-CpG) | state | GQ
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr1 | 1000 | 1001 | ... | 30 | 0 | 0 | 2 | 0 | 62 | 408 | CpG_hom | 62 |
| chr1 | 1001 | 1002 | ... | 28 | 0 | 0 | 0 | 0 | 77 | 412 | CpG_hom | 77 |
| chr1 | 1011 | 1012 | ... | 15 | 0 | 15 | 1 | 193 | 0 | 47 | CpG_het | 47 |
| chr1 | 1012 | 1013 | ... | 13 | 0 | 15 | 1 | 181 | 0 | 34 | CpG_het | 34 |
| chr1 | 1021 | 1022 | ... | 1 | 0 | 28 | 1 | 409 | 71 | 0 | non-CpG | 71 |
| chr1 | 1033 | 1034 | ... | 30 | 6 | 1 | 1 | 8 | 0 | 330 | CpG_het | 8 |

filter likelihood --CpG_hom/het --GQ 20 

| chrom | start | end | ... | N_valid_cov | N_delete | N_diff | N_nocall
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr1 | 1000 | 1001 | ... | 30 | 0 | 0 | 2 |
| chr1 | 1001 | 1002 | ... | 28 | 0 | 0 | 0 |
| chr1 | 1011 | 1012 | ... | 15 | 0 | 15 | 1 | 
| chr1 | 1012 | 1013 | ... | 13 | 0 | 15 | 1 | 
| chr1 | 1033 | 1034 | ... | 30 | 6 | 1 | 1 |

The PL columns can be retained using --keep-pl
