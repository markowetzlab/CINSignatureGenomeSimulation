# Simulating genomes with CIN-producing processes

## Quick start

If you just want to produce some copy number profiles with the five mutational processes provided, run the script `scripts/1_Example_CINGenomeSimulation.R` in R or RStudio. Change the number of samples with the variable `N` and perhaps edit the `NAME` variable into something more meaningful and then you are good to go.

## How to edit the mutational processes

This repo allows you to create genomes with copy number aberrations produced by different mutational processes. At the moment, five types of chromosomal instability which are described in the literature will introduce copy number changes into the human genome. These included:
1. chromosome missegregation via mitotic errors (CHR)**[1]**,
2. large-scale state transitions via homologous recombination deficiency (LST)**[2]**,
3. focal amplification via ecDNA circularisation and amplification (ecDNA)**[3]**,
4. early whole-genome duplication via cytokinesis failure (WGD early)**[4]**, and
5. late whole-genome duplication via endoreduplication (WGD late)**[4]**.

Their rate of mutation, segment sizes and copy number changes  are specified in the file `scripts/CINGenomeSimulation_Functions.R` in function `callSignature()`. Please see the section below for implementation details.


## Reproducing Drews et al. (2022)

For reproducing the simulated genomes as used in Drews et al. (2022), please run the script in R or RStudio `scripts/1_Example_CINGenomeSimulation.R`. It will plot an overview of the signatures created by the mutational processes, create a specified number of copy number profiles (default 240) and plot a few samples as examples. Then run the script `scripts/2_Example_Calculate_activities.R` which will extract copy number features, calculate the sum-of-posterior matrix and signature activities. The final script `scripts/3_Evaluate_Activities.R` takes the previously calculated activities and compares them to the originally simulated mutational processes.

If you want to extract novel copy number signatures, please use the code from the SignatureDiscovery repository (TBD: Add link to repo) on the sum-of-posterior matrix resulting from `scripts/2_Example_Calculate_activities.R`.


## Details on the implementation of mutational processes

This simulation was developed for Drews et al. (2022) (TBD: Add link once public) where we used the simulated copy number profiles to study copy number signatures and chromosomal instability. So most of the design choices are explained in the Methods section (TBD: Add link once public). This is a summary of how the mutational processes are implemented.

We simulated a total of 240 samples. We modelled WGD by setting a tetraploid background before inserting other CNAs (early-WGD) or by multiplying all copy number states after all CNAs were placed (late-WGD) in proportions roughly matching what is expected in cancer genomes **[5]**. For any given sample, three of the five mutational processes were active. Half of the samples had one dominant signature and the other half had two (Table 1 and 2). Given that most WGD events are early events, the ratio between early and late WGD samples was 4 to 1 (Nearly=16 and Nlate=4 per partition) **[6]**.

**Table 1:** Distribution of samples exposed to differently active mutational processes
x - Signature is active with low activity (background)
**D** - Signature is active with high activity (dominant)
| **Partition** | **Large-scale transitions (LST)** | **ecDNA(ecDNA)** | **Chromosome missegregation (CHR)** | **Whole-genome duplication (WGD early + late)** | **N** | **Sum** |
|---|---|---|---|---|---|---|
| 1 dominant signature | **D** | x | x |  | 20 | 20 |
|  | x | **D** | x |  | 20 | 40 |
|  | x | x | **D** |  | 20 | 60 |
| 2 dominant signatures | **D** | x | **D** |  | 20 | 80 |
|  | **D** | **D** | x |  | 20 | 100 |
|  | x | **D** | **D** |  | 20 | 120 |
| 2 dominant signatures<br>plus WGD | **D** | x |  | **D** | 16+4 | 140 |
|  | x | **D** |  | **D** | 16+4 | 160 |
|  |  | **D** | x | **D** | 16+4 | 180 |
|  |  | x | **D** | **D** | 16+4 | 200 |
|  | **D** |  | x | **D** | 16+4 | 220 |
|  | x |  | **D** | **D** | 16+4 | 240 |

We used 23 chromosomes (incl. X, excl. Y) of the hg19 reference genome and ignored gender-specific chromosome constellations. The average number of affected chromosomes, number of CNAs produced by a mutational process, as well as the segment size and induced copy number change, were taken from literature where possible (Table 3).

We could not find evidence in the literature for two parameters for our simulation: the changepoint for the LST/HRD process and the average number of CNAs per chromosome for the focal amplifications. Fortunately, our previous work on large numbers of high-grade serous ovarian cancer genomes (120 published **[7]**, ~2000 unpublished) allowed us to qualitatively choose these parameters as means of Poisson distributions.

Firstly, the changepoint for the LST process was not mentioned in the original paper **[2]**. In the ovarian genomes we had observed that most genomes coming from BRCA1/2 germline deficient samples tended to have single copy number changes. Therefore, we assumed a single copy loss / gain changepoint, which has also been qualitatively reported in the literature **[8]**.

Secondly, for the average number of CNAs per chromosome produced by focal amplifications, Deshpande et al. **[9]** stated that on average 5 rearrangements were involved per amplicon using deep WGS. However, in the ovarian genomes, where we used shallow WGS, we qualitatively observed that the number of events per amplicon was lower. This is expected given the lower resolution of the technology. The SNP6 array technology used in this study also has lower resolution than deep WGS, thus we reduced the number to 3.

**Table 2:** Mean affected chromosomes by simulated mutational processes
|  |  | **LST (HRD)** | **ecDNA** | **Chrom. misseg.** | **WGD** |
|---|---|---|---|---|---|
| **Affected chromosomes<br>(mean)** | **Dominant** | 8 **[2]** | 4 **[9]** | 8 **[10, 11]** | 23 |
|  | **Background** | 2 | 2 | 4 | - |

**Table 3:** Details of the four mutational processes
| **Full name** | **ID** | **Average CNAs per chromosome (Poisson)** | **Segment size** | **Changepoint** | **Background ploidy** |
|---|---|---|---|---|---|
| Large-scale transitions (LST) | LST | 7**[2]** | 5 - 30MB<br>**[2]** | 0.8 - 1.4 (loss and gain possible) | 2 |
| ecDNA | ECDNA | 3 | 1MB - 3MB**[12]** | 7-50**[3]** | 2 |
| Early whole-genome duplication | WGD (early) | 1 | Whole chromosome | 2 | 4 |
| Late whole-genome duplication | WGD (late) | 1 | Whole chromosome | x2 of existing segments | 4 |
| Chromosome missegregation | CHR | 1 | Whole chromosome | 0.8 - 1.2 (loss and gain possible) | 2 |

We allowed for multiple mutational processes to be active on the same chromosome by nesting CNAs (smaller CNAs fully overlapped by larger CNAs), but did not allow CNAs to overlap by only one breakpoint.

Homozygous deletions were rare and were only allowed on diploid genomes where an overlap of a HRD-related loss of DNA happened on top of a loss-of-heterozygosity (LOH) produced by a chromosome missegregation. Deletions with a size over 10MB were deemed artificial and were rescued to a subclonal copy number state around 0.5 **[13]**.

Activity of mutational processes were defined by randomly assigning them to individual chromosomes. The number of affected chromosomes was drawn from a Poisson distribution with the mean dependent on the strength (dominant or background) and the type of mutational process (Table 2). We allowed for multiple signatures to target the same chromosome. The number of CNAs per chromosome and mutational process were also drawn by a Poisson (see Table 3 for the means). If the number of CNAs drawn times the maximum allowed size of the segments was larger than the available space on the chromosome (length minus the centromeric regions) then first, the number was corrected to the maximum number of CNAs that can be accommodated if they all had the maximum size. If this did not resolve the issue, then the maximum allowed segment size was set to 2/3rds of the available space on the chromosome and the number of CNAs set to 1. The latter case could happen e.g. for large-scale transitions which can be up to 30MB but the smallest chromosome is around 50MB. If the number of CNAs drawn was zero, it got rescued to 1.

Then we drew the segment sizes and copy number for each CNA and sorted the CNAs by size. By looping over each CNA, we first determined all regions on the chromosome where this CNA could be placed. The rules were:
  - Not overlapping the ends of the chromosomes
  - Not ending or beginning in the centromeric region of the chromosome
  - Not ending or beginning in other CNAs but beginning and ending in the same CNA (nested CNAs) was allowed

This produced a list of allowed regions of the chromosomes. We sampled uniformly from the regions correcting by their size. This way each base pair position was equally likely to get drawn.

For large-scale state transitions (LST) and chromosome missegregation (CHR), the output copy number changes were around 1 (Table 3) and we allowed for gains or losses to happen. For chromosome missegregation the probability of a gain was 50%, and for LST it was 25% since losses were stronger associated with HRD in general **[8]**.

The resulting profiles contained a mixture of multiple mutational processes across the genome.

## Citations

**[1]:** Jallepalli, P. V. & Lengauer, C.; Nat. Rev. Cancer (2001)
**[2]:** Popova, T. et al.; Cancer Res. (2012)
**[3]:** Kim, H. et al.; Nat. Genet. (2020)
**[4]:** Davoli, T. & de Lange, T.; Annu. Rev. Cell Dev. Biol. (2011)
**[5]:** Dentro, S. C. et al.; Cell (2-21)
**[6]:** Bielski, C. M. et al.; Nat. Genet. (2018)
**[7]:** Macintyre, G.; Nat. Genet. (2018)
**[8]:** Abkevich, V. et al.; Br. J. Cancer (2012)
**[9]:** Deshpande, V. et al.; Nat. Commun. (2019)
**[10]:** Weaver, B. A. A. & Cleveland, D. W.; Curr. Opin. Cell Biol. (2006)
**[11]:** Taylor, A. M. et al.; Cancer Cell (2018)
**[12]:** Liao, Z. et al.; Biochim. Biophys. Acta Rev. Cancer (2020)
**[13]:** Cheng, J. et al.; Nat. Commun. (2017)


## Licence
The contents of this repository are copyright (c) 2022, University of Cambridge and Spanish National Cancer Research Centre (CNIO).

The contents of this repository are published and distributed under the GAP Available Source License v1.0 (ASL). 

The contents of this repository are distributed in the hope that it will be useful for non-commercial academic research, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the ASL for more details. 

The methods implemented in the code are the subject of pending patent application GB 2114203.9.

Any commercial use of this code is prohibited.

