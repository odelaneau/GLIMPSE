---
layout: default
title: phase
nav_order: 3
parent: Documentation
---
# phase
{: .no_toc .text-center }

{: .highlight }
Website under construction. A complete release of GLIMPSE2 will be available by the 7th of December 2022.

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

### Description
Tool for imputation and phasing.

### Usage
Simple run

<div class="code-example" markdown="1">
```bash
GLIMPSE2_phase --bam-list bams_1.0x.txt --reference binary_reference_panel_chr20_7702567_12266861.bin --output imputed_glimpse2_rp140k_1.0x_chr20_7702567_12266861.bcf --threads 4
```
</div>

---

### Command line options

#### Basic options

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-\-help             | NA      | NA       | Produces help message |
| \-\-seed             | INT     | 15052011 | Seed of the random number generator  |
| \-T \[ \-\-threads \]| INT     | 1        | Number of threads|


#### Input parameters

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-\-bam-file         | FILE    | NA       | Input BAM/CRAM file containing low-coverage sequencing reads. Only one of the following options can be declared: \-\-input-gl, \-\-bam-file, \-\-bam-list.|
| \-\-bam-list         | FILE    | NA       | List (.txt file) of input BAM/CRAM files containing BAM/CRAM file containing low-coverage sequencing reads. One file per line. A second column (space separated) can be used to specify the sample name, otherwise the name of the file is used. Only one of the following options can be declared: \-\-input-gl, \-\-bam-file, \-\-bam-list.|
| \-\-input-gl         | FILE    | NA       | VCF/BCF file containing the genotype likelihoods. Only one of the following options can be declared: \-\-input-gl, \-\-bam-file, \-\-bam-list.|
| \-R \[\-\-reference\]| FILE    | NA       | Haplotype reference in VCF/BCF or binary format |
| \-M \[ \-\-map \]    | FILE    | NA       | Genetic map |
| \-\-input-region     | STRING  | NA       | Imputation region with buffers |
| \-\-output-region    | STRING  | NA       | Imputation region without buffers |
| \-\-sparse-maf       | FLOAT   | 0.001   | **Expert setting.** Rare variant threshold |
| \-\-samples-file     | STRING  | NA       | File with sample names and ploidy information. One sample per line with a mandatory second column indicating ploidy (1 or 2). Sample names that are not present are assumed to have ploidy 2 (diploids). If the parameter is omitted, all samples are assumed to be diploid. GLIMPSE does NOT handle the use of sex (M/F) instead of ploidy.|
| \-\-ind-name         | STRING  | NA       | Only used together with \-\-bam-file. Name of the sample to be processed. If not specified the prefix of the BAM/CRAM file (\-\-bam-file) is used.|
| \-\-keep-monomorphic-ref-sites | NA       | NA       | **Expert setting.** Keeps monomorphic markers in the reference panel (removed by default) |
| \-\-impute-reference-only-variants | NA       | NA       | Only used together with \-\-input-gl. Allows imputation at variants only present in the reference panel (no GL called at these positions). The use of this option is intended only to allow imputation at sporadic missing variants. If the number of missing variants is non-sporadic, please re-run the genotype likelihood computation at all reference variants and avoid using this option, since data from the reads should be  used. A warning is thrown if reference-only variants are found.|
| \-\-input-field-gl   | NA       | NA       | Only used together with \-\-input-gl. Use FORMAT/GL field instead of FORMAT/PL to read genotyope likelihoods |

#### Model parameters

| Option name 	      | Argument|  Default  | Description |
|:--------------------|:--------|:----------|:-------------------------------------|
| \-\-burn-in         | INT     | 5         | **Expert setting.** Number of burn-in iterations of the Gibbs sampler
| \-\-main            | INT     | 15        | **Expert setting.** Number of main iterations of the Gibbs sampler
| \-\-ne              | INT     | 100000    | **Expert setting.** Effective diploid population size modelling recombination frequency
| \-\-min-gl          | FLOAT   | 1e-10     | **Expert setting.** Minimim haploid likelihood
| \-\-err-imp         | FLOAT   | 1e-12     | **Expert setting.** Imputation HMM error rate
| \-\-err-phase       | FLOAT   | 1e-4      | **Expert setting.** Phasing HMM error rate

#### Selection parameters

| Option name 	      | Argument|  Default  | Description |
|:--------------------|:--------|:----------|:-------------------------------------|
| \-\-pbwt-depth      | INT     | 12        | **Expert setting.** Number of neighbors in the sparse PBWT selection step (positive number).
| \-\-pbwt-modulo-cm  | FLOAT   | 5         | **Expert setting.** Frequency of PBWT selection in cM (positive number). This parameter is automatically adjusted in case of small imputation regions.
| \-\-Kinit           | INT     | 1000      | **Expert setting.** Number of states used for initialization (positive number). Can be set to zero only when --state-list is set, to skip the selection for the initialization step.
| \-\-Kpbwt           | INT     | 2000      | **Expert setting.** Maximum number of states selected from the sparse PBWT (positive number). Can be set to zero only when --state-list is set, to skip the selection for during the Gibbs iterations.
| \-\-state-list      | FILE    | 5         | **Expert setting.** List (.txt file) of haplotypes always present in the conditioning set, independent from state selection. Not affected by other selection parameters. Each row is a target haplotype (two lines per sample in case of diploid individuals) each column is a space separated list of reference haplotypes (in numberical order 0-(2N-1) ). Useful when prior knowledge of relatedness between the reference and target panel is known a priori.

#### BAM/CRAM options and filters

| Option name 	      | Argument|  Default  | Description |
|:--------------------|:--------|:----------|:-------------------------------------|
| \-\-call-model      | STRING  | standard  | Model to use to call the data. Only the standard model is available at the moment.
| \-\-call-indels     | NA      | NA        | **Expert setting.** Use the calling model to produce genotype likelihoods at indels. However the likelihoods from low-coverage data can be miscalibrated, therefore by default GLIMPSE does only imputation into the haplotype scaffold (assuming flat likelihoods)
| \-F \[ \-\-fasta \] | FILE    | NA        | Faidx-indexed reference sequence file in the appropriate genome build. Necessary for CRAM files
| \-\-mapq            | INT     | 10        | Minimum mapping quality for a read to be considered
| \-\-baseq           | INT     | 10        | Minimum phred-scalde based quality to be considered
| \-\-max-depth       | INT     | 40        | **Expert setting.** Max read depth at a site. If the number of reads exceeds this number, the reads at the sites are downsampled (e.g. to avoid artifactual coverage increases).
| \-\-keep-orphan-reads | NA      | NA        | **Expert setting.** Keep paired end reads where one of mates is unmapped
| \-\-ignore-orientation | NA      | NA        | **Expert setting.** Ignore the orientation of mate pairs
| \-\-check-proper-pairing | NA      | NA        | **Expert setting.** Discard reads that are not properly paired
| \-\-keep-failed-qc  | NA      | NA        | **Expert setting.** Keep reads that fail sequencing QC (as indicated by the sequencer).
| \-\-keep-duplicates | NA      | NA        | **Expert setting.** Keep duplicate sequencing reads in the process
| \-\-illumina13+     | NA      | NA        | **Expert setting.** Use illimina 1.3 encoding for the base quality (for older sequencing machines)

#### Output parameters

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-O \[\-\-output \]  | FILE    | NA       | Phased and imputed haplotypes in VCF/BCF/BGEN format |
| \-\-contigs-fai      | FILE    | NA       | If specified, header contig names and their lengths are copied from the  provided fasta index file (.fai). This allows to create imputed whole-genome files as contigs are present and can be easily merged by bcftools |
| \-\-bgen-bits        | INT     | 8        | **Expert setting.** Only used together when the output is in BGEN file format. Specifies the number of bits to be used for the encoding probabilities of the output BGEN file. If the output is in the .vcf\[.gz\]/.bcf format, this value is ignored. Accepted values: 1-32. |
| \-\-bgen-compr       | STRING  | zstd     | **Expert setting.** Only used together when the output is in BGEN file format. Specifies the compression of the output BGEN file. If the output is in the .vcf\[.gz\]/.bcf format, this value is ignored.  Accepted values: \[no,zlib,zstd\] |
| \-\-log              | FILE    | NA       | Log file  |


