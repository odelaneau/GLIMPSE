---
layout: default
title: concordance
nav_order: 5
parent: Documentation
---
# concordance
{: .no_toc .text-center }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

### Description
Program to compute the genotyping error rate at the sample or marker level.

### Usage
Simple run

<div class="code-example" markdown="1">
```bash
echo "chr20 ukb_white_british_sites.bcf validation_c20_snps.bcf imputed_ligated_c20_rp140k_1.0x.bcf" > concordance.txt

GLIMPSE2_concordance --gt-val --ac-bins 1 5 10 20 50 100 200 500 1000 2000 5000 10000  20000 50000 100000 140119 --threads 2 --input concordance.txt --output concordance_c20_rp140
```
</div>

#### Sample aliasing

By default, samples are matched between the imputed VCF and the validation VCF by exact ID. When the same biological sample appears under different IDs in the two files — for example, the same truth sample imputed at multiple downsampling depths — pass a `--samples` file with an optional second column that maps each imputed-side ID to the corresponding truth-side ID. Rows without a second column behave as before.

<div class="code-example" markdown="1">
```
# samples.txt — first column = imputed ID, optional second column = validation ID
NA12878.0_5x  NA12878
NA12878.1x    NA12878
NA12878.2x    NA12878
HG00096
HG00097
```
</div>

The three NA12878 rows alias to the single truth sample `NA12878`; `HG00096` and `HG00097` use the same ID in both files and need no alias. Per-sample output rows are labelled with the imputed-side ID.

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
| \-\-input            | FILE    | NA       | File with four columns listing in order: regions frequencies validation and imputed dataset. For genome-wide concordance, add more lines specifying different chromosomes. |
| \-\-samples          | FILE    | NA       | List of samples to process. One sample per line, giving the sample ID as it appears in the imputed VCF. Each row may optionally include a second whitespace-separated column giving the corresponding sample ID in the validation VCF (useful when the same biological sample appears under different IDs, e.g. multiple downsampled replicates of the same truth sample). The second column is optional on a per-row basis — a single file may freely mix rows with and without an alias. When omitted, the truth-side ID is assumed identical to the imputed-side ID. See "Sample aliasing" below for an example. |
| \-\-gt-val           | NA      | NA       | Uses hard called genotypes rather than phread-scaled likelihoods for the validation dataset, reading them from FORMAT/GT field. |
| \-\-gt-tar           | NA      | NA       | Uses FORMAT/GT field to determine the best-guess genotype rather than the FORMAT/GP (default). FORMAT/DS are FORMAT/GP fields are still required for calibration and rsquared calculations. |


#### Other parameters

| Option name 	      | Argument|  Default  | Description |
|:--------------------|:--------|:----------|:-------------------------------------|
| \-\-af-tag          | STRING  | AF        | Allele frequency INFO tag to use for binning. By default the allele frequency is estimated from the INFO/AF tag.  |
| \-\-use-alt-af      | NA      | NA        | If specified, the metrics work on the ALT allele frequency (range \[0,1\]), rather than minor allele frequency (range \[0,0.5\]).  |
| \-\-bins            | VECTOR  | NA        | Allele frequency bins used for rsquared computations. By default they should as MAF bins \[0-0.5\], while they should take the full range \[0-1\] if --use-ref-alt is used.  |
| \-\-ac-bins         | VECTOR  | NA        | User-defined allele count bins used for rsquared computations. |
| \-\-allele-counts   | VECTOR  | NA        | Default allele count bins used for rsquared computations. AN field must be defined in the frequency file. |
| \-\-min-val-gl      | FLOAT   | NA        | Minimum genotype likelihood probability P(G\|R) in validation data \[set to zero to have no filter of if using --gt-validation\] |
| \-\-min-val-dp      | INT     | NA        | Minimum coverage in validation data. If FORMAT/DP is missing and --minDP > 0, the program exits with an error. \[set to zero to have no filter of if using --gt-validation\] |
| \-\-min-tar-gp      | VECTOR  | NA        | Minimum GP probabilities to be used as a filter. By default it looks at the GP field to specify the filter, but will try to use FORMAT/PL if gt-tar option is specified. Leave empty if no filter is used. |
| \-\-out-r2-per-site | NA      | NA        | Output r2 at each site. |
| \-\-out-rej-sites   | NA      | NA        | Output sites where that cannot be used for the concordance.  |
| \-\-out-conc-sites  | NA      | NA        | Output sites where all target genotypes are concordant with the truth.  |
| \-\-out-disc-sites  | NA      | NA        | Output sites where at least one target genotype is diconcordant with the truth.  |
| \-\-groups          | FILE    | NA        | Alternative to frequency bins: group bins are user defined, provided in a file.  |


#### Output files

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-O \[\-\-output \]  | STRING  | NA       | Prefix of the output files (extensions are automatically added) |
| \-\-log              | STRING  | NA       | Log file  |

