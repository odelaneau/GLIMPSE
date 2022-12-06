---
layout: default
title: split_reference
nav_order: 2
parent: Documentation
---
# phase_rare
{: .no_toc .text-center }

{: .highlight }
Website under construction. A complete release of GLIMPSE2 will be available by the 7th of December 2022.


## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

### Description
Tool to create a binary reference panel for quick reading time.

### Usage
Simple run

<div class="code-example" markdown="1">
```bash
GLIMPSE2_split_reference --input-region chr20:7702567-12266861 --output-region chr20:7952603-12016861 --output binary_reference_panel --reference reference_panel_full_chr20.bcf --map chr20.b38.gmap.gz --threads 4
```
</div>

---

### Command line options

#### Basic options

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-\-help             | NA      | NA       | Produces help message |
| \-\-seed             | INT     | 15052011 | Seed of the random number generator  |
| \-T \[ \-\-thread \] | INT     | 1        | Number of threads |



#### Input parameters

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-\-input-region     | STRING  | NA       | Imputation region with buffers |
| \-\-output-region    | STRING  | NA       | Imputation region without buffers |
| \-R \[ \-\-reference \] | FILE    | NA       | Reference panel of haplotypes in VCF/BCF format |
| \-M \[ \-\-map \]    | FILE    | NA       | Genetic map |
| \-\-sparse-maf       | FLOAT   | =0.001   | Expert setting: rare variant threshold |
| \-\-keep-monomorphic-ref-sites | NA       | NA       | Keeps monomorphic markers in the reference panel, that are removed by default |

#### Output files

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-O \[\-\-output \]  | STRING  | NA       | Prefix of the output file (region and extension are automatically added) |
| \-\-log              | STRING  | NA       | Log file  |

