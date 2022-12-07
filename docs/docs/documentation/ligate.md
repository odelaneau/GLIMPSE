---
layout: default
title: ligate
nav_order: 4
parent: Documentation
---
# ligate
{: .no_toc .text-center }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

### Description
Ligatation of multiple phased BCF/VCF files into a single whole chromosome file. GLIMPSE2 is run in chunks that are ligated into chromosome-wide files maintaining the phasing.

### Usage
Simple run

<div class="code-example" markdown="1">
```bash
#ls -1v in order to keep the order within the chromosome
ls -1v chr20/*.imputed.bcf > list_imputed_files_chr20.txt

GLIMPSE2_ligate --input list_imputed_files_chr20.txt --output ligated_chr20.bcf --threads 2
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

#### Input files

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-I \[\-\-input \]   | STRING  | NA       | Text file containing all VCF/BCF to ligate, one file per line |

#### Output files

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-O \[\-\-output \]  | STRING  | NA       | Output ligated (phased) file in VCF/BCF format |
| \-\-no-index         | STRING  | NA       | If specified, the ligated VCF/BCF is not indexed by GLIMPSE2 for random access to genomic regions |
| \-\-log              | STRING  | NA       | Log file  |

