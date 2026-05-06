---
layout: default
title: inspect
nav_order: 6
parent: Documentation
---
# inspect
{: .no_toc .text-center }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

### Description
Inspect a GLIMPSE2 binary reference panel (`.bin` file produced by `GLIMPSE2_split_reference`) and print summary statistics. Useful for debugging, validation, and understanding chunk characteristics.

Statistics reported include:

- File size, chromosome, input/output regions (bp and Mbp)
- Genetic map span for input and output regions (cM)
- Number of haplotypes
- Variant counts (total, common, rare, common HQ, low quality)
- Variant types (SNP, MNP, indel, other)
- Allele-frequency distribution (monomorphic, singletons, MAC 2-5, MAF <1%, 1-5%, 5-50%)
- Core (output region) vs buffer variant breakdown

### Usage

<div class="code-example" markdown="1">
```bash
GLIMPSE2_inspect --input reference_panel/split/1000GP.chr22.noNA12878_chr22_16050075_17084716.bin
```
</div>

---

### Command line options

#### Basic options

| Option name          | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-\-help             | NA      | NA       | Produces help message |

#### Input files

| Option name          | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-I \[\-\-input \]   | STRING  | NA       | Binary reference panel file (.bin) to inspect |

#### Output files

| Option name          | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-\-log              | STRING  | NA       | Log file  |
