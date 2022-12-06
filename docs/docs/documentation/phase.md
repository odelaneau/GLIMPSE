---
layout: default
title: phase
nav_order: 2
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



#### Model parameters

| Option name 	      | Argument|  Default  | Description |
|:--------------------|:--------|:----------|:-------------------------------------|


#### Selection parameters

| Option name 	      | Argument|  Default  | Description |
|:--------------------|:--------|:----------|:-------------------------------------|

#### BAM/CRAM options and filters

| Option name 	      | Argument|  Default  | Description |
|:--------------------|:--------|:----------|:-------------------------------------|


#### Output parameters

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-O \[\-\-output \]  | STRING  | NA       | Phased and imputed haplotypes in VCF/BCF/BGEN format |
| \-\-log              | STRING  | NA       | Log file  |


