---
layout: default
title: chunk
nav_order: 1
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
Tool to create imputation chunks.

### Usage
Simple run

<div class="code-example" markdown="1">
```bash
GLIMPSE2_chunk --input file_chr20.bcf --map chr20.b38.gmap.gz --region chr20 --sequential --output chunks_chr20.txt
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


#### Input files

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-I \[ \-\-input \]  | FILE    | NA       | Reference or target dataset at all variable positions in VCF/BCF format. The GT field is not required |
| \-\-region           | STRING  | NA       | Chromosome or region to be split |
| \-M \[ \-\-map \]    | FILE    | NA       | Genetic map |
| \-\-sparse-maf       | FLOAT   | 0.001   | **Expert setting.** Rare variant threshold |

#### Window Parameters

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-\-window-cm        | FLOAT   | 4.0     | Minimal Window size in cM |
| \-\-window-mb        | FLOAT   | 4.0     | Minimal Window size in Mb |
| \-\-window-count     | INT     | 30000   | Minimal window size in #variants |
| \-\-buffer-cm        | FLOAT   | 0.5     | Minimal buffer size in cM |
| \-\-buffer-mb        | FLOAT   | 0.5     | Minimal buffer size in Mb |
| \-\-buffer-count     | INT     | 3000    | Minimal buffer size in #variants |

#### Model Parameters

| Option name 	              | Argument|  Default  | Description |
|:----------------------------|:--------|:----------|:-------------------------------------|
| \-\-recursive               | NA      | NA        | Recursive algorithm |
| \-\-sequential              | NA      | NA        | **Recommended.** Sequential algorithm|
| \-\-uniform-number-variants | NA      | NA        | **Experimental.** Uniformize the number of variants in the sequential algorithm |

#### Output files

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-O \[\-\-output \]  | STRING  | NA       | Region file containing the chunks for phasing and imputation |
| \-\-log              | STRING  | NA       | Log file  |

