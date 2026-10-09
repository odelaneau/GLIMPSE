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
Inspect a GLIMPSE2 binary reference panel (`.bin` file produced by `GLIMPSE2_split_reference`), print summary statistics and, optionally, write the reference haplotypes it contains as VCF/BCF. Useful for debugging, validation, and understanding chunk characteristics, e.g. to see exactly which sites and haplotypes `GLIMPSE2_phase` was given for a chunk.

Statistics reported include:

- File size, chromosome, input/output regions (bp and Mbp)
- Genetic map span for input and output regions (cM), and whether the cM positions were interpolated from a genetic map or use the constant 1 cM/Mb fallback
- Number of haplotypes
- Variant counts (total, common, rare, common HQ, low quality)
- Variant types (SNP, MNP, indel, other)
- Allele-frequency distribution (monomorphic, singletons, MAC 2-5, MAF <1%, 1-5%, 5-50%)
- Core (output region) vs buffer variant breakdown

### Usage

<div class="code-example" markdown="1">
```bash
GLIMPSE2_inspect --input reference_panel/split/1000GP.chr22.noNA12878_chr22_15649863_19744259.bin
```
</div>

To also write the haplotypes of the chunk's core (output) region:

<div class="code-example" markdown="1">
```bash
GLIMPSE2_inspect --input reference_panel/split/1000GP.chr22.noNA12878_chr22_15649863_19744259.bin --output chunk_01.bcf
```
</div>

### Haplotype output

With `--output`, the reference haplotypes stored in the `.bin` file are written as VCF/BCF (format chosen from the file extension: `.bcf`, `.vcf.gz` or `.vcf`; compressed outputs are indexed with a `.csi` file).

- **Samples**: the `.bin` file stores neither sample names nor ploidy, so each haplotype is written as its own haploid sample, in the order of the original panel, named `h1`, `h2`, ... zero-padded to a fixed width (e.g. `h00001` to `h18000` for 18,000 haplotypes). For an all-diploid panel, haplotypes `2i-1` and `2i` are the two phased haplotypes of the `i`-th original sample.
- **Sites**: `CHROM`, `POS`, `ID`, `REF` and `ALT` as in the original panel. Only the sites `GLIMPSE2_split_reference` kept are present: multi-allelic records and (unless `--keep-monomorphic-ref-sites` was used) monomorphic records were dropped when the panel was split. `QUAL`, `FILTER`, other INFO and FORMAT fields (including `INFO/END`, so symbolic alleles lose their span), and unphased status are not stored.
- **INFO fields**: `AC`, `AN` and `AF` in the reference panel; `RARE` flags variants stored sparsely (minor allele frequency below `--sparse-maf`); `CM` is the genetic position relative to the first core variant of the chunk, not a position on the chromosome's genetic map (variants in the left buffer are negative). `CM` is a 32-bit float: text VCF shows 6 significant digits, so nearby variants can show the same value. `BUFFER` flags buffer variants when `--include-buffers` is used.
- **Region**: by default only the core region is written; `--include-buffers` also writes the variants in the buffers on either side.

---

### Command line options

#### Basic options

| Option name          | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-\-help             | NA      | NA       | Produces help message |
| \-T \[\-\-threads \] | INT     | 1        | Number of threads used to compress the haplotype output |

#### Input files

| Option name          | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-I \[\-\-input \]   | STRING  | NA       | Binary reference panel file (.bin) to inspect |

#### Output files

| Option name          | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-O \[\-\-output \]  | STRING  | NA       | Write the reference haplotypes of the core (output) region to this VCF/BCF file, one haploid sample per haplotype |
| \-\-include-buffers  | NA      | NA       | Also write the variants in the buffer regions, flagged with INFO/BUFFER |
| \-\-compression-level| INT     | 6        | Compression level for VCF/BCF output: 0 = none (still BGZF-framed and indexable), 1 = fastest, 9 = smallest. Ignored for plain .vcf output. |
| \-\-log              | STRING  | NA       | Log file  |
