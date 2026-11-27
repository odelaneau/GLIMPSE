---
title: CHANGELOG
layout: default
---

# CHANGELOG

All notable changes to this project are documented in this file.

## v2.1.0

GLIMPSE2 2.1 imputes substantially faster and with much less memory than 2.0.1 (see Performance), runs natively on ARM64, and fixes many bugs, several of which change imputation, ligation, chunking or concordance output. Changes below are relative to v2.0.1.

### Performance

`GLIMPSE2_phase` is significantly faster and uses much less memory in every case we measured, but by how much depends strongly on the reference panel (its size and variant density), the sequencing depth, and the invocation: one sample or many, BAM/CRAM or `--input-gl`, and the number of threads. The table gives a few examples: 2.1 against v2.0.1, imputing a 3.2 Mb region of chr22 from 1x data.

| Reference panel | Input | Samples, threads | x86-64 speed-up | x86-64 memory | ARM64 speed-up | ARM64 memory |
|---|---|---|---|---|---|---|
| HGDP+1KG, 8,176 haplotypes | BAM | 1, 1 | 2.3× | 5.3× less | 1.8× | 5.2× less |
| Simulated, 10,000 haplotypes | GLs | 1, 1 | 2.3× | 5.3× less | 1.7× | 5.2× less |
| Simulated, 150,000 haplotypes | GLs | 1, 1 | 5.1× | 2.1× less | 3.9× | 2.1× less |
| HGDP+1KG, 8,176 haplotypes | BAM | 100, 4 | 2.2× | 4.5× less | 1.5× | 4.6× less |
| Simulated, 10,000 haplotypes | GLs | 100, 4 | 2.2× | 3.5× less | 1.5× | 3.6× less |
| Simulated, 150,000 haplotypes | GLs | 100, 4 | 2.4× | 2.2× less | 1.7× | 2.3× less |

For example, one sample against the 150,000-haplotype panel took 83 s and 2.4 GB of memory with v2.0.1, and 16 s and 1.1 GB with 2.1, on x86-64. Memory is peak resident memory. Each figure is the median of two runs on AWS m8a (AMD EPYC, g++ 13) or m8g (Graviton4, clang 23). v2.0.1 does not build on ARM64, so the ARM64 baseline is v2.0.1 plus the changes needed to build there.

The gains come from a faster sparse-PBWT state selection (which dominates run time on large panels), faster construction of the conditioning states and HMM inner loops, and an imputation HMM that keeps only every 64th forward row and recomputes the rest (`--hmm-checkpoint-block`). None of these changes alters output.

### New

* **ARM64 support**: Linux on aarch64 (e.g. AWS Graviton) and macOS on Apple silicon. On aarch64 and macOS the default compiler is now `clang++`, which builds a much faster `GLIMPSE2_phase` on aarch64 than `g++` (about 2× end to end on Graviton4); set `CXX` to override.
* **`GLIMPSE2_inspect`**, a new tool that summarises a binary reference panel (`.bin`): region, genetic-map span, number of haplotypes, variant types and allele-frequency spectrum. With `--output` it also writes the panel's haplotypes as VCF/BCF.
* GLIMPSE2_phase: **checkpointing**. `--checkpoint-file-out` writes a checkpoint after each iteration and `--checkpoint-file-in` resumes from one, e.g. after a preempted cloud job. With one thread, a resumed run gives the same output as an uninterrupted one.
* GLIMPSE2_phase: **retries transient read errors** (with exponential backoff) when reading BAM/CRAM, `--input-gl` and binary reference files, e.g. when streaming from cloud storage, and stops with an error if a BAM/CRAM cannot be read to the end instead of imputing from the part it read. GLIMPSE2_ligate and GLIMPSE2_concordance now report read errors instead of stopping early without one.
* GLIMPSE2_phase: `--hmm-checkpoint-block` (default 64) sets the memory/recomputation trade-off of the imputation HMM; 0 stores every forward row, as 2.0.1 did.
* `--compression-level` for the VCF/BCF output of GLIMPSE2_phase, GLIMPSE2_ligate, GLIMPSE2_concordance and GLIMPSE2_inspect.
* `-T` as a short form of `--threads` in every tool, and `-I` for `GLIMPSE2_ligate --input`.
* GLIMPSE2_ligate: the input chunks no longer need to be listed in order; they are sorted by contig and position.
* GLIMPSE2_concordance: an optional second column in the `--samples` file names the validation sample an imputed sample is compared with, so one truth sample can score several imputed copies (e.g. the same sample at different coverages).
* Precomputed GRCh38 chunks (4 cM windows) in `resources/chunks/b38/4cM/`.

### Changes to output

Results from 2.1 differ from 2.0.1 for the same inputs and seed:

* GLIMPSE2_phase with BAM/CRAM input: genotype likelihoods are computed differently, fixing several bugs.
  * Split multi-allelic SNPs now get likelihoods; every ALT record after the first at a position was left flat.
  * Haploid samples use all their reads.
  * Bases above `--max-depth` are downsampled at random instead of keeping the first ones.
  * PLs are normalised by the best likelihood and rounded.
  * Overlapping read mates are counted once.
  * Base qualities are capped at Q40.
* GLIMPSE2_phase: two bugs fixed in the rephasing of heterozygous sites without reads.
  * Both sampled alleles were written to the first haplotype, leaving the second unchanged ([#297](https://github.com/odelaneau/GLIMPSE/issues/297)).
  * A running sum was not reset between calls.
* GLIMPSE2_phase: INFO and GP.
  * The INFO score is calculated correctly ([#144](https://github.com/odelaneau/GLIMPSE/issues/144)).
  * INFO is now also computed for haploid-only output (FPLOIDY=1, e.g. male chrX non-PAR), where it was 1 at every site, and for mixed-ploidy output, which had none.
  * FORMAT/GP is declared `Number=G`.
* GLIMPSE2_ligate:
  * Fixed switch errors introduced where chunks overlap.
  * Phase alignment starts afresh at each chromosome, so a chromosome's output no longer depends on the chromosomes before it. On later chromosomes the haplotype order of a sample can differ; genotypes do not.
  * All-haploid chunks (FPLOIDY=1) are concatenated without phase alignment; consecutive samples were treated as one diploid sample.
* GLIMPSE2_chunk: chunk boundaries change.
  * Non-integer `--buffer-cm` values are used as given; they were rounded up, giving larger buffers.
  * `--buffer-count` counts common variants, as documented; it counted all variants.
  * Recursive chunking no longer loses the chunks it had already made.
  * `--window-cm 0` and `--buffer-cm 0`, recommended without a genetic map, are accepted.
* GLIMPSE2_concordance:
  * `rsquared.spl` reports r² (it reported r).
  * Calibration bins hold the GP range they are labelled with; most GPs were reported one bin low.
  * The second site is no longer counted twice with `--ac-bins`/`--allele-counts` or a mixed-ploidy imputed file.
  * `--out-r2-per-site` covers every `--input` line instead of only the last ([#179](https://github.com/odelaneau/GLIMPSE/issues/179)). `--out-rej-sites`, `--out-conc-sites` and `--out-disc-sites`, which also kept only the last line, stop with an error when given more than one `--input` line.
  * Output headers gain the columns they were missing, and typos in column names are fixed.
  * Per-sample rows follow the imputed file's sample order.

### Fixes

* GLIMPSE2_phase:
  * No longer crashes on a region where the genetic map is flat (0 cM, e.g. the chr9 centromere) ([#304](https://github.com/odelaneau/GLIMPSE/issues/304)). A region with no common SNPs for the PBWT selection stops with an error.
  * With BAM/CRAM input, records with `*`, symbolic or `N` alleles are skipped instead of stopping the run.
  * With BAM/CRAM input and a VCF/BCF reference panel, a monomorphic record in the region no longer stops the run with "AC/AN INFO fields in VCF are inconsistent with GT field".
  * A reference panel without samples (e.g. a sites-only VCF) stops with a clear error instead of a floating point exception ([#220](https://github.com/odelaneau/GLIMPSE/issues/220)); GLIMPSE2_split_reference also gives the clear error.
  * Fixed runaway memory use and a crash from an integer overflow with panels of at most `--Kpbwt` haplotypes and many heterozygous sites without reads.
  * A `--samples-file` that cannot be opened or read to the end is an error; every sample was silently treated as diploid.
* GLIMPSE2_chunk: small regions no longer produce an empty chunk file ([#170](https://github.com/odelaneau/GLIMPSE/issues/170)).
* GLIMPSE2_split_reference and GLIMPSE2_chunk: a read error, e.g. from a truncated input file, is an error instead of being taken for the end of the file.
* GLIMPSE2_ligate: chunks nested inside another chunk, and chunks with different FPLOIDY (e.g. chrX PAR and non-PAR together), are rejected with an error.
* GLIMPSE2_concordance:
  * Validation and imputed records with missing or too few FORMAT values (GT, GP, PL) are rejected instead of crashing or reusing another sample's or record's values.
  * Haploid genotypes in a mixed-ploidy imputed file subset with `--samples` are no longer all rejected as NaN.
  * Fixed a hang when the first region has no record in the frequency or imputed file.
* Every tool exits with status 1 on a command-line error (unrecognised option, invalid value, missing argument); they exited 0.

### Compatibility

* GLIMPSE2_ligate: `--no-index` has been removed; the output is always indexed.
* GLIMPSE2_concordance: some output column names have changed (typo fixes).
* Command-line errors now exit with status 1 (see Fixes), so scripts that relied on exit 0 will stop there.
* Binary reference panels (`.bin`) made with 2.0.1 load in 2.1 unchanged.

### Build

* Build settings are shared in `common.mk`; `make` builds every tool, `make system` works across Linux distributions, and HTSlib and Boost locations can be set on the command line. The Dockerfile is updated.
* Every tool reports version 2.1.0 (the 2.0.1 binaries reported 2.0.0).
* The command-line reference in the documentation is brought up to date with the tools.

## v2.0.1

Released 2023-05-05 (tag `v2.0.1`, a GitHub pre-release; this is the version on Bioconda). Its binaries still report version 2.0.0.

* GLIMPSE2_chunk: new defaults give smaller chunks: `--window-cm 2.5` (was 4), `--window-mb 2` (was 4), `--window-count 20000` common variants (was 30000), `--buffer-mb 0.4` (was 0.5) and `--buffer-count 2000` (was 3000).
* GLIMPSE2_chunk: consecutive output regions meet with no gap between them, and on a whole chromosome the first starts at position 1 and the last ends at the contig length when the header gives it. Fixed regions smaller than one window. Without a genetic map it warns, and advises setting `--window-cm` and `--buffer-cm` to 0.
* GLIMPSE2_ligate: haploid samples are accepted; chunks with them stopped with "Non-diploid samples found in overlap".
* GLIMPSE2_phase: fixed "States for individual 0 are zero. Error during selection." with small reference panels ([#121](https://github.com/odelaneau/GLIMPSE/issues/121)).
* GLIMPSE2_phase: with `--input-gl`, likelihoods at indels and other non-SNP records are ignored (treated as flat) by default; the new `--use-gl-indels` uses them, as 2.0.0 did.
* GLIMPSE2_phase: fixed the output format for `--output` names ending in `.vcf.gz`, and gives a clearer error when a binary reference panel cannot be read.
* The binary reference panel format changed: panels made with 2.0.0's GLIMPSE2_split_reference must be rebuilt.
* Positions in genetic maps and chunking use 64-bit integers.
* Build: executables are written to each tool's `bin/`; Dockerfile fixes.

## v2.0.0
* Major release. Introduced speedups and accuracy improvements. Version used for the preprint. (https://doi.org/10.1101/2022.11.28.518213)

## v1.1.1
* Removed bug in haploid imputation
* Removed bug reading imputing regions with spanning indels (bug reported by @MySelvan)
* Removed bug in INFO score calculation
* Refactory model classes (cast)
* Added --inputGL option
* Added --ban-repeated-sample-names option

## v1.1.0
* Better reporting of missing genotype likelihoods and ploidy
* Added multi-threaded input for VCF/BCF files
* Added haploid/diploid/mixed ploidy imputation
* Added FPLOIDY parameter in output
* Small improvements in the PBWT selection

## v1.0.1
* Version used for the paper (https://doi.org/10.1038/s41588-020-00756-0)
* Change in the state selection for phasing
* Bugfix for compatibility with bcftools 1.10.x regarding --main

## v1.0.0
* First release. Version used for the preprint (https://doi.org/10.1101/2020.04.14.040329)
