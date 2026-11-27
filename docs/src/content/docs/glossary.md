---
title: "Glossary"
description: "Terms used in the GLIMPSE documentation."
---

## Genetic map

A map of genetic distance along a chromosome, in centimorgans (cM), which measures how often recombination separates two positions. GLIMPSE2 uses it to define chunk sizes and as part of the model, and `chunk` and `split_reference` both take one.

## Genotype dosage

The expected number of alternative alleles at a site, a value between 0 and 2 for a diploid sample. GLIMPSE2 writes it in `FORMAT/DS`, next to the best-guess genotype (`GT`) and the genotype probabilities (`GP`).

## Genotype imputation

Genotype imputation is a probabilistic inference of genotypes, typically used for SNP array data or low-coverage WGS. It uses a large reference panel of phased haplotypes to find long shared identity-by-descent segments, which carry the information used for inference. For SNP array data, see [IMPUTE5](https://jmarchini.org/software/); for low-coverage sequencing, GLIMPSE2.

## Genotype likelihood

The probability of the sequencing reads at a site given each possible genotype. At low coverage the reads leave the genotype uncertain, so GLIMPSE2 works with these likelihoods instead of hard genotype calls.

## Haplotype phasing

Haplotype phasing is the process of separating the two parentally inherited copies of each chromosome into haplotypes. [SHAPEIT5](https://odelaneau.github.io/shapeit5/) is the tool from the same family for phasing large sequencing cohorts.

## Low-coverage whole-genome sequencing

Whole-genome sequencing at low depth, usually between 0.1x and 8x. Unlike high-coverage WGS, genotype calls are uncertain at this depth, so the data is usually imputed with specialised methods such as GLIMPSE.

## Reference panel

A large set of deeply sequenced haplotypes, used for genotype imputation and haplotype phasing of cohorts typed on SNP arrays or sequenced at low coverage.

## Singleton

A rare variant carried by a single chromosome in the dataset (minor allele count of 1). There is no information in the data to phase it, so statistical phasing reports a switch error rate of 50% at these sites. Family data is the typical way to phase them.

## SNP array

A DNA microarray that detects common SNPs in a population. In GWAS, SNP array data is often imputed using a reference panel of haplotypes.

## Switch error rate

A common metric for the quality of phased data, used when the true maternal and paternal haplotypes are known. It is the number of switch errors divided by the number of opportunities for a switch error.

## Whole-genome sequencing (WGS)

Determining the DNA sequence of an entire genome. High-coverage WGS (30x or more) gives high-quality genotype calls.
