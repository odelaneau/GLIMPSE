---
title: CHANGELOG
layout: default
---

# CHANGELOG

All notable changes to this project are documented in this file.

{: .highlight }
Version 2.0.0 will be released on the 7th of December 2022.

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
	* Change in the state  selection for phasing
	* Bugfix for compatibility with bcftools 1.10.x regarding --main

## v1.0.0
	* First release. Version used for the preprint (https://doi.org/10.1101/2020.04.14.040329)
