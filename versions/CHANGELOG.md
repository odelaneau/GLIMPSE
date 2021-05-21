# GLIMPSE versioning changelog

### GLIMPSE phase
* v1.1.1
	* Removed bug in haploid imputation
	* Removed bug reading imputing regions with spanning indels (bug reported by @MySelvan)
	* Removed bug in INFO score calculation
	* Refactory model classes (cast)
	* Added --inputGL option
	* Added --ban-repeated-sample-names option

* v1.1.0
	* Better reporting of missing genotype likelihoods and ploidy
	* Added multi-threaded input for VCF/BCF files
	* Added haploid/diploid/mixed ploidy imputation
	* Added FPLOIDY parameter in output
	* Small improvements in the PBWT selection

* v1.0.1
	* Change in the state  selection for phasing
	* Bugfix for compatibility with bcftools 1.10.x regarding --main 

* v1.0.0
	* First release. Version used for the paper

### GLIMPSE chunk

* v1.1.1
	* Lifted version to v1.1.1
	* Removed warnings

* v1.1.0
	* Lifted version to v1.1.0
	* Added multithreading for reading files (htslib)

* v1.0.1
	* Removed option (--input) requiring input dataset and merged with --reference. Only reference panel is now required using the --input option.

* v1.0.0
	* First release. Version used for the paper

### GLIMPSE ligate

* v1.1.1
	* Lifted version to v1.1.1
	* Removed warnings

* v1.1.0
	* Added additional checks for FPLOIDY
	* Added haploid/diploid/mixed ploidy support
	* Performs indexing on input files if not indexed

* v1.0.0
	* First release. Version used for the paper

### GLIMPSE sample

* v1.1.1
	* Lifted version to v1.1.1
	* Removed warnings

* v1.1.0
	* Added haploid/diploid/mixed ploidy support
	* Added multithreading for reading files (htslib)

* v1.0.0
	* First release. Version used for the paper

### GLIMPSE concordance

* v1.1.1
	* Big refactory of the call and output format

* v1.1.0
	* Added ploidy 1 or 2 option

* v1.0.0
	* First release. Version used for the paper

### GLIMPSE stats

* v1.0.0
	* First release.
