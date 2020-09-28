# GLIMPSE versioning changelog

### GLIMPSE phase

* v1.1.0
	* Better reporting of missing genotype likelihoods and ploidy
	* Added multi-threaded input for VCF/BCF files
	* Added haploid/diploid/mixed ploidy imputation
	* Added FPLOIDY parameter in output
	* Small improvements in the PBWT selection
* v1.0.1
	* Bugfix for compatibility with bcftools 1.10.x regarding --main 
* v1.0.0
	* First release. Version used for the paper

### GLIMPSE chunk

* v1.0.1
	* Removed option (--input) requiring input dataset and merged with --reference. Only reference panel is now required using the --input option.
* v1.0.0
	* First release. Version used for the paper

### GLIMPSE ligate

* v1.1.0
	* Performs indexing on input files if not indexed (avoiding to run bcftools after GLIMPSE_phase)
	* Added additional checks for FPLOIDY

* v1.0.0
	* First release. Version used for the paper

### GLIMPSE sample

* v1.0.0
	* First release. Version used for the paper

### GLIMPSE stats

* v1.0.0
	* First release.
