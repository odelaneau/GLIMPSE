/*******************************************************************************
 * @class genotype_reader
 * @brief Reads and processes genotype data from VCF/BCF files, integrating reference and target haplotypes.
 * 
 * This class handles:
 * - Parsing and scanning of genotype likelihoods or genotype calls from input files.
 * - Management of ploidy information for reference and target samples.
 * - Classification of variants into common and rare based on minor allele frequency.
 * - Integration with associated containers for haplotypes, genotypes, variants, and mpileup data.
 * 
 * It supports flexible input formats including GL (genotype likelihoods) and PL (phred-scaled likelihoods),
 * handles reference-only imputation, and manages monomorphic sites optionally.
 * 
 * Dependencies:
 * - haplotype_set, genotype_set, variant_map, glimpse_mpileup for data storage and manipulation.
 * 
 * Configuration parameters:
 * - sparse_maf: Threshold for minor allele frequency to define common variants.
 * - inputGL: Whether input data contains GL fields.
 * - impute_refonly: Whether to impute only reference panel variants.
 * - keep_mono: Whether to retain monomorphic reference sites.
 * - use_gl_indels: Whether to process indels from GL fields.
 * 
 * @note Sample ploidy can be customized using ploidy_samples map.
 * 
 ******************************************************************************/

#ifndef _GENOTYPE_READER_H
#define _GENOTYPE_READER_H

#include "otools.h"

#include "variant_map.h"
#include "containers/haplotype_set.h"
#include "containers/glimpse_mpileup.h"

class genotype_reader
{
public:
    /// Reference to haplotype set container.
	haplotype_set & H;
    /// Reference to genotype set container.
	genotype_set & G;
    /// Reference to variant map container.
	variant_map & V;
    /// Reference to glimpse mpileup data.
	glimpse_mpileup & M;

    /// Minor allele frequency threshold for defining common variants.
	const float sparse_maf;
    /// Flag indicating if input contains genotype likelihoods (GL).
	const bool inputGL;
    /// Flag indicating imputation mode (reference only).
	const bool impute_refonly;
    /// Flag to keep monomorphic reference sites.
	const bool keep_mono;
    /// Flag to include indels from GL data.
	const bool use_gl_indels;

    /// Number of samples in the reference panel.
	int n_ref_samples;
    /// Set of samples currently initializing.
	std::set < std::string > initializing_samples;
    /// Map from sample name to ploidy value.
	std::map < std::string, int> ploidy_samples;
    /// Vector holding ploidy per reference sample.
	std::vector<int> ploidy_ref_samples;

    /**
     * @brief Constructor
     * 
     * Initializes genotype_reader with references to data containers and config parameters.
     * 
     * @param H Reference to haplotype_set instance.
     * @param G Reference to genotype_set instance.
     * @param V Reference to variant_map instance.
     * @param M Reference to glimpse_mpileup instance.
     * @param _sparse_maf Minor allele frequency threshold to classify common variants.
     * @param _inputGL Flag for presence of genotype likelihoods in input.
     * @param _impute_refonly Flag to enable reference-only imputation.
     * @param keep_mono Flag to keep monomorphic reference sites.
     * @param use_gl_indels Flag to allow indel processing from GL fields.
     */
	genotype_reader(haplotype_set &, genotype_set &, variant_map &, glimpse_mpileup & M, const float _sparse_maf, const bool _inputGL, const bool _impute_refonly, const bool keep_mono, const bool use_gl_indels);

    /// Destructor, releases resources if any.
	~genotype_reader();

	//IO
	//void readInitializingSamples(string);
	/**
	 * @brief Reads a samples file and assigns ploidy values to samples.
	 * 
	 * This function reads a text file where each line contains two columns:
	 * the sample name and its corresponding ploidy value. It validates
	 * the ploidy against accepted values and stores the mapping of samples
	 * to their ploidy in the member variable `ploidy_samples`.
	 * 
	 * If a sample appears multiple times, only the first occurrence is kept,
	 * and a warning is issued.
	 * 
	 * @param ftext Path to the samples file to read.
	 * 
	 * @throws std::runtime_error if the file format is incorrect
	 *         or if an unrecognized ploidy value is encountered.
	 */
	void readSamplesFilePloidy(std::string);

	/**
	 * @brief Sets the ploidy information for reference samples from a BCF/VCF reader.
	 * 
	 * This function reads genotype data from the specified reference reader in
	 * the BCF/VCF streaming structure (`bcf_srs_t`) and determines the maximum ploidy
	 * for the reference samples. It validates that the ploidy is either 1 or 2,
	 * then assigns ploidy values per sample accordingly. Haploid and diploid counts
	 * are computed, and internal haplotype structures are allocated.
	 * 
	 * The function outputs verbose information on the number of target and reference
	 * samples and their haploid/diploid breakdowns.
	 * 
	 * @param sr Pointer to the BCF streaming structure containing reference data.
	 * @param id_ref Index of the reference reader in the streaming structure.
	 * 
	 * @throws std::runtime_error if the maximum ploidy is not 1 or 2, or if genotype
	 *         data format is inconsistent.
	 */
	void set_ploidy_ref(bcf_srs_t * sr, int id_ref);

	/**
	 * @brief Initializes and sets ploidy information for target samples.
	 * 
	 * This function configures ploidy settings for all target samples (`M.n_tar_samples`):
	 * - If `ploidy_samples` contains sample-specific ploidy info, it assigns each sample's
	 *   ploidy accordingly; otherwise, all samples default to diploid (ploidy = 2).
	 * - Calculates the number of diploid (`M.n_tar_diploid`) and haploid (`M.n_tar_haploid`) samples.
	 * - Sets index mappings for genotype (`M.tar_ind2gt`) and ploidy (`M.tar_ind2pl`), which
	 *   are used to track sample haplotype and genotype positions in data structures.
	 * - Computes total haplotypes count (`M.n_tar_haps`) and determines max ploidy (`M.max_ploidy`).
	 * - Updates the haplotype set (`H`) with sample counts, ploidy info, and index mappings.
	 * - Allocates genotype objects (`G.vecG`) for each sample, initializing them with sample
	 *   names, indices, number of variant sites (`H.n_tot_sites`), and ploidy-related data.
	 * 
	 * Note:
	 * - The flag `M.fploidy` indicates ploidy structure: positive for uniform ploidy,
	 *   negative if a mixture of haploid and diploid samples.
	 * - The function ensures the genotype structures are ready for downstream genotype
	 *   processing and imputation steps.
	 */
	void set_ploidy_tar();

	/**
	 * @brief Reads and processes genotype data from main and reference files.
	 * 
	 * This function initializes BCF/VCF readers for the main and reference genotype files,
	 * then performs two passes:
	 * 1. Scanning genotypes to gather summary information or preliminary parsing.
	 * 2. Parsing genotypes in detail for downstream processing.
	 * 
	 * Both scanning and parsing use a multi-threaded reader initialized with `nthreads`.
	 * Resources are properly allocated and released after each phase.
	 * 
	 * @param fmain Path to the main genotype file (VCF/BCF).
	 * @param fref Path to the reference genotype file (VCF/BCF).
	 * @param nthreads Number of threads to use for parallel reading.
	 */
	void readGenotypes(std::string , std::string , int nthreads);

	/**
	 * @brief Initializes the BCF/VCF streaming reader for main and reference files.
	 * 
	 * This function configures the streaming reader (`bcf_srs_t * sr`) to:
	 * - Disable collapsing of multi-allelic variants (`COLLAPSE_NONE`).
	 * - Require an index for efficient region-based access.
	 * - Set the genomic region of interest from the global variable `V.input_gregion`.
	 * - Set the number of threads for parallel reading if specified.
	 * 
	 * It then attempts to open both the main (`fmain`) and reference (`fref`) genotype files,
	 * verifying that files and their indexes are accessible. Errors are raised if files
	 * cannot be opened or indexes are missing.
	 * 
	 * @param sr Pointer to an initialized `bcf_srs_t` structure for streaming reading.
	 * @param fmain Path to the main genotype file (VCF/BCF).
	 * @param fref Path to the reference genotype file (VCF/BCF).
	 * @param nthreads Number of threads to use for parallel reading; ignored if <= 1.
	 * 
	 * @throws std::runtime_error if region setting fails, or if file/index loading fails.
	 */
	void initReader(bcf_srs_t * sr, std::string& fmain, std::string& fref, int nthreads);

	/**
	 * @brief Initializes the BCF/VCF streaming reader for a single genotype file.
	 * 
	 * Configures the streaming reader (`bcf_srs_t * sr`) to:
	 * - Require an index for efficient random access.
	 * - Set the genomic region of interest from `V.input_gregion`.
	 * - Set the number of threads for parallel reading if `nthreads > 1`.
	 * 
	 * Attempts to open the given genotype file and verifies index availability.
	 * Throws an error if the file or index cannot be loaded or if region setting fails.
	 * 
	 * @param sr Pointer to an initialized `bcf_srs_t` streaming reader structure.
	 * @param file Path to the genotype file (VCF/BCF) to open.
	 * @param nthreads Number of threads to use for parallel reading; ignored if <= 1.
	 * 
	 * @throws std::runtime_error on file or index loading failure, or if region setting fails.
	 */
	void initReader(bcf_srs_t * sr, std::string& file, int nthreads);

	/**
	 * @brief Scans genotype data from the primary BCF reader and initializes target samples.
	 * 
	 * This function performs the following steps:
	 * 1. Retrieves the number of target samples (`M.n_tar_samples`) from the header of the first
	 *    BCF reader (`sr->readers[0]`).
	 * 2. Extracts and stores the sample names in `M.tar_sample_names`.
	 * 3. Calls `scanGenotypesCommon(sr, 1)`, which performs additional common scanning tasks
	 *    (e.g., reading variant sites, initializing data structures).
	 * 4. Calls `set_ploidy_tar()` to assign ploidy values to each target sample based on
	 *    available ploidy information or defaults.
	 * 
	 * This function prepares internal data structures for downstream genotype processing
	 * and imputation.
	 * 
	 * @param sr Pointer to an initialized BCF streaming reader structure.
	 * 
	 * @see scanGenotypesCommon()
	 * @see set_ploidy_tar()
	 */
	void scanGenotypes(bcf_srs_t * sr);

	/**
	 * @brief Reads and parses target genotypes from a VCF/BCF file using a two-pass approach.
	 * 
	 * This function performs the following workflow:
	 * 
	 * - **Pass 1: Scanning**  
	 *   - Initialize reader to iterate variants and collect sample names.  
	 *   - Assign ploidy per sample.  
	 *   - For each variant, identify and store pointers to known variants or nullptr if unknown.  
	 * 
	 * - **Pass 2: Parsing**  
	 *   - Re-initialize reader to re-iterate variants.  
	 *   - Skip variants that are not biallelic or not present in known variants.  
	 *   - If configured, skip indels when `use_gl_indels` is false.  
	 *   - Read genotype likelihoods (GL) or phred-scaled likelihoods (PL) from FORMAT fields.  
	 *   - For each sample, convert likelihoods to internal encoded values as:  
	 *     @f[
	 *       \text{encoded} = \min\left(255, \left\lfloor -10 \times \text{likelihood} \right\rceil \right)
	 *     @f]
	 *   - Mark genotype sites as "flat" if all genotype likelihoods indicate homozygous calls.
	 * 
	 * Memory allocations for temporary arrays and streaming readers are managed and freed properly.
	 * 
	 * @param fmain Path to target genotype file (VCF/BCF).
	 * @param nthreads Number of threads used for parallel reading and parsing.
	 * 
	 * @note Processes only biallelic variants.
	 * @note Uses class/global structures:  
	 *       - `V`: Variant database  
	 *       - `M`: Sample and ploidy metadata  
	 *       - `H`: Haplotype structures  
	 *       - `G`: Genotype container  
	 *       - `vrb`: Verbose logger  
	 *       - `tac`: Timer utility  
	 * @note `inputGL` controls whether GL or PL likelihoods are parsed.
	 * @todo Review and test all TODO/FIXME comments.
	 */
	void readTarGenotypes(std::string , int);

	/**
	 * @brief Parses genotypes from a bcf_srs_t streaming reader for reference and target samples.
	 * 
	 * This function iterates through variant records in the streaming reader `sr` which contains
	 * both reference and target sample data. It performs genotype parsing, likelihood decoding,
	 * and updates internal data structures accordingly.
	 * 
	 * Workflow:
	 * - For each variant site:
	 *   - Check that the variant is biallelic and that required fields (AC, AN) are present in the reference.
	 *   - Skip monomorphic sites unless `keep_mono` is true.
	 *   - On the first variant, initialize reference sample ploidy using `set_ploidy_ref`.
	 *   - Parse reference genotypes to count alternate and reference allele counts, and track haplotype data.
	 *   - Validate that allele counts (AC/AN) match genotype data.
	 *   - If the target genotype data is available:
	 *     - Read genotype likelihoods (GL) if `inputGL` is true, or phred-scaled likelihoods (PL) otherwise.
	 *     - For each target sample and ploidy:
	 *       - Convert likelihood values to encoded format using:  
	 *         @f[
	 *           \text{encoded} = \min\left(255, \text{round}\left(-10 \times \text{likelihood}\right)\right)
	 *         @f]
	 *       - Assign encoded likelihoods to genotype objects.  
	 *       - Mark genotype site as "flat" if likelihoods indicate homozygous genotype calls.
	 * 
	 * @param sr Pointer to an initialized and populated bcf_srs_t streaming reader with at least
	 *           two readers: 0 for target and 1 for reference panel.
	 * 
	 * @throws Throws error if input files are null or required VCF INFO fields are missing or inconsistent.
	 * 
	 * @note Only supports biallelic variants (2 alleles).
	 * @note Uses class/global objects such as `H`, `G`, `M`, `V`, `vrb`, and `tac`.
	 * @note Relies on global flags: `inputGL`, `keep_mono`, and `impute_refonly`.
	 */
	void parseGenotypes(bcf_srs_t * sr);

	/**
	 * @brief Reads reference genotypes and associated BAM files with a two-pass approach.
	 * 
	 * This function processes the reference genotype file in two stages:
	 * 
	 * 1. **Scanning pass:**  
	 *    - Initializes a BCF streaming reader for the reference file.  
	 *    - Sets target sample ploidy information.  
	 *    - Performs a common genotype scan on the reference panel without target samples.  
	 *    - Cleans up the scanning reader.  
	 * 
	 * 2. **Parsing pass:**  
	 *    - Re-initializes a BCF streaming reader on the reference file.  
	 *    - Parses detailed reference genotypes using the dedicated parsing method.  
	 *    - Cleans up the parsing reader.  
	 * 
	 * @param fref Path to the reference genotype file (VCF/BCF).  
	 * @param nthreads Number of threads to use for file reading and processing.  
	 * 
	 * @note Assumes class/global members and helper functions:  
	 *       - `initReader()` to initialize BCF readers  
	 *       - `set_ploidy_tar()` to set ploidy for target samples  
	 *       - `scanGenotypesCommon()` to scan genotype data  
	 *       - `parseRefGenotypes()` to parse reference genotype details  
	 *       - `bcf_sr_destroy()` to free BCF streaming resources  
	 */
	void readGenotypesAndBAMs(std::string funphased, int nthreads);

	/**
	 * @brief Parses reference genotypes from a BCF streaming reader and updates internal haplotype structures.
	 * 
	 * This function iterates over variants in the reference genotype file to:
	 * - Skip variants that are not biallelic (only variants with 2 alleles are processed).
	 * - On the first variant, determine and set the ploidy of reference samples.
	 * - For each reference sample at each variant site:
	 *   - Extract genotype alleles.
	 *   - Update haplotype variant presence in either common or rare variant data structures:
	 *     - If variant is common (`H.flag_common[i_site]`), update haplotype variant matrix `HvarRef`.
	 *     - Otherwise, store variant index in `ShapRef` for rare variant haplotypes.
	 *   - Count reference (`cref`) and alternate (`calt`) alleles for validation.
	 * - Verify that the allele counts from genotypes match AC/AN INFO fields in the VCF to ensure data consistency.
	 * - Report parsing progress and timing.
	 * 
	 * @param sr Pointer to an initialized BCF streaming reader for the reference genotype file.
	 * 
	 * @throws Throws an error if `sr` is null or if AC/AN fields do not match genotype counts.
	 * 
	 * @note Uses global structures/classes:
	 *   - `H`: haplotype structures and variant flags
	 *   - `V`: variant positions and allele counts
	 *   - `vrb`: verbose logger
	 *   - `tac`: timer utility
	 */
	void parseRefGenotypes(bcf_srs_t * sr);

	/**
	 * @brief Scans variants from a BCF/VCF reference reader to classify sites and collect variant information.
	 * 
	 * This function iterates through variants in the reference file to:
	 * - Skip non-biallelic variants or incomplete data lines.
	 * - Extract allele count (AC) and allele number (AN) INFO fields to calculate allele frequencies.
	 * - Classify variants as common or rare based on minor allele frequency (MAF) threshold.
	 * 
	 * The minor allele frequency is calculated as:
	 * \[
	 * \text{MAF} = \min \left( \frac{C_{\text{ref}}}{C_{\text{ref}} + C_{\text{alt}}}, \frac{C_{\text{alt}}}{C_{\text{ref}} + C_{\text{alt}}} \right)
	 * \]
	 * where:
	 * - \( C_{\text{ref}} \) = count of reference alleles,
	 * - \( C_{\text{alt}} \) = count of alternate alleles.
	 * 
	 * - Variants with MAF ≥ \c sparse_maf are flagged as common, otherwise rare.
	 * - Monomorphic sites (with zero minor allele count) are skipped unless the \c keep_mono option is enabled.
	 * - Counts and flags for common and rare variants are maintained.
	 * - A sparse index mapping (\c H.common2tot) is built for common variants.
	 * - Variant details are stored in the \c V container.
	 * - Provides warnings and reports statistics on scanned variants.
	 * 
	 * @param sr Pointer to an initialized BCF streaming reader.
	 * @param ref_sr_n The index of the reader corresponding to the reference panel in the streaming reader.
	 * 
	 * @throws Throws error if AC/AN INFO fields are missing or if no variants are found.
	 * 
	 * @note Uses global variables/structures:
	 * - \c H for haplotype and variant summary statistics.
	 * - \c V for variant storage.
	 * - \c vrb for verbose logging.
	 * - \c stb for string conversions.
	 * 
	 * @see variant
	 */
	void scanGenotypesCommon(bcf_srs_t * sr, int ref_sr_n);
};

#endif
