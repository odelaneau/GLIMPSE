/**
 * @class caller
 * @brief Main class handling genotype calling, phasing, and imputation operations.
 *
 * This class encapsulates all the data structures, options, and methods required
 * to process genomic data, including reading input files, performing multi-threaded
 * computations, and outputting results.
 *
 * **Responsibilities:**
 * - Parse and store command-line options.
 * - Manage checkpointing using checksums.
 * - Store and process haplotype/genotype/variant data.
 * - Perform multi-threaded computations for imputation and phasing.
 * - Track and manage computation stages and statistics.
 *
 * **Modules Used:**
 * - `haplotype_set`, `genotype_set`, `variant_map` for genomic data storage.
 * - `phasing_hmm`, `imputation_hmm` for statistical modeling.
 * - `genotype_bam_caller` for BAM file genotype calling.
 *
 * @note Thread safety is partially ensured using `pthread_mutex_t` for worker coordination.
 */
#ifndef _CALLER_H
#define _CALLER_H

#include "otools.h"
#include "checksum_utils.h"

#include "containers/genotype_set.h"
#include "containers/haplotype_set.h"
#include "variant_map.h"
#include "io/genotype_reader.h"
#include "io/genotype_bam_caller.h"

#include "models/phasing_hmm.h"
#include "models/imputation_hmm.h"

class caller {
public:
	/** @name Command-line Options */
	///@{
	bpo::options_description descriptions; /**< Parsed command-line option descriptions. */
	bpo::variables_map options;            /**< Parsed command-line options and their values. */
	///@}

	/** @name Checkpointing */
	///@{
	checksum crc; /**< Checksum object used for verifying data consistency between runs. */
	///@}

	/** @name Internal Data */
	///@{
	haplotype_set H;      /**< Storage for haplotypes. */
	genotype_set G;       /**< Storage for genotypes. */
	variant_map V;        /**< Variant mapping and indexing. */
	glimpse_mpileup M;    /**< MPileup data handler for read-level processing. */

	InputFormat input_fmt;           /**< Format of input genomic data. */
	OutputFormat output_fmt;         /**< Format for output data. */
	OutputCompression output_compr;  /**< Compression type for output files. */
	int bgen_bits;                   /**< Number of bits for BGEN genotype encoding. */
	///@}

	/** @name Multi-threading */
	///@{
	int i_workers;                         /**< Number of worker threads. */
	int i_jobs;                            /**< Number of jobs to be processed. */
	std::vector<pthread_t> id_workers;     /**< Thread IDs for workers. */
	pthread_mutex_t mutex_workers;         /**< Mutex for synchronizing worker threads. */
	///@}

	/** @name Filtering Parameters */
	///@{
	float min_gl; /**< Minimum genotype likelihood threshold for inclusion. */
	///@}

	/** @name Computation Tracking */
	///@{
	int current_stage;                /**< Current stage of computation pipeline. */
	int current_iteration;            /**< Current iteration in the current stage. */
	int iterations_per_stage[3];      /**< Number of iterations for each computation stage. */
	stats1D statH;                     /**< Statistics for haplotypes. */
	stats1D statC;                     /**< Statistics for conditioning sets. */
	std::vector<std::vector<float>> HP0; /**< Haplotype posterior probabilities for state 0. */
	std::vector<std::vector<float>> HP1; /**< Haplotype posterior probabilities for state 1. */
	std::vector<std::vector<float>> HLC; /**< Conditional haplotype likelihoods. */
	std::vector<conditioning_set *> COND; /**< Conditioning states. */
	std::vector<imputation_hmm *> HMM;    /**< Imputation HMM objects. */
	std::vector<phasing_hmm *> DMM;       /**< Phasing HMM objects. */
	std::vector<genotype_bam_caller *> READER_BAM; /**< BAM genotype caller readers. */
	///@}

	/** @name Constructor & Destructor */
	///@{
	caller();  /**< Default constructor. */
	~caller(); /**< Destructor. */
	///@}

	//METHODS

	/**
	 * @brief Phases haplotypes and calculates genotype posteriors for a single individual.
	 *
	 * Workflow:
	 * 1. Select the individual for processing based on current stage.
	 * 2. Initialize or compute haplotype likelihoods depending on ploidy and stage.
	 * 3. Use HMM to compute posterior probabilities for haplotypes.
	 * 4. Sample haplotypes (H0 and H1 if diploid) and optionally rephase them.
	 * 5. Store genotype posteriors and haplotypes in main stage.
	 * 6. Update statistics (state count and polymorphic site percentage) with thread safety.
	 *
	 * Handles both haploid and diploid genomes.
	 *
	 * @param id_worker Thread ID running the phasing.
	 * @param id_job    ID of the individual to phase.
	 */
	void phase_individual(const int, const int);

	/**
	 * @brief Executes one iteration of the phasing process for genetic data.
	 *
	 * This function performs haplotype phasing on the current dataset, either in
	 * the initialization stage or in subsequent processing stages. It handles
	 * multi-threading, progress reporting, and checkpoint saving.
	 *
	 * @details
	 * - If in the initialization stage (`STAGE_INIT`), the function initializes rare
	 *   target haplotypes and performs initial rare variant selection.
	 * - Otherwise, it updates haplotypes, transposes rare target data, and matches
	 *   haplotypes using PBWT-based methods.
	 * - Supports both single-threaded and multi-threaded execution.
	 * - Tracks and reports progress, average number of states, and percentage of polymorphic sites.
	 * - Saves a checkpoint at the end of each iteration.
	 * - Clears initial states after the initialization stage to free memory.
	 *
	 * @note This method modifies the internal state of the `caller` object,
	 *       including haplotype storage, statistics, and stage progress.
	 *
	 * @see phase_individual()
	 * @see write_checkpoint()
	 *
	 * @warning Requires that input data `G` and `V` are already loaded and initialized.
	 *
	 * @dot
	 * digraph phase_iteration_flow {
	 *     node [shape=box, style="rounded,filled", fillcolor=lightgray];
	 *     start [label="Start"];
	 *     check_stage [label="current_stage == STAGE_INIT?"];
	 *     init_stage [label="H.initRareTar(G,V)\nH.performSelection_RARE_INIT_GL(V)"];
	 *     main_stage [label="H.updateHaplotypes(G)\nH.transposeRareTar()\nH.matchHapsFromCompressedPBWTSmall(V, stage==STAGE_MAIN)"];
	 *     threading [label="Multi-threaded or single-threaded\nphase_individual() calls"];
	 *     stats [label="Report statistics & progress"];
	 *     checkpoint [label="write_checkpoint()"];
	 *     cleanup [label="If STAGE_INIT:\nClear init_states"];
	 *     end [label="End"];
	 *
	 *     start -> check_stage;
	 *     check_stage -> init_stage [label="Yes"];
	 *     check_stage -> main_stage [label="No"];
	 *     init_stage -> threading;
	 *     main_stage -> threading;
	 *     threading -> stats;
	 *     stats -> checkpoint;
	 *     checkpoint -> cleanup;
	 *     cleanup -> end;
	 * }
	 * @enddot
	 *
	 * Example usage:
	 * @code
	 * caller c;
	 * c.phase_iteration();
	 * @endcode
	 */
	void phase_iteration();

	/**
	 * @brief Executes the main phasing loop across all stages.
	 *
	 * This method coordinates the phasing process by:
	 * - Initializing iteration counters and stage trackers.
	 * - Attempting to load a saved checkpoint to resume from a previous run.
	 * - Iterating through the defined phasing stages (`STAGE_INIT` to `STAGE_MAIN`).
	 * - Calling `phase_iteration()` for each iteration in the loop.
	 * - Incrementing iteration counters after each call.
	 * - Performing final post-processing steps on the genotype data after all stages are complete.
	 *
	 * ## Workflow:
	 * 1. **Initialization**
	 *    - Sets the current stage to `STAGE_INIT`.
	 *    - Sets iteration counter to `-1` so the first increment starts at 0.
	 *    - Reads checkpoint data if available.
	 *
	 * 2. **Main Loop**
	 *    - Repeatedly calls `phase_iteration()` for the current stage.
	 *    - Moves to the next iteration until all stages are processed.
	 *
	 * 3. **Finalization**
	 *    - For each genotype dataset in `G.vecG`, applies:
	 *        - Sorting
	 *        - Normalization
	 *        - Genotype inference
	 *
	 * @note
	 * - The checkpoint mechanism allows the loop to resume from a saved state.
	 * - This function assumes that `G.vecG` contains valid genotype objects.
	 *
	 * @see phase_iteration()
	 * @see read_checkpoint_if_available()
	 * @see increment_iteration()
	 */
	void phase_loop();

	/**
	 * @brief Advances the iteration counter and transitions between stages if needed.
	 *
	 * This function increments the current iteration number within the current stage.
	 * If the iteration count exceeds the number of iterations allocated for the stage,
	 * the stage counter is incremented and the iteration counter is reset to zero.
	 *
	 * ## Workflow:
	 * 1. Increment `current_iteration` by one.
	 * 2. If `current_iteration` reaches or exceeds `iterations_per_stage[current_stage]`:
	 *    - Increment `current_stage` by one (up to `STAGE_MAIN`).
	 *    - Reset `current_iteration` to zero.
	 *
	 * @note
	 * - This function assumes that `iterations_per_stage` contains valid entries for all stages.
	 * - The maximum stage is `STAGE_MAIN`; after that, no further stage transitions occur.
	 *
	 * @see phase_loop()
	 */
	void increment_iteration();

	//PARAMETERS
	/**
	 * @brief Declares and organizes command-line options for the program.
	 *
	 * This method sets up all command-line arguments using Boost.Program_options,
	 * grouping them into logical categories. Each group corresponds to a set of
	 * related configuration parameters (e.g., input files, algorithm settings,
	 * filtering options, and output controls).
	 *
	 * The options are organized into the following categories:
	 * 
	 * - <b>Basic parameters (`opt_base`)</b>
	 *   - Help, RNG seed, number of threads.
	 *
	 * - <b>Input parameters (`opt_input`)</b>
	 *   - Paths to BAM/CRAM/GL files, reference panels, genetic maps, region selectors,
	 *     rare variant thresholds, sample ploidy definitions, and checkpoint imports.
	 *   - Includes mutually exclusive input modes (`--input-gl`, `--bam-file`, `--bam-list`).
	 *
	 * - <b>VCF/BCF genotype likelihood input parameters (`opt_vcf_input`)</b>
	 *   - Flags controlling the handling of reference-only variants, GL/PL field usage,
	 *     and indel likelihood usage.
	 *
	 * - <b>Model parameters (`opt_algo`)</b>
	 *   - Number of burn-in and main iterations, population size, likelihood thresholds,
	 *     and HMM error rates for imputation and phasing.
	 *
	 * - <b>Selection parameters (`opt_selection`)</b>
	 *   - PBWT search depth, PBWT interval in centimorgans, initial and per-iteration
	 *     state limits, and optional state-list files for fixed haplotypes.
	 *
	 * - <b>BAM/CRAM options and filters (`opt_filters`)</b>
	 *   - Read calling model, indel calling, mapping/base quality thresholds,
	 *     maximum depth limits, and toggles for advanced read filtering behavior.
	 *
	 * - <b>Output parameters (`opt_output`)</b>
	 *   - Output file format, contig copying from FASTA index, BGEN encoding and compression
	 *     parameters, log file path, and checkpoint export file.
	 *
	 * Once defined, all groups are aggregated into the member variable `descriptions`.
	 *
	 * @note
	 * - Default values are set for many expert options.
	 * - Some options are mutually exclusive and must be validated at runtime.
	 *
	 * @see increment_iteration()
	 * @see phase_loop()
	 */
	void declare_options();

	/**
	 * @brief Parse and process command-line arguments for the GLIMPSE2 phase tool.
	 *
	 * This function uses Boost.Program_options to parse command-line arguments,
	 * configure runtime options, initialize logging, and display program metadata.
	 *
	 * **Main responsibilities:**
	 * - Parse provided arguments using predefined option descriptions.
	 * - Handle parsing errors and terminate if arguments are invalid.
	 * - Open a log file if the `--log` option is specified.
	 * - Display software metadata, authorship, version, commit information, and citations.
	 * - Show the help message if the `--help` flag is provided.
	 *
	 * **Options handled:**
	 * - `--log <filename>`: Write program output to the specified log file.
	 * - `--help`: Display usage information and exit.
	 *
	 * **Metadata displayed:**
	 * - Tool name and purpose.
	 * - Authors and contact details.
	 * - Version, commit ID, and commit date.
	 * - Citation references.
	 * - Current run date.
	 *
	 * @param args A vector of strings containing the command-line arguments
	 *        (including the program name as `args[0]`).
	 *
	 * @note If `--help` is specified, the function will display the help text
	 *       and immediately terminate the program.
	 * @note If invalid arguments are passed, an error message will be printed
	 *       to `stderr` and the program will exit.
	 *
	 * @throws boost::program_options::error If command-line parsing fails.
	 *
	 * @see descriptions
	 * @see vrb
	 * @see tac
	 */
	void parse_command_line(std::vector < std::string > &);

	/**
	 * @brief Validate and check consistency of parsed command-line options.
	 *
	 * This method performs an extensive set of validations on command-line
	 * parameters stored in the `options` variable map (parsed earlier by
	 * parse_command_line()). It ensures that required parameters are present,
	 * mutually exclusive options are respected, and parameter values are within
	 * valid ranges. It also determines the correct input and output file formats
	 * based on file extensions, and applies constraints specific to each format.
	 *
	 * **Main validation steps:**
	 *  - **Required input files:**
	 *      - Exactly one of `--bam`, `--bam-list`, or `--input-gl` must be provided.
	 *  - **Required output:**
	 *      - `--output` must be specified.
	 *  - **Reference panel:**
	 *      - `--reference` must be specified and must be in a supported format:
	 *        `.bcf`, `.vcf`, `.vcf.gz`, or binary GLIMPSE format.
	 *  - **Random seed:**
	 *      - If provided via `--seed`, it must be positive.
	 *  - **Thread count:**
	 *      - `--threads` must be ≥ 1.
	 *  - **Iteration constraints:**
	 *      - `--main` cannot exceed 15 (to extend beyond, `--burn` must be increased).
	 *
	 * **File format detection:**
	 *  - **Reference panel input:**
	 *      - `.bcf` or `.vcf` → InputFormat::BCF
	 *      - `.vcf.gz` → InputFormat::BCF
	 *      - Any other extension → InputFormat::GLIMPSE
	 *  - **Output:**
	 *      - `.vcf` → OutputFormat::VCF, OutputCompression::NONE
	 *      - `.bcf` → OutputFormat::BCF, OutputCompression::ZLIB
	 *      - `.vcf.gz` → OutputFormat::VCF, OutputCompression::ZLIB
	 *      - `.bgen` → OutputFormat::BGEN with compression type (`no`, `zlib`, `zstd`) and bits (`--bgen-bits` in 1..32)
	 *      - Other formats are rejected.
	 *      - If compiled without BGEN support, `.bgen` output triggers an error.
	 *
	 * **FASTA index validation:**
	 *  - If `--contigs-fai` is specified, verifies that the `.fai` index is readable and valid.
	 *
	 * **Region and MAF constraints:**
	 *  - For non-GLIMPSE input formats:
	 *      - Must specify `--input-region` and `--output-region`.
	 *      - Must specify `--sparse-maf` within (0, 0.5).
	 *  - For GLIMPSE binary input:
	 *      - `--input-region`, `--output-region`, `--sparse-maf`, and `--map` must NOT be specified manually.
	 *
	 * **State list and PBWT parameters:**
	 *  - If `--state-list` is provided:
	 *      - The file must be readable.
	 *      - `--Kinit` and `--Kpbwt` must be ≥ 0.
	 *      - If `Kpbwt > 0`, `--pbwt-depth` and `--pbwt-modulo-cm` must be > 0.
	 *  - If `--state-list` is NOT provided:
	 *      - `--Kinit`, `--Kpbwt`, `--pbwt-depth`, and `--pbwt-modulo-cm` must be > 0.
	 *
	 * **Depth constraints:**
	 *  - `--max-depth` must be ≥ 10.
	 *
	 * **Behavior on failure:**
	 *  - Any validation failure immediately calls `vrb.error(...)` which reports
	 *    the error and terminates execution.
	 *
	 * @note This method must be called immediately after parsing command-line options
	 *       to ensure all parameters are valid before any processing begins.
	 *
	 * @warning Many of these constraints are program-specific and must remain aligned
	 *          with the internal assumptions of the GLIMPSE phasing engine.
	 *
	 * @see parse_command_line() for the parsing stage before validation.
	 */
	void check_options();

	/**
	 * @brief Print all parsed and resolved GLIMPSE phase parameters in a human-readable format.
	 *
	 * This method outputs a structured summary of all important configuration
	 * parameters after parsing and validation have completed. It is primarily used
	 * for logging, debugging, and reproducibility, ensuring that all user-provided
	 * and default values are visible before computation begins.
	 *
	 * The printed information is grouped into logical sections, covering:
	 *
	 * **1. Input and reference configuration**
	 *  - **Input region**: Either from `--input-region` (if `InputFormat::BCF`)
	 *    or reported as "Given by binary reference panel".
	 *  - **Output region**: Similar logic to input region.
	 *  - **Sparse MAF**: Printed as a float, unless provided by the binary panel.
	 *  - **Recombination rates**: From `--map` if specified, otherwise a constant
	 *    1cM/Mb rate; for binary reference panels, this is predefined.
	 *  - **Ploidy**: If `--samples-file` is given, prints "Given by samples file";
	 *    otherwise, defaults to "Only diploid samples in region".
	 *  - **Keep monomorphic reference sites**: Printed as "YES"/"NO".
	 *
	 * **2. Model parameters**
	 *  - Burn-in and main iteration counts (`--burnin`, `--main`).
	 *  - Effective population size (`--ne`).
	 *  - Phase and imputation error rates, clamped to acceptable ranges.
	 *  - Minimum genotype likelihood value (`--min-gl`).
	 *
	 * **3. Selection parameters**
	 *  - K-init and K-pbwt selection sizes (with "not enabled" message if ≤ 0).
	 *  - PBWT depth and modulo-cM values, conditional on PBWT being enabled.
	 *  - State list path, or "No list provided".
	 *
	 * **4. Genotype calling (if not using GL input)**
	 *  - Calling model and indel handling mode.
	 *  - BAM/CRAM filtering thresholds: mapping quality, base quality, maximum depth.
	 *  - Flags for keeping QC-failed reads, orphan reads, duplicates, supplementary
	 *    alignments, and other BAM-specific behaviors.
	 *  - Additional sequencing data handling options:
	 *    pairing checks, orientation ignoring, Illumina-1.3+ flag.
	 *
	 * **5. Genotype likelihood input mode (if `--input-gl` is used)**
	 *  - Whether GL-indels are used (via PLs) or haplotype scaffolding.
	 *
	 * **6. Miscellaneous**
	 *  - Random seed value (`--seed`).
	 *  - Thread count (`--threads`).
	 *
	 * @note Many printed values are conditionally dependent on `input_fmt`
	 *       (`InputFormat::BCF` vs. binary reference panel) and on whether
	 *       GL inputs are used.
	 *
	 * @warning This function assumes that `check_options()` has already validated
	 *          all input parameters and that all required options are present.
	 *
	 * @see parse_command_line() for option parsing.
	 * @see check_options() for validation logic before printing.
	 */
	void verbose_options();

	/**
	 * @brief Prints detailed information about input/output files and formats used by the caller.
	 *
	 * This function outputs the key file paths, formats, and compression settings
	 * for the current run to the verbose logger (`vrb`). It uses a human-readable
	 * format for display, converting internal enum values to descriptive strings.
	 *
	 * The printed details include:
	 * - Input file (BAM/CRAM list, BAM/CRAM single file, or VCF/BCF likelihoods)
	 * - Reference sequence FASTA
	 * - Reference genotype file (either VCF or binary format)
	 * - Genetic map file (if provided)
	 * - Output file name, format, and compression method
	 * - Log file path (if provided)
	 *
	 * The output format and compression type are displayed using:
	 * - `fmt2string`: maps `OutputFormat` enum values to string names (VCF, BCF, BGEN)
	 * - `compr2string`: maps compression type enums to string names (NO, ZLIB, ZSTD)
	 * - For BGEN format, the number of bits used for storing probabilities is appended
	 *
	 * @note This method only reports file paths and configurations—it does not
	 *       perform any validation or file I/O.
	 *
	 * Example output:
	 * @code
	 * Files:
	 *   List BAM/CRAM        : [bam_list.txt]
	 *   Reference seq. fasta : [ref.fa]
	 *   Reference VCF        : [ref.vcf.gz]
	 *   Genetic Map          : [map.txt]
	 *   Output file          : [result.bgen]
	 *   Output format        : [BGEN format - 8 bits | ZSTD compression]
	 *   Output LOG           : [run.log]
	 * @endcode
	 */
	void verbose_files();

	//FILE I/O
	/**
	 * @brief Prints detailed information about the reference panel being used.
	 *
	 * This function logs key statistics about the reference panel, including:
	 * - Number of reference haplotypes
	 * - Total number of variant sites
	 * - Distribution of rare and common variants (with percentages)
	 * - Genomic input/output regions being processed
	 * - Sparse minor allele frequency (MAF) threshold
	 *
	 * The output is intended for verbose logging to help the user
	 * verify the reference panel and genomic regions being analyzed.
	 *
	 * @param ref_string A label or description of the reference panel 
	 *                   (e.g., `"1000 Genomes"` or `"Custom Panel"`).
	 *
	 * @note
	 * - Rare sites are defined based on a MAF threshold (implementation-specific).
	 * - Percentages for rare/common sites are computed relative to the total number of sites.
	 * - `H` contains reference panel statistics.
	 * - `V` contains variant/region information.
	 * - `vrb` handles formatted verbose logging.
	 *
	 * **Example Output:**
	 * @code
	 * > Custom reference panel [Nrh=5008] [L=1000000] [Lrare= 25000 (2.5%) - Lcommon= 975000 (97.5%)]
	 * > Input region         : [chr1:10000-50000]
	 * > Output region        : [chr1:15000-45000]
	 * > Sparse MAF           : [0.01]
	 * @endcode
	 */
	void print_ref_panel_info(const std::string ref_string);

	/**
	 * @brief Reads all input files and initializes internal data structures for imputation/phasing.
	 *
	 * This function performs the full initialization process, including:
	 * - Setting random seeds and thread configuration
	 * - Loading reference panel data (VCF/BCF or binary format)
	 * - Reading sample files and genotype likelihoods or BAMs
	 * - Initializing PBWT structures and HMM models
	 * - Loading genetic maps
	 * - Allocating buffers for imputation and phasing
	 * - Performing consistency checks and optional state-list loading
	 *
	 * The initialization is essential before running the imputation/phasing pipeline.
	 *
	 * @details
	 * **Main steps:**
	 * 1. **Seed and threading setup**
	 *    - Sets RNG seed from `--seed`
	 *    - Reads thread count from `--threads` and initializes worker threads if > 1
	 *    - Stores iteration counts for initialization, burn-in, and main phases
	 *
	 * 2. **Reference panel loading**
	 *    - If `input_fmt` is `BCF`, calls `buildCoordinates()` and uses `genotype_reader`
	 *      to load reference and target genotypes from VCF/BCF
	 *    - If `input_fmt` is binary, deserializes `H` and `V` using `boost::archive::binary_iarchive`
	 *    - Calls `print_ref_panel_info()` to log panel stats
	 *
	 * 3. **Sample and genotype likelihood handling**
	 *    - If `--samples-file` is given, reads sample IDs and ploidy
	 *    - If `--input-gl` is provided, loads target genotype likelihoods
	 *    - If `--bam-list` or `--bam-file` is provided, sets up mpileup processing and BAM reading
	 *
	 * 4. **Reference PBWT and rare-site handling**
	 *    - Calls `H.transposeRareRef()` to rearrange reference haplotypes
	 *
	 * 5. **Genetic map loading**
	 *    - Reads genetic map from `--map` if provided, else sets default uniform map
	 *    - Logs genomic span in base pairs and centimorgans
	 *
	 * 6. **Buffer allocation for HMMs**
	 *    - Allocates `HP0`, `HP1` for haplotype probabilities
	 *    - Allocates `HLC` for likelihoods, and `HMM` / `DMM` for imputation and phasing HMMs
	 *    - Creates `conditioning_set` objects for each thread
	 *
	 * 7. **PBWT allocation**
	 *    - Calls `H.allocatePBWT()` with depth and modulo parameters
	 *
	 * 8. **State-list loading**
	 *    - If `--state-list` is given, loads precomputed state lists
	 *
	 * 9. **Checksum update**
	 *    - If `--checkpoint-file-in` or `--checkpoint-file-out` is provided, updates CRC checksums
	 *
	 * @throws std::runtime_error if:
	 * - Thread count <= 0
	 * - Input files are missing or unreadable
	 * - Binary reference panel is empty or corrupted
	 * - No valid genotype likelihood or BAM input is provided
	 *
	 * @note
	 * - `H`, `G`, `V`, `M` are major data structures holding haplotypes, genotypes, variants, and metadata.
	 * - PBWT = Positional Burrows–Wheeler Transform.
	 * - This function must be called before any imputation/phasing stages.
	 *
	 * **Example:**
	 * @code
	 * caller c;
	 * c.read_files_and_initialise();
	 * @endcode
	 */
	void read_files_and_initialise();

	/**
	 * @brief Configures mpileup parameters and loads BAM/CRAM/FASTA inputs.
	 *
	 * This function initializes all necessary settings for the mpileup stage,
	 * including reference genome loading, read filtering options, quality thresholds,
	 * and the list of BAM/CRAM files to process.
	 *
	 * It reads settings from the `options` object, configures the `M` (mpileup state)
	 * structure, and performs input validation checks to ensure the run is correctly set up.
	 *
	 * ### Major Steps
	 * 1. **Prevent unintended reference sequence downloads**  
	 *    If the `--download-fasta-ref` option is not specified, the function sets  
	 *    environment variables `REF_CACHE` and `REF_PATH` to prevent `htslib`
	 *    from attempting to download sequence data, even when `-f` reference is provided.
	 *
	 * 2. **Reference genome loading**  
	 *    - If `--fasta` is specified, the reference genome file is loaded via `fai_load`.
	 *    - Requires a `.fai` index; errors if missing or unreadable.
	 *
	 * 3. **Read filtering flags**  
	 *    Sets `M.fflag` to filter unwanted reads:
	 *    - Always filters **unmapped** (`BAM_FUNMAP`) and **secondary** (`BAM_FSECONDARY`) reads.
	 *    - Filters **QC-fail** reads unless `--keep-failed-qc` is set.
	 *    - Filters **duplicates** unless `--keep-duplicates` is set.
	 *    - (Optionally) supplementary alignments can be filtered (commented out).
	 *
	 * 4. **Other filtering/validation parameters**  
	 *    - `M.keep_orphan`: keep orphan reads if `--keep-orphan-reads`.
	 *    - `M.check_orientation`: check read orientation unless `--ignore-orientation`.
	 *    - `M.check_proper_pair`: check proper pairing if `--check-proper-pairing`.
	 *    - Minimum mapping quality (`--mapq`), base quality (`--baseq`), and maximum depth (`--max-depth`)
	 *      are validated to ensure they are non-negative.
	 *    - Illumina 1.3+ quality encoding handled via `--illumina13+`.
	 *
	 * 5. **Call model restriction**  
	 *    - Only `"standard"` call model is currently supported.
	 *    - Any other value triggers an error.
	 *
	 * 6. **BAM/CRAM input handling**  
	 *    - If `--bam-list` is given:  
	 *      Reads one or two columns per line:  
	 *        - Column 1: BAM/CRAM filename (must be unique).  
	 *        - Column 2 (optional): Sample name (must be unique).  
	 *      If no sample name is provided, it is derived from the filename (without extension).
	 *    - If no list is provided, requires a single `--bam-file` with optional `--ind-name`.
	 *    - Errors if zero BAM files are given.
	 *
	 * ### Errors and validation
	 * The function will stop with an error if:
	 * - Reference FASTA cannot be loaded or lacks `.fai` index.
	 * - Invalid negative values for `--mapq`, `--baseq`, or `--max-depth`.
	 * - Duplicate BAM file paths or sample names in `--bam-list`.
	 * - No BAM/CRAM input files are provided.
	 *
	 * @note This function modifies the global `M` mpileup state object with all configured parameters.
	 *
	 * @warning `--download-fasta-ref` should be explicitly set if you intend to allow remote sequence fetching.
	 *
	 * @see fai_load
	 * @see BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
	 */
	void setup_mpileup();

	/**
	 * @brief Reads and processes BAM/CRAM files for genotype calling and coverage statistics.
	 *
	 * This function orchestrates reading input BAM/CRAM files, performing genotype calling
	 * via `genotype_bam_caller` instances (one per thread), and collecting sequencing depth statistics
	 * for each target sample. The results include per-sample coverage, missing data percentage,
	 * and observed versus expected counts from a Poisson distribution model.
	 *
	 * The function supports multi-threaded reading when more than one thread is specified.
	 * After processing, coverage statistics are written to a compressed output file
	 * (`*_stats_coverage.txt.gz`), and memory used for coverage data is released.
	 *
	 * ### Workflow:
	 * 1. Initialize BAM readers (one per thread).
	 * 2. Allocate statistics containers for coverage and depth counts.
	 * 3. Read BAM/CRAM files either in parallel (multi-threaded) or sequentially.
	 * 4. Delete BAM reader objects and clean up.
	 * 5. Output detailed per-sample coverage statistics to file.
	 * 6. Display average sequencing coverage in logs.
	 * 7. Release memory for statistics data.
	 *
	 * ### Coverage Statistics:
	 * For each target sample, the following are reported:
	 * - Mean coverage
	 * - Percentage of missing data (depth = 0)
	 * - Depth counts in ranges:
	 *   - Depth 0, 1, 2, 3, 4
	 *   - Depth 5–10, 11–30, and >30
	 * - Observed counts vs. expected counts from Poisson distribution with mean coverage.
	 *
	 * @note If mean coverage is `0` or extremely low, warnings are logged for possible BAM/CRAM file issues.
	 * @note Progress updates are displayed in the log if single-threaded processing is used.
	 *
	 * @warning This function assumes that `M.n_tar_samples`, `M.max_dp`, `H`, `G`, `V`, and `options`
	 *          are correctly initialized before calling.
	 *
	 * @thread_safety Uses pthreads for multi-threaded BAM reading. Statistics arrays are separate per sample.
	 *
	 * @file_output Writes coverage statistics to a gzipped text file named:
	 *              `<output>_stats_coverage.txt.gz`
	 *
	 * @exception May throw exceptions if:
	 *            - File writing fails
	 *            - BAM reading encounters unrecoverable errors
	 *
	 * @see genotype_bam_caller
	 *
	 * @return void
	 */
	void read_BAMs();
	
	/**
	 * @brief Main entry point for the phasing workflow.
	 *
	 * This function serves as the high-level orchestration method for running
	 * the complete phasing pipeline. It parses command-line arguments, validates
	 * user-provided options, initializes input data structures, executes the
	 * phasing algorithm, and writes the final results to output files.
	 *
	 * ### Workflow Steps:
	 * 1. **declare_options()**  
	 *    Registers all configurable command-line options.
	 *
	 * 2. **parse_command_line(args)**  
	 *    Parses the command-line arguments provided in `args` and stores
	 *    them in the internal options structure.
	 *
	 * 3. **check_options()**  
	 *    Validates the parsed options for logical consistency and required values.
	 *    Terminates with an error if invalid parameters are found.
	 *
	 * 4. **verbose_files()**  
	 *    Logs information about the input/output files that will be used.
	 *
	 * 5. **verbose_options()**  
	 *    Logs the full list of configured options for reproducibility.
	 *
	 * 6. **read_files_and_initialise()**  
	 *    Reads input files (reference panels, variant data, etc.) and initializes
	 *    the required internal data structures for phasing.
	 *
	 * 7. **phase_loop()**  
	 *    Executes the main iterative phasing loop until convergence or
	 *    until the configured number of iterations/stages is reached.
	 *
	 * 8. **write_files_and_finalise()**  
	 *    Writes the phased haplotypes, statistics, and any auxiliary output files.
	 *    Cleans up resources before program exit.
	 *
	 * ### Usage:
	 * This method is typically called once per program run, after preparing
	 * the command-line arguments vector.
	 *
	 * @param args Vector of command-line arguments, where `args[0]` is typically
	 *             the program name and subsequent entries are user-specified options.
	 *
	 * @pre `args` must contain valid arguments expected by `declare_options()`.
	 * @post Output files are written to disk as specified in the options.
	 *
	 * @exception May terminate execution if:
	 *            - Command-line parsing fails
	 *            - Required options are missing
	 *            - File I/O errors occur during initialization or finalization
	 *
	 * @note This function does not return a value; it controls the complete
	 *       phasing workflow from start to finish.
	 *
	 * @see declare_options()
	 * @see parse_command_line()
	 * @see check_options()
	 * @see read_files_and_initialise()
	 * @see phase_loop()
	 * @see write_files_and_finalise()
	 *
	 * @return void
	 */
	void phase(std::vector < std::string > &);

	/**
	 * @brief Finalizes the phasing process by writing output files and releasing resources.
	 *
	 * This method performs the last step of the phasing pipeline. It handles:
	 * - Destroying synchronization primitives if multi-threading was used.
	 * - Writing genotype data to disk in the requested format (VCF, BCF, or BGEN).
	 * - Displaying the output file path.
	 * - Reporting the total runtime of the program.
	 *
	 * ## Workflow:
	 * 1. **Thread cleanup**:
	 *    - If more than one thread was used, the function destroys the global `mutex_workers` mutex.
	 *
	 * 2. **Output file writing**:
	 *    - If the output format is **not** BGEN:
	 *      - Writes best-guess haplotypes using `genotype_writer::writeGenotypes()` to the specified VCF/BCF file.
	 *    - If the output format is **BGEN**:
	 *      - If BGEN support is compiled in (`__BGEN__` defined), writes the output using `writeGenotypesBgen()`.
	 *      - If BGEN is **not** supported, issues a warning, falls back to BCF output, and appends `.bcf` to the output file name.
	 *
	 * 3. **Output confirmation**:
	 *    - Prints the final output file path to the verbose log.
	 *
	 * 4. **Runtime reporting**:
	 *    - Prints the total elapsed execution time in both human-readable and raw seconds formats.
	 *
	 * @note
	 * - The function assumes that `options` contains all necessary parameters (e.g., `"output"`, `"threads"`, `"main"`, `"contigs-fai"`).
	 * - BGEN support must be explicitly enabled at compile time; otherwise, BCF will be used as a fallback.
	 * - If BGEN output is requested but not supported, no error will be thrown here; the user will be warned instead.
	 *
	 * @warning
	 * - Destroying the mutex without ensuring all worker threads have finished may lead to undefined behavior.
	 * - Ensure that the output directory exists and is writable before calling this function.
	 *
	 * @see genotype_writer
	 */
	void write_files_and_finalise();

	//REGION
	/**
	 * @brief Parses and validates genomic coordinate strings for input and output regions.
	 *
	 * This method reads the `input-region` and `output-region` options, parses them into
	 * chromosome IDs and start/stop positions, validates their formats and compatibility,
	 * and stores the resulting coordinates in the `V` structure.  
	 *
	 * Expected format for both input and output regions:
	 * @code
	 * chrX:Y-Z
	 * @endcode
	 * where:
	 * - `chrX`  = Chromosome identifier (e.g., `chr1`, `chr2`, `chrX`)
	 * - `Y`     = Start position (positive integer)
	 * - `Z`     = Stop position (positive integer greater than Y)
	 *
	 * ### Processing steps
	 * 1. Retrieves `input-region` and `output-region` strings from the program options.
	 * 2. Splits each region into chromosome ID and position range using `":"` as delimiter.
	 * 3. Validates that both have exactly two components (chromosome, range).
	 * 4. Checks that both input and output regions have the same chromosome ID.
	 * 5. Splits position ranges into start and stop coordinates using `"-"` as delimiter.
	 * 6. Validates numerical correctness of coordinates (start < stop, non-negative).
	 * 7. Ensures input region fully contains the output region.
	 * 8. Stores parsed and validated values in:
	 *    - `V.chrid`          → Chromosome ID
	 *    - `V.input_start`    → Input start coordinate
	 *    - `V.input_stop`     → Input stop coordinate
	 *    - `V.output_start`   → Output start coordinate
	 *    - `V.output_stop`    → Output stop coordinate
	 *    - `V.input_gregion`  → String representation of input region
	 *    - `V.output_gregion` → String representation of output region
	 * 9. Logs the parsed regions for user confirmation.
	 *
	 * @throws std::runtime_error if:
	 * - The region format is invalid.
	 * - Chromosome IDs differ between input and output.
	 * - Start/stop coordinates are illogical (e.g., start >= stop).
	 * - Input and output ranges are incompatible.
	 * - Any coordinate is negative.
	 *
	 * @note This function assumes that:
	 * - `options` contains valid entries for `"input-region"` and `"output-region"`.
	 * - `stb.split()` returns the number of tokens obtained.
	 * - `vrb` is a logging/error-reporting utility.
	 *
	 * @see stb::split(), stb::str(), vrb::error(), vrb::bullet()
	 */
	void buildCoordinates();

	//CHECKPOINTING
	/**
	 * @brief Writes the current state of the caller to a checkpoint file.
	 *
	 * This function creates a checkpoint file containing the current state of the program's
	 * execution so that it can be resumed later without starting from scratch.
	 * The checkpoint includes various runtime parameters, configuration options, and
	 * serialized data from the `G` object. The output file is first written to a temporary
	 * file to avoid partial writes in case of failure, and then atomically renamed to the
	 * final checkpoint file name.
	 *
	 * @details
	 * Steps performed:
	 * 1. Checks whether the `checkpoint-file-out` option is set in `options`.
	 * 2. If set, constructs the checkpoint filename and a corresponding temporary filename.
	 * 3. Opens the temporary file in binary mode for output.
	 * 4. Uses Boost's `binary_oarchive` to serialize:
	 *    - CRC checksum value (`crc.get_value()`).
	 *    - Current stage of execution (`current_stage`).
	 *    - Current iteration (`current_iteration`).
	 *    - Number of iterations per burn-in stage.
	 *    - Several configuration parameters from `options` (e.g., `ne`, `min-gl`, `err-imp`, etc.).
	 *    - Checkpoint data from the genotype container `G` via `G.serialize_checkpoint_data()`.
	 * 5. Renames the temporary file to the final checkpoint filename using `std::filesystem::rename()`.
	 * 6. Logs the completion time of the checkpoint operation.
	 *
	 * @note This function will do nothing if the `checkpoint-file-out` option is not specified.
	 *
	 * @pre The `options` map must contain valid entries for all serialized parameters if the checkpoint option is set.
	 * @pre The object `G` must implement the method `serialize_checkpoint_data(boost::archive::binary_oarchive&)`.
	 *
	 * @post A valid checkpoint file will be created at the location specified by `checkpoint-file-out`.
	 *
	 * @warning If the program crashes before `std::filesystem::rename()` is called,
	 *          the temporary checkpoint file (`.tmp`) may remain on disk.
	 *
	 * @exception May throw exceptions from:
	 * - `std::ofstream` if file writing fails.
	 * - `boost::archive::binary_oarchive` if serialization fails.
	 * - `std::filesystem::rename` if file renaming fails.
	 *
	 * @see G::serialize_checkpoint_data()
	 *
	 * @code
	 * caller myCaller;
	 * myCaller.write_checkpoint(); // Creates a checkpoint if option is set
	 * @endcode
	 */
	void write_checkpoint();

	/**
	 * @brief Reads and restores program state from a checkpoint file, if specified.
	 *
	 * This method checks if a checkpoint input file (`checkpoint-file-in`) is provided
	 * via program options. If so, it attempts to restore the computation state from
	 * the file, ensuring data integrity and parameter consistency.
	 *
	 * The checkpoint file stores:
	 * - Input data CRC for validation
	 * - Current stage and iteration
	 * - Burn-in iterations count
	 * - Various run parameters (e.g., `ne`, `min-gl`, `err-imp`, etc.)
	 * - Serialized data from object `G`
	 *
	 * @throws vrb.error if:
	 *  - The CRC in the checkpoint does not match the current input data
	 *  - The checkpoint indicates more iterations than the current run is configured for
	 *  - The burn-in iteration counts do not match between the checkpoint and current run
	 *  - Any stored parameter does not match the current run configuration
	 *
	 * @note Uses Boost serialization (`boost::archive::binary_iarchive`) for data I/O.
	 * @note Requires that `confirm_checkpoint_param<T>` correctly verifies each parameter.
	 *
	 * @see write_checkpoint()
	 */
	void read_checkpoint_if_available();
	
	/**
	 * @brief Confirms that a parameter stored in a checkpoint file matches the current run configuration.
	 *
	 * This function reads a parameter value of type `T` from a checkpoint archive and compares it
	 * against the corresponding parameter in the current run's options. If the values do not match,
	 * the function will throw an error with a descriptive message, indicating that the run must
	 * be configured with the same parameter value as in the checkpoint in order to proceed.
	 *
	 * @tparam T The data type of the parameter being checked (e.g., int, float).
	 * @tparam Archive The type of the archive used for reading the checkpoint data 
	 *                 (e.g., `boost::archive::binary_iarchive`).
	 *
	 * @param ar The input archive stream from which the checkpoint parameter value will be read.
	 * @param param_name The name of the parameter (as stored in `options`) to check.
	 *
	 * @throws std::runtime_error If the parameter value in the checkpoint does not match
	 *         the value in the current run's options.
	 *
	 * @note This function is typically used to ensure compatibility between the current run's settings
	 *       and the parameters used in generating the checkpoint file. If a mismatch is detected,
	 *       it prevents loading the checkpoint to avoid inconsistent results.
	 *
	 * @see read_checkpoint_if_available()
	 */
	template <typename T, class Archive>
	void confirm_checkpoint_param(Archive &ar, std::string param_name) {
		T checkpoint_value;
		ar >> checkpoint_value;
		T current_value = options[param_name].as<T>();
		if (checkpoint_value != current_value) {
			std::stringstream err_str;
			err_str<<"Checkpoint value was run with "<<param_name<<" set to "<<checkpoint_value<<" and "
			"this run has "<<param_name<<" set to "<<current_value<<".  You must set "<<param_name<<" "
			"to "<<checkpoint_value<<" in order to use this checkpoint file.";
			vrb.error(err_str.str());
		}
	}
};


#endif


