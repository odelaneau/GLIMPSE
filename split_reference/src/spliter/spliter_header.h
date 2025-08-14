/*******************************************************************************
 * Copyright (C) 2022-2023 Simone Rubinacci
 * Copyright (C) 2022-2023 Olivier Delaneau
 *
 * MIT Licence
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#ifndef _SPLITER_H
#define _SPLITER_H

#include "otools.h"

#include "variant_map.h"
#include "ref_genotype_reader.h"
#include "ref_haplotype_set.h"

#define STAGE_INIT	0
#define STAGE_BURN	1
#define STAGE_MAIN	2


/**
 * @class spliter
 * @brief Splits a genomic reference dataset into user-defined regions.
 *
 * The `spliter` class implements the workflow for reading user parameters,
 * parsing genomic coordinates, processing reference data, and writing
 * the resulting subset to disk. This module is part of the GLIMPSE toolkit.
 *
 * ## Workflow Overview
 *
 * The main execution flow is initiated by `spliter::phase()`:
 *
 * 1. **declare_options()**
 *    - Registers command-line options with default values and descriptions.
 *
 * 2. **parse_command_line(args)**
 *    - Parses command-line arguments into the `options` map.
 *
 * 3. **check_options()**
 *    - Ensures all required options are present and values are valid.
 *
 * 4. **verbose_files()**
 *    - Logs input/output file paths and required parameters.
 *
 * 5. **verbose_options()**
 *    - Logs additional run parameters (e.g., seed, MAF thresholds, threads).
 *
 * 6. **read_files_and_initialise()**
 *    - Reads the reference VCF and prepares in-memory data structures.
 *    - Calls `buildCoordinates()` to parse and validate genomic regions.
 *
 * 7. **write_files_and_finalise()**
 *    - Writes processed binary output and optional logs.
 *    - Destroys any allocated threading resources.
 *    - Prints total runtime.
 *
 * ---
 *
 * ### Detailed Function Descriptions
 *
 * #### `phase(std::vector<std::string>& args)`
 * @brief Orchestrates the complete `spliter` pipeline from argument parsing to finalization.
 *
 * #### `verbose_files()`
 * @brief Logs file-related parameters and checks required input/output regions.
 *
 * #### `verbose_options()`
 * @brief Logs processing options (MAF, monomorphic site handling, seed, threads).
 *
 * #### `buildCoordinates()`
 * @brief Parses and validates the `input-region` and `output-region` strings
 * into chromosome IDs and genomic coordinate ranges.
 * Ensures compatibility between input and output.
 *
 * #### `write_files_and_finalise()`
 * @brief Cleans up multi-threading resources, writes final logs, and reports runtime.
 *
 * ---
 *
 * ## Example
 * ```bash
 * GLIMPSE_split_reference \
 *   --reference ref.vcf.gz \
 *   --input-region chr20:1000000-2000000 \
 *   --output-region chr20:1100000-1900000 \
 *   --output subset.bcf \
 *   --threads 4
 * ```
 */
class spliter {
public:
	// =========================
    // COMMAND LINE OPTIONS
    // =========================

    /**
     * @brief Stores the description of available command-line options.
     *
     * Contains human-readable descriptions, defaults, and parameter types for all
     * supported command-line arguments.
     */
    bpo::options_description descriptions;

    /**
     * @brief Parsed command-line arguments and their values.
     *
     * Populated via `boost::program_options::store()` and accessed throughout
     * the splitting workflow.
     */
    bpo::variables_map options;

    /**
     * @brief Chromosome identifiers for parsed genomic regions.
     *
     * Each element corresponds to an entry in the `input_start`/`input_stop`
     * and `output_start`/`output_stop` vectors.
     */
    std::vector<std::string> chrid;

    /**
     * @brief Start coordinates (1-based) for each input genomic region.
     */
    std::vector<int> input_start;

    /**
     * @brief End coordinates (1-based, inclusive) for each input genomic region.
     */
    std::vector<int> input_stop;

    /**
     * @brief Input genomic regions in `chr:start-stop` format.
     */
    std::vector<std::string> input_gregion;

    /**
     * @brief Start coordinates (1-based) for each output genomic region.
     */
    std::vector<int> output_start;

    /**
     * @brief End coordinates (1-based, inclusive) for each output genomic region.
     */
    std::vector<int> output_stop;

    /**
     * @brief Output genomic regions in `chr:start-stop` format.
     */
    std::vector<std::string> output_gregion;

    // =========================
    // MULTI-THREADING
    // =========================

    /**
     * @brief Number of worker threads available for processing.
     */
    int i_workers;

    /**
     * @brief Number of processing jobs to be handled by threads.
     */
    int i_jobs;

    /**
     * @brief Identifiers for worker threads.
     */
    std::vector<pthread_t> id_workers;

    /**
     * @brief Mutex to synchronize access between worker threads.
     */
    pthread_mutex_t mutex_workers;

	//CONSTRUCTOR
	/**
	 * @brief Constructs a new spliter object with default values.
	 *
	 * Initializes all member variables to their default states:
	 * - `input_start` = 0
	 * - `input_stop` = 0
	 * - `i_workers` = 0
	 * - `i_jobs` = 0
	 *
	 * This constructor does not perform any file I/O or job scheduling; it
	 * simply sets up an empty spliter instance ready for configuration.
	 */
	spliter();

	/**
	 * @brief Destroys the spliter object.
	 *
	 * This destructor currently performs no explicit cleanup,
	 * as all resources are managed automatically.
	 * 
	 * @note If dynamic memory or file handles are added to the class
	 * in the future, this destructor should be updated to release them.
	 */
	~spliter();

	//METHODS
	/**
	 * @brief Performs haplotype phasing for a single individual in a specific job and worker thread.
	 *
	 * This method executes the phasing process for the given `id_job` using the resources
	 * assigned to the `id_worker`. The workflow depends on the current processing stage
	 * and the ploidy of the target sample.
	 *
	 * **Main steps:**
	 * 1. **Conditioning set selection** — Selects the conditioning haplotypes for the target individual.
	 * 2. **Haplotype likelihood initialization**  
	 *    - In `STAGE_INIT`, initializes haplotype likelihoods.  
	 *    - In other stages:
	 *      - For diploid or polyploid (`ploidy > 1`): builds haplotype likelihoods.  
	 *      - For haploid: initializes haplotype likelihoods directly.
	 * 3. **HMM posterior computation and sampling**  
	 *    - Computes posterior probabilities using the HMM.  
	 *    - Samples haplotype H0.  
	 *    - If polyploid, also recomputes likelihoods for H1, computes posteriors, samples H1,
	 *      and performs haplotype re-phasing via the diploid model matcher (DMM).
	 * 4. **Results storage** — In `STAGE_MAIN`, stores genotype posteriors and haplotypes
	 *    for later output.
	 * 5. **Statistics update** — Updates state-space size and polymorphism statistics.
	 *
	 * Thread safety:  
	 * Uses `pthread_mutex_lock` / `pthread_mutex_unlock` to ensure safe access to
	 * shared statistics when running with multiple threads.
	 *
	 * @param id_worker Index of the worker thread handling this job.
	 * @param id_job Index of the job (target individual/sample) being processed.
	 *
	 * @note This function assumes that all relevant data structures (`G`, `COND`, `HLC`, `HMM`,
	 * `DMM`, `HP0`, `HP1`, etc.) have been pre-allocated and initialized for each worker.
	 * @see caller::phase_all()
	 */
	void phase_individual(const int, const int);

	/**
	 * @brief Executes one full phasing iteration over all individuals.
	 *
	 * This function controls the main iteration loop for the phasing process.
	 * The behavior depends on the current processing stage (`current_stage`),
	 * and can run either in single-threaded or multi-threaded mode.
	 *
	 * **Stage-dependent operations:**
	 * - **STAGE_INIT**:
	 *   1. Logs initialization message.
	 *   2. Initializes rare target haplotypes (`H.initRareTar()`).
	 *   3. Performs rare haplotype selection (`H.performSelection_RARE_INIT_GL()`), non-parallel.
	 * - **Other stages**:
	 *   1. Logs stage and iteration number.
	 *   2. Updates haplotypes (`H.updateHaplotypes()`).
	 *   3. Transposes rare targets (`H.transposeRareTar()`).
	 *   4. Matches haplotypes using compressed PBWT (`H.matchHapsFromCompressedPBWTSmall()`),
	 *      with a boolean flag for `STAGE_MAIN` mode.
	 *
	 * **Phasing execution:**
	 * - Resets timers, counters (`i_workers`, `i_jobs`), and statistics (`statH`, `statC`).
	 * - Computes progress bar step size based on number of individuals (`G.n_ind`).
	 * - If multi-threaded (`n_thread > 1`):
	 *   - Spawns worker threads (`pthread_create`) to run `phase_callback`.
	 *   - Waits for all threads to complete (`pthread_join`).
	 * - Otherwise:
	 *   - Iterates over all individuals sequentially, calling `phase_individual()`
	 *     for each job.
	 *   - Updates progress bar after each individual.
	 *
	 * **Post-processing:**
	 * - Logs mean HMM state size and polymorphism percentage, with elapsed time.
	 * - Writes a checkpoint of current state (`write_checkpoint()`).
	 * - If `STAGE_INIT`, releases memory from `H.init_states`.
	 *
	 * @note This function orchestrates the phasing of all individuals
	 *       and should be called repeatedly for each iteration in the
	 *       phasing schedule.
	 *
	 * @see caller::phase_individual()
	 * @see phase_callback()
	 * @see write_checkpoint()
	 *
	 * @threadsafe Uses pthreads for parallel execution, with shared
	 *             statistics collected across threads.
	 */
	void phase_iteration();

	/**
	 * @brief Executes the complete phasing loop across all processing stages.
	 *
	 * This function is the top-level driver for the phasing process. It iterates
	 * through all phasing stages (`STAGE_INIT` to `STAGE_MAIN`), calling
	 * `phase_iteration()` for each iteration, until all scheduled iterations are complete.
	 *
	 * **Workflow:**
	 * 1. **Initialize iteration counters**
	 *    - Sets `current_stage` to `STAGE_INIT`.
	 *    - Sets `current_iteration` to -1 so the first `increment_iteration()`
	 *      call starts from iteration 0 (unless resuming from a checkpoint).
	 * 2. **Checkpoint recovery**
	 *    - Calls `read_checkpoint_if_available()` to restore saved iteration
	 *      state if a checkpoint file exists.
	 * 3. **Main phasing loop**
	 *    - Calls `increment_iteration()` to start the first iteration.
	 *    - While `current_stage <= STAGE_MAIN`:
	 *      - Runs one full iteration via `phase_iteration()`.
	 *      - Advances iteration counters with `increment_iteration()`.
	 * 4. **Finalization**
	 *    - After completing all stages, calls
	 *      `sortAndNormAndInferGenotype()` on each sample (`G.vecG[i]`)
	 *      to finalize genotype calls.
	 *
	 * @note This function handles stage progression and checkpoint recovery,
	 *       but delegates actual per-iteration work to `phase_iteration()`.
	 *
	 * @see caller::phase_iteration()
	 * @see caller::increment_iteration()
	 * @see caller::read_checkpoint_if_available()
	 * @see caller::sortAndNormAndInferGenotype()
	 */
	void phase_loop();

	//PARAMETERS
	/**
	 * @brief Declares all command-line options for the spliter module.
	 *
	 * This function defines and groups all supported CLI parameters using
	 * Boost Program Options (`bpo::options_description`), then registers them
	 * in the `descriptions` container for later parsing.
	 *
	 * **Option groups:**
	 *
	 * - **Basic parameters** (`opt_base`):
	 *   - `--help` — Prints the help message.
	 *   - `--seed` *(int, default: 15052011)* — Seed for the random number generator.
	 *   - `--threads` *(int, default: 1)* — Number of threads to use.
	 *
	 * - **Input parameters** (`opt_input`):
	 *   - `--reference, -R` *(string)* — Haplotype reference panel in VCF/BCF format.
	 *   - `--map, -M` *(string)* — Genetic map file.
	 *   - `--input-region` *(string)* — Imputation region with buffers.
	 *   - `--output-region` *(string)* — Imputation region without buffers.
	 *   - `--sparse-maf` *(float, default: 0.001)* — Rare variant threshold (expert setting).
	 *   - `--keep-monomorphic-ref-sites` — Keeps monomorphic markers in the reference panel (expert setting, disabled by default).
	 *
	 * - **Output parameters** (`opt_output`):
	 *   - `--output, -O` *(string)* — Output file prefix; region and extension are appended automatically.
	 *   - `--log` *(string)* — Log file path.
	 *
	 * @note This function only *declares* options; parsing and validation
	 *       are handled elsewhere in the program.
	 *
	 * @see spliter::parse_command_line()
	 * @see boost::program_options::options_description
	 */
	void declare_options();

	/**
	 * @brief Parses and processes the command-line arguments for the spliter module.
	 *
	 * This function takes a list of command-line arguments, parses them using
	 * Boost Program Options (`bpo::command_line_parser`), and stores the results
	 * in the `options` variable. It also initializes logging, prints program
	 * metadata, and handles the `--help` option.
	 *
	 * **Workflow:**
	 * 1. **Parsing and notification**
	 *    - Runs the Boost command-line parser with the previously declared
	 *      option descriptions (`descriptions`).
	 *    - Calls `bpo::store()` to save parsed values into `options`.
	 *    - Calls `bpo::notify()` to trigger Boost's type conversion and
	 *      validation.
	 *    - If parsing fails, prints the error and exits.
	 * 2. **Logging setup**
	 *    - If `--log` is provided, attempts to open the specified log file
	 *      via `vrb.open_log()`. If this fails, prints an error and exits.
	 * 3. **Program banner**
	 *    - Prints program title, authors, contact information, version,
	 *      commit ID/date, and relevant citations.
	 *    - Displays the run date from the timer/clock utility (`tac.date()`).
	 * 4. **Help handling**
	 *    - If `--help` is present, prints all declared option descriptions
	 *      and exits.
	 *
	 * @param args Vector of command-line argument strings (excluding `argv[0]`).
	 *
	 * @note This function must be called after `spliter::declare_options()`
	 *       so that the `descriptions` object is populated before parsing.
	 * @warning The program exits immediately on parse errors, log file
	 *          creation failures, or after displaying help.
	 *
	 * @see spliter::declare_options()
	 * @see boost::program_options::command_line_parser
	 * @see vrb.open_log()
	 */
	void parse_command_line(std::vector < std::string > &);

	/**
	 * @brief Validates the command-line options for the spliter module.
	 *
	 * This function checks that all required CLI parameters have been provided
	 * and that their values are within valid ranges. If any check fails, the
	 * program terminates with an error message via `vrb.error()`.
	 *
	 * **Validation rules:**
	 * - **Required parameters:**
	 *   - `--input-region` **and** `--output-region` must both be specified.
	 *   - `--output` must be provided.
	 *   - `--reference` must be provided.
	 * - **Value checks:**
	 *   - `--seed` (if provided) must be a non-negative integer.
	 *   - `--threads` must be a strictly positive integer.
	 *   - `--sparse-maf` must be in the range `[0, 0.5)`. Recommended range:
	 *     `[0.001%, 1%]`, with 0.1% as the suggested default.
	 *
	 * @note This function assumes that `spliter::parse_command_line()` has
	 *       already populated the `options` variable.
	 * @warning Any violation results in immediate program termination.
	 *
	 * @see spliter::declare_options()
	 * @see spliter::parse_command_line()
	 */
	void check_options();

	/**
	 * @brief Prints the parsed command-line options in a human-readable format.
	 *
	 * This function outputs the main runtime parameters to the verbose logger (`vrb`)
	 * so that the user can confirm the settings before execution begins.
	 *
	 * **Displayed parameters:**
	 * - **Sparse MAF** — Percentage value of `--sparse-maf` (`float × 100`).
	 * - **Keep monomorphic reference sites** — `YES` if `--keep-monomorphic-ref-sites`
	 *   was provided, otherwise `NO`.
	 * - **Seed** — Random number generator seed (`--seed` value).
	 * - **#Threads** — Number of threads to use (`--threads` value).
	 *
	 * @note This is intended for user feedback and reproducibility logging.
	 *       It assumes that `spliter::check_options()` has already verified
	 *       that all relevant options are valid.
	 *
	 * @see spliter::declare_options()
	 * @see spliter::parse_command_line()
	 * @see spliter::check_options()
	 */
	void verbose_options();

	/**
	 * @brief Prints the list of input and output files used by the spliter module.
	 *
	 * This function logs all relevant file paths to the verbose logger (`vrb`)
	 * so the user can verify them before execution. It also performs a minimal
	 * check to ensure required file-related options are present.
	 *
	 * **Displayed files:**
	 * - **Input region** (`--input-region`)
	 * - **Output region** (`--output-region`)
	 * - **Reference VCF** (`--reference`)
	 * - **Output binary** (`--output`)
	 * - **Genetic map** (`--map`) — optional, only shown if provided
	 * - **Log file** (`--log`) — optional, only shown if provided
	 *
	 * **Validation:**
	 * - Both `--input-region` and `--output-region` must be present;
	 *   otherwise, the function terminates with an error via `vrb.error()`.
	 *
	 * @note This is intended for user feedback and reproducibility purposes.
	 *       It assumes that `spliter::check_options()` has already verified
	 *       all file-related options.
	 *
	 * @see spliter::declare_options()
	 * @see spliter::parse_command_line()
	 * @see spliter::check_options()
	 * @see spliter::verbose_options()
	 */
	void verbose_files();

	//FILE I/O
	/**
	 * @brief Reads all required input files, prepares data structures, and splits the reference panel.
	 *
	 * This function performs the full initialization and data loading workflow for the
	 * `spliter` module. It reads the reference panel and optional genetic map, prepares
	 * per-region variant and haplotype data, and writes the processed regions into
	 * serialized GLIMPSE2 binary files.
	 *
	 * **Workflow:**
	 * 1. **Initialization**
	 *    - Prints module header.
	 *    - Sets random number generator seed (`rng.setSeed()`).
	 *    - If using multiple threads, initializes worker/job counters,
	 *      resizes `id_workers`, and sets up `mutex_workers`.
	 *    - Calls `buildCoordinates()` to prepare genomic region coordinates.
	 *
	 * 2. **Genetic map loading**
	 *    - If `--map` is provided, loads the genetic map via `gmap_reader::readGeneticMapFile()`.
	 *
	 * 3. **Per-region processing** *(loop over `input_gregion`)*
	 *    - Creates new `variant_map` (`V`) and `ref_haplotype_set` (`H`) objects.
	 *    - Assigns chromosome ID, coordinate ranges, and region strings to `V`.
	 *    - Generates an output-safe region string (`reg_out`) by replacing `:` and `-` with `_`.
	 *    - Instantiates a `ref_genotype_reader` with MAF threshold and monomorphic site settings.
	 *    - Reads the reference panel (`readRefPanel()`).
	 *    - Assigns the genetic map to `V` (from file if loaded, otherwise default).
	 *    - Logs region length in base pairs and centiMorgans.
	 *    - Builds a sparse PBWT index (`H.build_sparsePBWT()`).
	 *
	 * 4. **Output**
	 *    - Opens a binary output stream with filename:
	 *      `"<output_prefix>_<reg_out>.bin"`.
	 *    - Serializes `H` and `V` using `boost::archive::binary_oarchive`.
	 *    - Logs file writing completion time.
	 *
	 * @note The output `.bin` files contain both the reference haplotypes (`H`)
	 *       and variant map (`V`) for each region.
	 * @warning Assumes that all CLI options have been validated with `check_options()`.
	 *
	 * @see spliter::buildCoordinates()
	 * @see ref_genotype_reader
	 * @see variant_map
	 * @see ref_haplotype_set
	 */
	void read_files_and_initialise();

	/**
	 * @brief Configure and initialize parameters for HTSlib mpileup processing.
	 *
	 * This function prepares the reference genome, read filtering flags, quality thresholds, 
	 * and BAM file list before performing mpileup operations. It prevents unintended remote 
	 * reference downloads, validates user input, and ensures all required files are accessible.
	 *
	 * Dataflow:
	 * 1. **Reference Download Control**  
	 *    - If `--download-fasta-ref` is not set, disables HTSlib remote reference fetching
	 *      by setting environment variables `REF_CACHE` and `REF_PATH`.
	 *
	 * 2. **Reference Genome Loading**  
	 *    - If `--fasta` is provided, loads the reference FASTA file and its `.fai` index using `fai_load()`.  
	 *    - Terminates with an error if the reference cannot be loaded.
	 *
	 * 3. **BAM Read Filtering**  
	 *    - Initializes BAM filtering flags to exclude unmapped and secondary alignments.  
	 *    - Optionally excludes failed QC reads and duplicates unless explicitly kept via options.
	 *    - Other read filtering and pairing/orientation checks are also configured here.
	 *
	 * 4. **Quality and Depth Thresholds**  
	 *    - Reads minimum mapping quality (`--mapq`), minimum base quality (`--baseq`), 
	 *      and maximum depth (`--max-depth`) from options.  
	 *    - Validates that these values are non-negative.
	 *
	 * 5. **Calling Model Validation**  
	 *    - Checks `--call-model` option, currently only allows `"standard"`.
	 *
	 * 6. **BAM Input List Handling**  
	 *    - If `--bam-list` is provided:  
	 *      - Reads a list of BAM filenames (and optional sample names) from a file.  
	 *      - Ensures no duplicate BAM files or sample names are present.  
	 *      - Populates `M.bam_fnames` and `M.tar_sample_names` vectors.  
	 *    - If `--bam-list` is not provided:  
	 *      - Expects `--bam-file` to specify a single BAM file, and optionally `--ind-name` for the sample name.
	 *
	 * 7. **Final Checks**  
	 *    - Verifies that at least one BAM file is provided before proceeding.
	 *
	 * @throws runtime_error if required files or parameters are missing, invalid, or duplicated.
	 */
	void setup_mpileup();

	/**
	 * @brief Reads BAM/CRAM files, calls mpileup, and generates coverage statistics per sample.
	 *
	 * This function performs parallel or single-threaded reading of BAM/CRAM input files
	 * to compute coverage metrics, store per-sample depth distributions, and generate
	 * summary statistics for downstream processing.
	 *
	 * ## Workflow:
	 * 1. **Initialization**
	 *    - Starts a timer and prints a log message indicating BAM file reading.
	 *    - Allocates one `genotype_bam_caller` instance per thread to handle mpileup calls.
	 *    - Initializes coverage (`G.stats.cov_ind`) and depth count (`G.stats.depth_count`) containers.
	 *
	 * 2. **Parallel or Single-threaded Reading**
	 *    - If multiple threads are requested (`--threads`), uses `pthread_create()` to launch
	 *      `read_BAMs_callback` workers.
	 *    - If single-threaded, iterates through all target samples, calling
	 *      `READER_BAM[0]->call_mpileup()` for each.
	 *    - Updates a progress bar during single-threaded execution.
	 *
	 * 3. **Cleanup**
	 *    - Deletes all `genotype_bam_caller` objects and clears `READER_BAM` vector.
	 *
	 * 4. **Coverage Statistics Generation**
	 *    - Creates an output filename `<output>_stats_coverage.txt.gz` based on the user-provided
	 *      `--output` path.
	 *    - Writes a header describing coverage statistics for each sample.
	 *    - For each target sample:
	 *        - Computes mean coverage.
	 *        - Issues warnings for unmapped or high-missingness samples.
	 *        - Groups depth counts into bins: 0–4, 5–10, 11–30, and 31+.
	 *        - Compares observed counts with expected counts from a Poisson model.
	 *        - Writes results in tab-separated format to the output file.
	 *
	 * 5. **Finalization**
	 *    - Logs location of coverage statistics file and average sequencing coverage.
	 *    - Clears and deallocates coverage and depth count vectors.
	 *
	 * @note This function assumes that `setup_mpileup()` has already been called
	 *       to configure `M` (mpileup parameters) and BAM file lists.
	 *
	 * @throws runtime_error if there are issues with BAM file reading, allocation, or coverage computation.
	 */
	void read_BAMs();

	/**
	 * @brief Executes the main workflow of the spliter module.
	 *
	 * This function coordinates the complete processing pipeline for the `spliter` tool.
	 * It parses user-provided arguments, validates settings, loads input data, and
	 * produces the required output files.
	 *
	 * ## Workflow:
	 * 1. **Declare Options**
	 *    - Calls `declare_options()` to register all command-line options
	 *      supported by this module.
	 *
	 * 2. **Parse Command Line**
	 *    - Parses `args` using `parse_command_line(args)` and stores
	 *      values into the internal options structure.
	 *
	 * 3. **Check Options**
	 *    - Validates option values and required parameters via `check_options()`.
	 *
	 * 4. **Verbose Logging**
	 *    - Calls `verbose_files()` to print details of input/output file paths.
	 *    - Calls `verbose_options()` to print all active parameter values.
	 *
	 * 5. **Main Processing**
	 *    - `read_files_and_initialise()` loads input data, prepares internal
	 *      structures, and initializes any state needed for processing.
	 *    - `write_files_and_finalise()` writes processed results to output files
	 *      and releases allocated resources.
	 *
	 * @param args Vector of command-line arguments for the `spliter` module.
	 *
	 * @note This is the top-level entry point for the `spliter` workflow.
	 *       All helper functions must be implemented and called in the correct order.
	 */
	void phase(std::vector < std::string > &);

	/**
	 * @brief Finalizes the `spliter` processing and writes completion logs.
	 *
	 * This function performs the final cleanup and reporting steps after
	 * the main processing has completed.
	 *
	 * ## Actions:
	 * 1. **Multi-threading Cleanup**
	 *    - If more than one thread was used (`threads` > 1), the
	 *      `mutex_workers` mutex is destroyed to free system resources.
	 *
	 * 2. **Runtime Reporting**
	 *    - Calculates the total elapsed time since the start of execution
	 *      using `tac.abs_time()` and logs it via `vrb.bullet()`.
	 *
	 * @note This function does not write any output data itself; all
	 *       data writing should occur before this point. It focuses on
	 *       cleanup and performance reporting.
	 */
	void write_files_and_finalise();

	//REGION
	/**
	 * @brief Parses and validates user-specified genomic regions.
	 *
	 * This function reads the `input-region` and `output-region` parameters from
	 * the program options, parses them into chromosome identifiers and genomic
	 * coordinate ranges, and validates their compatibility.
	 *
	 * ## Expected format:
	 * - Both input and output regions must be in the form:
	 *   @code
	 *   chrID:start-end
	 *   @endcode
	 *   e.g., `chr1:100000-200000`
	 *
	 * ## Actions:
	 * 1. **Option Reading**
	 *    - Reads `input-region` and `output-region` values if both are provided.
	 *
	 * 2. **String Parsing**
	 *    - Splits each region string into chromosome ID and coordinate range.
	 *    - Further splits coordinate range into `start` and `stop` integers.
	 *
	 * 3. **Validation Checks**
	 *    - Chromosome IDs must match between input and output regions.
	 *    - Coordinates must satisfy:
	 *      - `start < stop`
	 *      - Non-negative positions
	 *      - Input region fully covers output region
	 *
	 * 4. **Data Storage**
	 *    - Appends parsed values to:
	 *      - `chrid`
	 *      - `input_start`, `input_stop`
	 *      - `output_start`, `output_stop`
	 *      - `input_gregion`, `output_gregion`
	 *
	 * 5. **Logging**
	 *    - Prints parsed and validated regions to verbose output.
	 *
	 * @throws vrb.error if format validation or coordinate compatibility fails.
	 */
	void buildCoordinates();
};


#endif


