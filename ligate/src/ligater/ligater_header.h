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

#ifndef _LIGATER_H
#define _LIGATER_H

#include "otools.h"

/**
 * @class ligater
 * @brief Class to ligate multiple phased VCF/BCF chunks into chromosome-wide files.
 *
 * The `ligater` class is responsible for combining multiple input VCF/BCF files
 * (produced by chunked phasing or imputation pipelines) into a single, chromosome-wide
 * output file while maintaining consistent phasing and genotype information.
 *
 * ## Workflow Overview:
 * 1. Parse command-line arguments to get input/output files, thread number, and seed.
 * 2. Validate input/output options and parameters.
 * 3. Read the list of input filenames and initialize internal structures.
 * 4. Ligate input files:
 *    - Open a small number of input files at a time to manage memory usage.
 *    - Iterate through variants in each file in genomic order.
 *    - Handle overlapping regions (chunk overlaps) and update genotype phasing.
 *    - Write processed variants to the output file.
 * 5. Finalize output and report total running time.
 *
 * ## Multi-threading:
 * Supports multi-threaded reading and writing through htslib when more than one thread is specified.
 *
 * ## Input/Output:
 * - Input: List of VCF/BCF files (one per line in a text file).
 * - Output: Single VCF/BCF file covering the entire chromosome.
 * - Optional log file for progress and debug messages.
 *
 * @note
 * The ligation ensures:
 * - Samples are consistent across all input files.
 * - Variants are written in ascending genomic order.
 * - Overlapping regions are correctly phased based on overlap statistics.
 */
class ligater {
public:
	/** 
     * @brief Command-line option descriptions (Boost.Program_options)
     * 
     * Stores all possible options that can be passed via command-line, including input/output files,
     * number of threads, seed for random generator, and help flag.
     */
    bpo::options_description descriptions;

    /**
     * @brief Parsed command-line option values.
     *
     * Populated after parsing with Boost.Program_options. Allows accessing user-specified
     * options like input filenames, output path, threads, seed, etc.
     */
    bpo::variables_map options;

    // FILE DATA

    /**
     * @brief Number of input files to ligate.
     */
    int nfiles;

    /**
     * @brief List of input filenames (VCF/BCF) read from the input file list.
     */
    std::vector<std::string> filenames;

    /**
     * @brief Stores previous reader indices used during phased merging.
     *
     * Helps track which readers were open in overlapping chunks.
     */
    std::vector<int> prev_readers;

    // SAMPLE DATA

    /**
     * @brief Number of samples across all input files.
     */
    int nsamples;

    /**
     * @brief Number of swaps applied during phasing for each half of an overlap.
     *
     * nswap[0]: swaps in the previous half, nswap[1]: swaps in the current half
     */
    std::array<int,2> nswap;

    /**
     * @brief Per-sample phasing swap flags for each half of the overlap.
     *
     * swap_phase[0][i]: whether sample i was swapped in previous half.
     * swap_phase[1][i]: whether sample i was swapped in current half.
     */
    std::array<std::vector<bool>,2> swap_phase;

    /**
     * @brief Count of genotype matches per sample in the overlap region.
     *
     * Used to determine whether phasing swap is needed for a sample.
     */
    std::vector<int> nmatch;

    /**
     * @brief Count of genotype mismatches per sample in the overlap region.
     *
     * Used together with nmatch to compute swap decisions.
     */
    std::vector<int> nmism;

    /**
     * @brief Half-lengths of overlap buffers.
     *
     * Each entry represents nsites_buff/2 for a given overlap chunk.
     */
    std::vector<int> nsites_buff_d2;

    /**
     * @brief Temporary genotype buffers for two overlapping chunks.
     *
     * GTa: genotypes from first chunk, GTb: genotypes from second chunk.
     * mGTa/mGTb store the allocated sizes for each buffer.
     */
    int32_t *GTa, *GTb;
    int32_t mGTa, mGTb;

	//CONSTRUCTOR
	/**
	 * @brief Constructs a new ligater instance.
	 *
	 * Initializes the ligater object with default values.
	 * This constructor does not perform any setup or parameter parsing;
	 * it simply creates an empty instance ready for configuration
	 * before execution.
	 */
	ligater();

	/**
	 * @brief Destroys the ligater instance.
	 *
	 * Cleans up resources associated with the ligater object.
	 * Since no dynamic memory or special cleanup is required
	 * in the current implementation, this destructor is empty.
	 */
	~ligater();

	//PARAMETERS
	/**
	 * @brief Declares command-line options for the ligater program.
	 *
	 * This method defines and groups the available command-line options
	 * into logical categories using Boost Program Options.
	 *
	 * The declared option groups are:
	 * - **Basic options**:
	 *   - `--help` : Displays a help message.
	 *   - `--seed` : Seed for the random number generator (default: `15052011`).
	 *   - `--threads` : Number of threads to use (default: `1`).
	 *
	 * - **Input files**:
	 *   - `--input` : Path to a text file listing all VCF/BCF files to be ligated, one file per line.
	 *
	 * - **Output files**:
	 *   - `--output, -O` : Path for the output ligated (phased) VCF/BCF file.
	 *   - `--log` : Path for the log file.
	 *
	 * The resulting option descriptions are stored in the `descriptions` member
	 * for later parsing and validation.
	 */
	void declare_options();

	/**
	 * @brief Parses and processes command-line arguments for the ligater program.
	 *
	 * This method:
	 * 1. Parses the provided `args` vector using the previously declared
	 *    Boost Program Options `descriptions`.
	 * 2. Stores the parsed values into the `options` member and triggers
	 *    validation with `boost::program_options::notify`.
	 * 3. Displays program metadata (title, authors, contact, version,
	 *    commit information, citation references, and current date).
	 * 4. If `--help` is specified, prints the available command-line options
	 *    and exits.
	 * 5. If the `--log` option is provided, attempts to open the log file.
	 *    If the file cannot be created, logs an error and terminates.
	 *
	 * @param args A vector of strings representing command-line arguments
	 *             (excluding the executable name).
	 *
	 * @throws Exits the program with `exit(0)` if:
	 *         - There is an error parsing arguments.
	 *         - The user requests `--help`.
	 *         - The specified log file cannot be created.
	 */
	void parse_command_line(std::vector < std::string > &);

	/**
	 * @brief Validates the presence and correctness of required command-line options.
	 *
	 * This method performs several checks on the parsed command-line arguments
	 * stored in the `options` member to ensure the ligater program has all
	 * necessary parameters before execution:
	 *
	 * 1. Ensures the `--input` option is specified, indicating the file
	 *    containing the list of VCF/BCF files to ligate.
	 * 2. Ensures the `--output` option is specified, defining the destination
	 *    VCF/BCF file for the ligated output.
	 * 3. If the `--seed` option is provided, verifies that its value is a
	 *    non-negative integer (≥ 0).
	 * 4. Ensures that the `--threads` option specifies a strictly positive integer
	 *    (≥ 1).
	 *
	 * @note If any of these checks fail, the method logs an error message via
	 *       `vrb.error()` and terminates the program.
	 */
	void check_options();

	/**
	 * @brief Prints the key runtime parameters to the verbose output stream.
	 *
	 * This method outputs a summary of important configuration parameters
	 * obtained from the parsed command-line options. The parameters are
	 * displayed in a structured, human-readable format for logging or
	 * user reference.
	 *
	 * Specifically, the method prints:
	 * - The random number generator seed (`--seed` option).
	 * - The number of threads to be used (`--threads` option).
	 *
	 * @note The output is formatted using the `vrb` verbose logging utility
	 *       and numerical values are converted to strings via the `stb.str()` helper.
	 */
	void verbose_options();

	/**
	 * @brief Displays input and output file paths.
	 *
	 * This method outputs a summary of the file-related parameters
	 * provided via command-line options. It is used to confirm to
	 * the user which files will be processed and generated.
	 *
	 * Specifically, the method prints:
	 * - The path to the input list file (`--input`), containing VCF/BCF files to be ligated.
	 * - The output VCF/BCF file path (`--output`).
	 * - The optional log file path (`--log`), if specified.
	 *
	 * @note The output is formatted using the `vrb` verbose logging utility.
	 */
	void verbose_files();

	//FILE I/O
	/**
	 * @brief Reads the list of VCF/BCF files and initializes resources.
	 *
	 * This method performs the initial setup required before file ligation:
	 * 1. Initializes the random number generator seed from the `--seed` option.
	 * 2. Reads the file paths from the list file specified by `--input`.
	 * 3. Stores all file paths into the `filenames` vector.
	 * 4. Counts the total number of files and stores it in `nfiles`.
	 *
	 * @throws std::runtime_error If the input file list is empty.
	 *
	 * @note
	 * - The input list file must contain one VCF/BCF file path per line.
	 * - The verbose logger (`vrb`) outputs the number of files read.
	 */
	void read_files_and_initialise();

	/**
	 * @brief Main entry point for the ligation process with command-line arguments.
	 *
	 * This high-level method orchestrates the entire ligation workflow:
	 * 1. Declares available command-line options (`declare_options()`).
	 * 2. Parses the provided arguments (`parse_command_line()`).
	 * 3. Validates all required parameters (`check_options()`).
	 * 4. Prints details about input/output files (`verbose_files()`).
	 * 5. Prints selected runtime parameters (`verbose_options()`).
	 * 6. Reads the input file list and prepares resources (`read_files_and_initialise()`).
	 * 7. Executes the ligation process (`ligate()` - internal implementation).
	 * 8. Writes the resulting files and finalizes (`write_files_and_finalise()`).
	 *
	 * @param args Vector of command-line arguments (excluding program name).
	 *
	 * @note
	 * This function is typically called once from `main()` after collecting CLI arguments.
	 * The inner call to `ligate()` refers to the core ligation logic, not this wrapper method.
	 */
	void ligate(std::vector < std::string > &);

	/**
	 * @brief Executes the full ligation pipeline for merging chunked genomic data into chromosome-wide VCF/BCF files.
	 *
	 * This method is the main entry point for the ligater module. It manages the complete workflow for reading,
	 * validating, and combining multiple phased VCF/BCF chunks into a single continuous output file per chromosome.
	 *
	 * ## Workflow
	 * 1. **Option Declaration**  
	 *    Calls declare_options() to set up Boost.Program_options parameters for:
	 *      - Random number generator seed
	 *      - Number of threads
	 *      - Input file list path
	 *      - Output file paths (VCF/BCF, log)
	 *
	 * 2. **Command-line Parsing & Validation**  
	 *    - Parses arguments via parse_command_line().
	 *    - Ensures required parameters (`--input`, `--output`) are present using check_options().
	 *    - Validates numerical parameters (positive seed, positive thread count).
	 *
	 * 3. **Verbose Reporting**  
	 *    - Displays file paths, seed, and thread count using verbose_files() and verbose_options().
	 *
	 * 4. **Initialization**  
	 *    - Sets the RNG seed.
	 *    - Reads the list of chunked VCF/BCF file paths from the specified input list file.
	 *    - Stores all file paths in `filenames` and counts them.
	 *    - Aborts if no files are found.
	 *
	 * 5. **Core Ligation**  
	 *    - Invokes ligate() (internal, no-argument version) to:
	 *      - Iterate through each chunk file in `filenames`
	 *      - Load and decode genomic records
	 *      - Merge records in coordinate order
	 *      - Handle overlapping variants and ensure continuity between chunks
	 *      - Maintain per-sample phasing consistency
	 *
	 * 6. **Finalization**  
	 *    - Calls write_files_and_finalise() to:
	 *        - Close output files
	 *        - Release resources
	 *        - Report total running time
	 *
	 * ## Algorithmic Notes
	 * - **Sorting & Merging**: Files are processed in genomic order; variant records are merged using position-based ordering.
	 * - **Overlap Handling**: Overlapping records between chunks are resolved to prevent duplicate output variants.
	 * - **Phasing Integrity**: The haplotype phase is preserved across chunk boundaries.
	 * - **Threading**: Multi-threading is supported for file parsing and record processing; mutexes protect shared structures.
	 *
	 * @param args Vector of command-line arguments, typically from main().
	 * @throws std::runtime_error if required parameters are missing or invalid, or if input files cannot be read.
	 * @note The internal ligate() method must be implemented to perform the actual merging logic.
	 */
	void ligate();

	/**
	 * @brief Finalizes the ligation process and reports runtime statistics.
	 *
	 * This method is called at the end of the ligation workflow to perform final
	 * reporting and cleanup steps. In the current implementation, it primarily
	 * handles verbose output of execution time.
	 *
	 * ## Workflow
	 * 1. **Title Output**  
	 *    Prints a "Finalization:" title block to indicate the end of processing.
	 *
	 * 2. **Runtime Measurement**  
	 *    - Uses the `tac` (timer/chronometer) instance to retrieve the total
	 *      elapsed wall-clock time since the start of the program.
	 *    - Converts the elapsed time to a string using `stb.str()`.
	 *    - Outputs the runtime in seconds via the `vrb.bullet()` logging method.
	 *
	 * ## Notes
	 * - This function could be extended in the future to close output file handles,
	 *   flush buffers, and release allocated resources.
	 * - The measurement includes all stages of processing:
	 *   command-line parsing, option validation, file reading, ligation, and writing.
	 *
	 * @pre The `tac` timer must have been started before the ligation process began.
	 * @post Outputs a human-readable runtime message to the configured verbose/log system.
	 */
	void write_files_and_finalise();
	//void scan_chunks();

	/**
	 * @brief Scans and compares overlapping variant sites between two VCF/BCF files.
	 *
	 * This function opens two indexed variant files (specified by their indices in `filenames`)
	 * and iterates through their overlapping sites starting from a given chromosome and position.
	 * For each overlapping bi-allelic site, genotype concordance and switch phase statistics
	 * are computed to evaluate phasing consistency between the two datasets.
	 *
	 * @param ifname
	 *   The index of the second file in the `filenames` list (1-based).  
	 *   This function will load files at positions `ifname - 2` and `ifname - 1` in `filenames`.
	 *
	 * @param seek_chr
	 *   Chromosome name (as C-string) at which to begin scanning.
	 *
	 * @param seek_pos
	 *   0-based genomic position at which to start scanning.
	 *
	 * ## Algorithm Overview
	 * 1. **Reader Setup**  
	 *    - Initializes a `bcf_srs_t` synchronized reader (`sr`) with indexing required.
	 *    - Configures multithreading based on the `threads` option.
	 *    - Adds two readers: the previous file and the current file.
	 *    - Errors if index creation/loading fails.
	 *
	 * 2. **Seek and Iterate**  
	 *    - Seeks to the specified chromosome and position (`bcf_sr_seek()`).
	 *    - Iterates with `bcf_sr_next_line()` until the overlapping region ends.
	 *
	 * 3. **Filtering & Genotype Extraction**  
	 *    - Ignores non-bi-allelic sites (`n_allele != 2`).
	 *    - Retrieves genotypes (`bcf_get_genotypes`) from both files.
	 *    - If GT is missing, raises an error.
	 *
	 * 4. **Distance & Match/Mismatch Counting**  
	 *    - Calls `update_distances()` to update per-sample match/mismatch statistics.
	 *
	 * 5. **Switch Phase Computation**  
	 *    - For each sample, determines if phasing needs to be swapped
	 *      (`swap_phase[1][i] = nmatch[i] < nmism[i]`).
	 *    - Tracks the number of swaps and collects per-sample het counts.
	 *
	 * 6. **Phase Quality Estimation**  
	 *    - Computes an entropy-based phase quality metric (`phaseQ`) for each sample:
	 *      \f[
	 *        Q = 99 \times \frac{0.7 + p \log(p) + (1-p) \log(1-p)}{0.7}
	 *      \f]
	 *      where \( p = \frac{\text{nmatch}}{\text{nmatch} + \text{nmism}} \).
	 *
	 * 7. **Buffer & Logging**  
	 *    - Stores half of the number of overlapping sites in `nsites_buff_d2`.
	 *    - Logs buffer index, genomic range, overlap length, total sites scanned,
	 *      average heterozygous sites per sample, switch rate, and average phaseQ.
	 *
	 * @throws std::runtime_error
	 *   If input files cannot be opened/indexed or if GT fields are missing at overlapping sites.
	 *
	 * @pre Both input VCF/BCF files must be indexed (`.csi` or `.tbi` present).
	 * @post Updates `nsites_buff_d2`, `swap_phase[1]`, and logs overlap statistics.
	 */
	void scan_overlap(const int ifname,const char* seek_chr, int seek_pos);


	//FUNCTIONS
	void updateHS(int *);
	int update_switching();

	/**
	 * @brief Update per-sample genotype phase match/mismatch counts for overlapping variants.
	 *
	 * This method compares phased genotype calls from two overlapping variant datasets
	 * (previously loaded into `GTa` and `GTb` arrays) to determine whether the phase
	 * between the two datasets agrees or disagrees for each sample.
	 *
	 * ## Algorithm
	 * 1. Iterate over all samples (`nsamples`).
	 * 2. For each sample:
	 *    - Extract the two alleles from the first dataset (`gta`) and second dataset (`gtb`).
	 *    - Skip this sample if:
	 *      - One or both alleles are missing (`bcf_gt_is_missing()`).
	 *      - Genotypes are not phased in either dataset (`bcf_gt_is_phased()`).
	 *      - The genotype is homozygous in either dataset (only heterozygotes are informative for phase comparison).
	 *    - Compare alleles between datasets:
	 *      - **Match case:** `(a0 == b0 && a1 == b1)` — The haplotype phase is identical.
	 *      - **Swap case:** `(a0 == b1 && a1 == b0)` — The haplotype phase is reversed.
	 *    - Depending on the current `swap_phase[0][i]` flag:
	 *      - If `swap_phase[0][i]` is `false`:
	 *         - Match case increments `nmatch[i]`.
	 *         - Swap case increments `nmism[i]`.
	 *      - If `swap_phase[0][i]` is `true`:
	 *         - Match case increments `nmism[i]`.
	 *         - Swap case increments `nmatch[i]`.
	 *
	 * ## Purpose
	 * This function is used during overlap scanning (`scan_overlap()`) to accumulate
	 * per-sample statistics on how well the phasing agrees between two overlapping blocks.
	 * These counts are later used to determine whether to flip the phase in one block
	 * when joining them together during ligation.
	 *
	 * @note
	 * - `GTa` and `GTb` are assumed to be populated with genotype data for the current variant.
	 * - Only heterozygous sites with phased data contribute to phase matching.
	 * - `swap_phase` array is used to track if a previous decision has already flipped the phase for a sample.
	 *
	 * @see scan_overlap()
	 */
	void update_distances();

	/**
	 * @brief Apply phase flips to a VCF/BCF record for selected samples.
	 *
	 * This function updates the genotype (GT) field of a single variant record based
	 * on the previously determined `swap_phase` decisions. It is called during ligation
	 * to ensure that overlapping blocks have consistent phasing across samples.
	 *
	 * ## Algorithm
	 * 1. Extract genotypes for all samples from the variant record using `bcf_get_genotypes()`.
	 *    - If the GT field is absent (`nGTs <= 0`), return immediately.
	 * 2. Iterate over all samples:
	 *    - Skip samples where `swap_phase[uphalf][i]` is false (no phase flip needed).
	 *    - Skip if the genotype is missing or unphased.
	 *    - Otherwise, swap the alleles while preserving the phased bit:
	 *       - `gt[0]` becomes the phased encoding of the second allele.
	 *       - `gt[1]` becomes the phased encoding of the first allele.
	 * 3. Update the VCF/BCF record with the modified genotype array using `bcf_update_genotypes()`.
	 *
	 * ## Purpose
	 * Ensures phase consistency across overlapping variant blocks when merging multiple
	 * VCF/BCF files into chromosome-wide outputs. This step is essential to correct
	 * switch errors introduced in previous blocks.
	 *
	 * @param hdr Pointer to the BCF/VCF header corresponding to `line`.
	 * @param line Pointer to the BCF/VCF record to update.
	 * @param uphalf Boolean flag indicating whether the current record belongs to the
	 *        “upper half” block (overlap buffer) for which phase flipping may apply.
	 *
	 * @note
	 * - `GTa` is reused as a temporary buffer for the genotypes of the current record.
	 * - Only phased, heterozygous genotypes are modified.
	 * - Relies on `swap_phase[uphalf][i]` being correctly set during `scan_overlap()`.
	 *
	 * @see scan_overlap(), update_distances()
	 */
	void phase_update(bcf_hdr_t *hdr, bcf1_t *line, const bool uphalf);

	/**
	 * @brief Remove most INFO fields from a VCF/BCF record, preserving only AC and AN.
	 *
	 * This function iterates over all INFO fields in a variant record and frees their
	 * memory, effectively clearing them from the record. Only the essential fields
	 * `AC` (Allele Count) and `AN` (Allele Number) are retained.
	 *
	 * ## Algorithm
	 * 1. Ensure the INFO fields are unpacked (`BCF_UN_INFO`) for access.
	 * 2. Loop through all INFO fields:
	 *    - Retrieve the field key using `bcf_hdr_int2id()`.
	 *    - Skip fields with key `"AC"` or `"AN"`.
	 *    - For other fields:
	 *       - Free the memory buffer if allocated (`vptr_free`).
	 *       - Mark the field as cleared by setting `vptr = NULL` and resetting offsets.
	 *       - Flag the record as having modified INFO fields (`shared_dirty |= BCF1_DIRTY_INF`).
	 *
	 * ## Purpose
	 * Reduces memory usage and file size when merging or ligating multiple VCF/BCF files
	 * by removing non-essential INFO annotations. Essential allele count information
	 * is preserved for downstream analyses.
	 *
	 * @param hdr Pointer to the BCF/VCF header corresponding to `line`.
	 * @param line Pointer to the BCF/VCF record to update.
	 *
	 * @note
	 * - Only the INFO fields `"AC"` and `"AN"` are preserved; all others are removed.
	 * - Must be called before writing or processing records where INFO is not needed.
	 *
	 * @see ligater::ligate(), ligater::phase_update()
	 */
	void remove_info(bcf_hdr_t *hdr, bcf1_t *line);

	/**
	 * @brief Remove all FORMAT fields from a VCF/BCF record except the genotype (GT) field.
	 *
	 * This function iterates over all FORMAT fields in a variant record and frees their
	 * memory, effectively clearing them from the record. The genotype field (`GT`) is
	 * preserved as it is essential for phasing and downstream analyses.
	 *
	 * ## Algorithm
	 * 1. Ensure the FORMAT fields are unpacked (`BCF_UN_FMT`) for access.
	 * 2. Loop through all FORMAT fields:
	 *    - Retrieve the field key using `bcf_hdr_int2id()`.
	 *    - Skip the `GT` field (key `"GT"`).
	 *    - For all other fields:
	 *       - Free the memory buffer if allocated (`p_free`).
	 *       - Mark the record as having modified individual-level data (`indiv_dirty = 1`).
	 *       - Set the pointer to `NULL` to indicate removal.
	 *
	 * ## Purpose
	 * Reduces memory usage and file size when merging or ligating multiple VCF/BCF files
	 * by removing non-essential FORMAT annotations while retaining the genotype information.
	 *
	 * @param hdr Pointer to the BCF/VCF header corresponding to `line`.
	 * @param line Pointer to the BCF/VCF record to update.
	 *
	 * @note
	 * - Only the FORMAT field `"GT"` is preserved; all other FORMAT fields are removed.
	 * - Must be called before writing or processing records where FORMAT fields are not needed.
	 *
	 * @see ligater::ligate(), ligater::remove_info(), ligater::phase_update()
	 */
	void remove_format(bcf_hdr_t *hdr, bcf1_t *line);

	/**
	 * @brief Write a single VCF/BCF record to the output file, applying phasing updates if needed.
	 *
	 * This function handles the final processing of a variant record before writing it
	 * to the output VCF/BCF file. It optionally updates the genotype phasing, translates
	 * the record to match the output header, and writes it to disk.
	 *
	 * ## Algorithm
	 * 1. Check if phasing swaps are needed for the current half (`nswap[uphalf]`):
	 *    - If yes, call `phase_update()` to adjust the genotype (`GT`) field accordingly.
	 * 2. Optionally remove INFO and non-GT FORMAT fields (currently commented out):
	 *    - `remove_info()` can clear unnecessary INFO annotations.
	 *    - `remove_format()` can clear all FORMAT fields except `GT`.
	 * 3. Translate the record to the output header using `bcf_translate()`:
	 *    - Ensures the variant record matches the structure and sample order of `out_hdr`.
	 * 4. Write the record to the output file using `bcf_write()`:
	 *    - Throws an error via `vrb.error()` if writing fails.
	 *
	 * ## Purpose
	 * Centralizes per-record output handling during the ligation process, including:
	 * - Maintaining consistent headers across multiple input files.
	 * - Applying phasing corrections for overlapping chunks.
	 * - Ensuring only the relevant data is written to the output file.
	 *
	 * @param fd Pointer to the output VCF/BCF file (`htsFile`) where the record will be written.
	 * @param out_hdr Pointer to the output header (`bcf_hdr_t*`) that defines the file format and samples.
	 * @param hdr_in Pointer to the input header corresponding to the current record.
	 * @param line Pointer to the variant record to be processed and written.
	 * @param uphalf Boolean flag indicating which “half” of overlapping chunks is being processed.
	 *
	 * @note
	 * - Phasing updates are applied only if swaps are required for the current half.
	 * - INFO and FORMAT field removal is currently commented out; can be enabled to reduce file size.
	 * - The function ensures that the output record conforms to the output header.
	 *
	 * @see ligater::phase_update(), ligater::remove_info(), ligater::remove_format()
	 */
	void write_record(htsFile *, bcf_hdr_t * ,  bcf_hdr_t * ,bcf1_t *, const bool uphalf);

};

#endif


