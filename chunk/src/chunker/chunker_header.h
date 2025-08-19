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

#ifndef _CHUNKER_H
#define _CHUNKER_H

#include "otools.h"
#include "gmap_reader.h"

/**
 * @struct chunk_info
 * @brief Holds information about genomic chunks and their buffered regions.
 * 
 * This struct stores start and stop positions of buffered and chunk regions,
 * along with properties such as size in centiMorgans (cM), size in base pairs (Mb),
 * and counts of variants within each chunk. It provides functionality to reset data,
 * add new chunks, and output chunk information to a file.
 */
struct chunk_info
{
	/// Chromosome identifier (e.g., "chr1")
	std::string chr;

	/// Start positions of buffered regions in base pairs
	std::vector<long int> buf_start;

	/// Stop positions of buffered regions in base pairs
	std::vector<long int> buf_stop;

	/// Start positions of chunk regions (phased output region) in base pairs
	std::vector<long int> cnk_start;

	/// Stop positions of chunk regions in base pairs
	std::vector<long int> cnk_stop;

	/// Size of each chunk window in centiMorgans (cM)
	std::vector<float> curr_window_cm_size;

	/// Size of each chunk window in base pairs (Mb)
	std::vector<long int> curr_window_mb_size;

	/// Number of total variants in each chunk window
	std::vector<long int> curr_window_count;

	/// Number of common variants in each chunk window
	std::vector<long int> curr_window_common_count;

	/**
     * @brief Default constructor.
     * 
     * Initializes an empty chunk_info object.
     */
	chunk_info(){};

	/**
     * @brief Constructor that sets the chromosome identifier.
     * 
     * @param _chr Chromosome ID (e.g., "chr1")
     */
	chunk_info(std::string _chr) : chr(_chr){};

	/**
     * @brief Clears all stored chunk information.
     */
	void reset()
	{
		buf_start.clear();
		buf_stop.clear();
		cnk_start.clear();
		cnk_stop.clear();
		curr_window_cm_size.clear();
		curr_window_mb_size.clear();
		curr_window_count.clear();
		curr_window_common_count.clear();
	}

	/**
     * @brief Adds a new chunk with buffered and chunk boundaries and statistics.
     * 
     * @param _buf_start Start position of the buffered region (bp).
     * @param _buf_stop Stop position of the buffered region (bp).
     * @param _cnk_start Start position of the chunk region (bp).
     * @param _cnk_stop Stop position of the chunk region (bp).
     * @param _curr_window_cm_size Size of the chunk window in centiMorgans (cM).
     * @param _curr_window_mb_size Size of the chunk window in base pairs.
     * @param _curr_window_count Total number of variants in the chunk.
     * @param _curr_window_common_count Number of common variants in the chunk.
     * 
     * @note This function asserts that buffered regions enclose chunk regions properly.
     */
	void add_chunk(long int _buf_start, long int _buf_stop, long int _cnk_start, long int _cnk_stop, float _curr_window_cm_size, long int _curr_window_mb_size, long int _curr_window_count, long int _curr_window_common_count)
	{
		assert(_buf_start < _buf_stop);
		assert(_cnk_start < _cnk_stop);
		assert(_buf_start <= _cnk_start);
		assert(_buf_stop >= _buf_stop);

		buf_start.push_back(_buf_start);
		buf_stop.push_back(_buf_stop);
		cnk_start.push_back(_cnk_start);
		cnk_stop.push_back(_cnk_stop);
		curr_window_cm_size.push_back(_curr_window_cm_size);
		curr_window_mb_size.push_back(_curr_window_mb_size);
		curr_window_count.push_back(_curr_window_count);
		curr_window_common_count.push_back(_curr_window_common_count);
	}

	/**
     * @brief Outputs the stored chunk information to a file stream.
     * 
     * This function also adjusts chunk boundaries by averaging overlapping edges
     * between consecutive chunks before writing.
     * 
     * @param fd Output file stream to write the chunk information.
     */
	void output_to_file(output_file& fd)
	{
		for (long int i=0; i<buf_start.size(); ++i)
		{
			if (i<buf_start.size()-1)
			{
				long int mean_curr_stop_next_start = (cnk_stop[i] + cnk_start[i+1]) /2;
				cnk_stop[i] = mean_curr_stop_next_start;
				cnk_start[i + 1] = mean_curr_stop_next_start+1;
			}
			fd << i << "\t" << chr << "\t" << chr << ":"<< buf_start[i] << "-" << buf_stop[i] << "\t" << chr << ":" << cnk_start[i] << "-" << cnk_stop[i] << "\t" << curr_window_cm_size[i] << "\t" << curr_window_mb_size[i] << "\t" << curr_window_count[i] <<"\t" << curr_window_common_count[i] << std::endl;
		}
	}
};


/**
 * @class chunker
 * @brief A class to split genomic regions into chunks based on genetic and physical coordinates.
 * 
 * This class is designed to process variant data within a specified genomic region
 * and split it into smaller chunks suitable for downstream analysis, such as imputation.
 * It supports command-line options, reads variant positions, applies genetic maps,
 * and uses recursive or sequential algorithms to define chunk boundaries with buffers.
 * 
 * The chunking criteria include window size in centiMorgans (cM), physical distance (Mb),
 * and variant counts, with optional buffering around chunks.
 * 
 * Major functionalities:
 * - Command-line parsing and option validation.
 * - Region parsing and coordinate management.
 * - Reading variant data and genetic maps.
 * - Recursive and sequential chunk splitting algorithms.
 * - Chunk metadata storage and output.
 * 
 * @note Uses Boost.Program_options for argument parsing.
 * 
 * @author Simone Rubinacci
 * @author Olivier Delaneau
 * @version GLIMPSE2_chunk (commit and release IDs)
 */
class chunker {
public:
	//COMMAND LINE OPTIONS
	/**
     * @brief Description of command line options, used for argument parsing and help messages.
     */
	bpo::options_description descriptions;

	/**
     * @brief Map storing the parsed command line option values.
     */
	bpo::variables_map options;

	/**
     * @brief Minor Allele Frequency threshold for sparsifying variants.
     * 
     * Used to filter out variants below this MAF in chunking.
     */
	float sparse_maf;

	/**
     * @brief Chromosome identifier string (e.g., "chr1").
     */
	std::string chrID;

	/**
     * @brief Positions of all variants in centiMorgans (cM).
     * 
     * This vector holds genetic map positions for all variants in the region.
     */
	std::vector < float > positions_all_cm;

	/**
     * @brief Positions of all variants in base pairs (Mb).
     * 
     * Positions are stored as integers representing base pair locations.
     */
	std::vector < long int > positions_all_mb;

	/**
     * @brief Positions of common variants in centiMorgans (cM).
     * 
     * Subset of variants passing common variant filtering.
     */
	std::vector < float > positions_common_cm;

	/**
     * @brief Positions of common variants in base pairs (Mb).
     */
	std::vector < long int > positions_common_mb;

	/**
     * @brief Mapping from all variants to common variants indices.
     * 
     * For each variant index in `positions_all_mb`, provides the corresponding
     * index in the common variant set.
     */
	std::vector < long int > all2common;

	/**
     * @brief Mapping from common variants to all variants indices.
     * 
     * For each common variant index, provides the corresponding index in the
     * full variant set.
     */
	std::vector < long int > common2all;

	/**
     * @brief Multimap of variant base pair positions to their indices.
     * 
     * This associative container maps base pair positions (keys) to
     * variant indices (values) in the full variant set.
     */
	std::multimap < long int, long int > map_positions_all;	//associative container of variant with position in bp

	/**
     * @brief Vector storing lengths of chunks in centiMorgans (cM).
     * 
     * Populated during chunk splitting.
     */
	std::vector<float> chunk_cm_length;

	/**
     * @brief Vector storing lengths of chunks in base pairs (Mb).
     * 
     * Populated during chunk splitting.
     */
	std::vector<long int> chunk_mb_length;

	/**
     * @brief Vector storing counts of common variants per chunk.
     * 
     * Populated during chunk splitting.
     */
	std::vector<long int> chunk_common_count;

	 /**
     * @brief Information object holding details about chunk boundaries and metadata.
     */
	chunk_info cnk_info;


	//PARAMETERS
	
	 /**
     * @brief Genomic region string provided by the user (e.g., "chr1:1000-2000").
     */
	std::string gregion;

	/**
     * @brief Start position of the region (in base pairs).
     */
	long int start;

	/**
     * @brief Stop position of the region (in base pairs).
     */
	long int stop;

	/**
     * @brief Flag indicating whether the entire chromosome is being processed.
     */
	bool whole_chr;

	/**
     * @brief Length of the contig or chromosome (in base pairs).
     */
	long int contig_len;

	/**
     * @brief Minimum window size in centiMorgans for chunk splitting.
     */
	float window_cm;

	/**
     * @brief Minimum window size in base pairs (Mb) for chunk splitting.
     */
	float window_mb;

	/**
     * @brief Minimum number of variants per window for chunk splitting.
     */
	long int window_count;

	/**
     * @brief Buffer size in centiMorgans added around chunks.
     */
	float buffer_cm;

	/**
     * @brief Buffer size in base pairs (Mb) added around chunks.
     */
	float buffer_mb;

	/**
     * @brief Buffer size in number of variants added around chunks.
     */
	long int buffer_count;

	//CONSTRUCTOR
	chunker();
	~chunker();

		/**
	 * @brief Reads variant data from a VCF/BCF file into chunker data structures.
	 *
	 * This function opens the specified VCF/BCF file, restricts reading to the given genomic
	 * region, and iterates through all biallelic variants. For each variant, it extracts allele
	 * counts (AC) and allele numbers (AN) from the INFO fields, computes the minor allele frequency (MAF),
	 * and classifies the variant as common or rare based on the @ref sparse_maf threshold.
	 *
	 * Variant positions, along with mapping tables between all variants and common variants,
	 * are stored in the chunker's internal member variables for later use in imputation/phasing.
	 *
	 * If @ref whole_chr is enabled, the function will also retrieve and store the contig length
	 * for the given chromosome from the VCF/BCF header or index. Requires the input file to have
	 * an accompanying index (.tbi/.csi) file.
	 *
	 * @param fmain    Path to the input VCF or BCF file (must be indexed).
	 * @param region   Genomic region string in format "CHR" or "CHR:START-END".
	 * @param nthreads Number of threads to use for reading (if > 1, enables multi-threaded reading).
	 *
	 * @throws std::runtime_error If:
	 *         - The file cannot be opened.
	 *         - The index cannot be loaded.
	 *         - The region cannot be set.
	 *         - The AC or AN INFO fields are missing or invalid.
	 *         - The file type is not recognized as VCF or BCF.
	 *
	 * @note Only biallelic variants (n_allele == 2) are processed. Multiallelic variants are skipped.
	 * @note AC and AN must be present in the INFO field and be single integer values.
	 * @note Positions are stored in 1-based coordinates (VCF format is 0-based).
	 * @note Updates the following member variables:
	 *       - positions_all_mb
	 *       - positions_common_mb
	 *       - map_positions_all
	 *       - all2common
	 *       - common2all
	 *       - contig_len (if whole_chr is true)
	 *       - chrID (chromosome ID of the first common variant)
	 *
	 * @details
	 * **Algorithm workflow:**
	 * 1. **Initialization**
	 *    - Start a timer for performance logging.
	 *    - Initialize a synchronous VCF/BCF reader (`bcf_srs_t`) from htslib.
	 *    - Configure reader to not collapse multi-allelic sites and to require an index.
	 *    - Enable multi-threading if `nthreads > 1`.
	 *
	 * 2. **Set region filter**
	 *    - Restrict reading to the given `region` (chromosome or range).
	 *    - Abort with an error if the region cannot be set.
	 *
	 * 3. **Open input file**
	 *    - Add the file as a reader.
	 *    - If opening or index loading fails, terminate with an error message.
	 *
	 * 4. **Optional: Get contig length** (only if `whole_chr == true`)
	 *    - Detect if the file is VCF or BCF and load the appropriate index (`.tbi` or `.csi`).
	 *    - Retrieve sequence names from the index.
	 *    - Match the `region` chromosome and attempt to read its length from the VCF header.
	 *    - If unavailable, assign a fallback contig length.
	 *
	 * 5. **Prepare buffers**
	 *    - Initialize counters for total and common variants.
	 *    - Allocate arrays for reading AC and AN INFO field values.
	 *
	 * 6. **Read and process variants**
	 *    - Loop over all lines in the file/region:
	 *      - Keep only **biallelic** variants (`n_allele == 2`).
	 *      - For the first common variant, store the chromosome ID (`chrID`).
	 *      - Retrieve AC and AN values and validate them.
	 *      - Compute **minor allele frequency (MAF)**:
	 *        \f[
	 *        MAF = \min\left(\frac{AN - AC}{AN}, \frac{AC}{AN}\right)
	 *        \f]
	 *      - Determine if the variant is **common** (`MAF >= sparse_maf`).
	 *      - Store the variant position (converted from 0-based to 1-based).
	 *      - Update:
	 *        - `positions_all_mb` with all variant positions.
	 *        - `map_positions_all` mapping position → index.
	 *        - `all2common` mapping all variants to common variant indices.
	 *      - If common, also update:
	 *        - `positions_common_mb` with positions of common variants.
	 *        - `common2all` mapping common variant indices back to all variant indices.
	 *
	 * 7. **Cleanup**
	 *    - Destroy the BCF/VCF reader.
	 *    - Free AC/AN buffers.
	 *
	 * 8. **Final reporting**
	 *    - Log total, rare, and common variant counts.
	 *    - Log elapsed processing time.
	 *
	 * @see sparse_maf
	 * @see whole_chr
	 */
	void readData(std::string fmain, std::string region, long int nthreads);

	//PARAMETERS

	/**
	 * @brief Declares the command-line options for the chunking process.
	 *
	 * This function uses Boost.Program_options to define:
	 * - **Basic options**: Help message, random seed, number of threads.
	 * - **Input parameters**: Input VCF/BCF file, target region, genetic map, sparse MAF threshold.
	 * - **Chunking parameters**: Minimum window and buffer sizes (in cM, Mb, or variant count).
	 * - **Algorithm parameters**: Recursive or sequential algorithms, uniform variant distribution.
	 * - **Output files**: Chunk output file, log file.
	 *
	 * All defined option groups are added to the `descriptions` member
	 * for later parsing.
	 */
	void declare_options();

	/**
	 * @brief Parses and processes the command-line arguments for the chunker module.
	 *
	 * This function uses Boost.Program_options to:
	 * 1. Parse the provided command-line arguments.
	 * 2. Store the results in the `options` member.
	 * 3. Trigger notifications for any registered option callbacks.
	 *
	 * Additional behavior:
	 * - If a log file is specified (`--log`), attempts to open it; 
	 *   reports an error and exits on failure.
	 * - Prints program metadata such as version, authors, contact info, 
	 *   citation references, and run date.
	 * - If the `--help` flag is present, displays the list of available options and exits.
	 *
	 * @param args Reference to a vector of strings containing the command-line arguments.
	 *
	 * @throws boost::program_options::error If argument parsing fails.
	 *         The function catches this internally, prints an error message, and exits the program.
	 *
	 * @note This function terminates the program (`exit(0)`) in several cases:
	 *       - Parsing errors
	 *       - Failure to open the specified log file
	 *       - `--help` flag is provided
	 */
	void parse_command_line(std::vector < std::string > &);

	/**
	 * @brief Validates command-line options for the chunking process.
	 *
	 * This method checks that all required options are provided and that 
	 * their values are within acceptable ranges. If any validation fails, 
	 * an error is reported and execution is terminated.
	 *
	 * **Validations performed:**
	 * - Ensures that the following required options are specified:
	 *   - `--input` : Path to the input file.
	 *   - `--region` : Genomic region to be split (ideally a chromosome).
	 *   - `--output` : Path to the phased output file.
	 * - Ensures that `--seed` (if provided) is a positive integer.
	 * - Ensures that `--threads` is a strictly positive integer.
	 * - Ensures that all window size parameters (`--window-cm`, `--window-mb`, `--window-count`)
	 *   are strictly positive.
	 * - Ensures that all buffer size parameters (`--buffer-cm`, `--buffer-mb`, `--buffer-count`)
	 *   are strictly positive.
	 * - Validates `--sparse-maf` to be in the range [0, 0.5) with recommended values
	 *   between 0.001 and 0.01.
	 * - Ensures that exactly one of the two algorithms is selected: `--recursive` or `--sequential`.
	 *
	 * @throws vrb.error if:
	 *   - A required option is missing.
	 *   - An option has an invalid value (e.g., non-positive sizes, out-of-range parameters).
	 *   - Both or neither of `--recursive` and `--sequential` are selected.
	 *
	 * @note This function must be called after parsing command-line arguments.
	 */
	void check_options();

	/**
	 * @brief Prints the current chunking configuration to the verbose output.
	 *
	 * This function outputs the values of all relevant parameters 
	 * in a structured, human-readable format. It is intended to be 
	 * called after command-line arguments are parsed and validated 
	 * (e.g., by `check_options()`), so that the user can verify 
	 * the configuration before processing begins.
	 *
	 * **Displayed parameters include:**
	 * - **Sparse MAF threshold** (`--sparse-maf`)
	 * - **Algorithm**: "Sequential" if `--sequential` is specified, otherwise "Recursive".
	 * - **Recombination rate mode**: Either from a genetic map (`--map`) or constant rate (1 cM/Mb).
	 * - **Minimum window size** in centiMorgans, megabases, and variant count.
	 * - **Minimum buffer size** in centiMorgans, megabases, and variant count.
	 * - **Seed** value for randomization (`--seed`).
	 * - **Number of threads** (`--threads`).
	 *
	 * @note This function produces output only for logging/debugging 
	 *       purposes; it does not modify any state or perform validation.
	 */
	void verbose_options();

	/**
	 * @brief Displays information about input/output files and genetic map settings.
	 *
	 * This function logs to the console (via `vrb`) details of the file-related 
	 * parameters provided through `this->options`. It includes the input VCF path, 
	 * analysis region, genetic map file (if available), output file path, and 
	 * optional log file path.  
	 *
	 * If a genetic map is not specified (`--map` missing), the function issues a 
	 * warning indicating that recombination rates will default to 1 cM/Mb and 
	 * advises setting `--window-cm` and `--buffer-cm` to zero.
	 *
	 * @note The function assumes `this->options` contains all necessary keys 
	 *       (`input`, `region`, `output`, and optionally `map`, `log`).
	 *
	 * Example output:
	 * @code
	 * Files:
	 *   Input VCF            : [example.vcf.gz]
	 *   Region               : [chr1:1-1000000]
	 *   Genetic Map          : [map.txt]
	 *   Output file          : [output.vcf.gz]
	 *   Output LOG           : [log.txt]
	 * @endcode
	 */
	void verbose_files();

	//FILE I/O

	/**
	 * @brief Splits a genomic region recursively with chunk information reset.
	 *
	 * This function initializes the chunk information state (`cnk_info.reset()`)
	 * before delegating the splitting operation to `split_recursive_no_reset()`.
	 * It is typically used when starting a new chunking process to ensure that 
	 * no residual state from previous operations affects the results.
	 *
	 * @param fd        Reference to the output file handler used for writing results.
	 * @param cidx      Reference to the current chunk index counter, updated during processing.
	 * @param chr       Chromosome identifier (e.g., "chr1") for the region to split.
	 * @param start_idx Start position (index) of the region to process.
	 * @param stop_idx  End position (index) of the region to process.
	 *
	 * @see split_recursive_no_reset()
	 */
	void split_recursive(output_file &, long int &, std::string &, long int, long int);

	/**
	 * @brief Recursively splits a genomic region into chunks without resetting chunk state.
	 *
	 * This function divides a genomic interval into smaller windows based on genetic
	 * (cM), physical (Mb), and variant count constraints. It evaluates whether each
	 * half of the region meets the minimum size requirements and, if so, recurses
	 * further. Otherwise, it finalizes the chunk, applies a buffer, and stores
	 * chunk metadata.
	 *
	 * The method differs from `split_recursive()` in that it does not reset the
	 * chunk information (`cnk_info`), allowing continued accumulation of chunk data
	 * across recursive calls.
	 *
	 * @param fd         Reference to the output file handler for writing results.
	 * @param cidx       Reference to the current chunk index counter, incremented 
	 *                   when terminal chunks are created.
	 * @param chr        Chromosome identifier (e.g., "chr1") for the region to split.
	 * @param start_idx  Index of the first variant in the region.
	 * @param stop_idx   Index of the last variant in the region.
	 *
	 * @throws If variant indices are invalid or inconsistent with common variant
	 *         mappings (`all2common`), the function calls `vrb.error()` and stops.
	 *
	 * @note
	 * - Splitting is based on:
	 *   - Minimum genetic length (`window_cm`)
	 *   - Minimum physical length (`window_mb`)
	 *   - Minimum common variant count (`window_count`)
	 * - Buffers are added around terminal windows via `add_buffer()`.
	 * - Information about terminal chunks is stored in `cnk_info`.
	 *
	 * @see split_recursive(), add_buffer(), cnk_info.add_chunk()
	 */
	void split_recursive_no_reset(output_file &, long int &, std::string &, long int, long int);

	/**
	 * @brief Splits genomic positions into sequential chunks based on window size constraints.
	 *
	 * This function iterates through genomic position indices (`positions_all_mb` and `positions_all_cm`)
	 * and creates sequential chunks that satisfy one or more of the following criteria:
	 * - Minimum number of variants (`window_count`)
	 * - Minimum physical distance in megabases (`window_mb`)
	 * - Minimum genetic distance in centimorgans (`window_cm`)
	 *
	 * It accounts for an additional buffer region on each side of the chunk (`add_buffer()`), 
	 * adjusts the boundaries for whole-chromosome processing (`whole_chr`), 
	 * and records chunk metadata such as physical length, genetic length, and variant count.
	 *
	 * The function updates an internal `cnk_info` object with the computed chunk boundaries 
	 * and increments a chunk index counter (`cidx`) for each chunk generated.
	 *
	 * @param fd               Output file handle (type `output_file&`) where chunk information 
	 *                         may be written if `output_to_file` is `true`.
	 * @param cidx             Reference to the current chunk index counter. 
	 *                         Incremented after each chunk is generated.
	 * @param chr              Chromosome name or identifier for the current region being processed.
	 * @param start_idx        Start index (0-based) in `positions_all_mb` for chunk generation.
	 * @param stop_idx         End index (0-based, exclusive) in `positions_all_mb` for chunk generation.
	 *                         Must be greater than `start_idx`.
	 * @param output_to_file   If `true`, chunk information will be written to `fd`. 
	 *                         If `false`, chunks are only recorded internally.
	 *
	 * @note
	 * - **Window stopping logic:**  
	 *   The algorithm determines `curr_window_stop_idx` using `window_count` in the common variant index space, 
	 *   then extends the window until it meets the minimum cm/mb/count thresholds.
	 * - **Buffer addition:**  
	 *   `add_buffer()` is called to extend the window to include flanking variants within the configured buffer size.
	 * - **Whole chromosome special case:**  
	 *   If `whole_chr` is true, the first chunk starts at position 1, and the last chunk extends to `contig_len`.
	 * - **Index mapping:**  
	 *   `all2common` and `common2all` map between the "all variants" space and the "common variants" space.
	 * - **Recording:**  
	 *   Metadata is pushed to `chunk_cm_length`, `chunk_mb_length`, and `chunk_common_count` vectors.
	 *
	 * @warning
	 * - `stop_idx` must be > `start_idx`; otherwise, the assertion fails.
	 * - This function assumes `positions_all_mb` and `positions_all_cm` are non-empty and sorted.
	 * - `window_count`, `window_mb`, `window_cm`, and `buffer_mb` must be set before calling this function.
	 *
	 * @internal
	 * Complexity is approximately O(N) in the number of positions processed between `start_idx` and `stop_idx`.
	 *
	 * **Example:**
	 * @code
	 * long int chunk_index = 0;
	 * std::string chromosome = "chr1";
	 * chunker my_chunker;
	 * my_chunker.split_sequential(output_file_handle, chunk_index, chromosome, 0, positions_all_mb.size(), true);
	 * @endcode
	 */
	void split_sequential(output_file &, long int &, std::string &, long int, long int, const bool);

	/**
	 * @brief Adds buffer regions to the left and right boundaries of a genomic interval.
	 *
	 * This function extends the given genomic interval `[start_idx, stop_idx]` by adding 
	 * buffer regions on both sides based on constraints in centimorgans (cM), megabases (Mb),
	 * and variant count. The extended boundaries are returned via `left_idx` and `right_idx`.
	 *
	 * The buffering is done as follows:
	 * - **Left Buffer**: Starts from `start_idx` and moves left until the buffer size in Mb,
	 *   cM, or variant count meets the specified thresholds.
	 * - **Right Buffer**: Starts from `stop_idx` and moves right until the buffer size in Mb,
	 *   cM, or variant count meets the specified thresholds.
	 * 
	 * If the boundaries reach the start or end of the chromosome, buffering stops.
	 *
	 * @param[in] start_idx    Index of the first variant in the target region (0-based).
	 * @param[in] stop_idx     Index of the last variant in the target region (0-based).
	 * @param[out] left_idx    Updated index marking the new left boundary after buffering.
	 * @param[out] right_idx   Updated index marking the new right boundary after buffering.
	 *
	 * @note The function uses the following class members:
	 * - `positions_all_mb` : Vector of variant positions in Mb.
	 * - `positions_all_cm` : Vector of variant positions in cM.
	 * - `buffer_count`     : Minimum number of variants to include in the buffer.
	 * - `buffer_mb`        : Minimum buffer size in Mb.
	 * - `buffer_cm`        : Minimum buffer size in cM.
	 *
	 * @warning The function assumes `positions_all_mb` and `positions_all_cm` are aligned
	 * and have at least `stop_idx + 1` elements.
	 *
	 * @complexity O(n) in worst case, where `n` is the number of variants in the buffer region.
	 *
	 * @see chunker::positions_all_mb, chunker::positions_all_cm
	 */
	void add_buffer(const long int, const long int, long int&, long int&);

	/**
	 * @brief Main chunking procedure for splitting genomic data into analysis windows.
	 *
	 * This method drives the chunking process by:
	 *  - Initializing random seeds and key chunking parameters from the command-line options.
	 *  - Building genomic coordinates and reading input variant data (optionally filtered by region).
	 *  - Loading or generating genetic map positions.
	 *  - Splitting the data into chunks using either recursive or sequential algorithms.
	 *  - Optionally performing a second pass to create chunks with a uniform number of variants.
	 *  - Writing chunk metadata to the specified output file.
	 *
	 * **Workflow Diagram:**
	 * @dot
	 * digraph ChunkerFlow {
	 *     node [shape=box, fontname="Helvetica", fontsize=10];
	 *
	 *     start [label="Start: Load parameters & init random seed"];
	 *     coords [label="Build genomic coordinates"];
	 *     readData [label="Read input variant data"];
	 *     checkEmpty [label="Positions empty?\nAbort with error"];
	 *     loadMap [label="Load or compute genetic map"];
	 *
	 *     chooseAlg [label="Choose chunking method"];
	 *     recursive [label="Recursive Split"];
	 *     sequential [label="Sequential Split\n(up to 10,000 tries)"];
	 *     fallback [label="Fallback to recursive\nif no sequential solution"];
	 *     tooSmall [label="If only one chunk:\nWarn user"];
	 *
	 *     uniformMode [label="Uniform number of variants mode?"];
	 *     uniformSeq [label="Sequential Split (Uniform mode)\n(up to 10,000 tries)"];
	 *
	 *     output [label="Write chunks to output file"];
	 *     end [label="End"];
	 *
	 *     start -> coords -> readData -> checkEmpty;
	 *     checkEmpty -> loadMap;
	 *     loadMap -> chooseAlg;
	 *     chooseAlg -> recursive [label="--recursive"];
	 *     chooseAlg -> sequential [label="--sequential"];
	 *     sequential -> fallback [label="No solution"];
	 *     sequential -> tooSmall [label="1 chunk only"];
	 *     recursive -> uniformMode;
	 *     fallback -> uniformMode;
	 *     tooSmall -> uniformMode;
	 *
	 *     uniformMode -> uniformSeq [label="Yes"];
	 *     uniformMode -> output [label="No"];
	 *     uniformSeq -> output;
	 *     output -> end;
	 * }
	 * @enddot
	 *
	 * @note
	 * - Requires `options` to be pre-populated via `parse_command_line()` before calling.
	 * - Uses `vrb` for logging and progress reporting.
	 *
	 * @warning
	 * - Empty position list after `readData()` causes immediate termination.
	 * - Sequential splitting may not converge within the retry limit.
	 *
	 * @complexity
	 * - Recursive: O(n log n) typical.
	 * - Sequential: O(n * attempts) worst case.
	 *
	 * @see chunker::split_recursive(), chunker::split_sequential(), chunker::add_buffer()
	 */
	void chunk();

	/**
	 * @brief Executes the chunking workflow.
	 *
	 * This function runs the chunking process by:
	 * 1. Declaring available command-line options (`declare_options()`).
	 * 2. Parsing the provided arguments (`parse_command_line()`).
	 * 3. Validating the parsed options (`check_options()`).
	 * 4. Printing verbose file information (`verbose_files()`).
	 * 5. Printing verbose option settings (`verbose_options()`).
	 * 6. Performing the actual chunking (`chunk()`).
	 *
	 * @param args Reference to a vector of strings containing 
	 *             the command-line arguments for the chunking process.
	 */
	void chunk(std::vector < std::string > &);

	/**
	 * @brief Builds the genomic coordinates for chunk processing.
	 *
	 * This method retrieves the genomic region specified in the command-line
	 * options, stores it in the `gregion` member variable, and parses it into
	 * a usable coordinate range for downstream chunking operations.
	 *
	 * The region string is extracted from the `--region` option, which is expected
	 * to follow the format:
	 * ```
	 * <chromosome>[:<start>-<end>]
	 * ```
	 * Examples:
	 * - `"chr1"` (entire chromosome 1)
	 * - `"chr1:1000000-2000000"` (from position 1,000,000 to 2,000,000 on chr1)
	 *
	 * @note This method must be called before any operation that requires
	 *       chromosome coordinates, as it initializes internal state.
	 *
	 * @pre The `options` map must contain a valid `"region"` key.
	 * @post `gregion` is populated and `parseRegion()` has been executed.
	 *
	 * @see parseRegion()
	 */
	void buildCoordinates();

	/**
	 * @brief Parses the genomic region string into chromosome, start, and stop coordinates.
	 *
	 * This function interprets the `gregion` string, expected to be in one of the formats:
	 * - "chrN" (whole chromosome)
	 * - "chrN:Y-" (chromosome from position Y to end)
	 * - "chrN:Y-Z" (chromosome from position Y to Z)
	 *
	 * The function extracts the chromosome ID and genomic coordinates, validating them
	 * and assigning default values if needed:
	 * - `start` is set to 0 if not specified.
	 * - `stop` is set to the maximum 32-bit integer value if not specified.
	 *
	 * If the whole chromosome is specified (no coordinates), the `whole_chr` flag is set.
	 *
	 * @throws Runtime error via `vrb.error` if the region string format is invalid, 
	 *         or if coordinates are out of logical bounds (e.g., negative start, start > stop).
	 *
	 * @pre `gregion` must be a non-empty string with a valid region specification.
	 * @post Member variables `start`, `stop`, and `whole_chr` are set accordingly.
	 */
	void parseRegion();

	/**
	 * @brief Retrieves all values associated with a given key from a multimap.
	 *
	 * This function searches the input `map_positions` multimap for all entries with the
	 * key equal to `pos` and collects their corresponding values into a vector.
	 *
	 * @param pos The key position to search for in the multimap.
	 * @param map_positions A multimap mapping keys of type `long int` to values of type `long int`.
	 *
	 * @return A vector of `long int` containing all values from `map_positions` that have the key `pos`.
	 *
	 * @note If no entries are found for the given key, the returned vector will be empty.
	 */
	std::vector < long int > getByPos(const long int, const std::multimap < long int, long int >& map_positions);
	
	/**
	 * @brief Sets the genetic map positions in centiMorgans (cM) for genomic markers.
	 *
	 * This function initializes and fills the `positions_cm` vector, which holds the
	 * genetic map positions in cM corresponding to the physical positions in `positions_mb`.
	 * It uses a `gmap_reader` object containing reference genetic map data and a mapping
	 * from physical to genetic positions.
	 *
	 * The process includes:
	 * - Setting exact cM positions from the genetic map using `setCentiMorgan`.
	 * - Interpolating cM values for positions not directly found in the genetic map using `interpolateCentiMorgan`.
	 *
	 * @param readerGM Reference to a genetic map reader containing base pair positions (`pos_bp`)
	 *                 and genetic map positions (`pos_cm`).
	 * @param map_positions A multimap associating physical positions (base pairs) to indices in the genetic map.
	 * @param positions_mb Vector of marker physical positions in megabases (Mb).
	 * @param[out] positions_cm Vector to be filled with corresponding genetic map positions in centiMorgans (cM).
	 *
	 * @throws Throws an error via `vrb.error` if no markers are found in the region (empty `positions_mb`).
	 *
	 * @note The function logs the number of positions directly set and the number of positions interpolated,
	 *       along with the elapsed time for the interpolation step.
	 */
	void setGeneticMap(const gmap_reader&, const std::multimap < long int, long int >&, const std::vector<long int>& positions_mb, std::vector<float>& positions_cm);
	
	/**
	 * @brief Sets a constant-rate genetic map when no genetic map file is provided.
	 *
	 * This function assigns genetic map positions (in centiMorgans, cM) to markers based on
	 * a constant recombination rate assumption of 1 cM per 1 Mb physical distance.
	 * It initializes the `positions_cm` vector with scaled values derived from physical
	 * positions in megabases (Mb), adjusting the baseline to zero.
	 *
	 * @param map_positions A multimap of physical positions to indices (not used in this implementation,
	 *                      but kept for interface consistency).
	 * @param positions_mb Vector of marker physical positions in base pairs.
	 * @param[out] positions_cm Vector to be filled with corresponding genetic map positions in centiMorgans (cM).
	 *
	 * @throws Throws an error via `vrb.error` if no markers are found in the region (empty `positions_mb`).
	 *
	 * @note This function logs that a constant map was used and reports the elapsed time.
	 */
	void setGeneticMap(const std::multimap < long int, long int >&, const std::vector<long int>& positions_mb, std::vector<float>& positions_cm);
	
	/**
	 * @brief Assigns genetic positions in centiMorgans (cM) to markers based on exact matches.
	 *
	 * This function updates the `positions_cm` vector by mapping physical base-pair positions
	 * (`pos_bp`) to their corresponding genetic positions (`pos_cM`) using the provided
	 * multimap `map_positions`. For each physical position in `pos_bp`, all corresponding
	 * marker indices from `map_positions` are assigned the genetic position from `pos_cM`.
	 *
	 * @param pos_bp Vector of physical base-pair positions from the genetic map.
	 * @param pos_cM Vector of corresponding genetic positions in centiMorgans (cM).
	 * @param map_positions Multimap linking physical positions (keys) to marker indices (values).
	 * @param[out] positions_cm Vector where genetic map positions (cM) for each marker will be stored.
	 *
	 * @return The total count of markers successfully assigned a genetic position.
	 *
	 * @note It is assumed that `pos_bp` and `pos_cM` have the same length.
	 */
	long int setCentiMorgan(const std::vector < long int > & pos_bp, const std::vector < double > & pos_cM, const std::multimap < long int, long int >& map_positions, std::vector<float> & positions_cm);
	
	/**
	 * @brief Interpolates genetic map positions (in centiMorgans) for markers lacking exact map positions.
	 * 
	 * This function fills in missing genetic positions in `positions_cm` by performing linear interpolation
	 * between known genetic map points (`pos_bp` and `pos_cM`). For markers outside the bounds of the
	 * known map, it extrapolates based on the average recombination rate.
	 * 
	 * The interpolation is done for markers whose genetic position is currently set to -1.
	 * 
	 * @param pos_bp Vector of known physical base-pair positions from the genetic map.
	 * @param pos_cM Vector of known genetic positions in centiMorgans corresponding to `pos_bp`.
	 * @param map_positions Multimap linking physical positions to marker indices (not directly used in interpolation).
	 * @param positions_mb Vector of physical positions of all markers to be assigned genetic positions.
	 * @param[out] positions_cm Vector of genetic map positions to be updated by interpolation.
	 *                        Positions initially set to -1 will be assigned interpolated values.
	 * 
	 * @return Number of markers for which the genetic map position was interpolated or extrapolated.
	 * 
	 * @note Assumes `pos_bp` and `pos_cM` are sorted in ascending order.
	 * @note Linear interpolation is used between the closest known genetic map positions.
	 *       Extrapolation outside known map bounds uses the average recombination rate.
	 */
	long int interpolateCentiMorgan(const std::vector < long int > & pos_bp, const std::vector < double > & pos_cM, const std::multimap < long int, long int >&, const std::vector<long int> & positions_mb, std::vector<float>& positions_cm);
	
	/**
	 * @brief Computes the length of the variant map region in base pairs.
	 * 
	 * This method calculates the physical length covered by the variant map
	 * as the difference between the base-pair position of the last variant and
	 * the first variant in `vec_pos`, plus one to account for inclusive coordinates.
	 * 
	 * @return The length of the region in base pairs (unsigned int).
	 * 
	 * @note Assumes `vec_pos` is a non-empty vector sorted by base-pair positions,
	 *       and that `vec_pos.front()` and `vec_pos.back()` are valid pointers to variant objects
	 *       with accessible `bp` (base-pair position) members.
	 */
	unsigned long int length() const;

	/**
	 * @brief Computes the genetic length of the variant map region in centiMorgans (cM).
	 * 
	 * This method calculates the genetic distance covered by the variant map
	 * as the difference between the centiMorgan (cM) position of the last variant
	 * and the first variant in `vec_pos`.
	 * 
	 * @return The genetic length of the region in centiMorgans (double).
	 * 
	 * @note Assumes `vec_pos` is a non-empty vector sorted by base-pair positions,
	 *       and that `vec_pos.front()` and `vec_pos.back()` are valid pointers to variant objects
	 *       with accessible `cm` (centiMorgan position) members.
	 */
	double lengthcM() const;

};


#endif


