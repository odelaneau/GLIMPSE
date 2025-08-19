/*******************************************************************************
 * @file checker.h
 * @brief Class definition for the GLIMPSE checker module.
 *
 * The `checker` class orchestrates the validation and summary steps for
 * genetic data analysis in GLIMPSE. It handles:
 * - Parsing and validating command-line options
 * - Reading and initializing input datasets
 * - Executing quality control checks
 * - Writing the processed results
 * 
 ******************************************************************************/

#ifndef _CHECKER_H
#define _CHECKER_H

#include "otools.h"
#include "call_set_header.h"

/**
 * @class checker
 * @brief Coordinates input parsing, validation, processing, and output
 *        for genotype checking workflows.
 *
 * This class acts as the main controller for the "checker" tool.
 * It contains:
 * - Command-line parsing and validation
 * - Management of the internal `call_set` object (`D`)
 * - Input/output handling
 *
 * **Typical Workflow:**
 * 1. `declare_options()` – Define allowed command-line arguments.
 * 2. `parse_command_line()` – Read and store arguments into `options`.
 * 3. `check_options()` – Ensure arguments are logically consistent.
 * 4. `read_files_and_initialise()` – Load datasets.
 * 5. Perform computations.
 * 6. `write_files_and_finalise()` – Save results.
 */
class checker {
public:
	/** @name Command-line options */
    ///@{
    bpo::options_description descriptions;  ///< CLI option descriptions.
    bpo::variables_map options;              ///< Parsed CLI option values.
    ///@}

    /** @name Internal data */
    ///@{
    call_set D;  ///< Main dataset container.
    ///@}

    /** @name Lifecycle */
    ///@{
    checker();   ///< Constructor: sets defaults and allocates structures.
    ~checker();  ///< Destructor: cleans up allocated resources.
    ///@}

	//PARAMETERS
	/**
	 * @brief Declares and configures all command-line options for the checker tool.
	 *
	 * This method uses Boost.Program_options (`bpo`) to define, group, and register 
	 * all the command-line options required by the application. The options are organized 
	 * into four categories: Basic options, Input parameters, Algorithm parameters, and 
	 * Output files. Once declared, all options are aggregated into the member 
	 * variable `descriptions` for later parsing.
	 *
	 * **Option Groups:**
	 * 
	 * - <b>Basic options (`opt_base`)</b>
	 *   - `--help` : Displays the help message.
	 *   - `--seed` : Integer seed for the random number generator (default: `15052011`).
	 *   - `--threads` : Number of threads to use (default: `1`).
	 *
	 * - <b>Input parameters (`opt_input`)</b>
	 *   - `--input` : Path to file listing regions, frequencies, validation, and imputed dataset.
	 *                 Multiple lines may be used for genome-wide concordance (per chromosome).
	 *   - `--samples` : File containing one sample ID per line to process.
	 *   - `--gt-val` : Use hard-called genotypes from `FORMAT/GT` in the validation dataset
	 *                  instead of phred-scaled likelihoods.
	 *   - `--gt-tar` : Use `FORMAT/GT` for best-guess genotypes in the target dataset instead 
	 *                  of `FORMAT/GP` (default). Still requires `FORMAT/DS` and `FORMAT/GP` 
	 *                  for calibration and r² calculations.
	 *
	 * - <b>Algorithm parameters (`opt_algo`)</b>
	 *   - `--af-tag` : INFO tag name for allele frequency binning (default: `"AF"`).
	 *   - `--use-alt-af` : Interpret allele frequencies as ALT allele frequency `[0,1]` 
	 *                      instead of minor allele frequency `[0,0.5]`.
	 *   - `--bins` : Custom allele frequency bins for r² computation (MAF range `[0-0.5]` by default,
	 *                `[0-1]` if `--use-ref-alt` is specified).
	 *   - `--ac-bins` : Custom allele count bins for r² computation.
	 *   - `--allele-counts` : Use default allele count bins (requires `AN` field in frequency file).
	 *   - `--min-val-gl` : Minimum genotype likelihood probability P(G|R) for validation data 
	 *                      (0 disables filtering; applies to `--gt-validation`).
	 *   - `--min-val-dp` : Minimum coverage in validation data (0 disables filtering; 
	 *                      program exits if `FORMAT/DP` missing and min DP > 0).
	 *   - `--min-tar-gp` : Minimum GP probabilities for filtering target genotypes. Uses `FORMAT/PL` 
	 *                      if `--gt-tar` is specified. Empty disables filtering.
	 *   - `--out-r2-per-site` : Output r² value for each site.
	 *   - `--out-rej-sites` : Output sites rejected from concordance analysis.
	 *   - `--out-conc-sites` : Output sites with complete concordance between target and truth genotypes.
	 *   - `--out-disc-sites` : Output sites with at least one discordant genotype.
	 *   - `--groups` : Path to file defining user-defined grouping bins (alternative to `--bins`).
	 *
	 * - <b>Output files (`opt_output`)</b>
	 *   - `--output` / `-O` : Prefix for all output files (extensions added automatically).
	 *   - `--log` : Path to log file.
	 *
	 * @note This function only declares available options; it does not parse them.
	 *       Parsing and handling of options must be done elsewhere in the program.
	 *
	 * @see boost::program_options::options_description
	 */
	void declare_options();

	/**
	 * @brief Parses and processes the command-line arguments for the checker tool.
	 *
	 * This method takes a list of command-line arguments, parses them using 
	 * Boost.Program_options (`bpo`), and stores the resulting values in the 
	 * internal `options` variable. It also initializes the log file (if specified), 
	 * and prints program metadata such as title, authors, version, and citation 
	 * information.
	 *
	 * **Main Steps:**
	 * 1. **Parsing arguments:**
	 *    - Uses `bpo::command_line_parser` with the previously declared `descriptions`
	 *      to parse options from the given `args` vector.
	 *    - Stores parsed options in `options` and applies them with `bpo::notify`.
	 *    - If a parsing error occurs, prints the error message and exits.
	 *
	 * 2. **Log file initialization:**
	 *    - If the `--log` option is provided, attempts to open the specified log file 
	 *      via `vrb.open_log()`.
	 *    - If the log file cannot be created, calls `vrb.error()` and terminates execution.
	 *
	 * 3. **Metadata output:**
	 *    - Prints program title, authors, contact information, version, commit ID/date, 
	 *      and relevant citations to the log and/or standard output.
	 *    - Displays the current run date via `tac.date()`.
	 *
	 * 4. **Help message:**
	 *    - If the `--help` option is provided, prints all option descriptions 
	 *      (`descriptions`) to standard output and exits.
	 *
	 * @param args Reference to a vector of strings containing the command-line arguments.
	 *
	 * @note This function assumes that `declare_options()` has already been called 
	 *       to populate `descriptions` with the full set of available options.
	 * @note If `--help` is specified, no further processing occurs after displaying the help message.
	 *
	 * @throws std::exit The function terminates the program if parsing fails, 
	 *         the log file cannot be created, or help is requested.
	 *
	 * @see checker::declare_options()
	 * @see boost::program_options::command_line_parser
	 */
	void parse_command_line(std::vector < std::string > &);

	/**
	 * @brief Validates the command-line options after parsing.
	 *
	 * This method checks the presence, mutual exclusivity, and value constraints 
	 * of the options parsed into the `options` variable. If any requirement is 
	 * not met, it calls `vrb.error()` with a descriptive message and terminates 
	 * execution.
	 *
	 * **Validation Rules:**
	 * 1. **Required options:**
	 *    - `--input` must be provided.
	 *    - `--output` must be provided.
	 *
	 * 2. **Mutually exclusive binning parameters:**
	 *    - Exactly one of the following must be specified:
	 *      - `--bins`
	 *      - `--ac-bins`
	 *      - `--allele-counts`
	 *      - `--groups`
	 *
	 * 3. **Thread count:**
	 *    - `--threads` must be a strictly positive integer (`>= 1`).
	 *
	 * 4. **Genotype filtering constraints:**
	 *    - If either `--min-val-gl` or `--min-val-dp` is set, `--gt-val` **must not** be set.
	 *      This ensures that genotype-based validation filters are not mixed with 
	 *      hard-called genotype validation mode.
	 *
	 * 5. **Validation dataset filtering:**
	 *    - If `--gt-val` is **not** set:
	 *      - `--min-val-gl` must be provided.
	 *      - `--min-val-dp` must be provided.
	 *    - This ensures explicit filtering parameters are set for non-GT validation mode.
	 *
	 * @note The commented-out check between `--gt-tar` and `--min-tar-gp` suggests 
	 *       an intended restriction that is currently disabled.
	 *
	 * @throws std::exit The program terminates if any of the above checks fail.
	 *
	 * @see checker::parse_command_line()
	 * @see checker::declare_options()
	 */
	void check_options();

	/**
	 * @brief Prints the parsed and validated command-line options in a human-readable format.
	 *
	 * This method outputs the current configuration parameters stored in the `options` 
	 * variable to the verbose logger (`vrb`). It is intended to provide users with a 
	 * summary of key parameters before the main computation begins.
	 *
	 * **Displayed Information:**
	 * 1. **Basic parameters:**
	 *    - Output prefix (`--output`)
	 *    - Random seed (`--seed`)
	 *    - Number of threads (`--threads`)
	 *
	 * 2. **Validation dataset filters:**
	 *    - If `--gt-val` is **not** set:
	 *      - Minimum validation genotype likelihood (`--min-val-gl`)
	 *      - Minimum validation read depth (`--min-val-dp`)
	 *    - If `--gt-val` **is** set:
	 *      - Message indicating all validation genotypes are used.
	 *
	 * 3. **Allele frequency binning method:**
	 *    - INFO tag used for allele frequency (`--af-tag`, default `"AF"`).
	 *
	 * 4. **Binning/grouping details:**
	 *    - If `--bins` is set: number of allele frequency bins.
	 *    - If `--ac-bins` is set: number of allele count bins.
	 *    - If `--groups` is set: indicates user-defined groups.
	 *    - Otherwise: indicates fixed allele frequency bins.
	 *
	 * 5. **Target genotype filtering:**
	 *    - If `--min-tar-gp` is set and contains values: number of GP thresholds.
	 *
	 * 6. **Output site-level metrics:**
	 *    - Whether to output:
	 *      - r² per site (`--out-r2-per-site`)
	 *      - rejected sites (`--out-rej-sites`)
	 *      - fully concordant sites (`--out-conc-sites`)
	 *      - discordant sites (`--out-disc-sites`)
	 *    - Shown as `"YES"` or `"NO"`.
	 *
	 * @note This method does not modify any program state—it only reports the current settings.
	 * @note Assumes that `parse_command_line()` and `check_options()` have already been executed successfully.
	 *
	 * @see checker::parse_command_line()
	 * @see checker::check_options()
	 */
	void verbose_options();

	/**
	 * @brief Prints information about input and auxiliary files specified in the command-line options.
	 *
	 * This method logs the paths to key files involved in the analysis using the 
	 * verbose logger (`vrb`). It provides a quick summary of which files are being 
	 * used for the run.
	 *
	 * **Displayed Information:**
	 * 1. **Primary input file:**
	 *    - The file specified by `--input` (lists regions, frequencies, validation, and imputed datasets).
	 *
	 * 2. **Log file (optional):**
	 *    - If `--log` is specified, prints the path to the log output file.
	 *
	 * 3. **Groups file (optional):**
	 *    - If `--groups` is specified, prints the path to the file containing user-defined group bins.
	 *
	 * @note This function is for reporting purposes only; it does not check whether 
	 *       the specified files exist or are readable.
	 * @note Should be called after `parse_command_line()` and `check_options()` to 
	 *       ensure the options have been validated.
	 *
	 * @see checker::verbose_options()
	 */
	void verbose_files();

	/**
	 * @brief Reads required input files, initializes parameters, and runs the concordance analysis.
	 *
	 * This method is the main initialization and execution entry point after command-line 
	 * parsing and option validation. It sets up the random number generator, determines 
	 * binning/grouping configurations, loads sample and dataset files, validates GP filters, 
	 * and triggers the concordance computation through the `D` object.
	 *
	 * **Execution Steps:**
	 * 
	 * 1. **Initialize RNG and threading:**
	 *    - Sets the random seed from `--seed`.
	 *
	 * 2. **Process bin/group configuration and validation thresholds:**
	 *    - Reads `--min-val-gl` and `--min-val-dp` if present.
	 *    - Initializes `D` using:
	 *      - Allele frequency bins from `--bins`, **or**
	 *      - Group definitions from `--groups`, **or**
	 *      - Default bins if neither is specified.
	 *
	 * 3. **Load target samples:**
	 *    - If `--samples` is specified, calls `D.setTargets()` with the provided sample list file.
	 *
	 * 4. **Read the input file list:**
	 *    - Opens the file specified by `--input`.
	 *    - Expects each line to contain **four columns**:
	 *      1. Regions file
	 *      2. Frequencies file
	 *      3. Validation dataset
	 *      4. Imputed dataset
	 *    - Populates `vec_regi`, `vec_freq`, `vec_true`, and `vec_esti` with these values.
	 *    - Throws an error if no valid sets are found.
	 *
	 * 5. **Run concordance analysis:**
	 *    - If `--min-tar-gp` is set:
	 *      - Loops over each GP filter value (valid range `[0,1]`).
	 *      - For each filter:
	 *        - Calls `D.readData()` with the given filter.
	 *        - Writes results with `write_files_and_finalise()` (suffix `_GPfilt_<value>`).
	 *    - If `--min-tar-gp` is not set:
	 *      - Runs concordance without GP filtering.
	 *
	 * 6. **Report execution time:**
	 *    - Logs running time for each concordance run and total elapsed time.
	 *
	 * @note This function assumes that:
	 *       - `parse_command_line()` and `check_options()` have already validated options.
	 *       - Input files are correctly formatted and accessible.
	 * @note All output is reported via the verbose logger (`vrb`).
	 *
	 * @throws vrb.error() Terminates execution if:
	 *         - Input files are missing or unreadable.
	 *         - GP filter values are outside `[0,1]`.
	 *         - No dataset sets are detected in `--input`.
	 *
	 * @see checker::check_options()
	 * @see checker::write_files_and_finalise()
	 */
	void read_files_and_initialise();

	/**
	 * @brief Main entry point for executing the checker tool.
	 *
	 * This method coordinates the complete workflow for running the concordance analysis,
	 * from declaring available options to reading input files and producing results.
	 * It sequentially calls all the major steps required for initialization and execution.
	 *
	 * **Execution Flow:**
	 * 1. **Declare options** — Calls `declare_options()` to define all available
	 *    command-line arguments.
	 * 2. **Parse command line** — Calls `parse_command_line()` to read arguments
	 *    from `args` and store them in the internal `options` container.
	 * 3. **Validate options** — Calls `check_options()` to ensure that all required
	 *    parameters are present and consistent.
	 * 4. **Report files** — Calls `verbose_files()` to log the paths of key input/output files.
	 * 5. **Report parameters** — Calls `verbose_options()` to display the current
	 *    run configuration.
	 * 6. **Initialize and run analysis** — Calls `read_files_and_initialise()` to load data,
	 *    apply filters, and execute the concordance computation.
	 *
	 * @param args Reference to a vector of strings containing the command-line arguments.
	 *
	 * @note This method acts as the high-level controller for the checker’s execution.
	 *       It assumes no prior initialization other than passing in the argument list.
	 *
	 * @see checker::declare_options()
	 * @see checker::parse_command_line()
	 * @see checker::check_options()
	 * @see checker::verbose_files()
	 * @see checker::verbose_options()
	 * @see checker::read_files_and_initialise()
	 */
	void check(std::vector < std::string > &);

	/**
	 * @brief Finalizes the analysis by writing output files.
	 *
	 * This method performs the final output step of the checker workflow.
	 * It writes the processed results to the file specified in the `--output`
	 * option, optionally appending a genotype probability filter suffix to
	 * the file name.
	 *
	 * **Execution Steps:**
	 * 1. Construct the output file path by concatenating:
	 *    - The base path from `options["output"]`
	 *    - The `gp_filter` string (if provided)
	 * 2. Call `D.writeData()` to write the results to disk.
	 *
	 * @param gp_filter Optional string appended to the output filename
	 *        to indicate genotype probability filter settings.
	 *
	 * @note This function assumes that all necessary data structures have
	 *       been populated prior to calling.
	 *
	 * @see checker::check()
	 * @see Data::writeData()
	 */
	void write_files_and_finalise(const std::string gp_filter = "");
};


#endif


