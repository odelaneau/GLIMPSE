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

#ifndef _CALLSET_H
#define _CALLSET_H

#include "otools.h"

#define N_BIN_CAL 100

class call_set {
public:
	// Sample IDs [N] & Bins [L]
	int N, L, D;
	float T;
	bool use_subset_samples;
	std::set < std::string > subset_samples_set;
	std::vector < std::string > subset_samples;

	std::vector < std::string > samples;
	std::vector < double > bins;

	int fploidy;
	std::vector< int > ploidy_samples;
	std::vector< int > ind2gpos;
	int n_haploid = 0, 	n_diploid = 0, ploidy = 0;

	// Per sample concordance [3xN]
	std::vector < unsigned long int > genotype_spl_errors_all;
	std::vector < unsigned long int > genotype_spl_totals_all;
	std::vector < unsigned long int > genotype_val_spl_totals_all;

	std::vector < unsigned long int > genotype_spl_errors_snps;
	std::vector < unsigned long int > genotype_spl_totals_snps;
	std::vector < unsigned long int > genotype_val_spl_totals_snps;

	std::vector < unsigned long int > genotype_spl_errors_indels;
	std::vector < unsigned long int > genotype_spl_totals_indels;
	std::vector < unsigned long int > genotype_val_spl_totals_indels;

	std::vector < unsigned long int > filtered_gp_spl_all;
	std::vector < unsigned long int > filtered_gp_spl_snps;
	std::vector < unsigned long int > filtered_gp_spl_indels;

	// Per bin concordance [3xL]
	std::vector < unsigned long int > genotype_bin_errors_all;
	std::vector < unsigned long int > genotype_bin_totals_all;
	std::vector < unsigned long int > genotype_val_bin_totals_all;

	std::vector < unsigned long int > genotype_bin_errors_snps;
	std::vector < unsigned long int > genotype_bin_totals_snps;
	std::vector < unsigned long int > genotype_val_bin_totals_snps;

	std::vector < unsigned long int > genotype_bin_errors_indels;
	std::vector < unsigned long int > genotype_bin_totals_indels;
	std::vector < unsigned long int > genotype_val_bin_totals_indels;

	std::vector < unsigned long int > filtered_gp_bin_all;
	std::vector < unsigned long int > filtered_gp_bin_snps;
	std::vector < unsigned long int > filtered_gp_bin_indels;

	// Concordance for calibration [100]
	std::vector < unsigned long int > genotype_cal_errors_all;
	std::vector < unsigned long int > genotype_cal_totals_all;

	std::vector < unsigned long int > genotype_cal_errors_snps;
	std::vector < unsigned long int > genotype_cal_totals_snps;

	std::vector < unsigned long int > genotype_cal_errors_indels;
	std::vector < unsigned long int > genotype_cal_totals_indels;

	// R2 per bin
	std::vector < std::string > rsquared_str;

	std::vector < stats2D > rsquared_bin_ds_all;
	std::vector < stats2D > rsquared_bin_gt_all;
	std::vector < stats1D > frequency_bin_all;

	std::vector < stats2D > rsquared_bin_ds_snps;
	std::vector < stats2D > rsquared_bin_gt_snps;
	std::vector < stats1D > frequency_bin_snps;

	std::vector < stats2D > rsquared_bin_ds_indels;
	std::vector < stats2D > rsquared_bin_gt_indels;
	std::vector < stats1D > frequency_bin_indels;

	// R2 per group
	int n_fields_in_group_files;
	std::map < std::string, std::pair < int, bool > > site2grp;

	// R2 per sample: from DS and GT
	std::vector < stats2D > rsquared_spl_ds_all;
	std::vector < stats2D > rsquared_spl_gt_all;

	std::vector < stats2D > rsquared_spl_ds_snps;
	std::vector < stats2D > rsquared_spl_gt_snps;

	std::vector < stats2D > rsquared_spl_ds_indels;
	std::vector < stats2D > rsquared_spl_gt_indels;

	//
	call_set ();
	~call_set();

	/**
	 * @brief Write a single VCF/BCF record to the output file.
	 *
	 * This method attempts to write the given `line` (variant record) to the
	 * provided HTSlib output file handle (`out_fd`) using the given output
	 * header (`out_hdr`).
	 *
	 * @param[out] out_fd   Pointer to an open HTSlib file object for writing.
	 * @param[out] out_hdr  Header structure corresponding to the output file.
	 * @param[in]  hdr_in   Header of the input file (currently unused).
	 * @param[in]  line     The variant record to be written.
	 *
	 * @details
	 * - Uses `bcf_write()` from HTSlib to serialize the record.
	 * - If writing fails, an error is thrown via `vrb.error()`.
	 * - The commented-out `bcf_update_format()` line suggests that genotype
	 *   format fields could be stripped or updated before writing.
	 *
	 * @throws std::runtime_error If the record cannot be written successfully.
	 *
	 * @note
	 * - `hdr_in` is not used in this function, but may be relevant for
	 *   transformations in other contexts.
	 * - This method assumes that the `line` record is already validated and
	 *   consistent with `out_hdr`.
	 */
	void write_record(htsFile *out_fd, bcf_hdr_t * out_hdr, bcf_hdr_t * hdr_in, bcf1_t *line);

	/**
	 * @brief Initialize the call set processing engine with allele frequency bins and validation thresholds.
	 *
	 * This method sets up the internal state of the `call_set` object based on user-specified
	 * allele frequency bin boundaries, genotype likelihood threshold (`T`), and read depth threshold (`D`).
	 *
	 * @param _maf_bins Vector of minor allele frequency (MAF) bin boundaries in ascending order.
	 *                  Must contain at least two elements (lower and upper bounds).
	 * @param _T        Minimum genotype likelihood probability for validation (GL ≥ T).
	 * @param _D        Minimum read depth for validation (Depth ≥ D).
	 *
	 * @details
	 * - Disables sample subset usage (`use_subset_samples = false`).
	 * - Stores the MAF bins and computes `L` = number of bins - 1.
	 * - Clears `site2grp` and `rsquared_str` containers to reset previous state.
	 * - Logs initialization details, including thresholds and bin configuration.
	 *
	 * Example:
	 * @code
	 * std::vector<double> bins = {0.0, 0.05, 0.1, 0.5};
	 * call_set cs;
	 * cs.initialize(bins, 0.9, 10);
	 * @endcode
	 *
	 * This example sets:
	 * - 3 MAF bins: [0.0–0.05], [0.05–0.1], [0.1–0.5]
	 * - GL threshold = 0.9
	 * - Depth threshold = 10
	 *
	 * @note
	 * The function assumes that `_maf_bins` is sorted and contains valid frequency values
	 * between 0 and 1.
	 */
	void initialize(std::vector < double >,double, int);

	/**
	 * @brief Initialize the call set processing engine with validation thresholds only.
	 *
	 * This overload configures the internal state using only the genotype likelihood
	 * threshold (`T`) and read depth threshold (`D`), without setting allele frequency bins.
	 *
	 * @param _T  Minimum genotype likelihood probability for validation (GL ≥ T).
	 * @param _D  Minimum read depth for validation (Depth ≥ D).
	 *
	 * @details
	 * - Disables sample subset usage (`use_subset_samples = false`).
	 * - Stores the thresholds in `T` and `D`.
	 * - Clears `site2grp` and `rsquared_str` to reset previous state.
	 * - Logs initialization details for thresholds.
	 *
	 * Example:
	 * @code
	 * call_set cs;
	 * cs.initialize(0.9, 10);
	 * @endcode
	 *
	 * This example sets:
	 * - GL threshold = 0.9
	 * - Depth threshold = 10
	 *
	 * @note
	 * Use this version when allele frequency binning is not required
	 * or is handled separately.
	 */
	void initialize(double, int); //allele-counts

	/**
	 * @brief Initialize the call set processing engine based on site-to-group mapping.
	 *
	 * This overload configures thresholds (`T`, `D`) and reads a group file (`fgrps`)
	 * mapping genomic sites to group IDs for downstream statistics (e.g., r² per group).
	 *
	 * @param fgrps Path to the group file.
	 *   - Format must be either:
	 *     - **3 columns**: `CHR POS GroupID`
	 *     - **5 columns**: `CHR POS REF ALT GroupID`
	 * @param _T  Minimum genotype likelihood probability for validation (GL ≥ T).
	 * @param _D  Minimum read depth for validation (Depth ≥ D).
	 *
	 * @details
	 * - Disables sample subset usage (`use_subset_samples = false`).
	 * - Clears and reinitializes `site2grp`, `rsquared_str`, and `bins`.
	 * - Parses each line of the group file:
	 *   1. Generates a unique site identifier (`uuid`) from either 2 or 4 site fields.
	 *   2. Finds or assigns a numeric group index for the `GroupID`.
	 *   3. Stores mapping in `site2grp` as `(group_index, false)`.
	 * - Verifies consistent column formatting throughout the file.
	 * - Logs summary information: thresholds, number of sites, and number of groups.
	 *
	 * @throws std::runtime_error If the group file format is invalid or inconsistent.
	 *
	 * Example:
	 * @code
	 * call_set cs;
	 * cs.initialize("groups.txt", 0.95, 8);
	 * @endcode
	 *
	 * This example:
	 * - Reads group mappings from `groups.txt`.
	 * - Sets GL threshold to 0.95 and depth threshold to 8.
	 *
	 * @note
	 * The `GroupID` strings are stored in `rsquared_str` in the order encountered.
	 * Group indices are assigned sequentially starting from 0.
	 */
	void initialize(std::vector < int >,double, int); //bins

	/**
	 * @brief Initializes the call set based on a group definition file.
	 *
	 * This method configures the call set to validate and group variant sites
	 * according to an external group file (`fgrps`). Each site in the file
	 * is assigned to a group index, and validation thresholds for genotype
	 * likelihood probabilities (`T`) and read depth (`D`) are set.
	 *
	 * The group file must be tab-delimited and contain either:
	 * - **3 columns**: CHR, POS, GroupID
	 * - **5 columns**: CHR, POS, REF, ALT, GroupID
	 *
	 * Sites are identified by a unique key constructed from the provided fields.
	 * All sites are stored in `site2grp`, with their corresponding group index
	 * recorded in `rsquared_str`. If the group already exists, the site is mapped
	 * to the existing group index; otherwise, a new group index is created.
	 *
	 * @param fgrps Path to the group definition file.
	 * @param _T Minimum genotype likelihood probability for validation.
	 * @param _D Minimum read depth for validation.
	 *
	 * @throws std::runtime_error If the group file format is invalid or empty.
	 *
	 * **Example group file (3 columns):**
	 * @code
	 * chr1    10583   Group1
	 * chr1    10611   Group1
	 * chr1    13302   Group2
	 * @endcode
	 *
	 * **Example group file (5 columns):**
	 * @code
	 * chr1    10583   G   A   Group1
	 * chr1    10611   C   T   Group1
	 * chr1    13302   T   G   Group2
	 * @endcode
	 */
	void initialize(std::string, double, int); //groups

	/**
	 * @brief Reads and stores a subset of target samples for analysis.
	 *
	 * This method loads a list of sample IDs from a file and restricts the
	 * analysis to only those samples. Each sample ID is stored in both a
	 * `std::set` (for quick lookup) and a `std::vector` (to preserve input order).
	 *
	 * **File format:**
	 * - One sample ID per line.
	 * - Extra whitespace after the sample ID is ignored.
	 * - Lines must contain at least one token; empty lines cause an error.
	 *
	 * @param fsamples Path to the file containing sample IDs to be included.
	 *
	 * @note
	 * - If the file contains duplicate sample IDs, only the first occurrence
	 *   is added to the `subset_samples` list; duplicates are ignored.
	 * - Sets `use_subset_samples` to `true`.
	 *
	 * @throws std::runtime_error If:
	 * - The file contains empty lines.
	 * - No samples are found in the file.
	 *
	 * **Example sample file:**
	 * @code
	 * HG00096
	 * HG00097
	 * HG00099
	 * @endcode
	 *
	 * **Side effects:**
	 * - Updates `subset_samples_set` and `subset_samples`.
	 * - Prints the number of loaded samples.
	 */
	void setTargets(std::string fsamples);

	/**
	 * @brief Determines the most likely genotype truth based on two probability-like values.
	 *
	 * This method takes two values (`pl0` and `pl1`) representing unnormalized likelihoods or probabilities
	 * of two genotype states (e.g., homozygous reference vs. heterozygous). It normalizes them and
	 * applies a threshold `T` to decide if the truth call is confident.
	 *
	 * @param pl0 Unnormalized likelihood/probability for genotype 0.
	 * @param pl1 Unnormalized likelihood/probability for genotype 1.
	 *
	 * @return
	 * - Returns `0` if genotype 0 is confidently more likely.
	 * - Returns `1` if genotype 1 is confidently more likely.
	 * - Returns `-1` if both genotypes have equal normalized probability.
	 * - Returns `-3` if any input likelihood is negative (invalid).
	 * - Returns `-4` if neither genotype passes the confidence threshold `T`.
	 *
	 * @note The confidence threshold `T` is assumed to be a member variable of the class.
	 */
	int getTruth(float, float);

	/**
	 * @brief Computes the calibration bin index based on genotype probabilities.
	 *
	 * This function calculates which calibration bin a genotype probability
	 * falls into by taking the maximum probability between `gp0` and `gp1`,
	 * scaling it to the range `[0, N_BIN_CAL-1]`, and truncating to an integer.
	 *
	 * @param gp0 Probability of genotype 0.
	 * @param gp1 Probability of genotype 1.
	 *
	 * @return The calibration bin index in the range `[0, N_BIN_CAL-1]`.
	 *
	 * @note Assumes `N_BIN_CAL` is a class constant or macro defining the number of calibration bins.
	 */
	int getCalibrationBin(float , float );

	/**
	 * @brief Finds the allele frequency bin index for a given probability.
	 *
	 * This function searches through the `bins` vector (assumed sorted in ascending order)
	 * to find the bin that contains the input `prob`. It returns the index of the bin where
	 * `prob` lies between `bins[b-1]` (exclusive) and `bins[b]` (inclusive).
	 *
	 * @param prob Allele frequency probability to bin.
	 *
	 * @return
	 * - The index of the bin containing `prob`.
	 * - `-1` if `prob` does not fall into any bin.
	 *
	 * @note `bins` should contain at least two values defining the bin edges.
	 */
	int getFrequencyBin(float);

	/**
	 * @brief Determines the most likely genotype truth for a multi-allelic site.
	 *
	 * This method evaluates genotype likelihood probabilities for three possible genotypes
	 * (e.g., diploid calls: 0, 1, 2) along with sequencing depth (`dp`) and ploidy information.
	 * It returns the most probable genotype if confident, or error codes for various edge cases.
	 *
	 * @param pl0 Likelihood/probability of genotype 0.
	 * @param pl1 Likelihood/probability of genotype 1.
	 * @param pl2 Likelihood/probability of genotype 2.
	 * @param dp Sequencing depth at the site. Compared to threshold `D`.
	 * @param ploidy Ploidy of the site (typically 1 or 2).
	 *
	 * @return
	 * - `0`, `1`, or `2` for the most confident genotype call.
	 * - `-1` if no genotype is clearly most likely (tie or uncertainty).
	 * - `-2` if depth is below threshold `D`.
	 * - `-3` if any likelihood is negative (invalid).
	 * - `-4` if all genotype likelihoods are below confidence threshold `T`.
	 * - `-5` if depth is missing (`bcf_int32_missing`).
	 *
	 * @note
	 * - For haploid sites (`ploidy == 1`), this function delegates to the 2-genotype version.
	 * - The thresholds `T` (confidence) and `D` (depth) are class members.
	 */
	int getTruth(float, float, float, int, int);

	/**
	 * @brief Determines the most likely genotype based on genotype probabilities and a filter.
	 *
	 * This function selects the most likely genotype (0 or 1) from two genotype probabilities (`gp0`, `gp1`),
	 * subject to a genotype probability filter (`gp_filter`) and a dosage value (`ds`).
	 *
	 * @param gp0 Probability of genotype 0.
	 * @param gp1 Probability of genotype 1.
	 * @param ds Dosage value, expected to be in the range [0,1].
	 * @param gp_filter Minimum genotype probability threshold required to make a call.
	 *
	 * @return
	 * - `0` if genotype 0 is most likely and passes the filter.
	 * - `1` if genotype 1 is most likely and passes the filter.
	 * - `-1` if probabilities are equal (tie).
	 * - `-2` if dosage `ds` is out of the expected range [0,1].
	 * - `-3` if either genotype probability is out of the range [0,1].
	 * - `-4` if the maximum genotype probability is below the filter threshold.
	 */
	int getMostLikely(const float , const float , const float, const float);

	/**
	 * @brief Determines the most likely genotype among three possibilities, with filtering.
	 *
	 * This function selects the most likely genotype (0, 1, or 2) from three genotype probabilities (`gp0`, `gp1`, `gp2`),
	 * considering the ploidy, dosage (`ds`), and a genotype probability filter (`gp_filter`).
	 *
	 * @param gp0 Probability of genotype 0.
	 * @param gp1 Probability of genotype 1.
	 * @param gp2 Probability of genotype 2.
	 * @param ds Dosage value, expected to be in the range [0, ploidy].
	 * @param ploidy Ploidy of the genotype (1 or 2).
	 * @param gp_filter Minimum genotype probability threshold required to make a call.
	 *
	 * @return
	 * - `0`, `1`, or `2` for the most likely genotype passing the filter.
	 * - `-1` if no genotype is clearly most likely (currently not implemented, see TODO).
	 * - `-2` if dosage `ds` is out of the valid range [0, ploidy].
	 * - `-3` if any genotype probability is out of the range [0,1].
	 * - `-4` if the maximum genotype probability is below the filter threshold.
	 *
	 * @note
	 * - For haploid genotypes (`ploidy == 1`), this method delegates to the 2-genotype version.
	 * - The TODO comment indicates returning `-1` for ties is planned but currently returns 0 as a placeholder.
	 */
	int getMostLikely(const float , const float , const float , const float, const int, const float);

	/**
	 * @brief Validates and returns the most likely genotype integer.
	 *
	 * This function checks that the dosage (`ds`) is within the expected range for the given ploidy,
	 * and that the genotype (`gt`) is non-negative. If validations fail, error codes are returned.
	 *
	 * @param gt The genotype call as an integer.
	 * @param ds Dosage value, expected to be in the range [0, ploidy].
	 * @param ploidy Ploidy of the genotype.
	 *
	 * @return
	 * - `gt` if valid.
	 * - `-2` if dosage `ds` is out of the valid range [0, ploidy].
	 * - `-3` if genotype `gt` is negative (invalid).
	 *
	 * @note The dosage check condition likely needs fixing: currently uses `&&` which will never be true.
	 *       Should probably be `if (ds < 0.0f || ds > ((float) ploidy) + 1e-7)`.
	 */
	int getMostLikelyGT(int, float, int);

	/**
	 * @brief Calculates the calibration bin index based on the maximum genotype probability.
	 *
	 * This method computes the calibration bin by finding the maximum genotype probability
	 * among the three provided probabilities (`gp0`, `gp1`, `gp2`) and scaling it to a bin index.
	 * The calibration bins are used to group genotype calls by confidence level.
	 *
	 * @param gp0 Probability of genotype 0.
	 * @param gp1 Probability of genotype 1.
	 * @param gp2 Probability of genotype 2.
	 *
	 * @return The calibration bin index as an integer, ranging from 0 to `N_BIN_CAL-1`.
	 *
	 * @note This method should not be affected by ploidy; for haploid (ploidy=1), `gp1` equals `gp2`.
	 */
	int getCalibrationBin(float , float , float );

	/**
	 * @brief Reads and processes variant call data for concordance analysis.
	 *
	 * This function manages the entire data flow of reading and processing variant data:
	 *  - Loads truth genotype call files, estimated genotype call files, and allele frequency files.
	 *  - Intersects sample lists across truth and estimated datasets to identify common samples.
	 *  - Reads variant sites and associated genotype likelihoods or hard calls, depending on options.
	 *  - Applies filters on genotype probabilities, read depths, and other quality metrics.
	 *  - Assigns variants to allele frequency bins or user-defined groups.
	 *  - Computes concordance metrics such as dosage r² between truth and estimated genotypes.
	 *  - Supports multi-threaded processing based on configured thread counts.
	 *  - Outputs detailed concordance statistics, optionally per-site or grouped.
	 *  - Handles logging and error reporting through a verbosity interface.
	 *
	 * @param ftruth Vector of file paths containing truth genotype data.
	 * @param festimated Vector of file paths containing imputed/estimated genotype data.
	 * @param ffrequencies Vector of file paths with allele frequency information for variants.
	 * @param region Vector of genomic regions to restrict analysis (e.g., chromosome segments).
	 * @param options Map of command-line options configuring filters, output, and processing behaviors.
	 * @param gp_filter Threshold on genotype probability below which genotypes are ignored.
	 * @param out_filename Prefix for output files where results and logs will be written.
	 *
	 * @note The function expects input files to be properly formatted VCF/BCF files or equivalent,
	 *       and that sample IDs can be matched between truth and estimated datasets.
	 *       It performs internal consistency checks and will throw errors if key files or samples are missing.
	 */
	void readData(std::vector < std::string > &, std::vector < std::string > &, std::vector < std::string > &, std::vector < std::string > &, bpo::variables_map&, const float gp_filter, const std::string out_filename);
	
	/**
	 * @brief Writes genotype concordance statistics to output files.
	 *
	 * This function generates and writes detailed genotype concordance statistics to output files.
	 * It computes various metrics such as genotype concordance by sample, error rates, and correlation coefficients.
	 * The data is organized and written in a structured format for further analysis.
	 *
	 * @param fout The base filename for the output files. The function will append appropriate extensions
	 *             to create the final output filenames.
	 *
	 * @throws std::ios_base::failure If there is an error opening or writing to the output files.
	 * @throws std::invalid_argument If the provided filename is empty or invalid.
	 *
	 * @note The function assumes that the necessary data structures and variables are properly initialized
	 *       and populated before calling this function. It also assumes that the input data is valid and
	 *       consistent.
	 */
	void writeData(std::string);

	void computeRsquaredPerBin(std::string output);
	void computeRsquaredPerBinPerSample(std::string output);

	void computeConcordancePerBIN(std::string output);
	void concordanceOverall(std::string output);
	void concordancePerIndividual(std::string output);
	void computeCalibration(std::string output);
};

inline
int call_set::getCalibrationBin(float gp0, float gp1) {
	float maxv = 0.0f;
	if (gp0 > maxv) { maxv = gp0; }
	if (gp1 > maxv) { maxv = gp1; }
	return (int)trunc(maxv * (N_BIN_CAL-1));
}

inline
int call_set::getCalibrationBin(float gp0, float gp1, float gp2) {
	//should not be affected by ploidy: gp1=gp2 with ploidy=1
	/*
	float maxv = 0.0f;
	if (gp0 > maxv) maxv = gp0;
	if (gp1 > maxv) maxv = gp1;
	if (gp2 > maxv) maxv = gp2;
	*/
	return (int)trunc(fmaxf(gp0, fmaxf(gp1, gp2)) * (N_BIN_CAL-1));
}

inline
int call_set::getTruth(float pl0, float pl1) {
	if (pl0 < 0.0f || pl1 < 0.0f) return -3;
	float sc = 1.0f / (pl0 + pl1);
	float p0 = pl0 * sc;
	float p1 = pl1 * sc;
	// Not certain enough about truth:
	if (p0 < T && p1 < T) return -4;
	// Certain enough about it:
	if (p0 > p1) return 0;
	if (p1 > p0) return 1;
	return -1;
}

inline
int call_set::getTruth(float pl0, float pl1, float pl2, int dp, int ploidy) {
	if (dp == bcf_int32_missing) return -5;
	if (dp < D) return -2;

	if (ploidy==1) return getTruth(pl0, pl1);


	if (pl0 < 0.0f || pl1 < 0.0f || pl2 < 0.0f) return -3;
	float sc = 1.0 / (pl0 + pl1 + pl2);
	float p0 = pl0 * sc;
	float p1 = pl1 * sc;
	float p2 = pl2 * sc;
	// Not certain enough about truth:
	if (p0 < T && p1 < T && p2 < T) return -4;
	// Certain enough about it:
	if (p0 > p1 && p0 > p2) return 0;
	if (p1 > p0 && p1 > p2) return 1;
	if (p2 > p0 && p2 > p1) return 2;
	return -1;
}

inline
int call_set::getFrequencyBin(float prob) {
	for (int b = 1 ; b < bins.size() ; b ++) {
		if (prob > bins[b-1] && prob <= bins[b]) return b-1;
	}
	return -1;
}

inline
int call_set::getMostLikely(const float gp0, const float gp1, const float ds, const float gp_filter) {
	if (ds < 0.0f && ds > 1.0f + 1e-7) return -2;
	if (gp0 < 0.0f || gp0 > 1.0f + 1e-7) return -3;
	if (gp1 < 0.0f || gp1 > 1.0f + 1e-7) return -3;

	if (std::max(gp0, gp1) < gp_filter) return -4;

	if (gp0 > gp1) return 0;
	if (gp1 > gp0) return 1;

	return -1;
}

inline
int call_set::getMostLikely(const float gp0, const float gp1, const float gp2, float ds, const int ploidy, const float gp_filter) {
	if (ploidy == 1) return getMostLikely(gp0,gp1,ds, gp_filter);

	if (ds < 0.0f && ds > 2.0f + 1e-7) return -2;
	if (gp0 < 0.0f || gp0 > 1.0f + 1e-7) return -3;
	if (gp1 < 0.0f || gp1 > 1.0f + 1e-7) return -3;
	if (gp2 < 0.0f || gp2 > 1.0f + 1e-7) return -3;

	if (std::max({gp0, gp1, gp2}) < gp_filter) return -4;

	if (gp0 > gp1 && gp0 > gp2) return 0;
	if (gp1 > gp0 && gp1 > gp2) return 1;
	if (gp2 > gp0 && gp2 > gp1) return 2;

	//TODO return -1. THIS IS JUST FOR TEST
	return 0;
}

inline
int call_set::getMostLikelyGT(int gt, float ds, int ploidy)
{
	if (ds < 0.0f && ds > ((float) ploidy) + 1e-7) return -2;
	if (gt < 0) return -3;
	return gt;

}

#endif