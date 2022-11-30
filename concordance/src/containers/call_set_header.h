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

#include <utils/otools.h>

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

	void write_record(htsFile *out_fd, bcf_hdr_t * out_hdr, bcf_hdr_t * hdr_in, bcf1_t *line);

	//
	void initialize(std::vector < double >,double, int); //bins
	void initialize(double, int); //allele-counts
	void initialize(std::vector < int >,double, int); //bins
	void initialize(std::string, double, int); //groups

	void setTargets(std::string fsamples);
	int getTruth(float, float);
	int getCalibrationBin(float , float );
	int getFrequencyBin(float);
	int getTruth(float, float, float, int, int);
	int getMostLikely(const float , const float , const float, const float);
	int getMostLikely(const float , const float , const float , const float, const int, const float);
	int getMostLikelyGT(int, float, int);
	int getCalibrationBin(float , float , float );

	//
	void readData(std::vector < std::string > &, std::vector < std::string > &, std::vector < std::string > &, std::vector < std::string > &, bpo::variables_map&, const float gp_filter, const std::string out_filename);
	void writeData(std::string);

	//
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
