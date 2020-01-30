#ifndef _CALLSET_H
#define _CALLSET_H

#include <utils/otools.h>

#define N_BIN_CAL 100

const float unphred[256] = { 1.000000e+00, 7.943282e-01, 6.309573e-01, 5.011872e-01, 3.981072e-01, 3.162278e-01, 2.511886e-01, 1.995262e-01, 1.584893e-01, 1.258925e-01, 1.000000e-01, 7.943282e-02, 6.309573e-02, 5.011872e-02, 3.981072e-02, 3.162278e-02, 2.511886e-02, 1.995262e-02, 1.584893e-02, 1.258925e-02, 1.000000e-02, 7.943282e-03, 6.309573e-03, 5.011872e-03, 3.981072e-03, 3.162278e-03, 2.511886e-03, 1.995262e-03, 1.584893e-03, 1.258925e-03, 1.000000e-03, 7.943282e-04, 6.309573e-04, 5.011872e-04, 3.981072e-04, 3.162278e-04, 2.511886e-04, 1.995262e-04, 1.584893e-04, 1.258925e-04, 1.000000e-04, 7.943282e-05, 6.309573e-05, 5.011872e-05, 3.981072e-05, 3.162278e-05, 2.511886e-05, 1.995262e-05, 1.584893e-05, 1.258925e-05, 1.000000e-05, 7.943282e-06, 6.309573e-06, 5.011872e-06, 3.981072e-06, 3.162278e-06, 2.511886e-06, 1.995262e-06, 1.584893e-06, 1.258925e-06, 1.000000e-06, 7.943282e-07, 6.309573e-07, 5.011872e-07, 3.981072e-07, 3.162278e-07, 2.511886e-07, 1.995262e-07, 1.584893e-07, 1.258925e-07, 1.000000e-07, 7.943282e-08, 6.309573e-08, 5.011872e-08, 3.981072e-08, 3.162278e-08, 2.511886e-08, 1.995262e-08, 1.584893e-08, 1.258925e-08, 1.000000e-08, 7.943282e-09, 6.309573e-09, 5.011872e-09, 3.981072e-09, 3.162278e-09, 2.511886e-09, 1.995262e-09, 1.584893e-09, 1.258925e-09, 1.000000e-09, 7.943282e-10, 6.309573e-10, 5.011872e-10, 3.981072e-10, 3.162278e-10, 2.511886e-10, 1.995262e-10, 1.584893e-10, 1.258925e-10, 1.000000e-10, 7.943282e-11, 6.309573e-11, 5.011872e-11, 3.981072e-11, 3.162278e-11, 2.511886e-11, 1.995262e-11, 1.584893e-11, 1.258925e-11, 1.000000e-11, 7.943282e-12, 6.309573e-12, 5.011872e-12, 3.981072e-12, 3.162278e-12, 2.511886e-12, 1.995262e-12, 1.584893e-12, 1.258925e-12, 1.000000e-12, 7.943282e-13, 6.309573e-13, 5.011872e-13, 3.981072e-13, 3.162278e-13, 2.511886e-13, 1.995262e-13, 1.584893e-13, 1.258925e-13, 1.000000e-13, 7.943282e-14, 6.309573e-14, 5.011872e-14, 3.981072e-14, 3.162278e-14, 2.511886e-14, 1.995262e-14, 1.584893e-14, 1.258925e-14, 1.000000e-14, 7.943282e-15, 6.309573e-15, 5.011872e-15, 3.981072e-15, 3.162278e-15, 2.511886e-15, 1.995262e-15, 1.584893e-15, 1.258925e-15, 1.000000e-15, 7.943282e-16, 6.309573e-16, 5.011872e-16, 3.981072e-16, 3.162278e-16, 2.511886e-16, 1.995262e-16, 1.584893e-16, 1.258925e-16, 1.000000e-16, 7.943282e-17, 6.309573e-17, 5.011872e-17, 3.981072e-17, 3.162278e-17, 2.511886e-17, 1.995262e-17, 1.584893e-17, 1.258925e-17, 1.000000e-17, 7.943282e-18, 6.309573e-18, 5.011872e-18, 3.981072e-18, 3.162278e-18, 2.511886e-18, 1.995262e-18, 1.584893e-18, 1.258925e-18, 1.000000e-18, 7.943282e-19, 6.309573e-19, 5.011872e-19, 3.981072e-19, 3.162278e-19, 2.511886e-19, 1.995262e-19, 1.584893e-19, 1.258925e-19, 1.000000e-19, 7.943282e-20, 6.309573e-20, 5.011872e-20, 3.981072e-20, 3.162278e-20, 2.511886e-20, 1.995262e-20, 1.584893e-20, 1.258925e-20, 1.000000e-20, 7.943282e-21, 6.309573e-21, 5.011872e-21, 3.981072e-21, 3.162278e-21, 2.511886e-21, 1.995262e-21, 1.584893e-21, 1.258925e-21, 1.000000e-21, 7.943282e-22, 6.309573e-22, 5.011872e-22, 3.981072e-22, 3.162278e-22, 2.511886e-22, 1.995262e-22, 1.584893e-22, 1.258925e-22, 1.000000e-22, 7.943282e-23, 6.309573e-23, 5.011872e-23, 3.981072e-23, 3.162278e-23, 2.511886e-23, 1.995262e-23, 1.584893e-23, 1.258925e-23, 1.000000e-23, 7.943282e-24, 6.309573e-24, 5.011872e-24, 3.981072e-24, 3.162278e-24, 2.511886e-24, 1.995262e-24, 1.584893e-24, 1.258925e-24, 1.000000e-24, 7.943282e-25, 6.309573e-25, 5.011872e-25, 3.981072e-25, 3.162278e-25, 2.511886e-25, 1.995262e-25, 1.584893e-25, 1.258925e-25, 1.000000e-25, 7.943282e-26, 6.309573e-26, 5.011872e-26, 3.981072e-26, 3.162278e-26};

class call_set {
public:
	// Sample IDs [N] & Bins [L]
	int N, L;
	double T;
	vector < string > samples;
	vector < double > bins;

	// Per sample concordance [3xN]
	vector < unsigned long int > genotype_spl_errors;
	vector < unsigned long int > genotype_spl_totals;

	// Per bin concordance [3xL]
	vector < unsigned long int > genotype_bin_errors;
	vector < unsigned long int > genotype_bin_totals;

	// Concordance for calibration [100]
	vector < unsigned long int > genotype_cal_errors;
	vector < unsigned long int > genotype_cal_totals;

	// R2 per bin
	vector < stats2D > rsquared_bin;
	vector < stats1D > frequency_bin;

	// R2 per sample
	vector < stats2D > rsquared_spl;

	//
	call_set ();
	~call_set();

	//
	void initialize(vector < double >, double);
	int getTruth(double, double, double);
	int getFrequencyBin(float);
	int getMostLikely(float , float , float );
	int getCalibrationBin(float , float , float );

	//
	void readData(vector < string > &, vector < string > &, vector < string > &, vector < string > &);
	void writeData(string);

	//
	void computeRsquaredPerBin(string output);
	void computeRsquaredPerBinPerSample(string output);

	void computeConcordancePerBIN(string output);
	void concordanceOverall(string output);
	void concordancePerIndividual(string output);
	void computeCalibration(string output);
};


inline
int call_set::getMostLikely(float gp0, float gp1, float gp2) {
	if (gp0 > gp1 && gp0 > gp2) return 0;
	if (gp1 > gp0 && gp1 > gp2) return 1;
	if (gp2 > gp0 && gp2 > gp1) return 2;
	return 0;
}

inline
int call_set::getCalibrationBin(float gp0, float gp1, float gp2) {
	float maxv = 0.0f;
	if (gp0 > maxv) { maxv = gp0; }
	if (gp1 > maxv) { maxv = gp1; }
	if (gp2 > maxv) { maxv = gp2; }
	return (int)trunc(maxv * (N_BIN_CAL-1));
}

inline
int call_set::getTruth(double pl0, double pl1, double pl2) {
	if (pl0 < 0.0f || pl1 < 0.0f || pl2 < 0.0f) return -1;
	double sc = 1.0 / (pl0 + pl1 + pl2);
	double p0 = pl0 * sc;
	double p1 = pl1 * sc;
	double p2 = pl2 * sc;
	// Not certain enough about truth:
	if (p0 < T && p1 < T && p2 < T) return -1;
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

#endif
