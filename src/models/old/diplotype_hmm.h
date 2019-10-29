#ifndef _DIPLOTYPE_HMM_H
#define _DIPLOTYPE_HMM_H

#include <utils/otools.h>
#include <containers/bitmatrix.h>
#include <objects/hmm_parameters.h>

#define HAP_NUMBER 8
#define HAP_SCALE 50

#define ALLELE(hap , pos) ((hap) & (1<<(pos)))

class diplotype_hmm {
private:
	//EXTERNAL DATA
	bitmatrix & H;
	hmm_parameters & M;
	vector < char > HET;
	vector < bool > ALT;

	//COORDINATES & CONSTANTS
	unsigned int n_haps;
	unsigned int n_vars;

	//SEGMENTATION
	vector < int > segments;

	//CURSORS
	int curr_locus;
	int curr_segment_index;
	int curr_segment_locus;

	//DYNAMIC ARRAYS
	double probSumT1;
	double probSumT2;
	vector < double > prob1;
	vector < double > prob2;
	vector < double > probSumK1;
	vector < double > probSumK2;
	vector < double > probSumH1;
	vector < double > probSumH2;
	vector < vector < double > > Alpha;
	vector < vector < double > > Beta;
	vector < vector < double > > AlphaSum;
	vector < double > BetaSum;

	//STATIC ARRAYS
	double EMIT0[3][HAP_NUMBER];
	double EMIT1[3][HAP_NUMBER];
	double HPROBS [HAP_NUMBER * HAP_NUMBER];
	double sumHProbs, sumDProbs;

	//INLINED AND UNROLLED ROUTINES
	void INIT1();
	void INIT2();
	void SUM1();
	void SUM2();
	void SUMK1();
	void SUMK2();
	void COLLAPSE1(bool);
	void COLLAPSE2(bool);
	void RUN1(bool);
	void RUN2(bool);
	bool TRANSH();

public:
	//CONSTRUCTOR/DESTRUCTOR
	diplotype_hmm(bitmatrix &, hmm_parameters &, int, int);
	~diplotype_hmm();

	void reallocate(vector < bool > & H0, vector < bool > & H1);
	void forward();
	void backward();
	void sampleHaplotypes(vector < bool > & H0, vector < bool > & H1);
};

inline
void diplotype_hmm::INIT2() {
	if (HET[curr_locus] == -1) {
		bool ag = ALT[curr_locus];
		for(int k = 0, i = 0 ; k != n_haps ; ++k, i += HAP_NUMBER) {
			bool ah = H.get(curr_locus, k);
			if (ag != ah) fill(prob2.begin() + i, prob2.begin() + i + HAP_NUMBER, M.ed);
			else fill(prob2.begin() + i, prob2.begin() + i + HAP_NUMBER, M.ee);
		}
	} else {
		for(int k = 0, i = 0 ; k != n_haps ; ++k, i += HAP_NUMBER) {
			bool ah = H.get(curr_locus, k);
			if (ah) memcpy(&prob2[i], &(EMIT1[HET[curr_locus]][0]), HAP_NUMBER*sizeof(double));
			else memcpy(&prob2[i], &(EMIT0[HET[curr_locus]][0]), HAP_NUMBER*sizeof(double));
		}
	}
}

inline
void diplotype_hmm::INIT1() {
	if (HET[curr_locus] == -1) {
		bool ag = ALT[curr_locus];
		for(int k = 0, i = 0 ; k != n_haps ; ++k, i += HAP_NUMBER) {
			bool ah = H.get(curr_locus, k);
			if (ag != ah) fill(prob1.begin() + i, prob1.begin() + i + HAP_NUMBER, M.ed);
			else fill(prob1.begin() + i, prob1.begin() + i + HAP_NUMBER, M.ee);
		}
	} else {
		for(int k = 0, i = 0 ; k != n_haps ; ++k, i += HAP_NUMBER) {
			bool ah = H.get(curr_locus, k);
			if (ah) memcpy(&prob1[i], &(EMIT1[HET[curr_locus]][0]), HAP_NUMBER*sizeof(double));
			else memcpy(&prob1[i], &(EMIT0[HET[curr_locus]][0]), HAP_NUMBER*sizeof(double));
		}
	}
}

inline
void diplotype_hmm::SUM2() {
	double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0, sum5 = 0.0, sum6 = 0.0, sum7 = 0.0;
	for(int k = 0, i = 0 ; k != n_haps ; ++k, i += HAP_NUMBER) {
		sum0 += prob2[i + 0];
		sum1 += prob2[i + 1];
		sum2 += prob2[i + 2];
		sum3 += prob2[i + 3];
		sum4 += prob2[i + 4];
		sum5 += prob2[i + 5];
		sum6 += prob2[i + 6];
		sum7 += prob2[i + 7];
	}
	probSumH2[0] = sum0;
	probSumH2[1] = sum1;
	probSumH2[2] = sum2;
	probSumH2[3] = sum3;
	probSumH2[4] = sum4;
	probSumH2[5] = sum5;
	probSumH2[6] = sum6;
	probSumH2[7] = sum7;
	probSumT2 = sum0 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7;
}

inline
void diplotype_hmm::SUM1() {
	double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0, sum5 = 0.0, sum6 = 0.0, sum7 = 0.0;
	for(int k = 0, i = 0 ; k != n_haps ; ++k, i += HAP_NUMBER) {
		sum0 += prob1[i + 0];
		sum1 += prob1[i + 1];
		sum2 += prob1[i + 2];
		sum3 += prob1[i + 3];
		sum4 += prob1[i + 4];
		sum5 += prob1[i + 5];
		sum6 += prob1[i + 6];
		sum7 += prob1[i + 7];
	}
	probSumH1[0] = sum0;
	probSumH1[1] = sum1;
	probSumH1[2] = sum2;
	probSumH1[3] = sum3;
	probSumH1[4] = sum4;
	probSumH1[5] = sum5;
	probSumH1[6] = sum6;
	probSumH1[7] = sum7;
	probSumT1 = sum0 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7;
}

inline
void diplotype_hmm::SUMK2() {
	for(int k = 0, i = 0 ; k != n_haps ; ++k, i += HAP_NUMBER) {
		probSumK2[k] = prob2[i+0] + prob2[i+1] + prob2[i+2] + prob2[i+3] + prob2[i+4] + prob2[i+5] + prob2[i+6] + prob2[i+7];
	}
}

inline
void diplotype_hmm::SUMK1() {
	for(int k = 0, i = 0 ; k != n_haps ; ++k, i += HAP_NUMBER) {
		probSumK1[k] = prob1[i+0] + prob1[i+1] + prob1[i+2] + prob1[i+3] + prob1[i+4] + prob1[i+5] + prob1[i+6] + prob1[i+7];
	}
}

inline
void diplotype_hmm::COLLAPSE2(bool forward) {
	double tmp_prob0 = M.nt[curr_locus-forward] / probSumT1;
	double tmp_prob1 = M.t[curr_locus-forward] / n_haps;
	for(int k = 0, i = 0 ; k != n_haps ; ++k, i += HAP_NUMBER) {
		double factor = probSumK1[k] * tmp_prob0 + tmp_prob1;
		prob2[i + 0] *= factor;
		prob2[i + 1] *= factor;
		prob2[i + 2] *= factor;
		prob2[i + 3] *= factor;
		prob2[i + 4] *= factor;
		prob2[i + 5] *= factor;
		prob2[i + 6] *= factor;
		prob2[i + 7] *= factor;
	}
}

inline
void diplotype_hmm::COLLAPSE1(bool forward) {
	double tmp_prob0 = M.nt[curr_locus-forward] / probSumT2;
	double tmp_prob1 = M.t[curr_locus-forward] / n_haps;
	for(int k = 0, i = 0 ; k != n_haps ; ++k, i += HAP_NUMBER) {
		double factor = probSumK2[k] * tmp_prob0 + tmp_prob1;
		prob1[i + 0] *= factor;
		prob1[i + 1] *= factor;
		prob1[i + 2] *= factor;
		prob1[i + 3] *= factor;
		prob1[i + 4] *= factor;
		prob1[i + 5] *= factor;
		prob1[i + 6] *= factor;
		prob1[i + 7] *= factor;
	}
}

inline
void diplotype_hmm::RUN2(bool forward) {
	double nt = M.nt[curr_locus-forward] / probSumT1;
	double tfreq = M.t[curr_locus-forward] / (n_haps * probSumT1);
	double tFreq0 = probSumH1[0] * tfreq;
	double tFreq1 = probSumH1[1] * tfreq;
	double tFreq2 = probSumH1[2] * tfreq;
	double tFreq3 = probSumH1[3] * tfreq;
	double tFreq4 = probSumH1[4] * tfreq;
	double tFreq5 = probSumH1[5] * tfreq;
	double tFreq6 = probSumH1[6] * tfreq;
	double tFreq7 = probSumH1[7] * tfreq;
	for(int k = 0, i = 0 ; k != n_haps ; ++k, i += HAP_NUMBER) {
		prob2[i + 0] *= prob1[i + 0] * nt + tFreq0;
		prob2[i + 1] *= prob1[i + 1] * nt + tFreq1;
		prob2[i + 2] *= prob1[i + 2] * nt + tFreq2;
		prob2[i + 3] *= prob1[i + 3] * nt + tFreq3;
		prob2[i + 4] *= prob1[i + 4] * nt + tFreq4;
		prob2[i + 5] *= prob1[i + 5] * nt + tFreq5;
		prob2[i + 6] *= prob1[i + 6] * nt + tFreq6;
		prob2[i + 7] *= prob1[i + 7] * nt + tFreq7;
	}
}


inline
void diplotype_hmm::RUN1(bool forward) {
	double nt = M.nt[curr_locus-forward] / probSumT2;
	double tfreq = M.t[curr_locus-forward] / (n_haps * probSumT2);
	double tFreq0 = probSumH2[0] * tfreq;
	double tFreq1 = probSumH2[1] * tfreq;
	double tFreq2 = probSumH2[2] * tfreq;
	double tFreq3 = probSumH2[3] * tfreq;
	double tFreq4 = probSumH2[4] * tfreq;
	double tFreq5 = probSumH2[5] * tfreq;
	double tFreq6 = probSumH2[6] * tfreq;
	double tFreq7 = probSumH2[7] * tfreq;
	for(int k = 0, i = 0 ; k != n_haps ; ++k, i += HAP_NUMBER) {
		prob1[i + 0] *= prob2[i + 0] * nt + tFreq0;
		prob1[i + 1] *= prob2[i + 1] * nt + tFreq1;
		prob1[i + 2] *= prob2[i + 2] * nt + tFreq2;
		prob1[i + 3] *= prob2[i + 3] * nt + tFreq3;
		prob1[i + 4] *= prob2[i + 4] * nt + tFreq4;
		prob1[i + 5] *= prob2[i + 5] * nt + tFreq5;
		prob1[i + 6] *= prob2[i + 6] * nt + tFreq6;
		prob1[i + 7] *= prob2[i + 7] * nt + tFreq7;
	}
}

inline
bool diplotype_hmm::TRANSH() {
	sumHProbs = 0.0;
	for (int h1 = 0 ; h1 < HAP_NUMBER ; h1++) {
		double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0, sum5 = 0.0, sum6 = 0.0, sum7 = 0.0;
		for (int k = 0 ; k < n_haps ; k ++) {
			double alpha = Alpha[curr_segment_index-1][k*HAP_NUMBER + h1] * M.nt[curr_locus-1] + AlphaSum[curr_segment_index - 1][h1] * M.t[curr_locus-1] / n_haps;
			sum0 += alpha * Beta[curr_segment_index][k*HAP_NUMBER + 0];
			sum1 += alpha * Beta[curr_segment_index][k*HAP_NUMBER + 1];
			sum2 += alpha * Beta[curr_segment_index][k*HAP_NUMBER + 2];
			sum3 += alpha * Beta[curr_segment_index][k*HAP_NUMBER + 3];
			sum4 += alpha * Beta[curr_segment_index][k*HAP_NUMBER + 4];
			sum5 += alpha * Beta[curr_segment_index][k*HAP_NUMBER + 5];
			sum6 += alpha * Beta[curr_segment_index][k*HAP_NUMBER + 6];
			sum7 += alpha * Beta[curr_segment_index][k*HAP_NUMBER + 7];
		}
		HPROBS[h1*HAP_NUMBER+0] = sum0;
		HPROBS[h1*HAP_NUMBER+1] = sum1;
		HPROBS[h1*HAP_NUMBER+2] = sum2;
		HPROBS[h1*HAP_NUMBER+3] = sum3;
		HPROBS[h1*HAP_NUMBER+4] = sum4;
		HPROBS[h1*HAP_NUMBER+5] = sum5;
		HPROBS[h1*HAP_NUMBER+6] = sum6;
		HPROBS[h1*HAP_NUMBER+7] = sum7;
		sumHProbs += sum0 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7;
	}
	return (isnan(sumHProbs) || sumHProbs < numeric_limits<double>::min());
}

#endif
