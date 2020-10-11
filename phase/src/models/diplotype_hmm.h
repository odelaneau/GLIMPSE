/*******************************************************************************
 * Copyright (C) 2020 Olivier Delaneau, University of Lausanne
 * Copyright (C) 2020 Simone Rubinacci, University of Lausanne
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

#ifndef _DIPLOTYPE_HMM_H
#define _DIPLOTYPE_HMM_H

#include <utils/otools.h>
#include <containers/conditioning_set.h>

#define HAP_NUMBER 8
#define HAP_SCALE 50

#define ALLELE(hap , pos) ((hap) & (1<<(pos)))

class diplotype_hmm {
private:

	const conditioning_set * C;

	//EXTERNAL DATA
	vector < char > HET;
	vector < bool > ALT;

	//COORDINATES & CONSTANTS
	unsigned int n_segs;

	//SEGMENTATION
	vector < int > segments;

	//CURSORS
	int curr_locus;
	int curr_segment_index;
	int curr_segment_locus;

	//DYNAMIC ARRAYS
	float probSumT1;
	float probSumT2;
	vector < float > prob1;
	vector < float > prob2;
	vector < float > probSumK1;
	vector < float > probSumK2;
	vector < float > probSumH1;
	vector < float > probSumH2;
	vector < float > Alpha;
	vector < float > Beta;
	vector < float > AlphaSum;
	vector < float > BetaSum;

	//STATIC ARRAYS
	float EMIT0[3][HAP_NUMBER];
	float EMIT1[3][HAP_NUMBER];
	float HPROBS [HAP_NUMBER * HAP_NUMBER];
	float sumHProbs, sumDProbs;

	//INLINED AND UNROLLED ROUTINES
	void INIT(bool);
	void SUM(bool);
	void SUMK(bool);
	void COLLAPSE(bool, bool);
	void RUN(bool, bool);
	bool TRANSH();

public:
	//CONSTRUCTOR/DESTRUCTOR
	diplotype_hmm(const conditioning_set * C);

	~diplotype_hmm();

	void reallocate(const vector < bool > &, const vector < bool > &);
	void forward();
	void backward();
	void rephaseHaplotypes(vector < bool > &, vector < bool > &);
};

inline
void diplotype_hmm::INIT(bool paired) {
	if (paired) {
		if (HET[curr_locus] == -1) {
			bool ag = ALT[curr_locus];
			for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER) {
				bool ah = _GET8(C->Hpoly[curr_locus][k/8], k%8);
				if (ag != ah) fill(prob2.begin() + i, prob2.begin() + i + HAP_NUMBER, C->ed);
				else fill(prob2.begin() + i, prob2.begin() + i + HAP_NUMBER, C->ee);
			}
		} else {
			for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER) {
				//bool ah = C->Hpoly[curr_locus*C->n_states+k];
				bool ah = _GET8(C->Hpoly[curr_locus][k/8], k%8);
				if (ah) memcpy(&prob2[i], &(EMIT1[HET[curr_locus]][0]), HAP_NUMBER*sizeof(float));
				else memcpy(&prob2[i], &(EMIT0[HET[curr_locus]][0]), HAP_NUMBER*sizeof(float));
			}
		}
	} else {
		if (HET[curr_locus] == -1) {
			bool ag = ALT[curr_locus];
			for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER) {
				//bool ah = C->Hpoly[curr_locus*C->n_states+k];
				bool ah = _GET8(C->Hpoly[curr_locus][k/8], k%8);
				if (ag != ah) fill(prob1.begin() + i, prob1.begin() + i + HAP_NUMBER, C->ed);
				else fill(prob1.begin() + i, prob1.begin() + i + HAP_NUMBER, C->ee);
			}
		} else {
			for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER) {
				//bool ah = C->Hpoly[curr_locus*C->n_states+k];
				bool ah = _GET8(C->Hpoly[curr_locus][k/8], k%8);
				if (ah) memcpy(&prob1[i], &(EMIT1[HET[curr_locus]][0]), HAP_NUMBER*sizeof(float));
				else memcpy(&prob1[i], &(EMIT0[HET[curr_locus]][0]), HAP_NUMBER*sizeof(float));
			}
		}
	}
}

inline
void diplotype_hmm::SUM(bool paired) {
	float sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0, sum5 = 0.0, sum6 = 0.0, sum7 = 0.0;
	if (paired) {
		for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER) {
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
	} else {
		for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER) {
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
}

inline
void diplotype_hmm::SUMK(bool paired) {
	if (paired) {
		for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER) {
			probSumK2[k] = prob2[i+0] + prob2[i+1] + prob2[i+2] + prob2[i+3] + prob2[i+4] + prob2[i+5] + prob2[i+6] + prob2[i+7];
		}
	} else {
		for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER) {
			probSumK1[k] = prob1[i+0] + prob1[i+1] + prob1[i+2] + prob1[i+3] + prob1[i+4] + prob1[i+5] + prob1[i+6] + prob1[i+7];
		}
	}
}

inline
void diplotype_hmm::COLLAPSE(bool forward, bool paired) {
	if (paired) {
		float tmp_prob0 = C->nt[curr_locus-forward] / probSumT1;
		float tmp_prob1 = C->t[curr_locus-forward] / C->n_states;
		for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER) {
			float factor = probSumK1[k] * tmp_prob0 + tmp_prob1;
			prob2[i + 0] *= factor;
			prob2[i + 1] *= factor;
			prob2[i + 2] *= factor;
			prob2[i + 3] *= factor;
			prob2[i + 4] *= factor;
			prob2[i + 5] *= factor;
			prob2[i + 6] *= factor;
			prob2[i + 7] *= factor;
		}
	} else {
		float tmp_prob0 = C->nt[curr_locus-forward] / probSumT2;
		float tmp_prob1 = C->t[curr_locus-forward] / C->n_states;
		for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER) {
			float factor = probSumK2[k] * tmp_prob0 + tmp_prob1;
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
}

inline
void diplotype_hmm::RUN(bool forward, bool paired) {
	if (paired) {
		float nont = C->nt[curr_locus-forward] / probSumT1;
		float tfreq = C->t[curr_locus-forward] / (C->n_states * probSumT1);
		float tFreq0 = probSumH1[0] * tfreq;
		float tFreq1 = probSumH1[1] * tfreq;
		float tFreq2 = probSumH1[2] * tfreq;
		float tFreq3 = probSumH1[3] * tfreq;
		float tFreq4 = probSumH1[4] * tfreq;
		float tFreq5 = probSumH1[5] * tfreq;
		float tFreq6 = probSumH1[6] * tfreq;
		float tFreq7 = probSumH1[7] * tfreq;
		for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER) {
			prob2[i + 0] *= prob1[i + 0] * nont + tFreq0;
			prob2[i + 1] *= prob1[i + 1] * nont + tFreq1;
			prob2[i + 2] *= prob1[i + 2] * nont + tFreq2;
			prob2[i + 3] *= prob1[i + 3] * nont + tFreq3;
			prob2[i + 4] *= prob1[i + 4] * nont + tFreq4;
			prob2[i + 5] *= prob1[i + 5] * nont + tFreq5;
			prob2[i + 6] *= prob1[i + 6] * nont + tFreq6;
			prob2[i + 7] *= prob1[i + 7] * nont + tFreq7;
		}
	} else {
		float nont = C->nt[curr_locus-forward] / probSumT2;
		float tfreq = C->t[curr_locus-forward] / (C->n_states * probSumT2);
		float tFreq0 = probSumH2[0] * tfreq;
		float tFreq1 = probSumH2[1] * tfreq;
		float tFreq2 = probSumH2[2] * tfreq;
		float tFreq3 = probSumH2[3] * tfreq;
		float tFreq4 = probSumH2[4] * tfreq;
		float tFreq5 = probSumH2[5] * tfreq;
		float tFreq6 = probSumH2[6] * tfreq;
		float tFreq7 = probSumH2[7] * tfreq;
		for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER) {
			prob1[i + 0] *= prob2[i + 0] * nont + tFreq0;
			prob1[i + 1] *= prob2[i + 1] * nont + tFreq1;
			prob1[i + 2] *= prob2[i + 2] * nont + tFreq2;
			prob1[i + 3] *= prob2[i + 3] * nont + tFreq3;
			prob1[i + 4] *= prob2[i + 4] * nont + tFreq4;
			prob1[i + 5] *= prob2[i + 5] * nont + tFreq5;
			prob1[i + 6] *= prob2[i + 6] * nont + tFreq6;
			prob1[i + 7] *= prob2[i + 7] * nont + tFreq7;
		}
	}
}

inline
bool diplotype_hmm::TRANSH() {
	sumHProbs = 0.0;
	for (int h1 = 0 ; h1 < HAP_NUMBER ; h1++) {
		float sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0, sum5 = 0.0, sum6 = 0.0, sum7 = 0.0;
		unsigned long prev_offset = (curr_segment_index-1) * C->n_states * HAP_NUMBER;
		unsigned long curr_offset = curr_segment_index * C->n_states * HAP_NUMBER;
		for (int k = 0 ; k < C->n_states ; k ++) {
			float alpha = Alpha[prev_offset + k*HAP_NUMBER + h1] * C->nt[curr_locus-1] + AlphaSum[(curr_segment_index-1)*HAP_NUMBER + h1] * C->t[curr_locus-1] / C->n_states;
			sum0 += alpha * Beta[curr_offset + k*HAP_NUMBER + 0];
			sum1 += alpha * Beta[curr_offset + k*HAP_NUMBER + 1];
			sum2 += alpha * Beta[curr_offset + k*HAP_NUMBER + 2];
			sum3 += alpha * Beta[curr_offset + k*HAP_NUMBER + 3];
			sum4 += alpha * Beta[curr_offset + k*HAP_NUMBER + 4];
			sum5 += alpha * Beta[curr_offset + k*HAP_NUMBER + 5];
			sum6 += alpha * Beta[curr_offset + k*HAP_NUMBER + 6];
			sum7 += alpha * Beta[curr_offset + k*HAP_NUMBER + 7];
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
	return (isnan(sumHProbs) || sumHProbs < numeric_limits<float>::min());
}

#endif
