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

#include <models/haplotype_hmm.h>

haplotype_hmm::haplotype_hmm(const conditioning_set * _C) {
	C = _C;
	Emissions = vector < float > (2*C->n_vars, 0.0);
}

haplotype_hmm::~haplotype_hmm() {
	Alpha.clear();
	AlphaSum.clear();
	Emissions.clear();
}

void haplotype_hmm::resize() {
	AlphaSum.resize(C->n_sites);
	Alpha.resize(C->n_sites * C->n_states);
}

void haplotype_hmm::init(const vector < float > & HL) {
	double p0, p1;
	for (int l = 0 ; l < C->n_vars ; l ++) {
		p0 = HL[2*l+0] * C->ee + HL[2*l+1] * C->ed;
		p1 = HL[2*l+0] * C->ed + HL[2*l+1] * C->ee;
		Emissions[2*l+0] = p0 / (p0+p1);
		Emissions[2*l+1] = p1 / (p0+p1);
	}
}

void haplotype_hmm::computePosteriors(const vector < float > & HL, vector < float > & HP) {
	resize();
	init(HL);
	forward();
	backward(HL, HP);
}

void haplotype_hmm::forward() {
	double fact1, fact2;
	for (int l = 0 ; l < C->n_sites ; l ++) {
		AlphaSum[l] = 0.0;
		if (l == 0) {
			fact1 = 1.0 / C->n_states;
			for (int k = 0 ; k < C->n_states ; k ++) {
				Alpha[k] = Emissions[2*C->Vpoly[l]+_GET8(C->Hpoly[l][k/8], k%8)] * fact1;
				AlphaSum[l] += Alpha[k];
			}
		} else {
			fact1 = C->t[l-1] / C->n_states;
			fact2 = C->nt[l-1] / AlphaSum[l-1];
			for (int k = 0 ; k < C->n_states ; k ++) {
				Alpha[l*C->n_states+k] = (Alpha[(l-1)*C->n_states+k] * fact2 + fact1) * Emissions[2*C->Vpoly[l]+_GET8(C->Hpoly[l][k/8], k%8)];
				AlphaSum[l] += Alpha[l*C->n_states+k];
			}
		}
		//loglik += log (AlphaSum[l]);
	}
}

void haplotype_hmm::backward(const vector < float > & HL, vector < float > & HP) {
	double betaSumNext, betaSumCurr;
	double prob0 = 0.0, prob1 = 0.0;
	double hemit[2][2], pcopy[2], prob[2], betaSumTmp[2];
	vector < float > beta = vector < float > (C->n_states, 1.0);
	for (int l = C->n_sites-1 ; l >= 0 ; l --) {
		// Initilialization
		prob[0]=0.0;prob[1]=0.0;
		betaSumTmp[0]=0.0;betaSumTmp[1]=0.0;
		hemit[0][0] = C->ee * HL[2*C->Vpoly[l]+0] / Emissions[2*C->Vpoly[l] + 0];
		hemit[0][1] = C->ed * HL[2*C->Vpoly[l]+0] / Emissions[2*C->Vpoly[l] + 1];
		hemit[1][0] = C->ed * HL[2*C->Vpoly[l]+1] / Emissions[2*C->Vpoly[l] + 0];
		hemit[1][1] = C->ee * HL[2*C->Vpoly[l]+1] / Emissions[2*C->Vpoly[l] + 1];
		// Backward pass
		betaSumCurr = 0.0;
		if (l == C->n_sites - 1) {
			for (int k = 0 ; k < C->n_states ; k ++) {
				prob[_GET8(C->Hpoly[l][k/8], k%8)] += Alpha[l*C->n_states+k] * beta[k];
				betaSumTmp[_GET8(C->Hpoly[l][k/8], k%8)]++;
			}
		} else {
			pcopy[0] = C->nt[l] * Emissions[2*C->Vpoly[l+1] + 0] / betaSumNext;
			pcopy[1] = C->nt[l] * Emissions[2*C->Vpoly[l+1] + 1] / betaSumNext;

			for (int k = 0 ; k < C->n_states ; k ++) {
				beta[k] = beta[k] * pcopy[_GET8(C->Hpoly[l+1][k/8], k%8)] + C->t[l];
				prob[_GET8(C->Hpoly[l][k/8], k%8)] += Alpha[l*C->n_states+k] * beta[k];
				betaSumTmp[_GET8(C->Hpoly[l][k/8], k%8)]+=beta[k];
			}
		}
		// Expectation pass
		prob0 = prob[0]*hemit[0][0] + prob[1]*hemit[0][1];
		prob1 = prob[0]*hemit[1][0] + prob[1]*hemit[1][1];
		HP[2*C->Vpoly[l]+0] = prob0 / (prob0 + prob1);
		HP[2*C->Vpoly[l]+1] = prob1 / (prob0 + prob1);
		betaSumCurr = betaSumTmp[0]*Emissions[2*C->Vpoly[l] + 0] + betaSumTmp[1]*Emissions[2*C->Vpoly[l] + 1];
		betaSumNext = betaSumCurr / C->n_states;
		//loglik += log (betaSumPrev);
	}
	// Monomorphic sites
	 for (int l = 0 ; l < C->Vmono.size() ; l ++) {
		 prob0 = 0.0f; prob1 = 0.0f;
		 if (!C->Hmono[C->Vmono[l]]) {
			 prob0 = C->ee * HL[2*C->Vmono[l]+0];
	         prob1 = C->ed * HL[2*C->Vmono[l]+1];
		 } else {
			 prob0 = C->ed * HL[2*C->Vmono[l]+0];
	         prob1 = C->ee * HL[2*C->Vmono[l]+1];
		 }
		 HP[2*C->Vmono[l]+0] = prob0 / (prob0 + prob1);
	     HP[2*C->Vmono[l]+1] = prob1 / (prob0 + prob1);
	 }
}
