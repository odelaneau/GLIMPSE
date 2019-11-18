////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2018 Olivier Delaneau, University of Lausanne
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////
#include <models/haplotype_hmm.h>

haplotype_hmm::haplotype_hmm(conditioning_set * _C) {
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

void haplotype_hmm::init(vector < float > & HL) {
	double p0, p1;
	for (int l = 0 ; l < C->n_vars ; l ++) {
		p0 = HL[2*l+0] * C->ee + HL[2*l+1] * C->ed;
		p1 = HL[2*l+0] * C->ed + HL[2*l+1] * C->ee;

		Emissions[2*l+0] = p0 / (p0+p1);
		Emissions[2*l+1] = p1 / (p0+p1);
		//Emissions[2*l+0] = 1.0;
		//Emissions[2*l+1] = (p1/(p0+p1))*((p0+p1)/p0);

	}
}

void haplotype_hmm::computePosteriors(vector < float > & HL, vector < float > & HP) {
	resize();
	//vrb.bullet("HAP Resize (" + stb.str(tac.rel_time()*1.0, 1) + "ms)");tac.clock();
	init(HL);
	//vrb.bullet("HAP Init (" + stb.str(tac.rel_time()*1.0, 1) + "ms)");tac.clock();
	forward();
	//vrb.bullet("HAP Forward (" + stb.str(tac.rel_time()*1.0, 1) + "ms)");tac.clock();
	backward(HL, HP);
	//vrb.bullet("HAP Backward (" + stb.str(tac.rel_time()*1.0, 1) + "ms)");tac.clock();
}

void haplotype_hmm::forward() {
	//double loglik = 0.0, fact1, fact2;
	double fact1, fact2;

	for (int l = 0 ; l < C->n_sites ; l ++) {
		AlphaSum[l] = 0.0;
		if (l == 0) {
			fact1 = 1.0 / C->n_states;
			for (int k = 0 ; k < C->n_states ; k ++) {
				Alpha[k] = Emissions[2*C->Vpoly[l]+C->Hpoly[l*C->n_states+k]] * fact1;
				AlphaSum[l] += Alpha[k];
			}
		} else {
			fact1 = C->t[l-1] / C->n_states;
			fact2 = C->nt[l-1] / AlphaSum[l-1];

			/*
			 * To accelerate it:
			 * 		For closeby variants , assume zero recombination rates, so that fact2=1.0 and fact1=0
			 * 		Then, rescale likelihoods so that P(G=0) = 1.0
			 * 		So you just have to iterate over non-zero alleles and multiply alpha values by P(G=1)
			 * 		This should greatly speed up things as zeros values are everywhere, especially when the reference panel is large!!!!
			 */

			for (int k = 0 ; k < C->n_states ; k ++) {
				Alpha[l*C->n_states+k] = (Alpha[(l-1)*C->n_states+k] * fact2 + fact1) * Emissions[2*C->Vpoly[l]+C->Hpoly[l*C->n_states+k]];
				AlphaSum[l] += Alpha[l*C->n_states+k];
			}
		}
		//loglik += log (AlphaSum[l]);
	}
	//cout << "F = " << loglik << endl;
}

void haplotype_hmm::backward(vector < float > & HL, vector < float > & HP) {
	//double loglik = 0.0, fact1, betaSumPrev, betaSumCurr, prob0, prob1;
	double betaSumNext, betaSumCurr;
	double prob0 = 0.0, prob1 = 0.0;

	double hemit[2][2], pcopy[2], prob[2], betaSumTmp[2];
	vector < float > beta = vector < float > (C->n_states, 1.0);
	for (int l = C->n_sites-1 ; l >= 0 ; l --) {
		//Set up values
		//prob0 = 0.0; prob1 = 0.0;
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
				prob[C->Hpoly[l*C->n_states+k]] += Alpha[l*C->n_states+k] * beta[k];
				betaSumTmp[C->Hpoly[l*C->n_states+k]]++;
			}
		} else {
			pcopy[0] = C->nt[l] * Emissions[2*C->Vpoly[l+1] + 0] / betaSumNext;
			pcopy[1] = C->nt[l] * Emissions[2*C->Vpoly[l+1] + 1] / betaSumNext;

			for (int k = 0 ; k < C->n_states ; k ++) {
				beta[k] = beta[k] * pcopy[C->Hpoly[(l+1)*C->n_states+k]] + C->t[l];
				prob[C->Hpoly[l*C->n_states+k]] += Alpha[l*C->n_states+k] * beta[k];
				betaSumTmp[C->Hpoly[l*C->n_states+k]]+=beta[k];
			}
		}

		prob0 = prob[0]*hemit[0][0] + prob[1]*hemit[0][1];
		prob1 = prob[0]*hemit[1][0] + prob[1]*hemit[1][1];
		HP[2*C->Vpoly[l]+0] = prob0 / (prob0 + prob1);
		HP[2*C->Vpoly[l]+1] = prob1 / (prob0 + prob1);
		betaSumCurr = betaSumTmp[0]*Emissions[2*C->Vpoly[l] + 0] + betaSumTmp[1]*Emissions[2*C->Vpoly[l] + 1];
		betaSumNext = betaSumCurr / C->n_states;
		//loglik += log (betaSumPrev);
	}
	//cout << "B = " << loglik << endl;
	for (int l = 0 ; l < C->Vmono.size() ; l ++) {
		prob0 = 0.0; prob1 = 0.0;
		if (!C->Hmono[C->Vmono[l]]) {
			prob0 = C->ee * HL[2*C->Vmono[l]+0];
			prob1 = C->ed * HL[2*C->Vmono[l]+1];
		} else {
			prob0 = C->ed * HL[2*C->Vmono[l]+0];
			prob1 = C->ee * HL[2*C->Vmono[l]+1];
		}
		//cout << prob0 / (prob0+prob1) << " " << prob1 / (prob0+prob1) << endl;
		HP[2*C->Vmono[l]+0] = prob0 / (prob0 + prob1);
		HP[2*C->Vmono[l]+1] = prob1 / (prob0 + prob1);
	}
}



/*
void haplotype_hmm::forward() {
	//double loglik = 0.0, fact1, fact2;
	double fact1, fact2;

	for (int l = 0 ; l < C->n_sites ; l ++) {
		AlphaSum[l] = 0.0;
		if (l == 0) {
			fact1 = 1.0 / C->n_states;
			for (int k = 0 ; k < C->n_states ; k ++) {
				Alpha[k] = Emissions[2*C->Vpoly[l]+C->Hpoly[l*C->n_states+k]] * fact1;
				AlphaSum[l] += Alpha[k];
			}
		} else {
			fact1 = C->t[l-1] / C->n_states;
			fact2 = C->nt[l-1] / AlphaSum[l-1];

			int i = 0;
			int repeat = (C->n_states/4);
			int left = (C->n_states%4);
			double sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;

			while (repeat --) {
				Alpha[l*C->n_states+i+0] = (Alpha[(l-1)*C->n_states+i+0] * fact2 + fact1) * Emissions[2*C->Vpoly[l]+C->Hpoly[l*C->n_states+i+0]];
				Alpha[l*C->n_states+i+1] = (Alpha[(l-1)*C->n_states+i+1] * fact2 + fact1) * Emissions[2*C->Vpoly[l]+C->Hpoly[l*C->n_states+i+1]];
				Alpha[l*C->n_states+i+2] = (Alpha[(l-1)*C->n_states+i+2] * fact2 + fact1) * Emissions[2*C->Vpoly[l]+C->Hpoly[l*C->n_states+i+2]];
				Alpha[l*C->n_states+i+3] = (Alpha[(l-1)*C->n_states+i+3] * fact2 + fact1) * Emissions[2*C->Vpoly[l]+C->Hpoly[l*C->n_states+i+3]];
				sum0 += Alpha[l*C->n_states+i+0];
				sum1 += Alpha[l*C->n_states+i+1];
				sum2 += Alpha[l*C->n_states+i+2];
				sum3 += Alpha[l*C->n_states+i+3];
				i += 4;
			}


			switch (left) {
			case 3: sum2 +=(Alpha[l*C->n_states+i+2] = (Alpha[l*C->n_states+i+2-C->n_states] * fact2 + fact1) * Emissions[2*C->Vpoly[l]+C->Hpoly[l*C->n_states+i+2]]);
			case 2: sum1 +=(Alpha[l*C->n_states+i+1] = (Alpha[l*C->n_states+i+1-C->n_states] * fact2 + fact1) * Emissions[2*C->Vpoly[l]+C->Hpoly[l*C->n_states+i+1]]);
			case 1: sum0 +=(Alpha[l*C->n_states+i+0] = (Alpha[l*C->n_states+i+0-C->n_states] * fact2 + fact1) * Emissions[2*C->Vpoly[l]+C->Hpoly[l*C->n_states+i+0]]);
			}

			AlphaSum[l] = sum0 + sum1 + sum2 + sum3;

			//for (int k = 0 ; k < C->n_states ; k ++) {
			//	Alpha[l*C->n_states+k] = (Alpha[(l-1)*C->n_states+k] * fact2 + fact1) * Emissions[2*C->Vpoly[l]+C->Hpoly[l*C->n_states+k]];
			//	AlphaSum[l] += Alpha[l*C->n_states+k];
			//}
		}
		//loglik += log (AlphaSum[l]);
	}
	//cout << "F = " << loglik << endl;
}

void haplotype_hmm::backward(vector < float > & HL, vector < float > & HP) {
	//double loglik = 0.0, fact1, betaSumPrev, betaSumCurr, prob0, prob1;
	double betaSumPrev, betaSumCurr, betaSumCurr0, betaSumCurr1, betaSumCurr2, betaSumCurr3;
	double prob0_0, prob0_1, prob0_2, prob0_3;
	double prob1_0, prob1_1, prob1_2, prob1_3;

	double emit[2], hemit[2][2], nextemit[2], pcopy[2];
	vector < float > beta = vector < float > (C->n_states, 1.0);
	for (int l = C->n_sites-1 ; l >= 0 ; l --) {
		//Set up values
		prob0_0 = 0.0;prob0_1 = 0.0;prob0_2 = 0.0;prob0_3 = 0.0;
		prob1_0 = 0.0;prob1_1 = 0.0;prob1_2 = 0.0;prob1_3 = 0.0;
		emit[0] = Emissions[2*C->Vpoly[l] + 0];
		emit[1] = Emissions[2*C->Vpoly[l] + 1];
		hemit[0][0] = C->ee * HL[2*C->Vpoly[l]+0] / emit[0];
		hemit[0][1] = C->ed * HL[2*C->Vpoly[l]+0] / emit[1];
		hemit[1][0] = C->ed * HL[2*C->Vpoly[l]+1] / emit[0];
		hemit[1][1] = C->ee * HL[2*C->Vpoly[l]+1] / emit[1];

		// Backward pass
		betaSumCurr = 0.0;betaSumCurr0 = 0.0;betaSumCurr1 = 0.0;betaSumCurr2 = 0.0;betaSumCurr3 = 0.0;
		if (l == C->n_sites - 1) {
			for (int k = 0 ; k < C->n_states ; k ++) {
				prob0_0 += Alpha[l*C->n_states+k] * beta[k] * hemit[0][C->Hpoly[l*C->n_states+k]];
				prob1_0 += Alpha[l*C->n_states+k] * beta[k] * hemit[1][C->Hpoly[l*C->n_states+k]];
				betaSumCurr += emit[C->Hpoly[l*C->n_states+k]];
			}
		} else {
			pcopy[0] = C->nt[l] * nextemit[0] / betaSumPrev;
			pcopy[1] = C->nt[l] * nextemit[1] / betaSumPrev;

			int i = 0;
			int repeat = (C->n_states/4);
			int left = (C->n_states%4);

			while (repeat --) {
				beta[i+0] = beta[i+0] * pcopy[C->Hpoly[(l+1)*C->n_states+i+0]] + C->t[l];
				beta[i+1] = beta[i+1] * pcopy[C->Hpoly[(l+1)*C->n_states+i+1]] + C->t[l];
				beta[i+2] = beta[i+2] * pcopy[C->Hpoly[(l+1)*C->n_states+i+2]] + C->t[l];
				beta[i+3] = beta[i+3] * pcopy[C->Hpoly[(l+1)*C->n_states+i+3]] + C->t[l];
				prob0_0 += Alpha[l*C->n_states+i+0] * beta[i+0] * hemit[0][C->Hpoly[l*C->n_states+i+0]];
				prob1_0 += Alpha[l*C->n_states+i+0] * beta[i+0] * hemit[1][C->Hpoly[l*C->n_states+i+0]];
				prob0_1 += Alpha[l*C->n_states+i+1] * beta[i+1] * hemit[0][C->Hpoly[l*C->n_states+i+1]];
				prob1_1 += Alpha[l*C->n_states+i+1] * beta[i+1] * hemit[1][C->Hpoly[l*C->n_states+i+1]];
				prob0_2 += Alpha[l*C->n_states+i+2] * beta[i+2] * hemit[0][C->Hpoly[l*C->n_states+i+2]];
				prob1_2 += Alpha[l*C->n_states+i+2] * beta[i+2] * hemit[1][C->Hpoly[l*C->n_states+i+2]];
				prob0_3 += Alpha[l*C->n_states+i+3] * beta[i+3] * hemit[0][C->Hpoly[l*C->n_states+i+3]];
				prob1_3 += Alpha[l*C->n_states+i+3] * beta[i+3] * hemit[1][C->Hpoly[l*C->n_states+i+3]];
				betaSumCurr0 += emit[C->Hpoly[l*C->n_states+i+0]] * beta[i+0];
				betaSumCurr1 += emit[C->Hpoly[l*C->n_states+i+1]] * beta[i+1];
				betaSumCurr2 += emit[C->Hpoly[l*C->n_states+i+2]] * beta[i+2];
				betaSumCurr3 += emit[C->Hpoly[l*C->n_states+i+3]] * beta[i+3];
				i += 4;
			}

			switch (left) {
			case 3:	beta[i+2] = beta[i+2] * pcopy[C->Hpoly[(l+1)*C->n_states+i+2]] + C->t[l];
					prob0_2 += Alpha[l*C->n_states+i+2] * beta[i+2] * hemit[0][C->Hpoly[l*C->n_states+i+2]];
					prob1_2 += Alpha[l*C->n_states+i+2] * beta[i+2] * hemit[1][C->Hpoly[l*C->n_states+i+2]];
					betaSumCurr2 += emit[C->Hpoly[l*C->n_states+i+2]] * beta[i+2];
			case 2:	beta[i+1] = beta[i+1] * pcopy[C->Hpoly[(l+1)*C->n_states+i+1]] + C->t[l];
					prob0_1 += Alpha[l*C->n_states+i+1] * beta[i+1] * hemit[0][C->Hpoly[l*C->n_states+i+1]];
					prob1_1 += Alpha[l*C->n_states+i+1] * beta[i+1] * hemit[1][C->Hpoly[l*C->n_states+i+1]];
					betaSumCurr1 += emit[C->Hpoly[l*C->n_states+i+1]] * beta[i+1];
			case 1:	beta[i+0] = beta[i+0] * pcopy[C->Hpoly[(l+1)*C->n_states+i+0]] + C->t[l];
					prob0_0 += Alpha[l*C->n_states+i+0] * beta[i+0] * hemit[0][C->Hpoly[l*C->n_states+i+0]];
					prob1_0 += Alpha[l*C->n_states+i+0] * beta[i+0] * hemit[1][C->Hpoly[l*C->n_states+i+0]];
					betaSumCurr0 += emit[C->Hpoly[l*C->n_states+i+0]] * beta[i+0];
			}

			betaSumCurr = betaSumCurr0 + betaSumCurr1 + betaSumCurr2 + betaSumCurr3;

			//for (int k = 0 ; k < C->n_states ; k ++) {
			//	beta[k] = beta[k] * pcopy[C->Hpoly[(l+1)*C->n_states+k]] + C->t[l];
			//	prob0 += Alpha[l*C->n_states+k] * beta[k] * hemit[0][C->Hpoly[l*C->n_states+k]];
			//	prob1 += Alpha[l*C->n_states+k] * beta[k] * hemit[1][C->Hpoly[l*C->n_states+k]];
			//	betaSumCurr += emit[C->Hpoly[l*C->n_states+k]] * beta[k];
			//}
		}

		double sum0 = prob0_0 + prob0_1 + prob0_2 + prob0_3;
		double sum1 = prob1_0 + prob1_1 + prob1_2 + prob1_3;
		HP[2*C->Vpoly[l]+0] = sum0 / (sum0 + sum1);
		HP[2*C->Vpoly[l]+1] = sum1 / (sum0 + sum1);

		betaSumPrev = betaSumCurr / C->n_states;
		nextemit[0] = emit[0];
		nextemit[1] = emit[1];
		//loglik += log (betaSumPrev);
	}
	//cout << "B = " << loglik << endl;
	for (int l = 0 ; l < C->Vmono.size() ; l ++) {
		double prob0 = 0.0, prob1 = 0.0;
		if (!C->Hmono[C->Vmono[l]]) {
			prob0 = C->ee * HL[2*C->Vmono[l]+0];
			prob1 = C->ed * HL[2*C->Vmono[l]+1];
		} else {
			prob0 = C->ed * HL[2*C->Vmono[l]+0];
			prob1 = C->ee * HL[2*C->Vmono[l]+1];
		}
		//cout << prob0 / (prob0+prob1) << " " << prob1 / (prob0+prob1) << endl;
		HP[2*C->Vmono[l]+0] = prob0 / (prob0 + prob1);
		HP[2*C->Vmono[l]+1] = prob1 / (prob0 + prob1);
	}
}
*/
