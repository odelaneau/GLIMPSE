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

haplotype_hmm::haplotype_hmm(haplotype_set * _H, conditioning_set * _C, probability_set* _P) {
	H = _H;
	C = _C;
	P = _P;
	//P = _P;
	Emissions = std::vector < float > (2*C->n_vars, 0.0);
	store_posteriors = false;
	posterior_threshold = 0.0f;
	hap=0;
}

haplotype_hmm::~haplotype_hmm() {
	Alpha.clear();
	AlphaSum.clear();
	Emissions.clear();
}

void haplotype_hmm::resize() {
	AlphaSum.resize(C->n_sites);
	Alpha.resize(C->n_sites * C->n_states);
	//posterior_threshold = std::min(0.005f, (float) 0.9999/C->n_states);
	posterior_threshold = 0.0f;
}

void haplotype_hmm::init(std::vector < float > & HL) {
	double p0, p1;
	for (int l = 0 ; l < C->n_vars ; l ++) {
		p0 = HL[2*l+0] * C->ee + HL[2*l+1] * C->ed;
		p1 = HL[2*l+0] * C->ed + HL[2*l+1] * C->ee;

		Emissions[2*l+0] = p0 / (p0+p1);
		Emissions[2*l+1] = p1 / (p0+p1);
	}
}

void haplotype_hmm::computePosteriors(std::vector < float > & HL, std::vector < float > & HP, int _hap, bool _store_posteriors) {
	store_posteriors = _store_posteriors;
	hap = _hap;
	resize();
	vrb.bullet("HAP Resize (" + stb.str(tac.rel_time()*1.0, 1) + "ms)");tac.clock();
	init(HL);
	vrb.bullet("HAP Init (" + stb.str(tac.rel_time()*1.0, 1) + "ms)");tac.clock();
	forward();
	vrb.bullet("HAP Forward (" + stb.str(tac.rel_time()*1.0, 1) + "ms)");tac.clock();
	backward(HL, HP);
	vrb.bullet("HAP Backward (" + stb.str(tac.rel_time()*1.0, 1) + "ms)");tac.clock();
}

void haplotype_hmm::forward() {
	//double loglik = 0.0, fact1, fact2;
	double fact1, fact2;

	for (int l = 0 ; l < C->n_sites ; l ++) {
		AlphaSum[l] = 0.0;
		if (l == 0) {
			fact1 = 1.0f / C->n_states;
			for (int k = 0 ; k < C->n_states ; k ++) {
				Alpha[k] = Emissions[2*l+H->H_opt_var.get(l,C->idxH[k])] * fact1;
				AlphaSum[l] += Alpha[k];
			}
		} else {
			fact1 = C->t[l-1] / C->n_states;
			fact2 = C->nt[l-1] / AlphaSum[l-1];

			/*for (int k = 0 ; k < idxH.size() ; k ++) ac += H_opt_var.get(l,idxH[k]);
			 * To accelerate it:
			 * 		For closeby variants , assume zero recombination rates, so that fact2=1.0 and fact1=0
			 * 		Then, rescale likelihoods so that P(G=0) = 1.0
			 * 		So you just have to iterate over non-zero alleles and multiply alpha values by P(G=1)
			 * 		This should greatly speed up things as zeros values are everywhere, especially when the reference panel is large!!!!
			 */

			for (int k = 0 ; k < C->n_states ; k ++) {
				Alpha[l*C->n_states+k] = (Alpha[(l-1)*C->n_states+k] * fact2 + fact1) * Emissions[2*l+H->H_opt_var.get(l,C->idxH[k])];
				AlphaSum[l] += Alpha[l*C->n_states+k];
			}
		}
		//loglik += log (AlphaSum[l]);
	}
	//cout << "F = " << loglik << endl;
}

void haplotype_hmm::backward(std::vector < float > & HL, std::vector < float > & HP) {
	//double loglik = 0.0, fact1, betaSumPrev, betaSumCurr, prob0, prob1;
	double betaSumNext, betaSumCurr;
	double prob0 = 0.0, prob1 = 0.0;

	double hemit[2][2], pcopy[2], prob[2], betaSumTmp[2];
	std::vector < float > beta = std::vector < float > (C->n_states, 1.0f);
	for (int l = C->n_sites-1 ; l >= 0 ; l --) {
		//Set up values
		//prob0 = 0.0; prob1 = 0.0;
		prob[0]=0.0;prob[1]=0.0;
		betaSumTmp[0]=0.0;betaSumTmp[1]=0.0;
		hemit[0][0] = C->ee * HL[2*l+0] / Emissions[2*l + 0];
		hemit[0][1] = C->ed * HL[2*l+0] / Emissions[2*l + 1];
		hemit[1][0] = C->ed * HL[2*l+1] / Emissions[2*l + 0];
		hemit[1][1] = C->ee * HL[2*l+1] / Emissions[2*l + 1];

		// Backward pass
		betaSumCurr = 0.0;
		if (l == C->n_sites - 1) {
			for (int k = 0 ; k < C->n_states ; k ++) {
				Alpha[l*C->n_states+k] *= beta[k];
				prob[H->H_opt_var.get(l,C->idxH[k])] += Alpha[l*C->n_states+k];
				betaSumTmp[H->H_opt_var.get(l,C->idxH[k])]++;
			}
		} else {
			pcopy[0] = C->nt[l] * Emissions[2*(l+1) + 0] / betaSumNext;
			pcopy[1] = C->nt[l] * Emissions[2*(l+1) + 1] / betaSumNext;

			for (int k = 0 ; k < C->n_states ; k ++) {
				beta[k] = beta[k] * pcopy[H->H_opt_var.get(l+1,C->idxH[k])] + C->t[l];
				Alpha[l*C->n_states+k] *= beta[k];
				prob[H->H_opt_var.get(l,C->idxH[k])] += Alpha[l*C->n_states+k];
				betaSumTmp[H->H_opt_var.get(l,C->idxH[k])]+=beta[k];
			}
		}

		AlphaSum[l]=prob[0]+prob[1];
		prob0 = prob[0]*hemit[0][0] + prob[1]*hemit[0][1];
		prob1 = prob[0]*hemit[1][0] + prob[1]*hemit[1][1];
		HP[2*l+0] = prob0 / (prob0 + prob1);
		HP[2*l+1] = prob1 / (prob0 + prob1);
		betaSumCurr = betaSumTmp[0]*Emissions[2*l + 0] + betaSumTmp[1]*Emissions[2*l + 1];
		betaSumNext = betaSumCurr / C->n_states;
		//loglik += log (betaSumPrev);

		if (store_posteriors) applyThreshold(l);
	}
	//cout << "B = " << loglik << endl;
	/*
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
	*/
}

void haplotype_hmm::applyThreshold(int l)
{
	if (P==nullptr) return;

	std::vector<int> newids;
	std::vector<float> newpL;
	std::vector<float> newpR;

	std::vector<int>& oldids = P->reference_haps[hap][l];
	std::vector<float>& oldpL = P->prob_stateL[hap][l];
	std::vector<float>& oldpR = P->prob_stateR[hap][l];

	//posterior_threshold
	double thrL, thrR, sumProbL, sumProbR,sumL, sumR;

    int lP1 = (l < C->n_sites-1) ? (l + 1) : l;

	thrL = posterior_threshold*AlphaSum[l];
	thrR = posterior_threshold*AlphaSum[lP1];
	sumL = AlphaSum[l];
	sumR = AlphaSum[lP1];

    int i = 0, k = 0;

    // Traverse both array
    while (i<oldids.size() && k < C->n_states)
    {
        if (oldids[i] < C->idxH[k])
        {
			newids.push_back(oldids[i]);
			newpL.push_back(oldpL[i]);
			newpR.push_back(oldpL[i]);
			i++;
		}
        else if (oldids[i] == C->idxH[k])
		{
			newids.push_back(oldids[i]);
			newpL.push_back(oldpL[i]+Alpha[l*C->n_states+k]/sumL);
			newpR.push_back(oldpR[i]+Alpha[lP1*C->n_states+k]/sumR);
			i++;
			k++;
		}
        else
        {
            if (thrL < Alpha[l*C->n_states+k] || thrR < Alpha[lP1*C->n_states+k])
            {
         		newids.push_back(C->idxH[k]);
    			newpL.push_back(Alpha[l*C->n_states+k]/sumL);
    			newpR.push_back(Alpha[lP1*C->n_states+k]/sumR);
            }
            k++;
        }
    }

    // Store remaining elements of first array
    while (i<oldids.size())
    {
    	newids.push_back(oldids[i]);
		newpL.push_back(oldpL[i]);
		newpR.push_back(oldpL[i]);
		i++;
    }

    // Store remaining elements of second array
    while (k < C->n_states)
    {
    	if (thrL < Alpha[l*C->n_states+k] || thrR < Alpha[lP1*C->n_states+k])
		{
			newids.push_back(C->idxH[k]);
			newpL.push_back(Alpha[l*C->n_states+k]/sumL);
			newpR.push_back(Alpha[lP1*C->n_states+k]/sumR);
		}
    	k++;
    }
    // swap vectors
	P->reference_haps[hap][l].swap(newids);
	P->prob_stateL[hap][l].swap(newpL);
	P->prob_stateR[hap][l].swap(newpR);
}
