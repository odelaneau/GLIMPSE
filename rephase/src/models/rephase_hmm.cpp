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

#include <models/rephase_hmm.h>

rephase_hmm::rephase_hmm(haplotype_set & _H, hmm_parameters & _M) : H(_H), M(_M) {
	target_ind = 0;
	idxK.clear();
}

rephase_hmm::~rephase_hmm() {
	Alpha.clear();
	AlphaSum.clear();
	idxK.clear();
}

void rephase_hmm::updateAndResize(unsigned int _ti) {
	target_ind = _ti;

	//Update idxK
	idxK.clear();
	for (int h = 0 ; h < H.cond_states[2*target_ind+0].size() ; h++) idxK.push_back(H.cond_states[2*target_ind+0][h]);
	for (int h = 0 ; h < H.cond_states[2*target_ind+1].size() ; h++) idxK.push_back(H.cond_states[2*target_ind+1][h]);
	sort(idxK.begin(), idxK.end());
	idxK.erase(unique(idxK.begin(), idxK.end()), idxK.end());

	//Update idxV
	idxV.clear();
	n_scaffolded = 0;
	typeV = vector < unsigned char > (H.CVflag.size(), VAR_HOM);
	for (int l = 0 ; l < H.CVflag.size() ; l++) {
		if (H.CVflag[l]) {
			typeV[l] = VAR_SCA;
			n_scaffolded++;
		} else if (H.RHflag[l][target_ind]) {
			typeV[l] = VAR_HET;
		}
	}

	//Resize memory allocations
	Alpha.resize(idxK.size() * n_scaffolded);
	AlphaSum.resize(n_scaffolded);
	hapProb = vector < vector < float > > (2, vector < float > (typeV.size(), 0.0f));
}

void rephase_hmm::forward(bool second) {
	double fact1, fact2;
	float emit [2];
	for (int al = 0, rl = 0 ; al < typeV.size() ; al ++) {
		if (idxV[l] == VAR_SCA) {
			AlphaSum[rl] = 0.0;
			bool tar_a = H.H_opt_var.get(al, 2*target_ind + second);
			emit[0] = tar_a?M.ed:M.ee;
			emit[1] = tar_a?M.ee:M.ed;

			if (rl == 0) {
				fact1 = 1.0f / idxK.size();
				for (int k = 0 ; k < idxK.size() ; k ++) {
					bool hid_a = H.H_opt_var.get(l, idxK[k]);
					Alpha[k] = emit[hid_a] * fact1;
					AlphaSum[rl] += Alpha[k];
				}
			} else {
				fact1 = C->t[rl-1] / C->n_states;
				fact2 = C->nt[rl-1] / AlphaSum[rl-1];
				for (int k = 0 ; k < C->n_states ; k ++) {
					bool hid_a = H.H_opt_var.get(al, idxK[k]);
					Alpha[rl*idxK.size()+k] = (Alpha[(rl-1)*idxK.size()+k] * fact2 + fact1) * emit[hid_a];
					AlphaSum[rl] += Alpha[rl*idxK.size()+k];
				}
			}
			rl ++;
		}
	}
}

void rephase_hmm::backward(bool second) {
	double fact1, fact2;
	float emit [2];
	vector < float > beta = vector < float > (idxK.size(), 1.0f);
	float betaSum = idxK.ize();
	for (int al = typeV.size() - 1 , sl = n_scaffolded - 1 ; al >= 0 ; al --) {
		if (idxV[al] == VAR_HET) {

		} else if (idxV[al] == VAR_SCA) {


			sl --;
		}
	}
}



		if (idxV[l] == VAR_SCA) {
			AlphaSum[rl] = 0.0;
			bool tar_a = H.H_opt_var.get(al, 2*target_ind + second);
			emit[0] = tar_a?M.ed:M.ee;
			emit[1] = tar_a?M.ee:M.ed;

			if (rl == 0) {
				fact1 = 1.0f / idxK.size();
				for (int k = 0 ; k < idxK.size() ; k ++) {
					bool hid_a = H.H_opt_var.get(l, idxK[k]);
					Alpha[k] = emit[hid_a] * fact1;
					AlphaSum[rl] += Alpha[k];
				}
			} else {
				fact1 = C->t[rl-1] / C->n_states;
				fact2 = C->nt[rl-1] / AlphaSum[rl-1];
				for (int k = 0 ; k < C->n_states ; k ++) {
					bool hid_a = H.H_opt_var.get(al, idxK[k]);
					Alpha[rl*idxK.size()+k] = (Alpha[(rl-1)*idxK.size()+k] * fact2 + fact1) * emit[hid_a];
					AlphaSum[rl] += Alpha[rl*idxK.size()+k];
				}
			}
			rl ++;
		}
	}
}






/*
 * I AM HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 */

void rephase_hmm::backward(vector < float > & HL, vector < float > & HP) {
	double emit[2][2];
	vector < float > beta0 = vector < float > (idxK.size(), 1.0);
	vector < float > beta1 = vector < float > (idxK.size(), 1.0);
	for (int l = idxV.size()-1 ; l >= 0 ; l --) {
		if (idxV[l] == VAR_SCA) {

		} else if (idxV[l] == VAR_SCA) {

		bool tar_a0 = H.H_opt_var.get(idxV[l], 2*target_ind + 0);
		bool tar_a1 = H.H_opt_var.get(idxV[l], 2*target_ind + 1);
		emit[0][0] = tar_a0?M.ed:M.ee;
		emit[0][1] = tar_a0?M.ee:M.ed;
		emit[1][0] = tar_a1?M.ed:M.ee;
		emit[1][1] = tar_a1?M.ee:M.ed;

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
