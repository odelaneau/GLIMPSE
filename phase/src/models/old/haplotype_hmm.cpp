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

haplotype_hmm::haplotype_hmm(bitmatrix & _H, hmm_parameters & _M, int _n_vars, int _n_haps) : H(_H), M(_M) {
	n_haps = _n_haps;
	n_vars = _n_vars;
	Alpha = vector < vector < float > > (n_vars, vector < float > (n_haps, 0.0));
	AlphaSum = vector < float > (n_vars, 0.0);
}

haplotype_hmm::~haplotype_hmm() {
	n_haps = 0;
	n_vars = 0;
	Alpha.clear();
	AlphaSum.clear();
}

double haplotype_hmm::forward(vector < double > & HL) {
	double loglik = 0.0, fact1, fact2;
	for (int l = 0 ; l < n_vars ; l ++) {
		float sum0 = 0.0, sum1 = 0.0;
		if (l == 0) {
			fact1 = 1.0 / n_haps;
			for (int k = 0 ; k < n_haps ; k += 2) {
				Alpha[l][k+0] = HL[2*l + H.get(l, k+0)] * fact1;
				Alpha[l][k+1] = HL[2*l + H.get(l, k+1)] * fact1;
				sum0 += Alpha[l][k+0];
				sum1 += Alpha[l][k+1];
			}
		} else {
			fact1 = M.t[l-1] / n_haps;
			fact2 = M.nt[l-1] / AlphaSum[l-1];
			for (int k = 0 ; k < n_haps ; k += 2) {
				Alpha[l][k+0] = (Alpha[l-1][k+0] * fact2 + fact1) * HL[2*l + H.get(l, k+0)];
				Alpha[l][k+1] = (Alpha[l-1][k+1] * fact2 + fact1) * HL[2*l + H.get(l, k+1)];
				sum0 += Alpha[l][k+0];
				sum1 += Alpha[l][k+1];
			}
		}

		//cout << l << " " << AlphaSum[l] << endl;

		AlphaSum[l] = sum0 + sum1;
		loglik += log (AlphaSum[l]);
	}
	return loglik;
}

double haplotype_hmm::backward(vector < double > & HL, vector < double > & HP) {
	double loglik = 0.0, fact1, fact2, betaSumPrev;
	vector < float > beta = vector < float > (n_haps, 1.0);
	for (int l = n_vars-1 ; l >= 0 ; l --) {
		// Backward pass
		double sum0 = 0.0, sum1 = 0.0;
		if (l == n_vars - 1) {
			for (int k = 0 ; k < n_haps ; k += 2) {
				sum0 += HL[2*l + H.get(l, k+0)];
				sum1 += HL[2*l + H.get(l, k+1)];
			}
		} else {
			fact1 = M.nt[l] / betaSumPrev;
			for (int k = 0 ; k < n_haps ; k += 2) {
				beta[k+0] = beta[k+0] * fact1 * HL[2*(l+1) + H.get(l+1, k+0)] + M.t[l];
				beta[k+1] = beta[k+1] * fact1 * HL[2*(l+1) + H.get(l+1, k+1)] + M.t[l];
				sum0 += HL[2*l + H.get(l, k+0)] * beta[k+0];
				sum1 += HL[2*l + H.get(l, k+1)] * beta[k+1];
			}
		}

		// Expectation pass
		double prob0 = 0.0;
		double prob1 = 0.0;
		for (int k = 0 ; k < n_haps ; k ++) {
			if (H.get(l, k)) prob1 += Alpha[l][k] * beta[k];
			else prob0 += Alpha[l][k] * beta[k];
		}
		HP[2*l+0] = prob0 / (prob0 + prob1);
		HP[2*l+1] = prob1 / (prob0 + prob1);
		betaSumPrev = (sum0 + sum1) / n_haps;
		loglik += log ((sum0 + sum1) / n_haps);
	}
	return loglik;
}
