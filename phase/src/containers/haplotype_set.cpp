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
#include <containers/haplotype_set.h>

haplotype_set::haplotype_set() {
	n_site = 0;
	n_hap = 0;
}

haplotype_set::~haplotype_set() {
	n_site = 0;
	n_hap = 0;
}

void haplotype_set::updateHaplotypes(genotype_set & G) {
	tac.clock();
	for (int i = 0 ; i < G.n_ind ; i ++) {
		for (int v = 0 ; v < n_site ; v ++) {
			bool a0 = G.vecG[i]->H0[v];
			bool a1 = G.vecG[i]->H1[v];
			H_opt_var.set(v, n_ref+2*i+0, a0);
			H_opt_var.set(v, n_ref+2*i+1, a1);
		}
	}
	vrb.bullet("HAP update (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}


void haplotype_set::initPositionalBurrowWheelerTransform(int _pbwt_depth, int _pbwt_modulo) {
	pbwt_depth = _pbwt_depth;
	pbwt_modulo = _pbwt_modulo;
	pbwt_arrays = vector < vector < int > > ((n_site/pbwt_modulo) + 1, vector < int > (n_hap, 0));
	pbwt_indexes = vector < vector < int > > ((n_site/pbwt_modulo) + 1, vector < int > (n_hap, 0));
}

void haplotype_set::updatePositionalBurrowWheelerTransform() {
	tac.clock();
	vector < int > A = vector < int >(n_hap, 0);
	vector < int > B = vector < int >(n_hap, 0);
	for (int l = 0, idx_store = 0 ; l < n_site ; l ++) {
		int u = 0, v = 0;
		for (int h = 0 ; h < n_hap ; h ++) {
			int alookup = l?A[h]:h;
			if (!H_opt_var.get(l,alookup)) A[u++] = alookup;
			else B[v++] = alookup;
		}
		std::copy(B.begin(), B.begin()+v, A.begin()+u);
		if ((l%pbwt_modulo) == 0) {
			std::copy(A.begin(), A.end(), pbwt_arrays[idx_store].begin());
			for (int h = 0 ; h < n_hap ; h ++) pbwt_indexes[idx_store][A[h]] = h;
			idx_store++;
		}
	}
	vrb.bullet("PBWT building (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::selectRandom(int K, conditioning_set * C) {
	vector < int > idxH;
	for (int h = 0 ; h < n_ref ; h ++) idxH.push_back(h);
	if (n_ref > K) {
		random_shuffle(idxH.begin(), idxH.end());
		idxH.erase(idxH.begin() + K, idxH.end());
		sort(idxH.begin(), idxH.end());
	}

	C->clear();
	C->n_states = idxH.size();
	for (int l = 0 ; l < n_site ; l ++) {
		unsigned int ac = 0;
		C->Hmono.push_back(H_opt_var.get(l,idxH[0]));
		for (int k = 0 ; k < idxH.size() ; k ++) ac += H_opt_var.get(l,idxH[k]);
		if (ac > 0 && ac < idxH.size()) {
			C->Vpoly.push_back(l);
			for (int k = 0 ; k < idxH.size() ; k ++) C->Hpoly.push_back(H_opt_var.get(l,idxH[k]));
		} else C->Vmono.push_back(l);
	}
	C->n_sites = C->Vpoly.size();
	C->updateTransitions();
}

void haplotype_set::selectPositionalBurrowWheelerTransform(int ind, conditioning_set * C) {
	vector < bool > hapFlag = vector < bool > (n_hap, false);
	for (int l = 0, idx_store = 0 ; l < n_site ; l += pbwt_modulo, idx_store++) {
		for (int p = 0 ; p < 2 ; p ++) {
			int h = 2*ind + p + n_ref;
			int i = pbwt_indexes[idx_store][h];

			//Backward
			bool ac;
            int o = 1, c = 0, hc = 0;
            bool a = H_opt_var.get(l,h);
            for (;;) {
            	if (i-o >= 0) hc = pbwt_arrays[idx_store][i-o];
            	else break;
            	ac = H_opt_var.get(l,hc);
            	if (ac != a) break;
            	hapFlag[hc] = (hc/2 != h/2);
            	c += (hc/2 != h/2);
            	if (c >= pbwt_depth) break;
            	o++;
            }

            //Forward
            o = 1; c = 0;
            for (;;) {
            	if (i+o < n_hap) hc = pbwt_arrays[idx_store][i+o];
            	else break;
            	ac = H_opt_var.get(l,hc);
            	if (ac != a) break;
            	hapFlag[hc] = (hc/2 != h/2);
            	c += (hc/2 != h/2);
            	if (c >= pbwt_depth) break;
            	o++;
            }
		}
	}
	vector < int > idxH;
	for (int h = 0 ; h < n_hap ; h ++) if (hapFlag[h]) idxH.push_back(h);

	C->clear();
	C->n_states = idxH.size();
	//vector < int > freq = vector < int > (10, 0);
	for (int l = 0 ; l < n_site ; l ++) {
		unsigned int ac = H_opt_var.get(l, 2*ind + n_ref + 0) + H_opt_var.get(l, 2*ind + n_ref + 1);
		unsigned int ac2 = H_opt_var.get(l, 2*ind + n_ref + 0) + H_opt_var.get(l, 2*ind + n_ref + 1);
		C->Hmono.push_back(H_opt_var.get(l,idxH[0]));
		for (int k = 0 ; k < idxH.size() ; k ++) ac2 += H_opt_var.get(l,idxH[k]);
		ac += ac2;
		//freq[(ac2>4)?4:ac2]++;
		if (ac > 0 && ac < (idxH.size()+2)) {
			C->Vpoly.push_back(l);
			for (int k = 0 ; k < idxH.size() ; k ++) C->Hpoly.push_back(H_opt_var.get(l,idxH[k]));
		} else C->Vmono.push_back(l);
	}
	//cout << freq[0] << " " << freq[1] << " " << freq[2] << " " << freq[3] << " " << freq[4] << endl;
	C->n_sites = C->Vpoly.size();
	C->updateTransitions();
}
