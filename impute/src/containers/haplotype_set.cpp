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
	n_ref_haps = 0;
	n_main_haps = 0;
	pbwt_modulo =1;
	pbwt_depth =1;
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
			H_opt_var.set(v, n_ref_haps+2*i+0, a0);
			H_opt_var.set(v, n_ref_haps+2*i+1, a1);
		}
	}
	vrb.bullet("HAP update (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}


void haplotype_set::initPositionalBurrowWheelerTransform(int _pbwt_depth, int _pbwt_modulo) {
	pbwt_depth = _pbwt_depth;
	pbwt_modulo = _pbwt_modulo;
	pbwt_arrays = vector < vector < int > > ((n_site/pbwt_modulo) + 1, vector < int > (n_ref_haps, 0));
	pbwt_indexes = vector < vector < int > > ((n_site/pbwt_modulo) + 1, vector < int > (n_main_haps, 0));
}

void haplotype_set::updatePositionalBurrowWheelerTransform() {
	tac.clock();
	vector < int > A_all(n_hap, 0);
	vector < int > B_all(n_hap, 0);
	vector < int > A_ref(n_ref_haps, 0);
	vector < int > B_ref(n_ref_haps, 0);
	vector < int > idx(n_main_haps, 0);
	vector < bool > update_idx(n_main_haps, false);

	for (int l = 0, idx_store = 0 ; l < n_site ; l ++) {
		int u = 0, v = 0;
		int us = 0, vs = 0;
		bool is_copy_site = (l%pbwt_modulo) == 0;
		for (int h = 0 ; h < n_hap ; h ++) {
			int alookup = l?A_all[h]:h;
			if (!H_opt_var.get(l,alookup))
			{
				A_all[u++] = alookup;
				if (is_copy_site)
				{
					if (alookup < n_ref_haps) A_ref[us++] = alookup;
					else
					{
						idx[alookup-n_ref_haps] = us;
						update_idx[alookup-n_ref_haps] = false;
					}
				}
			}
			else
			{
				B_all[v++] = alookup;
				if (is_copy_site)
				{
					if (alookup < n_ref_haps) B_ref[vs++] = alookup;
					else
					{
						idx[alookup-n_ref_haps] = vs;
						update_idx[alookup-n_ref_haps] = true;
					}
				}
			}


		}
		std::copy(B_all.begin(), B_all.begin()+v, A_all.begin()+u);

		if (is_copy_site) {
			//std::copy(B_ref.begin(), B_ref.begin()+vs, A_ref.begin()+us);
			//std::copy(A_ref.begin(), A_ref.end(), pbwt_arrays[idx_store].begin());
			std::copy(A_ref.begin(), A_ref.begin() +us, pbwt_arrays[idx_store].begin());
			std::copy(B_ref.begin(), B_ref.begin() +vs, pbwt_arrays[idx_store].begin()+us);
			//for (int h = 0 ; h < n_hap ; h ++) pbwt_indexes[idx_store][A[h]] = h;
			//std::copy(idx.begin(), idx.end(), pbwt_arrays[idx_store].begin());
			for (int h = 0 ; h < n_main_haps ; h ++) pbwt_indexes[idx_store][h] = update_idx[h] ? idx[h] + u : idx[h];
			idx_store++;
		}
	}
	vrb.bullet("PBWT building (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::selectRandom(int K, conditioning_set * C) {
	vector < int > idxH;
	for (int h = 0 ; h < n_ref_haps ; h ++) if (initializing_haps[h]) idxH.push_back(h);
	if (n_ref_haps > K) {
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

void haplotype_set::selectPositionalBurrowWheelerTransform(int ind, int maxK, conditioning_set * C) {
	vector < bool > hapFlag = vector < bool > (n_ref_haps, false);
	for (int l = 0, idx_store = 0 ; l < n_site ; l += pbwt_modulo, idx_store++) {
		for (int p = 0 ; p < 2 ; p ++) {
			int h = 2*ind + p + n_ref_haps;
			int i = pbwt_indexes[idx_store][2*ind + p];

			//Backward
			bool ac;
            int o = 0, c = 0, hc = 0;
            bool a = H_opt_var.get(l,h);
            for (;;) {
            	if (i-o >= 0) hc = pbwt_arrays[idx_store][i-o];
            	else break;
            	ac = H_opt_var.get(l,hc);
            	if (ac != a) break;
            	hapFlag[hc] = (hc/2 != h/2);
            	c += (hc/2 != h/2);
            	if (c >= pbwt_depth-1) break;
            	o++;
            }

            //Forward
            o = 1; c = 0;
            for (;;) {
            	if (i+o < n_ref_haps) hc = pbwt_arrays[idx_store][i+o];
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
	for (int h = 0 ; h < n_ref_haps ; h ++) if (hapFlag[h]) idxH.push_back(h);

	if (idxH.size() > maxK) {
		random_shuffle(idxH.begin(), idxH.end());
		idxH.erase(idxH.begin() + maxK, idxH.end());
		sort(idxH.begin(), idxH.end());
	}

	C->clear();
	C->n_states = idxH.size();
	for (int l = 0 ; l < n_site ; l ++) {
		unsigned int ac = H_opt_var.get(l, 2*ind + n_ref_haps + 0) + H_opt_var.get(l, 2*ind + n_ref_haps + 1);
		C->Hmono.push_back(H_opt_var.get(l,idxH[0]));
		for (int k = 0 ; k < idxH.size() ; k ++) ac += H_opt_var.get(l,idxH[k]);
		if (ac > 0 && ac < (idxH.size()+2)) {
			C->Vpoly.push_back(l);
			for (int k = 0 ; k < idxH.size() ; k ++) C->Hpoly.push_back(H_opt_var.get(l,idxH[k]));
		} else C->Vmono.push_back(l);
	}
	C->n_sites = C->Vpoly.size();
	C->updateTransitions();
}
