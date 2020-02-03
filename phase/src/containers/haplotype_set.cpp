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
	cond_states = vector < vector < int > > (n_hap - n_ref);		// Storage of states for each target_hap
	pbwt_array = vector < int > (n_hap, 0);
	pbwt_indexes = vector < int > (n_hap, 0);
}

void haplotype_set::updatePositionalBurrowWheelerTransform() {
	tac.clock();
	vector < int > A = vector < int >(n_hap, 0);
	vector < int > last_selected = vector < int > ((n_hap - n_ref)*pbwt_depth*2, -1);	//to kept track of the last ones that were recently push_backed into cond_states
	for (int e = 0 ; e < cond_states.size() ; e ++) cond_states[e].clear();
	for (int l = 0 ; l < n_site ; l ++) {
		// Building PBWT arrays
		int u = 0, v = 0;
		for (int h = 0 ; h < n_hap ; h ++) {
			int alookup = l?pbwt_array[h]:h;
			if (!H_opt_var.get(l,alookup)) pbwt_array[u++] = alookup;
			else A[v++] = alookup;
		}
		std::copy(A.begin(), A.begin()+v, pbwt_array.begin()+u);

		// Selecting using PBWT array
		// It might be better here to iterate over indexes as sorted in pbwt_arrays and to test is we hit a target haplotype.
		// Instead of iterating over target haplotypes (I guess it depends on the ratio #target/#reference).
		if ((l%pbwt_modulo) == 0) {
			// Build reverse indexing
			for (int h = 0 ; h < n_hap ; h ++) pbwt_indexes[pbwt_array[h]] = h;

			// Selecting conditioning haplotypes
			for (int ht = n_ref, htr = 0 ; ht < n_hap ; ht ++, htr ++) {
				int ac, o, c, hc = 0, a = H_opt_var.get(l,ht);
				int pbwt_idx = pbwt_indexes[ht];

				//Haplotypes that are *BEFORE* in PBWT array
				o = 1; c = 0;
	            for (;;) {
	            	if (pbwt_idx-o >= 0) hc = pbwt_array[pbwt_idx-o];	// Check if we hit boundaries
	            	else break;
	            	ac = H_opt_var.get(l,hc);
	            	if (ac != a) break;		//exit if alleles are different
	            	if (ht/2 != hc/2) {
	            		if (last_selected[htr * 2 * pbwt_depth + c] != hc) {	// Storage happens only when the new conditioning hap is different of the previously selected for this depth
	            			last_selected[htr * 2 * pbwt_depth + c] = hc;		// Update last found
	            			cond_states[htr].push_back(hc);						// Store new guess
	            		}
	            		c++;
	            		if (c >= pbwt_depth) break; 	// Check if we found enough states
	            	}
	            	o++;
	            }

	            //Haplotypes that are *AFTER* in PBWT array
	            o = 1; c = 0;
	            for (;;) {
	            	if (pbwt_idx+o < n_hap) hc = pbwt_array[pbwt_idx+o];	// Check if we hit boundaries
	            	else break;
	            	ac = H_opt_var.get(l,hc);
	            	if (ac != a) break;		//exit if alleles are different
	            	if (ht/2 != hc/2) {
	            		if (last_selected[htr * 2 * pbwt_depth + pbwt_depth + c] != hc) { 	// Storage happens only when the new conditioning hap is different of the previously selected for this depth
	            			last_selected[htr * 2 * pbwt_depth + pbwt_depth + c] = hc;		// Update last found (saved after those seen before in the PBWT array, i.e. index+pbwt_depth)
	            			cond_states[htr].push_back(hc);									// Store new guess
	            		}
	            		c++;
	            		if (c >= pbwt_depth) break; 	// Check if we found enough states
	            	}
	            	o++;
	            }
			}
		}
	}
	vrb.bullet("PBWT building & selection (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

/*

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

	if (idxH.size() > maxK) {
		random_shuffle(idxH.begin(), idxH.end());
		idxH.erase(idxH.begin() + maxK, idxH.end());
		sort(idxH.begin(), idxH.end());
	}

	C->clear();
	C->n_states = idxH.size();
	for (int l = 0 ; l < n_site ; l ++) {
		unsigned int ac = H_opt_var.get(l, 2*ind + n_ref + 0) + H_opt_var.get(l, 2*ind + n_ref + 1);
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
*/
