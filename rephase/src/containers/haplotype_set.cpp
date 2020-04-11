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

#include <containers/haplotype_set.h>

haplotype_set::haplotype_set() {
	n_site = 0;
	n_hap = 0;
}

haplotype_set::~haplotype_set() {
	n_site = 0;
	n_hap = 0;
}

void haplotype_set::mapRareHets(float maf) {
	RHflag = vector < vector < bool > > (n_site, vector < bool > (n_hap/2, false));
	RHprob = vector < vector < float > > (n_hap/2);
	for (int l = 0 ; l  < V.vec_pos.size() ; l ++) {
		if (V.vec_pos[l]->getMAF() < maf) {
			CVflag.push_back(false);
			for (int h = 0 ; h < n_hap ; h += 2) {
				RHflag[l][h/2] = (H_opt_var.get(l,h+0) != H_opt_var.get(l,h+1));
				if (RHflag[l][h/2]) RHprob[h/2].push_back(0.0f);
			}
		} else {
			CVflag.push_back(true);
			CVidx.push_back(l);
		}
	}
}

void haplotype_set::initPositionalBurrowWheelerTransform(int _pbwt_depth, int _pbwt_modulo) {
	pbwt_depth = _pbwt_depth;
	pbwt_modulo = _pbwt_modulo;
	cond_states = vector < vector < int > > (n_hap);
	pbwt_array = vector < int > (n_hap, 0);
	pbwt_indexes = vector < int > (n_hap, 0);
}

void haplotype_set::updatePositionalBurrowWheelerTransform() {
	tac.clock();
	vector < int > A = vector < int >(n_hap, 0);
	vector < int > last_selected = vector < int > (n_hap*pbwt_depth*2, -1);	//to kept track of the last ones that were recently push_backed into cond_states
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
			for (int h = 0 ; h < n_hap ; h ++) {
				int ac, o, c, hc = 0, a = H_opt_var.get(l,h);
				int pbwt_idx = pbwt_indexes[h];

				//Haplotypes that are *BEFORE* in PBWT array
				o = 1; c = 0;
	            for (;;) {
	            	if (pbwt_idx-o >= 0) hc = pbwt_array[pbwt_idx-o];	// Check if we hit boundaries
	            	else break;
	            	ac = H_opt_var.get(l,hc);
	            	if (ac != a) break;		//exit if alleles are different
	            	if (h/2 != hc/2) {
	            		if (last_selected[h * 2 * pbwt_depth + c] != hc) {	// Storage happens only when the new conditioning hap is different of the previously selected for this depth
	            			last_selected[h * 2 * pbwt_depth + c] = hc;		// Update last found
	            			cond_states[h].push_back(hc);						// Store new guess
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
	            	if (h/2 != hc/2) {
	            		if (last_selected[h * 2 * pbwt_depth + pbwt_depth + c] != hc) { 	// Storage happens only when the new conditioning hap is different of the previously selected for this depth
	            			last_selected[h * 2 * pbwt_depth + pbwt_depth + c] = hc;		// Update last found (saved after those seen before in the PBWT array, i.e. index+pbwt_depth)
	            			cond_states[h].push_back(hc);									// Store new guess
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

