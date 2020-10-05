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
	n_ref_haps = 0;
	pbwt_depth = 0;
	pbwt_modulo = 0;
	max_ploidy = 2;
	n_main_haps=0;
	fploidy=2;
}

haplotype_set::~haplotype_set() {
	n_site = 0;
	n_hap = 0;
}

void haplotype_set::updateHaplotypes(const genotype_set & G) {
	tac.clock();
	for (int i = 0 ; i < G.n_ind ; i ++) {
		const int ploidy = G.vecG[i]->ploidy;
		const int hapid = ind2hapid[i] + n_ref_haps;
		for (int v = 0 ; v < n_site ; v ++) {
			H_opt_var.set(v, hapid+0, G.vecG[i]->H0[v]);
			if (ploidy > 1) H_opt_var.set(v, hapid+1, G.vecG[i]->H1[v]);
		}
	}
	vrb.bullet("HAP update (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::initPositionalBurrowWheelerTransform(const int _pbwt_depth, const int _pbwt_modulo) {
	pbwt_depth = _pbwt_depth;
	pbwt_modulo = _pbwt_modulo;
	cond_states = vector < vector < int > > (n_main_haps);		// Storage of states for each target_hap
	pbwt_array = vector < int > (n_hap);
	pbwt_indexes = vector < int > (n_main_haps, 0);
}

void haplotype_set::updatePositionalBurrowWheelerTransform() {
	tac.clock();
	vector < int > A = vector < int >(n_hap, 0);

	int shift = rng.getInt(0,pbwt_modulo-1);

	std::iota(pbwt_array.begin(),pbwt_array.end(),0);
	random_shuffle(pbwt_array.begin(),pbwt_array.end());

	vector < int > last_selected = vector < int > (n_main_haps*pbwt_depth*2, -1);	//to kept track of the last ones that were recently push_backed into cond_states
	for (int e = 0 ; e < n_main_haps ; e ++) cond_states[e].clear();
	for (int l = 0 ; l < n_site ; l ++) {
		// Building PBWT arrays
		int u = 0, v = 0;
		for (int h = 0 ; h < n_hap ; h ++) H_opt_var.get(l,pbwt_array[h]) ? A[v++] = pbwt_array[h] : pbwt_array[u++] = pbwt_array[h];

		std::copy(A.begin(), A.begin()+v, pbwt_array.begin()+u);

		// Selecting using PBWT array
		// It might be better here to iterate over indexes as sorted in pbwt_arrays and to test is we hit a target haplotype.
		// Instead of iterating over target haplotypes (I guess it depends on the ratio #target/#reference).
		if ((l+shift)%pbwt_modulo == 0)
		{
			// Build reverse indexing
			for (int h = 0 ; h < n_hap ; h ++) if (pbwt_array[h] >= n_ref_haps) pbwt_indexes[pbwt_array[h]-n_ref_haps] = h;

			// Selecting conditioning haplotypes
			for (int htr = 0 ; htr < n_main_haps ; htr ++) {
				const int pbwt_idx = pbwt_indexes[htr];
				const int htrind = hapid2ind[htr+n_ref_haps];
				const bool a = pbwt_idx >= u;
				int* ptr_ls = last_selected.data() + (htr * 2 * pbwt_depth);

				//Haplotypes that are *BEFORE* in PBWT array
				for (int o = 1, c = 0; c < pbwt_depth; ++o) {
					const int idx = pbwt_idx-o;
					if (idx < 0 || ((idx >= u) != a)) break;	// Exit if we hit boundaries or if alleles are different
					const int hc = pbwt_array[idx];
					if (htrind != hapid2ind[hc])
	            	{
						if (ptr_ls[c] != hc) {	// Storage happens only when the new conditioning hap is different of the previously selected for this depth
							ptr_ls[c] = hc;		// Update last found
							cond_states[htr].push_back(hc);						// Store new guess
						}
						++c;
					}
				}

				//Haplotypes that are *AFTER* in PBWT array
				for (int o = 1, c = 0; c < pbwt_depth; ++o) {
					const int idx = pbwt_idx+o;
					if (idx >= n_hap || ((idx >= u) != a)) break;	// Exit if we hit boundaries or if alleles are different
					const int hc = pbwt_array[idx];
					if (htrind != hapid2ind[hc])
	            	{
	            		if (ptr_ls[c+pbwt_depth] != hc) { 	// Storage happens only when the new conditioning hap is different of the previously selected for this depth
	            			ptr_ls[c+pbwt_depth] = hc;		// Update last found (saved after those seen before in the PBWT array, i.e. index+pbwt_depth)
							cond_states[htr].push_back(hc);									// Store new guess
						}
						++c;
					}
				}
			}
		}
	}
	vrb.bullet("PBWT building & selection (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
