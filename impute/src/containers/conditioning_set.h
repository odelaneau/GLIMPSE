/*******************************************************************************
 * Copyright (C) 2018 Olivier Delaneau, University of Lausanne
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
#ifndef _CONDITIONING_SET_H
#define _CONDITIONING_SET_H

#include <utils/otools.h>
#include <containers/variant_map.h>
#include <containers/haplotype_set.h>

/*
 * Do not ask me why i did not put all this in a cpp file! Lazyness ...
 */

class conditioning_set {
public:
	//FIXED DATA
	variant_map & mapG;
	haplotype_set & H;
	unsigned int n_haps;
	unsigned int n_vars;
	unsigned int n_effective;

	//CONDITIONING STATES
	vector < int > idxH;
	vector < bool > Hpoly;
	vector < bool > Hmono;
	vector < unsigned int > Vpoly;
	vector < unsigned int > Vmono;
	unsigned int n_states;
	unsigned int n_sites;

	//TRANSITION PROBABILITIES
	vector < float > t;
	vector < float > nt;

	//EMISSION
	double ed;
	double ee;

	//CONSTRUCTOR/DESTRUCTOR/INITIALIZATION
	conditioning_set(variant_map & _mapG, haplotype_set & _H, unsigned int _n_haps, unsigned int _n_effective) : mapG(_mapG), H(_H) {
		Hmono.clear();
		Hpoly.clear();
		Vmono.clear();
		Vpoly.clear();
		t.clear();
		nt.clear();
		n_states = 0;
		n_sites = 0;
		n_haps = _n_haps;
		n_vars = mapG.size();
		n_effective= _n_effective;
		ed = 0.0001;
		ee = 0.9999;
	}

	~conditioning_set() {
		Hmono.clear();
		Hpoly.clear();
		Vmono.clear();
		Vpoly.clear();
		t.clear();
		nt.clear();
		n_states = 0;
		n_sites = 0;
	}

	void clear() {
		Hmono.clear();
		Hpoly.clear();
		Vmono.clear();
		Vpoly.clear();
		t.clear();
		nt.clear();
		n_states = 0;
		n_sites = 0;
	}

	// Update transition probabilities of the HMM for this particular conditioning set
	void updateTransitions() {
		t = vector < float > (Vpoly.size() - 1, 0.0);
		nt = vector < float > (Vpoly.size() - 1, 0.0);
		for (int l = 1 ; l < Vpoly.size() ; l ++) {
			float distcm = mapG.vec_pos[Vpoly[l]]->cm - mapG.vec_pos[Vpoly[l-1]]->cm;
			if (distcm < 0.00001f) distcm = 0.00001f;
			float rho = 0.04 * n_effective * distcm;
			t[l-1] = -1.0 * expm1(-1.0 * rho / n_haps);
			nt[l-1] = 1-t[l-1];
		}
	}

	// Build compact conditioning set, keeping track of polymorphic [need HMM compute] and monomorphic [no need of HMM compute] sites
	void compact(int ind) {
		clear();
		n_states = idxH.size();
		for (int l = 0 ; l < n_vars ; l ++) {
			unsigned int ac = H.H_opt_var.get(l, 2*ind + H.n_ref + 0) + H.H_opt_var.get(l, 2*ind + H.n_ref + 1);		// AC in target
			Hmono.push_back(H.H_opt_var.get(l,idxH[0]));																// Allele for monomorphic
			for (int k = 0 ; k < idxH.size() ; k ++) ac += H.H_opt_var.get(l,idxH[k]);									// AC in conditioning states
			if (ac > 0 && ac < (idxH.size()+2)) {																		// Is the variant polymorphic? (using both cond. and target haps)
				Vpoly.push_back(l);																						// If yes: store index
				for (int k = 0 ; k < idxH.size() ; k ++) Hpoly.push_back(H.H_opt_var.get(l,idxH[k]));					// If yes: store alleles
			} else Vmono.push_back(l);																					// If no: declare it as monomorphic
		}
		n_sites = Vpoly.size();		// Store number of polymorphic sites
	}

	// Select random set of reference haplotypes
	void selectRandom(int ind, int K) {
		//Selection
		idxH.clear();
		for (int h = 0 ; h < H.n_ref ; h ++) if (H.initializing_haps[h]) idxH.push_back(h);		// The if is there when pooling is used
		if (H.n_ref > K) {
			random_shuffle(idxH.begin(), idxH.end());
			idxH.erase(idxH.begin() + K, idxH.end());
			sort(idxH.begin(), idxH.end());
		}
		compact(ind);	//Compact it!
		updateTransitions();	// Update HMM parameters
	}

	// Select conditioning haplotypes using PBWT
	void selectPBWT(int ind, int K) {
		//Selection
		idxH.clear();
		for (int h0 = 0 ; h0 < H.cond_states[2*ind+0].size() ; h0++) idxH.push_back(H.cond_states[2*ind+0][h0]);
		for (int h1 = 0 ; h1 < H.cond_states[2*ind+1].size() ; h1++) idxH.push_back(H.cond_states[2*ind+1][h1]);
		sort(idxH.begin(), idxH.end());
		idxH.erase(unique(idxH.begin(), idxH.end()), idxH.end());
		if (idxH.size() > K) {
			random_shuffle(idxH.begin(), idxH.end());
			idxH.erase(idxH.begin() + K, idxH.end());
			sort(idxH.begin(), idxH.end());
		}
		compact(ind);	//Compact it!
		updateTransitions();	// Update HMM parameters
	}
};

#endif
