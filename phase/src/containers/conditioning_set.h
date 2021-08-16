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

#ifndef _CONDITIONING_SET_H
#define _CONDITIONING_SET_H

#include <utils/otools.h>
#include <containers/variant_map.h>
#include <containers/haplotype_set.h>

#define _SET8(n,i)	(n |= 1 << i)
#define _CLR8(n,i)	(n &= ~(1 << i))
#define _GET8(n,i)	((n >> i) & 1)

class conditioning_set {
public:
	//FIXED DATA
	const variant_map & mapG;
	const haplotype_set & H;
	const unsigned int n_haps;
	const unsigned int n_vars;
	const unsigned int n_effective;

	//CONDITIONING STATES
	vector < int > idxH;
	vector < vector < unsigned char > > Hpoly;
	vector < bool > Hmono;
	vector < unsigned int > Vpoly;
	vector < unsigned int > Vmono;
	unsigned int n_states;
	unsigned int n_sites;

	//TRANSITION PROBABILITIES
	vector < float > t;
	vector < float > nt;

	//EMISSION
	const float ed;
	const float ee;

	//CONSTRUCTOR/DESTRUCTOR/INITIALIZATION
	conditioning_set(const variant_map & _mapG, const haplotype_set & _H, const unsigned int _n_haps, const unsigned int _n_effective) :
		mapG(_mapG),
		H(_H),
		n_haps(_n_haps),
		n_vars(mapG.size()),
		n_effective(_n_effective),
		ed(0.0001f),
		ee(0.9999f) {
			Hmono.clear();
			Hpoly.clear();
			Vmono.clear();
			Vpoly.clear();
			t.clear();
			nt.clear();
			n_states = 0;
			n_sites = 0;
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
		if (Vpoly.size() == 0) return;
		t = vector < float > (Vpoly.size() - 1, 0.0);
		nt = vector < float > (Vpoly.size() - 1, 0.0);
		for (int l = 1 ; l < Vpoly.size() ; l ++) {
			float distcm = mapG.vec_pos[Vpoly[l]]->cm - mapG.vec_pos[Vpoly[l-1]]->cm;
			if (distcm < 0.00001f) distcm = 0.00001f;
			float rho = 0.04f * n_effective * distcm;
			t[l-1] = -1.0f * expm1(-1.0f * rho / n_haps);
			nt[l-1] = 1.0f-t[l-1];
		}
	}

	// Build compact conditioning set, keeping track of polymorphic [need HMM compute] and monomorphic [no need of HMM compute] sites
	void compact(int ind, const bool init) {
		clear();
		bool new_column = true;
		n_states = idxH.size();
		if (n_states == 0) {
			vrb.error("States for individual " + std::to_string(ind) + " are zero. Error during selection.  Try increasing the # samples in thef reference panel and/or the genomic region to consider.");
		}

		for (int l = 0 ; l < n_vars ; l ++) {
			unsigned int ac = H.H_opt_var.get(l, H.ind2hapid[ind] + H.n_ref_haps + 0);
			if (H.ploidy[ind] > 1) ac += H.H_opt_var.get(l, H.ind2hapid[ind] + H.n_ref_haps + 1);		// AC in target

			Hmono.push_back(H.H_opt_var.get(l,idxH[0]));																// Allele for monomorphic
			if (new_column) Hpoly.push_back(vector < unsigned char > (idxH.size()/8 + (idxH.size()%8>0), 0));
			else fill(Hpoly.back().begin(), Hpoly.back().end(), 0);

			for (int k = 0 ; k < idxH.size() ; k ++) if (H.H_opt_var.get(l,idxH[k])) { _SET8(Hpoly.back()[k/8], k%8); ac ++; };
			if (ac > 0 && ac < (idxH.size()+H.ploidy[ind])) {																		// Is the variant polymorphic? (using both cond. and target haps)
				Vpoly.push_back(l);											// If yes: store index
				new_column = true;
				//Hpoly.push_back(vector < unsigned char > (idxH.size()/8 + (idxH.size()%8>0), 0));
				//for (int k = 0 ; k < idxH.size() ; k ++) if (H.H_opt_var.get(l,idxH[k])) _SET(Hpoly.back()[k/8], k%8);
			} else {
				Vmono.push_back(l);																					// If no: declare it as monomorphic
				new_column = false;
			}
		}
		n_sites = Vpoly.size();		// Store number of polymorphic sites
	}

	// Select random set of reference haplotypes
	void selectRandom(int ind, int K) {
		//Selection
		idxH.clear();
		for (int h = 0 ; h < H.n_ref_haps ; h ++) if (H.initializing_haps[h]) idxH.push_back(h);		// The if is there when pooling is used
		if (H.n_ref_haps > K) {
			random_shuffle(idxH.begin(), idxH.end());
			idxH.erase(idxH.begin() + K, idxH.end());
			sort(idxH.begin(), idxH.end());
		}
		compact(ind, true);	//Compact it!
		updateTransitions();	// Update HMM parameters
	}

	// Select conditioning haplotypes using PBWT
	void selectPBWT(int ind, int K) {
		//Selection
		idxH.clear();
		//for (int h0 = 0 ; h0 < H.cond_states[H.ind2hapid[ind]+0].size() ; h0++) idxH.push_back(H.cond_states[H.ind2hapid[ind]+0][h0]);
		idxH.insert(idxH.end(), H.cond_states[H.ind2hapid[ind]+0].begin(), H.cond_states[H.ind2hapid[ind]+0].end());
		if (H.ploidy[ind] > 1)
			idxH.insert(idxH.end(), H.cond_states[H.ind2hapid[ind]+1].begin(), H.cond_states[H.ind2hapid[ind]+1].end());

		sort(idxH.begin(), idxH.end());
		idxH.erase(unique(idxH.begin(), idxH.end()), idxH.end());
		if (idxH.size() > K) {
			random_shuffle(idxH.begin(), idxH.end());
			idxH.erase(idxH.begin() + K, idxH.end());
			sort(idxH.begin(), idxH.end());
		}
		compact(ind, false);	//Compact it!
		updateTransitions();	// Update HMM parameters
	}
};

#endif

