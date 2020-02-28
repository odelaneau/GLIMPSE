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
#ifndef _GENOTYPE_H
#define _GENOTYPE_H

#include <utils/otools.h>

#define SET(n,i)	(n |= 1UL << i)
#define CLR(n,i)	(n &= ~(1UL << i));
#define GET(n,i)	((n >> i) & 1U);

// 5 bytes per genotype, but sparse!
struct inferred_genotype {
	float gp0, gp1, gp2;
	int hs, idx;

	inferred_genotype(int _idx, float _gp0, float _gp1, float _gp2, int _hs) : idx(_idx), gp0(_gp0), gp1(_gp1), gp2(_gp2), hs(_hs) {
	}

	bool operator<(const inferred_genotype & g) const {
		return idx < g.idx;
	}

	int infer() {
		if (gp0 >= gp1 && gp0 >= gp2) return 0;
		if (gp1 >= gp0 && gp1 >= gp2) return 1;
		if (gp2 >= gp0 && gp2 >= gp1) return 2;
		return 0;
	}
};

class genotype {
public:
	// INTERNAL DATA
	string name;
	unsigned int index;					// Index in containers
	unsigned int n_variants;			// Number of variants	(to iterate over Variants)

	// This is using lot of memory.
	// Solution:	Use a dictionnary given that GL triplets are highly repetitive.
	vector < unsigned char > GL;		// Original Genotype Likelihoods

	// Sparse storage for GPs and HSs
	vector < inferred_genotype > stored_data;
	int stored_cnt;

	// This should be removed as it is already stored in haplotype_set
	vector < bool > H0;					// First haplotype
	vector < bool > H1;					// Second haplotype


	//CORE METHODS
	genotype(int, int);
	~genotype();
	void allocate();
	void free();

	//PROB METHODS
	void initHaplotypeLikelihoods(vector < float > &);
	void sampleHaplotypeH0(vector < float > &);
	void sampleHaplotypeH1(vector < float > &);
	void makeHaplotypeLikelihoods(vector < float > &, bool);
	void storeGenotypePosteriorsAndHaplotypes(vector < float > &, vector < float > &);
	void sortAndNormAndInferGenotype();
};

#endif
