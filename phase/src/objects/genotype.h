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

#ifndef _GENOTYPE_H
#define _GENOTYPE_H

#include <utils/otools.h>

#define _SET32(n,i)	((n) |= 1U << (i))
#define _CLR32(n,i)	((n) &= ~(1U << (i)));
#define _GET32(n,i)	(((n) >> (i)) & 1U);

// 5 bytes per genotype, but sparse!
struct inferred_genotype {
	float gp0, gp1;
	int32_t hs, idx;

	inferred_genotype(const int _idx, const float _gp0, const float _gp1, const int _hs) : idx(_idx), gp0(_gp0), gp1(_gp1), hs(_hs) {
	}

	bool operator<(const inferred_genotype & g) const {
		return idx < g.idx;
	}

	int infer() const {
		float gp2 = 1.0f - gp1 - gp0;
		if (gp0 > gp1 && gp0 > gp2) return 0;
		if (gp1 > gp0 && gp1 > gp2) return 1;
		if (gp2 > gp0 && gp2 > gp1) return 2;
		return 0;
	}

	float getGp2() const {
		return std::max(1.0f - gp1 - gp0, 0.0f);
	}

	bool infer_haploid() const {
		return gp1 > gp0;
	}
};


class genotype {
public:
	// INTERNAL DATA
	const string name;
	const int index;					// Index in containers
	const int n_variants;			// Number of variants	(to iterate over Variants)
	const int ploidy;				//ploidy can be 1,2 (2 is default)
	const int hapid;
	int stored_cnt;

	// This is using lot of memory.
	// Solution:	Use a dictionnary given that GL triplets are highly repetitive.
	vector < unsigned char > GL;		// Original Genotype Likelihoods

	// Sparse storage for GPs and HSs
	vector < inferred_genotype > stored_data;

	// This should be removed as it is already stored in haplotype_set
	vector < bool > H0;					// First haplotype
	vector < bool > H1;					// Second haplotype


	//CORE METHODS
	genotype(const string _name, const int _index, const int _n_variants, const int _ploidy, const int _hapid);
	~genotype();
	void allocate();
	void free();

	//PROB METHODS
	void initHaplotypeLikelihoods(vector < float > &) const;
	void sampleHaplotypeH0(const vector < float > &);
	void sampleHaplotypeH1(const vector < float > &);
	void makeHaplotypeLikelihoods(vector < float > &, bool) const;
	void storeGenotypePosteriorsAndHaplotypes(const vector < float > &);
	void storeGenotypePosteriorsAndHaplotypes(const vector < float > &, const vector < float > &);
	void sortAndNormAndInferGenotype();
};

#endif
