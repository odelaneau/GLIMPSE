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

class genotype {
public:
	// INTERNAL DATA
	string name;
	unsigned int index;					// Index in containers
	unsigned int n_variants;			// Number of variants	(to iterate over Variants)

	// This is using lot of memory.
	// Solution:	Let's only store non-zero values (almost 33% memory reduction).
	//				Use a vector < bool > to flag zeros and non-zeros
	vector < unsigned char > GL;		// Original Genotype Likelihoods

	// These are using lot of memory.
	// Solution:	We can greatly compress these using sparse representation (massive memory reduction)
	// 				GP is mostly filled with triplets (1.0,0.0,0.0) and HAP with 0s (= all sampled haps are hom REF).
	// 				I propose to have a vector < bool > that flag sites at which storage is done.
	//				The big issue is that when a new site requires storage, this involves a memory block shift of all the next ones (we may use burnin-in iterations to spot these in advance)...
	//				In addition, we only have to store HAP at sites that are not GP=(1.0,0.0,0.0), as it's only zeros in this case.
	vector < float > GP;				// Output Genotype posteriors
	vector < int > HAP;					// Storing haplotypes
	int nHAPstored;

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
	void storeGenotypePosteriors(vector < float > &, vector < float > &);
	void storeSampledHaplotypes();

	void normGenotypePosteriors(int);
	void inferGenotypes();
};

#endif
