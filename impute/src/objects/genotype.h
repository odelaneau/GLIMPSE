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

class genotype {
public:
	// INTERNAL DATA
	std::string name;
	unsigned int index;					// Index in containers
	unsigned int n_variants;			// Number of variants	(to iterate over Variants)
	std::vector < unsigned char > GL;		// Original Genotype Likelihoods
	std::vector < float > GP;				// Output Genotype posteriors
	std::vector < bool > H0;					// First haplotype
	std::vector < bool > H1;					// Second haplotype

	//CORE METHODS
	genotype(int, int);
	~genotype();
	void allocate();
	void free();

	//PROB METHODS
	void initHaplotypeLikelihoods(std::vector < float > &);
	void sampleHaplotypeH0(std::vector < float > &);
	void sampleHaplotypeH1(std::vector < float > &);
	void makeHaplotypeLikelihoods(std::vector < float > &, bool);
	void storeGenotypePosteriors(std::vector < float > &, std::vector < float > &);

	void normGenotypePosteriors(int);
	void inferGenotypes();
};

#endif
