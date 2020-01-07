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
#ifndef _HAPLOTYPE_SET_H
#define _HAPLOTYPE_SET_H

#include <utils/otools.h>

#include <containers/bitmatrix.h>
#include <containers/genotype_set.h>
#include <containers/conditioning_set.h>


class haplotype_set {
public:
	//DATA
	unsigned long n_site;							// #variants
	unsigned long n_hap;							// #haplotypes
	unsigned long n_targ;							// #haplotypes
	unsigned long n_ref;							// #reference haplotypes
	unsigned long n_haps_pbwt;
	bitmatrix H_opt_var;							// Bit matrix of haplotypes (variant first)

	//PBWT
	bool refonly_pbwt;
	bool pbwt_built;
	int pbwt_depth;
	int pbwt_modulo;
	std::vector < std::vector < int > > pbwt_arrays;
	std::vector < std::vector < int > > pbwt_indexes;

	int n_iterations_main;

	//CONSTRUCTOR/DESTRUCTOR/INITIALIZATION
	haplotype_set();
	~haplotype_set();

	//ROUTINES
	void updateHaplotypes(genotype_set &);

	//PBWT ROUTINES
	void initPositionalBurrowWheelerTransform(int, int);
	void updatePositionalBurrowWheelerTransform();

	//STATE SELECTION
	void selectPositionalBurrowWheelerTransform(int, conditioning_set *);
	void selectRandomRefOnly(const int, conditioning_set *);
	void selectFmIndexRefOnly(const int ind, std::vector< int >& idxH);
	void selectStandardFull(const int ind, std::vector< int >& idxH);

	//HS iteration records
	void recordHS(genotype_set & G, int iter);

	//UTILITY
	int getC_bounded(const int l) const;
};

#endif
