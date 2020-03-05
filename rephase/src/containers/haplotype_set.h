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

struct rare_het {
	unsigned int ind_idx;
	unsigned int var_idx;
	float prob;
	bool original_a0;
	bool sampled_a0;

	rare_het(unsigned int ii, unsigned int vi, bool a0) {
		ind_idx = ii;
		var_idx = vi;
		prob = 0.0f;
		original_a0 = a0;
		sampled_a0 = a0;
	}
};

class haplotype_set {
public:
	//DATA
	unsigned long n_site;							// #variants
	unsigned long n_hap;							// #haplotypes
	bitmatrix H_opt_var;							// Bit matrix of haplotypes (variant first)

	//PBbWT
	int pbwt_depth, pbwt_modulo;
	vector < int > pbwt_array;
	vector < int > pbwt_indexes;
	vector < vector < int > > cond_states;

	//RARE HETS REPHASING
	vector < rare_het > RH;
	vector < bool > RHflag;

	//CONSTRUCTOR/DESTRUCTOR/INITIALIZATION
	haplotype_set();
	~haplotype_set();

	//PBWT ROUTINES
	void initPositionalBurrowWheelerTransform(int, int);
	void updatePositionalBurrowWheelerTransform();

	//RARE HETS ROUTINES
	void mapRareHets(float max_maf);
	void updateapsGivenRareHets();
};

#endif
