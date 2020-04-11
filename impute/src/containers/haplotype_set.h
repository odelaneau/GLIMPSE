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

#ifndef _HAPLOTYPE_SET_H
#define _HAPLOTYPE_SET_H

#include <utils/otools.h>

#include <containers/bitmatrix.h>
#include <containers/genotype_set.h>

class haplotype_set {
public:
	//DATA
	unsigned long n_site;							// #variants
	unsigned long n_hap;							// #reference+main haplotypes
	unsigned long n_ref;						// #reference haplotypes
	unsigned long n_main_haps;						// #main haplotypes
	bitmatrix H_opt_var;							// Bit matrix of haplotypes (variant first)
	vector < bool > initializing_haps;				// Reference haplotype pool (subset) used for initializing iteration

	//PBWT
	int pbwt_depth, pbwt_modulo;
	vector < int > pbwt_array;
	vector < int > pbwt_array_ref;
	vector < int > pbwt_indexes;
	vector < vector < int > > cond_states;

	//CONSTRUCTOR/DESTRUCTOR/INITIALIZATION
	haplotype_set();
	~haplotype_set();

	//ROUTINES
	void updateHaplotypes(genotype_set &);

	//PBWT ROUTINES
	void initPositionalBurrowWheelerTransform(int, int);
	void updatePositionalBurrowWheelerTransform();
};

#endif
