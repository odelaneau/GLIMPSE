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
#ifndef _HAPLOTYPE_HMM_H
#define _HAPLOTYPE_HMM_H

#include <utils/otools.h>
#include <containers/bitmatrix.h>
#include <objects/hmm_parameters.h>

class haplotype_hmm {
private:
	//EXTERNAL DATA
	bitmatrix & H;
	hmm_parameters & M;

	//COORDINATES & CONSTANTS
	unsigned int n_haps;
	unsigned int n_vars;

	//DYNAMIC ARRAYS
	vector < vector < float > > Alpha;
	vector < float > AlphaSum;

public:
	//CONSTRUCTOR/DESTRUCTOR
	haplotype_hmm(bitmatrix &, hmm_parameters &, int, int);
	~haplotype_hmm();

	double forward(vector < double > & HL);
	double backward(vector < double > & HL, vector < double > & HP);
};

#endif
