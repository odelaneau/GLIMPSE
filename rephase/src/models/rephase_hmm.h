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

#ifndef _REPHASE_HMM_H
#define _REPHASE_HMM_H

#include <utils/otools.h>
#include <containers/haplotype_set.h>
#include <objects/hmm_parameters.h>

#define VAR_SCA	0
#define VAR_HET	1
#define VAR_HOM	2

class rephase_hmm {
private:
	haplotype_set  & H;
	hmm_parameters & M;

	unsigned int target_ind;
	unsigned int n_scaffolded;
	vector < unsigned int > idxK;
	vector < unsigned char > typeV;

	vector < float > Alpha, AlphaSum;
	vector < vector < float > > hapProb;

public:
	//CONSTRUCTOR/DESTRUCTOR
	rephase_hmm(haplotype_set &, hmm_parameters &);
	~rephase_hmm();

	void updateAndResize(unsigned int);
	void forward();
	void backward();
};


#endif
