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
#include <containers/conditioning_set.h>
#include <containers/haplotype_set.h>
#include <containers/probability_set.h>

class haplotype_hmm {
private:
	haplotype_set * H;
	conditioning_set * C;
	probability_set * P;

	int hap;

	//DYNAMIC ARRAYS
	std::vector < float > Emissions;
	std::vector < float > Alpha;
	std::vector < float > AlphaSum;

	bool store_posteriors;
	float posterior_threshold;

public:
	//CONSTRUCTOR/DESTRUCTOR
	haplotype_hmm(haplotype_set * H,conditioning_set *, probability_set *);
	~haplotype_hmm();

	void resize();
	void init(std::vector < float > &);
	void forward();
	void backward(std::vector < float > &, std::vector < float > &);
	void computePosteriors(std::vector < float > &, std::vector < float > &, int, bool);
	void applyThreshold(int l);
};

#endif
