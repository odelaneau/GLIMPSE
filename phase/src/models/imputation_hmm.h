/*******************************************************************************
 * Copyright (C) 2022-2023 Simone Rubinacci
 * Copyright (C) 2022-2023 Olivier Delaneau
 *
 * MIT Licence
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
#define SIMDE_ENABLE_NATIVE_ALIASES
#include <simde/x86/avx2.h>
#include <simde/x86/fma.h>
#include <boost/align/aligned_allocator.hpp>

template <typename T>
using aligned_vector32 = std::vector<T, boost::alignment::aligned_allocator < T, 32 > >;

class imputation_hmm {
private:
	conditioning_set * C;
	unsigned int modK;

	// The backward pass needs every forward row, polymorphic sites x states floats in total. With
	// checkpoint_block > 0 only the first row of each block of that many sites is kept, and backward
	// recomputes each block when it reaches it: roughly one extra forward pass for far less memory.
	// 0 (or any value at least the number of sites) keeps every row in one block, so nothing is recomputed.
	const unsigned int checkpoint_block;
	unsigned int alpha_block;			// Sites per block for the current sample; set by resize()

	//DYNAMIC ARRAYS
	aligned_vector32 < float > Emissions;
	aligned_vector32 < float > AlphaCheckpoints;	// Forward row at the start of each block
	aligned_vector32 < float > AlphaBlock;			// Forward rows of the block currently being used
	aligned_vector32 < float > AlphaSum;
	aligned_vector32 < float > Beta;

	void forwardSite(const int, const std::vector < bool > &, const float *, float *);
	void recomputeBlock(const int, const std::vector < bool > &);

public:
	//CONSTRUCTOR/DESTRUCTOR
	imputation_hmm(conditioning_set *, const unsigned int);
	~imputation_hmm();

	void resize();
	void init(const std::vector < float > &);
	void forward(std::vector < bool > &);
	void backward(const std::vector < float > &, std::vector < bool > &, std::vector < float > &);
	void computePosteriors(const std::vector < float > &, std::vector < bool > &, std::vector < float > &);

};

#endif
