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

#ifndef _REF_GENOTYPE_READER_H
#define _REF_GENOTYPE_READER_H

#include <utils/otools.h>

#include <containers/variant_map.h>
#include "../containers/ref_haplotype_set.h"

class ref_genotype_reader
{
public:
	//DATA
	ref_haplotype_set& H;
	variant_map& V;

	const std::string region;
	const float sparse_maf;
	const bool keep_mono;

	int n_ref_samples;
	std::vector<int> ploidy_ref_samples;

	//CONSTRUCTORS/DESCTRUCTORS
	ref_genotype_reader(ref_haplotype_set &, variant_map &, const std::string regions, const float _sparse_maf, const bool _keep_mono);
	~ref_genotype_reader();

	//IO
	void set_ploidy_ref(const int ngt_ref, const int* gt_arr_ref, const int ngt_arr_ref);
	void readRefPanel(std::string fref,int nthreads);
	void initReader(bcf_srs_t * sr, const std::string fref, int nthreads);
	void scanGenotypesCommon(bcf_srs_t * sr, const std::string fref, int ref_sr_n);
	void parseRefGenotypes(bcf_srs_t * sr, const std::string fref);
};

#endif
