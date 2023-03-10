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

#ifndef _GENOTYPE_READER_H
#define _GENOTYPE_READER_H

#include <utils/otools.h>

#include <containers/variant_map.h>
#include <containers/haplotype_set.h>
#include <containers/glimpse_mpileup.h>

class genotype_reader
{
public:
	//DATA
	haplotype_set & H;
	genotype_set & G;
	variant_map & V;
	glimpse_mpileup & M;

	const float sparse_maf;
	const bool inputGL;
	const bool impute_refonly;
	const bool keep_mono;
	const bool use_gl_indels;

	int n_ref_samples;
	std::set < std::string > initializing_samples;
	std::map < std::string, int> ploidy_samples;
	std::vector<int> ploidy_ref_samples;

	//CONSTRUCTORS/DESCTRUCTORS
	genotype_reader(haplotype_set &, genotype_set &, variant_map &, glimpse_mpileup & M, const float _sparse_maf,const bool _inputGL, const bool _impute_refonly, const bool keep_mono, const bool use_gl_indels);
	~genotype_reader();

	//IO
	//void readInitializingSamples(string);
	void readSamplesFilePloidy(std::string);
	void set_ploidy_ref(bcf_srs_t * sr, int id_ref);
	void set_ploidy_tar();

	void readGenotypes(std::string , std::string , int nthreads);

	void initReader(bcf_srs_t * sr, std::string& fmain, std::string& fref, int nthreads);
	void initReader(bcf_srs_t * sr, std::string& file, int nthreads);

	void scanGenotypes(bcf_srs_t * sr);
	void readTarGenotypes(std::string , int);
	//void readTarGenotypesValidation(string, string, int);

	void parseGenotypes(bcf_srs_t * sr);

	//void initReaderAndBAMs(bcf_srs_t * sr, string& fref, int nthreads);
	void readGenotypesAndBAMs(std::string funphased, int nthreads);
	void parseRefGenotypes(bcf_srs_t * sr);
	void scanGenotypesCommon(bcf_srs_t * sr, int ref_sr_n);
};

#endif
