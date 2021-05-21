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

#ifndef _GENOTYPE_READER_H
#define _GENOTYPE_READER_H

#include <utils/otools.h>

#include <containers/variant_map.h>
#include <containers/haplotype_set.h>

class genotype_reader {
public:
	//DATA
	haplotype_set & H;
	genotype_set & G;
	variant_map & V;
	const string region;
	const bool impute_missing;
	const bool inputGL;
	const bool exclude_repeated_samples;

	std::vector<string> main_sample_names;
	set < string > initializing_samples;
	map <string, int> ploidy_samples;
	vector<int> ploidy_ref_samples;
	int n_haploid;
	int n_diploid;

	bcf_hdr_t *hdr_ref;
	int* imap;

	//COUNTS
	unsigned long n_variants;
	unsigned long n_included;
	unsigned long n_missing;
	unsigned long n_main_samples;
	unsigned long n_ref_samples;

	//CONSTRUCTORS/DESCTRUCTORS
	genotype_reader(haplotype_set &, genotype_set &, variant_map &, string regions, const bool _impute_missing, const bool _inputGL, const bool _exclude_repeated_samples);
	~genotype_reader();
	void allocateGenotypes();

	//IO
	void readInitializingSamples(string);
	void readSamplesFilePloidy(string);
	void initReader(bcf_srs_t * sr, string& fmain, string& fref, int nthreads);
	void readGenotypes(string funphased, string fphased, int nthreads);
	void scanGenotypes(bcf_srs_t * sr);
	void parseGenotypes(bcf_srs_t * sr);
	void set_ploidy_ref(bcf1_t * line_ref);

};

#endif

