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
#ifndef _GENOTYPE_WRITER_H
#define _GENOTYPE_WRITER_H

#include <utils/otools.h>

#include <containers/variant_map.h>
#include <containers/haplotype_set.h>
#include <containers/genotype_set.h>
#include <containers/probability_set.h>

#include <io/genotype_stream.h>



class genotype_writer {
public:
	//DATA
	haplotype_set & H;
	genotype_set & G;
	variant_map & V;
	std::vector < probability_set * >& P;		// Thresholded posterior probabilities
	gmap_reader & readerGM;

	int * genotypes;
	float * dosages;
	float * posteriors;

	std::string region;
	const float maf_common;

	//for info score
	double esum;
	double e2sum;
	double fsum;

	//COUNTS
	unsigned long n_variants;
	unsigned long n_main_samples;
	unsigned long n_ref_samples;

	//CONSTRUCTORS/DESCTRUCTORS
	genotype_writer(haplotype_set &, genotype_set &, variant_map &, std::vector < probability_set * >&, gmap_reader & readerGM, std::string regions, float _maf_common);
	~genotype_writer();

	void writeGenotypesAndImpute(std::string funphased, std::string freference, std::string fname, int start, int stop);

	void phase_marker(const int l, const genotype_stream& input_stream);
	void impute_marker(const int l, const float w,  const genotype_stream& input_stream);
	void impute_marker_border(const int l, const genotype_stream& input_stream, const bool left);

	float get_linear_interp_weight(const int l, const int current_position);
	int max3gt(const float& a,const float& b, const float& c) const;

};

#endif
