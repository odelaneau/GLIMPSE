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

#ifndef _GENOTYPE_WRITER_H
#define _GENOTYPE_WRITER_H

#include <utils/otools.h>

#include <containers/variant_map.h>
#include <containers/haplotype_set.h>
#include <containers/genotype_set.h>
#include <caller/caller_header.h>

#ifdef __BGEN__
#include "genfile/bgen/bgen.hpp"
#endif

class genotype_writer {
public:
	//DATA
	const haplotype_set & H;
	const genotype_set & G;
	const variant_map & V;
	const caller & C;

	//CONSTRUCTORS/DESCTRUCTORS
	genotype_writer(const haplotype_set &, const genotype_set &, const variant_map &, const caller & C);
	~genotype_writer();

	void writeGenotypes(const std::string fname, OutputFormat oformat, OutputCompression ocompr, const int n_bits_bgen, const int n_main, const int n_threads, const std::string fai_fname) const;
	void update_header_from_fai(bcf_hdr_t * hdr, const std::string fai_fname) const;

#ifdef __BGEN__
	void writeGenotypesBgen(const std::string fname, OutputFormat oformat, OutputCompression ocompr, const int n_bits_bgen, const int n_main, const int n_threads) const;
	void initialiseBGEN(genfile::bgen::Context& context,OutputCompression ocompr, const size_t n_sites_unbuf) const;
	void initialiseStream(genfile::bgen::Context& context,std::ostream& oStream) const ;
	void set_sample_names_impl(genfile::bgen::Context& context,std::ostream& oStream) const;
	void update_offset_and_header_block(genfile::bgen::Context& context,std::ostream& oStream, uint32_t offset) const;
#endif
};

#endif
