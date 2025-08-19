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

#define _DECLARE_TOOLBOX_HERE
#include "otools.h"
#include "variant_map.h"
#include "ref_haplotype_set.h"

/**
 * @class ref_genotype_reader
 * @brief Class for reading and processing reference genotype data.
 * 
 * This class manages reference haplotype sets and variant maps, allowing
 * initialization and parsing of reference genotype data from files and streams.
 */
class ref_genotype_reader
{
public:
	/// Reference to the haplotype set container
	ref_haplotype_set& H;

	/// Reference to the variant map container
	variant_map& V;

	/// Genomic region string (e.g., chromosome or region coordinates)
	const std::string region;

	/// Minor allele frequency threshold for sparse representation
	const float sparse_maf;

	/// Flag to keep monomorphic variants
	const bool keep_mono;

	/// Number of reference samples
	int n_ref_samples;

	/// Vector storing ploidy values for each reference sample
	std::vector<int> ploidy_ref_samples;

	/**
	 * @brief Constructor
	 * 
	 * Initializes the genotype reader with haplotype set, variant map, genomic region,
	 * sparse minor allele frequency threshold, and flag to keep monomorphic variants.
	 * 
	 * @param H Reference to haplotype set
	 * @param V Reference to variant map
	 * @param regions Genomic region string
	 * @param _sparse_maf Minor allele frequency threshold for sparse representation
	 * @param _keep_mono Whether to keep monomorphic variants (true/false)
	 */
	ref_genotype_reader(ref_haplotype_set &H, variant_map &V, const std::string regions, const float _sparse_maf, const bool _keep_mono);

	/// Destructor
	~ref_genotype_reader();

	/**
	 * @brief Set ploidy for reference samples
	 * 
	 * @param ngt_ref Number of genotype calls for reference samples
	 * @param gt_arr_ref Pointer to genotype array for reference samples
	 * @param ngt_arr_ref Size of genotype array
	 */
	void set_ploidy_ref(const int ngt_ref, const int* gt_arr_ref, const int ngt_arr_ref);

	/**
	 * @brief Read reference panel from file
	 * 
	 * @param fref File path of reference panel
	 * @param nthreads Number of threads to use
	 */
	void readRefPanel(std::string fref, int nthreads);

	/**
	 * @brief Initialize reader for parsing reference data
	 * 
	 * @param sr BCF stream reader pointer
	 * @param fref File path of reference panel
	 * @param nthreads Number of threads
	 */
	void initReader(bcf_srs_t * sr, const std::string fref, int nthreads);

	/**
	 * @brief Scan common genotypes from reference stream
	 * 
	 * @param sr BCF stream reader pointer
	 * @param fref File path of reference panel
	 * @param ref_sr_n Reference stream index or count
	 */
	void scanGenotypesCommon(bcf_srs_t * sr, const std::string fref, int ref_sr_n);

	/**
	 * @brief Parse genotypes from reference stream
	 * 
	 * @param sr BCF stream reader pointer
	 * @param fref File path of reference panel
	 */
	void parseRefGenotypes(bcf_srs_t * sr, const std::string fref);
};

#endif
