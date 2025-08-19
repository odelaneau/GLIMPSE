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

#ifndef _REF_HAPLOTYPE_SET_H
#define _REF_HAPLOTYPE_SET_H

#include "otools.h"
#include "checksum_utils.h"

#include "variant_map.h"
#include "bitmatrix.h"
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/vector.hpp"

static int p3decode[128] ;
#define ENCODE_MAX1 64		     /* ~64 */
#define ENCODE_MAX2 ((95-63) << 6)   /* ~1k - is this 32 or 31?*/
#define ENCODE_MAX3 ((127-96) << 11) /* ~64k - ditto */

// Code by Richard Durbin [https://github.com/richarddurbin/pbwt]
// doi: 10.1093/bioinformatics/btu014

/**
 * @brief Three-level run length encoding (pack3) decoding lookup initialization.
 *
 * pack3 encoding compresses runs of bits efficiently using 3 levels of run length
 * encoding to allow runs up to ~64k length to be stored in a few bytes.
 */
static void pack3init (void)
{
  int n ;
  for (n = 0 ; n < 64 ; ++n) p3decode[n] = n ;
  for (n = 64 ; n < 96 ; ++n) p3decode[n] = (n-64) << 6 ;
  for (n = 96 ; n < 128 ; ++n) p3decode[n] = (n-96) << 11 ;
}

/**
 * @brief Helper function to add encoded run-length data to compressed vector.
 * 
 * @param y Symbol value (bit 1 in top bit)
 * @param n Run length count
 * @param compressed_y Vector to append encoded bytes
 *
 * This function encodes runs of bits using 3 levels of run length encoding for efficient storage.
 */
static void pack3Add (unsigned char y, uint32_t n, std::vector<unsigned char>& compressed_y)
/* local utility for pack3 */
{
	y <<= 7 ;			/* first move the actual symbol to the top bit */

	while (n >= ENCODE_MAX3)
	{
		compressed_y.push_back(y | 0x7f);
		n -= ENCODE_MAX3;
	}
	if (n >= ENCODE_MAX2)
	{
		compressed_y.push_back(y | 0x60 | (n >> 11));
		n &= 0x7ff;
	}
	if (n >= ENCODE_MAX1)
	{
		compressed_y.push_back(y | 0x40 | (n >> 6));
		n &= 0x3f;
	}
	if (n) compressed_y.push_back(y | n);
}

/**
 * @brief Compress a vector of bits using three-level run length encoding (pack3).
 * 
 * @param uncompressed_y Input vector of bits to compress (values 0 or 1)
 * @param n_chars Number of characters (bits) in input vector
 * @param compressed_y Output vector to store compressed bytes
 *
 * This method applies run length encoding to the input bit vector, encoding
 * runs of identical bits efficiently with variable-length codes.
 */
static void pack3(const std::vector<unsigned char>& uncompressed_y, int n_chars, std::vector<unsigned char>& compressed_y)
{
	int m = 0, m0 = 0, i=0;
	while (m < n_chars)
	{
		unsigned char y = uncompressed_y[i++]; //take a symbol [0 or 1]
		m0 = m++; //m itererates over all symbols of the same type
		while (i < n_chars && uncompressed_y[i] == y)
		{
			m++;
			i++;
		}
		pack3Add(y, m-m0, compressed_y);
	}
}

/**
 * @class ref_haplotype_set
 * @brief Container for storing reference haplotype data including sparse and common variant representations.
 *
 * This class maintains data on variants and haplotypes for a reference panel,
 * including bitmatrix representations for common variants, sparse representations
 * for rare variants, and PBWT auxiliary data structures.
 */
class ref_haplotype_set {
public:
	/** Minor allele frequency threshold for sparse representation */
	float sparse_maf;

	//--- Counts ---
	/** Total number of variant sites (sparse + bitmatrix) */
	unsigned int n_tot_sites;

	/** Number of rare variant sites stored in sparse matrix */
	unsigned int n_rar_sites;

	/** Number of common variant sites stored in plain matrix */
	unsigned int n_com_sites;

	/** Number of high-quality common variant sites */
	unsigned int n_com_sites_hq;

	/** Number of reference haplotypes (samples) */
	unsigned int n_ref_haps;

	//--- Haplotype data representations ---

	/** Flags indicating if variant is common (true) or rare (false) */
	std::vector<bool> flag_common;

	/** Flags indicating major alleles for variants */
	std::vector<bool> major_alleles;

	/** Reverse mapping from common variant index to total variant site index */
	std::vector<int> common2tot;

	/** Rare alleles per haplotype (sparse representation) */
	std::vector<std::vector<int>> ShapRef;

	/** Rare variant sites per haplotype (indices) */
	std::vector<std::vector<int>> SvarRef;

	/** Bitmatrix storing haplotypes with variants as rows (common variants) */
	bitmatrix HvarRef;

	/** Compressed bitvector for packed representation */
	std::vector<unsigned char> Ypacked;

	/** Indices for small auxiliary data arrays */
	std::vector<std::vector<int>> A_small_idx;

	//--- PBWT auxiliary arrays (not stored permanently) ---

	/** PBWT array A for full data */
	std::vector<int> pbwt_array_A;

	/** PBWT array B for full data */
	std::vector<int> pbwt_array_B;

	/** PBWT array A for small data */
	std::vector<int> pbwt_small_A;

	/** PBWT array B for small data */
	std::vector<int> pbwt_small_B;

	//--- Constructor/Destructor ---

	/** Default constructor */
	ref_haplotype_set();

	/** Virtual destructor */
	virtual ~ref_haplotype_set();

	//--- Member functions ---

	/** Allocate data structures for haplotype set */
	void allocate();

	/**
	 * @brief Build sparse PBWT representation for rare variants from a variant map.
	 * @param M Variant map with variant data
	 */
	void build_sparsePBWT(const variant_map & M);

	/**
	 * @brief Update full PBWT auxiliary arrays for a reference haplotype segment.
	 * @param ref_rac_l Reference allele count or other metric
	 * @param l Index of current locus
	 * @param pbwt_ref_idx Vector of PBWT indices to update
	 */
	void update_full_pbwt_ay(const int ref_rac_l, const int l, std::vector<int>& pbwt_ref_idx);

	/**
	 * @brief Update PBWT arrays for small rare haplotypes.
	 * @param ref_rac_l Reference allele count or other metric
	 * @param rare_small_haps Vector of boolean flags indicating small rare haplotypes
	 */
	void update_small_pbwt_ay(const int ref_rac_l, std::vector<bool> rare_small_haps);

	/**
	 * @brief Initialize small rare haplotype data structures.
	 * @param M Variant map
	 * @param k Index parameter
	 * @param l Index parameter
	 * @param pbwt_ref_idx PBWT reference indices
	 * @param ref_small_hap Vector indicating small haplotypes
	 * @param map_big_small Mapping of big to small haplotypes
	 */
	void init_small_rare(const variant_map & M, const int k, const int l, const std::vector<int>& pbwt_ref_idx, std::vector<bool>& ref_small_hap, std::vector<int>& map_big_small);

	/**
	 * @brief Build initial common variant data structures for a given locus.
	 * @param l Index of locus
	 */
	void build_init_common(const int l);

	//--- Serialization support ---

	friend class boost::serialization::access;

	/**
	 * @brief Serialize the ref_haplotype_set object using Boost Serialization.
	 * @tparam Archive Archive type
	 * @param ar Archive instance
	 * @param version Serialization version
	 */
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & n_tot_sites;
		ar & n_rar_sites;
		ar & n_com_sites;
		ar & n_com_sites_hq;
		ar & n_ref_haps;
		ar & flag_common;
		ar & major_alleles;
		ar & common2tot;
		ar & ShapRef;
		ar & HvarRef;
		ar & Ypacked;
		ar & A_small_idx;
	}

	/**
	 * @brief Update checksum with internal data to support data integrity checks.
	 * @param crc Checksum instance to update
	 */
	void update_checksum(checksum &crc) const
	{
		crc.process_data(n_tot_sites);
		crc.process_data(n_rar_sites);
		crc.process_data(n_com_sites);
		crc.process_data(n_com_sites_hq);
		crc.process_data(n_ref_haps);
		crc.process_data(flag_common);
		crc.process_data(major_alleles);
		crc.process_data(common2tot);
		crc.process_data(ShapRef);
		HvarRef.update_checksum(crc);
		crc.process_data(Ypacked);
		crc.process_data(A_small_idx);
	}
};

#endif
