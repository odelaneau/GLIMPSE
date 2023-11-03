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

#include <utils/otools.h>
#include <utils/checksum_utils.h>

#include <containers/variant_map.h>
#include <containers/bitmatrix.h>
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/vector.hpp"

static int p3decode[128] ;
#define ENCODE_MAX1 64		     /* ~64 */
#define ENCODE_MAX2 ((95-63) << 6)   /* ~1k - is this 32 or 31?*/
#define ENCODE_MAX3 ((127-96) << 11) /* ~64k - ditto */

// Code by Richard Durbin [https://github.com/richarddurbin/pbwt]
// doi: 10.1093/bioinformatics/btu014

/* pack3 is a three level run length encoding: n times value
   yp & 0x80 = value
   yp & 0x40 == 0 implies n = yp & 0x3f
   yp & 0x40 == 1 implies
     yp & 0x20 == 0 implies n = (yp & 0x1f) << 6
     yp & 0x20 == 1 implies n = (yp & 0x1f) << 11
   This allows coding runs of length up to 64 * 32 * 32 = 64k in 3 bytes.
   Potential factor of ~1000 over a bit array.
   Build a lookup to avoid conditional operation in uncompression.
*/

static void pack3init (void)
{
  int n ;
  for (n = 0 ; n < 64 ; ++n) p3decode[n] = n ;
  for (n = 64 ; n < 96 ; ++n) p3decode[n] = (n-64) << 6 ;
  for (n = 96 ; n < 128 ; ++n) p3decode[n] = (n-96) << 11 ;
}

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


class ref_haplotype_set {
public:
	float sparse_maf;

	//COUNTS
	unsigned int n_tot_sites;					// #variants [sparse + bitmatrix]
	unsigned int n_rar_sites;					// #variants [rare / sparse matrix]
	unsigned int n_com_sites;					// #variants [common / plain matrix]
	unsigned int n_com_sites_hq;					// #variants [common / plain matrix]
	unsigned int n_ref_haps;					// #haplotypes [reference samples]

	//HAPLOTYPE DATA [plain/sparse bitmatrix representations]
	std::vector < bool > flag_common;				// Is the variant common ? / stored as plain and not sparse?
	std::vector < bool > major_alleles;				// Is the variant common ? / stored as plain and not sparse?
	std::vector < int > common2tot;						//reverse mapping common idx -> i_site
	std::vector < std::vector < int > > ShapRef;				// Rare alleles per haplotype
	std::vector < std::vector <int> > SvarRef;
	bitmatrix HvarRef;								// Bitmatrix of haplotypes, variant first / transpose of Hhap;

	std::vector < unsigned char > Ypacked;
	std::vector < std::vector <int > > A_small_idx;

	//we declare pbwt arrays but we don't store them
	std::vector < int > pbwt_array_A;
	std::vector < int > pbwt_array_B;
	std::vector<int> pbwt_small_A;
	std::vector<int> pbwt_small_B;

	//CONSTRUCTOR/DESTRUCTOR/INITIALIZATION
	ref_haplotype_set();
	virtual ~ref_haplotype_set();

	void allocate();
	void build_sparsePBWT(const variant_map & M);
	void update_full_pbwt_ay(const int ref_rac_l, const int l, std::vector<int>& pbwt_ref_idx);
	void update_small_pbwt_ay(const int ref_rac_l, std::vector<bool> rare_small_haps);
	void init_small_rare(const variant_map & M, const int k, const int l, const std::vector<int>& pbwt_ref_idx, std::vector<bool>& ref_small_hap, std::vector<int>& map_big_small);
	void build_init_common(const int l);


	friend class boost::serialization::access;
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

	void update_checksum(checksum &crc)
	{
		crc.process_data(n_tot_sites);
		crc.process_data(n_rar_sites);
		crc.process_data(n_com_sites);
		crc.process_data(n_com_sites_hq);
		crc.process_data(n_ref_haps);
		crc.prcoess_data(flag_common);
		crc.process_data(major_alleles);
		crc.process_data(common2tot);
		crc.process_data(ShapRef);
		HvarRef.update_checksum(crc);
		crc.process_data(Ypacked);
		crc.process_data(A_small_idx);
	}
};

#endif
