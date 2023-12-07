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

#ifndef _GENOTYPE_SET_H
#define _GENOTYPE_SET_H

#include <utils/otools.h>
#include <utils/checksum_utils.h>

#include <objects/genotype.h>
#include <containers/variant_map.h>

struct stats_cov {
	std::vector<stats1D> cov_ind;
	std::vector<std::vector<int>> depth_count;

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & cov_ind;
		ar & depth_count;
	}

	void update_checksum(checksum &crc)
	{
		for (const stat1D cov : cov_ind) {
			cov.update_checksum(crc);
		}
		crc.process_data(depth_count);
	}
};

class genotype_set {
public:
	//DATA
	int n_site, n_ind;					//Number of variants, number of individuals
	int n_hap;

	std::vector < genotype * > vecG;			//Vector of genotypes
	stats_cov stats;

	//CONSTRUCTOR/DESTRUCTOR
	genotype_set();
	~genotype_set();

	template<class Archive>
	void serialize_original_data_to_archive(Archive &ar) const
	{
		ar << n_site;
		ar << n_ind;
		ar << n_hap;
		const size_t vec_size = vecG.size();
		ar << vec_size;
		for (int i=0; i<vec_size; i++) {
			vecG[i]->serialize_original_data_to_archive(ar);
		}
		ar & stats;
	}

	template<class Archive>
	void serialize_checkpoint_data(Archive &ar)
	{
		size_t vec_size = vecG.size();
		ar & vec_size;
		for (int i=0; i<vec_size; i++) {
			vecG[i]->serialize_checkpoint_data(ar);
		}
	}

	void update_checksum(checksum &crc)
	{
		crc.process_data(n_site);
		crc.process_data(n_ind);
		crc.process_data(n_hap);
		crc.process_data(vecG);
		stats.update_checksum(crc);
	}
};

#endif
