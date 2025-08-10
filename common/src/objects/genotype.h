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

#ifndef _GENOTYPE_H
#define _GENOTYPE_H

#include "otools.h"
#include "checksum_utils.h"
#include "boost/serialization/serialization.hpp"

/**
 * @def _SET32
 * @brief Set the bit at position i in integer n.
 */
#define _SET32(n,i)	((n) |= 1U << (i))

/**
 * @def _CLR32
 * @brief Clear the bit at position i in integer n.
 */
#define _CLR32(n,i)	((n) &= ~(1U << (i)))

/**
 * @def _GET32
 * @brief Get the bit value at position i in integer n.
 */
#define _GET32(n,i)	(((n) >> (i)) & 1U)

/**
 * @struct inferred_genotype
 * @brief Stores inferred genotype data compactly with genotype probabilities.
 *
 * Contains genotype probabilities for two alleles (gp0, gp1), a haploid status flag (hds),
 * and the genotype index (idx). Provides methods to infer the most likely genotype.
 */
struct inferred_genotype {
	float gp0; ///< Genotype probability for allele 0
	float gp1; ///< Genotype probability for allele 1
	bool hds;  ///< Haploid status flag
	int32_t idx; ///< Genotype index

	/**
	 * @brief Default constructor initializes members to zero/false.
	 */
	inferred_genotype() : idx(0), gp0(0), gp1(0), hds(0)
	{
	}

	/**
	 * @brief Parameterized constructor.
	 * @param _idx Genotype index
	 * @param _gp0 Genotype probability allele 0
	 * @param _gp1 Genotype probability allele 1
	 * @param _hds Haploid status flag
	 */
	inferred_genotype(const int _idx, const float _gp0, const float _gp1, const bool _hds) : idx(_idx), gp0(_gp0), gp1(_gp1), hds(_hds) {
	}

	/**
	 * @brief Less-than operator for sorting by genotype index.
	 * @param g Another inferred_genotype
	 * @return True if this idx is less than g.idx
	 */
	bool operator<(const inferred_genotype & g) const {
		return idx < g.idx;
	}

	/**
	 * @brief Infer genotype as the allele with highest genotype probability.
	 * @return 0, 1, or 2 indicating inferred genotype
	 */
	int infer() const {
		const float gp2 = 1.0f - gp1 - gp0;
		if (gp0 > gp1 && gp0 > gp2) return 0;
		if (gp1 > gp0 && gp1 > gp2) return 1;
		if (gp2 > gp0 && gp2 > gp1) return 2;
		return 0;
	}

	/**
	 * @brief Compute genotype probability for the third allele.
	 * @return Probability clamped between 0 and 1
	 */
	float getGp2() const {
		return std::clamp(1.0f - gp1 - gp0, 0.0f, 1.0f);
	}

	/**
	 * @brief Infer haploid genotype based on which allele has higher probability.
	 * @return True if allele 1 probability is greater than allele 0
	 */
	bool infer_haploid() const {
		return gp1 > gp0;
	}

	/**
	 * @brief Serialization method for inferred_genotype.
	 */
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & gp0;
		ar & gp1;
		ar & hds;
		ar & idx;
	}
};

/**
 * @class genotype
 * @brief Stores genotype data including likelihoods and inferred haplotypes.
 *
 * Manages genotype likelihoods, sparse storage of posterior probabilities,
 * and haplotype data for a single sample.
 */
class genotype {
public:
	const std::string name; ///< Sample or genotype name
	const int index; ///< Index in container
	const int n_variants; ///< Number of variants to iterate over
	const int ploidy; ///< Ploidy (usually 1 or 2)
	const int hapid; ///< Haplotype ID
	int stored_cnt; ///< Count of stored genotypes

	std::vector < unsigned char > GL; ///< Original genotype likelihoods
	std::vector <bool> flat; ///< Flat flags for genotypes

	/// Sparse storage of inferred genotype posteriors and haploid status
	std::vector < inferred_genotype > stored_data;

	/// First haplotype alleles (to be removed since stored in haplotype_set)
	std::vector < bool > H0;

	/// Second haplotype alleles (to be removed since stored in haplotype_set)
	std::vector < bool > H1;

	/**
	 * @brief Constructor
	 * @param _name Sample name
	 * @param _index Index of genotype in container
	 * @param _n_variants Number of variants
	 * @param _ploidy Ploidy (default 2)
	 * @param _hapid Haplotype ID
	 */
	genotype(const std::string _name, const int _index, const int _n_variants, const int _ploidy, const int _hapid);

	/// Destructor
	~genotype();

	/// Allocate memory for genotype data
	void allocate();

	/// Free allocated memory
	void free();

	/**
	 * @brief Initialize haplotype likelihoods with minimum likelihood threshold.
	 * @param likelihoods Vector to initialize
	 * @param min_gl Minimum genotype likelihood
	 */
	void initHaplotypeLikelihoods(std::vector < float > &, const float min_gl);

	/**
	 * @brief Sample first haplotype (H0) from likelihoods.
	 * @param likelihoods Vector of likelihoods
	 */
	void sampleHaplotypeH0(const std::vector < float > &);

	/**
	 * @brief Sample second haplotype (H1) from likelihoods.
	 * @param likelihoods Vector of likelihoods
	 */
	void sampleHaplotypeH1(const std::vector < float > &);

	/**
	 * @brief Construct haplotype likelihoods.
	 * @param likelihoods Vector to fill
	 * @param use_first_hap Flag whether to use first haplotype
	 * @param min_gl Minimum genotype likelihood
	 */
	void makeHaplotypeLikelihoods(std::vector < float > &, bool, const float min_gl) const;

	/**
	 * @brief Store genotype posteriors and haplotypes from likelihoods.
	 * @param likelihoods Vector of likelihoods
	 */
	void storeGenotypePosteriorsAndHaplotypes(const std::vector < float > &);

	/**
	 * @brief Store genotype posteriors and haplotypes from two likelihood vectors.
	 * @param likelihoods1 First vector of likelihoods
	 * @param likelihoods2 Second vector of likelihoods
	 */
	void storeGenotypePosteriorsAndHaplotypes(const std::vector < float > &, const std::vector < float > &);

	/// Sort, normalize and infer genotype calls
	void sortAndNormAndInferGenotype();

	/**
	 * @brief Serialize checkpoint data (stored genotype info).
	 * @param ar Archive
	 */
	template<class Archive>
	void serialize_checkpoint_data(Archive &ar)
	{
		ar & stored_cnt;
		ar & stored_data;
		ar & H0;
		ar & H1;
	}

	/**
	 * @brief Update checksum with genotype data for data integrity.
	 * @param crc Checksum object
	 */
	void update_checksum(checksum &crc) const
	{
		crc.process_data(name);
		crc.process_data(index);
		crc.process_data(n_variants);
		crc.process_data(ploidy);
		crc.process_data(hapid);
		crc.process_data(GL);
		crc.process_data(flat);
	}
};

#endif
