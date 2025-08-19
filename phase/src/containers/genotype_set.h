/*******************************************************************************
 * @file genotype_set.h
 * @brief Defines the genotype_set class and related structures for storing and
 *        managing genotype data and coverage statistics.
 *
 * This header contains:
 * - **stats_cov** — Structure for coverage and depth statistics
 * - **genotype_set** — Class for storing genotypes across variants and individuals
 *
 * The `genotype_set` class provides methods for serialization, checksum updates,
 * and management of genotype data for phasing, imputation, or downstream analysis.
 *
 * @note
 * - All variables retain their original names from the source implementation.
 * - Designed to work with Boost.Serialization for checkpointing.
 *
 * @see genotype
 * @see variant_map
 ******************************************************************************/

#ifndef _GENOTYPE_SET_H
#define _GENOTYPE_SET_H

#include "otools.h"
#include "checksum_utils.h"
#include "genotype.h"
#include "variant_map.h"

/**
 * @struct stats_cov
 * @brief Stores coverage statistics per individual and depth count distribution.
 *
 * ### Members
 * - `cov_ind` — Per-individual coverage statistics (`stats1D` objects)
 * - `depth_count` — 2D matrix storing depth counts per variant and individual
 *
 * ### Serialization
 * This struct is serializable using Boost.Serialization and can be included in
 * checkpoint files.
 */
struct stats_cov {
	std::vector<stats1D> cov_ind;               /**< Per-individual coverage stats */
	std::vector<std::vector<int>> depth_count;  /**< Depth count matrix */

	friend class boost::serialization::access;

	/**
	 * @brief Serialize the coverage data for checkpointing.
	 * @tparam Archive Boost.Serialization archive type
	 * @param ar Archive reference
	 * @param version Serialization version
	 */
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & cov_ind;
		ar & depth_count;
	}

	/**
	 * @brief Update checksum with coverage data.
	 * @param crc Reference to checksum object
	 */
	void update_checksum(checksum &crc) const
	{
		for (stats1D cov : cov_ind) {
			cov.update_checksum(crc);
		}
		crc.process_data(depth_count);
	}
};

/**
 * @class genotype_set
 * @brief Container for genotype data across multiple variants and individuals.
 *
 * The `genotype_set` class manages a collection of `genotype` objects, tracks
 * counts of variants and individuals, and provides methods for checkpointing
 * and checksum verification.
 *
 * ### Data Members
 * - **n_site** — Number of variants
 * - **n_ind** — Number of individuals
 * - **n_hap** — Number of haplotypes
 * - **vecG** — Vector of pointers to genotype objects
 * - **stats** — Coverage statistics
 *
 * ### Key Features
 * - Boost.Serialization support for checkpointing
 * - Data integrity verification via checksums
 * - Access to genotype objects for downstream analysis
 *
 * @see genotype
 * @see stats_cov
 */
class genotype_set {
public:
	// DATA
	int n_site;                      /**< Number of variants */
	int n_ind;                       /**< Number of individuals */
	int n_hap;                       /**< Number of haplotypes */

	std::vector<genotype*> vecG;     /**< Vector of genotype pointers */
	stats_cov stats;                 /**< Coverage statistics */

	// CONSTRUCTOR / DESTRUCTOR
	genotype_set();
	~genotype_set();

	/**
	 * @brief Serialize genotype set data for checkpointing.
	 *
	 * Ensures that the vector size matches expectations; if mismatched,
	 * raises an error indicating possible data corruption.
	 *
	 * @tparam Archive Boost.Serialization archive type
	 * @param ar Archive reference
	 */
	template<class Archive>
	void serialize_checkpoint_data(Archive &ar)
	{
		size_t vec_size = vecG.size();
		ar & vec_size;
		if (Archive::is_loading::value) {
			if (vec_size != vecG.size()) {
				// Should be caught earlier by checksum verification
				std::stringstream err_str;
				err_str << "Checkpoint file attempting to load with vecG size mismatch. Possible data corruption?";
				vrb.error(err_str.str());
			}
		}
		for (int i = 0; i < vec_size; i++) {
			vecG[i]->serialize_checkpoint_data(ar);
		}
	}

	/**
	 * @brief Update checksum with genotype set data.
	 * @param crc Reference to checksum object
	 */
	void update_checksum(checksum &crc) const
	{
		crc.process_data(n_site);
		crc.process_data(n_ind);
		crc.process_data(n_hap);
		for (genotype* G : vecG) {
			G->update_checksum(crc);
		}
		stats.update_checksum(crc);
	}
};

#endif
