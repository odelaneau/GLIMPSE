/*******************************************************************************
 * @file conditioning_set.h
 * @brief Defines the conditioning_set class for handling variant conditioning states.
 *
 * This class manages variant classification, selection routines, and transition
 * probabilities for phasing and imputation processes. It uses references to a
 * variant map and haplotype set, and organizes variants into categories based
 * on allele frequency and other criteria.
 *
 * ## Variant Types
 * - **TYPE_COMMON (0)** — Common variant
 * - **TYPE_RARE (1)** — Rare variant
 * - **TYPE_MONO (2)** — Monomorphic variant
 *
 * ## Main Responsibilities
 * - Maintain lists of polymorphic and monomorphic sites
 * - Store variant types and major alleles
 * - Handle conditioning state selection and compaction
 * - Compute and update transition probabilities
 *
 * @note
 * - Uses `std::vector` for storing variant-related information
 * - Probabilities are clamped to avoid underflow
 *
 * @see variant_map
 * @see haplotype_set
 ******************************************************************************/

#ifndef _CONDITIONING_SET_H
#define _CONDITIONING_SET_H

#include "otools.h"
#include "variant_map.h"
#include "containers/haplotype_set.h"

// Types included in vector <char> var_type
#define TYPE_COMMON 0  /**< Common variant */
#define TYPE_RARE   1  /**< Rare variant */
#define TYPE_MONO   2  /**< Monomorphic variant */

/**
 * @class conditioning_set
 * @brief Stores and manages conditioning states for phasing and imputation.
 *
 * The `conditioning_set` class groups variants into logical categories and
 * manages the selection of conditioning states for phasing algorithms.
 *
 * ### Static Data
 * - `mapG` — Reference to the variant map
 * - `H` — Reference to the haplotype set
 *
 * ### Counts
 * - `n_ref_haps` — Number of reference haplotypes
 * - `n_eff_haps` — Effective number of haplotypes
 * - `n_com_sites` — Number of common sites
 * - `n_tot_sites` — Total number of sites
 * - `n_states` — Number of states
 *
 * ### Variant Information
 * - `var_type` — Variant type: common, rare, or monomorphic
 * - `major_alleles` — Major allele flags for direct imputation
 * - `polymorphic_sites` — List of common sites (TYPE_COMMON, TYPE_RARE)
 * - `monomorphic_sites` — List of monomorphic sites (TYPE_MONO)
 * - `lq_flag` — Low-quality variant flags
 *
 * ### Conditioning States
 * - `idxHaps_ref` — Indexes of conditioning states in reference haplotype set
 * - `Svar` — Sparse bitmatrix (variant-first)
 * - `Hvar` — Plain bitmatrix (variant-first)
 *
 * ### Probabilities
 * - `t`, `nt` — Transition probabilities
 * - `ed_phs`, `ee_phs` — Error/detection probabilities for phasing
 * - `ed_imp`, `ee_imp` — Error/detection probabilities for imputation
 *
 * ### Selection Parameters
 * - `Kinit` — Initial number of states
 * - `Kpbwt` — PBWT parameter for state search
 *
 * @see compactSelection()
 * @see select()
 * @see updateTransitions()
 */
class conditioning_set {
public:
	// STATIC DATA
	const variant_map & mapG;       /**< Reference to the variant map */
	const haplotype_set & H;        /**< Reference to the haplotype set */

	// COUNTS
	unsigned int n_ref_haps;        /**< Number of reference haplotypes */
	unsigned int n_eff_haps;        /**< Effective haplotype count */
	unsigned int n_com_sites;       /**< Number of common sites */
	unsigned int n_tot_sites;       /**< Total number of sites */
	unsigned int n_states;          /**< Number of states */

	// CONST
	const double nrho;              /**< Recombination rate constant */
	const double one_l;             /**< Precomputed constant for limits */

	// VARIANT TYPES
	std::vector<unsigned char> var_type;       /**< Type of variant: TYPE_COMMON / TYPE_RARE / TYPE_MONO */
	std::vector<bool> major_alleles;           /**< Major alleles for direct imputation */
	std::vector<unsigned int> polymorphic_sites; /**< List of common sites */
	std::vector<unsigned int> monomorphic_sites; /**< List of monomorphic sites */
	std::vector<bool> lq_flag;                  /**< Low-quality variant flags */

	// CONDITIONING STATES
	std::vector<unsigned int> idxHaps_ref;               /**< Indexes of conditioning states in reference haplotype set */
	std::vector<std::vector<unsigned int>> Svar;         /**< Sparse bitmatrix: variant-first */
	bitmatrix Hvar;                                      /**< Plain bitmatrix: variant-first */

	// TRANSITION & EMISSION PROBABILITIES
	std::vector<float> t;            /**< Transition probabilities */
	std::vector<float> nt;           /**< Non-transition probabilities */
	const float ed_phs;              /**< Error probability for phasing */
	const float ee_phs;              /**< Emission probability for phasing */
	const float ed_imp;              /**< Error probability for imputation */
	const float ee_imp;              /**< Emission probability for imputation */

	int Kinit;                       /**< Initial number of states */
	int Kpbwt;                       /**< PBWT search depth parameter */

	std::vector<unsigned int> swap_ref; /**< Reference haplotype swap buffer */
	std::vector<unsigned int> swap_tar; /**< Target haplotype swap buffer */

	// CONSTRUCTOR / DESTRUCTOR / INITIALIZATION
	conditioning_set(const variant_map &, const haplotype_set &, const unsigned int, const unsigned int, const int, const int, const float, const float, const bool);
	~conditioning_set();
	void clear();

	// SELECTION ROUTINES
	void compactSelection(const int ind, const int iter);
	void select(const int ind, const int iter);

	// UPDATE TRANSITION PROBS
	void updateTransitions();

	/**
	 * @brief Compute the transition probability between two haplotype positions.
	 * @param prev_abs_idx Index of previous haplotype position
	 * @param next_abs_idx Index of next haplotype position
	 * @return Transition probability value (clamped between 1e-7 and one_l)
	 */
	inline
	float getTransition(const int prev_abs_idx, const int next_abs_idx) const {
		return std::clamp(
			-expm1(nrho * (mapG.vec_pos[next_abs_idx]->cm - mapG.vec_pos[prev_abs_idx]->cm)),
			1e-7,
			one_l
		);
	}
};

#endif
