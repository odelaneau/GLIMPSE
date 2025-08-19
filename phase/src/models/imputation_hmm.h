/**
 * @file haplotype_hmm.h
 * @brief Declaration of the imputation_hmm class for genotype imputation using a Hidden Markov Model.
 *
 * This header defines the imputation_hmm class which implements an HMM to perform
 * genotype imputation from haplotype data using forward-backward algorithms.
 * The class relies on aligned memory vectors for efficient SIMD operations and
 * interacts closely with the conditioning_set class representing the conditioning haplotypes.
 *
 */

#ifndef _HAPLOTYPE_HMM_H
#define _HAPLOTYPE_HMM_H

#include "otools.h"
#include "containers/conditioning_set.h"
#include <immintrin.h>
#include <boost/align/aligned_allocator.hpp>

template <typename T>
using aligned_vector32 = std::vector<T, boost::alignment::aligned_allocator < T, 32 > >;

/**
 * @class imputation_hmm
 * @brief Implements a Hidden Markov Model (HMM) for genotype imputation using haplotype data.
 *
 * This class performs genotype imputation leveraging an HMM framework conditioned on
 * a reference haplotype set. It uses SIMD vectorization (AVX) to accelerate forward-backward
 * algorithms computing posterior genotype probabilities.
 *
 * @note This implementation assumes diploid genotypes with biallelic sites.
 *
 * @tparam T Template parameter for aligned vector type (usually float).
 */
class imputation_hmm {
private:
	conditioning_set * C; ///< Pointer to the conditioning set providing haplotype data and model parameters
	unsigned int modK; ///< Number of states rounded up to multiple of 8 for SIMD alignment

	// Dynamic arrays aligned to 32 bytes for SIMD operations:
	aligned_vector32 < float > Emissions; ///< Emission probabilities for each site and allele (size: 2 * number_of_sites)
	aligned_vector32 < float > Alpha;     ///< Forward probabilities matrix
	aligned_vector32 < float > AlphaSum;  ///< Sum of forward probabilities per polymorphic site (for normalization)
	aligned_vector32 < float > Beta;      ///< Backward probabilities vector

public:
	//CONSTRUCTOR/DESTRUCTOR
	/**
	 * @brief Constructor for the imputation_hmm class.
	 * 
	 * Initializes the imputation Hidden Markov Model (HMM) with a given
	 * conditioning set pointer. Prepares the emissions vector sized to
	 * twice the total number of sites in the conditioning set.
	 * 
	 * @param _C Pointer to the conditioning_set object used for imputation.
	 */
	imputation_hmm(conditioning_set *);

	/**
	 * @brief Destructor for the imputation_hmm class.
	 * 
	 * Clears internal data vectors used in the HMM, including Alpha,
	 * AlphaSum, and Emissions, to free resources and reset state.
	 */
	~imputation_hmm();

	/**
	 * @brief Resize internal containers to match the number of states and polymorphic sites.
	 * 
	 * This function calculates `modK` as the smallest multiple of 8
	 * that is greater or equal to the total number of states (`C->n_states`).
	 * It then resizes and initializes the `AlphaSum` vector to the size of
	 * polymorphic sites with zeros, the `Alpha` vector to the product of
	 * polymorphic sites and `modK` with zeros, and the `Beta` vector to size `modK`.
	 */
	void resize();

	/**
	 * @brief Initialize emission probabilities for the HMM using genotype likelihoods and error rates.
	 * 
	 * This method computes the emission probabilities for each site in the conditioning set.
	 * It uses the provided genotype likelihoods (HL) combined with error probabilities to
	 * model the probability of observing the data given the true underlying haplotype state.
	 * 
	 * The error model assumes two error probabilities:
	 * - ee_imp: Probability of correctly observing the allele (no error).
	 * - ed_imp: Probability of observing an erroneous allele.
	 * 
	 * For each site \f$l\f$, the emissions for two allelic states are computed as:
	 * \f[
	 * \begin{aligned}
	 * p_0 &= HL[2l] \times ee\_imp + HL[2l + 1] \times ed\_imp \\
	 * p_1 &= HL[2l] \times ed\_imp + HL[2l + 1] \times ee\_imp
	 * \end{aligned}
	 * \f]
	 * 
	 * These are then normalized to sum to 1 and stored in the Emissions vector as:
	 * \f[
	 * \begin{aligned}
	 * Emissions[2l] &= \frac{p_0}{p_0 + p_1} \\
	 * Emissions[2l + 1] &= \frac{p_1}{p_0 + p_1}
	 * \end{aligned}
	 * \f]
	 * 
	 * This reflects the probability of observing each allelic state given the genotype likelihoods
	 * and the sequencing/imputation error model.
	 * 
	 * @param HL Vector of genotype likelihoods with size 2 * number_of_sites.
	 *           Each pair (HL[2*l], HL[2*l+1]) corresponds to the likelihoods of two allelic states at site l.
	 * 
	 * @note This method assumes that HL values are proportional to likelihoods and uses a simple error
	 *       model to compute emission probabilities suitable for the HMM.
	 * 
	 * @see Durbin R. (2014) "Efficient haplotype matching and storage using the positional Burrows-Wheeler transform (PBWT)" 
	 *      Bioinformatics, 30(9), 1266–1272.
	 *      https://doi.org/10.1093/bioinformatics/btu014
	 * 
	 * @see Li, N., Stephens, M. (2003) "Modeling linkage disequilibrium and identifying recombination hotspots using SNP data."
	 *      Genetics, 165(4), 2213-2233.
	 */
	void init(const std::vector < float > &);

	/**
	 * @brief Runs the forward algorithm for the HMM to compute forward probabilities (Alpha).
	 * 
	 * This method computes the forward probabilities over polymorphic sites using an optimized
	 * vectorized implementation with AVX2 intrinsics (__m256) for performance.
	 * 
	 * The forward probabilities (\f$\alpha_l(k)\f$) represent the probability of observing the data up to
	 * site \f$l\f$ given the HMM state \f$k\f$ at site \f$l\f$.
	 * 
	 * The algorithm handles two cases per site:
	 * - When the site is flagged as 'flat' or 'low quality' (from `flat` or `C->lq_flag`),
	 *   a simplified transition model is applied, and emissions are ignored.
	 * - Otherwise, emission probabilities (from `Emissions`) are incorporated.
	 * 
	 * The HMM transition model uses two factors:
	 * - \f$fact1 = \frac{t[l-1]}{nStates}\f$, representing the probability of switching states.
	 * - \f$fact2 = \frac{nt[l-1]}{\mathrm{AlphaSum}[l-1]}\f$, representing the probability of staying in the same state, normalized by the sum of previous alpha probabilities.
	 *   normalized by the sum of previous alpha probabilities.
	 * 
	 * Vectorization is done over states in multiples of 8 to utilize AVX2 256-bit registers:
	 * - State probabilities are loaded and updated in blocks of 8 floats.
	 * - Emission probabilities are blended using bit masks to select the correct allele emission.
	 * 
	 * @param flat Vector<bool> indicating whether sites are 'flat' (true) and should use simplified model.
	 * 
	 * @note 
	 * - `Alpha` is a float vector storing forward probabilities for each state and polymorphic site.
	 * - `AlphaSum` stores the sum of forward probabilities per site, used for normalization.
	 * - `modK` is the padded number of states to a multiple of 8 for AVX2 vectorization.
	 * - The bitwise allele state information is accessed through `C->Hvar`.
	 * 
	 * @see Rabiner, L. R. (1989). "A tutorial on hidden Markov models and selected applications in speech recognition."
	 *      Proceedings of the IEEE, 77(2), 257-286.
	 * 
	 * @see Durbin R. (2014) "Efficient haplotype matching and storage using the positional Burrows-Wheeler transform (PBWT)" 
	 *      Bioinformatics, 30(9), 1266–1272.
	 */
	void forward(std::vector < bool > &);

	/**
	 * @brief Runs the backward algorithm for the HMM to compute backward probabilities (Beta) and posterior genotype probabilities (HP).
	 * 
	 * This method computes the backward probabilities (\f$\beta_l(k)\f$) for the hidden Markov model over polymorphic sites
	 * using AVX2 SIMD intrinsics for vectorized computation to improve performance.
	 * 
	 * The backward probabilities represent the probability of observing the data from site \f$l+1\f$ to the end given state \f$k\f$ at site \f$l\f$.
	 * 
	 * After computing Beta, the function combines forward (Alpha) and backward (Beta) probabilities with emission and transition probabilities
	 * to estimate posterior genotype probabilities for each site, stored in `HP`.
	 * 
	 * Key aspects of the implementation:
	 * - Handles "flat" or low-quality sites by simplifying the transition and emission model.
	 * - Uses vectorized operations over blocks of 8 states to leverage AVX2 instructions.
	 * - Blends emission probabilities using bitmasks derived from haplotype variant information.
	 * - Normalizes posterior genotype probabilities per site.
	 * - Processes monomorphic sites separately after polymorphic sites.
	 * 
	 * Mathematical notation:
	 * - \f$\beta_l(k)\f$: backward probability of state \f$k\f$ at site \f$l\f$.
	 * - \f$\alpha_l(k)\f$: forward probability at site \f$l\f$.
	 * - \f$t_l, nt_l\f$: transition parameters from conditioning set C.
	 * - \f$e_{imp}, d_{imp}\f$: emission error parameters (\c C->ee_imp, \c C->ed_imp).
	 * - \f$HL\f$: haplotype likelihoods input vector.
	 * - \f$HP\f$: posterior probabilities output vector.
	 * 
	 * @param HL Vector of haplotype likelihoods per allele per site.
	 * @param flat Vector indicating whether sites are "flat" (simplified model) or not.
	 * @param HP Output vector where posterior genotype probabilities are stored.
	 * 
	 * @note 
	 * - `Beta` is resized to hold backward probabilities per state.
	 * - Uses horizontal sums for vector registers to accumulate probabilities.
	 * - SIMD intrinsics (_mm256_*) are used to parallelize arithmetic over 8 floats at a time.
	 * - Bitwise allele masks extracted via `C->Hvar.getByte` and shifted for blending.
	 * 
	 * @see Rabiner, L. R. (1989). "A tutorial on hidden Markov models and selected applications in speech recognition." Proceedings of the IEEE, 77(2), 257-286.
	 * @see Durbin R. (2014). "Efficient haplotype matching and storage using the positional Burrows-Wheeler transform (PBWT)." Bioinformatics, 30(9), 1266–1272.
	 */
	void backward(const std::vector < float > &, std::vector < bool > &, std::vector < float > &);

	/**
	 * @brief Computes posterior genotype probabilities for all sites using the HMM forward-backward algorithm.
	 * 
	 * This function orchestrates the full HMM inference procedure by:
	 *  - Resizing internal data structures to match the number of states and sites.
	 *  - Initializing emission probabilities from haplotype likelihoods (\p HL).
	 *  - Running the forward algorithm to compute forward probabilities.
	 *  - Running the backward algorithm to compute backward probabilities and finally posterior probabilities.
	 * 
	 * The posterior genotype probabilities are stored in \p HP, which is updated in-place.
	 * 
	 * @param HL Vector of haplotype likelihoods per allele per site, used to initialize emission probabilities.
	 * @param flat Boolean vector indicating sites where simplified modeling (flat emission) is used.
	 * @param HP Output vector for storing posterior genotype probabilities per allele per site.
	 * 
	 * @note This method must be called after the conditioning set (C) has been properly initialized.
	 */
	void computePosteriors(const std::vector < float > &, std::vector < bool > &, std::vector < float > &);

};

#endif
