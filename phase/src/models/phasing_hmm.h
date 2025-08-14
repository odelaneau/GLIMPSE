/*******************************************************************************
 * @file phasing_hmm.h
 * @brief Hidden Markov Model (HMM) implementation for haplotype phasing with low-pass sequencing data.
 *
 * This class implements a diploid haplotype phasing HMM leveraging a conditioning set,
 * vectorized SIMD instructions (AVX), and advanced probabilistic modeling to phase
 * common and rare variants from genotype data.
 *
 * @details
 * The model integrates observed genotypes with a conditioning reference panel and applies
 * a state-space HMM with:
 * - States representing haplotype configurations
 * - Transitions informed by recombination rates between loci
 * - Emission probabilities incorporating sequencing error and genotype likelihoods
 * 
 * It handles different variant classes:
 * - Peak heterozygous variants (common, confidently called)
 * - Peak homozygous variants
 * - Flat heterozygous variants (rare or uncertain)
 * 
 * The implementation utilizes:
 * - Aligned SIMD vector operations for performance (_mm256_ intrinsics)
 * - Multi-pass forward-backward algorithms for posterior inference
 * - Probabilistic diplotype sampling and imputation for missing/rare genotypes
 *
 * @note
 * - Assumes AVX2 instruction set availability.
 * - Constants and helper macros (e.g., HAP_NUMBER, ALLELE) are defined for clarity.
 ******************************************************************************/

#ifndef _DIPLOTYPE_HMM_H
#define _DIPLOTYPE_HMM_H

#include "otools.h"
#include "containers/conditioning_set.h"
#include <immintrin.h>
#include <boost/align/aligned_allocator.hpp>

/**
 * @brief Vector type with 32-byte alignment for AVX operations.
 * @tparam T Element type.
 */
template <typename T>
using aligned_vector32 = std::vector<T, boost::alignment::aligned_allocator < T, 32 > >;

/// Number of haplotypes modeled per individual (fixed at 8).
#define HAP_NUMBER 8

/// Variant type codes for phasing model
#define VAR_PEAK_HET 0   ///< Common heterozygous variant (peak).
#define VAR_PEAK_HOM -1  ///< Common homozygous variant (peak).
#define VAR_FLAT_HET -2  ///< Rare or uncertain heterozygous variant (flat).

/// Macro to test presence of allele at haplotype `hap` and position `pos`.
#define ALLELE(hap, pos) ((hap) & (1<<(pos)))

/**
 * @brief Compute horizontal sum of 8 packed floats in __m256.
 * @param a __m256 register containing 8 floats.
 * @return Sum of all 8 floats.
 */
inline
float horizontal_add (const __m256& a) //@simo: implemented horizontal add with instructions in the 4B minimum
{
    __m128 vlow = _mm256_castps256_ps128(a);
    __m128 vhigh = _mm256_extractf128_ps(a, 1); // high 128
   vlow = _mm_add_ps(vlow, vhigh);     // add the low 128
   __m128 shuf = _mm_movehdup_ps(vlow);        // broadcast elements 3,1 to 2,0
   __m128 sums = _mm_add_ps(vlow, shuf);
   shuf = _mm_movehl_ps(shuf, sums); // high half -> low half
   sums = _mm_add_ss(sums, shuf);
   return _mm_cvtss_f32(sums);
}

/**
 * @class phasing_hmm
 * @brief Implements a diploid HMM for haplotype phasing with AVX vectorization.
 *
 * This class supports phasing common and rare variants using a conditioning
 * set and incorporates emission and transition modeling informed by sequencing
 * error rates and genotype probabilities.
 *
 * @note
 * - Uses aligned vectors for SIMD-friendly memory layout.
 * - Includes data structures for forward-backward passes, diplotype sampling,
 *   and imputation of rare variants.
 */
class phasing_hmm {
private:

	conditioning_set * C;                ///< Pointer to conditioning set with reference haplotypes and parameters.

    // Variant classification and indexing arrays
    std::vector<signed char> VAR_TYP;   ///< Variant type per locus: peak het, peak hom, or flat het.
    std::vector<bool> VAR_ALT;           ///< Alternative allele indicator for variants.
    std::vector<int> VAR_ABS;             ///< Absolute genomic indices of variants.
    std::vector<int> VAR_REL;             ///< Relative indices for variants within the conditioning set.

    // Segmentation and locus counters
    unsigned int n_segs;                 ///< Number of segments for phasing.
    unsigned int n_miss;                 ///< Number of missing or rare variant loci.

    std::vector<int> segments;           ///< Segment sizes for HMM blocks.

    // Current processing indices
    int curr_idx_locus;                  ///< Current variant locus index.
    int curr_abs_locus;                  ///< Absolute locus coordinate.
    int curr_rel_locus;                  ///< Relative locus index within segment.
    int curr_segment_index;              ///< Current segment index.
    int curr_segment_locus;              ///< Current locus within segment.
    int curr_missing_locus;              ///< Current index for missing variant locus.

    // Probabilities and dynamic arrays for forward-backward computations
    float probSumT;                      ///< Total sum of probabilities at current locus.
    aligned_vector32<float> prob;        ///< Probability vector for each haplotype-state pair.
    aligned_vector32<float> probSumK;    ///< Summed probabilities over states.
    aligned_vector32<float> probSumH;    ///< Summed probabilities over haplotypes.

    aligned_vector32<float> phasingProb;      ///< Stores phased probabilities per segment.
    aligned_vector32<float> phasingProbSum;   ///< Summed phasing probabilities over haplotypes.
    std::vector<float> phasingProbSumSum;     ///< Total sums of phasing probabilities.

    aligned_vector32<float> imputeProb;        ///< Imputation probability vectors for rare variants.
    aligned_vector32<float> imputeProbSum;     ///< Summed imputation probabilities.
    aligned_vector32<float> imputeProbSumSum;  ///< Total sums of imputation probabilities.
    aligned_vector32<float> imputeProbOf1s;    ///< Probability of allele=1 for imputed variants.
    std::vector<int> dip_sampled;               ///< Sampled diplotype indices per segment.

    // Static arrays for emission probabilities and diplotype probabilities
    std::vector<float> DProbs;                 ///< Diplotype sampling probabilities.
    std::vector<aligned_vector32<float>> EMIT0; ///< Emission probabilities for allele 0.
    std::vector<aligned_vector32<float>> EMIT1; ///< Emission probabilities for allele 1.
    aligned_vector32<float> HProbs;             ///< Haploid probabilities.
    float sumHProbs;                            ///< Sum of haplotype probabilities.
    float sumDProbs;                            ///< Sum of diplotype probabilities.

    // Transition probabilities for homozygous and heterozygous states
    float nt; ///< Probability of no recombination.
    float yt; ///< Probability of recombination.

	//INLINED ROUTINES
	/**
	 * @brief Initialize peak heterozygous emission probabilities for the current heterozygous site.
	 * 
	 * This function loads emission probabilities for the current heterozygous site into SIMD registers,
	 * then iterates over all HMM states to accumulate and store emission probabilities.
	 * 
	 * Specifically:
	 * - Loads emission probabilities \f$EMIT0[curr\_het]\f$ and \f$EMIT1[curr\_het]\f$ into two 256-bit SIMD registers.
	 * - For each state \f$k\f$ (total \f$C->n\_states\f$ states):
	 *   - Determines the allele state \f$ah = C->Hvar.get(curr\_rel\_locus, k)\f$ (0 or 1).
	 *   - Adds the corresponding emission probability vector \f$emits[ah]\f$ to the cumulative sum \f$_sum\f$.
	 *   - Stores the emission probability vector into the probability buffer \f$prob[i]\f$.
	 * - After the loop, the cumulative emission probabilities vector \f$_sum\f$ is stored into \f$probSumH\f$,
	 *   and its horizontal sum (sum of all floats) is stored in \f$probSumT\f$.
	 * 
	 * @note Uses AVX intrinsics for SIMD operations:
	 * - \c _mm256_load_ps: loads 8 floats into a __m256 register.
	 * - \c _mm256_add_ps: performs element-wise addition on __m256 registers.
	 * - \c _mm256_store_ps: stores __m256 register contents back to memory.
	 * 
	 * @param curr_het Index of the current heterozygous site to process.
	 */
	void INIT_PEAK_HET(int);

	/**
	 * @brief Initialize peak homozygous emission probabilities for the current homozygous site.
	 * 
	 * This function sets emission probabilities based on whether the allele matches the given allele `ag`.
	 * It uses SIMD vector instructions for efficient processing across all HMM states.
	 * 
	 * Specifically:
	 * - Defines two emission probability vectors:
	 *   - For matching allele: constant 1.0 vector \f$(1,1,\dots,1)\f$
	 *   - For non-matching allele: constant \f$\frac{ed\_phs}{ee\_phs}\f$ vector
	 * - For each HMM state \f$k\f$ (total \f$C->n\_states\f$ states):
	 *   - Checks if the allele at state differs from `ag`: \f$ag\_ah = (Hvar.get(curr\_rel\_locus, k) \neq ag)\f$
	 *   - Adds the corresponding emission vector \f$emits[ag\_ah]\f$ to the cumulative sum \f$_sum\f$.
	 *   - Stores this emission vector into the probability buffer at position \f$prob[i]\f$.
	 * - After processing all states:
	 *   - Stores the cumulative sum vector \f$_sum\f$ in \f$probSumH\f$.
	 *   - Computes and stores the horizontal sum (sum of all elements in \f$_sum\f$) in \f$probSumT\f$.
	 * 
	 * @note SIMD intrinsics used:
	 * - \c _mm256_set1_ps: sets all 8 floats in the __m256 register to the same value.
	 * - \c _mm256_add_ps: element-wise addition of __m256 registers.
	 * - \c _mm256_store_ps: stores __m256 register to memory.
	 * 
	 * @param ag The allele (true or false) considered homozygous at the current locus.
	 */
	void INIT_PEAK_HOM(bool);

	/**
	 * @brief Initialize emission probabilities uniformly for heterozygous sites (flat initialization).
	 * 
	 * This method assigns equal probabilities across all haplotype-state combinations,
	 * representing a uniform prior belief before any observation or evidence.
	 * 
	 * Specifically:
	 * - Sets all entries in the \c prob vector to \f$\frac{1}{\text{HAP\_NUMBER} \times C->n\_states}\f$,
	 *   distributing probability evenly across all haplotypes and states.
	 * - Sets all entries in \c probSumH to \f$\frac{1}{\text{HAP\_NUMBER}}\f$,
	 *   representing the sum of probabilities per haplotype.
	 * - Sets \c probSumT to 1.0, representing the total probability mass.
	 * 
	 * This flat initialization is typically used as a non-informative prior for heterozygous sites.
	 */
	void INIT_FLAT_HET();

	/**
	 * @brief Perform one iteration of the HMM forward update for a heterozygous site (peak heterozygous state).
	 * 
	 * This function updates the probability distribution over haplotype states at the current heterozygous locus.
	 * It incorporates transition probabilities and emission probabilities in a vectorized manner using AVX intrinsics.
	 * 
	 * @param curr_het Index of the current heterozygous site to process.
	 * 
	 * Detailed explanation:
	 * - \c _tFreq calculates the scaled transition frequency from the previous step weighted by \c yt,
	 *   normalized by the number of states and total probability sum.
	 * - \c _nt is the scaled non-transition probability weighted by \c nt and normalized by total probability.
	 * - \c emits holds emission probabilities for the two possible allelic states at \c curr_het.
	 * - The loop iterates over all haplotype states:
	 *   - Retrieves the allele state \c ah from the haplotype variant matrix.
	 *   - Computes updated state probability by multiplying previous probability, transition/non-transition factors, and emission.
	 *   - Sums the updated probabilities and stores them.
	 * - Finally, the function updates the per-haplotype probability sums (\c probSumH) and total probability sum (\c probSumT).
	 * 
	 * This step is part of the forward algorithm updating posterior probabilities given current data.
	 */
	void RUN_PEAK_HET(int);

	/**
	 * @brief Perform one iteration of the HMM forward update for a homozygous site (peak homozygous state).
	 * 
	 * This function updates the probability distribution over haplotype states at the current homozygous locus.
	 * It incorporates transition probabilities, emission mismatches, and uses AVX intrinsics for vectorization.
	 * 
	 * @param ag The allele genotype (true or false) for the current homozygous site.
	 * 
	 * Detailed explanation:
	 * - \c _tFreq computes the scaled transition frequency weighted by \c yt,
	 *   normalized by the number of states and total probability sum.
	 * - \c _nt holds the scaled non-transition probability weighted by \c nt and normalized by total probability.
	 * - \c _mism is the emission mismatch factor representing the relative likelihood of observing
	 *   the non-matching allele (ed_phs / ee_phs).
	 * - The loop iterates over all haplotype states:
	 *   - Retrieves the allele state \c ah from the haplotype variant matrix.
	 *   - Calculates updated probability by combining previous probability, transition/non-transition factors.
	 *   - If the allele state does not match the genotype (\c ag != \c ah), the probability is scaled by the mismatch factor.
	 *   - Accumulates the updated probabilities and stores them.
	 * - Finally, updates the per-haplotype sums (\c probSumH) and total probability sum (\c probSumT).
	 * 
	 * This function is a critical step in the forward algorithm for haplotype phasing in homozygous loci.
	 */
	void RUN_PEAK_HOM(bool);

	/**
	 * @brief Perform one iteration of the HMM forward update for a heterozygous site with uniform emission probabilities.
	 * 
	 * This function updates the haplotype state probabilities assuming a "flat" emission model where no specific allele
	 * information is used (i.e., emission probabilities are uniform across states).
	 * 
	 * @details
	 * - \c _tFreq computes the transition frequency weighted by \c yt and normalized by the number of states and total sum.
	 * - \c _nt holds the scaled non-transition probability weighted by \c nt and normalized by total probability sum.
	 * - For each haplotype state, the new probability is calculated as:
	 *   \f[
	 *     \text{prob}_\text{curr} = \text{prob}_\text{prev} \times \text{_nt} + \text{_tFreq}
	 *   \f]
	 * - The updated probabilities are summed and stored back into the state probability vector.
	 * - Finally, per-haplotype sums (\c probSumH) and total probability sum (\c probSumT) are updated.
	 * 
	 * This function is used in the phasing HMM forward algorithm for heterozygous loci without allele-specific emission.
	 */
	void RUN_FLAT_HET();

	/**
	 * @brief Collapse haplotype probabilities for a heterozygous site by distributing collapsed state sums across haplotype states.
	 * 
	 * This function updates the per-haplotype state probabilities (\c prob) from the collapsed state probabilities (\c probSumK),
	 * applying transition probabilities and emission likelihoods at the current heterozygous locus.
	 * 
	 * @param curr_het Index of the current heterozygous site.
	 * 
	 * @details
	 * - \c _tFreq and \c _nt represent transition and non-transition scaled probabilities, normalized by total sums.
	 * - For each haplotype state \f$k\f$, the probability is computed as:
	 *   \f[
	 *     \text{prob}_\text{curr} = (\text{probSumK}[k] \times \text{\_nt} + \text{\_tFreq}) \times \text{emission}(k)
	 *   \f]
	 *   where \c emission(k) is selected from \c EMIT0 or \c EMIT1 based on allele state.
	 * - The updated probabilities are summed into \c probSumH and total sum \c probSumT is updated.
	 * 
	 * This step is part of the phasing HMM algorithm used to refine haplotype probabilities at heterozygous sites.
	 */
	void COLLAPSE_PEAK_HET(int);

	/**
	 * @brief Collapse haplotype probabilities for a homozygous site by distributing collapsed state sums across haplotype states.
	 * 
	 * This function updates the per-haplotype state probabilities (\c prob) from the collapsed state probabilities (\c probSumK),
	 * applying transition probabilities, emission mismatches for homozygous genotypes, and scaling.
	 * 
	 * @param ag The allele genotype at the current locus (true/false indicating allele state).
	 * 
	 * @details
	 * - Computes transition factors \c _tFreq and \c _nt based on global sums \c probSumH and \c probSumT.
	 * - For each haplotype state \f$k\f$, calculates:
	 *   \f[
	 *     \text{prob}_{curr} = \left(\text{probSumK}[k] \times \_nt + \_tFreq \right) \times
	 *     \begin{cases}
	 *       1 & \text{if } ag = ah \\
	 *       \frac{ed\_phs}{ee\_phs} & \text{if } ag \neq ah
	 *     \end{cases}
	 *   \f]
	 *   where \c ah is the allele state of haplotype \f$k\f$ at the current locus.
	 * - The probabilities are accumulated to update \c probSumH and \c probSumT.
	 * 
	 * This function is part of the phasing HMM and handles the case where the observed genotype is homozygous,
	 * incorporating mismatch penalties between expected and observed alleles.
	 */
	void COLLAPSE_PEAK_HOM(bool);

	/**
	 * @brief Collapse haplotype probabilities for a heterozygous site assuming a flat emission model.
	 * 
	 * This function distributes collapsed state probabilities (\c probSumK) uniformly across haplotype states
	 * without considering allele-specific emissions (i.e., assumes equal likelihood for both alleles).
	 * 
	 * @details
	 * - Computes transition factors \c _tFreq and \c _nt based on global sums \c probSumH and \c probSumT.
	 * - For each haplotype state \f$k\f$, updates the per-haplotype probability \c prob as:
	 *   \f[
	 *     \text{prob}_{curr} = \text{probSumK}[k] \times \_nt + \_tFreq
	 *   \f]
	 * - Accumulates these probabilities into \c probSumH and updates the total sum \c probSumT.
	 * 
	 * This method is used in the phasing HMM when emission probabilities are assumed uniform for heterozygous loci,
	 * simplifying calculations by ignoring allele-specific emission differences.
	 */
	void COLLAPSE_FLAT_HET();

	/**
	 * @brief Compute the sum of probabilities over all haplotypes for each state.
	 * 
	 * For each state \f$k\f$, this function sums the probabilities stored in the \c prob vector
	 * corresponding to all haplotypes of that state.
	 * 
	 * The sums are stored in the \c probSumK vector at index \f$k\f$.
	 * 
	 * @note
	 * - Uses AVX intrinsics to efficiently load and horizontally sum 8 floats at a time.
	 * - Assumes that the \c prob vector is organized such that probabilities for each state \f$k\f$
	 *   are stored contiguously in blocks of size \c HAP_NUMBER.
	 */
    void SUMK();

	/**
	 * @brief Perform the haplotype transition update step in the phasing HMM.
	 *
	 * This method updates the haplotype transition probabilities \f$HProbs\f$ for the current locus
	 * by combining previous probabilities and the transition model.
	 *
	 * The transition model uses:
	 * - \f$yt = P(\text{switch state})\f$, obtained from the conditioning set between the current
	 *   and next loci.
	 * - \f$nt = 1 - yt\f$, the probability of not switching states.
	 *
	 * For each haplotype \f$h_1\f$, the new probability is computed as:
	 * \[
	 * HProbs[h_1, :] = \sum_{k=1}^{nStates} \left( \frac{prob[k, h_1] \cdot nt}{probSumT} + \frac{probSumH[h_1] \cdot yt}{probSumT \cdot nStates} \right) \cdot phasingProb[next\_locus, k, :]
	 * \]
	 *
	 * @return bool Returns true if the resulting sum of haplotype probabilities contains NaN,
	 *              infinity, or is smaller than the minimum positive float (indicating an error).
	 */
    bool TRANS_HAP();

	/**
	 * @brief Sample the diplotype (pair of haplotypes) for the current segment based on haplotype probabilities.
	 *
	 * This method calculates the diplotype probabilities \f$DProbs\f$ by combining haplotype
	 * transition probabilities \f$HProbs\f$ for the two haplotypes in the diplotype.
	 *
	 * The diplotype probability for diplotype \f$d\f$ is computed as:
	 * \[
	 * DProbs[d] = \frac{HProbs[h_0, d]}{\sum HProbs} \times \frac{HProbs[h_1, HAP\_NUMBER - d - 1]}{\sum HProbs}
	 * \]
	 * where \f$h_0\f$ and \f$h_1\f$ are the previously sampled haplotypes for the current segment.
	 *
	 * The method samples a new diplotype index from the distribution \f$DProbs\f$ using the
	 * random number generator \c rng and stores it in \c dip_sampled at the next segment index.
	 *
	 * @return bool Returns true if the sum of diplotype probabilities is NaN, infinite,
	 *              or smaller than the minimum positive float, indicating a numerical error.
	 *              Returns false otherwise.
	 */
    bool SAMPLE_DIP();

	/**
	 * @brief Computes the posterior probability of the alternate allele (1) for a heterozygous site under a flat model.
	 *
	 * This function calculates the posterior probabilities by combining the left and right haplotype probabilities,
	 * normalizing and clamping the results between 0 and 1. It uses SIMD intrinsics (__m256) for vectorized operations
	 * on haplotype states to speed up computation.
	 *
	 * @note
	 * - `C->n_states` is the number of conditioning states (haplotype states).
	 * - `HAP_NUMBER` is the number of haplotypes per state (typically 8 for AVX256).
	 * - `curr_missing_locus` is the index of the missing locus to impute.
	 * - `curr_rel_locus` is the index of the relative locus used for haplotype masking.
	 *
	 * @details
	 * For each conditioning state k:
	 * - Load and scale the imputation probability from the right side (_p1).
	 * - Load and scale the haplotype probability from the left side (_p2).
	 * - Multiply the scaled probabilities and accumulate sums separately for allele 0 and allele 1,
	 *   according to the haplotype allele at locus `curr_rel_locus`.
	 *
	 * After the loop, compute the normalized probability of allele 1:
	 * \f[
	 * P(1) = \frac{\text{sums}[1]}{\text{sums}[0] + \text{sums}[1]}
	 * \f]
	 * Clamp the value to \[0,1\] and store the result in `imputeProbOf1s`.
	 *
	 * @warning This function assumes that the vectors `imputeProb`, `prob`, `imputeProbSum`, `probSumH`,
	 * and `imputeProbOf1s` have been correctly allocated and aligned for AVX operations.
	 */
    void IMPUTE_FLAT_HET();

public:
	//CONSTRUCTOR/DESTRUCTOR
	/**
	 * @brief Constructs a phasing Hidden Markov Model (HMM) instance.
	 * 
	 * This constructor initializes the phasing HMM by setting up internal
	 * states and emission probability vectors based on the provided conditioning set.
	 * It prepares emission probability arrays for three genotype contexts, 
	 * allocating aligned memory for SIMD-optimized processing.
	 * 
	 * @param _C Pointer to a conditioning_set object containing haplotype and genotype data.
	 *           This set provides transition and emission parameters for the phasing model.
	 * 
	 * @details
	 * The following member variables are initialized:
	 * - Counters and indices related to loci and segments:
	 *   - `n_segs`, `n_miss`, `curr_idx_locus`, `curr_abs_locus`, `curr_rel_locus`, 
	 *     `curr_segment_index`, `curr_segment_locus`, `curr_missing_locus` are set to zero.
	 * - Probability sums and transition parameters:
	 *   - `probSumT`, `sumHProbs`, `sumDProbs`, `nt`, `yt` are set to 0.0f.
	 * - The conditioning set pointer `C` is assigned `_C`.
	 * 
	 * The emission vectors `EMIT0` and `EMIT1` are initialized as vectors of size 3,
	 * each containing aligned vectors of floats with length `HAP_NUMBER` (usually 8).
	 * These vectors hold emission probabilities for the phasing model, representing
	 * mismatch and match likelihoods for each haplotype state under different genotype assumptions.
	 * 
	 * Emission values are explicitly set as follows (example values based on phasing error rates):
	 * - `C->ee_phs` for emission probabilities of matching alleles.
	 * - `C->ed_phs` for emission probabilities of mismatching alleles.
	 * 
	 * The vector `HProbs` is initialized as an aligned vector of floats with size `HAP_NUMBER * HAP_NUMBER`.
	 * It stores pairwise haplotype probabilities during phasing computations.
	 * 
	 * Additionally, `DProbs` is initialized as a standard vector of floats with length `HAP_NUMBER`,
	 * used for storing diplotype probabilities.
	 */
	phasing_hmm(conditioning_set * C);

	/**
	 * @brief Destructor for the phasing_hmm class.
	 *
	 * Cleans up resources used by the phasing_hmm instance.
	 * Currently, no explicit resource deallocation or cleanup is required
	 * since all members use RAII-compliant containers.
	 */
	~phasing_hmm();

	/**
	 * @brief Reallocates and initializes internal data structures for phasing based on input haplotype states.
	 * 
	 * This function updates the variant type and indexing vectors based on the input haplotypes \p H0 and \p H1,
	 * and the flat variant mask \p flat. It then segments the variant sites for efficient phasing computations,
	 * counts missing and heterozygous variants, and resizes internal buffers accordingly.
	 * 
	 * @param[in] H0 Vector of boolean alleles representing the first haplotype state for polymorphic sites.
	 * @param[in] H1 Vector of boolean alleles representing the second haplotype state for polymorphic sites.
	 * @param[in,out] flat Vector of boolean flags marking variants to be treated as flat heterozygous.
	 * 
	 * @details
	 * The function performs the following steps:
	 * - Clears and rebuilds the internal variant type vectors (\c VAR_TYP, \c VAR_ALT, \c VAR_ABS, \c VAR_REL) by scanning polymorphic sites:
	 *   - Assigns variant types based on allele differences between \p H0 and \p H1 and whether the variant is flagged as flat or low quality.
	 * - Splits variants into segments, where each segment contains up to 4 peak heterozygous variants, to facilitate efficient phasing.
	 * - Counts the number of missing (flat heterozygous) variants.
	 * - Resizes and zero-initializes internal buffers for probability computations, such as \c prob, \c probSumH, \c probSumK, \c phasingProb, and \c imputeProb.
	 * 
	 * This setup is critical before running the phasing HMM algorithms to ensure all structures are consistent with the current variant set.
	 */
	void reallocate(const std::vector < bool > &, const std::vector < bool > &, std::vector < bool > &);

	/**
	 * @brief Perform the forward algorithm pass of the phasing Hidden Markov Model (HMM).
	 *
	 * This function computes the forward probabilities across all variant loci to
	 * estimate haplotype phase probabilities in a genotype phasing model. It accounts
	 * for variant-specific emissions, transition probabilities between states, and
	 * segmentation of variants into blocks for efficient computation.
	 *
	 * @details
	 * The forward pass is based on a Hidden Markov Model where the hidden states represent
	 * haplotype configurations and observed data are genotype variants. The algorithm
	 * sequentially processes variant sites to update the forward probability distribution
	 * over the haplotype states.
	 *
	 * Let:
	 * - \\(L\\) be the total number of variant loci.
	 * - \\(\\alpha_l(s)\\) be the forward probability at locus \\(l\\) and state \\(s\\).
	 * - \\(t_l = P(s_l \\neq s_{l-1})\\) be the transition probability of switching states at locus \\(l\\).
	 * - \\(e_l(s)\\) be the emission probability of observing the data at locus \\(l\\) given state \\(s\\).
	 *
	 * The recursion for forward probabilities is:
	 * \\f[
	 *   \\alpha_l(s) = \\left( \\sum_{s'} \\alpha_{l-1}(s') \\cdot P(s \\mid s') \\right) \\cdot e_l(s)
	 * \\f]
	 *
	 * This implementation handles three variant types:
	 * - Peak heterozygous variants (\\c VAR_PEAK_HET)
	 * - Peak homozygous variants (\\c VAR_PEAK_HOM)
	 * - Flat heterozygous variants (\\c VAR_FLAT_HET)
	 *
	 * @dot
	 * digraph forward_hmm {
	 *   rankdir=TB;
	 *   node [shape=box, style=filled, fillcolor=lightgrey, fontsize=10];
  	 *	 edge [fontsize=8];
	 *
	 *   Start [label="Start: l = 0"];
	 *   VarType [label="Variant type at locus l", shape=diamond, fillcolor=white];
	 *
	 *   PeakHet [label="PEAK_HET"];
	 *   PeakHom [label="PEAK_HOM"];
	 *   FlatHet [label="FLAT_HET"];
	 *   Other [label="Other (error)", shape=diamond, fillcolor=red];
	 *
	 *   IsFirstHet [label="Is first locus?", shape=diamond, fillcolor=white];
	 *   IsFirstHom [label="Is first locus?", shape=diamond, fillcolor=white];
	 *   IsFirstFlat [label="Is first locus?", shape=diamond, fillcolor=white];
	 *
	 *   InitPeakHet [label="INIT_PEAK_HET"];
	 *   RunPeakHet [label="RUN_PEAK_HET"];
	 *   CollapsePeakHet [label="COLLAPSE_PEAK_HET"];
	 *
	 *   InitPeakHom [label="INIT_PEAK_HOM"];
	 *   RunPeakHom [label="RUN_PEAK_HOM"];
	 *   CollapsePeakHom [label="COLLAPSE_PEAK_HOM"];
	 *
	 *   InitFlatHet [label="INIT_FLAT_HET"];
	 *   RunFlatHet [label="RUN_FLAT_HET"];
	 *   CollapseFlatHet [label="COLLAPSE_FLAT_HET"];
	 *
	 *   SegmentEnd [label="Segment end?", shape=diamond, fillcolor=white];
	 *   SumK [label="SUMK()"];
	 *   LastLocus [label="Last locus?", shape=diamond, fillcolor=white];
	 *   TransitionSample [label="TRANS_HAP() & SAMPLE_DIP()"];
	 *   SkipTransition [label="Skip transition"];
	 *   IsFlatHet [label="Variant is FLAT_HET?", shape=diamond, fillcolor=white];
	 *   ImputeFlatHet [label="IMPUTE_FLAT_HET()"];
	 *   Continue [label="Continue"];
	 *   UpdateCounters [label="Update locus and segment counters"];
	 *   LoopCheck [label="l < L-1?", shape=diamond, fillcolor=white];
	 *   End [label="End"];
	 *
	 *   Start -> VarType;
	 *   VarType -> PeakHet [label="PEAK_HET"];
	 *   VarType -> PeakHom [label="PEAK_HOM"];
	 *   VarType -> FlatHet [label="FLAT_HET"];
	 *   VarType -> Other [label="Other"];
	 *
	 *   PeakHet -> IsFirstHet;
	 *   PeakHom -> IsFirstHom;
	 *   FlatHet -> IsFirstFlat;
	 *   Other -> End [label="Throw error", color=red];
	 *
	 *   IsFirstHet -> InitPeakHet [label="Yes"];
	 *   IsFirstHet -> RunPeakHet [label="No", constraint=false];
	 *   IsFirstHet -> CollapsePeakHet [label="No", style=dotted];
	 *
	 *   RunPeakHet -> SegmentEnd;
	 *   InitPeakHet -> SegmentEnd;
	 *   CollapsePeakHet -> SegmentEnd;
	 *
	 *   IsFirstHom -> InitPeakHom [label="Yes"];
	 *   IsFirstHom -> RunPeakHom [label="No", constraint=false];
	 *   IsFirstHom -> CollapsePeakHom [label="No", style=dotted];
	 *
	 *   RunPeakHom -> SegmentEnd;
	 *   InitPeakHom -> SegmentEnd;
	 *   CollapsePeakHom -> SegmentEnd;
	 *
	 *   IsFirstFlat -> InitFlatHet [label="Yes"];
	 *   IsFirstFlat -> RunFlatHet [label="No", constraint=false];
	 *   IsFirstFlat -> CollapseFlatHet [label="No", style=dotted];
	 *
	 *   RunFlatHet -> SegmentEnd;
	 *   InitFlatHet -> SegmentEnd;
	 *   CollapseFlatHet -> SegmentEnd;
	 *
	 *   SegmentEnd -> SumK [label="Yes"];
	 *   SegmentEnd -> UpdateCounters [label="No"];
	 *
	 *   SumK -> LastLocus;
	 *   LastLocus -> TransitionSample [label="No"];
	 *   LastLocus -> SkipTransition [label="Yes"];
	 *
	 *   TransitionSample -> IsFlatHet;
	 *   SkipTransition -> IsFlatHet;
	 *
	 *   IsFlatHet -> ImputeFlatHet [label="Yes"];
	 *   IsFlatHet -> Continue [label="No"];
	 *
	 *   ImputeFlatHet -> UpdateCounters;
	 *   Continue -> UpdateCounters;
	 *
	 *   UpdateCounters -> LoopCheck;
	 *   LoopCheck -> Start [label="Yes"];
	 *   LoopCheck -> End [label="No"];
	 * }
	 * @enddot
	 *
	 * @throws std::runtime_error if an unknown variant type is encountered.
	 *
	 * @see INIT_PEAK_HET(), RUN_PEAK_HET(), COLLAPSE_PEAK_HET()
	 * @see INIT_PEAK_HOM(), RUN_PEAK_HOM(), COLLAPSE_PEAK_HOM()
	 * @see INIT_FLAT_HET(), RUN_FLAT_HET(), COLLAPSE_FLAT_HET()
	 * @see SUMK(), TRANS_HAP(), SAMPLE_DIP(), IMPUTE_FLAT_HET()
	 *
	 * @par References
	 * - Browning BL, Browning SR. "Haplotype phasing: existing methods and new developments." Nat Rev Genet. 2011.
	 * - Durbin R. "Efficient haplotype matching and storage using the positional Burrows–Wheeler transform (PBWT)." Bioinformatics. 2014.
	 * - Li N, Stephens M. "Modeling linkage disequilibrium and identifying recombination hotspots using single-nucleotide polymorphism data." Genetics. 2003.
	 */
	void forward();

	/**
	 * @brief Perform the backward algorithm pass of the phasing Hidden Markov Model (HMM).
	 *
	 * This function computes the backward probabilities over variant loci to complement
	 * the forward pass, enabling full posterior inference of haplotype phases. It
	 * iterates backward through variants, updating the state probabilities according to
	 * the model's transition and emission probabilities.
	 *
	 * @details
	 * The backward pass is a key step in the forward-backward algorithm for HMMs, used
	 * here to estimate the probability of the observed genotype data given the hidden
	 * haplotype states.
	 *
	 * Let:
	 * - \\(L\\) be the total number of variant loci.
	 * - \\(\\beta_l(s)\\) be the backward probability at locus \\(l\\) and state \\(s\\), defined as:
	 * \f[
	 * \beta_l(s) = P(O_{l+1:L} \mid s_l = s)
	 * \f]
	 * where \\(O_{l+1:L}\\) are observations from loci \\(l+1\\) to \\(L\\).
	 *
	 * The recursion for backward probabilities is:
	 * \f[
	 * \beta_l(s) = \sum_{s'} P(s_{l+1} = s' \mid s_l = s) \cdot e_{l+1}(s') \cdot \beta_{l+1}(s')
	 * \f]
	 * where \\(P(s_{l+1} \mid s_l)\\) are the state transition probabilities,
	 * and \\(e_{l+1}(s')\\) are the emission probabilities at locus \\(l+1\\).
	 *
	 * This implementation supports three variant types:
	 * - Peak heterozygous variants (\c VAR_PEAK_HET)
	 * - Peak homozygous variants (\c VAR_PEAK_HOM)
	 * - Flat heterozygous variants (\c VAR_FLAT_HET)
	 *
	 * @dot
	 * digraph BackwardPass {
	 *   rankdir=TB;
	 *   node [shape=box, style=rounded, fontsize=10];
	 *
	 *   Start [label="Start: l = L-1"];
	 *   VarType [shape=diamond, label="Variant type at locus l"];
	 *
	 *   PEAK_HET_Last [shape=diamond, label="Is last locus?"];
	 *   PEAK_HOM_Last [shape=diamond, label="Is last locus?"];
	 *   FLAT_HET_Last [shape=diamond, label="Is last locus?"];
	 *   Error [label="Throw error", shape=box];
	 *
	 *   INIT_P_HET [label="INIT_PEAK_HET"];
	 *   RUN_P_HET [label="RUN_PEAK_HET"];
	 *   COLL_P_HET [label="COLLAPSE_PEAK_HET"];
	 *
	 *   INIT_P_HOM [label="INIT_PEAK_HOM"];
	 *   RUN_P_HOM [label="RUN_PEAK_HOM"];
	 *   COLL_P_HOM [label="COLLAPSE_PEAK_HOM"];
	 *
	 *   INIT_F_HET [label="INIT_FLAT_HET"];
	 *   RUN_F_HET [label="RUN_FLAT_HET"];
	 *   COLL_F_HET [label="COLLAPSE_FLAT_HET"];
	 *
	 *   InSeg1 [shape=diamond, label="Inside segment?"];
	 *   InSeg2 [shape=diamond, label="Inside segment?"];
	 *   InSeg3 [shape=diamond, label="Inside segment?"];
	 *
	 *   SegStart [shape=diamond, label="Start of segment?"];
	 *   SUMK [label="SUMK()"];
	 *   StorePhase [label="Store phasing probabilities"];
	 *
	 *   IsFlatHet [shape=diamond, label="Variant is FLAT_HET?"];
	 *   StoreImp [label="Store imputation probabilities"];
	 *   Continue [label="Continue"];
	 *
	 *   Update [label="Update locus and segment counters"];
	 *   Loop [shape=diamond, label="l > 0?"];
	 *   End [label="End"];
	 *
	 *   Start -> VarType;
	 *   VarType -> PEAK_HET_Last [label="PEAK_HET"];
	 *   VarType -> PEAK_HOM_Last [label="PEAK_HOM"];
	 *   VarType -> FLAT_HET_Last [label="FLAT_HET"];
	 *   VarType -> Error [label="Other"];
	 *
	 *   PEAK_HET_Last -> INIT_P_HET [label="Yes"];
	 *   PEAK_HET_Last -> InSeg1 [label="No"];
	 *   InSeg1 -> RUN_P_HET [label="Yes"];
	 *   InSeg1 -> COLL_P_HET [label="No"];
	 *
	 *   PEAK_HOM_Last -> INIT_P_HOM [label="Yes"];
	 *   PEAK_HOM_Last -> InSeg2 [label="No"];
	 *   InSeg2 -> RUN_P_HOM [label="Yes"];
	 *   InSeg2 -> COLL_P_HOM [label="No"];
	 *
	 *   FLAT_HET_Last -> INIT_F_HET [label="Yes"];
	 *   FLAT_HET_Last -> InSeg3 [label="No"];
	 *   InSeg3 -> RUN_F_HET [label="Yes"];
	 *   InSeg3 -> COLL_F_HET [label="No"];
	 *
	 *   INIT_P_HET -> SegStart;
	 *   RUN_P_HET -> SegStart;
	 *   COLL_P_HET -> SegStart;
	 *   INIT_P_HOM -> SegStart;
	 *   RUN_P_HOM -> SegStart;
	 *   COLL_P_HOM -> SegStart;
	 *   INIT_F_HET -> SegStart;
	 *   RUN_F_HET -> SegStart;
	 *   COLL_F_HET -> SegStart;
	 *
	 *   SegStart -> SUMK [label="Yes"];
	 *   SUMK -> StorePhase;
	 *   SegStart -> StorePhase [label="No"];
	 *
	 *   StorePhase -> IsFlatHet;
	 *   IsFlatHet -> StoreImp [label="Yes"];
	 *   IsFlatHet -> Continue [label="No"];
	 *
	 *   StoreImp -> Update;
	 *   Continue -> Update;
	 *
	 *   Update -> Loop;
	 *   Loop -> Start [label="Yes"];
	 *   Loop -> End [label="No"];
	 * }
	 * @enddot
	 *
	 * @throws std::runtime_error if an unknown variant type is encountered.
	 *
	 * @see INIT_PEAK_HET(), RUN_PEAK_HET(), COLLAPSE_PEAK_HET()
	 * @see INIT_PEAK_HOM(), RUN_PEAK_HOM(), COLLAPSE_PEAK_HOM()
	 * @see INIT_FLAT_HET(), RUN_FLAT_HET(), COLLAPSE_FLAT_HET()
	 * @see SUMK()
	 *
	 * **References**:
	 * - Rabiner LR. "A tutorial on hidden Markov models and selected applications in speech recognition." Proceedings of the IEEE, 1989.
	 * - Browning BL, Browning SR. "Haplotype phasing: existing methods and new developments." Nature Reviews Genetics, 2011.
	 * - Li N, Stephens M. "Modeling linkage disequilibrium and identifying recombination hotspots using single-nucleotide polymorphism data." Genetics, 2003.
	 */
	void backward();

	/**
	 * @brief Phasing algorithm loop for loci from L-1 down to 0.
	 *
	 * This loop processes loci in reverse order, computing transition probabilities,
	 * handling different variant types (PEAK_HET, PEAK_HOM, FLAT_HET), and storing
	 * phasing and imputation probabilities.
	 *
	 * @details
	 * **Algorithm Steps:**
	 * - Iterate from last locus (L-1) to first (0).
	 * - For each locus:
	 *   1. Update absolute and relative locus indices.
	 *   2. Compute transition probability \(t_l\) and complement \(1 - t_l\).
	 *   3. Handle the variant type and position in segment:
	 *      - **PEAK_HET:** INIT, RUN, COLLAPSE phases
	 *      - **PEAK_HOM:** INIT, RUN, COLLAPSE phases
	 *      - **FLAT_HET:** INIT, RUN, COLLAPSE phases
	 *   4. At the start of segment: run `SUMK()` and store phasing probabilities.
	 *   5. For FLAT_HET: store imputation probabilities.
	 *   6. Update locus and segment counters.
	 *
	 * **Graph:**
	 * @dot
	 * digraph Phasing {
	 *   node [shape=box, fontsize=10];
	 *   rankdir=TB;
	 *
	 *   start [label="Start: l = L-1"];
	 *   update_idx [label="Update locus indices"];
	 *   compute_t [label="Compute t_l and 1 - t_l"];
	 *   switch_var [label="switch(VAR_TYP[l])"];
	 *   peak_het [label="PEAK_HET\n(INIT / RUN / COLLAPSE)"];
	 *   peak_hom [label="PEAK_HOM\n(INIT / RUN / COLLAPSE)"];
	 *   flat_het [label="FLAT_HET\n(INIT / RUN / COLLAPSE)"];
	 *   error_case [label="Default: throw error", shape=diamond, color=red];
	 *   sumk [label="If start of segment:\nSUMK() + store phasing probs"];
	 *   impute [label="If flat het:\nStore imputation probs"];
	 *   update_counters [label="Update locus & segment counters"];
	 *   loop [label="l-- (next locus)"];
	 *   end [label="End"];
	 *
	 *   start -> update_idx -> compute_t -> switch_var;
	 *   switch_var -> peak_het;
	 *   switch_var -> peak_hom;
	 *   switch_var -> flat_het;
	 *   switch_var -> error_case;
	 *   peak_het -> sumk;
	 *   peak_hom -> sumk;
	 *   flat_het -> sumk -> impute -> update_counters -> loop -> update_idx;
	 *   error_case -> end;
	 *   loop -> end [label="if l < 0"];
	 * }
	 * @enddot
	 *
	 * **References:**
	 * - Browning SR, Browning BL. Haplotype phasing: existing methods and new developments. *Nat Rev Genet*. 2011;12(10):703–714.
	 * - Delaneau O, Zagury JF, Marchini J. Improved whole-chromosome phasing for disease and population genetic studies. *Nat Methods*. 2013;10(1):5–6.
	 */
	void rephaseHaplotypes(std::vector < bool > &, std::vector < bool > &, std::vector < bool > &);
};

inline
void phasing_hmm::INIT_PEAK_HET(int curr_het)
{
	const std::array <__m256, 2 > emits = {_mm256_load_ps(&EMIT0[curr_het][0]),_mm256_load_ps(&EMIT1[curr_het][0])};
    __m256 _sum = _mm256_set1_ps(0.0f);
	for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER)
	{
		const bool ah = C->Hvar.get(curr_rel_locus, k);
		_sum = _mm256_add_ps(_sum, emits[ah]);
		_mm256_store_ps(&prob[i], emits[ah]);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = horizontal_add(_sum);
}

inline
void phasing_hmm::INIT_PEAK_HOM(bool ag)
{
	const std::array <__m256, 2 > emits = {_mm256_set1_ps(1.0f),_mm256_set1_ps(C->ed_phs/C->ee_phs)};
    __m256 _sum = _mm256_set1_ps(0.0f);
	for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER)
	{
		const bool ag_ah = C->Hvar.get(curr_rel_locus, k)!=ag;
		_sum = _mm256_add_ps(_sum, emits[ag_ah]);
		_mm256_store_ps(&prob[i], emits[ag_ah]);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = horizontal_add(_sum);
}

inline
void phasing_hmm::INIT_FLAT_HET() {
	fill(prob.begin(), prob.end(), 1.0f/(HAP_NUMBER * C->n_states));
	fill(probSumH.begin(), probSumH.end(), 1.0f/HAP_NUMBER);
    probSumT = 1.0f;
}

inline
void phasing_hmm::RUN_PEAK_HET(int curr_het)
{
	const __m256 _tFreq = _mm256_mul_ps(_mm256_load_ps(&probSumH[0]), _mm256_set1_ps(yt / (C->n_states * probSumT)));
	const __m256 _nt = _mm256_set1_ps(nt / probSumT);
	const std::array <__m256, 2 > emits = {_mm256_load_ps(&EMIT0[curr_het][0]),_mm256_load_ps(&EMIT1[curr_het][0])};
    __m256 _sum = _mm256_set1_ps(0.0f);
	for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER)
	{
		const bool ah = C->Hvar.get(curr_rel_locus, k);
		const __m256 _prob_curr = _mm256_mul_ps(_mm256_fmadd_ps(_mm256_load_ps(&prob[i]), _nt, _tFreq), emits[ah]);
		_sum = _mm256_add_ps(_sum, _prob_curr);
		_mm256_store_ps(&prob[i], _prob_curr);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = horizontal_add(_sum);
}

inline
void phasing_hmm::RUN_PEAK_HOM(bool ag)
{
	const __m256 _tFreq = _mm256_mul_ps(_mm256_load_ps(&probSumH[0]), _mm256_set1_ps(yt / (C->n_states * probSumT)));
	const __m256 _nt = _mm256_set1_ps(nt / probSumT);
    const __m256 _mism = _mm256_set1_ps(C->ed_phs/C->ee_phs);
    __m256 _sum = _mm256_set1_ps(0.0f);
	for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER)
	{
		const bool ah = C->Hvar.get(curr_rel_locus, k);
		const __m256 _prob_prev = _mm256_load_ps(&prob[i]);
		__m256 _prob_curr = _mm256_fmadd_ps(_prob_prev, _nt, _tFreq);
		if (ag!=ah) _prob_curr = _mm256_mul_ps(_prob_curr, _mism);
		_sum = _mm256_add_ps(_sum, _prob_curr);
		_mm256_store_ps(&prob[i], _prob_curr);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = horizontal_add(_sum);
}

inline
void phasing_hmm::RUN_FLAT_HET()
{
	const __m256 _tFreq = _mm256_mul_ps(_mm256_load_ps(&probSumH[0]), _mm256_set1_ps(yt / (C->n_states * probSumT)));
	const __m256 _nt = _mm256_set1_ps(nt / probSumT);
    __m256 _sum = _mm256_set1_ps(0.0f);
	for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER)
	{
		const __m256 _prob_curr = _mm256_fmadd_ps(_mm256_load_ps(&prob[i]), _nt, _tFreq);
		_sum = _mm256_add_ps(_sum, _prob_curr);
		_mm256_store_ps(&prob[i], _prob_curr);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = horizontal_add(_sum);
}

inline
void phasing_hmm::COLLAPSE_PEAK_HET(int curr_het)
{
	const __m256 _tFreq = _mm256_mul_ps(_mm256_load_ps(&probSumH[0]), _mm256_set1_ps(yt / (C->n_states * probSumT)));
	const __m256 _nt = _mm256_set1_ps(nt / probSumT);
	const std::array <__m256, 2 > emits = {_mm256_load_ps(&EMIT0[curr_het][0]),_mm256_load_ps(&EMIT1[curr_het][0])};
    __m256 _sum = _mm256_set1_ps(0.0f);
	for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER)
	{
		const bool ah = C->Hvar.get(curr_rel_locus, k);
		const __m256 _prob_curr = _mm256_mul_ps(_mm256_fmadd_ps(_mm256_set1_ps(probSumK[k]), _nt, _tFreq), emits[ah]);
		_sum = _mm256_add_ps(_sum, _prob_curr);
		_mm256_store_ps(&prob[i], _prob_curr);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = horizontal_add(_sum);
}

inline
void phasing_hmm::COLLAPSE_PEAK_HOM(bool ag)
{
	const __m256 _tFreq = _mm256_mul_ps(_mm256_load_ps(&probSumH[0]), _mm256_set1_ps(yt / (C->n_states * probSumT)));
   	const __m256 _nt = _mm256_set1_ps(nt / probSumT);
    __m256 _sum = _mm256_set1_ps(0.0f);
    const __m256 _mism = _mm256_set1_ps(C->ed_phs/C->ee_phs);

   	for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER)
   	{
   		const bool ah = C->Hvar.get(curr_rel_locus, k);
   		__m256 _prob_curr = _mm256_fmadd_ps(_mm256_set1_ps(probSumK[k]), _nt, _tFreq);
		if (ag!=ah) _prob_curr = _mm256_mul_ps(_prob_curr, _mism);
   		_sum = _mm256_add_ps(_sum, _prob_curr);
   		_mm256_store_ps(&prob[i], _prob_curr);
   	}
   	_mm256_store_ps(&probSumH[0], _sum);
   	probSumT = horizontal_add(_sum);
}

inline
void phasing_hmm::COLLAPSE_FLAT_HET ()
{
	const __m256 _tFreq = _mm256_mul_ps(_mm256_load_ps(&probSumH[0]), _mm256_set1_ps(yt / (C->n_states * probSumT)));
	const __m256 _nt = _mm256_set1_ps(nt / probSumT);
    __m256 _sum = _mm256_set1_ps(0.0f);
	for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER)
	{
   		__m256 _prob_curr = _mm256_fmadd_ps(_mm256_set1_ps(probSumK[k]), _nt, _tFreq);
		_sum = _mm256_add_ps(_sum, _prob_curr);
		_mm256_store_ps(&prob[i], _prob_curr);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = horizontal_add(_sum);
}


inline
void phasing_hmm::SUMK() {
	for(int k = 0, i = 0 ; k != C->n_states ; ++k, i += HAP_NUMBER)
		probSumK[k] = horizontal_add(_mm256_load_ps(&prob[i]));
}

inline
bool phasing_hmm::TRANS_HAP()
{
	const int states_haps = C->n_states*HAP_NUMBER;
	sumHProbs = 0.0f;
	yt = C->getTransition(VAR_ABS[curr_idx_locus], VAR_ABS[curr_idx_locus+1]);
	nt = 1.0f - yt;
	const __m256 _fact2 = _mm256_set1_ps( nt / probSumT);
	int h1 = 0, j=0;
	for (; h1 < HAP_NUMBER ; h1++, j += HAP_NUMBER)
	{
		const __m256 _fact1 = _mm256_set1_ps((probSumH[h1]/probSumT) * yt / C->n_states);
		__m256 _sum = _mm256_set1_ps(0.0f);
		for(int k=0, i=0; k != C->n_states ; ++k, i += HAP_NUMBER)
		{
			const __m256 _prob0 = _mm256_fmadd_ps(_mm256_set1_ps(prob[i+h1]), _fact2, _fact1);
			_sum = _mm256_add_ps(_sum, _mm256_mul_ps(_prob0, _mm256_load_ps(&phasingProb[(curr_segment_index+1)*states_haps+i])));
		}
		_mm256_store_ps(&HProbs[j], _sum);
		sumHProbs += horizontal_add(_sum);
	}
	return (std::isnan(sumHProbs) || std::isinf(sumHProbs) || sumHProbs < std::numeric_limits<float>::min());
}

inline
bool phasing_hmm::SAMPLE_DIP() {
	sumDProbs = 0.0f;
	for (int d = 0 ; d < HAP_NUMBER ; d ++) {
		int prev_h0 = dip_sampled[curr_segment_index];
		int prev_h1 = HAP_NUMBER - dip_sampled[curr_segment_index] - 1;
		DProbs[d] = (HProbs[prev_h0 * HAP_NUMBER + d] / sumHProbs) * (HProbs[prev_h1 * HAP_NUMBER + (HAP_NUMBER - d - 1)] / sumHProbs);
		sumDProbs += DProbs[d];
	}
	if (std::isnan(sumDProbs) || std::isinf(sumDProbs) || sumDProbs < std::numeric_limits<float>::min()) return true;
	dip_sampled[curr_segment_index+1] = rng.sample(DProbs, sumDProbs);
	return false;
}

inline
void phasing_hmm::IMPUTE_FLAT_HET()
{
	const int states_haps = C->n_states*HAP_NUMBER;
	const __m256 _one = _mm256_set1_ps(1.0f);
	const __m256 _zero = _mm256_set1_ps(0.0f);
	__m256 _scaleR = _mm256_load_ps(&imputeProbSum[curr_missing_locus*HAP_NUMBER]);
	__m256 _scaleL = _mm256_load_ps(&probSumH[0]);
	_scaleR = _mm256_div_ps(_one, _scaleR);
	_scaleL = _mm256_div_ps(_one, _scaleL);
	std::array <__m256, 2 > sums = {_mm256_set1_ps(0.0f),_mm256_set1_ps(0.0f)};

	for(int k = 0, i = 0 ; k !=  C->n_states ; ++k, i += HAP_NUMBER) {
		const bool ah = C->Hvar.get(curr_rel_locus, k);
		const __m256 _p1 = _mm256_mul_ps(_mm256_load_ps(&imputeProb[curr_missing_locus*states_haps + i]), _scaleR);
		const __m256 _p2 = _mm256_mul_ps(_mm256_load_ps(&prob[i]), _scaleL);
		sums[ah] = _mm256_add_ps(sums[ah], _mm256_mul_ps(_p1,_p2));
	}
	const __m256 _norm_sum = _mm256_div_ps(sums[1], _mm256_add_ps(sums[0], sums[1]));
	const __m256 _clamp = _mm256_max_ps(_zero, _mm256_min_ps(_one, _norm_sum));
	_mm256_store_ps(&imputeProbOf1s[curr_missing_locus*HAP_NUMBER], _clamp);
}

#endif

