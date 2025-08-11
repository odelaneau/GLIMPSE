/**
 * @file haplotype_set.h
 * @brief Declaration of the haplotype_set class for managing haplotype data and PBWT structures.
 * 
 * This class extends ref_haplotype_set and contains data and methods for
 * handling haplotypes from target and reference samples, including rare allele
 * tracking, ploidy management, and PBWT-based haplotype matching.
 * 
 * @author Simone Rubinacci, Olivier Delaneau
 * @license MIT License
 */

#ifndef _HAPLOTYPE_SET_H
#define _HAPLOTYPE_SET_H

#include "otools.h"
#include "checksum_utils.h"

#include "bitmatrix.h"
#include "ref_haplotype_set.h"
#include "containers/genotype_set.h"

/**
 * @class haplotype_set
 * @brief Extends ref_haplotype_set to handle both reference and target haplotypes with PBWT.
 * 
 * Contains:
 * - Counts of haplotypes and samples
 * - Haplotype data including rare variant storage
 * - Ploidy information and mappings
 * - PBWT-related data structures and matching algorithms
 * - Statistics for selection lengths and counters for various PBWT events
 */
class haplotype_set : public ref_haplotype_set {
public:
    /// Total number of haplotypes (target + reference)
    unsigned int n_tot_haps;

    /// Number of haplotypes in target samples only
    unsigned int n_tar_haps;

    /// Number of target samples
    unsigned int n_tar_samples;

    /// Bitmatrix storing haplotype variants for target samples
    bitmatrix HvarTar;

    /// Rare alleles per haplotype
    std::vector<std::vector<int>> ShapTar;

    /// Rare alleles per variant
    std::vector<std::vector<int>> SvarTar;

    /// Rare alleles per individual from genotype likelihoods
    std::vector<std::vector<int>> SindTarGL;

    /// Genetic map positions in centimorgans
    std::vector<float> cm_pos;

    /// Format ploidy indicator (1=haploid, 2=diploid, -2=mixed)
    int fploidy;

    /// Maximum ploidy found among samples
    int max_ploidy;

    /// Vector storing ploidy per target sample
    std::vector<int> tar_ploidy;

    /// Mapping from target individuals to haplotype IDs
    std::vector<int> tar_ind2hapid;

    /// Mapping from haplotype IDs to target individual indices
    std::vector<int> tar_hapid2ind;

    /// PBWT depth parameter controlling selection size
    int pbwt_depth;

    /// PBWT modulo in centimorgans, controlling segment length
    float pbwt_modulo_cm;

    /// Vector related to PBWT variant arrays (size n_ref_haps+1)
    std::vector<int> pbwt_array_V;

    /// PBWT indices per target haplotype
    std::vector<int> pbwt_index;

    /// PBWT indices for smaller subsets (rare variants)
    std::vector<int> pbwt_small_index;

    /// Additional vectors for small PBWT variant arrays and rare haplotype tracking
    std::vector<int> pbwt_small_V;
    std::vector<int> rareTarHaps;

    /// PBWT lower and upper bounds for matching windows per haplotype
    std::vector<int> f_k;
    std::vector<int> g_k;

    /// Similar bounds for small PBWT (rare variants)
    std::vector<int> f_k_small;
    std::vector<int> g_k_small;

    /// Track last reset position for each haplotype
    std::vector<int> last_reset;

    /// Track last rare variant position for each haplotype
    std::vector<int> last_rare;

    /// Flags indicating if a site is stored in PBWT
    std::vector<bool> pbwt_stored;

    /// Group boundaries of PBWT positions
    std::vector<int> pbwt_grp;

    /// Initialization parameters for PBWT
    int Kinit;
    int Kpbwt;

    /// Current PBWT depth used
    int K;

    /// Number of stored PBWT groups
    int nstored;

    /// Counters for different PBWT matching and selection events
    int counter_gf;
    int counter_sel_gf;
    int counter_rare_restarts;

    /// Nested vectors storing matched haplotypes states per target sample and depth
    std::vector<std::vector<std::vector<int>>> pbwt_states;

    /// Set of initial states per target haplotype
    std::vector<std::set<int>> init_states;

    /// List of haplotype states per target haplotype
    std::vector<std::vector<int>> list_states;

    /// Vector storing alleles for target haplotypes (0 or 1)
    std::vector<unsigned char> tar_hap;

    /// Map from rare variant indices to internal IDs
    std::map<int, int> rare_idx_to_id;

    /// Statistical trackers for segment length selections (in cM)
    stats1D length_sel_mod;
    stats1D length_sel_gf;
    stats1D length_sel_gf_rare;



	//CONSTRUCTOR/DESTRUCTOR/INITIALIZATION
	/**
	 * @brief Default constructor for the haplotype_set class.
	 *
	 * Initializes all counters and parameters related to haplotype storage,
	 * PBWT search, and ploidy to their default values.
	 *
	 * @details
	 * The constructor sets:
	 * - Site counters:
	 *   - `n_tot_sites` = 0 — total number of sites.
	 *   - `n_rar_sites` = 0 — number of rare variant sites.
	 *   - `n_com_sites` = 0 — number of common variant sites.
	 *
	 * - Haplotype and sample counters:
	 *   - `n_tot_haps` = 0 — total haplotypes stored.
	 *   - `n_tar_haps` = 0 — target haplotypes count.
	 *   - `n_ref_haps` = 0 — reference haplotypes count.
	 *   - `n_tar_samples` = 0 — target samples count.
	 *
	 * - PBWT and algorithm parameters:
	 *   - `Kinit` = 1000 — initial number of PBWT states.
	 *   - `K` = 1000 — maximum number of states retained.
	 *   - `Kpbwt` = 2000 — number of PBWT states for lookup.
	 *   - `nstored` = 0 — stored haplotypes counter.
	 *   - `pbwt_depth` = 0 — PBWT search depth.
	 *   - `pbwt_modulo_cm` = 0.0 — PBWT interval in centimorgans.
	 *
	 * - Ploidy and counters:
	 *   - `max_ploidy` = 2 — maximum allowed ploidy.
	 *   - `fploidy` = 2 — default ploidy factor.
	 *   - `counter_sel_gf` = 0 — counter for selected genotype fields.
	 *   - `counter_gf` = 0 — general genotype field counter.
	 *   - `counter_rare_restarts` = 0 — counter for rare restart events.
	 */
	haplotype_set();

	/**
	 * @brief Destructor for the haplotype_set class.
	 *
	 * This destructor resets key statistical counters and PBWT configuration values
	 * to their default zero state. While no dynamic memory is explicitly released here
	 * (as memory management is handled elsewhere), this ensures internal counters
	 * and parameters are cleared when the object is destroyed.
	 *
	 * The following members are reset:
	 * - **Site counters**:
	 *   - `n_tot_sites` – Total number of variant sites.
	 *   - `n_rar_sites` – Number of rare variant sites.
	 *   - `n_com_sites` – Number of common variant sites.
	 *
	 * - **Haplotype counters**:
	 *   - `n_tot_haps` – Total number of haplotypes.
	 *   - `n_tar_haps` – Number of target haplotypes.
	 *   - `n_ref_haps` – Number of reference haplotypes.
	 *
	 * - **PBWT parameters**:
	 *   - `pbwt_depth` – PBWT search depth.
	 *   - `pbwt_modulo_cm` – PBWT interval in centimorgans.
	 *
	 * @note
	 * This function does not perform deep cleanup of dynamically allocated
	 * resources. Such cleanup, if required, must be handled by the caller or
	 * other class methods.
	 */
	virtual ~haplotype_set();

	/**
	 * @brief Allocates memory and initializes data structures for haplotype storage and PBWT state tracking.
	 *
	 * This method sets up haplotype-related containers based on the number of target and reference haplotypes,
	 * total sites, and sample counts. It ensures that all major internal vectors are sized appropriately and ready
	 * for subsequent data assignment.
	 *
	 * Steps performed:
	 * 1. **Total haplotype calculation**:
	 *    - `n_tot_haps` is computed as the sum of `n_ref_haps` and `n_tar_haps`.
	 *
	 * 2. **Reference haplotype structures**:
	 *    - `HvarRef` – Allocated with `n_com_sites` × `n_ref_haps`.
	 *    - `ShapRef` – Sized as a vector of length `n_ref_haps`, each containing an empty `std::vector<int>`.
	 *    - `SvarRef` – Sized as a vector of length `n_tot_sites`, each containing an empty `std::vector<int>`.
	 *
	 * 3. **Target haplotype structures**:
	 *    - `HvarTar` – Allocated with `n_com_sites` × `n_tar_haps`.
	 *    - `ShapTar` – Sized as a vector of length `n_tar_haps`, each containing an empty `std::vector<int>`.
	 *    - `SvarTar` – Sized as a vector of length `n_tot_sites`, each containing an empty `std::vector<int>`.
	 *    - `SindTarGL` – Sized as a vector of length `n_tar_samples`, each containing an empty `std::vector<int>`.
	 *
	 * 4. **PBWT state initialization**:
	 *    - `pbwt_states` – A 3D vector sized `n_tar_samples`, each entry initialized as an empty 2D vector of `int`.
	 *    - `init_states` – A vector of `std::set<int>` sized `n_tar_samples`, each initialized as empty.
	 *    - `list_states` – A vector of `std::vector<int>` sized `n_tar_haps`, each initialized as empty.
	 *
	 * @note This method assumes that class member variables such as `n_com_sites`, `n_ref_haps`,
	 *       `n_tar_haps`, `n_tot_sites`, and `n_tar_samples` have been set before calling.
	 */
	void allocate();

	/**
	 * @brief Allocates memory for haplotype-related data structures, excluding genotype data.
	 *
	 * This method initializes and allocates containers for storing reference haplotypes,
	 * target haplotypes, and PBWT (Positional Burrows–Wheeler Transform) states.  
	 * It is intended for scenarios where only haplotype data (and not genotype likelihoods) 
	 * are needed.
	 *
	 * **Key actions performed:**
	 * - Calculates the total number of haplotypes (`n_tot_haps`) from the sum of reference
	 *   and target haplotypes.
	 * - Allocates variant storage structures for reference (`SvarRef`) and target (`SvarTar`)
	 *   haplotypes.
	 * - Allocates storage for target haplotype sequences (`ShapTar`), target variant indices
	 *   (`SindTarGL`), and PBWT-related states (`pbwt_states`, `init_states`, `list_states`).
	 *
	 * @note `HvarTar` uses a custom allocation method, while the others are standard STL containers.
	 * @note PBWT states are initialized as empty for each target sample, and initial states are
	 *       stored as empty sets.
	 *
	 * @pre Member variables such as `n_ref_haps`, `n_tar_haps`, `n_tot_sites`, 
	 *      `n_com_sites`, `n_tar_samples` must be correctly initialized before calling this function.
	 * @post All haplotype-related storage containers are allocated and ready for population.
	 */
	void allocate_hap_only();


	//ROUTINES
	/**
	 * @brief Initializes the list of rare target variants for each individual.
	 *
	 * This function identifies rare variants (targets) based on genotype likelihoods (GLs) 
	 * for each individual in the given genotype set. It filters out common variants and 
	 * those flagged as low quality, then calculates posterior genotype probabilities 
	 * and determines if a site should be considered a rare target.
	 *
	 * @param G The genotype set containing individuals, genotype likelihoods, and other metadata.
	 * @param M The variant map containing variant positions, allele counts, and quality flags.
	 *
	 * The method performs the following steps:
	 * - Skips initialization if `Kinit == 0`.
	 * - For each individual, iterates through all variant sites:
	 *   - Skips sites flagged as "flat" (uninformative).
	 *   - Ignores common variants and low-quality (LQ) sites.
	 *   - Compares GL values to identify likely non-reference genotypes.
	 *   - Calculates normalized posterior probabilities from unphred GL scores.
	 *   - Adjusts probabilities according to allele frequency (AF).
	 *   - Determines rare targets based on the adjusted probabilities.
	 * - Stores rare variant site indices for each individual in `SindTarGL`.
	 * - Collects statistics on the number of rare targets per individual and logs the mean.
	 *
	 * @note
	 *  - The method assumes diploid data but contains partial handling for haploid (`ploidy=1`).
	 *  - Uses `unphred` to convert GL Phred scores to probabilities.
	 *  - Relies on `flag_common`, `major_alleles`, `SindTarGL`, and `n_ref_haps` as class members.
	 *
	 * @warning 
	 * The ploidy handling for haploid genotypes (`ploidy=1`) is marked as TODO/FIXME and may require additional checks.
	 *
	 * @see genotype_set, variant_map
	 */
	void initRareTar(const genotype_set & G, const variant_map& M);

	/**
	 * @brief Updates target haplotypes with current genotype calls.
	 *
	 * This method synchronizes the internal haplotype representation (`HvarTar` for common variants
	 * and `ShapTar` for rare variants) with the genotype calls from the given genotype set.
	 *
	 * @param G The genotype set containing phased haplotypes for all individuals.
	 *
	 * The update process:
	 * - Iterates over all individuals in the genotype set.
	 * - Retrieves each individual's ploidy and corresponding haplotype index (`hapid`).
	 * - Clears rare allele storage (`ShapTar`) for the current individual's haplotypes.
	 * - Iterates over all variant sites:
	 *   - **Common variants (`flag_common[v] == true`)**:
	 *     - Writes the haplotype allele values (`H0`, `H1`) into the compressed haplotype matrix `HvarTar`.
	 *   - **Rare variants (`flag_common[v] == false`)**:
	 *     - If the allele differs from the major allele, records the variant index in `ShapTar` for that haplotype.
	 * - Measures and reports the update time using `tac` and `vrb`.
	 *
	 * @note
	 *  - This method updates only the target haplotypes, not the reference haplotypes.
	 *  - `HvarTar` stores common variant alleles in a memory-efficient form.
	 *  - `ShapTar` stores the positions of non-reference rare alleles for each haplotype.
	 *
	 * @warning
	 *  - Assumes that `tar_ind2hapid` has been initialized correctly to map individuals to haplotype indices.
	 *  - Requires that `major_alleles` and `flag_common` are correctly populated before calling.
	 *
	 * @see initRareTar()
	 */
	void updateHaplotypes(const genotype_set &);

	/**
	 * @brief Transposes the representation of rare variants for reference haplotypes.
	 *
	 * This method reorganizes the rare variant data for reference haplotypes
	 * from a haplotype-centric structure (`ShapRef`) into a site-centric structure (`SvarRef`).
	 * In the resulting structure, each site contains a list of haplotypes that carry
	 * the rare allele at that site.
	 *
	 * **Workflow:**
	 * - Clears the existing `SvarRef` entries for all sites.
	 * - Iterates through each reference haplotype.
	 * - For each rare variant index in `ShapRef[h]`, appends the haplotype ID `h`
	 *   to the corresponding site's vector in `SvarRef`.
	 *
	 * @note This transformation improves lookup efficiency for operations that
	 *       need to know all haplotypes carrying a rare variant at a given site.
	 *
	 * @pre `ShapRef` must contain valid rare variant indices for each haplotype.
	 * @post `SvarRef` will be populated so that `SvarRef[l]` contains all
	 *       reference haplotypes with the rare allele at site `l`.
	 *
	 * @warning This method clears all existing data in `SvarRef`.
	 *
	 * @complexity
	 * - **Time:** O(N_ref_haps × Avg_rare_variants_per_hap)
	 * - **Space:** O(N_tot_sites + total_rare_variant_entries)
	 */
	void transposeRareRef();

	/**
	 * @brief Transposes the representation of rare variants for target haplotypes.
	 *
	 * This method reorganizes the rare variant data for target haplotypes
	 * from a haplotype-centric structure (`ShapTar`) into a site-centric structure (`SvarTar`).
	 * After transposition, each site will store a list of target haplotypes that carry
	 * the rare allele at that site.
	 *
	 * **Workflow:**
	 * - Clears the existing `SvarTar` entries for all sites.
	 * - Iterates over all target haplotypes.
	 * - For each rare variant index in `ShapTar[h]`, appends the haplotype ID `h`
	 *   to the corresponding site's vector in `SvarTar`.
	 *
	 * @note This transformation is used to optimize operations that require
	 *       quick access to all target haplotypes with a rare allele at a specific site.
	 *
	 * @pre `ShapTar` must contain valid rare variant indices for each haplotype.
	 * @post `SvarTar[l]` will contain all target haplotypes with the rare allele at site `l`.
	 *
	 * @warning This method clears all existing data in `SvarTar`.
	 *
	 * @complexity
	 * - **Time:** O(N_tar_haps × Avg_rare_variants_per_hap)
	 * - **Space:** O(N_tot_sites + total_rare_variant_entries)
	 */
	void transposeRareTar();

	//Init
	/**
	 * @brief Perform initial selection of reference haplotypes for rare variants based on genotype likelihoods (GL).
	 * 
	 * This method initializes the `init_states` sets for each target sample by selecting reference haplotypes
	 * associated with rare variants identified in the target sample's genotype likelihood data (`SindTarGL`).
	 * The goal is to select a subset of reference haplotypes related to rare alleles for the initial PBWT state.
	 * 
	 * Selection strategy:
	 * - If `Kinit == 0`, the function exits early without selecting.
	 * - For each target sample:
	 *   - If the number of rare variants (`SindTarGL[ind]`) exceeds 80% of `Kinit`, it sorts rare variants by minor allele count (MAC)
	 *     and samples reference haplotypes linked to the rare variants with the smallest MAC values until ~80% of `Kinit` is reached.
	 *   - Otherwise, it distributes the sampling evenly across rare variants.
	 * - If the `init_states` set size is still less than `Kinit`, fills the remainder by uniform random sampling from all reference haplotypes.
	 * 
	 * @param M The variant map containing minor allele counts and site information.
	 * 
	 * @pre
	 * - `Kinit` must be non-negative.
	 * - `SindTarGL` must contain rare variant indices per target sample.
	 * - `SvarRef` must be populated, mapping rare variants to reference haplotypes.
	 * 
	 * @post
	 * - `init_states[ind]` for each target sample will contain up to `Kinit` selected reference haplotype indices.
	 * - `SindTarGL` will be cleared and shrunk.
	 * 
	 * @note
	 * - Uses `std::sample` with fixed sample sizes to randomly select haplotypes.
	 * - This function assumes `rng.randomEngine` is a valid random engine instance.
	 * 
	 * @complexity
	 * - Roughly O(n_tar_samples * max_kinit * log(k_init)) due to sorting and sampling.
	 */
	void performSelection_RARE_INIT_GL(const variant_map & M);

	/**
	 * @brief Reads and loads PBWT state lists for each target haplotype from a specified file.
	 * 
	 * Each line in the file corresponds to a target haplotype's list of states, 
	 * which are integer identifiers of reference haplotypes associated with that target haplotype.
	 * 
	 * @param file_list Path to the input file containing lists of states per haplotype.
	 * 
	 * @throws std::runtime_error If the file cannot be opened.
	 * @throws std::runtime_error If any state value is out of valid range (negative or >= number of reference haplotypes).
	 * @throws std::runtime_error If any haplotype's state list is empty.
	 * @throws std::runtime_error If the number of lines (haplotypes) read is less than expected (`n_tar_haps`).
	 * 
	 * @pre
	 * - `n_tar_haps` must be set to the expected number of target haplotypes.
	 * - `n_ref_haps` must be set to the number of reference haplotypes.
	 * - `list_states` must be initialized with `n_tar_haps` elements.
	 * 
	 * @post
	 * - `list_states` is populated with the state lists for each target haplotype.
	 * 
	 * @note
	 * - Lines in the input file should contain whitespace-separated integers.
	 * - Each integer must be in the range [0, n_ref_haps - 1].
	 */
	void read_list_states(const std::string file_list);


	//PBWT ROUTINES
	/**
	 * @brief Allocates and initializes data structures for PBWT (Positional Burrows-Wheeler Transform) based haplotype matching.
	 * 
	 * This function sets up the PBWT arrays, parameters, and internal states required for subsequent PBWT operations.
	 * It adapts parameters such as `pbwt_modulo_cm` based on the genetic map length and the requested PBWT depth.
	 * 
	 * @param _pbwt_depth The depth of PBWT (number of conditioning haplotypes per target haplotype).
	 * @param _pbwt_modulo_cm The genetic distance interval in centiMorgans (cM) at which PBWT arrays are stored.
	 * @param M The variant map containing positional and quality information for variants.
	 * @param G The genotype set containing reference and target genotypes.
	 * @param _Kinit The initial number of conditioning haplotypes used for initialization.
	 * @param _Kpbwt The number of conditioning haplotypes used in PBWT-based imputation.
	 * 
	 * @note
	 * - If `_Kpbwt` is zero or greater than or equal to the number of reference haplotypes, PBWT allocation is skipped.
	 * - This function initializes multiple PBWT arrays and state tracking vectors.
	 * - The function adjusts `pbwt_modulo_cm` if the imputation region is too small to ensure meaningful PBWT grouping.
	 * - Memory reservation and vector initializations are done based on the sample ploidies and variant counts.
	 * 
	 * @post
	 * - PBWT data structures and state vectors are allocated and initialized.
	 * - Genetic map positions (`cm_pos`) are computed and adjusted.
	 * - PBWT grouping (`pbwt_grp`) is established for common high-quality sites.
	 * 
	 * @see build_sparsePBWT()
	 */
	void allocatePBWT(const int _pbwt_depth, const float _pbwt_modulo_cm, const variant_map & V, const genotype_set & G, const int _Kinit, const int _Kpbwt);
	
	/**
	 * @brief Performs haplotype matching from the compressed sparse PBWT data structures for target haplotypes.
	 * 
	 * This function implements haplotype selection based on the small (compressed) PBWT representation, 
	 * considering both common and rare variants. It updates PBWT state arrays used for genotype imputation.
	 * 
	 * @param M The variant map containing variant positions and related metadata.
	 * @param main_iteration Boolean flag indicating whether this is the main iteration of the algorithm.
	 *                       (Currently unused but may control behavior in extended versions.)
	 * 
	 * @details
	 * - The function early-exits if PBWT is disabled (Kpbwt=0 or Kpbwt >= n_ref_haps).
	 * - Uses random seeding to initialize PBWT indexes.
	 * - Distinguishes common vs rare variants to invoke appropriate selection and PBWT reading functions.
	 * - Tracks various statistics for debugging and performance monitoring (commented out).
	 * 
	 * @note The function relies on several internal helper functions such as:
	 * - init_common(), init_rare()
	 * - read_full_pbwt_av(), read_small_pbwt_av()
	 * - select_common_pd_fg(), select_rare_pd_fg()
	 * 
	 * These helpers encapsulate the detailed logic for PBWT state initialization, reading, and selection.
	 * 
	 * @post
	 * - PBWT states for all target haplotypes are updated according to the compressed PBWT data.
	 * - Internal selection statistics are accumulated for performance analysis.
	 * 
	 * @see init_common(), init_rare(), read_full_pbwt_av(), read_small_pbwt_av()
	 */
	void matchHapsFromCompressedPBWTSmall(const variant_map & V, const bool main_iteration);

	/**
	 * @brief Initializes PBWT states for common variants at a given site.
	 * 
	 * This method updates the PBWT interval and index arrays for target haplotypes
	 * based on the common variant position indexed by `k` and the local index `l`.
	 * It handles the transition between rare and common haplotype sets and performs
	 * haplotype selection when intervals collapse.
	 * 
	 * @param k The index of the current site in the total variant sites.
	 * @param l The index in the compressed PBWT structure corresponding to the current site.
	 * @param prev_ref_rac_l_com The previous recombination adjusted count (reference) for common variants, used for selection.
	 * 
	 * @details
	 * - The function first checks if the small index array `A_small_idx[l]` is empty; if so, it returns early.
	 * - Iterates over all target haplotypes and distinguishes between rare and common haplotypes.
	 * - For rare haplotypes, uses specialized small PBWT indexes and intervals.
	 * - For common haplotypes, adjusts intervals using `A_small_idx[l]` and performs boundary checks.
	 * - When the haplotype interval collapses (g_dash <= f_dash), it triggers the haplotype selection routine `selectK()`.
	 * - Updates PBWT arrays and indexes to maintain consistency for the next iteration.
	 * 
	 * @post
	 * - PBWT state arrays (`pbwt_index`, `f_k`, `g_k`) and `pbwt_array_A` are updated for the current site.
	 * - Counters for selection events (`counter_sel_gf`) and length statistics are updated.
	 * 
	 * @note This function relies on member variables such as:
	 * - `rareTarHaps`, `pbwt_small_index`, `f_k_small`, `g_k_small`, `last_reset`, `last_rare`, `tar_hap`
	 * - `pbwt_index`, `f_k`, `g_k`, `pbwt_array_A`, `A_small_idx`, `cm_pos`, `counter_sel_gf`, `length_sel_gf`
	 * - And calls the member function `selectK()`.
	 */
	void init_common(const int k, const int l, const int ref_rac_l_com);

	/**
	 * @brief Initializes PBWT data structures and intervals for rare variants at a given site.
	 * 
	 * This function prepares the small PBWT arrays and updates interval bounds for rare
	 * target haplotypes associated with the variant site `k`. It also identifies the set
	 * of rare haplotypes to be considered from the current site onwards until the next common variant.
	 * 
	 * @param M Reference to the variant map containing variant metadata.
	 * @param k The index of the current variant site in the total set of variant sites.
	 * @param l The local index in the compressed PBWT structure corresponding to the current site.
	 * 
	 * @details
	 * - Copies relevant elements from the global PBWT array to the small PBWT arrays for the current site.
	 * - Gathers the list of rare target haplotypes spanning from site `k` until the next common variant site.
	 * - Sorts and removes duplicates from the list of rare haplotypes.
	 * - Updates PBWT interval boundaries (`f_k_small`, `g_k_small`) and index (`pbwt_small_index`) for each rare haplotype.
	 * - Resets intervals if they collapse (i.e., upper bound <= lower bound), marking the haplotype for reinitialization.
	 * 
	 * @pre
	 * - `A_small_idx[l]` must be a valid index mapping of small PBWT sites at index `l`.
	 * 
	 * @post
	 * - Small PBWT arrays `pbwt_small_A`, `pbwt_small_B`, `pbwt_small_V` are resized and updated.
	 * - The list `rareTarHaps` is populated with the rare haplotypes active at the current region.
	 * - Interval bounds and indices for rare haplotypes are updated accordingly.
	 */
	void init_rare(const variant_map & M,const int k, const int l);

	/**
	 * @brief Decode and update the full PBWT array for reference haplotypes at a given site.
	 * 
	 * This function reads encoded PBWT data from the compressed buffer pointed by `pY` 
	 * and updates the PBWT arrays `pbwt_array_A`, `pbwt_array_B`, and auxiliary array `pbwt_array_V`.
	 * It implements a run-length decoding combined with PBWT update logic to reconstruct the
	 * permutation of haplotypes at the current variant site.
	 * 
	 * @param pY Pointer reference to the compressed PBWT data buffer; this pointer is advanced.
	 * @param ref_rac_l The reference allele count at the current variant site, used in interval computations.
	 * 
	 * @note
	 * - `pbwt_array_A` and `pbwt_array_B` represent alternating states of the PBWT permutation arrays.
	 * - `pbwt_array_V` stores prefix sums or counts needed to update FM-index style intervals.
	 * 
	 * @see haplotype_set::read_small_pbwt_av for the rare-variant PBWT equivalent.
	 */
	void read_full_pbwt_av(const unsigned char*& pY, const int ref_rac_l);

	/**
	 * @brief Decode and update the small PBWT array for rare haplotypes at a given site.
	 * 
	 * Similar to `read_full_pbwt_av` but operates on the smaller subset of haplotypes
	 * representing rare variants, using `pbwt_small_A`, `pbwt_small_B`, and `pbwt_small_V`.
	 * It conditionally updates the auxiliary array `pbwt_small_V` based on `update_v`.
	 * 
	 * @param pY Pointer reference to the compressed PBWT data buffer; this pointer is advanced.
	 * @param ref_rac_l The reference allele count at the current variant site.
	 * @param update_v Whether to update the auxiliary vector `pbwt_small_V` (default true).
	 * 
	 * @see haplotype_set::read_full_pbwt_av for the full PBWT equivalent.
	 */
	void read_small_pbwt_av(const unsigned char*& pY, const int ref_rac_l, const bool update_v);
	
	/**
	 * @brief Select the top matching haplotypes around a target haplotype in the PBWT ordering.
	 * 
	 * This function selects up to `k0` haplotypes from the PBWT arrays around the index of the
	 * target haplotype `htr` for further analysis. It selects haplotypes both "up" (before) and 
	 * "down" (after) the current index within the search interval `[f_k[htr], g_k[htr])`.
	 * 
	 * The matching is filtered by allele `a` at the current variant site and restricted to 
	 * a genomic region defined by `ref_rac_l`.
	 * 
	 * @param htr Index of the target haplotype.
	 * @param k Current variant site index.
	 * @param ref_rac_l Reference allele count at site `k`.
	 * @param pbwt_array PBWT ordering array to select haplotypes from.
	 * @param k0 Number of haplotypes to select; defaults to member variable `K` if zero or negative.
	 * @param a Allele state (0 or 1) for filtering matching haplotypes.
	 * 
	 * @note
	 * - The function stores the selected haplotype indices in `pbwt_states[sample][o]`.
	 * - `sample` is the sample index corresponding to the haplotype `htr`.
	 */
	void selectK(const int htr, const int k, const int ref_rac_l, const std::vector < int >& pbwt_array, const int k0, const unsigned char a);
	
	/**
	 * @brief Selects up to k rare haplotypes around a target haplotype index from a small PBWT array.
	 * 
	 * This function operates similarly to `selectK` but uses the "small" PBWT arrays and indices
	 * specialized for rare haplotypes. It selects haplotypes near the current PBWT small index
	 * `pbwt_small_index[htr]` for a given haplotype `htr`, balancing the number of haplotypes selected
	 * above and below the index according to `k0` or the default `K`.
	 * 
	 * Only haplotypes passing the allele state filter `a` relative to the reference cutoff `ref_rac_l`
	 * are included. Selected haplotype indices are added to `pbwt_states` for the corresponding sample.
	 * 
	 * @param[in] htr Index of the target haplotype.
	 * @param[in] k Current PBWT position or site index.
	 * @param[in] ref_rac_l Reference allele count or cutoff to separate haplotypes.
	 * @param[in] pbwt_array The PBWT permutation array of haplotype indices for rare haplotypes.
	 * @param[in] k0 Optional limit for number of haplotypes to select (if <= 0, uses default K).
	 * @param[in] a Allele state (0 or 1) to filter haplotypes to select.
	 */
	void selectKrare(const int htr, const int k, const int ref_rac_l, const std::vector < int >& pbwt_array, const int k0, const unsigned char a);
	
	/**
	 * @brief Select top matching rare haplotypes around a target rare haplotype in the small PBWT ordering.
	 * 
	 * Similar to `selectK` but operates on the smaller PBWT arrays (`pbwt_small_A`) and
	 * indices related to rare variant haplotypes. Selects haplotypes around `pbwt_small_index[htr]`.
	 * 
	 * @param htr Index of the target haplotype.
	 * @param k Current variant site index.
	 * @param ref_rac_l Reference allele count at site `k`.
	 * @param pbwt_array PBWT ordering array for rare haplotypes.
	 * @param k0 Number of haplotypes to select; defaults to member variable `K` if zero or negative.
	 * @param a Allele state (0 or 1) for filtering matching haplotypes.
	 * 
	 * @see haplotype_set::selectK
	 */
	void select_common_pd_fg(const int k, const int l_hq, const int l_all, const int ref_rac_l, const int prev_ref_rac_l);
	
	/**
	 * @brief Update PBWT intervals for rare haplotypes at a given variant site.
	 * 
	 * Updates the PBWT search intervals and indices for rare target haplotypes carrying
	 * rare variants at site `k`. Resets intervals if they collapse and tracks resets.
	 * 
	 * @param k Index of the current variant site.
	 * @param ref_rac_l Reference allele count at site `k`.
	 * 
	 * @note
	 * - Uses `SvarTar[k]` to identify haplotypes with rare variants at site `k`.
	 * - Updates `pbwt_small_index`, `f_k_small`, `g_k_small`, and `last_reset` arrays.
	 * - Increments `counter_rare_restarts` when intervals are reset due to collapse.
	 */
	void select_rare_pd_fg(const int k, const int ref_rac_l);

	/**
	 * @brief Update the checksum with the current state of the haplotype set.
	 * 
	 * This function updates the provided checksum object with all relevant internal
	 * data members to ensure data integrity and consistency checking.
	 * 
	 * It calls the base class `ref_haplotype_set` checksum update, then processes
	 * haplotype counts, variant data, ploidy information, and sample-to-haplotype mappings.
	 * 
	 * @param crc Reference to a checksum object to update.
	 */
	void update_checksum(checksum &crc) const
	{
		ref_haplotype_set::update_checksum(crc);
		crc.process_data(n_tot_haps);
		crc.process_data(n_tar_haps);
		crc.process_data(n_tar_samples);
		HvarTar.update_checksum(crc);
		crc.process_data(ShapTar);
		crc.process_data(SvarTar);
		crc.process_data(SindTarGL);
		crc.process_data(cm_pos);
		crc.process_data(fploidy);
		crc.process_data(max_ploidy);
		crc.process_data(tar_ploidy);
		crc.process_data(tar_ind2hapid);
		crc.process_data(tar_hapid2ind);
	}
};

#endif
