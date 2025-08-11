/*******************************************************************************
 * @file glimpse_mpileup.h
 * @brief Definitions and structures for GLIMPSE mpileup processing.
 *
 * This header contains data structures, constants, and class declarations
 * used for mpileup processing in GLIMPSE. It includes utilities for
 * managing BAM/CRAM inputs, pileup iteration, genotype likelihoods, 
 * and filtering options.
 *
 * The mpileup functionality is used for variant calling and genotype
 * likelihood computation from sequencing data.
 *
 * @details
 * **Main components:**
 * - **Constants and Macros** for output file types and allele counts.
 * - **Reference sequence container** (`ref_t`).
 * - **BAM auxiliary structures** for pileup iteration (`aux_t`, `glimpse_bam_mplp_iter`, `glimpse_bam_plp_iter`).
 * - **Genotype likelihood models** and their supporting structures.
 * - **Main mpileup class** (`glimpse_mpileup`) containing input/output parameters and filter settings.
 *
 ******************************************************************************/

#ifndef SRC_CONTAINERS_GLIMPSE_MPILEUP_H_
#define SRC_CONTAINERS_GLIMPSE_MPILEUP_H_

#include <math.h>
#include <htslib/faidx.h>
#include <htslib/regidx.h>
#include "htslib/hts.h"

#include "otools.h"
#include "variant_map.h"

/**
 * @name Output File Type Constants
 * @{
 */
#define OFILE_VCFU  0   /**< Uncompressed VCF output */
#define OFILE_VCFC  1   /**< Compressed VCF output (.vcf.gz) */
#define OFILE_BCFC  2   /**< Compressed BCF output (.bcf) */
/** @} */

/**
 * @brief Number of alleles supported.
 */
const int N_ALLELES = 4;

/**
 * @brief HTSlib position type definition (64-bit).
 */
typedef int64_t hts_pos_t;

/**
 * @brief Container for reference sequence metadata.
 */
typedef struct {
    char *ref;      /**< Pointer to reference sequence string. */
    int ref_id;     /**< Reference sequence ID (e.g., chromosome index). */
    int ref_len;    /**< Length of the reference sequence. */
} ref_t;

/**
 * @brief Auxiliary data structure for BAM/CRAM input.
 */
struct aux_t {
    samFile *fp;        /**< File handle for BAM/CRAM. */
    hts_idx_t *idx;     /**< BAM/CRAM index. */
    bam_hdr_t *hdr;     /**< File header. */
    hts_itr_t *iter;    /**< Iterator for region queries; NULL if region not specified. */
    int min_mapQ;       /**< Minimum mapping quality filter threshold. */

    bool keep_orphan;       /**< Keep orphan reads flag. */
    bool check_orientation; /**< Enforce correct read orientation flag. */
    bool check_proper_pair; /**< Enforce proper pairing flag. */
    int fflag;              /**< Filter flag for reads. */
};

/**
 * @brief Multi-sample pileup iterator container.
 */
struct glimpse_bam_mplp_iter {
    int n;                                  /**< Number of samples. */
    int32_t min_tid;                        /**< Minimum reference ID observed. */
    std::vector<int32_t> tid;               /**< Per-sample reference IDs. */
    hts_pos_t min_pos;                      /**< Minimum genomic position. */
    std::vector<int32_t> pos;               /**< Positions for each sample. */
    std::vector<bam_plp_t> iter;            /**< Pileup iterators. */
    std::vector<int> n_plp;                 /**< Number of pileup entries per sample. */
    std::vector<const bam_pileup1_t *> plp; /**< Pileup entries per sample. */
};

/**
 * @brief Single-sample pileup iterator container.
 */
struct glimpse_bam_plp_iter {
    int32_t tid;                            /**< Reference ID. */
    int pos;                                /**< Genomic position. */
    bam_plp_t iter;                         /**< Pileup iterator. */
    const bam_pileup1_t *plp;               /**< Pointer to pileup entry. */
};

/**
 * @brief Genotype calling model types.
 */
enum call_model {
    standard,       /**< Standard diploid model. */
    pseudohaploid   /**< Pseudohaploid model. */
};

/**
 * @brief Auxiliary structure for genotype likelihood calling.
 */
struct bcf_call_aux_t {
    int max_dp;                         /**< Maximum depth. */
    int min_bq;                         /**< Minimum base quality. */
    std::vector<uint16_t> bases;        /**< Encoded base calls (5 bits per base). */
};

/**
 * @brief Genotype likelihood result structure.
 */
struct bcf_call_ret1_t {
    std::array<float,16> p; /**< Phred-scaled likelihood for each genotype. */
};

/**
 * @brief Per-sample genotype likelihoods and read counts.
 */
struct sample_gl {
    int read_depth;                 /**< Total read depth. */
    int dp_ind;                     /**< Depth index. */
    int ref_read;                   /**< Number of reference-supporting reads. */
    std::array<float,3> llk;        /**< Log-likelihoods for genotypes. */
};

/**
 * @brief Genotype call container for a single variant/sample.
 */
struct call_t {
    int ploidy;                         /**< Sample ploidy (1 = haploid, 2 = diploid). */
    sample_gl gls;                      /**< Genotype likelihoods and depths. */
    bcf_call_aux_t bca;                 /**< Auxiliary likelihood data. */
    bcf_call_ret1_t bcr;                 /**< Likelihood results. */
    bam_plp_t s_plp;                    /**< Pileup iterator. */
    const bam_pileup1_t* v_plp;         /**< Variant pileup entry. */
    int n_plp;                          /**< Number of pileup entries. */
    bool snp_called;                    /**< Whether a SNP was called. */
};

/**
 * @brief Main class for mpileup processing in GLIMPSE.
 */
class glimpse_mpileup {
public:
    /**
     * @brief Default constructor.
     */
    glimpse_mpileup();

    /**
     * @brief Copy constructor.
     * @param other Instance to copy.
     */
    glimpse_mpileup(const glimpse_mpileup& other);

    /**
     * @brief Destructor.
     */
    virtual ~glimpse_mpileup();

    // Input BAM files and sample names
    std::vector<std::string> bam_fnames;       /**< BAM/CRAM filenames. */
    std::vector<std::string> tar_sample_names; /**< Target sample names. */

    // Sample-related data
    int n_tar_samples;                         /**< Number of target samples. */
    int n_tar_haps;                            /**< Number of haplotypes. */
    std::vector<int> tar_ploidy;               /**< Ploidy per sample. */
    std::vector<int> tar_ind2gt;               /**< Sample-to-genotype index mapping. */
    std::vector<int> tar_ind2pl;               /**< Sample-to-PL index mapping. */

    int n_tar_diploid;                         /**< Number of diploid samples. */
    int n_tar_haploid;                         /**< Number of haploid samples. */
    int max_ploidy;                            /**< Maximum ploidy among samples. */
    int fploidy;                               /**< Forced ploidy (if set). */

    // Reference FASTA
    std::string fai_fname;                     /**< FASTA index filename. */
    faidx_t* fai;                              /**< FASTA index pointer. */

    // Filtering parameters
    bool keep_orphan;                          /**< Keep orphan reads flag. */
    bool check_orientation;                    /**< Enforce correct read orientation. */
    bool check_proper_pair;                    /**< Enforce proper pairing. */
    int fflag;                                 /**< Filter flag for reads. */
    bool illumina13;                           /**< Illumina 1.3+ base quality encoding. */
    int min_mq;                                /**< Minimum mapping quality. */
    int min_bq;                                /**< Minimum base quality. */
    int max_dp;                                /**< Maximum depth. */
};

#endif /* SRC_CONTAINERS_GLIMPSE_MPILEUP_H_ */
