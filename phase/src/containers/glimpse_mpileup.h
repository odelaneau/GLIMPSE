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

#ifndef SRC_CONTAINERS_GLIMPSE_MPILEUP_H_
#define SRC_CONTAINERS_GLIMPSE_MPILEUP_H_

#include <math.h>
#include <htslib/faidx.h>
#include <htslib/regidx.h>

#include "htslib/hts.h"

#include <utils/otools.h>
#include <containers/variant_map.h>

#define OFILE_VCFU	0
#define OFILE_VCFC	1
#define OFILE_BCFC	2

const int N_ALLELES = 4;
typedef int64_t hts_pos_t;

typedef struct {
    //char *ref[2];
    //int ref_id[2];
    //int ref_len[2];
	char *ref;
	int ref_id;
	int ref_len;
} ref_t;

struct aux_t {     			// auxiliary data structure
	samFile * fp;				// the file handle
	hts_idx_t *idx;
	bam_hdr_t * hdr;			// the file header
	hts_itr_t * iter;			// NULL if a region not specified
	int min_mapQ;				// mapQ threshold for filter
	//ref_t *ref;
	//faidx_t* fai;

	bool keep_orphan,check_orientation,check_proper_pair;	// parameters for filtering
	int fflag;												// Filter flag
	//unsigned int mtlen;										// Maximum observed fragment length
};

struct glimpse_bam_mplp_iter {
    int n;
    int32_t min_tid;
    std::vector<int32_t> tid;
    hts_pos_t min_pos;
    std::vector<int32_t> pos;
    std::vector<bam_plp_t> iter;
    std::vector<int> n_plp;
    std::vector<const bam_pileup1_t *> plp;
};

struct glimpse_bam_plp_iter {
    int32_t tid;
    int pos;
    bam_plp_t iter;
    const bam_pileup1_t * plp;
};


enum call_model { standard, pseudohaploid};

struct bcf_call_aux_t {
    int max_dp;
    int min_bq;
    std::vector<uint16_t> bases;  // 5bit: unused, 6:quality, 1:is_rev, 4:2-bit base or indel allele (index to bcf_callaux_t.indel_types)
};

struct bcf_call_ret1_t{
    //uint32_t ori_depth;     // ori_depth = anno[0..3] but before --min-BQ is applied
    std::array<float,16> p;        // phred-scaled likelihood of each genotype
};

struct sample_gl
{
	int read_depth;
	int dp_ind;
	int ref_read;
	std::array<float,3> llk;

};

struct call_t
{
	int ploidy;
	sample_gl gls;
	bcf_call_aux_t bca;
	bcf_call_ret1_t bcr;
	bam_plp_t s_plp;
	const bam_pileup1_t* v_plp;
	int n_plp;
	bool snp_called;
};

class glimpse_mpileup {
public:
	glimpse_mpileup();
	glimpse_mpileup(const glimpse_mpileup& );

	virtual ~glimpse_mpileup();

	std::vector<std::string> bam_fnames;
	std::vector<std::string> tar_sample_names;

	//samples
	int n_tar_samples;
	int n_tar_haps;
	std::vector<int> tar_ploidy;
	std::vector<int> tar_ind2gt;
	std::vector<int> tar_ind2pl;

	int n_tar_diploid;
	int n_tar_haploid;
	int max_ploidy;
	int fploidy;

	//fasta
	std::string fai_fname;
	faidx_t* fai;

	///filters
	bool keep_orphan,check_orientation,check_proper_pair;	// parameters for filtering
	int fflag;												// Filter flag
	bool illumina13;
	int min_mq;
	int min_bq;
	int max_dp;
};

#endif /* SRC_CONTAINERS_GLIMPSE_MPILEUP_H_ */
