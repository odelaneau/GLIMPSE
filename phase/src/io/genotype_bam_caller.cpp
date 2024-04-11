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

#include <io/genotype_bam_caller.h>

static int read_bam(void *data, bam1_t *b)
{
	aux_t * aux = (aux_t*) data;
	int ret, skip = 0;
	do{
		ret =  (aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b));
		if (ret < 0) break;
		else if (b->core.tid < 0 || (b->core.flag & aux->fflag) || b->core.qual < aux->min_mapQ) {skip = 1; continue;} // unmapped
		else if(b->core.flag & BAM_FPAIRED)
		{ //paired read
			if (b->core.flag & BAM_FMUNMAP) //mate is unmapped
			{ //mate unmapped
				if (!aux->keep_orphan) { skip = 1;continue;}
			}else{
				if (aux->check_orientation && ((b->core.tid != b->core.mtid) ||  //must be on the same chr
					((b->core.flag & BAM_FREVERSE) && (b->core.flag & BAM_FMREVERSE)) || // read 1 and 2 cannot be on the same strand
					((b->core.pos < b->core.mpos) && (b->core.flag & BAM_FREVERSE)) || //read 1 must be on the +ve strand
					((b->core.pos > b->core.mpos) && !(b->core.flag & BAM_FREVERSE)) //read 2 must be on the -ve strand
					))
				{
					skip = 1; continue;
				}
				if (aux->check_proper_pair && !(b->core.flag & BAM_FPROPER_PAIR))
				{
					skip = 1; continue;
				} //not properly paired
			}
		}
		skip = 0;
	} while(skip);
	return ret;
}

genotype_bam_caller::genotype_bam_caller(
		haplotype_set &_H,
		genotype_set & _G,
		const variant_map & _V,
		const glimpse_mpileup& _mpileup_data,
		const std::string _region,
		const uint32_t _beg,
		const uint32_t _end,
		std::string _callmodel,
		bool call_indels) :
				H(_H),G(_G),V(_V), mpileup_data(_mpileup_data), region(_region), beg(_beg),end(_end)
{
	caller.gls.read_depth = 0;
	caller.gls.dp_ind = 0;

	caller.bca.min_bq = mpileup_data.min_bq;
	caller.bca.max_dp = mpileup_data.max_dp;
	caller.bca.bases.reserve(((mpileup_data.max_dp + 31) & -32));

	//model
	group_reads = std::vector<std::function<void(call_t& caller, const variant* variant)>>(3);
	compute_llk = std::vector<std::function<void(call_t& caller, const variant* variant)>>(2);

	group_reads[0] = &bcf_call_glfgen_min;
	group_reads[1] = call_indels ? &bcf_call_glfgen_min_ins : &bcf_call_void;
	group_reads[2] = call_indels ? &bcf_call_glfgen_min_del : &bcf_call_void;

	if (_callmodel=="standard")
	{
		model = standard;

		compute_llk[0] = &standard_errmod;
		compute_llk[1] = call_indels ? &standard_errmod_indels_v : &flat_errmod_indels_v;
	}
	/*
	else if(_callmodel=="pseudohaploid")
	{
		model = pseudohaploid;

		compute_llk[0] = &pseudohaploid_errmod;
		compute_llk[1] = &pseudohaploid_errmod_indel;
	}
	*/
	else vrb.error("Model used is not supported: [" + _callmodel + "]");

	aux_data.check_orientation = mpileup_data.check_orientation;
	aux_data.check_proper_pair = mpileup_data.check_proper_pair;
	aux_data.fflag = mpileup_data.fflag;
	aux_data.keep_orphan = mpileup_data.keep_orphan;
	aux_data.min_mapQ = mpileup_data.min_mq;


	aux_data.iter=nullptr;
	aux_data.idx=nullptr;
	aux_data.hdr=nullptr;
	aux_data.fp=nullptr;
}

void genotype_bam_caller::clean()
{
	if ( aux_data.iter) hts_itr_destroy(aux_data.iter);
	aux_data.iter=nullptr;
	if (aux_data.idx) hts_idx_destroy(aux_data.idx);
	aux_data.idx=nullptr;
	if (aux_data.hdr) bam_hdr_destroy(aux_data.hdr);
	aux_data.hdr=nullptr;
	if (aux_data.fp)sam_close(aux_data.fp);
	aux_data.fp=nullptr;

	caller.gls.read_depth = 0;
	caller.gls.dp_ind = 0;
	std::fill(caller.gls.llk.begin(), caller.gls.llk.end(), 0.0f);

	//stats.cov_ind.clear();
	//std::fill(stats.depth_count.begin(), stats.depth_count.end(), 0);

	caller.ploidy=0;
}

genotype_bam_caller::~genotype_bam_caller()
{
	if ( aux_data.iter) hts_itr_destroy(aux_data.iter);
	aux_data.iter=nullptr;
	if (aux_data.idx) hts_idx_destroy(aux_data.idx);
	aux_data.idx=nullptr;
	if (aux_data.hdr) bam_hdr_destroy(aux_data.hdr);
	aux_data.hdr=nullptr;
	if (aux_data.fp)sam_close(aux_data.fp);
	aux_data.fp=nullptr;
}

void genotype_bam_caller::call_mpileup(int i)
{
	caller.ploidy = mpileup_data.tar_ploidy[i];

	aux_data.fp = sam_open(mpileup_data.bam_fnames[i].c_str(), "r");
	if (aux_data.fp == nullptr) vrb.error("Cannot open BAM/CRAM file: [" + mpileup_data.bam_fnames[i] + "]");

	if (hts_set_opt(aux_data.fp, CRAM_OPT_DECODE_MD, 0)) vrb.error("Failed to set CRAM_OPT_DECODE_MD value");
	if (!mpileup_data.fai_fname.empty()) if (hts_set_fai_filename(aux_data.fp, mpileup_data.fai_fname.c_str()) != 0)  vrb.error("Failed to process fasta");

	aux_data.hdr  = sam_hdr_read(aux_data.fp);
	if ( !aux_data.hdr ) vrb.error("Failed to read header: " + mpileup_data.bam_fnames[i] + ".");

	if (mpileup_data.bai_fnames.empty()) {
		aux_data.idx = sam_index_load2(aux_data.fp, mpileup_data.bam_fnames[i].c_str(), mpileup_data.bai_fnames[i].c_str());
	} else {
		aux_data.idx = sam_index_load(aux_data.fp, mpileup_data.bam_fnames[i].c_str());
	}
	if (aux_data.idx == NULL) vrb.error("Failed to load index for [" + mpileup_data.bam_fnames[i] + "]");
	aux_data.iter = sam_itr_querys(aux_data.idx, aux_data.hdr, region.c_str());
	if ( aux_data.iter==NULL  ) vrb.error("Failed to load region: [" + region + "] for file: [" + mpileup_data.bam_fnames[i] + "]");

	glimpse_mpileup_reg(i);
	clean();
}

int genotype_bam_caller::glimpse_mpileup_reg(int i)
{
    int tid=0, pos=0;
    const size_t n_sites_reg = V.vec_pos.size();
    const double n_sites_reg_d = n_sites_reg*1.0;

	size_t i_site = 0;
	stats1D coverage_reg;
	std::vector<unsigned char>& GLs  = G.vecG[i]->GL;
	const int sample_ploidy = mpileup_data.tar_ploidy[i];
	const int sample_ploidyP1 = mpileup_data.tar_ploidy[i]+1;

	caller.n_plp = 0;
	caller.s_plp = bam_plp_init(read_bam, (void*)&aux_data);
	unsigned long linecount = 0, prevstart = 0, prevend = 0;
	const bam_pileup1_t * v_plp;

	while (((caller.v_plp = bam_plp_auto(caller.s_plp, &tid, &pos, &caller.n_plp)) != 0) && i_site < n_sites_reg)
	{
		if (pos+1 < beg) continue;
		if (pos+1 > end) break;

		while (i_site < n_sites_reg && pos+1 > V.vec_pos[i_site]->bp)
		{
			++G.stats.depth_count[i][0];
			G.stats.cov_ind[i].push(0);
			++i_site;
		}
		caller.snp_called=false;
		while (i_site < n_sites_reg && pos+1 ==V.vec_pos[i_site]->bp)
		{

			const variant * variant_site = V.vec_pos[i_site];
			const bool is_indel = variant_site->type == VCF_INDEL;
			const bool is_delet = is_indel ? variant_site->alt.size() == 1 : 0;

			group_reads[is_indel + is_delet] (caller, variant_site);
			compute_llk[is_indel] (caller, variant_site);
			if (caller.gls.dp_ind > 0)
			{
				int max_llk = 0;
				for (int j=1;j<sample_ploidyP1;++j) if (caller.gls.llk[j] > caller.gls.llk[max_llk]) max_llk=j;
				for (int j=0; j<sample_ploidyP1; ++j) GLs[sample_ploidyP1 * i_site+j] = (unsigned char) std::min(static_cast<int>(-10*(caller.gls.llk[j]- max_llk)),255);
				if (!(GLs[sample_ploidyP1 * i_site+0] == GLs[sample_ploidyP1 * i_site+1] && GLs[sample_ploidyP1 * i_site+0] == GLs[sample_ploidyP1 * i_site+sample_ploidy]))
					G.vecG[i]->flat[i_site] = false;
			}
			++G.stats.depth_count[i][caller.gls.dp_ind];
			G.stats.cov_ind[i].push(caller.gls.dp_ind);

			caller.snp_called = !is_indel;
			++i_site;
		}
	}
	while (i_site < n_sites_reg)
	{
		++G.stats.depth_count[i][0];
		G.stats.cov_ind[i].push(0);
		++i_site;
	}

	caller.gls.dp_ind =0;

	bam_plp_reset(caller.s_plp);
	bam_plp_destroy(caller.s_plp);
    return 0;
}


void bcf_call_glfgen_min(call_t& caller, const variant* variant)
{
	if (caller.snp_called) return;

	const int _n = caller.n_plp;
	const bam_pileup1_t *pl = caller.v_plp;

	caller.bca.bases.clear();
	caller.gls.read_depth = 0;
	if (_n == 0) return;

	for (int i = 0; i < _n; ++i) {
		const bam_pileup1_t *p = pl + i;
		if (p->is_refskip) continue;
		if (p->is_del) continue;

		int b = seq_nt16_int[bam_seqi(bam_get_seq(p->b), p->qpos)]; // base
		int q = (int)bam_get_qual(p->b)[p->qpos];
		if (b > 3 || q < caller.bca.min_bq) continue;
		caller.bca.bases.push_back(q<<5 | (int)bam_is_rev(p->b)<<4 | b);
	}
	caller.gls.read_depth = caller.bca.bases.size();
}

void bcf_call_glfgen_min_ins(call_t& caller, const variant* variant)
{
	const int _n = caller.n_plp;
	const bam_pileup1_t *pl = caller.v_plp;

	caller.bca.bases.clear();
	caller.gls.read_depth = 0;
	if (_n == 0) return;

    const std::string& ref = variant->ref;
    const std::string& alt = variant->alt;

    const int ref_length = ref.size();
    const int indel_length = alt.size();
    if (ref_length > 1) return;
    if (indel_length < 2) return;

    for (int i = 0; i < _n; ++i)
	{
		const bam_pileup1_t *p = pl + i;
		int k=1;

		if (p->is_refskip || p->is_del) continue;

		if (seq_nt16_int[seq_nt16_table[ref[0]]] != seq_nt16_int[bam_seqi(bam_get_seq(p->b), p->qpos)]) continue;
		int q = (int)bam_get_qual(p->b)[p->qpos];

		if (p->indel)
		{
			//if (p->indel > 0) j=1;
			//else continue;
			const int read_length = p->b->core.l_qseq;
			const int length1 = std::min(indel_length, read_length - p->qpos);
			if (p->indel == length1-1) //the read has an indel of the right size
			{
				uint8_t *seq = bam_get_seq(p->b);
				for (k = 1; k < p->indel && seq_nt16_int[bam_seqi(seq, p->qpos + k)] == seq_nt16_int[seq_nt16_table[ref[k]]]; ++k) //we look for full matching, partial matchings will reinforce REF
						q+= (int)bam_get_qual(p->b)[p->qpos + k];

				if (k != p->indel) continue;
				q/=p->indel;
			} else continue;//indel but not supporting our indel


		}
		if (q < caller.bca.min_bq) continue;
		caller.bca.bases.push_back(q<<5 | (int)bam_is_rev(p->b)<<4 | p->indel!=0); //different coding, j=0 for ref, j=1 for indel match
	}
    caller.gls.read_depth = caller.bca.bases.size();
}

void bcf_call_void(call_t& caller, const variant* variant)
{
	caller.bca.bases.clear();
	caller.gls.read_depth = 0;
}

void bcf_call_glfgen_min_del(call_t& caller, const variant* variant)
{
	const int _n = caller.n_plp;
	const bam_pileup1_t *pl = caller.v_plp;

	caller.bca.bases.clear();
	caller.gls.read_depth = 0;
	if (_n == 0) return;

    const std::string& ref = variant->ref;
    const std::string& alt = variant->alt;

    const int ref_length = ref.size();
    const int del_length = -ref_length-1;
    if (ref_length <= 1) return;
    if (alt.size() != 1) return;

    for (int i = 0; i < _n; ++i)
	{
		const bam_pileup1_t *p = pl + i;
		int k=1;

		if (p->is_refskip || p->is_del) continue;

		if (seq_nt16_int[seq_nt16_table[ref[0]]] != seq_nt16_int[bam_seqi(bam_get_seq(p->b), p->qpos)]) continue;
		int q = (int)bam_get_qual(p->b)[p->qpos];

		if (p->indel == 0)
		{
			const int read_length = p->b->core.l_qseq;
			uint8_t *seq = bam_get_seq(p->b);
 			const int length = std::min(ref_length, read_length - p->qpos);
			for (k = 1; k < length && seq_nt16_int[bam_seqi(seq, p->qpos + k)] == seq_nt16_int[seq_nt16_table[ref[k]]]; ++k) //we look for full matching, partial matchings will reinforce REF
					q+= (int)bam_get_qual(p->b)[p->qpos + k];

			if (k != length) continue;
			q/=length;
		}
		else if (p->indel > del_length) continue;

		if (q < caller.bca.min_bq) continue;
		caller.bca.bases.push_back(q<<5 | (int)bam_is_rev(p->b)<<4 | p->indel!=0); //different coding, j=0 for ref, j=1 for indel match
	}
	//r.ori_depth=bca.bases.size();
    caller.gls.read_depth = caller.bca.bases.size();
}

void standard_errmod(call_t& caller, const variant* variant)
{
    int ref = seq_nt16_int[seq_nt16_table[variant->ref[0]]];
	int alt = seq_nt16_int[seq_nt16_table[variant->alt[0]]];

	if (ref >= N_ALLELES|| alt >= N_ALLELES) vrb.error("ERROR base sequences.");

	const int sploidy = caller.ploidy;

	caller.gls.dp_ind = 0;
	if (caller.gls.read_depth == 0) return;

	if (!caller.snp_called)
	{
		int ii,k,j=0;
	    if (caller.gls.read_depth > caller.bca.max_dp)
	    {
	    	// if we exceed bcf.maxbases (255) bases observed then shuffle them to sample and only keep the first bca.max_bases
	    	shuffle(caller.bca.bases.begin(), caller.bca.bases.end(), rng.randomEngine);
	    	caller.gls.read_depth = caller.bca.max_dp;
	    }
	    std::fill(caller.bcr.p.begin(), caller.bcr.p.end(), 0.0f);
	    std::sort(caller.bca.bases.begin(), caller.bca.bases.end());

		while (j<caller.gls.read_depth)
		{
	    	uint16_t b = caller.bca.bases[j];

	    	const int n_identical_reads =
	    			std::distance(caller.bca.bases.begin()+j, std::upper_bound(caller.bca.bases.begin()+j, caller.bca.bases.begin()+caller.gls.read_depth, b));

	        // extract quality
	        const int qual = b>>5 < 4? 4 : b>>5; //TODO check
	        // extract base
	        const int base = b&0xf;

	        if (sploidy > 1)
	        {
				const float errM1 = calling_utils::log10_one_minus_err[qual];
				const float errDivThree = calling_utils::log10_err_div_three[qual];
				const float oneHalfMinusErrDivThree= calling_utils::log10_one_div_two_minus_err_div_three[qual];

				caller.bcr.p[base*N_ALLELES+base] += errM1*n_identical_reads; //HOM_BB

				for (k = 0; k < base; ++k)
					caller.bcr.p[k*N_ALLELES+base] = caller.bcr.p[base*N_ALLELES+k] += oneHalfMinusErrDivThree*n_identical_reads; //HET_B

				for (k = base+1; k < N_ALLELES; ++k)
					caller.bcr.p[k*N_ALLELES+base] = caller.bcr.p[base*N_ALLELES+k] += oneHalfMinusErrDivThree*n_identical_reads; //HET_B

				for (ii = 0; ii < N_ALLELES; ++ii)
				{
					if (ii==base) continue;

					for (k=ii; k<N_ALLELES;++k)
					{
						if (k==base) continue;
						caller.bcr.p[k*N_ALLELES+ii] = caller.bcr.p[ii*N_ALLELES+k] += errDivThree*n_identical_reads;
					}
				}
	        }
	        else
	        {
	        	const float one_minus_err = calling_utils::log10_one_minus_err[qual];
	        	const float err = calling_utils::log10_err[qual];

	        	caller.bcr.p[base*N_ALLELES+base] += one_minus_err;

				for (ii = 0; ii < N_ALLELES; ++ii)
				{
					if (ii==base) continue;
					caller.bcr.p[ii*N_ALLELES+ii] += err;
				}
	        }
	        caller.gls.dp_ind+=n_identical_reads;
	        if (base==ref) caller.gls.ref_read += n_identical_reads;
			j+=n_identical_reads;
		}
	}
    if (sploidy > 1)
    {
    	caller.gls.llk[0] = caller.bcr.p[ref*N_ALLELES+ref];
    	caller.gls.llk[1] = caller.bcr.p[ref*N_ALLELES+alt];
    	caller.gls.llk[2] = caller.bcr.p[alt*N_ALLELES+alt];
    }
    else
    {
    	caller.gls.llk[0] = caller.bcr.p[ref*N_ALLELES+ref];
    	caller.gls.llk[1] = caller.bcr.p[alt*N_ALLELES+alt];
    	caller.gls.llk[2] = -255.5f;
    }

}

void pseudohaploid_errmod(call_t& caller, const variant* variant)
{
	int ref = seq_nt16_int[seq_nt16_table[variant->ref[0]]];
	int alt = seq_nt16_int[seq_nt16_table[variant->alt[0]]];

	//int n = bca.bases.size();
	caller.gls.dp_ind = 0;
	caller.gls.ref_read = 0;

	if (caller.gls.read_depth == 0) return;

	const int sploidy = caller.ploidy;

	// extract base
	if (!caller.snp_called)
	{
		const int base = (caller.bca.bases[rng.getInt(caller.gls.read_depth)])&0xf;

		caller.bcr.p[base*N_ALLELES+base] = 0.0f;
		for (int ii = 0; ii < N_ALLELES; ++ii)
		{
			if (ii==base) continue;
			caller.bcr.p[ii*N_ALLELES+ii] = -25.5f;
		}
	}
	if (caller.bcr.p[ref*N_ALLELES+ref] <0.0f && caller.bcr.p[alt*N_ALLELES+alt] < 0.0f) return;

	caller.gls.dp_ind = 1;
	if (caller.bcr.p[ref*N_ALLELES+ref]<1e-7) caller.gls.ref_read = 1;
	if (sploidy > 1)
	{
		caller.gls.llk[0] = caller.bcr.p[ref*N_ALLELES+ref];
		caller.gls.llk[2] = caller.bcr.p[alt*N_ALLELES+alt];
		caller.gls.llk[1] = -25.5f;
	}
	else
	{
		caller.gls.llk[0] = caller.bcr.p[ref*N_ALLELES+ref];
		caller.gls.llk[1] = caller.bcr.p[alt*N_ALLELES+alt];
		caller.gls.llk[2] = -25.5f;
	}
}

void pseudohaploid_errmod_indel(call_t& caller, const variant* variant)
{
	caller.gls.dp_ind = 0;
	caller.gls.ref_read = 0;
	if (caller.gls.read_depth == 0) return;

	const int sploidy = caller.ploidy;

	//sample a read
	uint16_t b = caller.bca.bases[rng.getInt(caller.gls.read_depth)];
	caller.gls.dp_ind=1;

	// extract base
	const int base = b&0xf;
	if (base==0) caller.gls.ref_read = 1;

	for (int j=0;j<3; ++j) caller.gls.llk[j] = -25.5f;
	caller.gls.llk[base*sploidy] = 0.0f;
}

void flat_errmod_indels_v(call_t& caller, const variant* variant)
{
	caller.gls.dp_ind=0;
    caller.gls.ref_read = 0;
    caller.gls.read_depth = 0;
    caller.gls.llk[0] = 0.0f;
    caller.gls.llk[1] = 0.0f;
    caller.gls.llk[2] = 0.0f;
}

void standard_errmod_indels_v(call_t& caller, const variant* variant)
{
    int j=0;
    caller.gls.dp_ind = 0;
    caller.gls.ref_read = 0;
	if (caller.gls.read_depth == 0) return;

	const int sploidy = caller.ploidy;

    if (caller.gls.read_depth > caller.bca.max_dp) {
    	shuffle(caller.bca.bases.begin(), caller.bca.bases.end(), rng.randomEngine);
    	caller.gls.read_depth = caller.bca.max_dp;
    }
    std::sort(caller.bca.bases.begin(), caller.bca.bases.end());

    caller.gls.llk[0] = 0.0f;
    caller.gls.llk[1] = 0.0f;
    caller.gls.llk[2] = 0.0f;

	while (j<caller.gls.read_depth)
	{
    	uint16_t b = caller.bca.bases[j];

    	const int n_identical_reads =
    			std::distance(caller.bca.bases.begin()+j, std::upper_bound(caller.bca.bases.begin()+j, caller.bca.bases.begin()+caller.gls.read_depth, b));
        const int qual = b>>5 < 4? 4 : b>>5;
        const int base = b&0xf; //zero ref, 1 alt

        if (sploidy > 1)
		{
            const float errM1 = calling_utils::log10_one_minus_err[qual];
    		const float errDivThree = calling_utils::log10_err_div_three[qual];
    		const float oneHalfMinusErrDivThree= calling_utils::log10_one_div_two_minus_err_div_three[qual];
    		caller.gls.llk[2*base] += errM1*n_identical_reads; //HOM_BB
    		caller.gls.llk[1] += oneHalfMinusErrDivThree*n_identical_reads; //HET_B;
    		caller.gls.llk[2-2*base] += errDivThree*n_identical_reads;
		}
        else
        {
        	const float one_minus_err = calling_utils::log10_one_minus_err[qual];
			const float err = calling_utils::log10_err[qual];

			caller.gls.llk[sploidy*base] += one_minus_err;
			caller.gls.llk[sploidy*(sploidy-base)] += err;
        }

        caller.gls.dp_ind+=n_identical_reads;
        if (base==0) caller.gls.ref_read += n_identical_reads;
		j+=n_identical_reads;
	}

    if (sploidy <= 1) caller.gls.llk[2] = -25.5f;
}
