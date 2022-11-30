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

#include "glimpse_mpileup.h"

glimpse_mpileup::glimpse_mpileup() :
	n_tar_samples(0), n_tar_haps(0),
	n_tar_diploid(0),n_tar_haploid(0), max_ploidy(0),fploidy(0),
	fai(NULL),
	keep_orphan(false),check_orientation(false),check_proper_pair(true),fflag(0),illumina13(false),min_mq(0),min_bq(0),max_dp(255)
{
}

//copy constructor
glimpse_mpileup::glimpse_mpileup(const glimpse_mpileup& mpileup_data) :
	n_tar_samples(mpileup_data.n_tar_samples), n_tar_haps(mpileup_data.n_tar_haps),
	n_tar_diploid(mpileup_data.n_tar_diploid),n_tar_haploid(mpileup_data.n_tar_haploid), max_ploidy(mpileup_data.max_ploidy),fploidy(mpileup_data.fploidy),
	fai(mpileup_data.fai),
	keep_orphan(mpileup_data.keep_orphan),check_orientation(mpileup_data.check_orientation),check_proper_pair(mpileup_data.check_proper_pair),fflag(mpileup_data.fflag),illumina13(mpileup_data.illumina13),min_mq(mpileup_data.min_mq),min_bq(mpileup_data.min_bq),max_dp(mpileup_data.max_dp)
{
	tar_ploidy = mpileup_data.tar_ploidy;
	tar_ind2gt = mpileup_data.tar_ind2gt;
	tar_ind2pl = mpileup_data.tar_ind2pl;
	fai_fname = mpileup_data.fai_fname;
}

glimpse_mpileup::~glimpse_mpileup()
{
	/*
	for (int i = 0; i < iter.n; ++i)
	{
		bam_plp_destroy(iter.iter[i]);
		iter.iter[i]=NULL;

		if ( aux_data[i].iter) hts_itr_destroy(aux_data[i].iter);
		aux_data[i].iter=NULL;

		hts_idx_destroy(aux_data[i].idx);
		aux_data[i].idx=NULL;
		sam_hdr_destroy(aux_data[i].hdr);
		aux_data[i].hdr=NULL;
		sam_close(aux_data[i].fp);
		aux_data[i].fp=NULL;

	}
	*/
	if (fai) fai_destroy(fai);
}

