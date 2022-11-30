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


#include <io/ref_genotype_reader.h>

#ifdef __XSI__
#include "c_api.h"
#endif

std::map<std::string, int> mapPloidy = {
		{"1",1},
		{"2",2}
};

ref_genotype_reader::ref_genotype_reader(
		ref_haplotype_set & _H,
		variant_map & _V,
		const std::string _region,
		const float _sparse_maf,
		const bool _keep_mono) :
				H(_H), V(_V), region(_region),
				sparse_maf(_sparse_maf), n_ref_samples(0), keep_mono(_keep_mono)
{
	H.sparse_maf = _sparse_maf;
}

ref_genotype_reader::~ref_genotype_reader()
{
}

void ref_genotype_reader::set_ploidy_ref(const int ngt_ref, const int* gt_arr_ref, const int ngt_arr_ref)
{
	int max_ploidy_ref = 0;
	int n_ref_haploid = 0;

	if (gt_arr_ref == nullptr || ngt_ref==0) vrb.error("Error setting ploidy while reading the first record.");

	max_ploidy_ref = ngt_ref/n_ref_samples;
	if (max_ploidy_ref < 1 || max_ploidy_ref > 2 || ngt_ref%n_ref_samples != 0) vrb.error("Max ploidy of the reference panel cannot be set to neither 1 or 2.");

	ploidy_ref_samples = std::vector<int> (n_ref_samples,2);

	if (max_ploidy_ref == 1)
	{
		n_ref_haploid=n_ref_samples;
		for(int i = 0 ; i < n_ref_samples ; ++i) ploidy_ref_samples[i] = 1;
	}
	else
	{
		for(int i = 0 ; i < n_ref_samples ; ++i)
		{
			if (gt_arr_ref[max_ploidy_ref*i+1] == bcf_int32_vector_end)
			{
				ploidy_ref_samples[i]=1;
				++n_ref_haploid;
			}
		}
	}
	H.n_ref_haps = n_ref_haploid + 2*(n_ref_samples-n_ref_haploid);

	//Allocate some parts of haplotype set
	H.allocate();

	//Verbose
	vrb.bullet("VCF/BCF scanning [Nr=" + stb.str(n_ref_samples) + " (Nrh=" + stb.str(n_ref_haploid) + " - Nrd=" + stb.str(n_ref_samples-n_ref_haploid) + ")" + "]");
	vrb.bullet("VCF/BCF scanning [Reg=" + V.input_gregion + " / L=" + stb.str(H.n_tot_sites) + "] [" + stb.str(tac.rel_time()*1.0/1000, 2) + "s]");
	vrb.wait("  * Reference panel parsing");
	tac.clock();
}


void ref_genotype_reader::initReader(bcf_srs_t *sr, const std::string fref, int nthreads)
{
	sr->require_index = 1;
	if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + region + "] in [" + fref + "]");
	if (bcf_sr_set_targets(sr, region.c_str(), 0, 0) == -1) vrb.error("Impossible to set target region [" + region + "] in [" + fref + "]");

	if (nthreads>1) bcf_sr_set_threads(sr, nthreads);

	//Opening files
	if(!(bcf_sr_add_reader (sr, fref.c_str())))
	{
		//we do not build an index here, as the target and reference panel could be accessed in parallel
		if (sr->errnum != idx_load_failed) vrb.error("Failed to open file: " + fref + "");
		else vrb.error("Failed to load index of the file: " + fref + "");
	}
}

void ref_genotype_reader::readRefPanel(std::string fref,int nthreads)
{
	if (stb.get_extension(fref) == "xsi") fref += "_var.bcf";

	bcf_srs_t * sr_scan =  bcf_sr_init();
	initReader(sr_scan, fref, nthreads);

	scanGenotypesCommon(sr_scan,fref, 0);
	bcf_sr_destroy(sr_scan);

	bcf_srs_t * sr_parse =  bcf_sr_init();
	initReader(sr_parse, fref, nthreads);
	parseRefGenotypes(sr_parse,fref);
	bcf_sr_destroy(sr_parse);
}

void ref_genotype_reader::scanGenotypesCommon(bcf_srs_t * sr,const std::string fref, int ref_sr_n /* Reference reader number*/) {
	if(sr==nullptr) vrb.error("Error reading input files");

	vrb.wait("  * VCF/BCF scanning");
	tac.clock();

	//Scan Ref file
	#ifdef __XSI__
		n_ref_samples = c_xcf_nsamples(fref.c_str());
	#else
		n_ref_samples = bcf_hdr_nsamples(sr->readers[ref_sr_n].header);
	#endif
	//n_ref_samples = bcf_hdr_nsamples(sr->readers[ref_sr_n].header);
	H.n_tot_sites = 0;

	//Scan file
	int nset;
	bcf1_t * line_ref = NULL;
	int rAC=0, nAC=0, *vAC = NULL;
	int rAN=0, nAN=0, *vAN = NULL;
	unsigned int cref = 0, calt = 0;
	int n_inc_sites = 0;
	int n_mis_sites = 0;

	sr->max_unpack = BCF_UN_INFO;
	H.n_com_sites_hq=0;
	bool warning_out=false;
	int prev_pos=-1;

	while ((nset = bcf_sr_next_line (sr)))
	{
		if (bcf_sr_get_line(sr, ref_sr_n)->n_allele != 2) continue; //always ref
		line_ref =  bcf_sr_get_line(sr, ref_sr_n);

		//Get AC / AN
		rAC = bcf_get_info_int32(sr->readers[ref_sr_n].header, line_ref, "AC", &vAC, &nAC);
		rAN = bcf_get_info_int32(sr->readers[ref_sr_n].header, line_ref, "AN", &vAN, &nAN);
		if ((nAC!=1)||(nAN!=1)) vrb.error("VCF for reference panel needs AC/AN INFO fields to be present");
		calt = vAC[0]; cref = (vAN[0] - vAC[0]);

		if (std::min(calt, cref) == 0 && !keep_mono) { if (!warning_out) vrb.warning("Monomorphic site found [AC field] in reference panel at position: " + std::to_string(line_ref->pos + 1) + ". ALL monomorphic variants will be skipped. Please check your reference panel file. Use the --keep-monomorphic-ref-sites option to force GLIMPSE to use monomorphic sites in the reference panel. This warning is shown only once."); warning_out=true; continue;}

		//Classify variant
		double MAF = std::min(cref*1.0f/(cref+calt), calt*1.0f/(cref+calt));
		const bool is_common = MAF >= sparse_maf;
		H.flag_common.push_back(is_common);
		H.major_alleles.push_back(calt > cref);
		H.n_com_sites += is_common;
		H.n_rar_sites += !is_common;
		if (is_common) H.common2tot.push_back(H.n_tot_sites);//build sparse idx mapping
		H.n_tot_sites ++;

		int8_t line_type = bcf_get_variant_types(line_ref);
		H.n_com_sites_hq += is_common && line_type==VCF_SNP && (line_ref->pos != prev_pos);
		//Push new variant into map
		V.push(new variant (line_ref->pos + 1, std::string(line_ref->d.id), line_ref->d.allele[0], line_ref->d.allele[1], line_type, V.size(), cref, calt, line_type==VCF_SNP && (line_ref->pos != prev_pos)));
		prev_pos=line_ref->pos;
	}
	free(vAC);
	free(vAN);
	if (H.n_tot_sites == 0) vrb.error("No variants to be imputed in files");

	if (ref_sr_n > 0 /* warning relies on reader number */) {
		n_mis_sites=H.n_tot_sites-n_inc_sites;
		vrb.bullet("VCF/BCF scanning [Reg=" + region + " / L=" + stb.str(H.n_tot_sites) + " (Li= " + stb.str(n_inc_sites) + " (" + stb.str(((float)n_inc_sites/H.n_tot_sites)*100.0) + "%) - Lm= " + stb.str(n_mis_sites)  + " (" + stb.str(((float)n_mis_sites/H.n_tot_sites)*100.0) + "%))]");
		if (n_mis_sites > 0) vrb.warning("There are variants of the reference panel with no genotype likelihood in the target dataset. Please compute the likelihoods at all variants in the reference panel in order to use all the information of the sequencing reads.");
	}

	vrb.bullet("VCF/BCF scanning [L=" + stb.str(H.n_tot_sites) + " (Lrare= " + stb.str(H.n_rar_sites) + " (" + stb.str(((float)H.n_rar_sites/H.n_tot_sites)*100.0, 1) + "%) - Lcommon= " + stb.str(H.n_com_sites) + " - hq= " + stb.str(H.n_com_sites)  + " (" + stb.str(((float)H.n_com_sites/H.n_tot_sites)*100.0, 1) + "%))]");
}

void ref_genotype_reader::parseRefGenotypes(bcf_srs_t * sr,const std::string fref) {
	if(sr==nullptr) vrb.error("Error reading reference file");

	unsigned int i_site = 0, i_common =0, n_ref_unphased = 0;
	int ngt_ref, *gt_arr_ref = NULL, ngt_arr_ref = 0;
	bcf1_t * line_ref;
	float prog_step = 1.0/H.n_tot_sites;
	float prog_bar = 0.0;
	sr->max_unpack=BCF_UN_ALL;

	unsigned int cref = 0, calt = 0;
	int line_max_ploidy=0;
	int idx_ref_hap=0;
	bool a = false;
	int *ptr;
	float *ptr_f;

	#ifdef __XSI__
		c_xcf* c_xcf_p = c_xcf_new();
		c_xcf_add_readers(c_xcf_p, sr);
	#endif

	int rAC=0, nAC=0, *vAC = NULL;
	int rAN=0, nAN=0, *vAN = NULL;
	while (bcf_sr_next_line (sr))
	{
		if (bcf_sr_get_line(sr, 0)->n_allele != 2) continue; //always ref

		line_ref =  bcf_sr_get_line(sr, 0);

		rAC = bcf_get_info_int32(sr->readers[0].header, line_ref, "AC", &vAC, &nAC);
		rAN = bcf_get_info_int32(sr->readers[0].header, line_ref, "AN", &vAN, &nAN);
		if ((nAC!=1)||(nAN!=1)) vrb.error("VCF for reference panel needs AC/AN INFO fields to be present");
		calt = vAC[0]; cref = (vAN[0] - vAC[0]);

		if (std::min(calt, cref) == 0 && !keep_mono) continue;

		#ifdef __XSI__
			ngt_ref = c_xcf_get_genotypes(c_xcf_p, 0, sr->readers[0].header, line_ref, &gt_arr_ref, &ngt_arr_ref);
		#else
			ngt_ref = bcf_get_genotypes(sr->readers[0].header, line_ref, &gt_arr_ref, &ngt_arr_ref);
		#endif
			//		ngt_ref = bcf_get_genotypes(sr->readers[0].header, line_ref, &gt_arr_ref, &ngt_arr_ref);

		if (i_site==0) set_ploidy_ref(ngt_ref, gt_arr_ref, ngt_arr_ref);
		line_max_ploidy = ngt_ref/n_ref_samples;
		idx_ref_hap=0, cref = 0, calt = 0;

		if (H.flag_common[i_site])
		{
			for(int i = 0 ; i < n_ref_samples ; ++i) {
				ptr = gt_arr_ref + i*line_max_ploidy;
				for (int j=0; j<ploidy_ref_samples[i]; j++) {
					a = (bcf_gt_allele(ptr[j])==1);
					if (a) H.HvarRef.set(i_common, idx_ref_hap, a);
					++idx_ref_hap;
					a?++calt:++cref;
				}
			}
		}
		else
		{
			for(int i = 0 ; i < n_ref_samples ; ++i) {
				ptr = gt_arr_ref + i*line_max_ploidy;
				for (int j=0; j<ploidy_ref_samples[i]; j++) {
					a = (bcf_gt_allele(ptr[j])==1);
					if (a != H.major_alleles[i_site]) H.ShapRef[idx_ref_hap].push_back(i_site);
					++idx_ref_hap;
					a?++calt:++cref;
				}
			}
		}
		//Check that AC and AN used for classification were actually correct
		if ((cref!=V.vec_pos[i_site]->cref)||(calt!=V.vec_pos[i_site]->calt))
			vrb.error("AC/AN INFO fields in VCF are inconsistent with GT field, update the values in the VCF");
		i_common+=H.flag_common[i_site];
		i_site++;
		prog_bar+=prog_step;
		vrb.progress("  * Reference panel parsing ", prog_bar);
	}
	#ifdef __XSI__
		c_xcf_delete(c_xcf_p);
	#endif
	free(gt_arr_ref);

	// Report
	vrb.bullet("Reference panel parsing done ("+stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

