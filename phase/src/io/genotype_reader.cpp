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

#include <io/genotype_reader.h>

std::map<std::string, int> mapPloidy = {
		{"1",1},
		{"2",2}
};

genotype_reader::genotype_reader(
		haplotype_set & _H,
		genotype_set & _G,
		variant_map & _V,
		glimpse_mpileup & _M,
		const float _sparse_maf,
		const bool _inputGL,
		const bool _impute_refonly,
		const bool _keep_mono,
		const bool _use_gl_indels) :
				H(_H), G(_G), V(_V), M(_M),
				sparse_maf(_sparse_maf),
				inputGL(_inputGL), impute_refonly(_impute_refonly), n_ref_samples(0), keep_mono(_keep_mono), use_gl_indels(_use_gl_indels)
{
	H.sparse_maf = _sparse_maf;
}

genotype_reader::~genotype_reader() {
}
/*
void genotype_reader::readInitializingSamples(string ftext) {
	string buffer;
	input_file fd (ftext);
	while (getline(fd, buffer)) initializing_samples.insert(buffer);
	vrb.bullet("Subset of "+stb.str(initializing_samples.size()) + " samples used for initialization");
	fd.close();
}
*/
void genotype_reader::readSamplesFilePloidy(std::string ftext) {
	std::string buffer;
	input_file fd (ftext);

	bool repeated_samples=false;
	std::vector<std::string> tokens(2);
	while (getline(fd, buffer))
	{
		if (stb.split(buffer, tokens) != 2) vrb.error("Samples file should contain two columns.");
		if (mapPloidy.find(tokens[1]) == mapPloidy.end()) vrb.error("Unrecognized ploidy. Accepted values: [1,2]");
		if (ploidy_samples.find(tokens[0]) == ploidy_samples.end()) ploidy_samples.insert(std::make_pair(tokens[0],mapPloidy[tokens[1]]));
		else repeated_samples=true;
	}
	if (repeated_samples) vrb.warning("Repeated sample in the file. Kept first instance");
	fd.close();
}

void genotype_reader::set_ploidy_ref(bcf_srs_t * sr, int id_ref)
{
	int ngt_ref, *gt_arr_ref = NULL, ngt_arr_ref = 0;
	int max_ploidy_ref = 0;
	int n_ref_haploid = 0;

	ngt_ref = bcf_get_genotypes(sr->readers[id_ref].header, bcf_sr_get_line(sr, id_ref), &gt_arr_ref, &ngt_arr_ref);
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

	free(gt_arr_ref);

	//Allocate some parts of haplotype set
	H.allocate();

	//Verbose
	vrb.bullet("VCF/BCF scanning [Nm=" + stb.str(M.n_tar_samples) + " (Nmh=" + stb.str(M.n_tar_haploid) + " - Nmd=" + stb.str(M.n_tar_diploid) + ")" + " / Nr=" + stb.str(n_ref_samples) + " (Nrh=" + stb.str(n_ref_haploid) + " - Nrd=" + stb.str(n_ref_samples-n_ref_haploid) + ")" + "]");
	vrb.bullet("VCF/BCF scanning [Reg=" + V.input_gregion + " / L=" + stb.str(H.n_tot_sites) + "] [" + stb.str(tac.rel_time()*1.0/1000, 2) + "s]");
	vrb.wait("  * Reference panel parsing");
	tac.clock();
}



void genotype_reader::initReader(bcf_srs_t *sr, std::string& file, int nthreads)
{
	sr->require_index = 1;
	if (bcf_sr_set_regions(sr, V.input_gregion.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + V.input_gregion + "] in [" + file + "]");
	if (bcf_sr_set_targets(sr, V.input_gregion.c_str(), 0, 0) == -1) vrb.error("Impossible to set target region [" + V.input_gregion + "] in [" + file + "]");

	if (nthreads>1) bcf_sr_set_threads(sr, nthreads);

	//Opening files
	if(!(bcf_sr_add_reader (sr, file.c_str())))
	{
		//we do not build an index here, as the target and reference panel could be accessed in parallel
		if (sr->errnum != idx_load_failed) vrb.error("Failed to open file: " + file + "");
		else vrb.error("Failed to load index of the file: " + file + "");
	}
}

void genotype_reader::readGenotypesAndBAMs(std::string fref,int nthreads)
{
	bcf_srs_t * sr_scan =  bcf_sr_init();
	initReader(sr_scan, fref, nthreads);

	set_ploidy_tar();
	scanGenotypesCommon(sr_scan, 0);

	bcf_sr_destroy(sr_scan);

	bcf_srs_t * sr_parse =  bcf_sr_init();
	initReader(sr_parse, fref, nthreads);
	parseRefGenotypes(sr_parse);
	bcf_sr_destroy(sr_parse);
}

void genotype_reader::initReader(bcf_srs_t *sr, std::string& fmain, std::string& fref, int nthreads)
{
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	if (bcf_sr_set_regions(sr, V.input_gregion.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + V.input_gregion + "] in [" + fmain + "]");
	if (bcf_sr_set_targets(sr, V.input_gregion.c_str(), 0, 0) == -1) vrb.error("Impossible to set target region [" + V.input_gregion + "] in [" + fmain + "]");

	if (nthreads>1) bcf_sr_set_threads(sr, nthreads);

	//Opening files
	std::array<std::string,2> fnames = {fmain,fref};
	for (int reader_id=0; reader_id<2; ++reader_id)
	{
		if(!(bcf_sr_add_reader (sr, fnames[reader_id].c_str())))
		{
			//we do not build an index here, as the target and reference panel could be accessed in parallel
			if (sr->errnum != idx_load_failed) vrb.error("Failed to open file: " + fnames[reader_id] + "");
			else vrb.error("Failed to load index of the file: " + fnames[reader_id] + "");
		}
	}
}


void genotype_reader::readGenotypes(std::string fmain, std::string fref, int nthreads)
{
	bcf_srs_t * sr_scan =  bcf_sr_init();

	initReader(sr_scan, fmain, fref, nthreads);
	scanGenotypes(sr_scan);
	bcf_sr_destroy(sr_scan);

	bcf_srs_t * sr_parse =  bcf_sr_init();
	initReader(sr_parse, fmain, fref, nthreads);
	parseGenotypes(sr_parse);
	bcf_sr_destroy(sr_parse);
}

void genotype_reader::scanGenotypesCommon(bcf_srs_t * sr, int ref_sr_n /* Reference reader number*/) {
	if(sr==nullptr) vrb.error("Error reading input files");

	vrb.wait("  * VCF/BCF scanning");
	tac.clock();

	//Scan Ref file
	n_ref_samples = bcf_hdr_nsamples(sr->readers[ref_sr_n].header);
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
	H.n_com_sites_hq = 0;
	bool warning_out =false;
	int prev_pos=-1;

	while ((nset = bcf_sr_next_line (sr)))
	{
		if (!bcf_sr_has_line(sr, ref_sr_n) || bcf_sr_get_line(sr, ref_sr_n)->n_allele != 2 || (nset < (ref_sr_n+1 /* warning relies on reader number */) && !impute_refonly)) continue; //always ref
		line_ref =  bcf_sr_get_line(sr, ref_sr_n);

		//Get AC / AN
		rAC = bcf_get_info_int32(sr->readers[ref_sr_n].header, line_ref, "AC", &vAC, &nAC);
		rAN = bcf_get_info_int32(sr->readers[ref_sr_n].header, line_ref, "AN", &vAN, &nAN);
		if ((nAC!=1)||(nAN!=1)) vrb.error("VCF for reference panel needs AC/AN INFO fields to be present");
		calt = vAC[0]; cref = (vAN[0] - vAC[0]);

		//TODO: Check AN for missing data too
		if (std::min(calt, cref) == 0 && !keep_mono) { if (!warning_out) vrb.warning("Monomorphic site found [AC field] in reference panel at position: " + stb.str(line_ref->pos + 1) + ". ALL monomorphic variants will be skipped. Please check your reference panel file. Use the --keep-monomorphic-ref-sites option to force GLIMPSE to use monomorphic sites in the reference panel. This warning is shown only once."); warning_out=true; continue;}
		if (ref_sr_n > 0 ) { /* warning relies on reader number */
			n_inc_sites += bcf_sr_has_line(sr, ref_sr_n-1);
		}
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
		vrb.bullet("VCF/BCF scanning [Reg=" + V.input_gregion + " / L=" + stb.str(H.n_tot_sites) + " (Li= " + stb.str(n_inc_sites) + " (" + stb.str(((float)n_inc_sites/H.n_tot_sites)*100.0) + "%) - Lm= " + stb.str(n_mis_sites)  + " (" + stb.str(((float)n_mis_sites/H.n_tot_sites)*100.0) + "%))]");
		if (n_mis_sites > 0) vrb.warning("There are variants of the reference panel with no genotype likelihood in the target dataset. Please compute the likelihoods at all variants in the reference panel in order to use all the information of the sequencing reads.");
	}

	vrb.bullet("VCF/BCF scanning [L=" + stb.str(H.n_tot_sites) + " (Lrare= " + stb.str(H.n_rar_sites) + " (" + stb.str(((float)H.n_rar_sites/H.n_tot_sites)*100.0, 1) + "%) - Lcommon= " + stb.str(H.n_com_sites)  + " (" + stb.str(((float)H.n_com_sites/H.n_tot_sites)*100.0, 1) + "%))]");
}

void genotype_reader::scanGenotypes(bcf_srs_t * sr) {
	M.n_tar_samples = bcf_hdr_nsamples(sr->readers[0].header);
	M.tar_sample_names = std::vector<std::string> (M.n_tar_samples);

	for (int i = 0 ; i < M.n_tar_samples ; i ++)
		M.tar_sample_names[i] = std::string(sr->readers[0].header->samples[i]);

	scanGenotypesCommon(sr, 1);
	set_ploidy_tar();

}
void genotype_reader::readTarGenotypes(std::string fmain, int nthreads)
{
	//TODO //FIXME //TO TEST
	bcf_srs_t * sr_scan =  bcf_sr_init();

	initReader(sr_scan, fmain, nthreads);
	M.n_tar_samples = bcf_hdr_nsamples(sr_scan->readers[0].header);
	M.tar_sample_names = std::vector<std::string> (M.n_tar_samples);

	for (int i = 0 ; i < M.n_tar_samples ; i ++)
		M.tar_sample_names[i] = std::string(sr_scan->readers[0].header->samples[i]);

	set_ploidy_tar();

	vrb.wait("  * VCF/BCF scanning");
	tac.clock();

	//Scan file
	bcf1_t * line_main = NULL;
	int npl_main, npl_arr_main = 0, *pl_arr_main = NULL;
	int ngl_main, ngl_arr_main = 0;
	float *gl_arr_main = NULL;
	const int max_ploidyP1=H.max_ploidy+1;
	int *ptr;
	int main_file_max_ploidy=0;
	float *ptr_f;
	std::vector<variant*> vec_pos_tar;
	sr_scan->max_unpack = BCF_UN_INFO;

	while (bcf_sr_next_line (sr_scan))
	{
		if (!bcf_sr_has_line(sr_scan, 0) || bcf_sr_get_line(sr_scan, 0)->n_allele != 2) continue; //always ref
		line_main =  bcf_sr_get_line(sr_scan, 0);
		int8_t line_type = bcf_get_variant_types(line_main);
		std::vector < variant * > vecS = V.getByRef(line_main->pos + 1, line_main->d.allele[0], line_main->d.allele[1]);
		if (vecS.size() > 0) vec_pos_tar.push_back(vecS[0]);
		else vec_pos_tar.push_back(nullptr);
	}

	bcf_sr_destroy(sr_scan);

	vrb.wait("  * VCF/BCF parsing");
	tac.clock();

	bcf_srs_t * sr_parse =  bcf_sr_init();
	initReader(sr_parse, fmain, nthreads);
	int i_site=0;

	sr_parse->max_unpack = BCF_UN_ALL;
	const int n_ref_haps = H.n_ref_haps;

	while (bcf_sr_next_line (sr_parse))
	{
		if (!bcf_sr_has_line(sr_parse, 0) || bcf_sr_get_line(sr_parse, 0)->n_allele != 2) continue; //always ref
		if (vec_pos_tar[i_site] == nullptr)
		{
			++i_site;
			continue;
		}
		else
		{
			line_main = bcf_sr_get_line(sr_parse, 0);

			//request: switch on/off indel calling from gl file
			if (bcf_get_variant_types(line_main)!=VCF_SNP && !use_gl_indels) { ++i_site; continue;}

			if (inputGL)
			{
				ngl_main = bcf_get_format_float(sr_parse->readers[0].header, line_main, "GL", &gl_arr_main, &ngl_arr_main);
				main_file_max_ploidy = ngl_main/H.n_tar_samples;
				std::array<long int, 3> llk;

				if (max_ploidyP1 == main_file_max_ploidy)
				{
					for(int i = 0 ; i < H.n_tar_samples ; ++i)
					{
						ptr_f = gl_arr_main + i*max_ploidyP1;
						const int ploidy = G.vecG[i]->ploidy;
						unsigned char *gl = G.vecG[i]->GL.data() + (ploidy+1)*vec_pos_tar[i_site]->idx;

						if (bcf_float_is_missing(ptr_f[0]) || bcf_float_is_vector_end(ptr_f[0])) continue;
						if (bcf_float_is_missing(ptr_f[1]) || bcf_float_is_vector_end(ptr_f[1])) continue;
						if (ploidy > 1 && (bcf_float_is_missing(ptr_f[2]) || bcf_float_is_vector_end(ptr_f[2]))) continue;

						llk[0] = lroundf(-10 * ptr_f[0]);
						llk[1] = lroundf(-10 * ptr_f[1]);
						if (ploidy > 1) llk[2] = lroundf(-10 * ptr_f[2]);

						if (llk[0] < 0 || llk[1] < 0 || (ploidy >1 && llk[2] < 0))
							vrb.error("Found positive value in FORMAT/GL for individual [" +
									M.tar_sample_names[i] + "] at site [" +
									bcf_hdr_id2name(sr_parse->readers[0].header, line_main->rid) + ":" + std::to_string(line_main->pos + 1) + "_" +
									line_main->d.allele[0] + "_" + line_main->d.allele[1] + "]");

						gl[0] = (unsigned char) std::min(llk[0], 255l);// (ptr[0]<0 || ptr[0]>=256) ? (unsigned char) 255 : (unsigned char) ptr[0];
						gl[1] = (unsigned char) std::min(llk[1], 255l);//(ptr[1]<0 || ptr[1]>=256) ? (unsigned char) 255 : (unsigned char) ptr[1];

						if (ploidy > 1) gl[2] = (unsigned char) std::min(llk[2], 255l);//(ptr[2]<0 || ptr[2]>=256) ? (unsigned char) 255 : (unsigned char) ptr[2];
						if (!(gl[0] == gl[1] && gl[0] == gl[ploidy]))
							G.vecG[i]->flat[vec_pos_tar[i_site]->idx] = false;
					}
				}
			}
			else
			{
				npl_main = bcf_get_format_int32(sr_parse->readers[0].header, line_main, "PL", &pl_arr_main, &npl_arr_main);

				int main_file_max_ploidy = npl_main/H.n_tar_samples;
				if (max_ploidyP1 == main_file_max_ploidy)
				{
					for(int i = 0 ; i < H.n_tar_samples ; ++i)
					{
						ptr = pl_arr_main + i*max_ploidyP1;
						const int ploidy = G.vecG[i]->ploidy;
						unsigned char *gl = G.vecG[i]->GL.data() + (ploidy+1)*vec_pos_tar[i_site]->idx;

						if (ptr[0]==bcf_int32_missing || ptr[0]==bcf_int32_vector_end) continue;
						if (ptr[1]==bcf_int32_missing || ptr[1]==bcf_int32_vector_end) continue;
						if (ploidy > 1 && (ptr[2]==bcf_int32_missing || ptr[2]==bcf_int32_vector_end)) continue;

						if (ptr[0] < 0 || ptr[1] < 0 || (ploidy >1 && ptr[2] < 0))
							vrb.error("Found negative value in FORMAT/PL for individual [" +
									M.tar_sample_names[i] + "] at site [" +
									bcf_hdr_id2name(sr_parse->readers[0].header, line_main->rid) + ":" + stb.str(line_main->pos + 1) + "_" +
									line_main->d.allele[0] + "_" + line_main->d.allele[1] + "]");

						gl[0] = (unsigned char) std::min(ptr[0], 255);// (ptr[0]<0 || ptr[0]>=256) ? (unsigned char) 255 : (unsigned char) ptr[0];
						gl[1] = (unsigned char) std::min(ptr[1], 255);//(ptr[1]<0 || ptr[1]>=256) ? (unsigned char) 255 : (unsigned char) ptr[1];

						if (ploidy > 1) gl[2] = (unsigned char) std::min(ptr[2], 255);//(ptr[2]<0 || ptr[2]>=256) ? (unsigned char) 255 : (unsigned char) ptr[2];
						if (!(gl[0] == gl[1] && gl[0] == gl[ploidy]))
							G.vecG[i]->flat[vec_pos_tar[i_site]->idx] = false;

						//if (vec_pos_tar[i_site]->getMAC() / (2.0*n_ref_haps) < 1e-5)
						//	G.vecG[i]->flat[vec_pos_tar[i_site]->idx] = true;
					}
				}
			}
			++i_site;
		}
	}
	bcf_sr_destroy(sr_parse);
	free(pl_arr_main);
	free(gl_arr_main);
}

void genotype_reader::set_ploidy_tar()
{
	M.tar_ploidy = std::vector<int> (M.n_tar_samples,2);
	M.tar_ind2gt = std::vector<int> (M.n_tar_samples);
	M.tar_ind2pl = std::vector<int> (M.n_tar_samples);

	M.n_tar_diploid=0;
	M.n_tar_haploid=0;
	if (ploidy_samples.size() > 0)
	{
		for (int i = 0 ; i < M.n_tar_samples ; i ++)
		{
			M.tar_ind2gt[i] = 2*M.n_tar_diploid + M.n_tar_haploid;
			M.tar_ind2pl[i] = 3*M.n_tar_diploid + 2*M.n_tar_haploid;
			M.tar_ploidy[i] = (ploidy_samples.find(M.tar_sample_names[i])!=ploidy_samples.end()) ? ploidy_samples[M.tar_sample_names[i]] : 2;
			(M.tar_ploidy[i] > 1) ? ++M.n_tar_diploid : ++M.n_tar_haploid;
		}
		M.n_tar_haps = 2*M.n_tar_diploid + M.n_tar_haploid;
	}
	else
	{	//default 2 ploidy
		M.n_tar_haps = 2*M.n_tar_samples;
		M.n_tar_diploid = M.n_tar_samples;
		for (int i = 0 ; i < M.n_tar_samples ; i ++)
		{
			M.tar_ind2gt[i] = 2*i;
			M.tar_ind2pl[i] = 3*i;
		}
	}
	M.max_ploidy = 1 + (M.n_tar_diploid != 0);
	M.fploidy = M.max_ploidy;
	if (M.n_tar_diploid > 0 && M.n_tar_haploid > 0) M.fploidy=-M.fploidy;

	//This is crutial, as we have few dependencies
	H.n_tar_samples = M.n_tar_samples;
	H.n_tar_haps = M.n_tar_haps;
	H.max_ploidy = M.max_ploidy;
	H.fploidy = M.fploidy;
	H.tar_ploidy = M.tar_ploidy;
	H.tar_ind2hapid = M.tar_ind2gt;

	//Allocate genotype set
	G.vecG = std::vector < genotype * > (H.n_tar_samples, NULL);
	G.n_ind = H.n_tar_samples;
	for (unsigned int i = 0 ; i < M.n_tar_samples ; i ++) {
		G.vecG[i] = new genotype (M.tar_sample_names[i], i, H.n_tot_sites, H.tar_ploidy[i], H.tar_ind2hapid[i]);
		G.vecG[i]->allocate();
	}
}


void genotype_reader::parseGenotypes(bcf_srs_t * sr) {
	if(sr==nullptr) vrb.error("Error reading input files");

	unsigned int i_site = 0, nset = 0, n_ref_unphased = 0, i_common = 0;
	int npl_main, npl_arr_main = 0, *pl_arr_main = NULL;
	int ngl_main, ngl_arr_main = 0;
	float *gl_arr_main = NULL;
	int ngt_ref, *gt_arr_ref = NULL, ngt_arr_ref = 0;
	bcf1_t * line_main, * line_ref;
	float prog_step = 1.0/H.n_tot_sites;
	float prog_bar = 0.0;
	const int max_ploidyP1=H.max_ploidy+1;
	//sr->max_unpack=BCF_UN_FMT;
	sr->max_unpack=BCF_UN_ALL;

	unsigned int cref = 0, calt = 0;

	int line_max_ploidy=0;
	int idx_ref_hap=0;
	bool a = false;
	int *ptr;
	int main_file_max_ploidy=0;
	float *ptr_f;
	int rAC=0, nAC=0, *vAC = NULL;
	int rAN=0, nAN=0, *vAN = NULL;
	while ((nset = bcf_sr_next_line (sr)))
	{
		if (!bcf_sr_has_line(sr,1) || bcf_sr_get_line(sr, 1)->n_allele != 2 || (nset < 2 && !impute_refonly)) continue; //always ref

		line_ref =  bcf_sr_get_line(sr, 1);
		rAC = bcf_get_info_int32(sr->readers[1].header, line_ref, "AC", &vAC, &nAC);
		rAN = bcf_get_info_int32(sr->readers[1].header, line_ref, "AN", &vAN, &nAN);
		if ((nAC!=1)||(nAN!=1)) vrb.error("VCF for reference panel needs AC/AN INFO fields to be present");
		calt = vAC[0]; cref = (vAN[0] - vAC[0]);

		if (std::min(calt, cref) == 0 && !keep_mono) continue;

		if (i_site==0) set_ploidy_ref(sr,1);

		ngt_ref = bcf_get_genotypes(sr->readers[1].header, line_ref, &gt_arr_ref, &ngt_arr_ref);
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

		//Read target data
		if (nset == 2)
		{
			line_main = bcf_sr_get_line(sr, 0);

			if (inputGL)
			{
				ngl_main = bcf_get_format_float(sr->readers[0].header, line_main, "GL", &gl_arr_main, &ngl_arr_main);

				main_file_max_ploidy = ngl_main/H.n_tar_samples;
				std::array<long int, 3> llk;

				if (max_ploidyP1 == main_file_max_ploidy)
				{
					for(int i = 0 ; i < H.n_tar_samples ; ++i)
					{
						ptr_f = gl_arr_main + i*max_ploidyP1;
						const int ploidy = G.vecG[i]->ploidy;
						unsigned char *gl = G.vecG[i]->GL.data() + (ploidy+1)*i_site;

						if (bcf_float_is_missing(ptr_f[0]) || bcf_float_is_vector_end(ptr_f[0])) continue;
						if (bcf_float_is_missing(ptr_f[1]) || bcf_float_is_vector_end(ptr_f[1])) continue;
						if (ploidy > 1 && (bcf_float_is_missing(ptr_f[2]) || bcf_float_is_vector_end(ptr_f[2]))) continue;

						llk[0] = lroundf(-10 * ptr_f[0]);
						llk[1] = lroundf(-10 * ptr_f[1]);
						if (ploidy > 1) llk[2] = lroundf(-10 * ptr_f[2]);

						if (llk[0] < 0 || llk[1] < 0 || (ploidy >1 && llk[2] < 0))
							vrb.error("Found positive value in FORMAT/GL for individual [" +
									M.tar_sample_names[i] + "] at site [" +
									bcf_hdr_id2name(sr->readers[1].header, line_ref->rid) + ":" + std::to_string(line_ref->pos + 1) + "_" +
									line_ref->d.allele[0] + "_" + line_ref->d.allele[1] + "]");

						gl[0] = (unsigned char) std::min(llk[0], 255l);// (ptr[0]<0 || ptr[0]>=256) ? (unsigned char) 255 : (unsigned char) ptr[0];
						gl[1] = (unsigned char) std::min(llk[1], 255l);//(ptr[1]<0 || ptr[1]>=256) ? (unsigned char) 255 : (unsigned char) ptr[1];

						if (ploidy > 1) gl[2] = (unsigned char) std::min(llk[2], 255l);//(ptr[2]<0 || ptr[2]>=256) ? (unsigned char) 255 : (unsigned char) ptr[2];
						if (!(gl[0] == gl[1] && gl[0] == gl[ploidy])) G.vecG[i]->flat[i_site] = false;
						/*
						if (!(gl[0] == gl[1] && gl[0] == gl[ploidy]))
						{
							G.vecG[i]->flat[i_site] = false;
							int best_gt = 0;
							for (int j=1;j<ploidy+1;++j) if (gl[j] < gl[best_gt]) best_gt=j;

							if (!H.flag_common[i_site] && H.major_alleles[i_site]*2 == )
							{

							}
							if (H.major_alleles[i_site]*2 == best_gt) continue;
							const int f = (rng.getFloat() < 0.5);
							ShapTar[hapid+f].push_back(v);
							if (best_gt!=1) ShapTar[hapid+(1-f)].push_back(v);
						}*/

					}
				}
			}
			else
			{
				npl_main = bcf_get_format_int32(sr->readers[0].header, line_main, "PL", &pl_arr_main, &npl_arr_main);

				int main_file_max_ploidy = npl_main/H.n_tar_samples;
				if (max_ploidyP1 == main_file_max_ploidy)
				{
					for(int i = 0 ; i < H.n_tar_samples ; ++i)
					{
						ptr = pl_arr_main + i*max_ploidyP1;
						const int ploidy = G.vecG[i]->ploidy;
						unsigned char *gl = G.vecG[i]->GL.data() + (ploidy+1)*i_site;

						if (ptr[0]==bcf_int32_missing || ptr[0]==bcf_int32_vector_end) continue;
						if (ptr[1]==bcf_int32_missing || ptr[1]==bcf_int32_vector_end) continue;
						if (ploidy > 1 && (ptr[2]==bcf_int32_missing || ptr[2]==bcf_int32_vector_end)) continue;

						if (ptr[0] < 0 || ptr[1] < 0 || (ploidy >1 && ptr[2] < 0))
							vrb.error("Found negative value in FORMAT/PL for individual [" +
									M.tar_sample_names[i] + "] at site [" +
									bcf_hdr_id2name(sr->readers[1].header, line_ref->rid) + ":" + std::to_string(line_ref->pos + 1) + "_" +
									line_ref->d.allele[0] + "_" + line_ref->d.allele[1] + "]");

						gl[0] = (unsigned char) std::min(ptr[0], 255);// (ptr[0]<0 || ptr[0]>=256) ? (unsigned char) 255 : (unsigned char) ptr[0];
						gl[1] = (unsigned char) std::min(ptr[1], 255);//(ptr[1]<0 || ptr[1]>=256) ? (unsigned char) 255 : (unsigned char) ptr[1];

						if (ploidy > 1) gl[2] = (unsigned char) std::min(ptr[2], 255);//(ptr[2]<0 || ptr[2]>=256) ? (unsigned char) 255 : (unsigned char) ptr[2];
						if (!(gl[0] == gl[1] && gl[0] == gl[ploidy])) G.vecG[i]->flat[i_site] = false;
					}
				}
			}
		}

		//H.ref_mac[i_site] = cref;
		//V.push(new variant (bcf_hdr_id2name(hdr_ref, line_ref->rid), line_ref->pos + 1, string(line_ref->d.id), line_ref->d.allele[0], line_ref->d.allele[1], V.size(), cref, calt));

		//Increment
		//H.n_com_sites += H.flag_common[i_site];
		//H.n_rar_sites += !H.flag_common[i_site];
		i_common+=H.flag_common[i_site];
		i_site++;
		prog_bar+=prog_step;
		vrb.progress("  * VCF/BCF parsing", prog_bar);
	}
	free(pl_arr_main);
	free(gl_arr_main);
	free(gt_arr_ref);

	// Report
	vrb.bullet("VCF/BCF parsing done ("+stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void genotype_reader::parseRefGenotypes(bcf_srs_t * sr) {
	if(sr==nullptr) vrb.error("Error reading reference file");

	tac.clock();
	unsigned int i_site = 0, i_common =0, n_ref_unphased = 0;
	int ngt_ref, *gt_arr_ref = NULL, ngt_arr_ref = 0;
	bcf1_t * line_ref;
	float prog_step = 1.0/H.n_tot_sites;
	float prog_bar = 0.0;
	const int max_ploidyP1=H.max_ploidy+1;
	sr->max_unpack=BCF_UN_ALL;

	unsigned int cref = 0, calt = 0;
	int line_max_ploidy=0;
	int idx_ref_hap=0;
	bool a = false;
	int *ptr;
	float *ptr_f;

	while (bcf_sr_next_line (sr))
	{
		if (bcf_sr_get_line(sr, 0)->n_allele != 2) continue; //always ref

		line_ref =  bcf_sr_get_line(sr, 0);
		if (i_site==0) set_ploidy_ref(sr,0);

		ngt_ref = bcf_get_genotypes(sr->readers[0].header, line_ref, &gt_arr_ref, &ngt_arr_ref);
		line_max_ploidy = ngt_ref/n_ref_samples;
		idx_ref_hap=0, cref = 0, calt = 0;

		for(int i = 0 ; i < n_ref_samples ; ++i)
		{
			ptr = gt_arr_ref + i*line_max_ploidy;
			for (int j=0; j<ploidy_ref_samples[i]; j++) {
				a = (bcf_gt_allele(ptr[j])==1);
				if (H.flag_common[i_site]) H.HvarRef.set(i_common, idx_ref_hap, a);
				else if (a != H.major_alleles[i_site]) H.ShapRef[idx_ref_hap].push_back(i_site);
				++idx_ref_hap;
				calt+=a;
				cref+=!a;
			}
		}

		//Check that AC and AN used for classification were actually correct
		if ((cref!=V.vec_pos[i_site]->cref)||(calt!=V.vec_pos[i_site]->calt)) vrb.error("AC/AN INFO fields in VCF are inconsistent with GT field, update the values in the VCF");
		i_common+=H.flag_common[i_site];
		i_site++;
		prog_bar+=prog_step;
		vrb.progress("  * Reference panel parsing ", prog_bar);
	}
	free(gt_arr_ref);

	// Report
	vrb.bullet("Reference panel parsing done ("+stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}




