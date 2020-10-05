/*******************************************************************************
 * Copyright (C) 2020 Olivier Delaneau, University of Lausanne
 * Copyright (C) 2020 Simone Rubinacci, University of Lausanne
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

genotype_reader::genotype_reader(haplotype_set & _H, genotype_set & _G, variant_map & _V, const string _region) : H(_H), G(_G), V(_V), region(_region) {
	n_variants = 0;
	n_main_samples = 0;
	n_ref_samples = 0;
	n_glikelihoods = 0;
	n_missing=0;
}

genotype_reader::~genotype_reader() {
	n_variants = 0;
	n_main_samples = 0;
	n_ref_samples = 0;
}

void genotype_reader::readInitializingSamples(string ftext) {
	string buffer;
	input_file fd (ftext);
	while (getline(fd, buffer)) initializing_samples.insert(buffer);
	vrb.bullet("Subset of "+stb.str(initializing_samples.size()) + " samples used for initialization");
	fd.close();
}

void genotype_reader::readSamplesFilePloidy(string ftext) {
	string buffer;
	input_file fd (ftext);

	bool repeated_samples=false;
	std::vector<string> tokens(2);
	while (getline(fd, buffer))
	{
		if (stb.split(buffer, tokens) != 2) vrb.error("Samples file should contain two columns.");
		if (mapPloidy.find(tokens[1]) == mapPloidy.end()) vrb.error("Unrecognized ploidy. Accepted values: [1,2]");
		if (ploidy_samples.find(tokens[0]) == ploidy_samples.end()) ploidy_samples.insert(std::make_pair(tokens[0],mapPloidy[tokens[1]]));
		else repeated_samples=true;
	}
	if (repeated_samples) vrb.warning("Repeated sample in the file. Kept first instance");

	vrb.bullet("Set of "+stb.str(ploidy_samples.size()) + " samples with associated ploidy read from sample file.");
	fd.close();
}

void genotype_reader::allocateGenotypes() {
	//Genotypes
	G.vecG = vector < genotype * > (n_main_samples);
	for (unsigned int i = 0 ; i < n_main_samples ; i ++) {
		G.vecG[i] = new genotype (main_sample_names[i], i, n_variants, H.ploidy[i], H.ind2hapid[i]);
		G.vecG[i]->allocate();
	}
	G.n_ind = n_main_samples;
	G.n_hap = H.n_main_haps;
	G.n_site = n_variants;
	//Haplotypes
	H.n_site = n_variants;
	H.n_hap = H.n_main_haps + H.n_ref_haps;
	H.H_opt_var.allocate(H.n_site, H.n_hap);
}

void genotype_reader::scanGenotypes(string fmain, string fref, int nthreads, const bool impute_missing) {
	vrb.wait("  * VCF/BCF scanning");
	tac.clock();
	int n_haploid = 0, n_diploid=0;

	float prog_step = 1.0/n_variants;
	float prog_bar = 0.0;

	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	bcf_sr_set_threads(sr, nthreads);
	if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + region + "] in [" + fmain + "]");
	if(!(bcf_sr_add_reader (sr, fmain.c_str()))) vrb.error("Problem opening index file for [" + fmain + "]");
	if(!(bcf_sr_add_reader (sr, fref.c_str()))) vrb.error("Problem opening index file for [" + fref + "]");

	n_variants = 0;
	n_glikelihoods=0;
	n_missing=0;
	n_main_samples = bcf_hdr_nsamples(sr->readers[0].header);
	n_ref_samples = bcf_hdr_nsamples(sr->readers[1].header);

	//introduce a bit of initilization here.
	// Needed for ploidy
	ploidy_ref_samples = std::vector<int> (n_ref_samples,2);
	main_sample_names = std::vector<string> (n_main_samples);
	H.ploidy = std::vector<int> (n_main_samples,2);
	H.ind2hapid = std::vector<int> (n_main_samples);
	//H.hapid2ind = std::vector<int> (2*n_main_samples);
	for (int i = 0 ; i < n_main_samples ; i ++) main_sample_names[i] = string(sr->readers[0].header->samples[i]);

	if (ploidy_samples.size() > 0) {
		for (int i = 0 ; i < n_main_samples ; i ++)
		{
			H.ind2hapid[i] = 2*n_diploid + n_haploid;
			H.ploidy[i] = (ploidy_samples.find(main_sample_names[i])!=ploidy_samples.end()) ? ploidy_samples[main_sample_names[i]] : 2;
			(H.ploidy[i] > 1) ? ++n_diploid : ++n_haploid;
		}
		H.n_main_haps = 2*n_diploid + n_haploid;
	}
	else
	{	//default 2 ploidy
		H.n_main_haps = 2*n_main_samples;
		n_diploid = n_main_samples;
		for (int i = 0 ; i < n_main_samples ; i ++) H.ind2hapid[i] = 2*i;
	}
	H.max_ploidy = 1 + (n_diploid != 0);

	H.fploidy = H.max_ploidy;
	if (n_diploid > 0 && n_haploid > 0) H.fploidy=-H.fploidy;

	vrb.bullet("#samples = " + stb.str(n_diploid+n_haploid) + " ["  + stb.str(n_haploid) + " haploid" + n_haploid!=1? "s" : "" + "/ " + stb.str(n_diploid) + " diploid"+ n_haploid!=1? "s" : "" + "]");

	//Scan file
	int nset;
	bcf1_t * line_ref;
	bool to_set_ploidy_ref = true;

	int max_ploidy_ref = 0;
	int n_ref_haploid = 0;
	while ((nset = bcf_sr_next_line (sr))) {
		if (nset < 2 && !(impute_missing && bcf_sr_has_line(sr,1))) continue;

		line_ref =  bcf_sr_get_line(sr, 1);
		if (line_ref->n_allele == 2)
		{
			if (to_set_ploidy_ref) //read first line to set the ploidy of the reference panel
			{
				int ngt_ref, *gt_arr_ref = NULL, ngt_arr_ref = 0;
				ngt_ref = bcf_get_genotypes(sr->readers[1].header, line_ref, &gt_arr_ref, &ngt_arr_ref);
				max_ploidy_ref = ngt_ref/n_ref_samples;
				if (max_ploidy_ref < 1 || max_ploidy_ref > 2 || ngt_ref%n_ref_samples != 0) vrb.error("Max ploidy of the reference panel cannot be set to neither 1 or 2.");

				if (max_ploidy_ref == 1) n_ref_haploid=n_ref_samples;
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
				to_set_ploidy_ref=false;
			}
			nset == 2 ?	++n_glikelihoods : ++n_missing;
			++n_variants;
		}

	}
	bcf_sr_destroy(sr);

	if (n_variants == 0) vrb.error("No variants to be imputed in files");

	//hap2ind
	H.hapid2ind = std::vector<int> (H.n_main_haps + H.n_ref_haps);
	int idx_ref_hap=0;
	for(int i = 0 ; i < n_ref_samples ; ++i) for (int j =0; j<ploidy_ref_samples[i]; ++j) H.hapid2ind[idx_ref_hap++]=i;
	for(int i = 0 ; i < n_main_samples ; ++i) for (int j =0; j<H.ploidy[i]; ++j) H.hapid2ind[idx_ref_hap++]=n_ref_samples+i;

	if (initializing_samples.size() > 0) {
		H.initializing_haps = vector < bool > (H.n_ref_haps, false);
		int j = 0;
		for (int i = 0 ; i < n_ref_samples ; i ++) {
			if (initializing_samples.count(string(sr->readers[1].header->samples[i]))) {
				H.initializing_haps[j+0] = true;
				if (ploidy_ref_samples[i] > 1) H.initializing_haps[j+1] = true;
			}
			j+=ploidy_ref_samples[i];
		}
	} else H.initializing_haps = vector < bool > (H.n_ref_haps,true);

	vrb.bullet("VCF/BCF scanning [Reg=" + region + " / L=" + stb.str(n_variants) + " (Lgl= " + stb.str(n_glikelihoods) + " (" + stb.str(((float)n_glikelihoods/n_variants)*100.0) + "%) - Lnl= " + stb.str(n_missing)  + " (" + stb.str(((float)n_missing/n_variants)*100.0) + "%))] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	vrb.bullet("VCF/BCF scanning [Nm=" + stb.str(n_main_samples) + " (Nmh=" + stb.str(n_haploid) + " - Nmd=" + stb.str(n_diploid) + ")" + " / Nr=" + stb.str(n_ref_samples) + " (Nrh=" + stb.str(n_ref_haploid) + " - Nrd=" + stb.str(n_ref_samples-n_ref_haploid) + ")" + "]");

	if (n_missing > 0) vrb.warning("There are variants of the reference panel with no genotype likelihood in the target dataset. Please compute the likelihoods at all variants in the reference panel in order to use all the information of the sequencing reads.");
}

void genotype_reader::readGenotypes(string funphased, string freference, int nthreads, const bool impute_missing) {
	tac.clock();
	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	bcf_sr_set_threads(sr, nthreads);
	bcf_sr_set_regions(sr, region.c_str(), 0);
	bcf_sr_add_reader (sr, funphased.c_str());
	bcf_sr_add_reader (sr, freference.c_str());

	unsigned int i_variant = 0, nset = 0, n_ref_unphased = 0;
	int ngl_main, ngl_arr_main = 0, *gl_arr_main = NULL;
	int ngt_ref, *gt_arr_ref = NULL, ngt_arr_ref = 0;
	bcf1_t * line_main, * line_ref;
	float prog_step = 1.0/n_variants;
	float prog_bar = 0.0;

	while ((nset = bcf_sr_next_line (sr))) {
		if (nset < 2 && !(impute_missing && bcf_sr_has_line(sr,1))) continue;

		line_ref =  bcf_sr_get_line(sr, 1);
		if (line_ref->n_allele == 2)
		{
			bcf_unpack(line_ref, BCF_UN_STR);
			string chr = bcf_hdr_id2name(sr->readers[1].header, line_ref->rid);
			int pos = line_ref->pos + 1;
			string id = string(line_ref->d.id);
			string ref = string(line_ref->d.allele[0]);
			string alt = string(line_ref->d.allele[1]);
			unsigned int cref = 0, calt = 0;

			ngt_ref = bcf_get_genotypes(sr->readers[1].header, line_ref, &gt_arr_ref, &ngt_arr_ref);
			int line_max_ploidy = ngt_ref/n_ref_samples;
			int idx_ref_hap=0;
			for(int i = 0 ; i < n_ref_samples ; ++i)
			{
				int32_t *ptr = gt_arr_ref + i*line_max_ploidy;
				for (int j=0; j<ploidy_ref_samples[i]; j++)
				{
					bool a0 = (bcf_gt_allele(ptr[j])==1);
					H.H_opt_var.set(i_variant, idx_ref_hap, a0);
					++idx_ref_hap;
					a0?calt++:cref++;
				}
			}

			if (nset == 2)
			{
				line_main = bcf_sr_get_line(sr, 0);
				ngl_main = bcf_get_format_int32(sr->readers[0].header, line_main, "PL", &gl_arr_main, &ngl_arr_main);

				int main_file_max_ploidy = ngl_main/n_main_samples;
				std::array<float, 3> tmp;
				float sum=0.0;

				if (main_file_max_ploidy == 3 || main_file_max_ploidy==2)
				{
					for(int i = 0 ; i < n_main_samples ; ++i)
					{
						int32_t *ptr = gl_arr_main + i*main_file_max_ploidy;
						const int ploidy = G.vecG[i]->ploidy;
						unsigned char *gl = G.vecG[i]->GL.data() + (ploidy+1)*i_variant;

						if ( ptr[0]==bcf_int32_missing || ptr[1]==bcf_int32_missing || ploidy>1 ? ptr[2]==bcf_int32_missing : false) continue;
						if ( ptr[0]==bcf_int32_vector_end || ptr[1]==bcf_int32_vector_end || ploidy>1 ? ptr[2]==bcf_int32_vector_end : false) continue;

						gl[0] = (ptr[0]<0 || ptr[0]>=256) ? (unsigned char) 255 : ptr[0];
						gl[1] = (ptr[1]<0 || ptr[1]>=256) ? (unsigned char) 255 : ptr[1];
						if (ploidy > 1) gl[2] = (ptr[2]<0 || ptr[2]>=256) ? (unsigned char) 255 : ptr[2];
					}
				}
			}

			V.push(new variant (chr, pos, id, ref, alt, V.size(), cref, calt));
			++i_variant;
			prog_bar+=prog_step;
			vrb.progress("  * VCF/BCF parsing", prog_bar);
		}
	}
	free(gl_arr_main);
	free(gt_arr_ref);
	bcf_sr_destroy(sr);

	// Report
	vrb.bullet("VCF/BCF parsing done ("+stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
