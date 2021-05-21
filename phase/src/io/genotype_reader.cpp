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

genotype_reader::genotype_reader(
		haplotype_set & _H,
		genotype_set & _G,
		variant_map & _V,
		const string _region,
		const bool _impute_missing,
		const bool _inputGL,
		const bool _exclude_repeated_samples) :
				H(_H), G(_G), V(_V), region(_region),
				impute_missing(_impute_missing), inputGL(_inputGL),
				exclude_repeated_samples(_exclude_repeated_samples),
				n_haploid(0),n_diploid(0),
				hdr_ref(NULL), imap(NULL),
				n_variants(0), n_main_samples(0),
				n_ref_samples(0), n_included (0), n_missing(0)
{
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

	//H.ref_mac = vector<int> (n_variants);
	//H.targ_mac = vector<int> (n_variants);
}

void genotype_reader::initReader(bcf_srs_t *sr, string& fmain, string& fref, int nthreads)
{
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + region + "] in [" + fmain + "]");
	if (bcf_sr_set_targets(sr, region.c_str(), 0, 0) == -1) vrb.error("Impossible to set target region [" + region + "] in [" + fmain + "]");

	if (nthreads>1) bcf_sr_set_threads(sr, nthreads);

	//Opening files
	std::array<string,2> fnames = {fmain,fref};
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

void genotype_reader::readGenotypes(string fmain, string fref, int nthreads)
{
	bcf_srs_t * sr_scan =  bcf_sr_init();

	initReader(sr_scan, fmain, fref, nthreads);
	scanGenotypes(sr_scan);
	bcf_sr_destroy(sr_scan);

	bcf_srs_t * sr_parse =  bcf_sr_init();
	initReader(sr_parse, fmain, fref, nthreads);
	parseGenotypes(sr_parse);
	bcf_sr_destroy(sr_parse);

	bcf_hdr_destroy(hdr_ref);
	free(imap);
}

void genotype_reader::scanGenotypes(bcf_srs_t * sr) {
	if(sr==nullptr) vrb.error("Error reading input files");

	vrb.wait("  * VCF/BCF scanning");
	tac.clock();

	n_main_samples = bcf_hdr_nsamples(sr->readers[0].header);
	main_sample_names = std::vector<string> (n_main_samples);
	H.ploidy = std::vector<int> (n_main_samples,2);
	H.ind2hapid = std::vector<int> (n_main_samples);

	std::unordered_set<string> main_samples_names_set;
	for (int i = 0 ; i < n_main_samples ; i ++)
	{
		main_sample_names[i] = string(sr->readers[0].header->samples[i]);
		main_samples_names_set.insert(main_sample_names[i]);
	}

	n_ref_samples = bcf_hdr_nsamples(sr->readers[1].header);

	vector<char*> ptr_list;
	ptr_list.reserve(n_ref_samples);
	for (int i = 0; i < n_ref_samples ; i ++)
	{
		if (!exclude_repeated_samples || main_samples_names_set.find(string(sr->readers[1].header->samples[i])) == main_samples_names_set.end())
			ptr_list.push_back(sr->readers[1].header->samples[i]);
	}
	char** ref_samples_names = &ptr_list[0];
	int n_unique_ref_samples = ptr_list.size();

	if (n_unique_ref_samples == 0) vrb.error("Number of samples in the reference panel: " + std::to_string(n_ref_samples) + ", exclusive samples (not in the target panel): " + std::to_string(n_unique_ref_samples));

	imap = (int*)malloc(n_unique_ref_samples * sizeof(int));
	hdr_ref = bcf_hdr_subset(sr->readers[1].header, n_unique_ref_samples, ref_samples_names, imap);

	if ( !hdr_ref ) vrb.error("Error occurred while subsetting samples");
	if ( n_unique_ref_samples != bcf_hdr_nsamples(hdr_ref) )
	{
		for (int i=0; i<n_unique_ref_samples; i++) if ( imap[i]<0 ) vrb.error("No such sample: " + string(ref_samples_names[i]));
		vrb.error("Error occurred while subsetting samples [imap]");
	}
	if (n_ref_samples - n_unique_ref_samples) vrb.warning("Found " + std::to_string(n_ref_samples - n_unique_ref_samples) + " repeated sample IDs between target and reference panel. GLIMPSE excludes these samples from the reference panel.");
	n_ref_samples = n_unique_ref_samples;

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

	string pl1  = n_haploid!=1? "s" : "";
	string pl2  = n_diploid!=1? "s" : "";

	//Scan file
	int nset;
	n_variants = 0;
	n_included=0;
	n_missing=0;
	sr->max_unpack = BCF_UN_STR;

	while ((nset = bcf_sr_next_line (sr)))
	{
		if (!bcf_sr_has_line(sr,1) || bcf_sr_get_line(sr, 1)->n_allele != 2 || (nset < 2 && !impute_missing)) continue; //always ref
		n_included += bcf_sr_has_line(sr,0);
		++n_variants;
	}
	n_missing=n_variants-n_included;
	if (n_variants == 0) vrb.error("No variants to be imputed in files");

	vrb.bullet("#samples = " + stb.str(n_diploid+n_haploid) + " ["  + stb.str(n_haploid) + " haploid" + pl1 + "/ " + stb.str(n_diploid) + " diploid"+ pl2 + "]");
	vrb.bullet("VCF/BCF scanning [Reg=" + region + " / L=" + stb.str(n_variants) + " (Li= " + stb.str(n_included) + " (" + stb.str(((float)n_included/n_variants)*100.0) + "%) - Lm= " + stb.str(n_missing)  + " (" + stb.str(((float)n_missing/n_variants)*100.0) + "%))] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");

	if (n_missing > 0) vrb.warning("There are variants of the reference panel with no genotype likelihood in the target dataset. Please compute the likelihoods at all variants in the reference panel in order to use all the information of the sequencing reads.");
}

void genotype_reader::parseGenotypes(bcf_srs_t * sr) {
	if(sr==nullptr) vrb.error("Error reading input files");

	tac.clock();
	unsigned int i_variant = 0, nset = 0, n_ref_unphased = 0;
	int npl_main, npl_arr_main = 0, *pl_arr_main = NULL;
	int ngl_main, ngl_arr_main = 0;
	float *gl_arr_main = NULL;
	int ngt_ref, *gt_arr_ref = NULL, ngt_arr_ref = 0;
	bcf1_t * line_main, * line_ref;
	float prog_step = 1.0/n_variants;
	float prog_bar = 0.0;
	const int max_ploidyP1=H.max_ploidy+1;
	sr->max_unpack=BCF_UN_FMT;

	unsigned int cref = 0, calt = 0;
	int line_max_ploidy=0;
	int idx_ref_hap=0;
	bool a = false;
	int *ptr;
	int main_file_max_ploidy=0;
	float *ptr_f;

	while ((nset = bcf_sr_next_line (sr))) {
		if (!bcf_sr_has_line(sr,1) || bcf_sr_get_line(sr, 1)->n_allele != 2 || (nset < 2 && !impute_missing)) continue; //always ref

		line_ref =  bcf_sr_get_line(sr, 1);
		bcf_subset(hdr_ref, line_ref, n_ref_samples, imap);
		if (i_variant==0) set_ploidy_ref(line_ref);

		ngt_ref = bcf_get_genotypes(hdr_ref, line_ref, &gt_arr_ref, &ngt_arr_ref);
		line_max_ploidy = ngt_ref/n_ref_samples;
		idx_ref_hap=0, cref = 0, calt = 0;
		for(int i = 0 ; i < n_ref_samples ; ++i)
		{
			ptr = gt_arr_ref + i*line_max_ploidy;
			for (int j=0; j<ploidy_ref_samples[i]; j++)
			{
				a = (bcf_gt_allele(ptr[j])==1);
				H.H_opt_var.set(i_variant, idx_ref_hap, a);
				++idx_ref_hap;
				a?++calt:++cref;
			}
		}

		if (nset == 2)
		{
			line_main = bcf_sr_get_line(sr, 0);

			if (inputGL)
			{
				ngl_main = bcf_get_format_float(sr->readers[0].header, line_main, "GL", &gl_arr_main, &ngl_arr_main);

				main_file_max_ploidy = ngl_main/n_main_samples;
				std::array<long int, 3> llk;

				if (max_ploidyP1 == main_file_max_ploidy)
				{
					for(int i = 0 ; i < n_main_samples ; ++i)
					{
						ptr_f = gl_arr_main + i*max_ploidyP1;
						const int ploidy = G.vecG[i]->ploidy;
						unsigned char *gl = G.vecG[i]->GL.data() + (ploidy+1)*i_variant;

						if (bcf_float_is_missing(ptr_f[0]) || bcf_float_is_vector_end(ptr_f[0])) continue;
						if (bcf_float_is_missing(ptr_f[1]) || bcf_float_is_vector_end(ptr_f[1])) continue;
						if (ploidy > 1 && (bcf_float_is_missing(ptr_f[2]) || bcf_float_is_vector_end(ptr_f[2]))) continue;

						llk[0] = lroundf(-10 * ptr_f[0]);
						llk[1] = lroundf(-10 * ptr_f[1]);
						if (ploidy > 1) llk[2] = lroundf(-10 * ptr_f[2]);

						if (llk[0] < 0 || llk[1] < 0 || (ploidy >1 && llk[2] < 0))
							vrb.error("Found positive value in FORMAT/GL for individual [" +
									main_sample_names[i] + "] at site [" +
									bcf_hdr_id2name(hdr_ref, line_ref->rid) + ":" + to_string(line_ref->pos + 1) + "_" +
									line_ref->d.allele[0] + "_" + line_ref->d.allele[1] + "]");

						gl[0] = (unsigned char) std::min(llk[0], 255l);// (ptr[0]<0 || ptr[0]>=256) ? (unsigned char) 255 : (unsigned char) ptr[0];
						gl[1] = (unsigned char) std::min(llk[1], 255l);//(ptr[1]<0 || ptr[1]>=256) ? (unsigned char) 255 : (unsigned char) ptr[1];

						if (ploidy > 1) gl[2] = (unsigned char) std::min(llk[2], 255l);//(ptr[2]<0 || ptr[2]>=256) ? (unsigned char) 255 : (unsigned char) ptr[2];
					}
				}
				//else assuming missing data?
			}
			else
			{
				npl_main = bcf_get_format_int32(sr->readers[0].header, line_main, "PL", &pl_arr_main, &npl_arr_main);

				int main_file_max_ploidy = npl_main/n_main_samples;
				if (max_ploidyP1 == main_file_max_ploidy)
				{
					for(int i = 0 ; i < n_main_samples ; ++i)
					{
						ptr = pl_arr_main + i*max_ploidyP1;
						const int ploidy = G.vecG[i]->ploidy;
						unsigned char *gl = G.vecG[i]->GL.data() + (ploidy+1)*i_variant;

						if (ptr[0]==bcf_int32_missing || ptr[0]==bcf_int32_vector_end) continue;
						if (ptr[1]==bcf_int32_missing || ptr[1]==bcf_int32_vector_end) continue;
						if (ploidy > 1 && (ptr[2]==bcf_int32_missing || ptr[2]==bcf_int32_vector_end)) continue;

						if (ptr[0] < 0 || ptr[1] < 0 || (ploidy >1 && ptr[2] < 0))
							vrb.error("Found negative value in FORMAT/PL for individual [" +
									main_sample_names[i] + "] at site [" +
									bcf_hdr_id2name(hdr_ref, line_ref->rid) + ":" + to_string(line_ref->pos + 1) + "_" +
									line_ref->d.allele[0] + "_" + line_ref->d.allele[1] + "]");

						gl[0] = (unsigned char) std::min(ptr[0], 255);// (ptr[0]<0 || ptr[0]>=256) ? (unsigned char) 255 : (unsigned char) ptr[0];
						gl[1] = (unsigned char) std::min(ptr[1], 255);//(ptr[1]<0 || ptr[1]>=256) ? (unsigned char) 255 : (unsigned char) ptr[1];

						if (ploidy > 1) gl[2] = (unsigned char) std::min(ptr[2], 255);//(ptr[2]<0 || ptr[2]>=256) ? (unsigned char) 255 : (unsigned char) ptr[2];
					}
				}
				//else assuming missing data?
			}
		}

		//H.ref_mac[i_variant] = cref;
		V.push(new variant (bcf_hdr_id2name(hdr_ref, line_ref->rid), line_ref->pos + 1, string(line_ref->d.id), line_ref->d.allele[0], line_ref->d.allele[1], V.size(), cref, calt));
		++i_variant;
		prog_bar+=prog_step;
		vrb.progress("  * VCF/BCF parsing", prog_bar);
	}
	free(pl_arr_main);
	free(gl_arr_main);
	free(gt_arr_ref);

	// Report
	vrb.bullet("VCF/BCF parsing done ("+stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void genotype_reader::set_ploidy_ref(bcf1_t * line_ref)
{
	int ngt_ref, *gt_arr_ref = NULL, ngt_arr_ref = 0;
	int max_ploidy_ref = 0;
	int n_ref_haploid = 0;

	ngt_ref = bcf_get_genotypes(hdr_ref, line_ref, &gt_arr_ref, &ngt_arr_ref);
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

	H.hapid2ind = std::vector<int> (H.n_main_haps + H.n_ref_haps);
	int idx_ref_hap=0;
	for(int i = 0 ; i < n_ref_samples ; ++i) for (int j =0; j<ploidy_ref_samples[i]; ++j) H.hapid2ind[idx_ref_hap++]=i;
	for(int i = 0 ; i < n_main_samples ; ++i) for (int j =0; j<H.ploidy[i]; ++j) H.hapid2ind[idx_ref_hap++]=n_ref_samples+i;

	if (initializing_samples.size() > 0)
	{
		H.initializing_haps = vector < bool > (H.n_ref_haps, false);
		int j = 0;
		for (int i = 0 ; i < n_ref_samples ; i ++) {
			if (initializing_samples.count(string(hdr_ref->samples[i]))) {
				H.initializing_haps[j+0] = true;
				if (ploidy_ref_samples[i] > 1) H.initializing_haps[j+1] = true;
			}
			j+=ploidy_ref_samples[i];
		}
	}
	else H.initializing_haps = vector < bool > (H.n_ref_haps,true);
	free(gt_arr_ref);

	allocateGenotypes();

	vrb.bullet("VCF/BCF scanning [Nm=" + stb.str(n_main_samples) + " (Nmh=" + stb.str(n_haploid) + " - Nmd=" + stb.str(n_diploid) + ")" + " / Nr=" + stb.str(n_ref_samples) + " (Nrh=" + stb.str(n_ref_haploid) + " - Nrd=" + stb.str(n_ref_samples-n_ref_haploid) + ")" + "]");
	vrb.wait("  * VCF/BCF parsing");
}

