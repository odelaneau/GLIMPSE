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

#include <sampler/sampler_header.h>

#define OFILE_VCFU	0
#define OFILE_VCFC	1
#define OFILE_BCFC	2

#define GET(n,i)	(((n)>>(i))&1U)
#define TOG(n,i)	((n)^=(1UL<<(i)))

std::map<std::string, uint8_t> mapPloidy = {
		{"1",1},
		{"2",2}
};

std::map<int, string> fploidy_to_msg = {
		{-2, "Mixed haploid/diploid samples in the region"},
		{1,"Only haploid samples in the region"},
		{2,"Only diploid samples in the region"}
};

bool phaseHet(int nmain, int curr_hs, int prev_hs, bool prev_a0, bool prev_a1) {
	int phase0 = 0, phase1 = 0;
	for (int s = 0 ; s < nmain ; s ++) {
		bool prev_s0 = GET(prev_hs, 2*s+0);
		bool prev_s1 = GET(prev_hs, 2*s+1);
		bool curr_s0 = GET(curr_hs, 2*s+0);
		bool curr_s1 = GET(curr_hs, 2*s+1);
		if (prev_s0 == prev_a0 && prev_s1 == prev_a1 && curr_s0 != curr_s1) {
			phase0 += !curr_s0;
			phase1 += curr_s0;
		}
		if (prev_s0 == prev_a1 && prev_s1 == prev_a0 && curr_s0 != curr_s1) {
			phase0 += curr_s0;
			phase1 += !curr_s0;
		}
	}
	if (phase0>phase1) return false;
	else if (phase0<phase1) return true;
	else return rng.flipCoin();
}

void sampler::sample() {
	tac.clock();
	fploidy = 2;
	string filename = options["input"].as < string > ();
	vrb.title("Generating haplotype pairs for [" + filename + "]");

	//Create all input file descriptors
	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	int n_threads  = options["thread"].as < int > ();
	if (n_threads > 1) bcf_sr_set_threads(sr, n_threads);
	if(!(bcf_sr_add_reader (sr, filename.c_str())))
	{
		if (sr->errnum != idx_load_failed) vrb.error("Failed to open file: " + filename + "");
		//Error reading the index! Attempt to create the index now
		//vrb.warning("Index not found for [" + filenames[f] + "]. GLIMPSE will attempt to create an index for the file.");
		//Remove the reader, as it's broken
		bcf_sr_remove_reader (sr, 0);
		//create index using htslib (csi, using default bcftools option 14)
		int ret = bcf_index_build3(filename.c_str(), NULL, 14, options["thread"].as < int > ());

		if (ret != 0)
		{
			if (ret == -2)
				vrb.error("index: failed to open " + filename);
			else if (ret == -3)
				vrb.error("index: " + filename + " is in a format that cannot be usefully indexed");
			else
				vrb.error("index: failed to create index for + " + filename);
		}
		if(!(bcf_sr_add_reader (sr, filename.c_str()))) vrb.error("Problem opening/creating index file for [" + filename + "]");
		else vrb.bullet("Index file for [" + filename + "] has been successfully created.\n");
	}

	//Extract sample IDs
	vrb.bullet("Extract sample IDs");
	int nsamples = bcf_hdr_nsamples(sr->readers[0].header);
	sample_names = std::vector<string> (nsamples);
	for (int i = 0 ; i < nsamples ; i ++) sample_names[i] = string(sr->readers[0].header->samples[i]);
	//vrb.bullet("#samples = " + stb.str(nsamples));

	//Extract #pairs per sample
	bcf_hrec_t * header_record = bcf_hdr_get_hrec(sr->readers[0].header, BCF_HL_GEN, "NMAIN", NULL, NULL);
	if (header_record == NULL) vrb.error("Cannot retrieve NMAIN flag in VCF header");
	int nmain = atoi(header_record->value);
	if (nmain <1 || nmain > 15) vrb.error("NMAIN out of bounds : " + stb.str(nmain));
	vrb.bullet("#main_iterations = " + stb.str(nmain));

	header_record = bcf_hdr_get_hrec(sr->readers[0].header, BCF_HL_GEN, "FPLOIDY", NULL, NULL);
	if (header_record == NULL) vrb.warning("Cannot retrieve FPLOIDY flag in VCF header [" + filename + "], used a version of GLIMPSE < 1.1.0? Assuming diploid genotypes only [FPLOIDY=2].");
	else fploidy = atoi(header_record->value);
	if (fploidy_to_msg.find(fploidy) == fploidy_to_msg.end()) vrb.error("FPLOIDY out of bounds : " + stb.str(fploidy));
	vrb.bullet("FPLOIDY = "+ to_string(fploidy) + " [" + fploidy_to_msg[fploidy] + "]");

	int * buffer = NULL, nbuffer = 0, rbuffer = 0, *hs_fields = NULL, *gt_fields = NULL;
	int n_gt_fields = 0;

	bcf1_t * line;

	//Ploidy
	ploidy = std::vector<int> (nsamples);
	max_ploidy = std::abs(fploidy);
	n_haploid = 0;
	n_diploid = 0;
	if (fploidy > 0)
	{
		fploidy == 2 ? n_diploid = nsamples: n_haploid = nsamples;
		for (int i = 0 ; i < nsamples ; i ++) ploidy[i] = max_ploidy;
	}
	else
	{
		if (bcf_sr_next_line (sr)  == 0) vrb.error("No marker found in region");
		line =  bcf_sr_get_line(sr, 0);

		int ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_fields, &n_gt_fields);
		const int line_max_ploidy = ngt/nsamples;
		assert(line_max_ploidy==max_ploidy); //we do not allow missing data
		for(int i = 0 ; i < nsamples; ++i)
		{
			ploidy[i] = 2 - (gt_fields[max_ploidy*i+1] == bcf_int32_vector_end);
			ploidy[i] > 1? ++n_diploid : ++n_haploid;
		}
		if (n_diploid == 0 && n_haploid == 0) vrb.error("No sample found.");
		bcf_sr_seek (sr, NULL, 0);
	}
	string pl1  = n_haploid!=1? "s" : "";
	string pl2  = n_diploid!=1? "s" : "";
	vrb.bullet("#samples = " + stb.str(n_diploid+n_haploid) + " ["  + stb.str(n_haploid) + " haploid" + pl1 + "/ " + stb.str(n_diploid) + " diploid"+ pl2 + "]");

	//Generating task data
	bool first_het = true;
	bool maximize = !options.count("sample");
	vector < int > sampled_conf = vector < int > (nsamples, -1);
	vector < int > prev_conf = vector < int > (nsamples, -1);
	vector < bool > prev_haps = vector < bool > (2*nsamples, false);
	if (!maximize)  {
		for (int i = 0 ; i < nsamples ; i ++)
			sampled_conf[i] = rng.getInt(nmain);
	}

	//Create output file + header by duplication
	string file_format = "w";
	string fname = options["output"].as < string > ();
	unsigned int file_type = OFILE_VCFU;
	if (fname.size() > 6 && fname.substr(fname.size()-6) == "vcf.gz") { file_format = "wz"; file_type = OFILE_VCFC; }
	if (fname.size() > 3 && fname.substr(fname.size()-3) == "bcf") { file_format = "wb"; file_type = OFILE_BCFC; }
	htsFile * fp = hts_open(fname.c_str(),file_format.c_str());
	if (n_threads > 1) hts_set_threads(fp, n_threads);
	bcf_hdr_t * hdr = bcf_hdr_dup(sr->readers[0].header);
	bcf_hdr_write(fp, hdr);

	//Main loop
	int n_variants = 0;

	while (bcf_sr_next_line (sr)) {
		//Retrieve variant informations
		line =  bcf_sr_get_line(sr, 0);

		//Retrieve GT
		int ngt = bcf_get_genotypes(hdr, line, &gt_fields, &n_gt_fields);
		//assert(ngt == 2*nsamples && n_gt_fields == 2*nsamples);

		//Retrieve HS
		int n_hs_fields = 0;
		int nhs = bcf_get_format_int32(hdr, line, "HS", &hs_fields, &n_hs_fields);
		assert(nhs == nsamples && n_hs_fields == nsamples);

		const int line_max_ploidy = ngt/nsamples;
		assert(line_max_ploidy==max_ploidy); //we do not allow missing data
		//Process & update record
		if (maximize) {
			for (int i = 0 ; i  < nsamples ; i ++)
			{
				if (ploidy[i] > 1)
				{
					bool a0 = (bcf_gt_allele(gt_fields[max_ploidy*i+0])==1);
					bool a1 = (bcf_gt_allele(gt_fields[max_ploidy*i+1])==1);
					if (a0 != a1) {
						int curr_hs = hs_fields[i];
						int prev_hs = prev_conf[i];
						bool prev_a0 = prev_haps[max_ploidy*i+0];
						bool prev_a1 = prev_haps[max_ploidy*i+1];
						if (prev_hs < 0) {
							a0 = false; a1 = true;
							prev_haps[max_ploidy*i+0] = a0;
							prev_haps[max_ploidy*i+1] = a1;
							prev_conf[i] = curr_hs;
						} else if (phaseHet(nmain, curr_hs, prev_hs, prev_a0, prev_a1)) {
							a0 = true; a1 = false;
							prev_haps[max_ploidy*i+0] = a0;
							prev_haps[max_ploidy*i+1] = a1;
							prev_conf[i] = curr_hs;
						} else {
							a0 = false; a1 = true;
							prev_haps[max_ploidy*i+0] = a0;
							prev_haps[max_ploidy*i+1] = a1;
							prev_conf[i] = curr_hs;
						}
					}
					gt_fields[max_ploidy*i+0] = bcf_gt_phased(a0);
					gt_fields[max_ploidy*i+1] = bcf_gt_phased(a1);
				}
				else
				{
					bool a0 = (bcf_gt_allele(gt_fields[max_ploidy*i+0])==1);
					gt_fields[max_ploidy*i+0] = bcf_gt_phased(a0);
					if (max_ploidy > 1) gt_fields[max_ploidy*i+1] = bcf_int32_vector_end;
				}

			}
			bcf_update_genotypes(hdr, line, gt_fields, nsamples*max_ploidy);
		} else {
			for (int i = 0 ; i  < nsamples ; i ++)
			{
				bool a0 = GET(hs_fields[i], max_ploidy*sampled_conf[i]+0);
				gt_fields[max_ploidy*i+0] = bcf_gt_phased(a0);

				if (ploidy[i] > 1)
				{
					bool a1 = GET(hs_fields[i], max_ploidy*sampled_conf[i]+1);
					gt_fields[max_ploidy*i+1] = bcf_gt_phased(a1);
				} else if (max_ploidy > 1) gt_fields[max_ploidy*i+1] = bcf_int32_vector_end;
			}
			bcf_update_genotypes(hdr, line, gt_fields, nsamples*max_ploidy);
		}

		//Clear other unwanted format field (all but GT)
		bcf_update_format(hdr, line, "HS", NULL, 0, BCF_HT_INT);
		bcf_update_format(hdr, line, "GP", NULL, 0, BCF_HT_INT);
		bcf_update_format(hdr, line, "DS", NULL, 0, BCF_HT_INT);
		// Write updated record in file
		bcf_write1(fp, hdr, line);
		//
		n_variants++;
	}
	//Close file descriptors & free arrays
	free(gt_fields);
	free(hs_fields);
	bcf_hdr_destroy(hdr);
	bcf_sr_destroy(sr);
	if (hts_close(fp)) vrb.error("Non zero status when closing VCF/BCF file descriptor");
	//Last verbose
	if (n_variants == 0) vrb.error("No variants to be phased in files");
	vrb.title("Writing completed [L=" + stb.str(n_variants) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");

	vrb.title("Creating index");
	//create index using htslib (csi, using default bcftools option 14)
	if (!bcf_index_build3(options["output"].as < string > ().c_str(), NULL, 14, options["thread"].as < int > ())) vrb.print("Index successfully created");
	else vrb.warning("Problem building the index for the output file. Try to build it using tabix/bcftools.");
}

