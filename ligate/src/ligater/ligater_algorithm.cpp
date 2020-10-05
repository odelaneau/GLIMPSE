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

#include <ligater/ligater_header.h>
#include <sys/stat.h>

#define OFILE_VCFU	0
#define OFILE_VCFC	1
#define OFILE_BCFC	2

#define GET(n,i)	(((n)>>(i))&1U)
#define TOG(n,i)	((n)^=(1UL<<(i)))

std::map<int, string> fploidy_to_msg = {
		{-2, "Mixed haploid/diploid samples in the region"},
		{1,"Only haploid samples in the region"},
		{2,"Only diploid samples in the region"}
};


void ligater::updateHS(int * values) {
	for (int i = 0 ; i < n_diploid ; i++) {
		for (int m = 0 ; m < nmain ; m ++) {
			if (switching[i * nmain + m] && (GET(values[diploid_idx[i]], 2*m+0) != GET(values[diploid_idx[i]], 2*m+1))) {
				TOG(values[diploid_idx[i]], 2*m+0);
				TOG(values[diploid_idx[i]], 2*m+1);
			}
		}
	}
}

int ligater::update_switching() {
	//as switching and distances are defined on diploids, this
	//can be blind about the ploidy type

	int n_changes = 0;
	for (int d = 0 ; d < distances.size() ; d += 2) {
		if (switching[d/2]) {
			if (distances[d+0] < distances[d+1]) switching[d/2] = true;
			else {
				switching[d/2] = false;
				n_changes ++;
			}
		} else {
			if (distances[d+0] > distances[d+1]) {
				switching[d/2] = true;
				n_changes ++;
			} else switching[d/2] = false;
		}
	}
	fill(distances.begin(), distances.end(), 0);
	return n_changes;
}

void ligater::update_distances_and_write_record(htsFile * fd, bcf_hdr_t * hdr, bcf1_t * r_buffer, bcf1_t * r_body) {
	int n_buffer_hs_fields = 0, n_body_hs_fields = 0;
	int nhs0 = bcf_get_format_int32(hdr, r_body, "HS", &body_hs_fields, &n_body_hs_fields);
	int nhs1 = bcf_get_format_int32(hdr, r_buffer, "HS", &buffer_hs_fields, &n_buffer_hs_fields);
	assert(nhs0 == nsamples && n_body_hs_fields == nsamples);
	assert(nhs1 == nsamples && n_buffer_hs_fields == nsamples);

	for (int i = 0 ; i < n_diploid; i++) {
		for (int m = 0 ; m < nmain ; m ++)
		{
			bool buf0 = GET(buffer_hs_fields[diploid_idx[i]], 2*m+0);
			bool buf1 = GET(buffer_hs_fields[diploid_idx[i]], 2*m+1);
			bool bod0 = GET(body_hs_fields[diploid_idx[i]], 2*m+0);
			bool bod1 = GET(body_hs_fields[diploid_idx[i]], 2*m+1);

			distances[(i * nmain + m) * 2 + 0] += ((buf0 != bod0) + (buf1 != bod1));
			distances[(i * nmain + m) * 2 + 1] += ((buf0 != bod1) + (buf1 != bod0));
		}

	}

	updateHS(body_hs_fields);
	bcf_update_format_int32(hdr, r_body, "HS", body_hs_fields, nsamples);
	bcf_update_info(hdr, r_body, "BUF", NULL, 0, BCF_HT_INT);  // the type does not matter with n=0 / This removes INFO/BUF field
	bcf_write1(fd, hdr, r_body);
}

void ligater::write_record(htsFile *fd, bcf_hdr_t * hdr, bcf1_t *r_body) {
	int n_body_hs_fields = 0;
	int nhs = bcf_get_format_int32(hdr, r_body, "HS", &body_hs_fields, &n_body_hs_fields);
	assert(nhs == nsamples && n_body_hs_fields == nsamples);
	updateHS(body_hs_fields);
	bcf_update_format_int32(hdr, r_body, "HS", body_hs_fields, nsamples);
	bcf_update_info(hdr, r_body, "BUF", NULL, 0, BCF_HT_INT);  // the type does not matter with n=0 / This removes INFO/BUF field
	bcf_write1(fd, hdr, r_body);
}

void ligater::ligate() {
	tac.clock();
	vrb.title("Ligating chunks");

	//Create all input file descriptors
	vrb.bullet("Create all file descriptors");
	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	for (int f = 0 ; f < nfiles ; f ++)
	{
		//Attempt to read the file with the index
		if(!(bcf_sr_add_reader (sr, filenames[f].c_str())))
		{
			//Error reading the index! Attempt to create the index now
			vrb.warning("Index not found for [" + filenames[f] + "]. GLIMPSE will attempt to create an index for the file.");
			//Remove the reader, as it's broken
			bcf_sr_remove_reader (sr, f);
			//create index using htslib (csi, using default bcftools option 14)
			int ret = bcf_index_build3(filenames[f].c_str(), NULL, 14, options["thread"].as < int > ());

			if (ret != 0)
			{
				if (ret == -2)
					vrb.error("index: failed to open " + filenames[f]);
				else if (ret == -3)
					vrb.error("index: " + filenames[f] + " is in a format that cannot be usefully indexed");
				else
					vrb.error("index: failed to create index for + " + filenames[f]);
			}
			if(!(bcf_sr_add_reader (sr, filenames[f].c_str()))) vrb.error("Problem opening index file for [" + filenames[f] + "]");
			else vrb.bullet("Index file for [" + filenames[f] + "] has been successfully created.");
		}
	}

	//Extract sample IDs
	vrb.bullet("Extract sample IDs");
	nsamples = bcf_hdr_nsamples(sr->readers[0].header);
	for (int i = 0 ; i < nsamples ; i ++) sampleIDs.push_back(string(sr->readers[0].header->samples[i]));
	//vrb.bullet("#samples = " + stb.str(nsamples));

	bcf_hrec_t * header_record = bcf_hdr_get_hrec(sr->readers[0].header, BCF_HL_GEN, "FPLOIDY", NULL, NULL);
	if (header_record == NULL) vrb.warning("Cannot retrieve FPLOIDY flag in VCF header [" + filenames[0] + "], used a version of GLIMPSE < 1.1.0? Assuming diploid genotypes only [FPLOIDY=2].");
	else fploidy = atoi(header_record->value);
	if (fploidy_to_msg.find(fploidy) == fploidy_to_msg.end()) vrb.error("FPLOIDY out of bounds : " + stb.str(fploidy));
	vrb.bullet("FPLOIDY = "+ to_string(fploidy) + " [" + fploidy_to_msg[fploidy] + "]");


	//Check sample overlap between all files
	vrb.bullet("Checking samples and FPLOIDY consistency across all input files");
	for (int f = 1 ; f < nfiles ; f ++) {
		int cfploidy = 2;
		header_record = bcf_hdr_get_hrec(sr->readers[0].header, BCF_HL_GEN, "FPLOIDY", NULL, NULL);
		if (header_record == NULL) vrb.warning("Cannot retrieve FPLOIDY flag in VCF header [" + filenames[f] + "], used a version of GLIMPSE < 1.1.0? Assuming diploid genotypes only [FPLOIDY=2].");
		else cfploidy = atoi(header_record->value);
		if (fploidy != cfploidy) vrb.error("Plody format is different across files [" + filenames[0] + "] and [" + filenames[f] + "] fploidy=" + stb.str(fploidy) + " vs " + stb.str(cfploidy));

		int cnsamples = bcf_hdr_nsamples(sr->readers[f].header);
		if (cnsamples != nsamples) vrb.error("Numbers of samples are different between [" + filenames[0] + "] and [" + filenames[f] + "] n=" + stb.str(nsamples) + " vs " + stb.str(cnsamples));
		for (int i = 0 ; i < nsamples ; i ++) {
			string spl = string(sr->readers[f].header->samples[i]);
			if (spl != sampleIDs[i]) vrb.error("Sample overlap problem between [" + filenames[0] + "] and [" + filenames[f] + "] idx=" + stb.str(i) + " id1=" + sampleIDs[i] + " / id2=" + spl);
		}
	}
	vrb.bullet("Samples overlap and FPLOIDY are consistent across files.");

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
		if (bcf_sr_next_line (sr)  == 0) vrb.error("No marker found in files");

		int * gt_fields = NULL;
		int n_gt_fields = 0;

		bcf1_t * line =  bcf_sr_get_line(sr, 0);

		int ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_fields, &n_gt_fields);
		const int line_max_ploidy = ngt/nsamples;
		assert(line_max_ploidy==max_ploidy); //we do not allow missing data
		for(int i = 0 ; i < nsamples; ++i)
		{
			ploidy[i] = 2 - (gt_fields[max_ploidy*i+1] == bcf_int32_vector_end);

			if (ploidy[i] > 1)	++n_diploid;
			else ++n_haploid;
		}
		if (n_diploid == 0 && n_haploid == 0) vrb.error("No sample found.");
		bcf_sr_seek (sr, NULL, 0);

		//let's keep track of the diploids indeces:
		diploid_idx = std::vector<int>(n_diploid);

		for (int i=0, j=0; i<nsamples; ++i)
		{
			if (ploidy[i] == 2) diploid_idx[j++] = i;
		}
	}
	vrb.bullet("#samples = " + stb.str(n_diploid+n_haploid) + " ["  + stb.str(n_haploid) + " haploid" + n_haploid!=1? "s" : "" + "/ " + stb.str(n_diploid) + " diploid"+ n_haploid!=1? "s" : "" + "]");

	// Extract number of main iterations
	header_record = bcf_hdr_get_hrec(sr->readers[0].header, BCF_HL_GEN, "NMAIN", NULL, NULL);
	if (header_record == NULL) vrb.error("Cannot retrieve NMAIN flag in VCF header");
	nmain = atoi(header_record->value);
	if (nmain <1 || nmain > 15) vrb.error("NMAIN out of bounds : " + stb.str(nmain));

	switching = vector < bool > (nmain * n_diploid, false);
	distances = vector < int > (nmain * max_ploidy * n_diploid, 0);
	vrb.bullet("#main_iterations = " + stb.str(nmain));

	//Create output file + header by duplication
	string file_format = "w";
	string fname = options["output"].as < string > ();
	unsigned int file_type = OFILE_VCFU;
	if (fname.size() > 6 && fname.substr(fname.size()-6) == "vcf.gz") { file_format = "wz"; file_type = OFILE_VCFC; }
	if (fname.size() > 3 && fname.substr(fname.size()-3) == "bcf") { file_format = "wb"; file_type = OFILE_BCFC; }
	htsFile * fp = hts_open(fname.c_str(),file_format.c_str());
	bcf_hdr_t * hdr = bcf_hdr_dup(sr->readers[0].header);
	bcf_hdr_write(fp, hdr);

	//Main loop
	int n_variants = 0;
	int *buffer0=NULL, *buffer1=NULL, nbuffer0=0, nbuffer1=0, rbuffer0=0, rbuffer1=0;
	bcf1_t * line0, * line1;
	while (bcf_sr_next_line (sr)) {

		//DETERMINE READERS HAVING THE CURRENT RECORD
		vector < int > active_readers;
		vector < bool > file_has_record = vector < bool > (nfiles, false);
		for (int f = 0 ; f < nfiles ; f ++) {
			file_has_record[f] = bcf_sr_has_line(sr, f);
			if (file_has_record[f]) active_readers.push_back(f);
		}

		//Retrieve variant informations
		line0 =  bcf_sr_get_line(sr, active_readers[0]);
		bcf_unpack(line0, BCF_UN_STR);
		string chr = string(bcf_hdr_id2name(sr->readers[active_readers[0]].header, line0->rid));
		int position = line0->pos + 1;
		string rsid = string(line0->d.id);
		string ref = string(line0->d.allele[0]);
		string alt = string(line0->d.allele[1]);
		if (active_readers.size() >= 3) vrb.error("Too many files overlap at position [" + stb.str(position) + "]");
		if (active_readers.size() == 0) vrb.error("Not enough files overlap at position [" + stb.str(position) + "]");

		//CASE0: WE ARE NOT IN AN OVERLAPPING REGION [STANDARD CASE]
		if (active_readers.size() == 1) {
			// verbose
			if (prev_readers.size() == 0) vrb.title("First chunk of data from position [" + stb.str(position) + "]");
			//Write Output using current switching
			write_record(fp, hdr, line0);
			//Update current stage of active reader
			for (int r = 0 ; r < prev_readers.size() ; r++) current_stages[prev_readers[r]] = STAGE_DONE;
			current_stages[active_readers[0]] = STAGE_BODY;
		}

		//CASE1: WE ARE IN AN OVERLAPPING REGION [LIGATION CASE]
		if (active_readers.size() == 2) {
			line1 =  bcf_sr_get_line(sr, active_readers[1]);
			bcf_unpack(line0, BCF_UN_STR);
			//Retrieve in which stages the readers were for previous variants
			unsigned char prev_stage0 = current_stages[active_readers[0]];
			unsigned char prev_stage1 = current_stages[active_readers[1]];
			//Retrieve who is buffer, who is main
			rbuffer0 = bcf_get_info_int32(sr->readers[active_readers[0]].header,line0,"BUF",&buffer0, &nbuffer0);
			rbuffer1 = bcf_get_info_int32(sr->readers[active_readers[1]].header,line1,"BUF",&buffer1, &nbuffer1);
			if ((buffer0[0] + buffer1[0]) == 0) vrb.error("Overlap between chunk specific variants");
			if ((buffer0[0] + buffer1[0]) == 2) vrb.error("Overlap between buffer-specific variants");
			//Check in which stages the readers are now
			unsigned char curr_stage0, curr_stage1;
			if (buffer0[0]) {
				switch (prev_stage0) {
				case STAGE_NONE: curr_stage0 = STAGE_UBUF; break;
				case STAGE_UBUF: curr_stage0 = STAGE_UBUF; break;
				case STAGE_BODY: curr_stage0 = STAGE_DBUF; break;
				case STAGE_DBUF: curr_stage0 = STAGE_DBUF; break;
				}
			} else curr_stage0 = STAGE_BODY;
			if (buffer1[0]) {
				switch (prev_stage1) {
				case STAGE_NONE: curr_stage1 = STAGE_UBUF; break;
				case STAGE_UBUF: curr_stage1 = STAGE_UBUF; break;
				case STAGE_BODY: curr_stage1 = STAGE_DBUF; break;
				case STAGE_DBUF: curr_stage1 = STAGE_DBUF; break;
				}
			} else curr_stage1 = STAGE_BODY;
			//Check if stage has changed?
			bool changedStatus = (prev_stage0 != curr_stage0) + (prev_stage1 != curr_stage1);
			// f0 is upBUF / f1 is BODY
			if (curr_stage0 == STAGE_UBUF && curr_stage1 == STAGE_BODY) {
				// update hamming distances
				update_distances_and_write_record(fp, hdr, line0, line1);
			}
			// f1 is upBUF / f0 is BODY
			else if (curr_stage1 == STAGE_UBUF && curr_stage0 == STAGE_BODY) {
				update_distances_and_write_record(fp, hdr, line1, line0);
			}
			// f0 is dwBUF / f1 is BODY
			else if (curr_stage0 == STAGE_DBUF && curr_stage1 == STAGE_BODY) {
				// update switching if necessary
				if (changedStatus) {
					vrb.title("Next chunk of data from position [" + stb.str(position) + "]");
					int nc = update_switching();
					vrb.bullet("#switches=" + stb.str(nc));
				}
				// write F1 using current switching
				write_record(fp, hdr, line1);
			}
			// f1 is dwBUF / f0 is BODY
			else if (curr_stage1 == STAGE_DBUF && curr_stage0 == STAGE_BODY) {
				// update switching if necessary
				if (changedStatus) {
					vrb.title("Next chunk of data from position [" + stb.str(position) + "]");
					int nc = update_switching();
					vrb.bullet("#switches=" + stb.str(nc));
				}
				// write F0 using current switching
				write_record(fp, hdr, line0);
			}
			//IMPOSSIBLE CASE
			else vrb.error("Problematic overlapping between chunks): " + stb.str(curr_stage0) + " " + stb.str(curr_stage1));
			//Update the stages in which the readers are now
			current_stages[active_readers[0]] = curr_stage0;
			current_stages[active_readers[1]] = curr_stage1;
		}
		//Store current readers
		prev_readers = active_readers;
		//
		n_variants++;
	}
	//Close file descriptors & free arrays
	free(body_hs_fields);
	free(buffer_hs_fields);
	bcf_hdr_destroy(hdr);
	bcf_sr_destroy(sr);
	if (hts_close(fp)) vrb.error("Non zero status when closing VCF/BCF file descriptor");

	//Last verbose
	if (n_variants == 0) vrb.error("No variants to be phased in files");

	vrb.title("Writing completed [L=" + stb.str(n_variants) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	vrb.title("Creating index");
	//create index using htslib (csi, using default bcftools option 14)
	if (!bcf_index_build3(string(options["output"].as < string > ()).c_str(), NULL, 14, options["thread"].as < int > ())) vrb.print("Index successfully created");
	else vrb.warning("Problem building the index. Try to build using tabix/bcftools.");
}

