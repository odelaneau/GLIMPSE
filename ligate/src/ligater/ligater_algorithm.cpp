////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2018 Olivier Delaneau, University of Lausanne
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////
#include <ligater/ligater_header.h>

#define OFILE_VCFU	0
#define OFILE_VCFC	1
#define OFILE_BCFC	2

void ligater::update_distances_and_write_record(htsFile * fd, bcf1_t * r0, bcf1_t * r1) {

}

void ligater::write_record(htsFile *fd, bcf1_t *r) {

}

void ligater::ligate() {
	tac.clock();
	vrb.title("Ligating chunks");

	//Create all input file descriptors
	vrb.bullet("Create all file descriptors");
	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	for (int f = 0 ; f < nfiles ; f ++) {
		if(!(bcf_sr_add_reader (sr, filenames[f].c_str()))) vrb.error("Problem opening index file for [" + filenames[f] + "]");
	}

	//Extract sample IDs
	vrb.bullet("Extract sample IDs");
	nsamples = bcf_hdr_nsamples(sr->readers[0].header);
	for (int i = 0 ; i < nsamples ; i ++) sampleIDs.push_back(string(sr->readers[0].header->samples[i]));
	vrb.bullet("#samples = " + stb.str(nsamples));

	// Extract number of main iterations
	bcf_hrec_t * header_record = bcf_hdr_get_hrec(sr->readers[0].header, BCF_HL_GEN, "NMAIN", NULL, NULL);
	if (header_record == NULL) vrb.error("Cannot retrieve NMAIN flag in VCF header");
	int n_main = atoi(header_record->value);
	if (n_main <1 || n_main > 16) vrb.error("NMAIN out of bounds : " + stb.str(n_main));
	switching = vector < bool > (n_main * nsamples, false);
	distances = vector < int > (n_main * 2 * nsamples, 0);
	vrb.bullet("#main_iterations = " + stb.str(n_main));


	//Check sample overlap between all files
	vrb.bullet("Check sample consistency across all input files");
	for (int f = 1 ; f < nfiles ; f ++) {
		int cnsamples = bcf_hdr_nsamples(sr->readers[f].header);
		if (cnsamples != nsamples) vrb.error("Numbers of samples are different between [" + filenames[0] + "] and [" + filenames[f] + "] n=" + stb.str(nsamples) + " vs " + stb.str(cnsamples));
		for (int i = 0 ; i < nsamples ; i ++) {
			string spl = string(sr->readers[f].header->samples[i]);
			if (spl != sampleIDs[i]) vrb.error("Sample overlap problem between [" + filenames[0] + "] and [" + filenames[f] + "] idx=" + stb.str(i) + " id1=" + sampleIDs[i] + " / id2=" + spl);
		}
	}
	vrb.bullet("looks good!");

	//Create output file + header by duplication
	string file_format = "w";
	string fname = options["output"].as < string > ();
	unsigned int file_type = OFILE_VCFU;
	if (fname.size() > 6 && fname.substr(fname.size()-6) == "vcf.gz") { file_format = "wz"; file_type = OFILE_VCFC; }
	if (fname.size() > 3 && fname.substr(fname.size()-3) == "bcf") { file_format = "wb"; file_type = OFILE_BCFC; }
	htsFile * fp = hts_open(fname.c_str(),file_format.c_str());
	bcf_hdr_t * hdr = bcf_hdr_dup(sr->readers[0].header);

	//Main loop
	int n_variants = 0;
	int buffer0, buffer1, nbuffer0, nbuffer1, rbuffer0, rbuffer1;
	bcf1_t * line0, * line1;
	while (bcf_sr_next_line (sr)) {

		//DETERMINE READERS HAVING THE CURRENT RECORD
		vector < int > active_readers;
		vector < bool > file_has_record = vector < bool > (nfiles, false);
		for (int f = 0 ; f < nfiles ; f ++) {
			file_has_record[f] = bcf_sr_has_line(sr, f);
			active_readers.push_back(f);
		}

		if (active_readers.size() >= 3) vrb.error("Too much overlap at a given position!");
		if (active_readers.size() == 0) vrb.error("Not enough overlap at a given position!");

		//Retrieve variant informations
		line0 =  bcf_sr_get_line(sr, active_readers[0]);
		bcf_unpack(line0, BCF_UN_STR);
		string chr = string(bcf_hdr_id2name(sr->readers[active_readers[0]].header, line0->rid));
		int position = line0->pos + 1;
		string rsid = string(line0->d.id);
		string ref = string(line0->d.allele[0]);
		string alt = string(line0->d.allele[1]);

		//CASE0: WE ARE NOT IN AN OVERLAPPING REGION [STANDARD CASE]
		if (active_readers.size() == 1) {
			//Update switching
			if (prev_readers.size() == 2) {
				assert(switching.size() * 2 == distances.size());
				for (int d = 0 ; d < distances.size() ; d += 2) if (distances[d+0] > distances[d+1]) switching[d/2] = !switching[d/2];
				fill(distances.begin(), distances.end(), 0);
			}

			//Write Output using current switching
			write_record(fp, line0);

			//Update current stage of active reader
			for (int r = 0 ; r < prev_readers.size() ; r++) current_stages[prev_readers[r]] = STAGE_DONE;
			current_stages[active_readers[0]] = STAGE_BODY;
		}

		//CASE1: WE ARE IN AN OVERLAPPING REGION [LIGATION CASE]
		if (active_readers.size() == 2) {
			line1 =  bcf_sr_get_line(sr, active_readers[1]);

			//Retrieve in which stages the readers were for previous variants
			unsigned char prev_stage0 = current_stages[active_readers[0]];
			unsigned char prev_stage1 = current_stages[active_readers[1]];

			//Retrieve who is buffer, who is main
			rbuffer0 = bcf_get_info_int32(sr->readers[active_readers[0]].header,line0,"BUFFER",&buffer0, &nbuffer0);
			rbuffer1 = bcf_get_info_int32(sr->readers[active_readers[1]].header,line1,"BUFFER",&buffer1, &nbuffer1);
			if ((buffer0 + buffer1) == 0) vrb.error("Only main regions are overlapping");
			if ((buffer0 + buffer1) == 2) vrb.error("Only buffer are overlapping");

			//Check in which stages the readers are now
			unsigned char curr_stage0, curr_stage1;
			if (rbuffer0) {
				switch (prev_stage0) {
				case STAGE_NONE: curr_stage0 = STAGE_UBUF; break;
				case STAGE_UBUF: curr_stage0 = STAGE_UBUF; break;
				case STAGE_BODY: curr_stage0 = STAGE_DBUF; break;
				case STAGE_DBUF: curr_stage0 = STAGE_DBUF; break;
				}
			} else curr_stage0 = STAGE_BODY;
			if (rbuffer1) {
				switch (prev_stage1) {
				case STAGE_NONE: curr_stage1 = STAGE_UBUF; break;
				case STAGE_UBUF: curr_stage1 = STAGE_UBUF; break;
				case STAGE_BODY: curr_stage1 = STAGE_DBUF; break;
				case STAGE_DBUF: curr_stage1 = STAGE_DBUF; break;
				}
			} else curr_stage1 = STAGE_BODY;

			//Check is stage has changed?
			bool change0 = (prev_stage0 != curr_stage0);
			bool change1 = (prev_stage1 != curr_stage1);

			// f0 is upBUF / f1 is BODY
			if (curr_stage0 == STAGE_UBUF && curr_stage1 == STAGE_BODY) {
				// update hamming distances
				update_distances_and_write_record(fp, line0, line1);
			}
			// f1 is upBUF / f0 is BODY
			else if (curr_stage1 == STAGE_UBUF && curr_stage0 == STAGE_BODY) {
				update_distances_and_write_record(fp, line1, line0);
			}
			// f0 is dwBUF / f1 is BODY
			else if (curr_stage0 == STAGE_DBUF && curr_stage1 == STAGE_BODY) {
				// write F1 using current switching
				write_record(fp, line1);
			}
			// f1 is dwBUF / f0 is BODY
			else if (curr_stage1 == STAGE_DBUF && curr_stage0 == STAGE_BODY) {
				// write F0 using current switching
				write_record(fp, line0);
			}
			//IMPOSSIBLE CASE
			else vrb.error("Problematic situation: " + stb.str(curr_stage0) + " " + stb.str(curr_stage1));

			//Update the stages in which the readers are now
			current_stages[active_readers[0]] = curr_stage0;
			current_stages[active_readers[1]] = curr_stage1;
		}

		//Store current readers
		prev_readers = active_readers;

		//
		n_variants++;

	}
	bcf_sr_destroy(sr);
	if (n_variants == 0) vrb.error("No variants to be phased in files");
	vrb.bullet("Writing completed [L=" + stb.str(n_variants) + " (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");


}
