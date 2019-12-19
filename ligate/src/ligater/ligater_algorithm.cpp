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

void ligater::ligate() {
	tac.clock();
	vrb.title("Ligating chunks");

	//Create all file descriptors
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
	switching = vector < bool > (nsamples, false);
	haplotypes = vector < vector < bool > > (16 * nsamples);
	vrb.bullet("#samples = " + stb.str(nsamples));

	//Check sample overlap between all files
	vrb.bullet("Check sample consistency across all input files");
	for (int f = 1 ; f < nfiles ; f ++) {
		int cnsamples = bcf_hdr_nsamples(sr->readers[f].header);
		if (cnsamples != nsamples) vrb.error("Numbers of samples are different between [" + filenames[0] + "] and [" + filenames[f] + "] n=" + stb.str(nsamples) + " vs " + stb.str(cnsamples));

		for (int i = 0 ; i < nsamples ; i ++) {
			string spl = string(sr->readers[f].header->samples[i]);
			if (spl != samplesIDs[i]) vrb.error("Sample overlap problem between [" + filenames[0] + "] and [" + filenames[f] + "] idx=" + stb.str(i) + " id1=" + sampleIDs[i] + " / id2=" + spl);
		}
	}
	vrb.bullet("looks good!");

	//Main loop
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
		string rsid = string(line_t->d.id);
		string ref = string(line_t->d.allele[0]);
		string alt = string(line_t->d.allele[1]);

		//CASE0: WE ARE NOT IN AN OVERLAPPING REGION [STANDARD CASE]
		if (active_readers.size() == 1) {

			body_already_passed[active_readers[0]] = true;
		}

		//CASE1: WE ARE IN AN OVERLAPPING REGION [LIGATION CASE]
		if (active_readers.size() == 2) {
			line1 =  bcf_sr_get_line(sr, active_readers[1]);

			//Retrieve who is buffer, whois main
			rbuffer0 = bcf_get_info_int32(sr->readers[active_readers[0]].header,line_f,"BUFFER",&buffer0, &nbuffer0);
			rbuffer1 = bcf_get_info_int32(sr->readers[active_readers[1]].header,line_f,"BUFFER",&buffer1, &nbuffer1);
			if ((buffer0 + buffer1) == 0) vrb.error("Only main regions are overlapping");
			if ((buffer0 + buffer1) == 2) vrb.error("Only buffer are overlapping");

			//CASE1.0: Downstream buffer0
			if (buffer0 && body_already_passed[active_readers[1]]) {

			}
			//CASE1.1: Downstream buffer1
			if (buffer1 && body_already_passed[active_readers[0]]) {

			}
			//CASE1.2: Upstream buffer0
			if (buffer0 && !body_already_passed[active_readers[1]]) {

			}
			//CASE1.3: Upstream buffer1
			if (buffer1 && !body_already_passed[active_readers[0]]) {

			}











		//CASE: WE ARE IN AN OVERLAPPING REGION [LIGATION CASE]
		if (nset == 2) {
		}


			line_ref =  bcf_sr_get_line(sr, 1);
			if (line_main->n_allele == 2 && line_ref->n_allele == 2) n_variants ++;
		}
	}
	bcf_sr_destroy(sr);
	if (n_variants == 0) vrb.error("No variants to be phased in files");
	vrb.bullet("VCF/BCF scanning [Nm=" + stb.str(n_main_samples) + " / Nr=" + stb.str(n_ref_samples) + " / L=" + stb.str(n_variants) + " / Reg=" + region + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");


}

