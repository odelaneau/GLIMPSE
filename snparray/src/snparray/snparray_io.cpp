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
#include <snparray/snparray_header.h>

#define OFILE_VCFU	0
#define OFILE_VCFC	1
#define OFILE_BCFC	2

const float unphred[256] = { 1.000000e+00, 7.943282e-01, 6.309573e-01, 5.011872e-01, 3.981072e-01, 3.162278e-01, 2.511886e-01, 1.995262e-01, 1.584893e-01, 1.258925e-01, 1.000000e-01, 7.943282e-02, 6.309573e-02, 5.011872e-02, 3.981072e-02, 3.162278e-02, 2.511886e-02, 1.995262e-02, 1.584893e-02, 1.258925e-02, 1.000000e-02, 7.943282e-03, 6.309573e-03, 5.011872e-03, 3.981072e-03, 3.162278e-03, 2.511886e-03, 1.995262e-03, 1.584893e-03, 1.258925e-03, 1.000000e-03, 7.943282e-04, 6.309573e-04, 5.011872e-04, 3.981072e-04, 3.162278e-04, 2.511886e-04, 1.995262e-04, 1.584893e-04, 1.258925e-04, 1.000000e-04, 7.943282e-05, 6.309573e-05, 5.011872e-05, 3.981072e-05, 3.162278e-05, 2.511886e-05, 1.995262e-05, 1.584893e-05, 1.258925e-05, 1.000000e-05, 7.943282e-06, 6.309573e-06, 5.011872e-06, 3.981072e-06, 3.162278e-06, 2.511886e-06, 1.995262e-06, 1.584893e-06, 1.258925e-06, 1.000000e-06, 7.943282e-07, 6.309573e-07, 5.011872e-07, 3.981072e-07, 3.162278e-07, 2.511886e-07, 1.995262e-07, 1.584893e-07, 1.258925e-07, 1.000000e-07, 7.943282e-08, 6.309573e-08, 5.011872e-08, 3.981072e-08, 3.162278e-08, 2.511886e-08, 1.995262e-08, 1.584893e-08, 1.258925e-08, 1.000000e-08, 7.943282e-09, 6.309573e-09, 5.011872e-09, 3.981072e-09, 3.162278e-09, 2.511886e-09, 1.995262e-09, 1.584893e-09, 1.258925e-09, 1.000000e-09, 7.943282e-10, 6.309573e-10, 5.011872e-10, 3.981072e-10, 3.162278e-10, 2.511886e-10, 1.995262e-10, 1.584893e-10, 1.258925e-10, 1.000000e-10, 7.943282e-11, 6.309573e-11, 5.011872e-11, 3.981072e-11, 3.162278e-11, 2.511886e-11, 1.995262e-11, 1.584893e-11, 1.258925e-11, 1.000000e-11, 7.943282e-12, 6.309573e-12, 5.011872e-12, 3.981072e-12, 3.162278e-12, 2.511886e-12, 1.995262e-12, 1.584893e-12, 1.258925e-12, 1.000000e-12, 7.943282e-13, 6.309573e-13, 5.011872e-13, 3.981072e-13, 3.162278e-13, 2.511886e-13, 1.995262e-13, 1.584893e-13, 1.258925e-13, 1.000000e-13, 7.943282e-14, 6.309573e-14, 5.011872e-14, 3.981072e-14, 3.162278e-14, 2.511886e-14, 1.995262e-14, 1.584893e-14, 1.258925e-14, 1.000000e-14, 7.943282e-15, 6.309573e-15, 5.011872e-15, 3.981072e-15, 3.162278e-15, 2.511886e-15, 1.995262e-15, 1.584893e-15, 1.258925e-15, 1.000000e-15, 7.943282e-16, 6.309573e-16, 5.011872e-16, 3.981072e-16, 3.162278e-16, 2.511886e-16, 1.995262e-16, 1.584893e-16, 1.258925e-16, 1.000000e-16, 7.943282e-17, 6.309573e-17, 5.011872e-17, 3.981072e-17, 3.162278e-17, 2.511886e-17, 1.995262e-17, 1.584893e-17, 1.258925e-17, 1.000000e-17, 7.943282e-18, 6.309573e-18, 5.011872e-18, 3.981072e-18, 3.162278e-18, 2.511886e-18, 1.995262e-18, 1.584893e-18, 1.258925e-18, 1.000000e-18, 7.943282e-19, 6.309573e-19, 5.011872e-19, 3.981072e-19, 3.162278e-19, 2.511886e-19, 1.995262e-19, 1.584893e-19, 1.258925e-19, 1.000000e-19, 7.943282e-20, 6.309573e-20, 5.011872e-20, 3.981072e-20, 3.162278e-20, 2.511886e-20, 1.995262e-20, 1.584893e-20, 1.258925e-20, 1.000000e-20, 7.943282e-21, 6.309573e-21, 5.011872e-21, 3.981072e-21, 3.162278e-21, 2.511886e-21, 1.995262e-21, 1.584893e-21, 1.258925e-21, 1.000000e-21, 7.943282e-22, 6.309573e-22, 5.011872e-22, 3.981072e-22, 3.162278e-22, 2.511886e-22, 1.995262e-22, 1.584893e-22, 1.258925e-22, 1.000000e-22, 7.943282e-23, 6.309573e-23, 5.011872e-23, 3.981072e-23, 3.162278e-23, 2.511886e-23, 1.995262e-23, 1.584893e-23, 1.258925e-23, 1.000000e-23, 7.943282e-24, 6.309573e-24, 5.011872e-24, 3.981072e-24, 3.162278e-24, 2.511886e-24, 1.995262e-24, 1.584893e-24, 1.258925e-24, 1.000000e-24, 7.943282e-25, 6.309573e-25, 5.011872e-25, 3.981072e-25, 3.162278e-25, 2.511886e-25, 1.995262e-25, 1.584893e-25, 1.258925e-25, 1.000000e-25, 7.943282e-26, 6.309573e-26, 5.011872e-26, 3.981072e-26, 3.162278e-26};

void snparray::readAndWriteData(string finput, string foutput, string region) {
	tac.clock();
	vrb.title("Reading input & writing output by streaming");

	//Open input file
	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	bcf_sr_set_regions(sr, region.c_str(), 0);
	bcf_sr_add_reader (sr, finput.c_str());

	//Open output file
	string file_format = "w";
	unsigned int file_type = OFILE_VCFU;
	if (foutput.size() > 6 && foutput.substr(foutput.size()-6) == "vcf.gz") { file_format = "wz"; file_type = OFILE_VCFC; }
	if (foutput.size() > 3 && foutput.substr(foutput.size()-3) == "bcf") { file_format = "wb"; file_type = OFILE_BCFC; }
	htsFile * fp = hts_open(foutput.c_str(),file_format.c_str());
	bcf_hdr_t * hdr = bcf_hdr_dup(sr->readers[0].header);
	bcf_hdr_write(fp, hdr);

	//Get number of samples in input
	int n_samples = bcf_hdr_nsamples(sr->readers[0].header);
	vrb.bullet("#samples = " + stb.str(n_samples));

	//Main loop
	unsigned int nset = 0, nfound = 0, nunfound = 0;
	int ngl_t, ngl_arr_t = 0, *gl_arr_t = NULL;
	int ngt_t, ngt_arr_t = 0, *gt_arr_t = NULL;
	int ndp_t, ndp_arr_t = 0, *dp_arr_t = NULL;
	bcf1_t * line;
	while ((nset = bcf_sr_next_line (sr))) {
		line =  bcf_sr_get_line(sr, 0);
		if (line->n_allele == 2) {
			bcf_unpack(line, BCF_UN_STR);

			// Build variant ID
			string chr = bcf_hdr_id2name(sr->readers[0].header, line->rid);
			int pos = line->pos + 1;
			string currID = chr + ":" + stb.str(pos);

			if (snpIDs.count(currID) > 0) {
				// Get data for the variant
				ngl_t = bcf_get_format_int32(sr->readers[0].header, line, "PL", &gl_arr_t, &ngl_arr_t);
				ngt_t = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr_t, &ngt_arr_t);
				ndp_t = bcf_get_genotypes(sr->readers[0].header, line, &dp_arr_t, &ndp_arr_t);

				// Remove PL & DP from variant
				bcf_update_format(hdr, line, "PL", NULL, 0, BCF_HT_INT);
				bcf_update_format(hdr, line, "DP", NULL, 0, BCF_HT_INT);

				// Perform calling & update GT field
				int n_missing = 0;
				for (int i = 0 ; i  < n_samples ; i ++) {
					int genotype = -1;
					if (gl_arr_t[3*i+0] != bcf_int32_missing && gl_arr_t[3*i+1] != bcf_int32_missing && gl_arr_t[3*i+2] != bcf_int32_missing) {
						double p0 = unphred[(unsigned char )gl_arr_t[3*i+0]];
						double p1 = unphred[(unsigned char )gl_arr_t[3*i+1]];
						double p2 = unphred[(unsigned char )gl_arr_t[3*i+2]];
						double sc = p0+p1+p2;
						p0*=sc;p1*=sc;p2*=sc;
						int dp = 0;
						if (dp_arr_t[i] != bcf_int32_missing) dp = dp_arr_t[i];
						if (p0 > minPROB && dp >= minDP) { genotype=0; }
						if (p1 > minPROB && dp >= minDP) { genotype=1; }
						if (p2 > minPROB && dp >= minDP) { genotype=2; }
					}
					switch (genotype) {
					case 0:	gt_arr_t[2*i+0] = bcf_gt_unphased(false); gt_arr_t[2*i+1] = bcf_gt_unphased(false); break;
					case 1:	gt_arr_t[2*i+0] = bcf_gt_unphased(false); gt_arr_t[2*i+1] = bcf_gt_unphased(true); break;
					case 2:	gt_arr_t[2*i+0] = bcf_gt_unphased(true); gt_arr_t[2*i+1] = bcf_gt_unphased(true); break;
					default: gt_arr_t[2*i+0] = bcf_gt_missing; gt_arr_t[2*i+1] = bcf_gt_missing; break;
					}
					n_missing += (genotype>=0);
				}

				// Write updated record if enough data
				double missing_rate = n_missing*1.0f/n_samples;
				if (missing_rate <= minMISS) {
					bcf_update_genotypes(hdr, line, gt_arr_t, n_samples*2);
					bcf_write1(fp, hdr, line);
					nfound ++;
				} else nunfound++;
			} else nunfound ++;
			if ((nfound + nunfound) % 10000 == 0)
				vrb.bullet("#processed = " + stb.str(nfound+nunfound) + " / #found = " + stb.str(nfound) + " / # unfound = " + stb.str(nunfound) + " / rate = " + stb.str(nfound*100.0/(nfound+nunfound)) + "%");
		}
	}

	if (hts_close(fp)) vrb.error("Non zero status when closing VCF/BCF file descriptor");
	free(gt_arr_t);
	free(gl_arr_t);
	bcf_hdr_destroy(hdr);
	bcf_sr_destroy(sr);
	vrb.bullet("#processed = " + stb.str(nfound+nunfound) + " / #found = " + stb.str(nfound) + " / # unfound = " + stb.str(nunfound) + " / rate = " + stb.str(nfound*100.0/(nfound+nunfound)) + "%");
	vrb.bullet("Fraction of SNPs in list found = " + stb.str(nfound*100.0/snpIDs.size()) + "%");
	switch (file_type) {
	case OFILE_VCFU: vrb.bullet("VCF writing [Uncompressed / N=" + stb.str(n_samples) + " / L=" + stb.str(nfound) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_VCFC: vrb.bullet("VCF writing [Compressed / N=" + stb.str(n_samples) + " / L=" + stb.str(nfound) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_BCFC: vrb.bullet("BCF writing [Compressed / N=" + stb.str(n_samples) + " / L=" + stb.str(nfound) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	}
}
