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

#include <chunker/chunker_header.h>

void chunker::readData(std::string fmain, std::string region, int nthreads) {
	tac.clock();
	vrb.title("Reading input files");
	vrb.bullet("Reading file         : [" + fmain + "]");


	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;

	if (nthreads > 1) bcf_sr_set_threads(sr, nthreads);
	if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + region + "] in [" + fmain + "]");
	if(!(bcf_sr_add_reader (sr, fmain.c_str())))
	{
		//we do not build an index here, as the target and reference panel could be accessed in parallel
		if (sr->errnum != idx_load_failed) vrb.error("Failed to open file: " + fmain + "");
		else vrb.error("Failed to load index of the file: " + fmain + "");
	}

	int n_variants = 0;
	int n_comm_variants_cnt=0;
	float* af_ptr = nullptr, af;
	int *itmp=NULL, mitmp=0, tret=0;
	int naf_f, naf_arr_f;
	int rAC=0, nAC=0, *vAC = NULL;
	int rAN=0, nAN=0, *vAN = NULL;

	bcf1_t * line_main;
	while (bcf_sr_next_line (sr)) {
		line_main =  bcf_sr_get_line(sr, 0);
		if (line_main->n_allele == 2)
		{
			n_variants++;
			if (n_comm_variants_cnt < 1) chrID = bcf_hdr_id2name(sr->readers[0].header, line_main->rid);
			rAC = bcf_get_info_int32(sr->readers[0].header, line_main, "AC", &vAC, &nAC);
			rAN = bcf_get_info_int32(sr->readers[0].header, line_main, "AN", &vAN, &nAN);
			if ((nAC!=1)||(nAN!=1)) vrb.error("VCF/BCF needs AC/AN INFO fields to be present");
			//Classify variant
			double MAF = std::min(((vAN[0] - vAC[0])*1.0f)/(vAN[0]), (vAN[0]*1.0f)/(vAN[0]));
			const bool is_common = MAF >= sparse_maf;

			positions_all_mb.push_back(line_main->pos + 1);
			map_positions_all.insert(std::pair < int , int> (line_main->pos + 1, positions_all_mb.size()-1));
			if (is_common)
			{
				n_comm_variants_cnt++;
				positions_common_mb.push_back(line_main->pos + 1);
				common2all.push_back(positions_all_mb.size()-1);
				map_positions_common.insert(std::pair < int , int> (line_main->pos + 1, positions_common_mb.size()-1));
			}

		}
	}
	bcf_sr_destroy(sr);
	if (itmp) free(itmp);
	if (vAC) free(vAC);
	if (vAN) free(vAN);
	if (report_common_variants)
		vrb.bullet("#variants            : [ " + stb.str(n_variants) + " | #common variants = " + stb.str(n_comm_variants_cnt) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	else
		vrb.bullet("#variants            : [" + stb.str(n_variants) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
