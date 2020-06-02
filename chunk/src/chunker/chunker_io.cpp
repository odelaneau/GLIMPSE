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

#include <chunker/chunker_header.h>

void chunker::readData(string fmain, string region) {
	tac.clock();
	vrb.title("Reading input files");
	vrb.bullet("Main      : [" + fmain + "]");

	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	if(!(bcf_sr_add_reader (sr, fmain.c_str()))) vrb.error("Problem opening index file for [" + fmain + "]");
	int nset, n_variants = 0;
	bcf1_t * line_main;
	while ((nset = bcf_sr_next_line (sr))) {
		line_main =  bcf_sr_get_line(sr, 0);
		if (line_main->n_allele == 2) {
			bcf_unpack(line_main, BCF_UN_STR);
			chrID = bcf_hdr_id2name(sr->readers[0].header, line_main->rid);
			positions.push_back(line_main->pos + 1);
			n_variants++;
		}
	}
	bcf_sr_destroy(sr);
	vrb.bullet("#variants = " + stb.str(n_variants) + " (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
