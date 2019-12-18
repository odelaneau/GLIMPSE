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
#include <chunker/chunker_header.h>

void chunker::readData(vector < string > & fmain, vector < string > & fref) {
	assert(fmain.size() == fref.size());
	for (int f = 0 ; f < fmain.size() ; f ++) {
		vrb.title("Reading set of input files [" + stb.str(f+1) + "/" + stb.str(fmain.size()) + "]");
		vrb.bullet("Main      : [" + fmain[f] + "]");
		vrb.bullet("Reference : [" + fref[f] + "]");
		tac.clock();
		bcf_srs_t * sr =  bcf_sr_init();
		sr->collapse = COLLAPSE_NONE;
		sr->require_index = 1;
		if(!(bcf_sr_add_reader (sr, fmain[f].c_str()))) vrb.error("Problem opening index file for [" + fmain[f] + "]");
		if(!(bcf_sr_add_reader (sr, fref[f].c_str()))) vrb.error("Problem opening index file for [" + fref [f] + "]");
		int nset, n_variants = 0;
		bcf1_t * line_main;
		while ((nset = bcf_sr_next_line (sr))) {
			if (nset == 2) {
				line_main =  bcf_sr_get_line(sr, 0);
				if (line_main->n_allele == 2) {
					bcf_unpack(line_main, BCF_UN_STR);
					string chr = bcf_hdr_id2name(sr->readers[0].header, line_main->rid);
					if (C.size() == 0 || C.back() != chr) {
						C.push_back(chr);
						V.push_back(vector < int > ());
					}
					V.back().push_back(line_main->pos + 1);
					n_variants++;
				}
			}
		}
		bcf_sr_destroy(sr);
		vrb.bullet("C=" + stb.str(C.size()) + " L=" + stb.str(n_variants) + " (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	}
}
