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
#include <io/genotype_reader.h>

genotype_reader::genotype_reader(haplotype_set & _H, genotype_set & _G, variant_map & _V, string _region) : H(_H), G(_G), V(_V) {
	n_variants = 0;
	n_main_samples = 0;
	n_ref_samples = 0;
	region = _region;
}

genotype_reader::~genotype_reader() {
	n_variants = 0;
	n_main_samples = 0;
	n_ref_samples = 0;
	region = "";
}

void genotype_reader::allocateGenotypes() {
	//Genotypes
	G.vecG = vector < genotype * > (n_main_samples);
	for (unsigned int i = 0 ; i < n_main_samples ; i ++) {
		G.vecG[i] = new genotype (i, n_variants);
		G.vecG[i]->allocate();
	}
	G.n_ind = n_main_samples;
	G.n_site = n_variants;
	//Haplotypes
	H.n_hap = 2 * (n_ref_samples + n_main_samples);
	H.n_ref = 2 * n_ref_samples;
	H.n_site = n_variants;
	H.H_opt_var.allocate(H.n_site, H.n_hap);
}

void genotype_reader::scanGenotypes(string fmain, string fref) {
	vrb.wait("  * VCF/BCF scanning");
	tac.clock();
	validation = false;
	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + region + "] in [" + fmain + "]");
	if(!(bcf_sr_add_reader (sr, fmain.c_str()))) vrb.error("Problem opening index file for [" + fmain + "]");
	if(!(bcf_sr_add_reader (sr, fref.c_str()))) vrb.error("Problem opening index file for [" + fref + "]");
	n_variants = 0;
	n_main_samples = bcf_hdr_nsamples(sr->readers[0].header);
	n_ref_samples = bcf_hdr_nsamples(sr->readers[1].header);
	int nset;
	bcf1_t * line_main, * line_ref;
	while ((nset = bcf_sr_next_line (sr))) {
		if (nset == 2) {
			line_main =  bcf_sr_get_line(sr, 0);
			line_ref =  bcf_sr_get_line(sr, 1);
			if (line_main->n_allele == 2 && line_ref->n_allele == 2) n_variants ++;
		}
	}
	bcf_sr_destroy(sr);
	if (n_variants == 0) vrb.error("No variants to be phased in files");
	vrb.bullet("VCF/BCF scanning [Nm=" + stb.str(n_main_samples) + " / Nr=" + stb.str(n_ref_samples) + " / L=" + stb.str(n_variants) + " / Reg=" + region + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void genotype_reader::readGenotypes(string funphased, string freference) {
	tac.clock();
	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	bcf_sr_set_regions(sr, region.c_str(), 0);
	bcf_sr_add_reader (sr, funphased.c_str());
	bcf_sr_add_reader (sr, freference.c_str());
	for (int i = 0 ; i < n_main_samples ; i ++) G.vecG[i]->name = string(sr->readers[0].header->samples[i]);
	unsigned int i_variant = 0, nset = 0, n_ref_unphased = 0;
	int ngl_main, ngl_arr_main = 0, *gl_arr_main = NULL;
	int ngt_ref, *gt_arr_ref = NULL, ngt_arr_ref = 0;
	bcf1_t * line_main, * line_ref;
	while ((nset = bcf_sr_next_line (sr))) {
		if (nset == 2) {
			line_main =  bcf_sr_get_line(sr, 0);
			line_ref =  bcf_sr_get_line(sr, 1);
			if (line_main->n_allele == 2 && line_ref->n_allele == 2) {
				bcf_unpack(line_main, BCF_UN_STR);
				string chr = bcf_hdr_id2name(sr->readers[0].header, line_main->rid);
				int pos = line_main->pos + 1;
				string id = string(line_main->d.id);
				string ref = string(line_main->d.allele[0]);
				string alt = string(line_main->d.allele[1]);
				variant * newV = new variant (chr, pos, id, ref, alt, V.size());
				unsigned int cref = 0, calt = 0, cmis = 0;

				ngt_ref = bcf_get_genotypes(sr->readers[1].header, line_ref, &gt_arr_ref, &ngt_arr_ref); assert(ngt_ref == 2 * n_ref_samples);
				for(int i = 0 ; i < 2 * n_ref_samples ; i += 2) {
					bool a0 = (bcf_gt_allele(gt_arr_ref[i+0])==1);
					bool a1 = (bcf_gt_allele(gt_arr_ref[i+1])==1);
					H.H_opt_var.set(i_variant, i+0, a0);
					H.H_opt_var.set(i_variant, i+1, a1);
					a0?calt++:cref++;
					a1?calt++:cref++;
				}

				ngl_main = bcf_get_format_int32(sr->readers[0].header, line_main, "PL", &gl_arr_main, &ngl_arr_main);
				if (ngl_main == 3 * n_main_samples) {
					for(int i = 0 ; i < 3 * n_main_samples ; i += 3) {
						if (gl_arr_main[i+0] != bcf_float_missing && gl_arr_main[i+1] != bcf_float_missing && gl_arr_main[i+2] != bcf_float_missing) {
							G.vecG[i/3]->GL[3*i_variant+0] = (unsigned char)gl_arr_main[i+0];
							G.vecG[i/3]->GL[3*i_variant+1] = (unsigned char)gl_arr_main[i+1];
							G.vecG[i/3]->GL[3*i_variant+2] = (unsigned char)gl_arr_main[i+2];
							//cout << (unsigned int)gl_arr_main[i+0] << " " << (unsigned int)gl_arr_main[i+1] << " " << (unsigned int)gl_arr_main[i+2]<<  endl;
						}
					}
				}

				newV->cref = cref;newV->calt = calt;newV->cmis = cmis;
				i_variant ++;
				V.push(newV);
				vrb.progress("  * VCF/BCF parsing", i_variant*1.0/n_variants);
			}
		}
	}
	free(gl_arr_main);
	free(gt_arr_ref);
	bcf_sr_destroy(sr);

	// Report
	string str3 = stb.str(tac.rel_time()*1.0/1000, 2) + "s";
	vrb.bullet("VCF/BCF parsing done ("+str3+")");
}


