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

#include <stater/stater_header.h>

void stater::readDataAndComputeStats(vector < string > files, vector < string > regions){
	tac.clock();
	double sum_avgw = 0.0f;
	int n_main_samples = 0;
	basic_stats spl_rstat_hwe;
	vector < basic_stats > spl_rstat;

	for (int f = 0 ; f < files.size() ; f ++) {
		vrb.title("Reading input file");
		vrb.bullet("VCF/BCF file : [" + files[f] + "] / Region: [" + regions[f] + "]");

		bcf_srs_t * sr =  bcf_sr_init();
		sr->collapse = COLLAPSE_NONE;
		sr->require_index = 1;
		if(!(bcf_sr_add_reader (sr, files[f].c_str()))) vrb.error("Problem opening index file for [" + files[f] + "]");

		//Extract and allocate per-sample infos
		int ntmp_main_samples = bcf_hdr_nsamples(sr->readers[0].header);
		if (n_main_samples == 0) {
			n_main_samples = ntmp_main_samples;
			for (int i = 0 ; i < n_main_samples ; i ++) spl_ids.push_back(string(sr->readers[0].header->samples[i]));
			spl_info = vector < double > (spl_ids.size(), 0.0f);
			spl_avgp = vector < double > (spl_ids.size(), 0.0f);
			spl_avgw = vector < double > (spl_ids.size(), 0.0f);
			spl_rstat = vector < basic_stats > (spl_ids.size());
		} else {
			if (ntmp_main_samples != n_main_samples) vrb.error("Different number of samples between input files!");
			else {
				for (int i = 0 ; i < n_main_samples ; i ++) {
					string sid = string(sr->readers[0].header->samples[i]);
					if (sid != spl_ids[i]) vrb.error("Samples are not ordered the same way between input files!");
				}
			}
		}

		int nset, n_variants = 0;
		bcf1_t * line_main;
		int ngp_arr_main = 0;
		float *gp_arr_main = NULL;
		while ((nset = bcf_sr_next_line (sr))) {
			line_main =  bcf_sr_get_line(sr, 0);
			if (line_main->n_allele == 2) {
				bcf_unpack(line_main, BCF_UN_STR);

				var_chr.push_back(bcf_hdr_id2name(sr->readers[0].header, line_main->rid));
				var_pos.push_back(line_main->pos + 1);
				var_ids.push_back(line_main->d.id);
				var_ref.push_back(string(line_main->d.allele[0]));
				var_alt.push_back(string(line_main->d.allele[1]));

				var_freq.push_back(0.0f);
				var_info.push_back(0.0f);
				var_avgp.push_back(0.0f);

				float ds_sum=0.0, ds2_sum=0.0, ds4_sum=0.0, max_p = 0.0f;
				int ngp_main = bcf_get_format_float(sr->readers[0].header, line_main, "GP", &gp_arr_main, &ngp_arr_main);
				if (ngp_main == 3 * n_main_samples) {
					for(int i = 0 ; i < 3 * n_main_samples ; i += 3) {
						float gp0 = gp_arr_main[i+0];
						float gp1 = gp_arr_main[i+1];
						float gp2 = gp_arr_main[i+2];
						float ds = gp1 + 2*gp2;
						float mp = max(max(gp0,gp1),gp2);

						ds_sum += ds;
						ds2_sum += ds * ds;
						ds4_sum += gp1 + 4.0*gp2;
						max_p += mp;
						spl_rstat[i/3].push(ds);
						spl_avgp[i/3] += mp;
					}

					float freq_alt_main = ds_sum / (2 * n_main_samples);
					float infoscore = (freq_alt_main>0.0 && freq_alt_main<1.0) ? (float)(1.0 - (ds4_sum - ds2_sum) / (2 * n_main_samples * freq_alt_main * (1.0 - freq_alt_main))) : 1.0f;
					infoscore = (infoscore<0.0f)?0.0f:infoscore;
					infoscore = roundf(infoscore * 1000.0) / 1000.0;
					float maf = freq_alt_main;
					if (freq_alt_main > 0.5) maf = 1.0f - freq_alt_main;

					for(int i = 0 ; i < 3 * n_main_samples ; i += 3) {
						float gp0 = gp_arr_main[i+0];
						float gp1 = gp_arr_main[i+1];
						float gp2 = gp_arr_main[i+2];
						float mp = max(max(gp0,gp1),gp2);
						spl_avgw[i/3] += mp * 2 * maf;
					}
					sum_avgw += 2 * maf;

					var_freq.back() = freq_alt_main;
					var_info.back() = infoscore;
					var_avgp.back() = max_p * 1.0f/n_main_samples;
					spl_rstat_hwe.push(2.0f * freq_alt_main);
					n_variants++;
				}
			}
		}
		free(gp_arr_main);
		bcf_sr_destroy(sr);
		vrb.bullet("#variants = " + stb.str(n_variants) + " / #samples = " + stb.str(n_main_samples) + " (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	}

	//Finalization
	for (int i = 0 ; i < n_main_samples ; i ++) {
		//cout << i << " " << spl_rstat[i].variance() << " " << spl_rstat_hwe.variance() << endl;
		spl_info[i] = 1.0f - (spl_rstat[i].variance() / spl_rstat_hwe.variance());
		spl_avgp[i] /= var_pos.size();
		spl_avgw[i] /= sum_avgw;
	}
}
