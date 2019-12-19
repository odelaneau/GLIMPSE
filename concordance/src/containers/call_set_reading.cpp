#include "call_set_header.h"

void call_set::readData(vector < string > & ftruth, vector < string > & festimated, vector < string > & ffrequencies, vector < string > & region) {
	tac.clock();
	int n_true_samples, n_esti_samples, cpt_overlap;
	vector < int > mappingT, mappingE;
	for (int f = 0 ; f < ftruth.size() ; f ++) {
		vrb.title("Reading set of input files [" + stb.str(f+1) + "/" + stb.str(ftruth.size()) + "]");
		vrb.bullet("Truth       [" + ftruth[f] + "]");
		vrb.bullet("Estimates   [" + festimated[f] + "]");
		vrb.bullet("Frequencies [" + ffrequencies[f] + "]");
		vrb.bullet("Region      [" + region[f] + "]");

		bcf_srs_t * sr =  bcf_sr_init();
		sr->collapse = COLLAPSE_NONE;
		sr->require_index = 1;
		bcf_sr_set_regions(sr, region[f].c_str(), 0);
		bcf_sr_add_reader (sr, ftruth[f].c_str());
		bcf_sr_add_reader (sr, festimated[f].c_str());
		bcf_sr_add_reader (sr, ffrequencies[f].c_str());

		if (f == 0) {
			cpt_overlap = 0;
			n_true_samples = bcf_hdr_nsamples(sr->readers[0].header);
			n_esti_samples = bcf_hdr_nsamples(sr->readers[1].header);
			mappingT = vector < int > (n_true_samples, -1);
			mappingE = vector < int > (n_esti_samples, -1);
			for (int i = 0 ; i < n_true_samples ; i ++) {
				for (int j = 0 ; j < n_esti_samples ; j ++) {
					string ts = string(sr->readers[0].header->samples[i]);
					string es = string(sr->readers[1].header->samples[j]);
					if (ts == es) {
						mappingT[i] = cpt_overlap;
						mappingE[j] = cpt_overlap;
						samples.push_back(ts);
						cpt_overlap ++;
					}
				}
			}
			vrb.bullet("#overlapping samples = " + stb.str(cpt_overlap));
		}

		unsigned int i_variant = 0, nset = 0;
		int ngt_t, ngt_arr_t = 0, *gt_arr_t = NULL;
		int ngl_t, ngl_arr_t = 0, *gl_arr_t = NULL;
		int ndp_t, ndp_arr_t = 0, *dp_arr_t = NULL;
		int ngt_e, ngt_arr_e = 0, *gt_arr_e = NULL;
		int nds_e, nds_arr_e = 0;
		int naf_f, naf_arr_f = 0;
		int ngp_e, ngp_arr_e = 0;
		float *af_arr_f = NULL, *ds_arr_e = NULL, *gp_arr_e = NULL;

		bcf1_t * line_t, * line_e, * line_f;
		while ((nset = bcf_sr_next_line (sr))) {
			if (nset == 3) {
				line_t =  bcf_sr_get_line(sr, 0);
				line_e =  bcf_sr_get_line(sr, 1);
				line_f =  bcf_sr_get_line(sr, 2);
				if (line_t->n_allele == 2 && line_e->n_allele == 2 && line_f->n_allele == 2) {
					bcf_unpack(line_t, BCF_UN_STR);
					chromosomes.push_back(string(bcf_hdr_id2name(sr->readers[0].header, line_t->rid)));
					positions.push_back(line_t->pos + 1);
					rsid.push_back(string(line_t->d.id));
					ref.push_back(string(line_t->d.allele[0]));
					alt.push_back(string(line_t->d.allele[1]));

					naf_f = bcf_get_info_float(sr->readers[2].header,line_f,"AF",&af_arr_f, &naf_arr_f);
					ngt_t = bcf_get_genotypes(sr->readers[0].header, line_t, &gt_arr_t, &ngt_arr_t);
					ngl_t = bcf_get_format_int32(sr->readers[0].header, line_t, "PL", &gl_arr_t, &ngl_arr_t);
					ngt_e = bcf_get_genotypes(sr->readers[1].header, line_e, &gt_arr_e, &ngt_arr_e);
					ndp_t = bcf_get_format_int32(sr->readers[0].header, line_t, "DP", &dp_arr_t, &ndp_arr_t);
					nds_e = bcf_get_format_float(sr->readers[1].header, line_e, "DS", &ds_arr_e, &nds_arr_e);
					ngp_e = bcf_get_format_float(sr->readers[1].header, line_e, "GP", &gp_arr_e, &ngp_arr_e);
					assert(naf_f == 1);
					assert(ngt_t == 2*n_true_samples);
					//assert(ngl_t == 3*n_true_samples);
					assert(ngt_e == 2*n_esti_samples);
					//assert(ndp_t == n_true_samples);
					assert(nds_e == n_esti_samples);
					assert(ngp_e == 3*n_esti_samples);

					if (af_arr_f[0] <= 0.5) {
						frq.push_back(af_arr_f[0]);
						flip.push_back(false);
					} else {
						frq.push_back(1.0-af_arr_f[0]);
						flip.push_back(true);
					}
					GTtrue.push_back(vector < char > (cpt_overlap, -1));
					GLtrue.push_back(vector < unsigned char > (3*cpt_overlap, 0));
					DPtrue.push_back(vector < short > (cpt_overlap, -1));
					int count = 0;
					for(int i = 0 ; i < n_true_samples ; i ++) {
						int idx = mappingT[i];
						if (idx >= 0) {
							bool a0 = bcf_gt_allele(gt_arr_t[2*i+0]);
							bool a1 = bcf_gt_allele(gt_arr_t[2*i+1]);
							count += a0 + a1;
							GTtrue.back()[idx] = a0 + a1;
							if (ndp_t == n_true_samples) DPtrue.back()[idx] = dp_arr_t[i];
							else DPtrue.back()[idx] = 100;
							if (ngl_t > 0) {
								GLtrue.back()[3*idx+0] = gl_arr_t[3*i+0];
								GLtrue.back()[3*idx+1] = gl_arr_t[3*i+1];
								GLtrue.back()[3*idx+2] = gl_arr_t[3*i+2];
							} else {
								GLtrue.back()[3*idx+0] = 255;
								GLtrue.back()[3*idx+1] = 255;
								GLtrue.back()[3*idx+2] = 255;
								GLtrue.back()[3*idx+a0+a1] = 0;
							}
						}
					}

					GTesti.push_back(vector < char > (cpt_overlap, -1));
					DSesti.push_back(vector < float > (cpt_overlap, -1));
					GPesti.push_back(vector < float > ( 3 * cpt_overlap, -1.0));
					for(int i = 0 ; i < n_esti_samples ; i ++) {
						int idx = mappingE[i];
						if (idx >= 0) {
							GTesti.back()[idx] = (bcf_gt_allele(gt_arr_e[2*i+0])==1) + (bcf_gt_allele(gt_arr_e[2*i+1])==1);
							DSesti.back()[idx] = ds_arr_e[i];
							GPesti.back()[3*idx+0] = gp_arr_e[3*i+0];
							GPesti.back()[3*idx+1] = gp_arr_e[3*i+1];
							GPesti.back()[3*idx+2] = gp_arr_e[3*i+2];
						}
					}
					i_variant ++;
				}
			}
		}
		free(af_arr_f);
		free(gt_arr_t);
		free(gt_arr_e);
		free(dp_arr_t);
		free(ds_arr_e);
		free(gp_arr_e);
		bcf_sr_destroy(sr);
		vrb.bullet("L variants overlap = " + stb.str(i_variant));
	}
	vrb.bullet("Total #variants = " + stb.str(positions.size()));
}
