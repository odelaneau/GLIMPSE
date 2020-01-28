#include "call_set_header.h"

void call_set::readData(vector < string > & ftruth, vector < string > & festimated, vector < string > & ffrequencies, vector < string > & region) {
	tac.clock();
	int n_true_samples, n_esti_samples;
	unsigned long int n_variants_all_chromosomes = 0;
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

		// If first file; we initialize data structures
		if (f == 0) {
			// Processing samples + overlap
			N = 0;
			n_true_samples = bcf_hdr_nsamples(sr->readers[0].header);
			n_esti_samples = bcf_hdr_nsamples(sr->readers[1].header);
			mappingT = vector < int > (n_true_samples, -1);
			mappingE = vector < int > (n_esti_samples, -1);
			for (int i = 0 ; i < n_true_samples ; i ++) {
				for (int j = 0 ; j < n_esti_samples ; j ++) {
					string ts = string(sr->readers[0].header->samples[i]);
					string es = string(sr->readers[1].header->samples[j]);
					if (ts == es) {
						mappingT[i] = N;
						mappingE[j] = N;
						samples.push_back(ts);
						N ++;
					}
				}
			}

			// Allocating data structures
			genotype_spl_errors = vector < int > (3 * N, 0);
			genotype_spl_totals = vector < int > (3 * N, 0);
			genotype_bin_errors = vector < int > (3 * L, 0);
			genotype_bin_totals = vector < int > (3 * L, 0);
			genotype_cal_errors = vector < int > (3 * N_BIN_CAL, 0);
			genotype_cal_totals = vector < int > (3 * N_BIN_CAL, 0);
			rsquared_bin = vector < stats2D > (L);
			rsquared_spl = vector < stats2D > (N);
			rsquared_binspl = vector < stats2D > (N*L);
			frequency_bin = vector < stats1D > (L);
			frequency_binspl = vector < stats1D > (N*L);
			vrb.bullet("#overlapping samples = " + stb.str(N));
		}

		unsigned long nvarianttot = 0, nvariantval = 0, nset = 0, ngenoval = 0, n_errors = 0;
		int ngl_t, ngl_arr_t = 0, *gl_arr_t = NULL, int_swap;
		float *af_arr_f = NULL, *ds_arr_e = NULL, *gp_arr_e = NULL, float_swap;
		int nds_e, nds_arr_e = 0;
		int naf_f, naf_arr_f = 0;
		int ngp_e, ngp_arr_e = 0;
		vector < unsigned char > PLs = vector < unsigned char > (3*N, 0);
		vector < float > DSs = vector < float > (N, 0.0f);
		vector < float > GPs = vector < float > (3*N, 0.0f);


		bcf1_t * line_t, * line_e, * line_f;
		while ((nset = bcf_sr_next_line (sr))) {
			if (nset == 3) {
				line_t =  bcf_sr_get_line(sr, 0);
				line_e =  bcf_sr_get_line(sr, 1);
				line_f =  bcf_sr_get_line(sr, 2);
				if (line_t->n_allele == 2 && line_e->n_allele == 2 && line_f->n_allele == 2) {
					bcf_unpack(line_t, BCF_UN_STR);

					naf_f = bcf_get_info_float(sr->readers[2].header,line_f,"AF",&af_arr_f, &naf_arr_f);
					ngl_t = bcf_get_format_int32(sr->readers[0].header, line_t, "PL", &gl_arr_t, &ngl_arr_t);
					nds_e = bcf_get_format_float(sr->readers[1].header, line_e, "DS", &ds_arr_e, &nds_arr_e);
					ngp_e = bcf_get_format_float(sr->readers[1].header, line_e, "GP", &gp_arr_e, &ngp_arr_e);

					assert(naf_f == 1);
					assert(ngl_t == 3*n_true_samples);
					assert(nds_e == n_esti_samples);
					assert(ngp_e == 3*n_esti_samples);

					// Meta data for variant
					bool flip = (af_arr_f[0] > 0.5);
					bool maf = min(af_arr_f[0], 1.0f - af_arr_f[0]);
					int frq_bin = getFrequencyBin(maf);

					// Read Truth
					for(int i = 0 ; i < n_true_samples ; i ++) {
						int idx = mappingT[i];
						if (idx >= 0) {
							PLs[3*idx+0] = gl_arr_t[3*i+0];
							PLs[3*idx+1] = gl_arr_t[3*i+1];
							PLs[3*idx+2] = gl_arr_t[3*i+2];
							if (flip) { int_swap = PLs[3*idx+2]; PLs[3*idx+2] = PLs[3*idx+0]; PLs[3*idx+0] = int_swap; }
						}
					}

					// Read Estimates
					for(int i = 0 ; i < n_esti_samples ; i ++) {
						int idx = mappingE[i];
						if (idx >= 0) {
							DSs[idx] = ds_arr_e[i];
							GPs[3*idx+0] = gp_arr_e[3*i+0];
							GPs[3*idx+1] = gp_arr_e[3*i+1];
							GPs[3*idx+2] = gp_arr_e[3*i+2];
							if (flip) { float_swap = GPs[3*idx+2]; GPs[3*idx+2] = GPs[3*idx+0]; GPs[3*idx+0] = float_swap; DSs[idx] = 2.0f - DSs[idx]; }
						}
					}

					// Process variant
					bool has_validation = false;
					for (int i = 0 ; i < N ; i ++) {
						if (frq_bin >= 0) {	// Do this variant fall within a given frequency bin?
							int true_genotype = getTruth(PLs[3*i+0], PLs[3*i+1], PLs[3*i+2]);
							if (true_genotype >= 0) {
								int esti_genotype = getMostLikely(GPs[3*i+0], GPs[3*i+1], GPs[3*i+2]);
								int cal_bin = getCalibrationBin(GPs[3*i+0], GPs[3*i+1], GPs[3*i+2]);

								// [0] Overall concordance for verbose
								n_errors += (true_genotype != esti_genotype);

								// [1] Update concordance per sample
								switch (true_genotype) {
								case 0: genotype_spl_errors[3*i+0] += (esti_genotype != 0); genotype_spl_totals[3*i+0]++; break;
								case 1: genotype_spl_errors[3*i+1] += (esti_genotype != 1); genotype_spl_totals[3*i+1]++; break;
								case 2: genotype_spl_errors[3*i+2] += (esti_genotype != 2); genotype_spl_totals[3*i+2]++; break;
								}

								// [2] Update concordance per bin
								switch (true_genotype) {
								case 0:	genotype_bin_errors[3*frq_bin+0] += (esti_genotype != 0); genotype_bin_totals[3*i+0]++; break;
								case 1:	genotype_bin_errors[3*frq_bin+1] += (esti_genotype != 1); genotype_bin_totals[3*i+1]++; break;
								case 2:	genotype_bin_errors[3*frq_bin+2] += (esti_genotype != 2); genotype_bin_totals[3*i+2]++; break;
								}

								// [3] Update concordance per calibration bin
								switch (true_genotype) {
								case 0:	genotype_cal_errors[3*cal_bin+0] += (esti_genotype != 0); genotype_cal_totals[3*i+0]++; break;
								case 1:	genotype_cal_errors[3*cal_bin+1] += (esti_genotype != 1); genotype_cal_totals[3*i+1]++; break;
								case 2:	genotype_cal_errors[3*cal_bin+2] += (esti_genotype != 2); genotype_cal_totals[3*i+2]++; break;
								}

								// [4] Update Rsquare per bin
								rsquared_bin[frq_bin].push(DSs[i], (float)true_genotype);
								frequency_bin[frq_bin].push(maf);

								// [5] Update Rsquare per sample
								rsquared_spl[i].push(DSs[i], (float)true_genotype);

								// [6] Update Rsquare per bin x sample
								rsquared_binspl[frq_bin * N + i].push(DSs[i], (float)true_genotype);
								frequency_binspl[frq_bin * N + i].push(maf);

								// Increment counts
								has_validation = true;
								ngenoval ++;
							}
						}
					}

					// increment number of variants
					nvariantval += has_validation;
					nvarianttot ++;
				}
			}
		}
		free(af_arr_f);
		free(gl_arr_t);
		free(ds_arr_e);
		free(gp_arr_e);
		bcf_sr_destroy(sr);
		n_variants_all_chromosomes += nvariantval;
		vrb.bullet("#variants in the overlap = " + stb.str(nvarianttot));
		vrb.bullet("#variants with validation data = " + stb.str(nvariantval));
		vrb.bullet("#genotypes used in validation = " + stb.str(ngenoval));
		vrb.bullet("%error rate in this file = " + stb.str(n_errors * 100.0f / ngenoval));
	}
	vrb.bullet("Total #variants = " + stb.str(n_variants_all_chromosomes));
}
