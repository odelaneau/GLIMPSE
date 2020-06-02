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

#include "call_set_header.h"

float mean(int n, vector<float>& arr) {
    float sum = 0;
    for (int i = 0; i < n; i++)
        sum += arr[i];

    return sum / n;
}

float stdDev(int n, vector<float>& arr, float mean) {
	double sum = 0;
    for (int i = 0; i < n; i++)
        sum += pow(abs(arr[i] - mean), 2);

    return sqrt(sum / n);
}

float pearson(int n, vector<float>& x, vector<float>& y) {
    float sum = 0;
    float xMean = mean(n, x);
    float xStdDev = stdDev(n, x, xMean);
    float yMean  = mean(n, y);
    float yStdDev = stdDev(n, y, yMean);

    for (int i = 0; i < n; i++)
        sum += (x[i] - xMean) * (y[i] - yMean);

    return sum / (n * xStdDev * yStdDev);
}

float pearson_prec(int n, vector<float>& x, vector<float>& y, float xMean, float xStdDev, float yMean, float yStdDev) {
    float sum = 0;
    float xMean1 = mean(n, x);
    float yMean1 = mean(n, y);
    float yStdDev1 = stdDev(n, y, yMean1);
    float xStdDev1 = stdDev(n, x, xMean1);
    for (int i = 0; i < n; i++)
        sum += (x[i] - xMean1) * (y[i] - yMean1);

    return sum / (n * xStdDev1 * yStdDev1);
}

void call_set::readData(vector < string > & ftruth, vector < string > & festimated, vector < string > & ffrequencies, vector < string > & region, string info_af, int nthreads) {
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
		bcf_sr_set_threads(sr, nthreads);
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
				string ts = string(sr->readers[0].header->samples[i]);
				if (!use_subset_samples || subset_samples.count(ts)>0) {
					for (int j = 0 ; j < n_esti_samples ; j ++) {
						string es = string(sr->readers[1].header->samples[j]);
						if (ts == es) {
							mappingT[i] = N;
							mappingE[j] = N;
							samples.push_back(ts);
							N ++;
						}
					}
				}
			}

			// Allocating data structures
			genotype_spl_errors = vector < unsigned long int > (3 * N, 0);
			genotype_spl_totals = vector < unsigned long int > (3 * N, 0);
			genotype_cal_errors = vector < unsigned long int > (3 * N_BIN_CAL, 0);
			genotype_cal_totals = vector < unsigned long int > (3 * N_BIN_CAL, 0);
			rsquared_spl = vector < stats2D > (N);

			if (L > 0) {
				genotype_bin_errors = vector < unsigned long int > (3 * L, 0);
				genotype_bin_totals = vector < unsigned long int > (3 * L, 0);
				rsquared_bin = vector < stats2D > (L);
				frequency_bin = vector < stats1D > (L);
				avg_rsquared_bin = vector < stats1D > (L);
			} else {
				genotype_bin_errors = vector < unsigned long int > (3 * rsquared_str.size(), 0);
				genotype_bin_totals = vector < unsigned long int > (3 * rsquared_str.size(), 0);
				rsquared_bin = vector < stats2D > (rsquared_str.size());
				frequency_bin = vector < stats1D > (rsquared_str.size());
				avg_rsquared_bin = vector < stats1D > (rsquared_str.size());
			}
			vrb.bullet("#overlapping samples = " + stb.str(N));
		}

		unsigned long nvarianttot = 0, nvariantval = 0, nset = 0, ngenoval = 0, n_errors = 0;
		unsigned long ngenoval1 = 0;

        // counts of valid fields
		int ngl_t;
		int ndp_t;
		int ngt_t;
        int nds_e;
        int ngp_e;
        float naf_f;


		// pointers for the arrays holding retrieved field values
		int *gl_arr_t = nullptr;
		int *dp_arr_t = nullptr;
        char* *gt_arr_t = nullptr;
		float *af_ptr = nullptr;
		float *ds_arr_e = nullptr;
		float *gp_arr_e = nullptr;
		float float_swap;

		// to hold values for last arg to bcf_get_*
        int ngl_arr_t;
        int ndp_arr_t;
        int ngt_arr_t;
        int nds_arr_e;
		int naf_arr_f;
		int ngp_arr_e;


		// container vectors for Sample-level data like genotypes and such
		vector < double > PLs = vector < double > (3*N, 0.0f);
		vector < float > DSs = vector < float > (N, 0.0f);
		vector < float > GPs = vector < float > (3*N, 0.0f);
		vector < int > DPs = vector < int > (N, 0);
		vector <char> GTs = vector<char>(N, 0);

		vector < float > sample_dosages(N);
		vector < float > sample_truth(N);
		float sum_ds_e = 0.0f;
		float sum_ds_t = 0.0f;

		map < string, pair < int, bool > > :: iterator itG;
		bcf1_t * line_t, * line_e, * line_f;
		// loop through each variant
		while ((nset = bcf_sr_next_line (sr))) {
			if (nset == 3) {
			    // line of each the three vcfs
				line_t =  bcf_sr_get_line(sr, 0);
				line_e =  bcf_sr_get_line(sr, 1);
				line_f =  bcf_sr_get_line(sr, 2);

				// checks for biallelicity at all vcfs
				if (line_t->n_allele == 2 && line_e->n_allele == 2 && line_f->n_allele == 2) {
					bcf_unpack(line_t, BCF_UN_STR);

                    // build arrays for this variant
				    naf_f = bcf_get_info_float(sr->readers[2].header, line_f, info_af.c_str(), &af_ptr, &naf_arr_f);
					ngl_t = bcf_get_format_int32(sr->readers[0].header, line_t, "PL", &gl_arr_t, &ngl_arr_t);
					ndp_t = bcf_get_format_int32(sr->readers[0].header, line_t, "DP", &dp_arr_t, &ndp_arr_t);
                    // GT
                    ngt_t = bcf_get_format_string(sr->readers[0].header, line_t, "GT", &gt_arr_t, &ngt_arr_t);

                    nds_e = bcf_get_format_float(sr->readers[1].header, line_e, "DS", &ds_arr_e, &nds_arr_e);
					ngp_e = bcf_get_format_float(sr->readers[1].header, line_e, "GP", &gp_arr_e, &ngp_arr_e);

					sum_ds_e = 0.0f;
					sum_ds_t = 0.0f;
					bool has_validation = false;
					// check that the number of each of these fields is correct.
					// seemingly, naf_f is the number of AF fields with nonzero values. so it will skip 0 af sites
					//if ((naf_f==1)&&(ngl_t==3*n_true_samples)&&(nds_e==n_esti_samples)&&(ngp_e==3*n_esti_samples)&&(ndp_t==n_true_samples)) {
                    if ((ngl_t==3*n_true_samples)&&(nds_e==n_esti_samples)&&(ngp_e==3*n_esti_samples)&&(ndp_t==n_true_samples)) {
						// Meta data for variant
						float af =  af_ptr[0];
						bool flip = (af > 0.5);
						float maf = min(af, 1.0f - af);
						int grp_bin = -1;
						if (L>0) grp_bin = getFrequencyBin(maf);
						else {
							string chr = bcf_hdr_id2name(sr->readers[0].header, line_t->rid);
							int pos = line_t->pos + 1;
							string uuid = chr + "_" + stb.str(pos);
							itG = site2grp.find(uuid);
							if (itG != site2grp.end()) {
								grp_bin = itG->second.first;
								itG->second.second = true;
							}
						}

						// Read Truth
						for(int i = 0 ; i < n_true_samples ; i ++) {
							int idx = mappingT[i];
							if (idx >= 0) {
                                // if truth genotypes has
								if (gl_arr_t[3*i+0] != bcf_int32_missing && gl_arr_t[3*i+1] != bcf_int32_missing && gl_arr_t[3*i+2] != bcf_int32_missing) {
								    // compute PLs off GLs
								    PLs[3*idx+0] = unphred[min(gl_arr_t[3*i+0], 255)];
									PLs[3*idx+1] = unphred[min(gl_arr_t[3*i+1], 255)];
									PLs[3*idx+2] = unphred[min(gl_arr_t[3*i+2], 255)];

									// grab genotypes

									// fill in DPs
									if (dp_arr_t[i]!= bcf_int32_missing){
									    DPs[idx] = dp_arr_t[i];
									}
									else DPs[idx] = 0;
									if (flip) {
									    float_swap = PLs[3*idx+2];
									    PLs[3*idx+2] = PLs[3*idx+0];
									    PLs[3*idx+0] = float_swap;
									}
								} else { PLs[3*idx+0] = -1.0f; PLs[3*idx+1] = -1.0f; PLs[3*idx+2] = -1.0f; }
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
								if (flip) {
								    float_swap = GPs[3*idx+2];
								    GPs[3*idx+2] = GPs[3*idx+0];
								    GPs[3*idx+0] = float_swap;
								    DSs[idx] = 2.0f - DSs[idx];
								}
							}
						}

						// Process variant
						if (grp_bin >= 0) {	// Do this variant fall within a given frequency bin?
							for (int i = 0 ; i < N ; i ++) {
								int true_genotype = getTruth(PLs[3*i+0], PLs[3*i+1], PLs[3*i+2], DPs[i]);

								// if true_genotype is NOT VALID (i.e. getTruth returns -1), then it skips the sample variant
								if (true_genotype >= 0) { // <--- TODO: OK, here is where some are getting filtered out.
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
									case 0:	genotype_bin_errors[3*grp_bin+0] += (esti_genotype != 0); genotype_bin_totals[3*grp_bin+0]++; break;
									case 1:	genotype_bin_errors[3*grp_bin+1] += (esti_genotype != 1); genotype_bin_totals[3*grp_bin+1]++; break;
									case 2:	genotype_bin_errors[3*grp_bin+2] += (esti_genotype != 2); genotype_bin_totals[3*grp_bin+2]++; break;
									}

									// [3] Update concordance per calibration bin
									switch (true_genotype) {
									case 0:	genotype_cal_errors[3*cal_bin+0] += (esti_genotype != 0); genotype_cal_totals[3*cal_bin+0]++; break;
									case 1:	genotype_cal_errors[3*cal_bin+1] += (esti_genotype != 1); genotype_cal_totals[3*cal_bin+1]++; break;
									case 2:	genotype_cal_errors[3*cal_bin+2] += (esti_genotype != 2); genotype_cal_totals[3*cal_bin+2]++; break;
									}

									// [4] Update Rsquare per bin
									rsquared_bin[grp_bin].push(DSs[i], true_genotype*1.0f);
									frequency_bin[grp_bin].push(maf);

									// [5] Update Rsquare per sample
									rsquared_spl[i].push(DSs[i], true_genotype*1.0f);

									// [6] Update data for AVG Rsquare
									sample_truth[i] = true_genotype*1.0f;
									sum_ds_e += DSs[i];
									sum_ds_t += true_genotype*1.0f;

									// Increment counts
									has_validation = true;
									ngenoval ++;
									ngenoval1 ++;
								}
                                else {ngenoval1 ++;}
							}
							// [6] Compute Rsquare variant and update data
							float xMean = mean(N, DSs);
							float xStdDev = stdDev(N, DSs, xMean);
							float yMean  = mean(N, sample_truth);
							float yStdDev = stdDev(N, sample_truth, yMean);

							if (xStdDev > 0.0f && yStdDev > 0.0f) {
                                avg_rsquared_bin[grp_bin].push(
                                        std::pow(pearson_prec(N, DSs, sample_truth, xMean, xStdDev, yMean, yStdDev),
                                                 2));
                            }
						}

					}
					// increment number of variants
					nvariantval += has_validation;
					nvarianttot ++;
				}
			}
		}
		free(af_ptr);
		free(gl_arr_t);
		free(ds_arr_e);
		free(gp_arr_e);
		bcf_sr_destroy(sr);
		n_variants_all_chromosomes += nvariantval;
		vrb.bullet("#variants in the overlap = " + stb.str(nvarianttot));
		vrb.bullet("#variants with validation data = " + stb.str(nvariantval));
		vrb.bullet("#genotypes used in validation = " + stb.str(ngenoval));
        vrb.bullet("#jeremy genotypes used in validation = " + stb.str(ngenoval1));
		vrb.bullet("%error rate in this file = " + stb.str(n_errors * 100.0f / ngenoval));
	}
	vrb.bullet("Total #variants = " + stb.str(n_variants_all_chromosomes));

	if (L == 0) {
		unsigned int n_found_groups = 0;
		for (map < string, pair < int, bool > > :: iterator itG = site2grp.begin() ; itG != site2grp.end() ; ++ itG) n_found_groups += itG->second.second;
		vrb.bullet("Total #variants in groups found = " + stb.str(n_found_groups));
	}
}
