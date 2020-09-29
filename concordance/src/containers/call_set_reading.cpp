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

std::map<std::string, int> mapPloidy = {
		{"1",1},
		{"2",2}
};

std::map<int, string> fploidy_to_msg = {
		{-2, "Mixed haploid/diploid samples in the region"},
		{1,"Only haploid samples in the region"},
		{2,"Only diploid samples in the region"}
};

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
    for (int i = 0; i < n; i++)
        sum += (x[i] - xMean) * (y[i] - yMean);

    return sum / (n * xStdDev * yStdDev);
}

void call_set::readData(vector < string > & ftruth, vector < string > & festimated, vector < string > & ffrequencies, vector < string > & region, bpo::variables_map& options) {
	tac.clock();
	int n_true_samples, n_esti_samples;
	unsigned long int n_variants_all_chromosomes = 0;
	vector < int > mappingT, mappingE;

	string af_tag = options["af-tag"].as<string>();
	int nthreads = options["thread"].as < int > ();
	fploidy=2;

	for (int f = 0 ; f < ftruth.size() ; f ++) {
		vrb.title("Reading set of input files [" + stb.str(f+1) + "/" + stb.str(ftruth.size()) + "]");
		vrb.bullet("Truth       [" + ftruth[f] + "]");
		vrb.bullet("Estimates   [" + festimated[f] + "]");
		vrb.bullet("Frequencies [" + ffrequencies[f] + "]");
		vrb.bullet("Region      [" + region[f] + "]");

		bcf_srs_t * sr =  bcf_sr_init();

		sr->collapse = COLLAPSE_NONE;
		sr->require_index = 1;
		if (nthreads > 1) bcf_sr_set_threads(sr, nthreads);
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

			bcf_hrec_t * header_record = bcf_hdr_get_hrec(sr->readers[1].header, BCF_HL_GEN, "FPLOIDY", NULL, NULL);
			if (header_record == NULL) vrb.warning("Cannot retrieve FPLOIDY flag in VCF header [" + festimated[f] + "], used a version of GLIMPSE < 1.1.0? Assuming diploid genotypes only [FPLOIDY=2].");
			else fploidy = atoi(header_record->value);
			if (fploidy_to_msg.find(fploidy) == fploidy_to_msg.end()) vrb.error("FPLOIDY out of bounds : " + stb.str(fploidy));
			vrb.bullet("FPLOIDY = "+ to_string(fploidy) + " [" + fploidy_to_msg[fploidy] + "]");

			//Ploidy
			ploidy = std::vector<int> (N);
			ind2gpos = std::vector<int> (N);
			max_ploidy = std::abs(fploidy);
			n_haploid = 0;
			n_diploid = 0;
			if (fploidy > 0)
			{
				fploidy == 2 ? n_diploid = N: n_haploid = N;
				for (int i = 0 ; i < N ; i ++)
				{
					ploidy[i] = fploidy;
					ind2gpos[i] = (fploidy+1)*i;
				}
			}
			else
			{
				int nset=0;
				while (!bcf_sr_has_line(sr,1)) bcf_sr_next_line (sr);

				if (!bcf_sr_has_line(sr,1)) vrb.error("No marker found in the intersection for the imputed file.");

				int * gt_fields = NULL;
				int n_gt_fields = 0;

				bcf1_t * line =  bcf_sr_get_line(sr, 1);
				int ngt = bcf_get_genotypes(sr->readers[1].header, line, &gt_fields, &n_gt_fields);
				const int line_max_ploidy = ngt/n_esti_samples;
				assert(line_max_ploidy==max_ploidy); //we do not allow missing data
				int j=0;
				for(int i = 0 ; i < n_esti_samples; ++i)
				{
					if (mappingE[i] == -1) continue;
					ind2gpos[j] = 3*n_diploid + 2*n_haploid;
					ploidy[j] = 2 - (gt_fields[max_ploidy*i+1] == bcf_int32_vector_end);
					ploidy[j] > 1? ++n_diploid : ++n_haploid;
					++j;
				}
				//assert(n_diploid > 0 && n_haploid > 0);
				bcf_sr_seek (sr, NULL, 0);

				if (n_diploid == 0 && n_haploid == 0)
					vrb.error("No samples found.");
			}

			// Allocating data structures
			genotype_spl_errors = vector < unsigned long int > (3*n_diploid + 2*n_haploid, 0);
			genotype_spl_totals = vector < unsigned long int > (3*n_diploid + 2*n_haploid, 0);
			genotype_cal_errors = vector < unsigned long int > ((max_ploidy+1) * N_BIN_CAL, 0);
			genotype_cal_totals = vector < unsigned long int > ((max_ploidy+1) * N_BIN_CAL, 0);
			rsquared_spl = vector < stats2D > (N);

			if (L > 0) {
				genotype_bin_errors = vector < unsigned long int > ((max_ploidy+1) * L, 0);
				genotype_bin_totals = vector < unsigned long int > ((max_ploidy+1) * L, 0);
				rsquared_bin = vector < stats2D > (L);
				frequency_bin = vector < stats1D > (L);
				avg_rsquared_bin = vector < stats1D > (L);
			} else {
				genotype_bin_errors = vector < unsigned long int > ((max_ploidy+1) * rsquared_str.size(), 0);
				genotype_bin_totals = vector < unsigned long int > ((max_ploidy+1) * rsquared_str.size(), 0);
				rsquared_bin = vector < stats2D > (rsquared_str.size());
				frequency_bin = vector < stats1D > (rsquared_str.size());
				avg_rsquared_bin = vector < stats1D > (rsquared_str.size());
			}

			assert(n_diploid+n_haploid==N);
			vrb.bullet("#overlapping samples = " + stb.str(n_diploid+n_haploid) + " ["  + stb.str(n_haploid) + " haploids / " + stb.str(n_diploid) + " diploids]");

			if (n_diploid > 0 && n_haploid > 0)
					vrb.warning("Samples with mixed ploidy detected. In order to get a correct value of the concordance metrics it is NECESSARY to split the study samples by ploidy before running GLIMPSE_concordance or use the --samples option.");

		}

		unsigned long nvarianttot = 0, nvariantval = 0, nset = 0, ngenoval = 0, n_errors = 0;
		int ngl_t, ndp_t, ngl_arr_t = 0, ndp_arr_t = 0, *gl_arr_t = NULL, *dp_arr_t = NULL;
		float *af_ptr=NULL,*ds_arr_e = NULL, *gp_arr_e = NULL, float_swap;
		int nds_e, nds_arr_e = 0;
		float naf_f;
		int naf_arr_f=0, ngp_e, ngp_arr_e = 0;
		vector < double > PLs = vector < double > ((3*n_diploid + 2*n_haploid), 0.0f);
		vector < float > DSs = vector < float > (N, 0.0f);
		vector < float > GPs = vector < float > ((3*n_diploid + 2*n_haploid), 0.0f);
		vector < int > DPs = vector < int > (N, 0);

		vector < float > sample_dosages(N);
		vector < float > sample_truth(N);
		float sum_ds_e = 0.0f;
		float sum_ds_t = 0.0f;

		map < string, pair < int, bool > > :: iterator itG;
		bcf1_t * line_t, * line_e, * line_f;

		while ((nset = bcf_sr_next_line (sr))) {
			if (nset == 3) {
				line_f =  bcf_sr_get_line(sr, 2);
				if (line_f->n_allele == 2) {
					line_t =  bcf_sr_get_line(sr, 0);
					line_e =  bcf_sr_get_line(sr, 1);

				    naf_f = bcf_get_info_float(sr->readers[2].header, line_f, af_tag.c_str(), &af_ptr, &naf_arr_f);
					ngl_t = bcf_get_format_int32(sr->readers[0].header, line_t, "PL", &gl_arr_t, &ngl_arr_t);
					ndp_t = bcf_get_format_int32(sr->readers[0].header, line_t, "DP", &dp_arr_t, &ndp_arr_t);
					nds_e = bcf_get_format_float(sr->readers[1].header, line_e, "DS", &ds_arr_e, &nds_arr_e);
					ngp_e = bcf_get_format_float(sr->readers[1].header, line_e, "GP", &gp_arr_e, &ngp_arr_e);
					
					sum_ds_e = 0.0f;
					sum_ds_t = 0.0f;
					bool has_validation = false;
					if ((naf_f==1)&&(ngl_t%n_true_samples==0)&&(nds_e==n_esti_samples)&&(ngp_e%n_esti_samples==0)&&(ndp_t==n_true_samples)) {

						// Meta data for variant
						int ploidy_e_record=ngp_e/n_esti_samples;
						int ploidy_t_record=ngl_t/n_true_samples;

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
								if (ploidy[idx] == 1)
								{
									int gpos=ind2gpos[idx];
									if (gl_arr_t[ploidy_t_record*i+0] != bcf_int32_missing && gl_arr_t[ploidy_t_record*i+1] != bcf_int32_missing) {
										PLs[gpos+0] = unphred[min(gl_arr_t[ploidy_t_record*i+0], 255)];
										PLs[gpos+1] = unphred[min(gl_arr_t[ploidy_t_record*i+1], 255)];
										if (dp_arr_t[i]!= bcf_int32_missing) DPs[idx] = dp_arr_t[i];
										else DPs[idx] = 0;
										if (flip) { float_swap = PLs[gpos+1]; PLs[gpos+1] = PLs[gpos+0]; PLs[gpos+0] = float_swap; }
									} else { PLs[gpos+0] = -1.0f; PLs[gpos+1] = -1.0f; }
								}
								else
								{
									int gpos=ind2gpos[idx];
									if (gl_arr_t[ploidy_t_record*i+0] != bcf_int32_missing && gl_arr_t[ploidy_t_record*i+1] != bcf_int32_missing && gl_arr_t[ploidy_t_record*i+2] != bcf_int32_missing) {
										PLs[gpos+0] = unphred[min(gl_arr_t[ploidy_t_record*i+0], 255)];
										PLs[gpos+1] = unphred[min(gl_arr_t[ploidy_t_record*i+1], 255)];
										PLs[gpos+2] = unphred[min(gl_arr_t[ploidy_t_record*i+2], 255)];
										if (dp_arr_t[i]!= bcf_int32_missing) DPs[idx] = dp_arr_t[i];
										else DPs[idx] = 0;
										if (flip) { float_swap = PLs[gpos+2]; PLs[gpos+2] = PLs[gpos+0]; PLs[gpos+0] = float_swap; }
									} else { PLs[gpos+0] = -1.0f; PLs[gpos+1] = -1.0f; PLs[gpos+2] = -1.0f; }
								}
							}
						}

						// Read Estimates
						for(int i = 0 ; i < n_esti_samples ; i ++) {
							int idx = mappingE[i];
							if (idx >= 0) {
								int gpos=ind2gpos[idx];
								DSs[idx] = ds_arr_e[i];
								if (ploidy[idx] == 1)
								{
									GPs[gpos+0] = std::round(gp_arr_e[ploidy_e_record*i+0] * 100.0f) / 100.0f;
									GPs[gpos+1] = std::round(gp_arr_e[ploidy_e_record*i+1] * 100.0f) / 100.0f;
									if (flip) { float_swap = GPs[gpos+1]; GPs[gpos+1] = GPs[gpos+0]; GPs[gpos+0] = float_swap; DSs[idx] = 1.0f - DSs[idx]; }
								}
								else
								{
									GPs[gpos+0] = std::round(gp_arr_e[ploidy_e_record*i+0] * 100.0f) / 100.0f;
									GPs[gpos+1] = std::round(gp_arr_e[ploidy_e_record*i+1] * 100.0f) / 100.0f;
									GPs[gpos+2] = std::round(gp_arr_e[ploidy_e_record*i+2] * 100.0f) / 100.0f;

									if (flip) { float_swap = GPs[gpos+2]; GPs[gpos+2] = GPs[gpos+0]; GPs[gpos+0] = float_swap; DSs[idx] = 2.0f - DSs[idx]; }
								}
							}
						}

						// Process variant
						if (grp_bin >= 0) {	// Do this variant fall within a given frequency bin?
							for (int i = 0 ; i < N ; i ++) {
								int gpos=ind2gpos[i];
								int true_genotype = getTruth(PLs[gpos+0], PLs[gpos+1], ploidy[i] == 1 ? 0.0 : PLs[gpos+2], DPs[i]);
								if (true_genotype >= 0) {
									int esti_genotype = getMostLikely(GPs[gpos+0], GPs[gpos+1], ploidy[i] == 1? 0.0 : GPs[gpos+2]);
									int cal_bin = getCalibrationBin(GPs[gpos+0], GPs[gpos+1], ploidy[i] == 1 ? 0.0 : GPs[gpos+2]);

									// [0] Overall concordance for verbose
									n_errors += (true_genotype != esti_genotype);

									// [1] Update concordance per sample
									switch (true_genotype) {
									case 0: genotype_spl_errors[gpos+0] += (esti_genotype != 0); genotype_spl_totals[gpos+0]++; break;
									case 1: genotype_spl_errors[gpos+1] += (esti_genotype != 1); genotype_spl_totals[gpos+1]++; break;
									case 2: genotype_spl_errors[gpos+2] += (esti_genotype != 2); genotype_spl_totals[gpos+2]++; break;
									}

									// [2] Update concordance per bin
									switch (true_genotype) {
									case 0:	genotype_bin_errors[(max_ploidy+1)*grp_bin+0] += (esti_genotype != 0); genotype_bin_totals[(max_ploidy+1)*grp_bin+0]++; break;
									case 1:	genotype_bin_errors[(max_ploidy+1)*grp_bin+1] += (esti_genotype != 1); genotype_bin_totals[(max_ploidy+1)*grp_bin+1]++; break;
									case 2:	genotype_bin_errors[(max_ploidy+1)*grp_bin+2] += (esti_genotype != 2); genotype_bin_totals[(max_ploidy+1)*grp_bin+2]++; break;
									}

									// [3] Update concordance per calibration bin
									switch (true_genotype) {
									case 0:	genotype_cal_errors[(max_ploidy+1)*cal_bin+0] += (esti_genotype != 0); genotype_cal_totals[(max_ploidy+1)*cal_bin+0]++; break;
									case 1:	genotype_cal_errors[(max_ploidy+1)*cal_bin+1] += (esti_genotype != 1); genotype_cal_totals[(max_ploidy+1)*cal_bin+1]++; break;
									case 2:	genotype_cal_errors[(max_ploidy+1)*cal_bin+2] += (esti_genotype != 2); genotype_cal_totals[(max_ploidy+1)*cal_bin+2]++; break;
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
								}
							}
							// [6] Compute Rsquare variant and update data
							float xMean = mean(N, DSs);
							float xStdDev = stdDev(N, DSs, xMean);
							float yMean  = mean(N, sample_truth);
							float yStdDev = stdDev(N, sample_truth, yMean);

							if (xStdDev > 0.0f && yStdDev > 0.0f)
								avg_rsquared_bin[grp_bin].push(std::pow(pearson_prec(N,DSs,sample_truth, xMean, xStdDev, yMean, yStdDev),2));
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
		if (nvarianttot==0) vrb.error("No variant found in the intersection of files. Files are probably not aligned correctly. Please verify that chromosome names and regions are matching for the imputed, validation and allele frequency file.");
		if (nvariantval==0) vrb.error("No usuable validation variant has been found in the intersection of files. Verify that the validation file has PL and DP fields defined and the imputed file has DS and GP fields defined in the same region.");
		vrb.bullet("#genotypes used in validation = " + stb.str(ngenoval));
		vrb.bullet("%error rate in this file = " + stb.str(n_errors * 100.0f / ngenoval));
	}
	vrb.bullet("Total #variants = " + stb.str(n_variants_all_chromosomes));

	if (L == 0) {
		unsigned int n_found_groups = 0;
		for (map < string, pair < int, bool > > :: iterator itG = site2grp.begin() ; itG != site2grp.end() ; ++ itG) n_found_groups += itG->second.second;
		vrb.bullet("Total #variants in groups found = " + stb.str(n_found_groups));
	}
}
