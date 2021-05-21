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
	float sum = 0.0f;
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
	//vector < int > mappingT, mappingE;
	char** isec_samples_names = NULL;
	int* truth_imap;
	int* estimated_imap;
	int n_isec_samples_names = 0;

	string af_tag = options["af-tag"].as<string>();
	const bool gt_validation = options.count("gt-validation");
	const bool gt_target = options.count("gt-target");
	const bool use_alt_af = options.count("use-alt-af");
	int nthreads = options["thread"].as < int > ();
	fploidy=2;

	for (int f = 0 ; f < ftruth.size() ; f ++)
	{
		vrb.title("Reading set of input files [" + stb.str(f+1) + "/" + stb.str(ftruth.size()) + "]");
		vrb.bullet("Truth       [" + ftruth[f] + "]");
		vrb.bullet("Estimates   [" + festimated[f] + "]");
		vrb.bullet("Frequencies [" + ffrequencies[f] + "]");
		vrb.bullet("Region      [" + region[f] + "]");

		bcf_srs_t * sr =  bcf_sr_init();
		bcf_hdr_t *hdr_truth;
		bcf_hdr_t *hdr_estimated;

		sr->collapse = COLLAPSE_NONE;
		sr->require_index = 1;
		if (nthreads > 1) bcf_sr_set_threads(sr, nthreads-1);
		bcf_sr_set_regions(sr, region[f].c_str(), 0);

		//Opening files
		std::array<string,3> fnames = {ftruth[f],festimated[f],ffrequencies[f]};
		for (int reader_id=0; reader_id<3; ++reader_id)
		{
			if(!(bcf_sr_add_reader (sr, fnames[reader_id].c_str())))
			{
				if (sr->errnum != idx_load_failed) vrb.error("Failed to open file: " + fnames[reader_id] + "");
				bcf_sr_remove_reader (sr, reader_id);
				int ret = bcf_index_build3(fnames[reader_id].c_str(), NULL, 14, options["thread"].as < int > ());

				if (ret != 0)
				{
					if (ret == -2)
						vrb.error("index: failed to open " + fnames[reader_id]);
					else if (ret == -3)
						vrb.error("index: " + fnames[reader_id] + " is in a format that cannot be usefully indexed");
					else
						vrb.error("index: failed to create index for + " + fnames[reader_id]);
				}
				if(!(bcf_sr_add_reader (sr, fnames[reader_id].c_str()))) vrb.error("Problem opening/creating index file for [" + fnames[reader_id] + "]");
				else vrb.bullet("Index file for [" + fnames[reader_id] + "] has been successfully created.\n");
			}
		}

		// If first file; we initialize data structures
		if (f == 0)
		{
			std::unordered_set<string> truth_samples_names_set;
			std::unordered_set<string> isec_samples_names_set;

			int n_true_samples_file = bcf_hdr_nsamples(sr->readers[0].header);
			int n_esti_samples_file = bcf_hdr_nsamples(sr->readers[1].header);
			for (int i = 0 ; i < n_true_samples_file ; i ++) truth_samples_names_set.insert(sr->readers[0].header->samples[i]);
			for (int i = 0 ; i < n_esti_samples_file ; i ++)
				if (truth_samples_names_set.find(sr->readers[1].header->samples[i]) != truth_samples_names_set.end())
					isec_samples_names_set.insert(sr->readers[1].header->samples[i]);

			if (!options.count("samples")) subset_samples.assign(isec_samples_names_set.begin(),isec_samples_names_set.end());
			for (int i = 0; i < subset_samples.size() ; i ++)
			{
				if (isec_samples_names_set.find(subset_samples[i]) != isec_samples_names_set.end())
				{
					isec_samples_names = (char**) realloc(isec_samples_names, (n_isec_samples_names+1)*sizeof(const char*));
					isec_samples_names[n_isec_samples_names] = strdup(subset_samples[i].c_str());
					n_isec_samples_names++;
				}
			}
			if (n_isec_samples_names == 0) vrb.error("No sample in common between datasets");

			truth_imap = (int*)malloc(n_isec_samples_names * sizeof(int));
			estimated_imap = (int*)malloc(n_isec_samples_names * sizeof(int));

			hdr_truth = bcf_hdr_subset(sr->readers[0].header, n_isec_samples_names, isec_samples_names, truth_imap);
			hdr_estimated = bcf_hdr_subset(sr->readers[1].header, n_isec_samples_names, isec_samples_names, estimated_imap);

			if ( !hdr_truth || !hdr_estimated || bcf_hdr_nsamples(hdr_truth) != bcf_hdr_nsamples(hdr_estimated) || bcf_hdr_nsamples(hdr_truth)<=0) vrb.error("Error occurred while subsetting truth samples");

			n_true_samples = bcf_hdr_nsamples(hdr_truth);
			n_esti_samples = bcf_hdr_nsamples(hdr_estimated);

			N=n_true_samples;
			samples.reserve(N);
			for (int i=0; i<N;++i) samples.push_back(string(hdr_truth->samples[i]));

/*
			////////////////////////////////////////////////////////////////////////////////
			// Processing samples + overlap
			N = 0;
			n_true_samples = bcf_hdr_nsamples(hdr_truth);
			n_esti_samples = bcf_hdr_nsamples(hdr_estimated);

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
*/
			bcf_hrec_t * header_record = bcf_hdr_get_hrec(hdr_estimated, BCF_HL_GEN, "FPLOIDY", NULL, NULL);
			if (header_record == NULL) vrb.warning("Cannot retrieve FPLOIDY flag in VCF header [" + festimated[f] + "], used GLIMPSE version < 1.1.0? Assuming diploid genotypes [FPLOIDY=2].");
			else fploidy = atoi(header_record->value);
			if (fploidy_to_msg.find(fploidy) == fploidy_to_msg.end()) vrb.error("FPLOIDY out of bounds : " + stb.str(fploidy));
			vrb.bullet("FPLOIDY = "+ to_string(fploidy) + " [" + fploidy_to_msg[fploidy] + "]");

			//Ploidy
			ploidy_samples = std::vector<int> (N);
			ind2gpos = std::vector<int> (N);
			ploidy = std::abs(fploidy);
			n_haploid = 0;
			n_diploid = 0;
			const int ploidyP1 = ploidy + 1;
			if (fploidy > 0)
			{
				fploidy == 2 ? n_diploid = N: n_haploid = N;
				for (int i = 0 ; i < N ; i ++)
				{
					ploidy_samples[i] = fploidy;
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
				bcf_subset(hdr_estimated, line, n_isec_samples_names, estimated_imap);
				int ngt = bcf_get_genotypes(hdr_estimated, line, &gt_fields, &n_gt_fields);
				const int line_max_ploidy = ngt/n_esti_samples;
				assert(line_max_ploidy==ploidy); //we do not allow missing data
				int j=0;
				for(int i = 0 ; i < n_esti_samples; ++i)
				{
					//if (mappingE[i] == -1) continue;
					ind2gpos[j] = 3*n_diploid + 2*n_haploid;
					ploidy_samples[j] = 2 - (gt_fields[ploidy*i+1] == bcf_int32_vector_end);
					ploidy_samples[j] > 1? ++n_diploid : ++n_haploid;
					++j;
				}
				//assert(n_diploid > 0 && n_haploid > 0);
				bcf_sr_seek (sr, NULL, 0);

				if (n_diploid == 0 && n_haploid == 0) vrb.error("No samples found.");
				free(gt_fields);
			}

			// Allocating data structures
			genotype_spl_errors_all = vector < unsigned long int > (3*n_diploid + 2*n_haploid, 0);
			genotype_spl_totals_all = vector < unsigned long int > (3*n_diploid + 2*n_haploid, 0);
			genotype_cal_errors_all = vector < unsigned long int > (ploidyP1 * N_BIN_CAL, 0);
			genotype_cal_totals_all = vector < unsigned long int > (ploidyP1 * N_BIN_CAL, 0);
			rsquared_spl_ds_all = vector < stats2D > (N);
			rsquared_spl_gt_all = vector < stats2D > (N);

			//snps
			// Allocating data structures
			genotype_spl_errors_snps = vector < unsigned long int > (3*n_diploid + 2*n_haploid, 0);
			genotype_spl_totals_snps = vector < unsigned long int > (3*n_diploid + 2*n_haploid, 0);
			genotype_cal_errors_snps = vector < unsigned long int > (ploidyP1 * N_BIN_CAL, 0);
			genotype_cal_totals_snps = vector < unsigned long int > (ploidyP1 * N_BIN_CAL, 0);
			rsquared_spl_ds_snps = vector < stats2D > (N);
			rsquared_spl_gt_snps = vector < stats2D > (N);

			//indels
			// Allocating data structures
			genotype_spl_errors_indels = vector < unsigned long int > (3*n_diploid + 2*n_haploid, 0);
			genotype_spl_totals_indels = vector < unsigned long int > (3*n_diploid + 2*n_haploid, 0);
			genotype_cal_errors_indels = vector < unsigned long int > (ploidyP1 * N_BIN_CAL, 0);
			genotype_cal_totals_indels = vector < unsigned long int > (ploidyP1 * N_BIN_CAL, 0);
			rsquared_spl_ds_indels = vector < stats2D > (N);
			rsquared_spl_gt_indels = vector < stats2D > (N);

			if (L > 0) {
				genotype_bin_errors_all = vector < unsigned long int > (ploidyP1 * L, 0);
				genotype_bin_totals_all = vector < unsigned long int > (ploidyP1 * L, 0);
				rsquared_bin_ds_all = vector < stats2D > (L);
				rsquared_bin_gt_all = vector < stats2D > (L);
				rsquared_bin_ds_all = vector < stats2D > (L);
				rsquared_bin_gt_all = vector < stats2D > (L);
				frequency_bin_all = vector < stats1D > (L);

				genotype_bin_errors_snps = vector < unsigned long int > (ploidyP1 * L, 0);
				genotype_bin_totals_snps = vector < unsigned long int > (ploidyP1 * L, 0);
				rsquared_bin_ds_snps = vector < stats2D > (L);
				rsquared_bin_gt_snps = vector < stats2D > (L);
				frequency_bin_snps = vector < stats1D > (L);

				genotype_bin_errors_indels = vector < unsigned long int > (ploidyP1 * L, 0);
				genotype_bin_totals_indels = vector < unsigned long int > (ploidyP1 * L, 0);
				rsquared_bin_ds_indels = vector < stats2D > (L);
				rsquared_bin_gt_indels = vector < stats2D > (L);
				frequency_bin_indels = vector < stats1D > (L);

			} else {
				genotype_bin_errors_all = vector < unsigned long int > (ploidyP1 * rsquared_str.size(), 0);
				genotype_bin_totals_all = vector < unsigned long int > (ploidyP1 * rsquared_str.size(), 0);
				rsquared_bin_ds_all = vector < stats2D > (rsquared_str.size());
				rsquared_bin_gt_all = vector < stats2D > (rsquared_str.size());
				frequency_bin_all = vector < stats1D > (rsquared_str.size());

				genotype_bin_errors_snps = vector < unsigned long int > (ploidyP1 * rsquared_str.size(), 0);
				genotype_bin_totals_snps = vector < unsigned long int > (ploidyP1 * rsquared_str.size(), 0);
				rsquared_bin_ds_snps = vector < stats2D > (rsquared_str.size());
				rsquared_bin_gt_snps = vector < stats2D > (rsquared_str.size());
				frequency_bin_snps = vector < stats1D > (rsquared_str.size());

				genotype_bin_errors_indels = vector < unsigned long int > (ploidyP1 * rsquared_str.size(), 0);
				genotype_bin_totals_indels = vector < unsigned long int > (ploidyP1 * rsquared_str.size(), 0);
				rsquared_bin_ds_indels = vector < stats2D > (rsquared_str.size());
				rsquared_bin_gt_indels = vector < stats2D > (rsquared_str.size());
				frequency_bin_indels = vector < stats1D > (rsquared_str.size());
			}

			assert(n_diploid+n_haploid==N);
			string pl1  = n_haploid!=1? "s" : "";
			string pl2  = n_diploid!=1? "s" : "";
			vrb.bullet("#overlapping samples = " + stb.str(n_diploid+n_haploid) + " ["  + stb.str(n_haploid) + " haploid" + pl1 + "/ " + stb.str(n_diploid) + " diploid"+ pl2 + "]");
			if (n_diploid > 0 && n_haploid > 0)
					vrb.error("Samples with mixed ploidy detected. In order to get a correct value of the concordance metrics it is necessary to split the study samples by ploidy before running GLIMPSE_concordance or use the --samples option.");

			ploidy= n_haploid >0 ? 1 : 2;
		}
		else
		{
			hdr_truth = bcf_hdr_subset(sr->readers[0].header, n_isec_samples_names, isec_samples_names, truth_imap);
			hdr_estimated = bcf_hdr_subset(sr->readers[1].header, n_isec_samples_names, isec_samples_names, estimated_imap);
		}

		const int ploidyP1 = ploidy + 1;
		unsigned long nvarianttot = 0, nvariantvalall = 0, nvariantvalsnps = 0, nvariantvalindels = 0, nvariantvalfull=0, nset = 0, ngenoval_all = 0, n_errors_all = 0, n_errors_snps = 0, ngenoval_snps = 0, n_errors_indels = 0, ngenoval_indels = 0;
		unsigned long count_tg_err_2 = 0, count_tg_err_3 = 0, count_tg_err_4 = 0, count_tg_err_5 = 0, count_tg_err_1 = 0;
		unsigned long count_es_err_2 = 0, count_es_err_3 = 0, count_es_err_1 = 0;
		unsigned long afnobin = 0;

		int npl_t, ndp_t, npl_arr_t = 0, ndp_arr_t = 0, *pl_arr_t = NULL, *dp_arr_t = NULL;
		float *af_ptr=NULL,*ds_arr_e = NULL, *gp_arr_e = NULL, float_swap;
		int * gt_arr_e = NULL, ngt_arr_e = 0, ngt_e; //ngt_t is npl_t in case of --gt-validation
		int nds_e, nds_arr_e = 0;
		int naf_f, naf_arr_f=0;
		int ngp_e, ngp_arr_e = 0;
		int *itmp=NULL, mitmp=0, tret=0;

		vector < float > GLs = vector < float > ((3*n_diploid + 2*n_haploid), 0.0f);
		vector < float > DSs = vector < float > (N, 0.0f);
		vector < int > GTs = vector < int > (N, 0);
		vector < float > GPs = vector < float > ((3*n_diploid + 2*n_haploid), 0.0f);
		vector < int > DPs = vector < int > (N, 0);

		vector < float > sample_dosages(N);

		map < string, pair < int, bool > > :: iterator itG;
		bcf1_t * line_t, * line_e, * line_f;
		int nobi=0,noval=0;
		long int n_disc_var_truth = 0, n_disc_var_est = 0, n_disc_var_af = 0;

		while ((nset = bcf_sr_next_line (sr)))
		{
			if (nset == 3)
			{
				line_f =  bcf_sr_get_line(sr, 2);

				const int line_type = bcf_get_variant_types(line_f);
				/* LINE TYPES:
				 	 	line_type&VCF_REF -> No alts
				        line_type&VCF_SNP -> SNP
						line_type&VCF_INDEL -> INDEL
						line_type&VCF_MNP -> MNP
						line_type&VCF_OTHER -> other
				 */
				if (line_f->n_allele == 2 && ((line_type==VCF_SNP) || (line_type==VCF_INDEL)))
				{
					float af = 0.0f;

					line_t =  bcf_sr_get_line(sr, 0);
					line_e =  bcf_sr_get_line(sr, 1);

					bcf_subset(hdr_truth, line_t, n_isec_samples_names, truth_imap);
					bcf_subset(hdr_estimated, line_e, n_isec_samples_names, estimated_imap);

				    naf_f = bcf_get_info_float(sr->readers[2].header, line_f, af_tag.c_str(), &af_ptr, &naf_arr_f);

				    if (naf_f != 1)
				    {
			            int AC = -1, AN = 0;
			            tret = bcf_get_info_int32(sr->readers[2].header, line_f, "AN", &itmp, &mitmp);
			            if ( tret==1 )
			            {
			                AN = itmp[0];
			                tret = bcf_get_info_int32(sr->readers[2].header, line_f, "AC", &itmp, &mitmp);
			                if ( tret>0 )
			                    AC = itmp[0];
			            }
			            if ( AN>0 && AC>=0 )
			            {
			            	af = (double) AC / AN;
					    	naf_f=1;
					    }
				    } else af =  af_ptr[0];


				    if (gt_validation) npl_t = bcf_get_genotypes(hdr_truth, line_t, &pl_arr_t, &npl_arr_t);
				    else npl_t = bcf_get_format_int32(hdr_truth, line_t, "PL", &pl_arr_t, &npl_arr_t);
					ndp_t = bcf_get_format_int32(hdr_truth, line_t, "DP", &dp_arr_t, &ndp_arr_t);

					ngt_e = bcf_get_genotypes(hdr_estimated, line_e, &gt_arr_e, &ngt_arr_e);
					nds_e = bcf_get_format_float(hdr_estimated, line_e, "DS", &ds_arr_e, &nds_arr_e);
					ngp_e = bcf_get_format_float(hdr_estimated, line_e, "GP", &gp_arr_e, &ngp_arr_e);

					int n_validation_marker=0;
					if ((naf_f==1)&&(npl_t%n_true_samples==0)&&(gt_target||((nds_e==n_esti_samples)&&(ngp_e%n_esti_samples==0)))&&((ndp_t==n_true_samples) || D==0))
					{
						// Meta data for variant
						const int ploidy_gt_e_record=ngt_e/n_esti_samples;
						const int ploidy_e_record=ngp_e/n_esti_samples;
						const int ploidy_t_record=npl_t/n_true_samples;

						const bool flip = use_alt_af? false : (af > 0.5f);
						const float maf = use_alt_af? af : min(af, 1.0f - af);
						int grp_bin = -1;
						if (L>0) grp_bin = getFrequencyBin(maf);
						else
						{
							string chr = bcf_hdr_id2name(hdr_truth, line_t->rid);
							int pos = line_t->pos + 1;
							string uuid = chr + "_" + stb.str(pos);
							itG = site2grp.find(uuid);
							if (itG != site2grp.end()) {
								grp_bin = itG->second.first;
								itG->second.second = true;
							}
						}

						// Read Truth
						for(int i = 0 ; i < n_true_samples ; i ++)
						{
							const int index = ploidy_t_record*i;//ploidy
							const int gpos=ind2gpos[i];

							if (gt_validation)
							{
								const int pgt = (ploidy_samples[i]>1);

								if ( (pl_arr_t[index+0]==bcf_gt_missing || pl_arr_t[index+pgt]==bcf_gt_missing)
									|| (pl_arr_t[index+0]==bcf_int32_vector_end || pl_arr_t[index+pgt]==bcf_int32_vector_end))
								{
									GLs[gpos+0] = -1.0f;
									GLs[gpos+1] = -1.0f;
									GLs[gpos+pgt] = -1.0f;
								}
								else
								{
									int gt = (bcf_gt_allele(pl_arr_t[index+0])==1);
									if (ploidy_samples[i]>1) gt += (bcf_gt_allele(pl_arr_t[index+1])==1);
									GLs[gpos+0] = (gt == 0);
									GLs[gpos+1] = (gt == 1);
									GLs[gpos+ploidy] = gt == ploidy;
								}

							}
							else
							{
								if ( (pl_arr_t[index+0]==bcf_int32_missing || pl_arr_t[index+1]==bcf_int32_missing || pl_arr_t[index+ploidy]==bcf_int32_missing)
										|| (pl_arr_t[index+0]==bcf_int32_vector_end || pl_arr_t[index+1]==bcf_int32_vector_end || pl_arr_t[index+ploidy]==bcf_int32_vector_end)
										|| (pl_arr_t[index+0]<0) || (pl_arr_t[index+1]<0) || (pl_arr_t[index+ploidy]<0))
								{
									GLs[gpos+0] = -1.0f;
									GLs[gpos+1] = -1.0f;
									GLs[gpos+ploidy] = -1.0f;
								}
								else
								{
									GLs[gpos+0] = unphred[std::min(pl_arr_t[index+0],255)];
									GLs[gpos+1] = unphred[std::min(pl_arr_t[index+1],255)];
									GLs[gpos+ploidy] = unphred[std::min(pl_arr_t[index+ploidy],255)];
								}
							}
							DPs[i] = D > 0 ? dp_arr_t[i] : 0;
							if (flip) { float_swap = GLs[gpos+ploidy]; GLs[gpos+ploidy] = GLs[gpos+0]; GLs[gpos+0] = float_swap; }
						}

						// Read Estimates
						for(int i = 0 ; i < n_esti_samples ; i ++)
						{
							const int gpos=ind2gpos[i];
							const int pgt = (ploidy_samples[i]>1);
							const int index = ploidy_gt_e_record*i;//ploidy

							if ( (gt_arr_e[index+0]==bcf_gt_missing || gt_arr_e[index+pgt]==bcf_gt_missing)
								|| (gt_arr_e[index+0]==bcf_int32_vector_end || gt_arr_e[index+pgt]==bcf_int32_vector_end))
							{
								GTs[i] = -1;
							}
							else
							{
								int gt = (bcf_gt_allele(gt_arr_e[index+0])==1);
								if (ploidy_samples[i]>1) gt += (bcf_gt_allele(gt_arr_e[index+1])==1);
								GTs[i] = gt;
							}
							if (gt_target)
							{
								DSs[i] = GTs[i];

								GPs[gpos+0] = (GTs[i] == 0);
								GPs[gpos+1] = (GTs[i] == 1);
								GPs[gpos+ploidy] = (GTs[i] == ploidy);
							}
							else
							{
								DSs[i] = ds_arr_e[i];

								GPs[gpos+0] = gp_arr_e[ploidy_e_record*i+0];
								GPs[gpos+1] = gp_arr_e[ploidy_e_record*i+1];
								GPs[gpos+ploidy] = gp_arr_e[ploidy_e_record*i+ploidy];
							}
							if (flip) { float_swap = GPs[gpos+ploidy]; GPs[gpos+ploidy] = GPs[gpos+0]; GPs[gpos+0] = float_swap; DSs[i] = ((float) ploidy) - DSs[i]; GTs[i] = ploidy - GTs[i];}
						}

						// Process variant
						if (grp_bin >= 0)
						{	// Do this variant fall within a given frequency bin?
							for (int i = 0 ; i < N ; i ++)
							{
								const int gpos=ind2gpos[i];
								const int true_genotype = getTruth(GLs[gpos+0], GLs[gpos+1], GLs[gpos+ploidy], DPs[i], ploidy_samples[i]);
								if (true_genotype >= 0)
								{
									int esti_genotype;
									if (gt_target) esti_genotype = getMostLikelyGT(GTs[i],DSs[i], ploidy_samples[i]);
									else esti_genotype = getMostLikely(GPs[gpos+0], GPs[gpos+1], GPs[gpos+ploidy], DSs[i], ploidy_samples[i]);

									if (esti_genotype >= 0)
									{
										const int cal_bin = getCalibrationBin(GPs[gpos+0], GPs[gpos+1], GPs[gpos+ploidy]);
										const int is_error = (true_genotype != esti_genotype);
										// [0] Overall concordance for verbose
										n_errors_all += is_error;
										ngenoval_all ++;
										// [1] Update concordance per sample
										genotype_spl_errors_all[gpos+true_genotype] += is_error;
										genotype_spl_totals_all[gpos+true_genotype]++;
										// [2] Update concordance per bin
										genotype_bin_errors_all[ploidyP1*grp_bin+true_genotype] += is_error;
										genotype_bin_totals_all[ploidyP1*grp_bin+true_genotype]++;
										// [3] Update concordance per calibration bin
										genotype_cal_errors_all[ploidyP1*cal_bin+true_genotype] += is_error;
										genotype_cal_totals_all[ploidyP1*cal_bin+true_genotype]++;
										// [4] Update Rsquare per bin, DS and best-guess
										rsquared_bin_ds_all[grp_bin].push(DSs[i], true_genotype*1.0f);
										rsquared_bin_gt_all[grp_bin].push(esti_genotype*1.0f, true_genotype*1.0f);
										frequency_bin_all[grp_bin].push(maf);
										// [5] Update Rsquare per sample using DS
										rsquared_spl_ds_all[i].push(DSs[i], true_genotype*1.0f);
										// [6] Update Rsquare per sample using estimated genotype
										rsquared_spl_gt_all[i].push(esti_genotype*1.0f, true_genotype*1.0f);

										if (line_type==VCF_SNP)
										{
											n_errors_snps += is_error;
											ngenoval_snps ++;
											genotype_spl_errors_snps[gpos+true_genotype] += is_error;
											genotype_spl_totals_snps[gpos+true_genotype]++;
											genotype_bin_errors_snps[ploidyP1*grp_bin+true_genotype] += is_error;
											genotype_bin_totals_snps[ploidyP1*grp_bin+true_genotype]++;
											genotype_cal_errors_snps[ploidyP1*cal_bin+true_genotype] += is_error;
											genotype_cal_totals_snps[ploidyP1*cal_bin+true_genotype]++;
											rsquared_bin_ds_snps[grp_bin].push(DSs[i], true_genotype*1.0f);
											rsquared_bin_gt_snps[grp_bin].push(esti_genotype*1.0f, true_genotype*1.0f);
											frequency_bin_snps[grp_bin].push(maf);
											rsquared_spl_ds_snps[i].push(DSs[i], true_genotype*1.0f);
											rsquared_spl_gt_snps[i].push(esti_genotype*1.0f, true_genotype*1.0f);
										}
										if (line_type==VCF_INDEL)
										{
											n_errors_indels += is_error;
											ngenoval_indels ++;
											genotype_spl_errors_indels[gpos+true_genotype] += is_error;
											genotype_spl_totals_indels[gpos+true_genotype]++;
											genotype_bin_errors_indels[ploidyP1*grp_bin+true_genotype] += is_error;
											genotype_bin_totals_indels[ploidyP1*grp_bin+true_genotype]++;
											genotype_cal_errors_indels[ploidyP1*cal_bin+true_genotype] += is_error;
											genotype_cal_totals_indels[ploidyP1*cal_bin+true_genotype]++;
											rsquared_bin_ds_indels[grp_bin].push(DSs[i], true_genotype*1.0f);
											rsquared_bin_gt_indels[grp_bin].push(esti_genotype*1.0f, true_genotype*1.0f);
											frequency_bin_indels[grp_bin].push(maf);
											rsquared_spl_ds_indels[i].push(DSs[i], true_genotype*1.0f);
											rsquared_spl_gt_indels[i].push(esti_genotype*1.0f, true_genotype*1.0f);
										}
										// Increment counts
										n_validation_marker++;
									}
									else
									{
										switch (esti_genotype) {
										case -2:  count_es_err_2++; break;
										case -3:  count_es_err_3++; break;
										default:  count_es_err_1++; break;
										}
									}

								}
								else
								{
									switch (true_genotype) {
									case -2:  count_tg_err_2++; break;
									case -3:  count_tg_err_3++; break;
									case -4:  count_tg_err_4++; break;
									case -5:  count_tg_err_5++; break;
									default:
										count_tg_err_1++; break;
									}
								}
							}
						} else afnobin++;
					}
					else noval++;

					// increment number of variants
					nvariantvalall += n_validation_marker>0;
					nvariantvalfull += n_validation_marker==N;

					if (line_type==VCF_SNP && n_validation_marker) nvariantvalsnps++;
					if (line_type==VCF_INDEL && n_validation_marker) nvariantvalindels++;

					nvarianttot ++;
				}
				else nobi++;
			}
			else
			{
				if (bcf_sr_has_line(sr,0))
					++n_disc_var_truth;
				if (bcf_sr_has_line(sr,1))
					++n_disc_var_est;
				if (bcf_sr_has_line(sr,2))
					++n_disc_var_af;
			}
		}
		free(af_ptr);
		free(pl_arr_t);
		free(gt_arr_e);
		free(ds_arr_e);
		free(gp_arr_e);
		bcf_hdr_destroy(hdr_truth);
		bcf_hdr_destroy(hdr_estimated);
		bcf_sr_destroy(sr);
		n_variants_all_chromosomes += nvariantvalall;
		if (afnobin > 0) vrb.bullet("#variants discarded as no associated bin can be set (e.g. AF<=0 or AF>1) = " + stb.str(afnobin));
		if (nvarianttot==0) vrb.error("No variant found in the intersection of files. Files are probably not aligned correctly. Please verify that chromosome names and regions are matching for the imputed, validation and allele frequency file.");
		if (nvariantvalall==0) vrb.error("No usable validation variant has been found in the intersection of files. Verify that the validation file has FORMAT/PL and FORMAT/DP fields defined, the imputed file has FORMAT/DS and FORMAT/GP fields defined in the same region and the allele frequency file has the INFO/AF (or --af-tag) field defined at every marker.");

		vrb.print("");
		vrb.bullet("#variants in the overlap (biallelic SNPs and indels) = " + stb.str(nvarianttot));
		vrb.bullet("#variants with all genotypes in the validation data = " + stb.str(nvariantvalfull));
		vrb.bullet("#variants with at least one genotype in the validation data = " + stb.str(nvariantvalall) + " [SNPs = " + stb.str(nvariantvalsnps) + " (" + stb.str(nvariantvalsnps*100.0f/nvariantvalall) + "%), indels = " + stb.str(nvariantvalindels)  + " (" + stb.str(nvariantvalindels*100.0f/nvariantvalall) + "%)]");

		string ra,rs,ri;
		rs = ngenoval_all?stb.str(ngenoval_snps*100.0f/ngenoval_all):"-";
		ri = ngenoval_all?stb.str(ngenoval_indels*100.0f/ngenoval_all):"-";
		vrb.bullet("#genotypes used in validation = " + stb.str(ngenoval_all) + " [#SNPs = " + stb.str(ngenoval_snps) + " (" + rs + "%), indels = " + stb.str(ngenoval_indels)  + " (" + ri + "%)]");

		vrb.print("");
		vrb.bullet("Statistics on discarded true genoypes:");
		vrb.print("     #FORMAT/DP missing: " + stb.str(count_tg_err_5));
		vrb.print("     #FORMAT/DP < MinDP: " + stb.str(count_tg_err_2));
		if (gt_validation)
		{
			vrb.print("     #Missing/malformatted GTs: " + stb.str(count_tg_err_3));
		}
		else
		{
			vrb.print("     #Missing/malformatted PLs: " + stb.str(count_tg_err_3));
			vrb.print("     #Max(GLs) < MinPROB: " + stb.str(count_tg_err_4));
			vrb.print("     #Other (e.g. PLs have the same value): " + stb.str(count_tg_err_1));
		}

		vrb.print("");
		vrb.bullet("Statistics on discarded dosages / genotype probabilities:");
		vrb.print("     #Missing/malformatted DSs: " + stb.str(count_es_err_2));
		if (gt_target) vrb.print("     #Missing/malformatted GTs: " + stb.str(count_es_err_3));
		else vrb.print("     #Missing/malformatted GPs: " + stb.str(count_es_err_3));
		vrb.print("     #Other (e.g. no unique best-guess genotype): " + stb.str(count_es_err_1));

		vrb.print("");
		ra = ngenoval_all?stb.str(n_errors_all*100.0f/ngenoval_all):"-";
		rs = ngenoval_snps?stb.str(n_errors_snps*100.0f/ngenoval_snps):"-";
		ri = ngenoval_indels?stb.str(n_errors_indels*100.0f/ngenoval_indels):"-";

		vrb.bullet("Error rate in this file = " + ra + "% [SNPs = " + rs + "%, indels = " + ri + "%)]");

	}
	free(truth_imap);
	free(estimated_imap);
	for(int i = 0; i < n_isec_samples_names; i++) free(isec_samples_names[i]);

	vrb.print("");
	vrb.bullet("Total #variants = " + stb.str(n_variants_all_chromosomes));

	if (L == 0) {
		unsigned int n_found_groups = 0;
		for (map < string, pair < int, bool > > :: iterator itG = site2grp.begin() ; itG != site2grp.end() ; ++ itG) n_found_groups += itG->second.second;
		vrb.bullet("Total #variants in groups found = " + stb.str(n_found_groups));
	}
}
