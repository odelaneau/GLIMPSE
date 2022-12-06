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

#include "call_set_header.h"

std::map<std::string, int> mapPloidy = {
		{"1",1},
		{"2",2}
};

std::map<int, std::string> fploidy_to_msg = {
		{-2, "Mixed haploid/diploid samples in the region"},
		{1,"Only haploid samples in the region"},
		{2,"Only diploid samples in the region"}
};


void call_set::readData(std::vector < std::string > & ftruth, std::vector < std::string > & festimated, std::vector < std::string > & ffrequencies, std::vector < std::string > & region, bpo::variables_map& options, const float gp_filter, const std::string out_filename) {
	tac.clock();
	int n_true_samples, n_esti_samples;
	unsigned long int n_variants_all_chromosomes = 0;
	//std::vector < int > mappingT, mappingE;
	char** isec_samples_names = NULL;
	int* truth_imap;
	int* estimated_imap;
	int n_isec_samples_names = 0;

	std::string af_tag = options["af-tag"].as<std::string>();
	const bool gt_validation = options.count("gt-val");
	const bool gt_target = options.count("gt-tar");
	const bool use_gp_filter = options.count("min-tar-gp");
	const bool use_alt_af = options.count("use-alt-af");
	const bool out_r2_per_site = options.count("out-r2-per-site");
	const bool out_rej_sites = options.count("out-rej-sites");
	const bool out_conc_sites =  options.count("out-conc-sites");
	const bool out_disc_sites =  options.count("out-disc-sites");
	int nthreads = options["threads"].as < int > ();
	fploidy=2;

	stats2D r2_variant;
	output_file out_file_r2_sites("");

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
		std::array<std::string,3> fnames = {ftruth[f],festimated[f],ffrequencies[f]};
		for (int reader_id=0; reader_id<3; ++reader_id)
		{
			if(!(bcf_sr_add_reader (sr, fnames[reader_id].c_str())))
			{
				if (sr->errnum != idx_load_failed) vrb.error("Failed to open file: " + fnames[reader_id] + "");
				bcf_sr_remove_reader (sr, reader_id);
				int ret = bcf_index_build3(fnames[reader_id].c_str(), NULL, 14, options["threads"].as < int > ());

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
			std::unordered_set<std::string> truth_samples_names_set;
			std::unordered_set<std::string> isec_samples_names_set;

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
			for (int i=0; i<N;++i) samples.push_back(std::string(hdr_truth->samples[i]));

/*
			////////////////////////////////////////////////////////////////////////////////
			// Processing samples + overlap
			N = 0;
			n_true_samples = bcf_hdr_nsamples(hdr_truth);
			n_esti_samples = bcf_hdr_nsamples(hdr_estimated);

			mappingT = std::vector < int > (n_true_samples, -1);
			mappingE = std::vector < int > (n_esti_samples, -1);
			for (int i = 0 ; i < n_true_samples ; i ++) {
				std::string ts = std::string(sr->readers[0].header->samples[i]);
				if (!use_subset_samples || subset_samples.count(ts)>0) {
					for (int j = 0 ; j < n_esti_samples ; j ++) {
						std::string es = std::string(sr->readers[1].header->samples[j]);
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
			vrb.bullet("FPLOIDY = "+ std::to_string(fploidy) + " [" + fploidy_to_msg[fploidy] + "]");

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

			if (options.count("allele-counts") || options.count("ac-bins"))
			{
				std::vector<int> ac_bin = {4,16,80,400,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000,2000000,5000000,10000000,20000000,50000000,100000000};
				if (options.count("ac-bins")) ac_bin = options["ac-bins"].as < std::vector < int > > ();

				int nset=0;
				while (! bcf_sr_has_line(sr,2)) bcf_sr_next_line (sr);
				int *itmp=NULL, mitmp=0, tret=0;
				tret = bcf_get_info_int32(sr->readers[2].header, bcf_sr_get_line(sr,2), "AN", &itmp, &mitmp);
			    if ( tret==1 && itmp[0] > 0)
			    {
			    	int AN = itmp[0] / ploidy;
			    	double delta = 0.5/itmp[0];
			    	bins.clear();
			    	bins.push_back(0);
			    	for (int k=0; k<ac_bin.size() && ac_bin[k] <= AN;++k)
			    	{
			    		bins.push_back(ac_bin[k]*1.0/itmp[0]+delta);
			    	}
			    	if (0.5 - bins.back() >= 0.1) bins.push_back(0.5+delta);
			    	else bins.back()=0.5+delta;

			    	L=bins.size()-1;
			    } else vrb.error("ERROR reading AN field required for --allele-counts option");
				bcf_sr_seek (sr, NULL, 0);
				if (itmp) free(itmp);
			}

			// Allocating data structures
			genotype_spl_errors_all = std::vector < unsigned long int > (3*n_diploid + 2*n_haploid, 0);
			genotype_spl_totals_all = std::vector < unsigned long int > (3*n_diploid + 2*n_haploid, 0);
			genotype_val_spl_totals_all = std::vector < unsigned long int > (3*n_diploid + 2*n_haploid, 0);
			filtered_gp_spl_all = std::vector < unsigned long int > (n_diploid + n_haploid, 0);
			genotype_cal_errors_all = std::vector < unsigned long int > (ploidyP1 * N_BIN_CAL, 0);
			genotype_cal_totals_all = std::vector < unsigned long int > (ploidyP1 * N_BIN_CAL, 0);
			rsquared_spl_ds_all = std::vector < stats2D > (N);
			rsquared_spl_gt_all = std::vector < stats2D > (N);

			//snps
			// Allocating data structures
			genotype_spl_errors_snps = std::vector < unsigned long int > (3*n_diploid + 2*n_haploid, 0);
			genotype_spl_totals_snps = std::vector < unsigned long int > (3*n_diploid + 2*n_haploid, 0);
			genotype_val_spl_totals_snps = std::vector < unsigned long int > (3*n_diploid + 2*n_haploid, 0);
			filtered_gp_spl_snps = std::vector < unsigned long int > (n_diploid +n_haploid, 0);
			genotype_cal_errors_snps = std::vector < unsigned long int > (ploidyP1 * N_BIN_CAL, 0);
			genotype_cal_totals_snps = std::vector < unsigned long int > (ploidyP1 * N_BIN_CAL, 0);
			rsquared_spl_ds_snps = std::vector < stats2D > (N);
			rsquared_spl_gt_snps = std::vector < stats2D > (N);

			//indels
			// Allocating data structures
			genotype_spl_errors_indels = std::vector < unsigned long int > (3*n_diploid + 2*n_haploid, 0);
			genotype_spl_totals_indels = std::vector < unsigned long int > (3*n_diploid + 2*n_haploid, 0);
			genotype_val_spl_totals_indels = std::vector < unsigned long int > (3*n_diploid + 2*n_haploid, 0);
			filtered_gp_spl_indels = std::vector < unsigned long int > (n_diploid + n_haploid, 0);
			genotype_cal_errors_indels = std::vector < unsigned long int > (ploidyP1 * N_BIN_CAL, 0);
			genotype_cal_totals_indels = std::vector < unsigned long int > (ploidyP1 * N_BIN_CAL, 0);
			rsquared_spl_ds_indels = std::vector < stats2D > (N);
			rsquared_spl_gt_indels = std::vector < stats2D > (N);

			if (L > 0) {
				genotype_bin_errors_all = std::vector < unsigned long int > (ploidyP1 * L, 0);
				genotype_bin_totals_all = std::vector < unsigned long int > (ploidyP1 * L, 0);
				genotype_val_bin_totals_all = std::vector < unsigned long int > (ploidyP1 * L, 0);
				filtered_gp_bin_all = std::vector < unsigned long int > (L, 0);
				rsquared_bin_ds_all = std::vector < stats2D > (L);
				rsquared_bin_gt_all = std::vector < stats2D > (L);
				frequency_bin_all = std::vector < stats1D > (L);

				genotype_bin_errors_snps = std::vector < unsigned long int > (ploidyP1 * L, 0);
				genotype_bin_totals_snps = std::vector < unsigned long int > (ploidyP1 * L, 0);
				genotype_val_bin_totals_snps = std::vector < unsigned long int > (ploidyP1 * L, 0);
				filtered_gp_bin_snps = std::vector < unsigned long int > (L, 0);
				rsquared_bin_ds_snps = std::vector < stats2D > (L);
				rsquared_bin_gt_snps = std::vector < stats2D > (L);
				frequency_bin_snps = std::vector < stats1D > (L);

				genotype_bin_errors_indels = std::vector < unsigned long int > (ploidyP1 * L, 0);
				genotype_bin_totals_indels = std::vector < unsigned long int > (ploidyP1 * L, 0);
				genotype_val_bin_totals_indels = std::vector < unsigned long int > (ploidyP1 * L, 0);
				filtered_gp_bin_indels = std::vector < unsigned long int > (L, 0);
				rsquared_bin_ds_indels = std::vector < stats2D > (L);
				rsquared_bin_gt_indels = std::vector < stats2D > (L);
				frequency_bin_indels = std::vector < stats1D > (L);

			} else { //groups
				genotype_bin_errors_all = std::vector < unsigned long int > (ploidyP1 * rsquared_str.size(), 0);
				genotype_bin_totals_all = std::vector < unsigned long int > (ploidyP1 * rsquared_str.size(), 0);
				genotype_val_bin_totals_all = std::vector < unsigned long int > (ploidyP1 * rsquared_str.size(), 0);
				filtered_gp_bin_all =  std::vector < unsigned long int > (ploidyP1 * rsquared_str.size(), 0);
				rsquared_bin_ds_all = std::vector < stats2D > (rsquared_str.size());
				rsquared_bin_gt_all = std::vector < stats2D > (rsquared_str.size());
				frequency_bin_all = std::vector < stats1D > (rsquared_str.size());

				genotype_bin_errors_snps = std::vector < unsigned long int > (ploidyP1 * rsquared_str.size(), 0);
				genotype_bin_totals_snps = std::vector < unsigned long int > (ploidyP1 * rsquared_str.size(), 0);
				genotype_val_bin_totals_snps = std::vector < unsigned long int > (ploidyP1 * rsquared_str.size(), 0);
				filtered_gp_bin_snps =  std::vector < unsigned long int > (ploidyP1 * rsquared_str.size(), 0);
				rsquared_bin_ds_snps = std::vector < stats2D > (rsquared_str.size());
				rsquared_bin_gt_snps = std::vector < stats2D > (rsquared_str.size());
				frequency_bin_snps = std::vector < stats1D > (rsquared_str.size());

				genotype_bin_errors_indels = std::vector < unsigned long int > (ploidyP1 * rsquared_str.size(), 0);
				genotype_bin_totals_indels = std::vector < unsigned long int > (ploidyP1 * rsquared_str.size(), 0);
				genotype_val_bin_totals_indels = std::vector < unsigned long int > (ploidyP1 * rsquared_str.size(), 0);
				filtered_gp_bin_indels =  std::vector < unsigned long int > (ploidyP1 * rsquared_str.size(), 0);
				rsquared_bin_ds_indels = std::vector < stats2D > (rsquared_str.size());
				rsquared_bin_gt_indels = std::vector < stats2D > (rsquared_str.size());
				frequency_bin_indels = std::vector < stats1D > (rsquared_str.size());
			}

			assert(n_diploid+n_haploid==N);
			std::string pl1  = n_haploid!=1? "s" : "";
			std::string pl2  = n_diploid!=1? "s" : "";
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
		unsigned long count_es_err_2 = 0, count_es_err_3 = 0, count_es_err_4 = 0, count_es_err_1 = 0;
		unsigned long afnobin = 0;

		int npl_t, ndp_t, npl_arr_t = 0, ndp_arr_t = 0, *pl_arr_t = NULL, *pl_arr_e = NULL, *dp_arr_t = NULL;
		float *af_ptr=NULL,*ds_arr_e = NULL, *gp_arr_e = NULL, float_swap;
		int * gt_arr_e = NULL, ngt_arr_e = 0, ngt_e; //ngt_t is npl_t in case of --gt-validation
		int nds_e, nds_arr_e = 0;
		int naf_f, naf_arr_f=0;
		int ngp_e, ngp_arr_e = 0;
		int *itmp=NULL, mitmp=0, tret=0;

		std::vector < float > GLs = std::vector < float > ((3*n_diploid + 2*n_haploid), 0.0f);
		std::vector < float > DSs = std::vector < float > (N, 0.0f);
		std::vector < int > GTs = std::vector < int > (N, 0);
		std::vector < float > GPs = std::vector < float > ((3*n_diploid + 2*n_haploid), 0.0f);
		std::vector < int > DPs = std::vector < int > (N, 0);

		std::vector < float > sample_dosages(N);

		std::map < std::string, std::pair < int, bool > > :: iterator itG;
		bcf1_t * line_t, * line_e, * line_f;
		int nobi=0,noval=0;
		long int n_disc_var_truth = 0, n_disc_var_est = 0, n_disc_var_af = 0;
		bool nan_printed = false;

		htsFile * out_fp_rej_sites = nullptr;
		bcf_hdr_t * out_hdr_rej_sites = nullptr;
		htsFile * out_fp_conc_sites = nullptr;
		bcf_hdr_t * out_hdr_conc_sites = nullptr;
		htsFile * out_fp_disc_sites = nullptr;
		bcf_hdr_t * out_hdr_disc_sites = nullptr;

		if (out_r2_per_site)
		{
			std::string out_filename_full = out_filename + "_r2_sites.txt.gz";
			out_file_r2_sites.open(out_filename_full);
			out_file_r2_sites << "chr\tpos\trsid\told_rsid\tallele1\tallele2\tmaf\tinfo\tds_r2\n";
		}
		if (out_rej_sites)
		{
			std::string out_filename_full = out_filename + "_rej_sites.bcf";
			std::string out_file_format = "wb";
			out_fp_rej_sites = hts_open(out_filename_full.c_str(),out_file_format.c_str());
			if (nthreads > 1) hts_set_threads(out_fp_rej_sites, nthreads);
			out_hdr_rej_sites = bcf_hdr_dup(sr->readers[2].header);
			if (bcf_hdr_write(out_fp_rej_sites, out_hdr_rej_sites)) vrb.error("Failed to write rejected sites header to output file");
		}
		if (out_conc_sites)
		{
			std::string out_filename_full = out_filename + "_conc_sites.bcf";
			std::string out_file_format = "wb";
			out_fp_conc_sites = hts_open(out_filename_full.c_str(),out_file_format.c_str());
			if (nthreads > 1) hts_set_threads(out_fp_conc_sites, nthreads);
			out_hdr_conc_sites = bcf_hdr_dup(sr->readers[2].header);
			if (bcf_hdr_write(out_fp_conc_sites, out_hdr_conc_sites)) vrb.error("Failed to write concondant sites header to output file");
		}
		if (out_disc_sites)
		{
			std::string out_filename_full = out_filename + "_disc_sites.bcf";
			std::string out_file_format = "wb";
			out_fp_disc_sites = hts_open(out_filename_full.c_str(),out_file_format.c_str());
			if (nthreads > 1) hts_set_threads(out_fp_disc_sites, nthreads);
			out_hdr_disc_sites = bcf_hdr_dup(sr->readers[2].header);
			if (bcf_hdr_write(out_fp_disc_sites, out_hdr_disc_sites)) vrb.error("Failed to write discordant sites header to output file");
		}


		while ((nset = bcf_sr_next_line (sr)))
		{
			bool rej_site = false;
			long int conc_gts = 0;
			long int disc_gts = 0;

			if (nset == 3)
			{
				line_f =  bcf_sr_get_line(sr, 2);

				const int line_type = bcf_get_variant_types(line_f);
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
					if (gt_target && use_gp_filter)
					{
						ngp_e = bcf_get_format_int32(hdr_estimated, line_e, "PL", &pl_arr_e, &ngp_arr_e);//POTENTIAL BUG HERE int->float
						if (ngp_e <= 0 && ngp_e != -3)
							vrb.error("Cannot read PL field from the target required for --min-tar-gp. Check your files or remove --min-tar-gp option to just use the GT field.");
					}
					else ngp_e = bcf_get_format_float(hdr_estimated, line_e, "GP", &gp_arr_e, &ngp_arr_e);

					long int n_validation_marker=0;
					if ((naf_f==1)&&(npl_t%n_true_samples==0)&&(gt_target||((nds_e==n_esti_samples)&&(ngp_e%n_esti_samples==0)))&&((ndp_t==n_true_samples) || D==0))
					{
						// Meta data for variant
						const int ploidy_gt_e_record=ngt_e/n_esti_samples;
						const int ploidy_e_record=ngp_e/n_esti_samples;
						const int ploidy_t_record=npl_t/n_true_samples;

						const bool flip = use_alt_af? false : (af > 0.5f);
						const float maf = use_alt_af? af : std::min(af, 1.0f - af);
						int grp_bin = -1;
						if (L>0) grp_bin = getFrequencyBin(maf);
						else
						{
							std::string chr = std::string(bcf_hdr_id2name(hdr_truth, line_t->rid));
							int pos = line_t->pos + 1;
							std::string ref = line_t->d.allele[0];
							std::string alt = line_t->d.allele[1];

							std::string uuid = n_fields_in_group_files == 5 ? chr + "_" + stb.str(pos) + "_" + ref + "_" + alt : chr + "_" + stb.str(pos);
							itG = site2grp.find(uuid);
							if (itG != site2grp.end()) {
								grp_bin = itG->second.first;
								itG->second.second = true;
							}
							/*
							else if (n_fields_in_group_files==5)
							{
								//second attempt, switch ref alt
								std::string uuid = n_fields_in_group_files == 5 ? chr + "_" + stb.str(pos) + "_" + alt + "_" + ref : chr + "_" + stb.str(pos);
								itG = site2grp.find(uuid);
								if (itG != site2grp.end()) {
									grp_bin = itG->second.first;
									itG->second.second = true;
								}
							}
							*/
						}

						// Process variant
						if (grp_bin >= 0)
						{	// Do this variant fall within a given frequency bin?

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
									DPs[i] = D > 0 ? dp_arr_t[i] : 0;
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
									DPs[i] = D > 0 ? dp_arr_t[i] : 0;
								}
								if (flip) { float_swap = GLs[gpos+ploidy]; GLs[gpos+ploidy] = GLs[gpos+0]; GLs[gpos+0] = float_swap; }
							}

							// Read Estimates
							for(int i = 0 ; i < n_esti_samples ; i ++)
							{
								const int gpos=ind2gpos[i];
								const int pgt = (ploidy_samples[i]>1);
								const int index = ploidy_gt_e_record*i;//ploidy

								if ( (gt_arr_e[index+0]==bcf_gt_missing || gt_arr_e[index+pgt]==bcf_gt_missing)
									|| (gt_arr_e[index+0]==bcf_int32_vector_end || gt_arr_e[index+pgt]==bcf_int32_vector_end)
									|| std::isnan(gt_arr_e[index+0]) || std::isnan(gt_arr_e[index+pgt]))
								{
									GTs[i] = -1;
									GPs[gpos+0] = -1;
									GPs[gpos+1] = -1;
									GPs[gpos+ploidy] = -1;
									DSs[i] = -1;
								}
								else
								{
									int gt = (bcf_gt_allele(gt_arr_e[index+0])==1);
									if (ploidy_samples[i]>1) gt += (bcf_gt_allele(gt_arr_e[index+1])==1);
									GTs[i] = gt;

									if (gt_target)
									{
										if (use_gp_filter)
										{
											if ( (pl_arr_e[ploidy_e_record*i+0]==bcf_int32_missing || pl_arr_e[ploidy_e_record*i+1]==bcf_int32_missing || pl_arr_e[ploidy_e_record*i+ploidy]==bcf_int32_missing)
													|| (pl_arr_e[ploidy_e_record*i+0]==bcf_int32_vector_end || pl_arr_e[ploidy_e_record*i+1]==bcf_int32_vector_end || pl_arr_e[ploidy_e_record*i+ploidy]==bcf_int32_vector_end)
													|| (pl_arr_e[ploidy_e_record*i+0]<0) || (pl_arr_e[ploidy_e_record*i+1]<0) || (pl_arr_e[ploidy_e_record*i+ploidy]<0))
											{
												GPs[gpos+0]=GPs[gpos+1]=GPs[gpos+ploidy]=DSs[i]=-1;
											}
											else
											{
												const float gl0 = unphred[std::min(pl_arr_e[ploidy_e_record*i+0],255)];
												const float gl1 = unphred[std::min(pl_arr_e[ploidy_e_record*i+1],255)];
												const float gl2 = ploidy > 1 ? unphred[std::min(pl_arr_e[ploidy_e_record*i+2],255)] : 0;
												const float sum = gl0 + gl1 + gl2;
												if (sum<=0) vrb.error("Error reading PL field in target");
												GPs[gpos+0] = gl0 / sum;
												GPs[gpos+1] = gl1 / sum;
												if (ploidy > 1) GPs[gpos+2] = gl2 / sum;
												DSs[i] = ploidy > 1 ? 2*GPs[gpos+2]+GPs[gpos+1] : GPs[gpos+1];
											}
										}
										else
										{
											DSs[i] = GTs[i];
											GPs[gpos+0] = (GTs[i] == 0);
											GPs[gpos+1] = (GTs[i] == 1);
											GPs[gpos+ploidy] = (GTs[i] == ploidy);
										}
									}
									else
									{

										if (std::isnan( gp_arr_e[ploidy_e_record*i+0]) || std::isnan( gp_arr_e[ploidy_e_record*i+1]) || std::isnan( gp_arr_e[ploidy_e_record*i+2]))
										{
											if (!nan_printed)
											{
												vrb.warning("Found NAN values in the imputed data. Skipping site. Please check you files. This message is printed only once.");
												nan_printed = true;
											}
											GPs[gpos+0]=GPs[gpos+1]=GPs[gpos+ploidy]=DSs[i]=-1;
										}
										else
										{
											DSs[i] = ds_arr_e[i];

											GPs[gpos+0] = gp_arr_e[ploidy_e_record*i+0];
											GPs[gpos+1] = gp_arr_e[ploidy_e_record*i+1];
											GPs[gpos+ploidy] = gp_arr_e[ploidy_e_record*i+ploidy];
										}
									}
									if (flip) { float_swap = GPs[gpos+ploidy]; GPs[gpos+ploidy] = GPs[gpos+0]; GPs[gpos+0] = float_swap; DSs[i] = ((float) ploidy) - DSs[i]; GTs[i] = ploidy - GTs[i];}
								}
							}

							for (int i = 0 ; i < N ; i ++)
							{
								const int gpos=ind2gpos[i];
								const int true_genotype = getTruth(GLs[gpos+0], GLs[gpos+1], GLs[gpos+ploidy], DPs[i], ploidy_samples[i]);
								if (true_genotype >= 0)
								{
									genotype_val_spl_totals_all[gpos+true_genotype]++;
									genotype_val_bin_totals_all[ploidyP1*grp_bin+true_genotype]++;
									if (line_type==VCF_SNP)
									{
										genotype_val_spl_totals_snps[gpos+true_genotype]++;
										genotype_val_bin_totals_snps[ploidyP1*grp_bin+true_genotype]++;
									}
									if (line_type==VCF_INDEL)
									{
										genotype_val_spl_totals_indels[gpos+true_genotype]++;
										genotype_val_bin_totals_indels[ploidyP1*grp_bin+true_genotype]++;
									}

									int esti_genotype;
									//if (gt_target) esti_genotype = getMostLikelyGT(GTs[i],DSs[i], ploidy_samples[i]);
									//else
									esti_genotype = getMostLikely(GPs[gpos+0], GPs[gpos+1], GPs[gpos+ploidy], DSs[i], ploidy_samples[i], gp_filter);

									if (esti_genotype >= 0)
									{
										const int cal_bin = getCalibrationBin(GPs[gpos+0], GPs[gpos+1], GPs[gpos+ploidy]);
										const int is_error = (true_genotype != esti_genotype);
										conc_gts+=!is_error;
										disc_gts+=is_error;
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
										r2_variant.push(DSs[i], true_genotype*1.0f);
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
										rej_site = true;
										switch (esti_genotype)
										{
											case -2:  count_es_err_2++; break;
											case -3:  count_es_err_3++; break;
											case -4:
												filtered_gp_spl_all[i] ++;
												filtered_gp_bin_all[grp_bin] ++;

												if (line_type==VCF_SNP)
												{
													filtered_gp_spl_snps[i] ++;
													filtered_gp_bin_snps[grp_bin] ++;

												}
												else if (line_type==VCF_INDEL)
												{
													filtered_gp_spl_indels[i] ++;
													filtered_gp_bin_indels[grp_bin] ++;
												}

												count_es_err_4++;
												break;
											default:  count_es_err_1++; break;
										}
									}
								}
								else
								{
									rej_site = true;
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
							if (out_r2_per_site)
							{
								float min_all_freq = -1.0;
					            int AC = -1, AN = 0;
					            tret = bcf_get_info_int32(sr->readers[2].header, line_f, "AN", &itmp, &mitmp);
					            if ( tret==1 )
					            {
					                AN = itmp[0];
					                tret = bcf_get_info_int32(sr->readers[2].header, line_f, "AC", &itmp, &mitmp);
					                if ( tret>0 )
					                    AC = itmp[0];

					                if ( AN>0 && AC>=0 ) min_all_freq = (float) AC / AN;
					            }
					            std::string smin_all_freq = min_all_freq >= 0 ? stb.str(std::min(min_all_freq, 1.0f - min_all_freq),6) : "NA";
								out_file_r2_sites <<  std::string(bcf_hdr_id2name(sr->readers[2].header, line_f->rid)) + "\t" + stb.str(line_f->pos + 1) + "\t" + std::string(line_f->d.id) + "\t" + std::string(line_e->d.id) + "\t" + std::string(line_f->d.allele[0]) + "\t" + std::string(line_f->d.allele[1]) + "\t" + smin_all_freq + "\t" + stb.str(maf,6) + "\t" +  stb.str(std::pow(r2_variant.corrXY(),2),6) + "\n";
							}
							r2_variant.clear();
						} else { rej_site = true; afnobin++;}
					}
					else { rej_site = true; noval++;}

					// increment number of variants
					nvariantvalall += n_validation_marker>0;
					nvariantvalfull += n_validation_marker==N;

					if (line_type==VCF_SNP && n_validation_marker) nvariantvalsnps++;
					if (line_type==VCF_INDEL && n_validation_marker) nvariantvalindels++;

					nvarianttot ++;
				}
				else { rej_site = true; nobi++;}

				if (out_rej_sites && rej_site) write_record(out_fp_rej_sites, out_hdr_rej_sites, sr->readers[2].header, line_f);

				if (!rej_site)
				{
					if (out_conc_sites && disc_gts == 0 && conc_gts > 0) write_record(out_fp_conc_sites, out_hdr_conc_sites, sr->readers[2].header, line_f);
					if (out_disc_sites && disc_gts > 0) write_record(out_fp_disc_sites, out_hdr_disc_sites, sr->readers[2].header, line_f);
				}
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
		if (itmp) free(itmp);
		free(af_ptr);
		free(pl_arr_t);
		free(pl_arr_e);
		free(gt_arr_e);
		free(ds_arr_e);
		free(gp_arr_e);
		free(dp_arr_t);
		bcf_hdr_destroy(hdr_truth);
		bcf_hdr_destroy(hdr_estimated);
		bcf_sr_destroy(sr);

		if (out_r2_per_site) out_file_r2_sites.close();
		if (out_rej_sites)
		{
			bcf_hdr_destroy(out_hdr_rej_sites);
			if (hts_close(out_fp_rej_sites)) vrb.error("Non zero status when closing VCF/BCF file descriptor [rejected sites]");
		}
		if (out_conc_sites)
		{
			bcf_hdr_destroy(out_hdr_conc_sites);
			if (hts_close(out_fp_conc_sites)) vrb.error("Non zero status when closing VCF/BCF file descriptor [concordant sites]");
		}
		if (out_disc_sites)
		{
			bcf_hdr_destroy(out_hdr_disc_sites);
			if (hts_close(out_fp_disc_sites)) vrb.error("Non zero status when closing VCF/BCF file descriptor [discordant sites]");
		}
		n_variants_all_chromosomes += nvariantvalall;
		if (afnobin > 0) vrb.bullet("#variants discarded as no associated bin can be set (e.g. AF<=0 or AF>1 or no group) = " + stb.str(afnobin));
		if (nvarianttot==0) vrb.error("No variant found in the intersection of files. Files are probably not aligned correctly. Please verify that chromosome names and regions are matching for the imputed, validation and allele frequency file.");
		if (nvariantvalall==0) vrb.error("No usable validation variant has been found in the intersection of files. Verify that the validation file has FORMAT/PL and FORMAT/DP fields defined, the imputed file has FORMAT/DS and FORMAT/GP fields defined in the same region and the allele frequency file has the INFO/AF (or --af-tag) field defined at every marker.");

		vrb.print("");
		vrb.bullet("#variants in the overlap (biallelic SNPs and indels) = " + stb.str(nvarianttot));
		vrb.bullet("#variants with all genotypes in the validation data = " + stb.str(nvariantvalfull));
		vrb.bullet("#variants with at least one genotype in the validation data = " + stb.str(nvariantvalall) + " [SNPs = " + stb.str(nvariantvalsnps) + " (" + stb.str(nvariantvalsnps*100.0f/nvariantvalall) + "%), indels = " + stb.str(nvariantvalindels)  + " (" + stb.str(nvariantvalindels*100.0f/nvariantvalall) + "%)]");

		std::string ra,rs,ri;
		rs = ngenoval_all?stb.str(ngenoval_snps*100.0f/ngenoval_all):"-";
		ri = ngenoval_all?stb.str(ngenoval_indels*100.0f/ngenoval_all):"-";
		vrb.bullet("#genotypes used in validation = " + stb.str(ngenoval_all) + " [#SNPs = " + stb.str(ngenoval_snps) + " (" + rs + "%), indels = " + stb.str(ngenoval_indels)  + " (" + ri + "%)]");

		vrb.print("");
		vrb.bullet("Statistics on discarded true genoypes:");
		vrb.print("     #FORMAT/DP missing: " + stb.str(count_tg_err_5));
		vrb.print("     #FORMAT/DP < min-val-dp: " + stb.str(count_tg_err_2));
		if (gt_validation)
		{
			vrb.print("     #Missing/malformatted GTs: " + stb.str(count_tg_err_3));
		}
		else
		{
			vrb.print("     #Missing/malformatted PLs: " + stb.str(count_tg_err_3));
			vrb.print("     #Max(GLs) < min-val-gl: " + stb.str(count_tg_err_4));
			vrb.print("     #Other (e.g. PLs have the same value): " + stb.str(count_tg_err_1));
		}

		vrb.print("");
		vrb.bullet("Statistics on discarded dosages / genotype probabilities:");
		vrb.print("     #Missing/malformatted DSs: " + stb.str(count_es_err_2));
		if (gt_target) vrb.print("     #Missing/malformatted GTs: " + stb.str(count_es_err_3));
		else vrb.print("     #Missing/malformatted GPs: " + stb.str(count_es_err_3));
		vrb.print("     #Filtered target GPs [gp_filter: " + stb.str(gp_filter) + "]: " + stb.str(count_es_err_4));
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
	if (isec_samples_names) free(isec_samples_names);

	vrb.print("");
	vrb.bullet("Total #variants = " + stb.str(n_variants_all_chromosomes));

	if (L == 0) {
		unsigned int n_found_groups = 0;
		for (std::map < std::string, std::pair < int, bool > > :: iterator itG = site2grp.begin() ; itG != site2grp.end() ; ++ itG) n_found_groups += itG->second.second;
		vrb.bullet("Total #variants in groups found = " + stb.str(n_found_groups));
	}
}
