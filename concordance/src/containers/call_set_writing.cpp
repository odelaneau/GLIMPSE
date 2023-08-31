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

void call_set::writeData(std::string fout) {
	vrb.title("Writting output files");

	const int ploidyP1 = ploidy + 1;

	// [1] Concordance per sample
	vrb.bullet("Concordance per sample: [" + fout + ".error.spl.txt.gz]");
	output_file fd1 (fout + ".error.spl.txt.gz");

	fd1<<"#Genotype concordance by sample (SNPs)\n";
	fd1 << "#GCsS" << " ";
	if (ploidy > 1) fd1 << "id sample_name #val_gt_RR #val_gt_RA #val_gt_AA #filtered_gp RR_hom_matches RA_het_matches AA_hom_matches RR_hom_mismatches RA_het_mismatches AA_hom_mismatches RR_hom_mismatches_rate_percent RA_het_mismatches_rate_percent AA_hom_mimatches non_reference_discordanc_rate_percent best_gt_rsquared imputed_ds_rsquared\n";
	else fd1 << "id sample_name #val_gt_R #val_gt_A #filtered_gp R_matches A_matches R_mismatches A_mismatches R_mismatches_rate_percent A_mismatches_rate_percent non_reference_discordanc_rate_percent best_gt_rsquared imputed_ds_rsquared\n";

	for (int i = 0 ; i < N ; i++)
	{
		int gpos=ind2gpos[i];

		unsigned long nrd_mm = genotype_spl_errors_snps[gpos+0];
		unsigned long nrd_m = 0;
		for (int p=1; p<= ploidy; ++p)
		{
			nrd_mm += genotype_spl_errors_snps[gpos+p];
			nrd_m += (genotype_spl_totals_snps[gpos+p] - genotype_spl_errors_snps[gpos+p]);
		}
		double nrd = (nrd_m+nrd_mm) > 0 ? nrd_mm*100.0/(nrd_m+nrd_mm) : 0.0;

		nrd = (nrd_m+nrd_mm) > 0 ? nrd_mm*100.0/(nrd_m+nrd_mm) : 0.0;

		fd1 << "GCsS ";
		fd1 << stb.str(i) << " ";
		fd1 << samples[i] << " ";

		fd1 << (genotype_val_spl_totals_snps[gpos+0]) << " ";
		fd1 << (genotype_val_spl_totals_snps[gpos+1]) << " ";
		if (ploidy > 1) fd1 << (genotype_val_spl_totals_snps[gpos+2]) << " ";

		fd1 << (filtered_gp_spl_snps[i]) << " ";

		fd1 << (genotype_spl_totals_snps[gpos+0]-genotype_spl_errors_snps[gpos+0]) << " ";
		fd1 << (genotype_spl_totals_snps[gpos+1]-genotype_spl_errors_snps[gpos+1]) << " ";
		if (ploidy > 1) fd1 << (genotype_spl_totals_snps[gpos+2]-genotype_spl_errors_snps[gpos+2]) << " ";

		fd1  << genotype_spl_errors_snps[gpos+0] << " ";
		fd1 << genotype_spl_errors_snps[gpos+1] << " ";
		if (ploidy > 1) fd1 << genotype_spl_errors_snps[gpos+2] << " ";

		double d0 = genotype_spl_totals_snps[gpos+0]>0? genotype_spl_errors_snps[gpos+0] * 100.0 / genotype_spl_totals_snps[gpos+0] : 0;
		double d1 = genotype_spl_totals_snps[gpos+1]>0? genotype_spl_errors_snps[gpos+1] * 100.0 / genotype_spl_totals_snps[gpos+1] : 0;
		double d2 = (ploidy > 1 && genotype_spl_totals_snps[gpos+2])>0? genotype_spl_errors_snps[gpos+2] * 100.0 / genotype_spl_totals_snps[gpos+2] : 0;

		fd1 << std::setprecision(3) << std::fixed << d0 << " ";
		fd1 << std::setprecision(3) << std::fixed << d1 << " ";
		if (ploidy > 1) fd1 << std::setprecision(3) << std::fixed << d2 << " ";
		//nrd rate
		fd1 << std::setprecision(3) << std::fixed << std::round(nrd * 1000000.0)/1000000.0 << " ";

		double rsq0 = rsquared_spl_gt_snps[i].corrXY();
		fd1 << std::setprecision(6) << std::fixed << std::round(rsq0 * rsq0 * 1000000.0)/1000000.0 << " ";

		rsq0 = rsquared_spl_ds_snps[i].corrXY();
		fd1 << std::setprecision(6) << std::fixed << std::round(rsq0 * rsq0 * 1000000.0)/1000000.0;

		fd1 << std::endl;
	}

	fd1<<"#Genotype concordance by sample (indels)\n";
	fd1 << "#GCsI" << " ";
	if (ploidy > 1) fd1 << "id sample_name #val_gt_RR #val_gt_RA #val_gt_AA #filtered_gp RR_hom_matches RA_het_matches AA_hom_matches RR_hom_mismatches RA_het_mismatches AA_hom_mismatches RR_hom_mismatches_rate_percent RA_het_mismatches_rate_percent AA_hom_mimatches non_reference_discordanc_rate_percent best_gt_rsquared imputed_ds_rsquared\n";
	else fd1 << "id sample_name #val_gt_R #val_gt_A #filtered_gp R_matches A_matches R_mismatches A_mismatches R_mismatches_rate_percent A_mismatches_rate_percent non_reference_discordanc_rate_percent best_gt_rsquared imputed_ds_rsquared\n";

	for (int i = 0 ; i < N ; i++)
	{
		int gpos=ind2gpos[i];

		unsigned long nrd_mm = genotype_spl_errors_indels[gpos+0];
		unsigned long nrd_m = 0;
		for (int p=1; p<= ploidy; ++p)
		{
			nrd_mm += genotype_spl_errors_indels[gpos+p];
			nrd_m += (genotype_spl_totals_indels[gpos+p] - genotype_spl_errors_indels[gpos+p]);
		}
		double nrd = (nrd_m+nrd_mm) > 0 ? nrd_mm*100.0/(nrd_m+nrd_mm) : 0.0;

		nrd = (nrd_m+nrd_mm) > 0 ? nrd_mm*100.0/(nrd_m+nrd_mm) : 0.0;

		fd1 << "GCsI ";
		fd1 << stb.str(i) << " ";
		fd1 << samples[i] << " ";

		fd1 << (genotype_val_spl_totals_indels[gpos+0]) << " ";
		fd1 << (genotype_val_spl_totals_indels[gpos+1]) << " ";
		if (ploidy > 1) fd1 << (genotype_val_spl_totals_indels[gpos+2]) << " ";

		fd1 << (filtered_gp_spl_indels[i]) << " ";

		fd1 << (genotype_spl_totals_indels[gpos+0]-genotype_spl_errors_indels[gpos+0]) << " ";
		fd1 << (genotype_spl_totals_indels[gpos+1]-genotype_spl_errors_indels[gpos+1]) << " ";
		if (ploidy > 1) fd1 << (genotype_spl_totals_indels[gpos+2]-genotype_spl_errors_indels[gpos+2]) << " ";

		fd1  << genotype_spl_errors_indels[gpos+0] << " ";
		fd1 << genotype_spl_errors_indels[gpos+1] << " ";
		if (ploidy > 1) fd1 << genotype_spl_errors_indels[gpos+2] << " ";

		double d0 = genotype_spl_totals_indels[gpos+0]>0? genotype_spl_errors_indels[gpos+0] * 100.0 / genotype_spl_totals_indels[gpos+0] : 0;
		double d1 = genotype_spl_totals_indels[gpos+1]>0? genotype_spl_errors_indels[gpos+1] * 100.0 / genotype_spl_totals_indels[gpos+1] : 0;
		double d2 = (ploidy > 1 && genotype_spl_totals_indels[gpos+2])>0? genotype_spl_errors_indels[gpos+2] * 100.0 / genotype_spl_totals_indels[gpos+2] : 0;

		fd1 << std::setprecision(3) << std::fixed << d0 << " ";
		fd1 << std::setprecision(3) << std::fixed << d1 << " ";
		if (ploidy > 1) fd1 << std::setprecision(3) << std::fixed << d2 << " ";

		//nrd rate
		fd1 << std::setprecision(3) << std::fixed << std::round(nrd * 1000000.0)/1000000.0 << " ";

		double rsq0 = rsquared_spl_gt_indels[i].corrXY();
		fd1 << std::setprecision(6) << std::fixed << std::round(rsq0 * rsq0 * 1000000.0)/1000000.0 << " ";

		rsq0 = rsquared_spl_ds_indels[i].corrXY();
		fd1 << std::setprecision(6) << std::fixed << std::round(rsq0 * rsq0 * 1000000.0)/1000000.0;

		fd1 << std::endl;
	}

	fd1<<"#Genotype concordance by sample (Variants: SNPs + indels)\n";
	fd1 << "#GCsV" << " ";
	if (ploidy > 1) fd1 << "id sample_name #val_gt_RR #val_gt_RA #val_gt_AA #filtered_gp RR_hom_matches RA_het_matches AA_hom_matches RR_hom_mismatches RA_het_mismatches AA_hom_mismatches RR_hom_mismatches_rate_percent RA_het_mismatches_rate_percent AA_hom_mimatches non_reference_discordanc_rate_percent best_gt_rsquared imputed_ds_rsquared\n";
	else fd1 << "id sample_name #val_gt_R #val_gt_A #filtered_gp R_matches A_matches R_mismatches A_mismatches R_mismatches_rate_percent A_mismatches_rate_percent non_reference_discordanc_rate_percent best_gt_rsquared imputed_ds_rsquared\n";

	for (int i = 0 ; i < N ; i++)
	{
		int gpos=ind2gpos[i];

		unsigned long nrd_mm = genotype_spl_errors_all[gpos+0];
		unsigned long nrd_m = 0;
		for (int p=1; p<= ploidy; ++p)
		{
			nrd_mm += genotype_spl_errors_all[gpos+p];
			nrd_m += (genotype_spl_totals_all[gpos+p] - genotype_spl_errors_all[gpos+p]);
		}
		double nrd = (nrd_m+nrd_mm) > 0 ? nrd_mm*100.0/(nrd_m+nrd_mm) : 0.0;

		nrd = (nrd_m+nrd_mm) > 0 ? nrd_mm*100.0/(nrd_m+nrd_mm) : 0.0;

		fd1 << "GCsV" << " ";
		fd1 << stb.str(i) << " ";
		fd1 << samples[i] << " ";

		fd1 << (genotype_val_spl_totals_all[gpos+0]) << " ";
		fd1 << (genotype_val_spl_totals_all[gpos+1]) << " ";
		if (ploidy > 1) fd1 << (genotype_val_spl_totals_all[gpos+2]) << " ";

		fd1 << (filtered_gp_spl_all[i]) << " ";

		fd1 << (genotype_spl_totals_all[gpos+0]-genotype_spl_errors_all[gpos+0]) << " ";
		fd1 << (genotype_spl_totals_all[gpos+1]-genotype_spl_errors_all[gpos+1]) << " ";
		if (ploidy > 1) fd1 << (genotype_spl_totals_all[gpos+2]-genotype_spl_errors_all[gpos+2]) << " ";

		fd1  << genotype_spl_errors_all[gpos+0] << " ";
		fd1 << genotype_spl_errors_all[gpos+1] << " ";
		if (ploidy > 1) fd1 << genotype_spl_errors_all[gpos+2] << " ";

		double d0 = genotype_spl_totals_all[gpos+0]>0? genotype_spl_errors_all[gpos+0] * 100.0 / genotype_spl_totals_all[gpos+0] : 0;
		double d1 = genotype_spl_totals_all[gpos+1]>0? genotype_spl_errors_all[gpos+1] * 100.0 / genotype_spl_totals_all[gpos+1] : 0;
		double d2 = (ploidy > 1 && genotype_spl_totals_all[gpos+2])>0? genotype_spl_errors_all[gpos+2] * 100.0 / genotype_spl_totals_all[gpos+2] : 0;

		fd1 << std::setprecision(3) << std::fixed << d0 << " ";
		fd1 << std::setprecision(3) << std::fixed << d1 << " ";
		if (ploidy > 1) fd1 << std::setprecision(3) << std::fixed << d2 << " ";

		//nrd rate
		fd1 << std::setprecision(3) << std::fixed << std::round(nrd * 1000000.0)/1000000.0 << " ";

		double rsq0 = rsquared_spl_gt_all[i].corrXY();
		fd1 << std::setprecision(6) << std::fixed << std::round(rsq0 * rsq0 * 1000000.0)/1000000.0 << " ";

		rsq0 = rsquared_spl_ds_all[i].corrXY();
		fd1 << std::setprecision(6) << std::fixed << std::round(rsq0 * rsq0 * 1000000.0)/1000000.0;
		fd1 << std::endl;
	}

	fd1.close();

	// [2] Concordance per bin

	fd1.open(fout + ".error.grp.txt.gz");
	if (L > 0)
	{
		vrb.bullet("Concordance by frequency bin: [" + fout + ".error.grp.txt.gz]");
		fd1<<"#Genotype concordance by allele frequency bin (SNPs)\n";
		fd1 << "#GCsSAF" << " ";
		if (ploidy > 1) fd1 << "id n_genotypes mean_AF #val_gt_RR #val_gt_RA #val_gt_AA filtered_gp RR_hom_matches RA_het_matches AA_hom_matches RR_hom_mismatches RA_het_mismatches AA_hom_mismatches RR_hom_mismatches_rate_percent RA_het_mismatches_rate_percent AA_hom_mimatches_rate_percent best_gt_rsquared imputed_ds_rsquared\n";
		else fd1 << "id n_genotypes mean_AF #val_gt_R #val_gt_A filtered_gp R_matches A_matches R_mismatches A_mismatches R_mismatches_rate_percent A_mismatches_rate_percent best_gt_rsquared imputed_ds_rsquared\n";


		for (int b = 0 ; b < L ; b++)
		{
			fd1 << "GCsSAF" << " " << b << " " << frequency_bin_snps[b].size() << " " << frequency_bin_snps[b].mean() << " ";
			fd1 << (genotype_val_bin_totals_snps[ploidyP1*b+0]) << " ";
			fd1 << (genotype_val_bin_totals_snps[ploidyP1*b+1]) << " ";
			if (ploidy > 1) fd1 << (genotype_val_bin_totals_snps[ploidyP1*b+2]) << " ";
			fd1 << filtered_gp_bin_snps[b] << " ";
			fd1 << (genotype_bin_totals_snps[ploidyP1*b+0] - genotype_bin_errors_snps[ploidyP1*b+0]) << " ";
			fd1 << (genotype_bin_totals_snps[ploidyP1*b+1] - genotype_bin_errors_snps[ploidyP1*b+1]) << " ";
			if (ploidy > 1) fd1 << (genotype_bin_totals_snps[ploidyP1*b+2] - genotype_bin_errors_snps[ploidyP1*b+2]) << " ";
			fd1 << genotype_bin_errors_snps[ploidyP1*b+0] << " ";
			fd1 << genotype_bin_errors_snps[ploidyP1*b+1] << " ";
			if (ploidy > 1) fd1 << genotype_bin_errors_snps[ploidyP1*b+2] << " ";

			double d0 = genotype_bin_totals_snps[ploidyP1*b+0]>0? genotype_bin_errors_snps[ploidyP1*b+0] * 100.0 / genotype_bin_totals_snps[ploidyP1*b+0] : 0;
			double d1 = genotype_bin_totals_snps[ploidyP1*b+1]>0? genotype_bin_errors_snps[ploidyP1*b+1] * 100.0 / genotype_bin_totals_snps[ploidyP1*b+1] : 0;
			double d2 = (ploidy > 1 && genotype_bin_totals_snps[ploidyP1*b+2]>0)? genotype_bin_errors_snps[ploidyP1*b+2] * 100.0 / genotype_bin_totals_snps[ploidyP1*b+2] : 0;

			fd1 << std::setprecision(3) << std::fixed << d0 << " ";
			fd1 << std::setprecision(3) << std::fixed << d1 << " ";
			if (ploidy > 1) fd1 << std::setprecision(3) << std::fixed << d2 << " ";

			double rsq0 = rsquared_bin_gt_snps[b].corrXY();
			fd1 << std::setprecision(6) << std::fixed << std::round(rsq0* rsq0 * 1000000.0)/1000000.0 << " ";

			rsq0 = rsquared_bin_ds_snps[b].corrXY();
			fd1 << std::setprecision(6) << std::fixed << std::round(rsq0* rsq0 * 1000000.0)/1000000.0;

			fd1 << std::endl;
		}

		fd1<<"#Genotype concordance by allele frequency bin (indels)\n";
		fd1 << "#GCsIAF" << " ";
		if (ploidy > 1) fd1 << "id n_genotypes mean_AF #val_gt_RR #val_gt_RA #val_gt_AA filtered_gp RR_hom_matches RA_het_matches AA_hom_matches RR_hom_mismatches RA_het_mismatches AA_hom_mismatches RR_hom_mismatches_rate_percent RA_het_mismatches_rate_percent AA_hom_mimatches_rate_percent best_gt_rsquared imputed_ds_rsquared\n";
		else fd1 << "id n_genotypes mean_AF #val_gt_R #val_gt_A filtered_gp R_matches A_matches R_mismatches A_mismatches R_mismatches_rate_percent A_mismatches_rate_percent best_gt_rsquared imputed_ds_rsquared\n";

		for (int b = 0 ; b < L ; b++)
		{
			fd1 << "GCsIAF" << " " << b << " " << frequency_bin_indels[b].size() << " " << frequency_bin_indels[b].mean() << " ";
			fd1 << (genotype_val_bin_totals_indels[ploidyP1*b+0]) << " ";
			fd1 << (genotype_val_bin_totals_indels[ploidyP1*b+1]) << " ";
			if (ploidy > 1) fd1 << (genotype_val_bin_totals_indels[ploidyP1*b+2]) << " ";
			fd1 << filtered_gp_bin_indels[b] << " ";
			fd1 << (genotype_bin_totals_indels[ploidyP1*b+0] - genotype_bin_errors_indels[ploidyP1*b+0]) << " ";
			fd1 << (genotype_bin_totals_indels[ploidyP1*b+1] - genotype_bin_errors_indels[ploidyP1*b+1]) << " ";
			if (ploidy > 1) fd1 << (genotype_bin_totals_indels[ploidyP1*b+2] - genotype_bin_errors_indels[ploidyP1*b+2]) << " ";
			fd1 << genotype_bin_errors_indels[ploidyP1*b+0] << " ";
			fd1 << genotype_bin_errors_indels[ploidyP1*b+1] << " ";
			if (ploidy > 1) fd1 << genotype_bin_errors_indels[ploidyP1*b+2] << " ";

			double d0 = genotype_bin_totals_indels[ploidyP1*b+0]>0? genotype_bin_errors_indels[ploidyP1*b+0] * 100.0 / genotype_bin_totals_indels[ploidyP1*b+0] : 0;
			double d1 = genotype_bin_totals_indels[ploidyP1*b+1]>0? genotype_bin_errors_indels[ploidyP1*b+1] * 100.0 / genotype_bin_totals_indels[ploidyP1*b+1] : 0;
			double d2 = (ploidy > 1 && genotype_bin_totals_indels[ploidyP1*b+2]>0)? genotype_bin_errors_indels[ploidyP1*b+2] * 100.0 / genotype_bin_totals_indels[ploidyP1*b+2] : 0;

			fd1 << std::setprecision(3) << std::fixed << d0 << " ";
			fd1 << std::setprecision(3) << std::fixed << d1 << " ";
			if (ploidy > 1) fd1 << std::setprecision(3) << std::fixed << d2 << " ";

			double rsq0 = rsquared_bin_gt_indels[b].corrXY();
			fd1 << std::setprecision(6) << std::fixed << std::round(rsq0* rsq0 * 1000000.0)/1000000.0 << " ";

			rsq0 = rsquared_bin_ds_indels[b].corrXY();
			fd1 << std::setprecision(6) << std::fixed << std::round(rsq0* rsq0 * 1000000.0)/1000000.0;

			fd1 << std::endl;
		}

		fd1<<"#Genotype concordance by allele frequency bin (Variants: SNPs + indels)\n";
		fd1 << "#GCsVAF" << " ";
		if (ploidy > 1) fd1 << "id n_genotypes mean_AF #val_gt_RR #val_gt_RA #val_gt_AA filtered_gp RR_hom_matches RA_het_matches AA_hom_matches RR_hom_mismatches RA_het_mismatches AA_hom_mismatches RR_hom_mismatches_rate_percent RA_het_mismatches_rate_percent AA_hom_mimatches_rate_percent best_gt_rsquared imputed_ds_rsquared\n";
		else fd1 << "id n_genotypes mean_AF #val_gt_R #val_gt_A filtered_gp R_matches A_matches R_mismatches A_mismatches R_mismatches_rate_percent A_mismatches_rate_percent best_gt_rsquared imputed_ds_rsquared\n";

		for (int b = 0 ; b < L ; b++)
		{
			fd1 << "GCsVAF" << " " << b << " " << frequency_bin_all[b].size() << " " << std::setprecision(6) << std::fixed << frequency_bin_all[b].mean() << " ";
			fd1 << (genotype_val_bin_totals_all[ploidyP1*b+0]) << " ";
			fd1 << (genotype_val_bin_totals_all[ploidyP1*b+1]) << " ";
			if (ploidy > 1) fd1 << (genotype_val_bin_totals_all[ploidyP1*b+2]) << " ";
			fd1 << filtered_gp_bin_all[b] << " ";
			fd1 << (genotype_bin_totals_all[ploidyP1*b+0] - genotype_bin_errors_all[ploidyP1*b+0]) << " ";
			fd1 << (genotype_bin_totals_all[ploidyP1*b+1] - genotype_bin_errors_all[ploidyP1*b+1]) << " ";
			if (ploidy > 1) fd1 << (genotype_bin_totals_all[ploidyP1*b+2] - genotype_bin_errors_all[ploidyP1*b+2]) << " ";
			fd1 << genotype_bin_errors_all[ploidyP1*b+0] << " ";
			fd1 << genotype_bin_errors_all[ploidyP1*b+1] << " ";
			if (ploidy > 1) fd1 << genotype_bin_errors_all[ploidyP1*b+2] << " ";

			double d0 = genotype_bin_totals_all[ploidyP1*b+0]>0? genotype_bin_errors_all[ploidyP1*b+0] * 100.0 / genotype_bin_totals_all[ploidyP1*b+0] : 0;
			double d1 = genotype_bin_totals_all[ploidyP1*b+1]>0? genotype_bin_errors_all[ploidyP1*b+1] * 100.0 / genotype_bin_totals_all[ploidyP1*b+1] : 0;
			double d2 = (ploidy > 1 && genotype_bin_totals_all[ploidyP1*b+2]>0)? genotype_bin_errors_all[ploidyP1*b+2] * 100.0 / genotype_bin_totals_all[ploidyP1*b+2] : 0;

			fd1 << std::setprecision(3) << std::fixed << d0 << " ";
			fd1 << std::setprecision(3) << std::fixed << d1 << " ";
			if (ploidy > 1) fd1 << std::setprecision(3) << std::fixed << d2 << " ";

			double rsq0 = rsquared_bin_gt_all[b].corrXY();
			fd1 << std::setprecision(6) << std::fixed << std::round(rsq0* rsq0 * 1000000.0)/1000000.0 << " ";

			rsq0 = rsquared_bin_ds_all[b].corrXY();
			fd1 << std::setprecision(6) << std::fixed << std::round(rsq0* rsq0 * 1000000.0)/1000000.0;

			fd1 << std::endl;
		}
	}
	else
	{
		vrb.bullet("Concordance by group: [" + fout + ".error.grp.txt.gz]");
		fd1<<"#Genotype concordance by group bin (Variants: SNPs + indels)\n";
		for (int b = 0 ; b < rsquared_str.size() ; b++)
		{
			fd1 << b << " " << rsquared_str[b] << " " << frequency_bin_all[b].size() << " " << frequency_bin_all[b].mean();
			fd1 << " " << genotype_bin_errors_all[ploidyP1*b+0] << " " << genotype_bin_totals_all[ploidyP1*b+0];
			fd1 << " " << genotype_bin_errors_all[ploidyP1*b+1] << " " << genotype_bin_totals_all[ploidyP1*b+1];
			if (ploidy > 1) fd1 << " " << genotype_bin_errors_all[ploidyP1*b+2] << " " << genotype_bin_totals_all[ploidyP1*b+2];
			fd1 << " " << genotype_bin_errors_all[ploidyP1*b+0] * 100.0 / genotype_bin_totals_all[ploidyP1*b+0];
			fd1 << " " << genotype_bin_errors_all[ploidyP1*b+1] * 100.0 / genotype_bin_totals_all[ploidyP1*b+1];
			if (ploidy > 1) fd1 << " " << genotype_bin_errors_all[ploidyP1*b+2] * 100.0 / genotype_bin_totals_all[ploidyP1*b+2];
			fd1 << std::endl;
		}
	}
	fd1.close();

	// [3] Calibration
	vrb.bullet("Concordance per calibration bin: [" + fout + ".error.cal.txt.gz]");
	fd1.open(fout + ".error.cal.txt.gz");
	float step_size = 1.0 / N_BIN_CAL;

	fd1<<"#Genotype caoncordance by calibration bin (SNPs)\n";
	fd1 << "#GCsSC" << " ";
	if (ploidy > 1) fd1 << "id start_bin end_bin avg_bin RR_hom_matches RA_het_matches AA_hom_matches RR_hom_mismatches RA_het_mismatches AA_hom_mismatches RR_hom_mismatches_rate_percent RA_het_mismatches_rate_percent AA_hom_mimatches_rate_percent all_mismatches_rate_percent\n";
	else fd1 << "id start_bin end_bin avg_bin R_matches A_matches R_mismatches A_mismatches R_mismatches_rate_percent A_mismatches_rate_percent all_mismatches_rate_percent\n";

	for (int b = 0 ; b < N_BIN_CAL ; b++) {
		unsigned long int sumE = genotype_cal_errors_snps[ploidyP1*b+0] + genotype_cal_errors_snps[ploidyP1*b+1];
		if (ploidy > 1) sumE += genotype_cal_errors_snps[ploidyP1*b+2];
		unsigned long int sumT = genotype_cal_totals_snps[ploidyP1*b+0] + genotype_cal_totals_snps[ploidyP1*b+1];
		if (ploidy > 1) sumT += genotype_cal_totals_snps[ploidyP1*b+2];

		double errA = sumT? sumE * 100.0 / sumT : 0;

		fd1 << std::setprecision(2) << std::fixed << "GCsSC" << " " << b << " " << b*step_size << " " << (b+1)*step_size << " " << std::setprecision(3) << std::fixed << b*step_size + step_size/2.0 << " ";

		fd1 << (genotype_cal_totals_snps[ploidyP1*b+0] - genotype_cal_errors_snps[ploidyP1*b+0]) << " ";
		fd1 << (genotype_cal_totals_snps[ploidyP1*b+1] - genotype_cal_errors_snps[ploidyP1*b+1]) << " ";
		if (ploidy > 1) fd1 << (genotype_cal_totals_snps[ploidyP1*b+2] - genotype_cal_errors_snps[ploidyP1*b+2]) << " ";

		fd1 << genotype_cal_errors_snps[ploidyP1*b+0] << " ";
		fd1 << genotype_cal_errors_snps[ploidyP1*b+1] << " ";
		if (ploidy > 1) fd1 << std::setprecision(2) << std::fixed << genotype_cal_errors_snps[ploidyP1*b+2] << " ";

		double d0 = genotype_cal_totals_snps[ploidyP1*b+0]>0? genotype_cal_errors_snps[ploidyP1*b+0]*100.0 / genotype_cal_totals_snps[ploidyP1*b+0] : 0;
		double d1 = genotype_cal_totals_snps[ploidyP1*b+1]>0? genotype_cal_errors_snps[ploidyP1*b+1]*100.0 / genotype_cal_totals_snps[ploidyP1*b+1] : 0;
		double d2 = (ploidy > 1 && genotype_cal_totals_snps[ploidyP1*b+2]>0)? genotype_cal_errors_snps[ploidyP1*b+2]*100.0 / genotype_cal_totals_snps[ploidyP1*b+2] : 0;

		fd1 << std::setprecision(6) << std::fixed << d0 << " ";
		fd1 << std::setprecision(6) << std::fixed << d1 << " ";
		if (ploidy > 1) fd1 << std::setprecision(4) << std::fixed << d2 << " ";

		fd1 << std::setprecision(6) << std::fixed << errA;
		fd1 << std::endl;
	}

	fd1<<"#Genotype caoncordance by calibration bin (indels)\n";
	fd1 << "#GCsIC" << " ";
	if (ploidy > 1) fd1 << "id start_bin end_bin avg_bin RR_hom_matches RA_het_matches AA_hom_matches RR_hom_mismatches RA_het_mismatches AA_hom_mismatches RR_hom_mismatches_rate_percent RA_het_mismatches_rate_percent AA_hom_mimatches_rate_percent all_mismatches_rate_percent\n";
	else fd1 << "id start_bin end_bin avg_bin R_matches A_matches R_mismatches A_mismatches R_mismatches_rate_percent A_mismatches_rate_percent all_mismatches_rate_percent\n";

	for (int b = 0 ; b < N_BIN_CAL ; b++) {
		unsigned long int sumE = genotype_cal_errors_indels[ploidyP1*b+0] + genotype_cal_errors_indels[ploidyP1*b+1];
		if (ploidy > 1) sumE += genotype_cal_errors_indels[ploidyP1*b+2];
		unsigned long int sumT = genotype_cal_totals_indels[ploidyP1*b+0] + genotype_cal_totals_indels[ploidyP1*b+1];
		if (ploidy > 1) sumT += genotype_cal_totals_indels[ploidyP1*b+2];

		double errA = sumT? sumE * 100.0 / sumT : 0;

		fd1 << std::setprecision(2) << std::fixed << "GCsIC" << " " << b << " " << b*step_size << " " << (b+1)*step_size << " " << std::setprecision(3) << std::fixed << b*step_size + step_size/2.0 << " ";

		fd1 << (genotype_cal_totals_indels[ploidyP1*b+0] - genotype_cal_errors_indels[ploidyP1*b+0]) << " ";
		fd1 << (genotype_cal_totals_indels[ploidyP1*b+1] - genotype_cal_errors_indels[ploidyP1*b+1]) << " ";
		if (ploidy > 1) fd1 << (genotype_cal_totals_indels[ploidyP1*b+2] - genotype_cal_errors_indels[ploidyP1*b+2]) << " ";

		fd1 << genotype_cal_errors_indels[ploidyP1*b+0] << " ";
		fd1 << genotype_cal_errors_indels[ploidyP1*b+1] << " ";
		if (ploidy > 1) fd1 << std::setprecision(2) << std::fixed << genotype_cal_errors_indels[ploidyP1*b+2] << " ";

		double d0 = genotype_cal_totals_indels[ploidyP1*b+0]>0? genotype_cal_errors_indels[ploidyP1*b+0]*100.0 / genotype_cal_totals_indels[ploidyP1*b+0] : 0;
		double d1 = genotype_cal_totals_indels[ploidyP1*b+1]>0? genotype_cal_errors_indels[ploidyP1*b+1]*100.0 / genotype_cal_totals_indels[ploidyP1*b+1] : 0;
		double d2 = (ploidy > 1 && genotype_cal_totals_indels[ploidyP1*b+2]>0)? genotype_cal_errors_indels[ploidyP1*b+2]*100.0 / genotype_cal_totals_indels[ploidyP1*b+2] : 0;

		fd1 << std::setprecision(6) << std::fixed << d0 << " ";
		fd1 << std::setprecision(6) << std::fixed << d1 << " ";
		if (ploidy > 1) fd1 << std::setprecision(4) << std::fixed << d2 << " ";

		fd1 << std::setprecision(6) << std::fixed << errA;
		fd1 << std::endl;
	}

	fd1<<"#Genotype concordance by calibration bin (Variants: SNPs + indels)\n";
	fd1 << "#GCsVC" << " ";
	if (ploidy > 1) fd1 << "id start_bin end_bin avg_bin RR_hom_matches RA_het_matches AA_hom_matches RR_hom_mismatches RA_het_mismatches AA_hom_mismatches RR_hom_mismatches_rate_percent RA_het_mismatches_rate_percent AA_hom_mimatches_rate_percent all_mismatches_rate_percent\n";
	else fd1 << "id start_bin end_bin avg_bin R_matches A_matches R_mismatches A_mismatches R_mismatches_rate_percent A_mismatches_rate_percent all_mismatches_rate_percent\n";

	for (int b = 0 ; b < N_BIN_CAL ; b++) {
		unsigned long int sumE = genotype_cal_errors_all[ploidyP1*b+0] + genotype_cal_errors_all[ploidyP1*b+1];
		if (ploidy > 1) sumE += genotype_cal_errors_all[ploidyP1*b+2];
		unsigned long int sumT = genotype_cal_totals_all[ploidyP1*b+0] + genotype_cal_totals_all[ploidyP1*b+1];
		if (ploidy > 1) sumT += genotype_cal_totals_all[ploidyP1*b+2];

		double errA = sumT? sumE * 100.0 / sumT : 0;

		fd1 << std::setprecision(2) << std::fixed << "GCsVC" << " " << b << " " << b*step_size << " " << (b+1)*step_size << " " << std::setprecision(3) << std::fixed << b*step_size + step_size/2.0 << " ";

		fd1 << (genotype_cal_totals_all[ploidyP1*b+0] - genotype_cal_errors_all[ploidyP1*b+0]) << " ";
		fd1 << (genotype_cal_totals_all[ploidyP1*b+1] - genotype_cal_errors_all[ploidyP1*b+1]) << " ";
		if (ploidy > 1) fd1 << (genotype_cal_totals_all[ploidyP1*b+2] - genotype_cal_errors_all[ploidyP1*b+2]) << " ";

		fd1 << genotype_cal_errors_all[ploidyP1*b+0] << " ";
		fd1 << genotype_cal_errors_all[ploidyP1*b+1] << " ";
		if (ploidy > 1) fd1 << std::setprecision(2) << std::fixed << genotype_cal_errors_all[ploidyP1*b+2] << " ";

		double d0 = genotype_cal_totals_all[ploidyP1*b+0]>0? genotype_cal_errors_all[ploidyP1*b+0]*100.0 / genotype_cal_totals_all[ploidyP1*b+0] : 0;
		double d1 = genotype_cal_totals_all[ploidyP1*b+1]>0? genotype_cal_errors_all[ploidyP1*b+1]*100.0 / genotype_cal_totals_all[ploidyP1*b+1] : 0;
		double d2 = (ploidy > 1 && genotype_cal_totals_all[ploidyP1*b+2]>0)? genotype_cal_errors_all[ploidyP1*b+2]*100.0 / genotype_cal_totals_all[ploidyP1*b+2] : 0;

		fd1 << std::setprecision(6) << std::fixed << d0 << " ";
		fd1 << std::setprecision(6) << std::fixed << d1 << " ";
		if (ploidy > 1) fd1 << std::setprecision(4) << std::fixed << d2 << " ";

		fd1 << std::setprecision(6) << std::fixed << errA;
		fd1 << std::endl;
	}
	fd1.close();

	// [4] Rsquare per bin
	output_file fd4 (fout + ".rsquare.grp.txt.gz");
	if (L > 0)
	{
		vrb.bullet("Rsquare per frequency bin: [" + fout + ".rsquare.grp.txt.gz" + "]");
		for (int b = 0 ; b < L ; b++) {
			double rsq_gt = rsquared_bin_gt_all[b].corrXY();
			double rsq_ds = rsquared_bin_ds_all[b].corrXY();
			fd4 << b << " " << frequency_bin_all[b].size() << " " << frequency_bin_all[b].mean();
			fd4 << " " << rsq_gt*rsq_gt << " " << rsq_ds*rsq_ds << std::endl;
		}
	} else {
		vrb.bullet("Rsquare per frequency bin: [" + fout + ".rsquare.grp.txt.gz " + "]");
		for (int b = 0 ; b < rsquared_str.size() ; b++)
		{
			double rsq_gt = (rsquared_bin_gt_all[b].sdNAx() || rsquared_bin_gt_all[b].sdNAy()) ? std::numeric_limits<double>::quiet_NaN() : rsquared_bin_gt_all[b].corrXY();
			double rsq_ds = (rsquared_bin_ds_all[b].sdNAx() || rsquared_bin_ds_all[b].sdNAy()) ? std::numeric_limits<double>::quiet_NaN() : rsquared_bin_ds_all[b].corrXY();

			fd4 << rsquared_str[b] << " " << frequency_bin_all[b].size() << " " << frequency_bin_all[b].mean();
			fd4 << " " << rsq_gt*rsq_gt << " " << rsq_ds*rsq_ds << std::endl;
		}
	}
	fd4.close();

	// [5] Rsquare per sample
	output_file fd5 (fout + ".rsquare.spl.txt.gz");
	vrb.bullet("Rsquare per sample: [" + fout + ".rsquare.spl.txt.gz" + "]");

	for (int i = 0 ; i < N ; i++)
	{
		double rsq_gt = rsquared_spl_gt_all[i].corrXY();
		double rsq_ds = rsquared_spl_ds_all[i].corrXY();
		fd5 << samples[i] << " " << rsq_gt*rsq_gt << " " << rsq_ds*rsq_ds << std::endl;
	}
	fd5.close();
}
