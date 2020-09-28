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

void call_set::writeData(string fout) {
	tac.clock();
	vrb.title("Writting output files");

	// [1] Concordance per sample
	vrb.bullet("Concordance per sample");
	output_file fd1 (fout + ".error.spl.txt.gz");
	for (int i = 0 ; i < N ; i++) {
		int gpos=ind2gpos[i];
		fd1 << samples[i];
		fd1 << " " << genotype_spl_errors[gpos+0] << " " << genotype_spl_totals[gpos+0];
		fd1 << " " << genotype_spl_errors[gpos+1] << " " << genotype_spl_totals[gpos+1];
		if (ploidy[i] > 1) fd1 << " " << genotype_spl_errors[gpos+2] << " " << genotype_spl_totals[gpos+2];
		fd1 << " " << genotype_spl_errors[gpos+0] * 100.0 / genotype_spl_totals[gpos+0];
		fd1 << " " << genotype_spl_errors[gpos+1] * 100.0 / genotype_spl_totals[gpos+1];
		if (ploidy[i] > 1) fd1 << " " << genotype_spl_errors[gpos+2] * 100.0 / genotype_spl_totals[gpos+2];
		fd1 << endl;
	}
	fd1.close();

	// [2] Concordance per bin
	output_file fd2 (fout + ".error.grp.txt.gz");
	if (L > 0) {
		vrb.bullet("Concordance per frequency bin");
		for (int b = 0 ; b < L ; b++) {
			fd2 << b << " " << frequency_bin[b].size() << " " << frequency_bin[b].mean();
			fd2 << " " << genotype_bin_errors[(max_ploidy+1)*b+0] << " " << genotype_bin_totals[(max_ploidy+1)*b+0];
			fd2 << " " << genotype_bin_errors[(max_ploidy+1)*b+1] << " " << genotype_bin_totals[(max_ploidy+1)*b+1];
			if (max_ploidy > 1) fd2 << " " << genotype_bin_errors[(max_ploidy+1)*b+2] << " " << genotype_bin_totals[(max_ploidy+1)*b+2];
			fd2 << " " << genotype_bin_errors[(max_ploidy+1)*b+0] * 100.0 / genotype_bin_totals[(max_ploidy+1)*b+0];
			fd2 << " " << genotype_bin_errors[(max_ploidy+1)*b+1] * 100.0 / genotype_bin_totals[(max_ploidy+1)*b+1];
			if (max_ploidy > 1) fd2 << " " << genotype_bin_errors[(max_ploidy+1)*b+2] * 100.0 / genotype_bin_totals[(max_ploidy+1)*b+2];
			fd2 << endl;
		}
	} else {
		vrb.bullet("Concordance per group");
		for (int b = 0 ; b < rsquared_str.size() ; b++) {
			fd2 << b << " " << rsquared_str[b] << " " << frequency_bin[b].size() << " " << frequency_bin[b].mean();
			fd2 << " " << genotype_bin_errors[(max_ploidy+1)*b+0] << " " << genotype_bin_totals[(max_ploidy+1)*b+0];
			fd2 << " " << genotype_bin_errors[(max_ploidy+1)*b+1] << " " << genotype_bin_totals[(max_ploidy+1)*b+1];
			if (max_ploidy > 1) fd2 << " " << genotype_bin_errors[(max_ploidy+1)*b+2] << " " << genotype_bin_totals[(max_ploidy+1)*b+2];
			fd2 << " " << genotype_bin_errors[(max_ploidy+1)*b+0] * 100.0 / genotype_bin_totals[(max_ploidy+1)*b+0];
			fd2 << " " << genotype_bin_errors[(max_ploidy+1)*b+1] * 100.0 / genotype_bin_totals[(max_ploidy+1)*b+1];
			if (max_ploidy > 1) fd2 << " " << genotype_bin_errors[(max_ploidy+1)*b+2] * 100.0 / genotype_bin_totals[(max_ploidy+1)*b+2];
			fd2 << endl;
		}
	}
	fd2.close();

	// [3] Calibration
	vrb.bullet("Concordance per calibration bin");
	output_file fd3 (fout + ".error.cal.txt.gz");
	float step_size = 1.0 / N_BIN_CAL;
	for (int b = 0 ; b < N_BIN_CAL ; b++) {
		float err0 = genotype_cal_errors[(max_ploidy+1)*b+0] * 100.0 / genotype_cal_totals[(max_ploidy+1)*b+0];
		float err1 = genotype_cal_errors[(max_ploidy+1)*b+1] * 100.0 / genotype_cal_totals[(max_ploidy+1)*b+1];
		float err2 = genotype_cal_errors[(max_ploidy+1)*b+2] * 100.0 / genotype_cal_totals[(max_ploidy+1)*b+2];
		float errA = (genotype_cal_errors[(max_ploidy+1)*b+0]+genotype_cal_errors[(max_ploidy+1)*b+1]+genotype_cal_errors[(max_ploidy+1)*b+2]) * 100.0 / (genotype_cal_totals[(max_ploidy+1)*b+0]+genotype_cal_totals[(max_ploidy+1)*b+1]+genotype_cal_totals[(max_ploidy+1)*b+2]);
		fd3 << b << " " << b*step_size << " " << (b+1)*step_size << " " << (b*step_size + step_size/2);
		fd3 << " " << genotype_cal_errors[(max_ploidy+1)*b+0] << " " << genotype_cal_totals[(max_ploidy+1)*b+0];
		fd3 << " " << genotype_cal_errors[(max_ploidy+1)*b+1] << " " << genotype_cal_totals[(max_ploidy+1)*b+1];
		fd3 << " " << genotype_cal_errors[(max_ploidy+1)*b+2] << " " << genotype_cal_totals[(max_ploidy+1)*b+2];
		fd3 << " " << err0 << " " << err1 << " " << err2 << " " << errA << endl;
	}
	fd3.close();

	// [4] Rsquare per bin
	output_file fd4 (fout + ".rsquare.grp.txt.gz");
	if (L > 0) {
		vrb.bullet("Rsquare per frequency bin");
		for (int b = 0 ; b < L ; b++) {
			double rsq0 = rsquared_bin[b].corrXY();
			double rsq2 = rsq0*rsq0;
			fd4 << b << " " << frequency_bin[b].size() << " " << frequency_bin[b].mean();
			fd4 << " " << rsq0 << " " << rsq2 << endl;
		}
	} else {
		vrb.bullet("Rsquare per frequency bin");
		for (int b = 0 ; b < rsquared_str.size() ; b++) {
			double rsq0 = rsquared_bin[b].corrXY();
			double rsq2 = rsq0*rsq0;
			fd4 << b << " " << rsquared_str[b] << " " << frequency_bin[b].size() << " " << frequency_bin[b].mean();
			fd4 << " " << rsq0 << " " << rsq2 << endl;
		}
	}
	fd4.close();

	// [5] Rsquare per sample
	vrb.bullet("Rsquare per sample");
	output_file fd5 (fout + ".rsquare.spl.txt.gz");
	for (int i = 0 ; i < N ; i++) {
		double rsq0 = rsquared_spl[i].corrXY();
		double rsq2 = rsq0*rsq0;
		fd5 << samples[i] << " " << rsq0 << " " << rsq2 << endl;
	}
	fd5.close();

	// [6] AVG Rsquare per bin
	vrb.bullet("Average Rsquare per frequency bin");
	output_file fd6 (fout + ".avg.rsquare.grp.txt.gz");
	for (int b = 0 ; b < L ; b++) {
		fd6 << b << " " << avg_rsquared_bin[b].size() << " " << avg_rsquared_bin[b].variance();
		fd6 << " " << avg_rsquared_bin[b].mean() /*<< " " << rsq2 */<< endl;
	}
	fd6.close();
}
