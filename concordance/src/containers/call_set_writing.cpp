#include "call_set_header.h"

	int N, L, T;
	vector < string > samples;
	vector < double > bins;

	// Per sample concordance [3xN]
	vector < int > genotype_spl_errors;
	vector < int > genotype_spl_totals;

	// Per bin concordance [3xL]
	vector < int > genotype_bin_errors;
	vector < int > genotype_bin_totals;

	// Concordance for calibration [100]
	vector < int > genotype_cal_errors;
	vector < int > genotype_cal_totals;

	// R2 per bin
	vector < stats2D > rsquared_bin;
	vector < stats1D > frequency_bin;

	// R2 per sample
	vector < stats2D > rsquared_spl;

	// R2 per bin x sample
	vector < stats2D > rsquared_binspl;
	vector < stats1D > frequency_binspl;


void call_set::writeData(string fout) {
	tac.clock();
	vrb.title("Writting output files");

	// [1] Concordance per sample
	vrb.bullet("Concordance per sample");
	output_file fd1 (fout + ".error.spl.txt.gz");
	for (int i = 0 ; i < N ; i++) {
		float err0 = genotype_spl_errors[3*i+0] * 100.0 / genotype_spl_totals[3*i+0];
		float err1 = genotype_spl_errors[3*i+1] * 100.0 / genotype_spl_totals[3*i+1];
		float err2 = genotype_spl_errors[3*i+2] * 100.0 / genotype_spl_totals[3*i+2];
		fd1 << samples[i];
		fd1 << " " << genotype_spl_errors[3*i+0] << " " << genotype_spl_totals[3*i+0];
		fd1 << " " << genotype_spl_errors[3*i+1] << " " << genotype_spl_totals[3*i+1];
		fd1 << " " << genotype_spl_errors[3*i+2] << " " << genotype_spl_totals[3*i+2];
		fd1 << " " << err0 << " " << err1 << " " << err2 << endl;
	}
	fd1.close();

	// [2] Concordance per bin
	vrb.bullet("Concordance per frequency bin");
	output_file fd2 (fout + ".error.frq.txt.gz");
	for (int b = 0 ; b < L ; b++) {
		float err0 = genotype_bin_errors[3*b+0] * 100.0 / genotype_bin_totals[3*b+0];
		float err1 = genotype_bin_errors[3*b+1] * 100.0 / genotype_bin_totals[3*b+1];
		float err2 = genotype_bin_errors[3*b+2] * 100.0 / genotype_bin_totals[3*b+2];
		fd2 << b << " " << frequency_bin[b].size() << " " << frequency_bin[b].mean();
		fd2 << " " << genotype_bin_errors[3*b+0] << " " << genotype_bin_totals[3*b+0];
		fd2 << " " << genotype_bin_errors[3*b+1] << " " << genotype_bin_totals[3*b+1];
		fd2 << " " << genotype_bin_errors[3*b+2] << " " << genotype_bin_totals[3*b+2];
		fd2 << " " << err0 << " " << err1 << " " << err2 << endl;
	}
	fd2.close();

	// [3] Calibration
	vrb.bullet("Concordance per calibration bin");
	output_file fd3 (fout + ".error.cal.txt.gz");
	float step_size = 1.0 / N_BIN_CAL;
	for (int b = 0 ; b < N_BIN_CAL ; b++) {
		float err0 = genotype_cal_errors[3*b+0] * 100.0 / genotype_cal_totals[3*b+0];
		float err1 = genotype_cal_errors[3*b+1] * 100.0 / genotype_cal_totals[3*b+1];
		float err2 = genotype_cal_errors[3*b+2] * 100.0 / genotype_cal_totals[3*b+2];
		float errA = (genotype_cal_errors[3*b+0]+genotype_cal_errors[3*b+1]+genotype_cal_errors[3*b+2]) * 100.0 / (genotype_cal_totals[3*b+0]+genotype_cal_totals[3*b+1]+genotype_cal_totals[3*b+2]);
		fd3 << b << " " << b*step_size << " " << (b+1)*step_size << " " << (b*step_size + step_size/2);
		fd3 << " " << genotype_cal_errors[3*b+0] << " " << genotype_cal_totals[3*b+0];
		fd3 << " " << genotype_cal_errors[3*b+1] << " " << genotype_cal_totals[3*b+1];
		fd3 << " " << genotype_cal_errors[3*b+2] << " " << genotype_cal_totals[3*b+2];
		fd3 << " " << err0 << " " << err1 << " " << err2 << " " << errA << endl;
	}
	fd3.close();

	// [4] Rsquare per bin
	vrb.bullet("Rsquare per frequency bin");
	output_file fd4 (fout + ".rsquare.frq.txt.gz");
	for (int b = 0 ; b < L ; b++) {
		double rsq0 = rsquared_bin[b].corrXY();
		double rsq2 = rsq0*rsq0;
		fd4 << b << " " << frequency_bin[b].size() << " " << frequency_bin[b].mean();
		fd4 << " " << rsq0 << " " << rsq2 << endl;
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
}
