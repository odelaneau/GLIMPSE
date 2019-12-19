#include "call_set_header.h"

void call_set::computeConcordancePerBIN(string output) {
	vrb.title ("Computing Rsquared genotype dosages in [" + output + "]");
	bins = vector < int > (positions.size(), -1);

	bins = vector < int > (positions.size(), -1);
	for (int l = 0 ; l < frq.size() ; l ++) bins[l] = getBin(frq[l]);

	vector < float > bfrq = vector < float >(BIN.size()-1, 0.0);
	vector < int > cfrq = vector < int >(BIN.size()-1, 0);
	vector < int > err_RR = vector < int >(BIN.size()-1, 0);
	vector < int > err_RA = vector < int >(BIN.size()-1, 0);
	vector < int > err_AA = vector < int >(BIN.size()-1, 0);
	vector < int > tot_RR = vector < int >(BIN.size()-1, 0);
	vector < int > tot_RA = vector < int >(BIN.size()-1, 0);
	vector < int > tot_AA = vector < int >(BIN.size()-1, 0);
	for (int l = 0 ; l < GTtrue.size() ; l ++) {
		bool add = false;
		for (int i = 0 ; i < GTtrue[l].size() ; i ++) {
			if (goodForValidation(l, i)) {
				int dsT = (int)GTtrue[l][i];
				switch (dsT) {
				case 0:	err_RR[bins[l]] += (GTesti[l][i] != GTtrue[l][i]);
						tot_RR[bins[l]] ++;
						break;
				case 1:	err_RA[bins[l]] += (GTesti[l][i] != GTtrue[l][i]);
						tot_RA[bins[l]] ++;
						break;
				case 2:	err_AA[bins[l]] += (GTesti[l][i] != GTtrue[l][i]);
						tot_AA[bins[l]] ++;
						break;
				}
				add = true;
			}
		}
		if (add) {
			bfrq[bins[l]] += frq[l];
			cfrq[bins[l]] ++;
		}
	}
	output_file fd (output);
	for (int b = 0 ; b < BIN.size()-1 ; b ++)
		if (cfrq[b] > 0)
			fd << cfrq[b] << " " << bfrq[b] << " " << err_RR[b] * 100.0 / tot_RR[b] << " " << err_RA[b] * 100.0 / tot_RA[b] << " " << err_AA[b] * 100.0 / tot_AA[b] << endl;
	fd.close();
	vrb.bullet("done");
}

void call_set::concordanceOverall(string output) {
	vrb.title ("Computing overall genotype concordance in [" + output + "]");
	int err_RR = 0;
	int tot_RR = 0;
	int err_RA = 0;
	int tot_RA = 0;
	int err_AA = 0;
	int tot_AA = 0;
	for (int l = 0 ; l < GTtrue.size() ; l ++) {
		for (int i = 0 ; i < GTtrue[l].size() ; i ++) {
			if (goodForValidation(l, i)) {
				switch (GTtrue[l][i]) {
				case 0:	err_RR += (GTesti[l][i] != GTtrue[l][i]);
						tot_RR ++;
						break;
				case 1:	err_RA += (GTesti[l][i] != GTtrue[l][i]);
						tot_RA ++;
						break;
				case 2:	err_AA += (GTesti[l][i] != GTtrue[l][i]);
						tot_AA ++;
						break;
				}
			}
		}
	}
	output_file fd (output);
	fd << " " << tot_RR << " " << err_RR << " " << tot_RA << " " << err_RA << " " << tot_AA << " " << err_AA << endl;
	fd.close();

	vrb.bullet ("RR = " + stb.str(tot_RR) + " errors = " + stb.str(err_RR) + " rate = " + stb.str(err_RR * 100.0 / tot_RR, 3) + "%");
	vrb.bullet ("RA = " + stb.str(tot_RA) + " errors = " + stb.str(err_RA) + " rate = " + stb.str(err_RA * 100.0 / tot_RA, 3) + "%");
	vrb.bullet ("AA = " + stb.str(tot_AA) + " errors = " + stb.str(err_AA) + " rate = " + stb.str(err_AA * 100.0 / tot_AA, 3) + "%");
}

void call_set::concordancePerIndividual(string output) {
	vrb.title ("Computing overall genotype concordance in [" + output + "]");
	output_file fd (output);
	for (int i = 0 ; i < samples.size() ; i ++) {
		int err_RR = 0;
		int tot_RR = 0;
		int err_RA = 0;
		int tot_RA = 0;
		int err_AA = 0;
		int tot_AA = 0;
		for (int l = 0 ; l < GTtrue.size() ; l ++) {
			if (goodForValidation(l, i)) {
				switch (GTtrue[l][i]) {
				case 0:	err_RR += (GTesti[l][i] != GTtrue[l][i]);
						tot_RR ++;
						break;
				case 1:	err_RA += (GTesti[l][i] != GTtrue[l][i]);
						tot_RA ++;
						break;
				case 2:	err_AA += (GTesti[l][i] != GTtrue[l][i]);
						tot_AA ++;
						break;
				}
			}
		}
		fd << samples[i] << " " << tot_RR << " " << err_RR << " " << tot_RA << " " << err_RA << " " << tot_AA << " " << err_AA << endl;
	}
	fd.close();
}
