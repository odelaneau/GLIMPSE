#include "call_set_header.h"

float squaredCorrelationCoefficient(vector < float > & X, vector < float > & Y) {
	double vmeanX = 0.0, vmeanY = 0.0;
	int cmeanX = 0.0, cmeanY = 0.0;

	for (int i = 0; i < X.size(); i++) {
		vmeanX += X[i];
		vmeanY += Y[i];
		cmeanX ++;
		cmeanY ++;
	}

	vmeanX /= cmeanX;
	vmeanY /= cmeanY;

	double sumX = 0.0, sumY = 0.0, sumXY = 0.0;
	for (int i = 0; i < X.size(); i++) {
		sumX += (X[i] - vmeanX) * (X[i] - vmeanX);
		sumY += (Y[i] - vmeanY) * (Y[i] - vmeanY);
		sumXY += (X[i] - vmeanX) * (Y[i] - vmeanY);
	}

    float correlation = sumXY / (sqrt(sumX) * sqrt(sumY));
    return correlation * correlation;
}

void call_set::computeRsquaredPerBin(string output) {
	vrb.title ("Computing Rsquared genotype dosages in [" + output + "]");

	bins = vector < int > (positions.size(), -1);
	for (int l = 0 ; l < frq.size() ; l ++) bins[l] = getBin(frq[l]);

	vector < float > bfrq = vector < float >(BIN.size()-1, 0.0);
	vector < int > cfrq = vector < int >(BIN.size()-1, 0);
	vector < vector < float > > trueDS = vector < vector < float > > (BIN.size()-1);
	vector < vector < float > > estiDS = vector < vector < float > > (BIN.size()-1);
	for (int l = 0 ; l < GTtrue.size() ; l ++) {
		bool add = false;
		for (int i = 0 ; i < GTtrue[l].size() ; i ++) {
			if (goodForValidation(l, i)) {
				int dsT = (int)GTtrue[l][i];
				float dsP = DSesti[l][i];
				if (flip[l]) {
					dsT = 2.0 - dsT;
					dsP = 2.0 - dsP;
				}

				trueDS[bins[l]].push_back((float)dsT);
				estiDS[bins[l]].push_back(dsP);
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
			fd << cfrq[b] << " " << bfrq[b]/cfrq[b] << " " << squaredCorrelationCoefficient(trueDS[b], estiDS[b]) << endl;
	fd.close();
	vrb.bullet("done");
}

void call_set::computeRsquaredPerBinPerSample(string output) {
	vrb.title ("Computing Rsquared genotype dosages in [" + output + "]");

	bins = vector < int > (positions.size(), -1);
	for (int l = 0 ; l < frq.size() ; l ++) bins[l] = getBin(frq[l]);

	output_file fd (output);

	for (int i = 0 ; i < samples.size() ; i ++) {
		vector < float > bfrq = vector < float >(BIN.size()-1, 0.0);
		vector < int > cfrq = vector < int >(BIN.size()-1, 0);
		vector < vector < float > > trueDS = vector < vector < float > > (BIN.size()-1);
		vector < vector < float > > estiDS = vector < vector < float > > (BIN.size()-1);
		for (int l = 0 ; l < GTtrue.size() ; l ++) {
			bool add = false;
			if (goodForValidation(l, i)) {
				int dsT = (int)GTtrue[l][i];
				float dsP = DSesti[l][i];
				if (flip[l]) {
					dsT = 2.0 - dsT;
					dsP = 2.0 - dsP;
				}

				trueDS[bins[l]].push_back((float)dsT);
				estiDS[bins[l]].push_back(dsP);
				add = true;
			}
			if (add) {
				bfrq[bins[l]] += frq[l];
				cfrq[bins[l]] ++;
			}
		}
		for (int b = 0 ; b < BIN.size()-1 ; b ++)
			if (cfrq[b] > 0)
				fd << samples[i] << " " << cfrq[b] << " " << bfrq[b]/cfrq[b] << " " << squaredCorrelationCoefficient(trueDS[b], estiDS[b]) << endl;
	}
	fd.close();
	vrb.bullet("done");
}
