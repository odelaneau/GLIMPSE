#include "call_set_header.h"

float getMAX(float a, float b, float c) {
	float d = (a<b)?b:a;
	return ((d<c)?c:d);
}


void call_set::computeCalibration(string output) {
	vrb.title ("Computing calibration of genotype posteriors in [" + output + "]");
	int nBIN = 100;
	vector < int > errGT = vector < int > (nBIN, 0);
	vector < int > totGT = vector < int > (nBIN, 0);
	for (int l = 0 ; l < GTtrue.size() ; l ++) {
		for (int i = 0 ; i < GTtrue[l].size() ; i ++) {
			if (goodForValidation(l, i)) {
				float gp0 = GPesti[l][3*i+0];
				float gp1 = GPesti[l][3*i+1];
				float gp2 = GPesti[l][3*i+2];
				int dsT = (int)GTtrue[l][i];
				float mmm = getMAX(gp0, gp1, gp2);
				int bin = (int)trunc( mmm * (nBIN-1));

				switch (dsT) {
				case 0:	errGT[bin] += ((gp0 <= gp1) || (gp0 <= gp2));
						totGT[bin] ++;
						break;
				case 1:	errGT[bin] += ((gp1 <= gp0) || (gp1 <= gp2));
						totGT[bin] ++;
						break;
				case 2:	errGT[bin] += ((gp2 <= gp0) || (gp2 <= gp1));
						totGT[bin] ++;
						break;
				}
			}
		}
	}

	output_file fd (output);
	for (int b = 0 ; b < nBIN ; b ++) fd << b*1.0/nBIN << " " << (b+1)*1.0/nBIN << " " << errGT[b] << " " << totGT[b] << endl;
	fd.close();
	vrb.bullet("done");
}
