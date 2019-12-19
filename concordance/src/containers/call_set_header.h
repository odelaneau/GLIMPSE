
#include <utils/otools.h>

class call_set {
public:
	//META
	vector < string > samples;
	vector < int > positions;
	vector < string > chromosomes;
	vector < string > ref, alt;
	vector < string > rsid;
	vector < float > frq;
	vector < bool > flip;
	vector < int > bins;

	//TRUTH
	vector < vector < char > > GTtrue;
	vector < vector < unsigned char > > GLtrue;
	vector < vector < short > > DPtrue;

	//ESTIMATION
	vector < vector < char > > GTesti;
	vector < vector < float > > DSesti;
	vector < vector < float > > GPesti;

	//parameters
	int minDP;
	double minPROB;
	vector < double > BIN;

	//
	call_set ();
	~call_set();

	//
	void set(int, vector < double >, double);
	bool goodForValidation(int, int);
	int getBin(double f);

	//
	void readData(vector < string > &, vector < string > &, vector < string > &, vector < string > &);

	//
	void computeRsquaredPerBin(string output);
	void computeRsquaredPerBinPerSample(string output);

	void computeConcordancePerBIN(string output);
	void concordanceOverall(string output);
	void concordancePerIndividual(string output);
	void computeCalibration(string output);
};

inline bool call_set::goodForValidation(int l, int i) {
	bool isGood = true;
	if (minDP > 0) isGood = isGood && (DPtrue[l][i] >= minDP);
	if (minPROB > 0) {
		double p0 = pow(10, -0.1 * GLtrue[l][3*i+0]);
		double p1 = pow(10, -0.1 * GLtrue[l][3*i+1]);
		double p2 = pow(10, -0.1 * GLtrue[l][3*i+2]);
		double su = p0 + p1 + p2;
		p0 /= su; p1 /= su; p2 /= su;
		//cout << (int)GLtrue[l][3*i+0] << " " << (int)GLtrue[l][3*i+1] << " " << (int)GLtrue[l][3*i+2] << "\t" << p0 << " " << p1 << " " << p2 << endl;
		isGood = isGood && (p0 >= minPROB || p1 >= minPROB || p2 >= minPROB);
	}
	return isGood;
}

inline int call_set::getBin(double f) {
	for (int b = 1 ; b < BIN.size() ; b ++) {
		if (b == 1 && f >= BIN[b-1] && f <= BIN[b]) return b-1;
		if (b != 1 && f > BIN[b-1] && f <= BIN[b]) return b-1;
	}
}
