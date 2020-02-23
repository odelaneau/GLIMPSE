#include "call_set_header.h"

call_set::call_set () {
}

void call_set::initialize(vector < double > _bins, double _T, int _D) {
	vrb.title("Initializing engine based on frequency bins:");
	T = _T;
	D = _D;
	bins = _bins;
	L = bins.size() - 1;
	site2grp.clear();
	rsquared_str.clear();
	vrb.bullet("Probs >= " + stb.str(T));
	vrb.bullet("Depth >= " + stb.str(D));
	vrb.bullet("#intervals in bining = " + stb.str(L));
	vrb.bullet("bins [" + stb.str(bins) + "]");
}

void call_set::initialize(string fgrps, double _T, int _D) {
	vrb.title("Initializing engine based on groups:");
	T = _T;
	D = _D;
	bins.clear();
	L = 0;
	vrb.bullet("Probs >= " + stb.str(T));
	vrb.bullet("Depth >= " + stb.str(D));

	//Read site2groups
	string buffer;
	input_file fd (fgrps);
	vector < string > tokens;
	map < string, pair < int, bool > > :: iterator itG;
	while (getline(fd, buffer, '\n')) {
		if (stb.split(buffer, tokens) > 3) {
			string uuid = tokens[0] + "_" + tokens[1];
			itG = site2grp.find(uuid);
			if (itG == site2grp.end()) {
				//Search for group index
				int grp_idx = -1;
				for (int g  = 0 ; g < rsquared_str.size() && grp_idx < 0; g ++) if (rsquared_str[g] == tokens[2]) grp_idx = g;
				if (grp_idx < 0) {
					grp_idx = rsquared_str.size();
					rsquared_str.push_back(tokens[2]);
				}
				//Add new site of interest
				site2grp.insert(pair < string, pair < int, bool > > (uuid, pair < int, bool > (grp_idx, false)));
			}
		}
	}
	vrb.bullet("#sites  = " + stb.str(site2grp.size()));
	vrb.bullet("#groups = " + stb.str(rsquared_str.size()));
}

call_set::~call_set() {

}
