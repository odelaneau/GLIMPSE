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

call_set::call_set () {
}

void call_set::initialize(vector < double > _bins, double _T, int _D) {
	vrb.title("Initializing engine based on frequency bins:");
	use_subset_samples=false;
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
	use_subset_samples=false;
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
	fd.close();
}

void call_set::setTargets(string fsamples) {
	vrb.title("Reading the subset of samples to consider in the analysis");
	use_subset_samples=true;
	string buffer;
	input_file fd (fsamples);
	vector < string > tokens;

	bool repeated_samples=false;
	while (getline(fd, buffer))
	{
		if (stb.split(buffer, tokens) < 1) vrb.error("Empty line found in sample file.");
		subset_samples.insert(tokens[0]);
	}
	vrb.bullet("#samples  = " + stb.str(subset_samples.size()));
	fd.close();
}

call_set::~call_set() {

}


