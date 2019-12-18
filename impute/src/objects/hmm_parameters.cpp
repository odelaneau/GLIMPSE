////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2018 Olivier Delaneau, University of Lausanne
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////

#include <objects/hmm_parameters.h>

hmm_parameters::hmm_parameters() {
	ed = 0.0001;
	ee = 0.9999;
}

hmm_parameters::~hmm_parameters() {
	t.clear();
	nt.clear();
	tfreq.clear();
}

void hmm_parameters::initialise(variant_map & mapG, int Neff, int Nhap, int Nstate, bool gmap) {
	t = std::vector < double > (mapG.size() - 1, 0.0);
	nt = std::vector < double > (mapG.size() - 1, 0.0);
	tfreq = std::vector < double > (mapG.size() - 1, 0.0);
	tac.clock();
	for (int l = 1 ; l < mapG.size() ; l ++) {
		double rho;
		if (gmap) rho = 0.04 * Neff * (mapG.vec_pos[l]->cm - mapG.vec_pos[l-1]->cm);
		else rho = 0.0002 *  (mapG.vec_pos[l]->bp - mapG.vec_pos[l-1]->bp);
		if (rho == 0.0) rho = 0.00001;
		nt[l-1] = exp(-1.0 * rho / Nhap);
		t[l-1] = 1-nt[l-1];
		tfreq[l-1] = t[l-1] * 1.0 / Nstate;
		//t[l-1] = -1.0 * expm1(-1.0 * rho / Nhap);
		//nt[l-1] = 1.0 - t[l-1];
		//tfreq[l-1] = t[l-1] * 1.0 / Nstate;
	}
	vrb.bullet("Update map (" + stb.str(tac.rel_time()*1.0, 1) + "ms)");
	//emissions
	double theta = 0.0;
	for (int i = 1 ; i < (Nhap - 1) ; i ++) theta += 1.0/i;
	theta = 1.0 / theta;
	ed = ( theta / ( Nhap + theta ) ) * 0.5;
	ee = ed + ( Nhap / ( Nhap + theta ) );
}

