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

#include <objects/hmm_parameters.h>

hmm_parameters::hmm_parameters() {
	ed = 0.0001;
	ee = 0.9999;
}

hmm_parameters::~hmm_parameters() {
	t.clear();
	nt.clear();
}

void hmm_parameters::initialise(variant_map & mapG, haplotype_set & H, int Neff) {
	t = vector < double > (H.CVidx.size() - 1, 0.0f);
	nt = vector < double > (H.CVidx.size() - 1, 0.0f);
	for (int l = 1 ; l < H.CVidx.size() ; l ++) {
		float dist_cm = (mapG.vec_pos[H.CVidx[l]]->cm - mapG.vec_pos[H.CVidx[l-1]]->cm);
		if (dist_cm < 1e-7) dist_cm = 1e-7;
		double rho = 0.04 * Neff * dist_cm;
		t[l-1] = -1.0 * expm1(-1.0 * rho / H.n_hap);
		nt[l-1] = 1-t[l-1];
	}
}

