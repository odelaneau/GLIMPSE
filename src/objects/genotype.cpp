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

#include <objects/genotype.h>

const float unphred[256] = { 1.000000e+00, 7.943282e-01, 6.309573e-01, 5.011872e-01, 3.981072e-01, 3.162278e-01, 2.511886e-01, 1.995262e-01, 1.584893e-01, 1.258925e-01, 1.000000e-01, 7.943282e-02, 6.309573e-02, 5.011872e-02, 3.981072e-02, 3.162278e-02, 2.511886e-02, 1.995262e-02, 1.584893e-02, 1.258925e-02, 1.000000e-02, 7.943282e-03, 6.309573e-03, 5.011872e-03, 3.981072e-03, 3.162278e-03, 2.511886e-03, 1.995262e-03, 1.584893e-03, 1.258925e-03, 1.000000e-03, 7.943282e-04, 6.309573e-04, 5.011872e-04, 3.981072e-04, 3.162278e-04, 2.511886e-04, 1.995262e-04, 1.584893e-04, 1.258925e-04, 1.000000e-04, 7.943282e-05, 6.309573e-05, 5.011872e-05, 3.981072e-05, 3.162278e-05, 2.511886e-05, 1.995262e-05, 1.584893e-05, 1.258925e-05, 1.000000e-05, 7.943282e-06, 6.309573e-06, 5.011872e-06, 3.981072e-06, 3.162278e-06, 2.511886e-06, 1.995262e-06, 1.584893e-06, 1.258925e-06, 1.000000e-06, 7.943282e-07, 6.309573e-07, 5.011872e-07, 3.981072e-07, 3.162278e-07, 2.511886e-07, 1.995262e-07, 1.584893e-07, 1.258925e-07, 1.000000e-07, 7.943282e-08, 6.309573e-08, 5.011872e-08, 3.981072e-08, 3.162278e-08, 2.511886e-08, 1.995262e-08, 1.584893e-08, 1.258925e-08, 1.000000e-08, 7.943282e-09, 6.309573e-09, 5.011872e-09, 3.981072e-09, 3.162278e-09, 2.511886e-09, 1.995262e-09, 1.584893e-09, 1.258925e-09, 1.000000e-09, 7.943282e-10, 6.309573e-10, 5.011872e-10, 3.981072e-10, 3.162278e-10, 2.511886e-10, 1.995262e-10, 1.584893e-10, 1.258925e-10, 1.000000e-10, 7.943282e-11, 6.309573e-11, 5.011872e-11, 3.981072e-11, 3.162278e-11, 2.511886e-11, 1.995262e-11, 1.584893e-11, 1.258925e-11, 1.000000e-11, 7.943282e-12, 6.309573e-12, 5.011872e-12, 3.981072e-12, 3.162278e-12, 2.511886e-12, 1.995262e-12, 1.584893e-12, 1.258925e-12, 1.000000e-12, 7.943282e-13, 6.309573e-13, 5.011872e-13, 3.981072e-13, 3.162278e-13, 2.511886e-13, 1.995262e-13, 1.584893e-13, 1.258925e-13, 1.000000e-13, 7.943282e-14, 6.309573e-14, 5.011872e-14, 3.981072e-14, 3.162278e-14, 2.511886e-14, 1.995262e-14, 1.584893e-14, 1.258925e-14, 1.000000e-14, 7.943282e-15, 6.309573e-15, 5.011872e-15, 3.981072e-15, 3.162278e-15, 2.511886e-15, 1.995262e-15, 1.584893e-15, 1.258925e-15, 1.000000e-15, 7.943282e-16, 6.309573e-16, 5.011872e-16, 3.981072e-16, 3.162278e-16, 2.511886e-16, 1.995262e-16, 1.584893e-16, 1.258925e-16, 1.000000e-16, 7.943282e-17, 6.309573e-17, 5.011872e-17, 3.981072e-17, 3.162278e-17, 2.511886e-17, 1.995262e-17, 1.584893e-17, 1.258925e-17, 1.000000e-17, 7.943282e-18, 6.309573e-18, 5.011872e-18, 3.981072e-18, 3.162278e-18, 2.511886e-18, 1.995262e-18, 1.584893e-18, 1.258925e-18, 1.000000e-18, 7.943282e-19, 6.309573e-19, 5.011872e-19, 3.981072e-19, 3.162278e-19, 2.511886e-19, 1.995262e-19, 1.584893e-19, 1.258925e-19, 1.000000e-19, 7.943282e-20, 6.309573e-20, 5.011872e-20, 3.981072e-20, 3.162278e-20, 2.511886e-20, 1.995262e-20, 1.584893e-20, 1.258925e-20, 1.000000e-20, 7.943282e-21, 6.309573e-21, 5.011872e-21, 3.981072e-21, 3.162278e-21, 2.511886e-21, 1.995262e-21, 1.584893e-21, 1.258925e-21, 1.000000e-21, 7.943282e-22, 6.309573e-22, 5.011872e-22, 3.981072e-22, 3.162278e-22, 2.511886e-22, 1.995262e-22, 1.584893e-22, 1.258925e-22, 1.000000e-22, 7.943282e-23, 6.309573e-23, 5.011872e-23, 3.981072e-23, 3.162278e-23, 2.511886e-23, 1.995262e-23, 1.584893e-23, 1.258925e-23, 1.000000e-23, 7.943282e-24, 6.309573e-24, 5.011872e-24, 3.981072e-24, 3.162278e-24, 2.511886e-24, 1.995262e-24, 1.584893e-24, 1.258925e-24, 1.000000e-24, 7.943282e-25, 6.309573e-25, 5.011872e-25, 3.981072e-25, 3.162278e-25, 2.511886e-25, 1.995262e-25, 1.584893e-25, 1.258925e-25, 1.000000e-25, 7.943282e-26, 6.309573e-26, 5.011872e-26, 3.981072e-26, 3.162278e-26};

genotype::genotype(int _index, int _n_variants) {
	name = "";
	index = _index;
	n_variants = _n_variants;
}

genotype::~genotype() {
	name = "";
	index = -1;
	n_variants = 0;
	free();
}

void genotype::free() {
	GL.clear();
	GP.clear();
	H0.clear();
	H1.clear();
}

void genotype::allocate() {
	GL = vector < unsigned char > (3 * n_variants, 0);
	GP = vector < float > (2 * n_variants, 0.0);
	H0 = vector < bool > (n_variants, false);
	H1 = vector < bool > (n_variants, false);
}

void genotype::initHaplotypeLikelihoods(vector < float > & HL) {
	//double sum, ph0, ph1, lh0, lh1, lg0, lg1, lg2;
	double sum, lg0, lg1, lg2;
	for (int l = 0 ; l < n_variants ; l ++) {
		lg0 = unphred[GL[3*l+0]];
		lg1 = unphred[GL[3*l+1]];
		lg2 = unphred[GL[3*l+2]];
		sum = lg0+lg1+lg2;
		lg0 /= sum;
		lg1 /= sum;
		lg2 /= sum;
		HL[2*l+0] = lg0 + 0.5 * lg1;
		HL[2*l+1] = lg2 + 0.5 * lg1;
		/*
		lh0 = lg0 + 0.5 * lg1;
		lh1 = lg2 + 0.5 * lg1;
		ph0 = ee * lh0 + ed * lh1;
		ph1 = ed * lh0 + ee * lh1;
		HL[2*l+0] = ph0 / (ph0 + ph1);
		HL[2*l+1] = ph1 / (ph0 + ph1);
		*/
	}
}

void genotype::makeHaplotypeLikelihoods(vector < float > & HL, bool first) {
	//double sum, ph0, ph1, lh0, lh1, lg0, lg1, lg2;
	double sum, lg0, lg1, lg2;
	for (int l = 0 ; l < n_variants ; l ++) {
		bool condAllele = first?H1[l]:H0[l];
		lg0 = unphred[GL[3*l+0]];
		lg1 = unphred[GL[3*l+1]];
		lg2 = unphred[GL[3*l+2]];
		sum = lg0+lg1+lg2;
		lg0 /= sum;
		lg1 /= sum;
		lg2 /= sum;
		if (condAllele) {
			HL[2*l+0] = lg1 / (lg1 + lg2);
			HL[2*l+1] = lg2 / (lg1 + lg2);
		} else {
			HL[2*l+0] = lg0 / (lg0 + lg1);
			HL[2*l+1] = lg1 / (lg0 + lg1);
		}

		/*
		if (condAllele) {
			lh0 = lg1 / (lg1 + lg2);
			lh1 = lg2 / (lg1 + lg2);
		} else {
			lh0 = lg0 / (lg0 + lg1);
			lh1 = lg1 / (lg0 + lg1);
		}
		ph0 = ee * lh0 + ed * lh1;
		ph1 = ed * lh0 + ee * lh1;
		HL[2*l+0] = ph0 / (ph0 + ph1);
		HL[2*l+1] = ph1 / (ph0 + ph1);
		*/
	}
}

void genotype::sampleHaplotypeH0(vector < float > & HP0) {
	for (int l = 0 ; l < n_variants ; l ++) {
		if (rng.getDouble() < HP0[2*l+0]) H0[l] = false;
		else H0[l] = true;
	}
}

void genotype::sampleHaplotypeH1(vector < float > & HP1) {
	for (int l = 0 ; l < n_variants ; l ++) {
		if (rng.getDouble() < HP1[2*l+0]) H1[l] = false;
		else H1[l] = true;
	}
}

void genotype::storeGenotypePosteriors(vector < float > & HP0, vector < float > & HP1) {
	for (int l = 0 ; l < n_variants ; l ++) {
		double p0 = HP0[2*l+0] * HP1[2*l+0];
		double p1 = HP0[2*l+0] * HP1[2*l+1] + HP0[2*l+1] * HP1[2*l+0];
		double p2 = HP0[2*l+1] * HP1[2*l+1];
		GP[2*l+0] += p0 / (p0 + p1 + p2);
		GP[2*l+1] += p1 / (p0 + p1 + p2);
	}
}

void genotype::normGenotypePosteriors(int nStorages) {
	for (int l = 0 ; l < n_variants ; l ++) {
		GP[2*l+0] /= nStorages;
		GP[2*l+1] /= nStorages;
	}
}

void genotype::inferGenotypes() {
	for (int l = 0 ; l < n_variants ; l ++) {
		float gp0 = GP[2*l+0];
		float gp1 = GP[2*l+1];
		float gp2 = 1.0 - gp0 - gp1;
		if (gp0 >= gp1 && gp0 >= gp2) { H0[l] = false; H1[l] = false; }
		if (gp1 >= gp0 && gp1 >= gp2) { H0[l] = false; H1[l] = true; }
		if (gp2 >= gp0 && gp2 >= gp1) { H0[l] = true; H1[l] = true; }
	}
}
