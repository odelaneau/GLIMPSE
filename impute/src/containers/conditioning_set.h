/*******************************************************************************
 * Copyright (C) 2018 Olivier Delaneau, University of Lausanne
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
#ifndef _CONDITIONING_SET_H
#define _CONDITIONING_SET_H

#include <utils/otools.h>
#include <containers/variant_map.h>


class conditioning_set {
public:
	//FIXED DATA
	variant_map & mapG;
	unsigned int n_haps;
	unsigned int n_vars;
	unsigned int n_effective;

	//CONDITIONING STATES
	vector < bool > Hpoly;
	vector < bool > Hmono;
	vector < unsigned int > Vpoly;
	vector < unsigned int > Vmono;
	unsigned int n_states;
	unsigned int n_sites;

	//TRANSITION PRBABILITIES
	vector < double > t;
	vector < double > nt;

	//EMISSION
	double ed;
	double ee;

	//CONSTRUCTOR/DESTRUCTOR/INITIALIZATION
	conditioning_set(variant_map & _mapG, unsigned int _n_haps, unsigned int _n_effective) : mapG(_mapG) {
		Hmono.clear();
		Hpoly.clear();
		Vmono.clear();
		Vpoly.clear();
		t.clear();
		nt.clear();
		n_states = 0;
		n_sites = 0;
		n_haps = _n_haps;
		n_vars = mapG.size();
		n_effective= _n_effective;
		ed = 0.0001;
		ee = 0.9999;
	}

	~conditioning_set() {
		Hmono.clear();
		Hpoly.clear();
		Vmono.clear();
		Vpoly.clear();
		t.clear();
		nt.clear();
		n_states = 0;
		n_sites = 0;
	}

	void clear() {
		Hmono.clear();
		Hpoly.clear();
		Vmono.clear();
		Vpoly.clear();
		t.clear();
		nt.clear();
		n_states = 0;
		n_sites = 0;
	}

	void verbose() {
		cout << "\t\tCOND: " << Hmono.size() << " " << Hpoly.size() << " " << Vmono.size() << " " << Vpoly.size() << " " << t.size() << " " << nt.size() << " " << n_states << " " << n_sites << endl;
	}

	void updateTransitions() {
		t = vector < double > (Vpoly.size() - 1, 0.0);
		nt = vector < double > (Vpoly.size() - 1, 0.0);
		for (int l = 1 ; l < Vpoly.size() ; l ++) {
			double rho = 0.04 * n_effective * (mapG.vec_pos[Vpoly[l]]->cm - mapG.vec_pos[Vpoly[l-1]]->cm);
			if (rho == 0.0) rho = 0.00001;
			t[l-1] = -1.0 * expm1(-1.0 * rho / n_haps);
			nt[l-1] = 1-t[l-1];
		}
	}
};

#endif
