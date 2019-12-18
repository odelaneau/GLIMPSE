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
//$Id: hmm_parameters.h 597 2012-07-18 14:00:06Z koskos $

#ifndef _HMM_PARAMETERS_H
#define _HMM_PARAMETERS_H

#include <utils/otools.h>
#include <containers/variant_map.h>

class hmm_parameters {
public :
	//DATA
	std::vector < double > t;
	std::vector < double > tfreq;
	std::vector < double > nt;
	double ee;
	double ed;
	double efreq;
	double dfreq;

	//CONSTRUCTOR/DESTRUCTOR
	hmm_parameters();
	~hmm_parameters();

	//METHODS
	void initialise(variant_map &, int, int, int, bool);
};

#endif
