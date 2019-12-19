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
#include <checker/checker_header.h>

void checker::check() {
	int minDP = (options.count("minDP"))?options["minDP"].as < int > ():0;
	double minPROB = (options.count("minPROB"))? options["minPROB"].as < double > ():0.0;
	D.set(minDP, options["bins"].as < vector < double > > (), minPROB);

	D.computeRsquaredPerBin(options["output"].as < string > () + ".rsq.txt");
	D.computeRsquaredPerBinPerSample(options["output"].as < string > () + ".rsq.ind.txt");
	D.computeConcordancePerBIN(options["output"].as < string > () + ".bin.txt");
	D.concordancePerIndividual(options["output"].as < string > () + ".ind.txt");
	D.concordanceOverall(options["output"].as < string > () + ".all.txt");
	D.computeCalibration(options["output"].as < string > () + ".cal.txt");
}
