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

#include <checker/checker_header.h>

void checker::read_files_and_initialise() {
	//step0: Initialize seed and multi-threading
	rng.setSeed(options["seed"].as < int > ());

	//step1: process bins/groups & parameters
	if (options.count("bins")) {
		D.initialize(options["bins"].as < vector < double > > (), options["minPROB"].as < double > (), options["minDP"].as < int > ());
	} else {
		D.initialize(options["groups"].as < string > (), options["minPROB"].as < double > (), options["minDP"].as < int > ());
	}

	//step2: read list of samples
	if (options.count("samples")) D.setTargets(options["samples"].as < string > ());

	//step3: read input files
	vrb.title("Reading list of input files in [" + options["input"].as < string > () + "]");
	string buffer;
	vector < string > tokens;
	vector < string > vec_regi, vec_true, vec_esti, vec_freq;
	input_file fd (options["input"].as < string > ());
	while (getline(fd, buffer, '\n')) {
		if (stb.split(buffer, tokens) == 4) {
			vec_regi.push_back(tokens[0]);
			vec_freq.push_back(tokens[1]);
			vec_true.push_back(tokens[2]);
			vec_esti.push_back(tokens[3]);
		}
	}
	vrb.bullet(stb.str(vec_esti.size()) + " sets of files detected");
	D.readData(vec_true, vec_esti, vec_freq, vec_regi, options["info_af"].as<string>(),options["thread"].as < int > ());
}
