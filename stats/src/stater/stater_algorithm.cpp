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


#include <stater/stater_header.h>

void stater::stats() {
	//Initialize
	rng.setSeed(options["seed"].as < int > ());

	//Read input files and regions
	vrb.title("Reading list of input files in [" + options["input"].as < string > () + "]");
	string buffer;
	vector < string > files, regions, tokens;
	input_file fd(options["input"].as < string > ());
	while (getline(fd, buffer)) {
		if (stb.split(buffer, tokens) < 2) {
			regions.push_back(tokens[0]);
			files.push_back(tokens[1]);
		}
	}
	vrb.bullet(stb.str(files.size()) + " files detected");

	//Read and Process data
	readDataAndComputeStats(files, regions);

	//Output statistics per-variant
	vrb.title("Writing per-variant statistics in [" + options["output"].as < string > () + ".var.bed]");
	output_file fdo1(options["output"].as < string > () + ".var.bed");
	for (int l = 0 ; l  < var_chr.size() ; l ++) {
		fdo1 << var_chr[l];
		fdo1 << "\t" << var_pos[l]-1;
		fdo1 << "\t" << var_pos[l];
		fdo1 << "\t" << var_ids[l];
		fdo1 << "\t" << var_ref[l];
		fdo1 << "\t" << var_alt[l];
		fdo1 << "\t" << stb.str(var_freq[l], 5);
		fdo1 << "\t" << stb.str(min(var_freq[l], 1.0f - var_freq[l]), 5);
		fdo1 << "\t" << stb.str(var_info[l], 5);
		fdo1 << "\t" << stb.str(var_avgp[l], 5) << endl;
	}
	fdo1.close();

	//Output statistics per-sample
	vrb.title("Writing per-sample statistics in [" + options["output"].as < string > () + ".spl.txt]");
	output_file fdo2(options["output"].as < string > () + ".spl.txt");
	for (int i = 0 ; i  < spl_ids.size() ; i ++) {
		fdo2 << spl_ids[i];
		fdo2 << "\t" << stb.str(spl_info[i], 5);
		fdo2 << "\t" << stb.str(spl_avgp[i], 5) << endl;
	}
	fdo2.close();

	//Finalize
	vrb.title("Total running time = " + stb.str(tac.abs_time()) + " seconds");
}
