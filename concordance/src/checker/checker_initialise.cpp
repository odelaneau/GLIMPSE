/*******************************************************************************
 * Copyright (C) 2022-2023 Simone Rubinacci
 * Copyright (C) 2022-2023 Olivier Delaneau
 *
 * MIT Licence
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
	double min_val_gl = 0.0;
	int min_val_dp = 0;
	if (options.count("min-val-gl")) min_val_gl=options["min-val-gl"].as < double > ();
	if (options.count("min-val-dp")) min_val_dp=options["min-val-dp"].as < int > ();

	if (options.count("bins")) {
		D.initialize(options["bins"].as < std::vector < double > > (), min_val_gl, min_val_dp);
	} else if (options.count("groups"))
	{
		D.initialize(options["groups"].as < std::string > (),  min_val_gl, min_val_dp);
	} else D.initialize(min_val_gl, min_val_dp);

	//step2: read list of samples
	if (options.count("samples")) D.setTargets(options["samples"].as < std::string > ());

	//step3: read input files
	vrb.title("Reading list of input files in [" + options["input"].as < std::string > () + "]");
	std::string buffer;
	std::vector < std::string > tokens;
	std::vector < std::string > vec_regi, vec_true, vec_esti, vec_freq;
	input_file fd (options["input"].as < std::string > ());
	if (fd.fail()) vrb.error("Cannot open input file: " + options["input"].as < std::string > () + "");
	while (getline(fd, buffer, '\n'))
	{
		if (stb.split(buffer, tokens) == 4)
		{
			vec_regi.push_back(tokens[0]);
			vec_freq.push_back(tokens[1]);
			vec_true.push_back(tokens[2]);
			vec_esti.push_back(tokens[3]);
		}
	}
	if (vec_esti.size()) vrb.bullet(stb.str(vec_esti.size()) + " sets of files detected");
	else vrb.error("No set of files detected in: " + options["input"].as < std::string > () + "");

	std::string out_filename = options["output"].as < std::string > ();
	if (options.count("min-tar-gp"))
	{
		std::vector<float> gp_filters;
		if (options.count("min-tar-gp")) gp_filters = options["min-tar-gp"].as < std::vector <float> > ();
		if (gp_filters.size() == 0) gp_filters.push_back(0.0f);

		for (int i=0; i< gp_filters.size(); ++i)
		{
			float gp_filter = gp_filters[i];
			std::string gp_filter_str = stb.str(gp_filter);
			if (gp_filter < 0 || gp_filter > 1) vrb.error("Values in the --min-tar-gp format must be defined in the range [0-1]. Found value: " + gp_filter_str);

			vrb.title("--- Concordance with GP filter " + gp_filter_str + " started ---");
			D.readData(vec_true, vec_esti, vec_freq, vec_regi, options, gp_filter, out_filename + "_GPfilt_" + gp_filter_str);
			write_files_and_finalise("_GPfilt_" + gp_filter_str);
			vrb.title("--- Concordance with GP filter " + gp_filter_str + " finished - running time = " + stb.str(tac.rel_time()*1.0/1000, 2) + " seconds");
		}
	}
	else
	{
		vrb.title("--- Concordance started ---");
		D.readData(vec_true, vec_esti, vec_freq, vec_regi, options, 0, out_filename);
		write_files_and_finalise();
		vrb.title("--- Concordance finished - running time = " + stb.str(tac.rel_time()*1.0/1000, 2) + " seconds");
	}
	//step2: Measure overall running time
	vrb.title("Total running time = " + stb.str(tac.abs_time()) + " seconds");

}
