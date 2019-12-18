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
#include <caller/caller_header.h>

caller::caller() {
}

caller::~caller() {
}

void caller::phase(vector < string > & args) {
	declare_options();
	parse_command_line(args);
	check_options();
	verbose_files();
	verbose_options();
	read_files_and_initialise();
	phase_loop();
	write_files_and_finalise();
}

void caller::buildCoordinates() {
	string region = options["region"].as < string > ();
	int sbuffer = options["buffer"].as < int > () * 1000;
	vrb.title("Parsing specified genomic region");
	vector < string > t1, t2;
	int ret1 = stb.split(region, t1, ":");
	if (ret1 == 2) {
		chrid = t1[0];
		int ret2 = stb.split(t1[1], t2, "-");
		if (ret2 != 2) vrb.error("Genomic region incorrectly specified (case 1)");
		output_start = atoi(t2[0].c_str());
		output_stop = atoi(t2[1].c_str());
		if (output_start >= output_stop) vrb.error("Genomic region incorrectly specified (case 2)");
		input_start = output_start - sbuffer;
		input_stop = output_stop + sbuffer;
		if (input_start < 0) input_start = 0;
		if (input_start >= input_stop) vrb.error("Genomic region incorrectly specified (case 3)");
		gregion = chrid + ":" + stb.str(input_start) + "-" + stb.str(input_stop);
		vrb.bullet("Input region  [" + gregion + "]");
		vrb.bullet("Output region [" + options["region"].as < string > () + "]");
	} else if (ret1 == 1) {
		chrid = t1[0];
		input_start = 0;
		input_stop = 1000000000;
		output_start = 0;
		output_stop = 1000000000;
		gregion = chrid;
		vrb.bullet("Input region  [" + gregion + "]");
		vrb.bullet("Output region [" + gregion + "]");
	} else vrb.error("Genomic region incorrectly specified (case 4)");
}
