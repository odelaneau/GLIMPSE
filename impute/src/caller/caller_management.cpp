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
	string input_region = options["input-region"].as < string > ();
	string output_region = options["output-region"].as < string > ();
	vrb.title("Parsing specified genomic regions");
	vector < string > input_t1, input_t2;
	vector < string > output_t1, output_t2;
	int input_ret = stb.split(input_region, input_t1, ":");
	int output_ret = stb.split(output_region, output_t1, ":");
	if (input_ret != 2) vrb.error("Input region needs to be specificied as chrX:Y-Z (chromosome ID cannot be extracted)");
	if (output_ret != 2) vrb.error("Output region needs to be specificied as chrX:Y-Z (chromosome ID cannot be extracted)");
	chrid = input_t1[0];
	if (chrid != output_t1[0]) vrb.error("Chromosome IDs in input and output regions are different!");
	input_ret = stb.split(input_t1[1], input_t2, "-");
	output_ret = stb.split(output_t1[1], output_t2, "-");
	if (input_ret != 2) vrb.error("Input region needs to be specificied as chrX:Y-Z (genomic positions cannot be extracted)");
	if (output_ret != 2) vrb.error("Output region needs to be specificied as chrX:Y-Z (genomic positions cannot be extracted)");
	input_start = atoi(input_t2[0].c_str());
	input_stop = atoi(input_t2[1].c_str());
	output_start = atoi(output_t2[0].c_str());
	output_stop = atoi(output_t2[1].c_str());
	if (input_start >= input_stop) vrb.error("Input genomic region coordinates are incorrect (start >= stop)");
	if (output_start >= output_stop) vrb.error("Output genomic region coordinates are incorrect (start >= stop)");
	if (input_start > output_start) vrb.error("Input/Output genomic region coordinates are imcompatible (input_start > output_start)");
	if (input_stop < output_stop) vrb.error("Input/Output genomic region coordinates are incompatible (input_stop < output_stop)");
	if (input_start < 0) vrb.error("Input genomic region coordinates are incorrect (input_start < 0)");
	if (output_start < 0) vrb.error("Input genomic region coordinates are incorrect (output_start < 0)");
	input_gregion = chrid + ":" + stb.str(input_start) + "-" + stb.str(input_stop);
	output_gregion = chrid + ":" + stb.str(output_start) + "-" + stb.str(output_stop);
	vrb.bullet("Input region  [" + input_gregion + "]");
	vrb.bullet("Output region  [" + output_gregion + "]");
}
