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

#include <caller/caller_header.h>

caller::caller():
		i_workers(0), i_jobs(0), current_stage(0), min_gl(1e-10),
		input_fmt(InputFormat::BCF), output_fmt(OutputFormat::BCF), output_compr(OutputCompression::ZLIB), bgen_bits(8)
{
}

caller::~caller()
{
	for (int t = 0; t < HMM.size(); t++)
		delete HMM[t];
	HMM.clear();
	for (int t = 0; t < DMM.size(); t++)
		delete DMM[t];
	DMM.clear();
	for (int t = 0; t < COND.size(); t++)
		delete COND[t];
	COND.clear();
}

void caller::phase(std::vector < std::string > & args) {
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
	std::string input_region = options["input-region"].as < std::string > ();
	std::string output_region = options["output-region"].as < std::string > ();
	vrb.title("Parsing specified genomic regions");
	std::vector < std::string > input_t1, input_t2;
	std::vector < std::string > output_t1, output_t2;
	int input_ret = stb.split(input_region, input_t1, ":");
	int output_ret = stb.split(output_region, output_t1, ":");
	if (input_ret != 2) vrb.error("Input region needs to be specificied as chrX:Y-Z (chromosome ID cannot be extracted)");
	if (output_ret != 2) vrb.error("Output region needs to be specificied as chrX:Y-Z (chromosome ID cannot be extracted)");
	V.chrid = input_t1[0];
	if (V.chrid != output_t1[0]) vrb.error("Chromosome IDs in input and output regions are different!");
	input_ret = stb.split(input_t1[1], input_t2, "-");
	output_ret = stb.split(output_t1[1], output_t2, "-");
	if (input_ret != 2) vrb.error("Input region needs to be specificied as chrX:Y-Z (genomic positions cannot be extracted)");
	if (output_ret != 2) vrb.error("Output region needs to be specificied as chrX:Y-Z (genomic positions cannot be extracted)");
	V.input_start = atoi(input_t2[0].c_str());
	V.input_stop = atoi(input_t2[1].c_str());
	V.output_start = atoi(output_t2[0].c_str());
	V.output_stop = atoi(output_t2[1].c_str());
	if (V.input_start >= V.input_stop) vrb.error("Input genomic region coordinates are incorrect (start >= stop)");
	if (V.output_start >= V.output_stop) vrb.error("Output genomic region coordinates are incorrect (start >= stop)");
	if (V.input_start > V.output_start) vrb.error("Input/Output genomic region coordinates are imcompatible (input_start > output_start)");
	if (V.input_stop < V.output_stop) vrb.error("Input/Output genomic region coordinates are incompatible (input_stop < output_stop)");
	if (V.input_start < 0) vrb.error("Input genomic region coordinates are incorrect (input_start < 0)");
	if (V.output_start < 0) vrb.error("Input genomic region coordinates are incorrect (output_start < 0)");
	V.input_gregion = V.chrid + ":" + stb.str(V.input_start) + "-" + stb.str(V.input_stop);
	V.output_gregion = V.chrid + ":" + stb.str(V.output_start) + "-" + stb.str(V.output_stop);
	vrb.bullet("Input region  [" + V.input_gregion + "]");
	vrb.bullet("Output region  [" + V.output_gregion + "]");
}
