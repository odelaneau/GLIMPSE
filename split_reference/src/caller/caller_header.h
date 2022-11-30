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

#ifndef _CALLER_H
#define _CALLER_H

#include <utils/otools.h>

#include <containers/variant_map.h>
#include <io/ref_genotype_reader.h>
#include "../containers/ref_haplotype_set.h"

#define STAGE_INIT	0
#define STAGE_BURN	1
#define STAGE_MAIN	2

class caller {
public:
	//COMMAND LINE OPTIONS
	bpo::options_description descriptions;
	bpo::variables_map options;

	std::vector<std::string> chrid;
	std::vector<int> input_start;
	std::vector<int> input_stop;
	std::vector<std::string> input_gregion;
	std::vector<int> output_start;
	std::vector<int> output_stop;
	std::vector<std::string> output_gregion;

	//MULTI-THREADING
	int i_workers, i_jobs;
	std::vector < pthread_t > id_workers;
	pthread_mutex_t mutex_workers;

	//CONSTRUCTOR
	caller();
	~caller();

	//METHODS
	void phase_individual(const int, const int);
	void phase_iteration();
	void phase_loop();

	//PARAMETERS
	void declare_options();
	void parse_command_line(std::vector < std::string > &);
	void check_options();
	void verbose_options();
	void verbose_files();

	//FILE I/O
	void read_files_and_initialise();
	void setup_mpileup();
	void read_BAMs();

	void phase(std::vector < std::string > &);
	void write_files_and_finalise();

	//REGION
	void buildCoordinates();
};


#endif


