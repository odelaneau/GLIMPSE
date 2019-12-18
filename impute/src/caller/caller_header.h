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

#ifndef _CALLER_H
#define _CALLER_H

#include <utils/otools.h>
#include <objects/hmm_parameters.h>

#include <containers/genotype_set.h>
#include <containers/haplotype_set.h>
#include <containers/probability_set.h>
#include <containers/variant_map.h>

#include <models/haplotype_hmm.h>
#include <models/diplotype_hmm.h>
#include <models/switchphasing.h>
#include <models/switchandflipphasing.h>

#define STAGE_INIT	0
#define STAGE_BURN	1
#define STAGE_MAIN	2

class caller {
public:
	//COMMAND LINE OPTIONS
	bpo::options_description descriptions;
	bpo::variables_map options;

	//INTERNAL DATA
	haplotype_set H;
	genotype_set G;
	variant_map V;
	gmap_reader readerGM;
	bool perform_delayed_imputation;

	//GENOMIC REGION
	std::string chrid;
	int input_start;
	int input_stop;
	int output_start;
	int output_stop;
	std::string gregion;

	//MULTI-THREADING
	int i_workers, i_jobs;
	std::vector < pthread_t > id_workers;
	pthread_mutex_t mutex_workers;

	//COMPUTE DATA
	int current_stage;
	basic_stats statH;
	//basic_stats statC;
	std::vector < std::vector < float > > HP0;			// Haplotype posteriors 0
	std::vector < std::vector < float > > HP1;			// Haplotype posteriors 1
	std::vector < std::vector < float > > HLC;			// Conditional haplotype likelihoods
	std::vector < conditioning_set * > COND;		// Conditionning states
	std::vector < probability_set * > TPROB;		// Thresholded posterior probabilities
	std::vector < haplotype_hmm * > HMM;
	std::vector < diplotype_hmm * > DMM;
	std::vector < switchandflipphasing * > FMM;
	std::vector < switchphasing * > SMM;

	//CONSTRUCTOR
	caller();
	~caller();

	//METHODS
	void phase_individual(int, int);
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
	void phase(std::vector < std::string > &);
	void write_files_and_finalise();

	//REGION
	void buildCoordinates();
};


#endif


