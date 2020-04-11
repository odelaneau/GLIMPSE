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

#include "../../versions/versions.h"
#include <snparray/snparray_header.h>

void snparray::declare_options() {
	bpo::options_description opt_base ("Basic options");
	opt_base.add_options()
			("help", "Produce help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator");

	bpo::options_description opt_input ("Input files");
	opt_input.add_options()
			("input", bpo::value< string > (), "Input Genotypes")
			("sites", bpo::value< string > (), "Sites to be kept in output")
			("region", bpo::value< string >(), "Target regions");

	bpo::options_description opt_algo ("Parameters");
	opt_algo.add_options()
			("minPROB", bpo::value < double >(), "Prob for calling")
			("minDP", bpo::value < int >(), "Minimal depth of coverage for calling")
			("minMIS", bpo::value < double >(), "Minimal missing data rate post calling for variant inclusion");

	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("output", bpo::value< string >(), "Output file")
			("log", bpo::value< string >(), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_algo).add(opt_output);
}

void snparray::parse_command_line(vector < string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { cerr << "Error parsing command line arguments: " << string(e.what()) << endl; exit(0); }

	if (options.count("help")) { cout << descriptions << endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < string > () +"]");

	vrb.title("[GLIMPSE] Simulation of SNP array from high coverage sequencing data");
	vrb.bullet("Author        : Simone RUBINACCI & Olivier DELANEAU, University of Lausanne");
	vrb.bullet("Contact       : simone.rubinacci@unil.ch & olivier.delaneau@unil.ch");
	vrb.bullet("Version       : " + string(SNPARRAY_VERSION));
	vrb.bullet("Run date      : " + tac.date());
}

void snparray::check_options() {
	if (!options.count("input"))
		vrb.error("You must specify --input");

	if (!options.count("sites"))
		vrb.error("You must specify --sites");

	if (!options.count("output"))
		vrb.error("You must specify --output");

	if (!options.count("region"))
		vrb.error("You must specify --region");

	if (!options.count("minPROB"))
		vrb.error("You must specify --minPROB");

	if (!options.count("minDP"))
		vrb.error("You must specify --minDP");

	if (!options.count("minMISS"))
		vrb.error("You must specify --minMISS");
}

void snparray::verbose_files() {
	vrb.title("Files:");
	vrb.bullet("Input    : " + stb.str(options["input"].as < string > ()));
	vrb.bullet("Sites    : " + stb.str(options["sites"].as < string > ()));
	vrb.bullet("Output   : " + stb.str(options["output"].as < string > ()));
	vrb.bullet("Region   : " + stb.str(options["region"].as < string > ()));
}

void snparray::verbose_options() {
	vrb.title("Parameters:");
	vrb.bullet("Seed    : " + stb.str(options["seed"].as < int > ()));
	vrb.bullet("MinPROB : " + stb.str(options["minPROB"].as < double > ()));
	vrb.bullet("MinDP   : " + stb.str(options["minDP"].as < int > ()));
	vrb.bullet("MinMISS : " + stb.str(options["minMIS"].as < double > ()));
}
