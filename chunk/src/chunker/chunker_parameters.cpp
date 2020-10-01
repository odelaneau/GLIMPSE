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
#include <chunker/chunker_header.h>

void chunker::declare_options() {
	bpo::options_description opt_base ("Basic options");
	opt_base.add_options()
			("help", "Produce help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator");

	bpo::options_description opt_input ("Input files");
	opt_input.add_options()
			("input,I", bpo::value< string >(), "Reference or target dataset at all variable positions in VCF/BCF format. The GT field is not required.")
			("region", bpo::value< string >(), "Chromosome or region to be splitted");

	bpo::options_description opt_algo ("Parameters");
	opt_algo.add_options()
			("window-size", bpo::value<int>()->default_value(1000000), "Minimal Window size in bp")
			("window-count", bpo::value<int>()->default_value(1000), "Minimal window size in #variants")
			("buffer-size", bpo::value<int>()->default_value(200000), "Minimal buffer size in bp")
			("buffer-count", bpo::value<int>()->default_value(100), "Minimal buffer size in #variants");

	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("output,O", bpo::value< string >(), "Coordinate files for driving GLIMPSE phase runs")
			("log", bpo::value< string >(), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_algo).add(opt_output);
}

void chunker::parse_command_line(vector < string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { cerr << "Error parsing command line arguments: " << string(e.what()) << endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < string > () +"]");

	vrb.title("[GLIMPSE] Split chromosomes into chunks");
	vrb.bullet("Author        : Simone RUBINACCI & Olivier DELANEAU, University of Lausanne");
	vrb.bullet("Contact       : simone.rubinacci@unil.ch & olivier.delaneau@unil.ch");
	vrb.bullet("Version       : " + string(CHUNK_VERSION));
	vrb.bullet("Run date      : " + tac.date());

	if (options.count("help")) { cout << descriptions << endl; exit(0); }
}

void chunker::check_options() {
	if (!options.count("input"))
		vrb.error("You must specify one input file using --input");

	if (!options.count("region"))
		vrb.error("You must specify a region to be split using --region (ideally a chromosome)");

	if (!options.count("output"))
		vrb.error("You must specify a phased output file with --output");

	if (options.count("seed") && options["seed"].as < int > () < 0)
		vrb.error("Random number generator needs a positive seed value");
}

void chunker::verbose_files() {
	vrb.title("Files:");
	vrb.bullet("Input VCF      : [" + options["input"].as < string > () + "]");
	vrb.bullet("Chromosome     : [" + options["region"].as < string > () + "]");
	vrb.bullet("Output file    : [" + options["output"].as < string > () + "]");
	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < string > () + "]");
}

void chunker::verbose_options() {
	vrb.title("Parameters:");
	vrb.bullet("Seed             : " + stb.str(options["seed"].as < int > ()));
	vrb.bullet("Min. Window size : " + stb.str(options["window-size"].as < int > ()) + "bp / " + stb.str(options["window-count"].as < int > ()) + " variants");
	vrb.bullet("Min. Buffer size : " + stb.str(options["buffer-size"].as < int > ()) + "bp / " + stb.str(options["buffer-count"].as < int > ()) + " variants");
}
