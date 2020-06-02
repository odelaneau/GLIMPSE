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
#include <stater/stater_header.h>

void stater::declare_options() {
	bpo::options_description opt_base ("Basic options");
	opt_base.add_options()
			("help", "Produce help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator");

	bpo::options_description opt_input ("Input files");
	opt_input.add_options()
			("input,I", bpo::value< string >(), "Text file listing chromosome IDs (column 1) and input VCF/BCF filenames (column 2) with GP populated.");

	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("output,O", bpo::value< string >(), "Prefix of the output filenames")
			("log", bpo::value< string >(), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_output);
}

void stater::parse_command_line(vector < string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { cerr << "Error parsing command line arguments: " << string(e.what()) << endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < string > () +"]");

	vrb.title("[GLIMPSE] Compute QC metrics for GLIMPSE imputed files");
	vrb.bullet("Author        : Simone RUBINACCI & Olivier DELANEAU, University of Lausanne");
	vrb.bullet("Contact       : simone.rubinacci@unil.ch & olivier.delaneau@unil.ch");
	vrb.bullet("Version       : " + string(STATS_VERSION));
	vrb.bullet("Run date      : " + tac.date());

	if (options.count("help")) { cout << descriptions << endl; exit(0); }
}

void stater::check_options() {
	if (!options.count("input"))
		vrb.error("You must specify one input file using --input");

	if (!options.count("output"))
		vrb.error("You must specify a phased output file with --output");

	if (options.count("seed") && options["seed"].as < int > () < 0)
		vrb.error("Random number generator needs a positive seed value");
}

void stater::verbose_files() {
	vrb.title("Files:");
	vrb.bullet("Input file     : [" + options["input"].as < string > () + "]");
	vrb.bullet("Output file    : [" + options["output"].as < string > () + "]");
	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < string > () + "]");
}

void stater::verbose_options() {
	vrb.title("Parameters:");
	vrb.bullet("Seed             : " + stb.str(options["seed"].as < int > ()));
}
