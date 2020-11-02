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
#include <pir/pir_header.h>

void pir::declare_options() {
	bpo::options_description opt_base ("Basic options");
	opt_base.add_options()
			("help", "Produce help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator");

	bpo::options_description opt_input ("Input files");
	opt_input.add_options()
			("bam", bpo::value < string >(), "BAM file containing sequencing reads")
			("input", bpo::value < string >(), "VCF/BCF file variant sites at which extraction should happen")
			("region", bpo::value < string >(), "Genomic region to be processed");

	bpo::options_description opt_algo ("Parameters");
	opt_algo.add_options()
			("min-map-qual,q", boost::program_options::value< int >()->default_value(10), "Minimum mapping quality for a read to be considered. Set this to only include uniquely mapped reads. (REQUIRED)")
			("min-base-qual,Q", boost::program_options::value< int >()->default_value(10), "Minimum phred quality for a base to be considered.");

	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("output,O", bpo::value< string >(), "Output phase informative reads in gzipped TXT format")
			("log", bpo::value< string >(), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_algo).add(opt_output);
}

void pir::parse_command_line(vector < string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { cerr << "Error parsing command line arguments: " << string(e.what()) << endl; exit(0); }


	vrb.title("[GLIMPSE] Extract phase informative reads in a BAM file given variant sites coordinates");
	vrb.bullet("Author        : Simone RUBINACCI & Olivier DELANEAU, University of Lausanne");
	vrb.bullet("Contact       : simone.rubinacci@unil.ch & olivier.delaneau@unil.ch");
	//vrb.bullet("Version       : " + string(SAMPLE_VERSION));
	vrb.bullet("Run date      : " + tac.date());

	if (options.count("help")) { cout << descriptions << endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < string > () +"]");
}

void pir::check_options() {
	if (!options.count("input"))
		vrb.error("You must specify the an input VCF/BCF using --input");

	if (!options.count("bam"))
		vrb.error("You must specify the an input BAM using --bam");

	if (!options.count("region"))
		vrb.error("You must specify a genomic region using --region");

	if (!options.count("output"))
		vrb.error("You must specify an output file with --output");

	if (options.count("seed") && options["seed"].as < int > () < 0)
		vrb.error("Random number generator needs a positive seed value");
}

void pir::verbose_files() {
	vrb.title("Files:");
	vrb.bullet("Input VCF     : [" + options["input"].as < string > () + "]");
	vrb.bullet("Input BAM     : [" + options["bam"].as < string > () + "]");
	vrb.bullet("Input region  : [" + options["region"].as < string > () + "]");
	vrb.bullet("Output TXT    : [" + options["output"].as < string > () + "]");
	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < string > () + "]");
}

void pir::verbose_options() {
	vrb.title("Parameters:");
	vrb.bullet("Seed       : " + stb.str(options["seed"].as < int > ()));
	vrb.bullet("Min MAPQ   : " + stb.str(options["min-map-qual"].as < int > ()));
	vrb.bullet("Min BASEQ  : " + stb.str(options["min-base-qual"].as < int > ()));
}
