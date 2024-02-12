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

#include "../../versions/versions.h"
#include <ligater/ligater_header.h>

void ligater::declare_options() {
	bpo::options_description opt_base ("Basic options");
	opt_base.add_options()
			("help", "Produces help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator")
			("threads", bpo::value<int>()->default_value(1), "Number of threads");

	bpo::options_description opt_input ("Input files");
	opt_input.add_options()
			("input", bpo::value < std::string >(), "Text file containing all VCF/BCF to ligate, one file per line");

	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("output,O", bpo::value< std::string >(), "Output ligated (phased) file in VCF/BCF format")
			("log", bpo::value< std::string >(), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_output);
}

void ligater::parse_command_line(std::vector < std::string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { std::cerr << "Error parsing command line arguments: " << std::string(e.what()) << std::endl; exit(0); }

	vrb.title("[GLIMPSE2] Ligate multiple output files into chromosome-wide files");
	vrb.bullet("Authors              : Simone RUBINACCI & Olivier DELANEAU, University of Lausanne");
	vrb.bullet("Contact              : simone.rubinacci@unil.ch & olivier.delaneau@unil.ch");
	vrb.bullet("Version       	 : GLIMPSE2_ligate v" + std::string(LIGATE_VERSION) + " / commit = " + std::string(__COMMIT_ID__) + " / release = " + std::string (__COMMIT_DATE__));
	vrb.bullet("Citation	         : BiorXiv, (2022). DOI: https://doi.org/10.1101/2022.11.28.518213");
	vrb.bullet("        	         : Nature Genetics 53, 120–126 (2021). DOI: https://doi.org/10.1038/s41588-020-00756-0");
	vrb.bullet("Run date      	 : " + tac.date());

	if (options.count("help")) { std::cout << descriptions << std::endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < std::string > ()))
	vrb.error("Impossible to create log file [" + options["log"].as < std::string > () +"]");
}

void ligater::check_options() {
	if (!options.count("input"))
		vrb.error("You must specify the list of VCF/BCF to ligate using --input");

	if (!options.count("output"))
		vrb.error("You must specify an output VCF/BCF file with --output");

	if (options.count("seed") && options["seed"].as < int > () < 0)
		vrb.error("Random number generator needs a positive seed value");

	if (options["threads"].as < int > () < 1)
		vrb.error("Number of threads is a strictly positive number.");
}

void ligater::verbose_files() {
	std::array<std::string,2> no_yes = {"NO","YES"};

	vrb.title("Files:");
	vrb.bullet("Input LIST     : [" + options["input"].as < std::string > () + "]");
	vrb.bullet("Output VCF     : [" + options["output"].as < std::string > () + "]");
	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < std::string > () + "]");
}

void ligater::verbose_options() {
	vrb.title("Parameters:");
	vrb.bullet("Seed           : [" + stb.str(options["seed"].as < int > ()) + "]");
	vrb.bullet("#Threads       : [" + stb.str(options["threads"].as < int > ()) + "]");
}
