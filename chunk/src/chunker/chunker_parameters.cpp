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
#include <chunker/chunker_header.h>

void chunker::declare_options() {
	bpo::options_description opt_base ("Basic options");
	opt_base.add_options()
			("help", "Produces help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator")
			("threads", bpo::value<int>()->default_value(1), "Number of threads");

	bpo::options_description opt_input ("Input parameters");
	opt_input.add_options()
			("input,I", bpo::value< std::string >(), "Reference or target dataset at all variable positions in VCF/BCF format. The GT field is not required")
			("region", bpo::value< std::string >(), "Chromosome or region to be split")
			("map,M", bpo::value < std::string >(), "Genetic map")
			("sparse-maf", bpo::value<float>()->default_value(0.001f), "(Expert setting) Rare variant threshold");

	bpo::options_description opt_param ("Parameters");
	opt_param.add_options()
			("window-cm", bpo::value<float>()->default_value(4.0), "Minimal window size in cM")
			("window-mb", bpo::value<float>()->default_value(4.0), "Minimal window size in Mb")
			("window-count", bpo::value<int>()->default_value(30000), "Minimal window size in #variants")
			("buffer-cm", bpo::value<float>()->default_value(0.5), "Minimal buffer size in cM")
			("buffer-mb", bpo::value<float>()->default_value(0.5), "Minimal buffer size in Mb")
			("buffer-count", bpo::value<int>()->default_value(3000), "Minimal buffer size in #variants");

	bpo::options_description opt_algo ("Model parameters");
	opt_algo.add_options()
			("recursive", "Recursive algorithm")
			("sequential", "(Recommended). Sequential algorithm")
			("uniform-number-variants","(Experimental) Uniform the number of variants in the sequential algorithm");

	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("output,O", bpo::value< std::string >(), "File containing the chunks for phasing and imputation")
			("log", bpo::value< std::string >(), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_param).add(opt_algo).add(opt_output);
}

void chunker::parse_command_line(std::vector < std::string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { std::cerr << "Error parsing command line arguments: " << std::string(e.what()) << std::endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < std::string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < std::string > () +"]");

	vrb.title("[GLIMPSE2] Split chromosomes into chunks");
	vrb.bullet("Authors              : Simone RUBINACCI & Olivier DELANEAU, University of Lausanne");
	vrb.bullet("Contact              : simone.rubinacci@unil.ch & olivier.delaneau@unil.ch");
	vrb.bullet("Version       	 : GLIMPSE2_chunk v" + std::string(CHUNK_VERSION) + " / commit = " + std::string(__COMMIT_ID__) + " / release = " + std::string (__COMMIT_DATE__));
	vrb.bullet("Citation	         : BiorXiv, (2022). DOI: https://doi.org/10.1101/2022.11.28.518213");
	vrb.bullet("        	         : Nature Genetics 53, 120–126 (2021). DOI: https://doi.org/10.1038/s41588-020-00756-0");
	vrb.bullet("Run date      	 : " + tac.date());

	if (options.count("help")) { std::cout << descriptions << std::endl; exit(0); }
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

	if (options["threads"].as < int > () < 1)
		vrb.error("Number of threads is a strictly positive number.");

	if (options["window-cm"].as < float > () <= 0)
		vrb.error("Window size in cM must be positive");
	if (options["window-mb"].as < float > () <= 0)
		vrb.error("Window size in Mb must be positive");
	if (options["window-count"].as < int > () <= 0)
		vrb.error("Window size in number of markers must be positive");

	if (options["buffer-cm"].as < float > () <= 0)
		vrb.error("Buffer size in cM must be positive");
	if (options["buffer-mb"].as < float > () <= 0)
		vrb.error("Buffer size in Mb must be positive");
	if (options["buffer-count"].as < int > () <= 0)
		vrb.error("Buffer size in number of markers must be positive");

	float s_maf = options["sparse-maf"].as < float > ();
	if (s_maf >= 0.5 || s_maf < 0) vrb.error("The sparse MAF parameter should not be set too high or low. Ideally within the range [0.01-0.001] 0.001 MAF is the recommended setting]");

	if (options.count("recursive") && options.count("sequential"))
		vrb.error("One of the two algorithms must be selected. Please choose one between recursive and sequential");
	if (!options.count("recursive") && !options.count("sequential"))
		vrb.error("One of the two algorithms must be selected. Please choose one between recursive and sequential");

}

void chunker::verbose_files() {
	vrb.title("Files:");
	vrb.bullet("Input VCF            : [" + options["input"].as < std::string > () + "]");
	vrb.bullet("Region               : [" + options["region"].as < std::string > () + "]");
	if (options.count("map"))
		vrb.bullet("Genetic Map          : [" + options["map"].as < std::string > () + "]");
	vrb.bullet("Output file          : [" + options["output"].as < std::string > () + "]");
	if (options.count("log")) vrb.bullet("Output LOG           : [" + options["log"].as < std::string > () + "]");
}

void chunker::verbose_options() {
	std::string opt_algo = options.count("sequential") ? "Sequential" : "Recursive";
	std::string opt_map = options.count("map") ? "Given by genetic map" : "Constant rate of 1cM/Mb";

	vrb.title("Parameters:");
	vrb.bullet("Sparse MAF           : [" + stb.str(options["sparse-maf"].as < float > ()) + "]");
	vrb.bullet("Algorithm            : [" + opt_algo + "]");
	vrb.bullet("Recombination rates  : [" + opt_map + "]");
	vrb.bullet("Min. Window size     : [" + stb.str(options["window-cm"].as < float > ()) + "cM | " + stb.str(options["window-mb"].as < float > ()) + "Mb | " + stb.str(options["window-count"].as < int > ()) + " variants]");
	vrb.bullet("Min. Buffer size     : [" + stb.str(options["buffer-cm"].as < float > ()) + "cM | " + stb.str(options["buffer-mb"].as < float > ()) + "Mb | " + stb.str(options["buffer-count"].as < int > ()) + " variants]");

	vrb.title("Other parameters");
	vrb.bullet("Seed                 : [" + stb.str(options["seed"].as < int > ()) + "]");
	vrb.bullet("#Threads             : [" + stb.str(options["threads"].as < int > ()) + "]");
}
