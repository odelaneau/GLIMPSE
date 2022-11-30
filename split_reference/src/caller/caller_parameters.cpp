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
#include <caller/caller_header.h>


void caller::declare_options() {
	bpo::options_description opt_base ("Basic parameters");
	opt_base.add_options()
			("help", "Produce help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator")
			("thread", bpo::value<int>()->default_value(1), "Number of threads");

	bpo::options_description opt_input ("Input parameters");
	opt_input.add_options()
			("input-file", bpo::value < std::string >(), "Text file containing chunks (output of GLIMPSE_chunk)")
			("input-region", bpo::value < std::string >(), "Region to be scanned (buffers included)")
			("output-region", bpo::value < std::string >(), "Phased genomic region to output (same as --input-region but without buffers)")
			("reference,R", bpo::value < std::string >(), "Reference panel of haplotypes in VCF/BCF format")
			("map,M", bpo::value < std::string >(), "Genetic map")
			("sparse-maf", bpo::value<float>()->default_value(0.001f), "Expert setting: sites below this MAF are represented/processed using sparse data structures when using VCF/BCF reference panel. Please do not change if not necessary: performance of the software is highly dependent on this parameter.")
			("keep-monomorphic-ref-sites", "Keeps monomorphic markers in the reference panel, that are removed by default.");

	bpo::options_description opt_output ("Output parameters");
	opt_output.add_options()
			("output,O", bpo::value< std::string >(), "Prefix of the output file [region and extension is automatically added]")
			("log", bpo::value< std::string >(), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_output);
}

void caller::parse_command_line(std::vector < std::string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { std::cerr << "Error parsing command line arguments: " << std::string(e.what()) << std::endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < std::string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < std::string > () +"]");

	vrb.title("[GLIMPSE2] Split reference panel into binary GLIMPSE2 files");
	vrb.bullet("Authors              : Simone RUBINACCI & Olivier DELANEAU, University of Lausanne");
	vrb.bullet("Contact              : simone.rubinacci@unil.ch & olivier.delaneau@unil.ch");
	vrb.bullet("Version       	 : GLIMPSE2_split_reference v" + std::string(LIGATE_VERSION) + " / commit = " + std::string(__COMMIT_ID__) + " / release = " + std::string (__COMMIT_DATE__));
	vrb.bullet("Citation	         : BiorXiv, (2022). DOI: https://doi.org/10.1101/2022.11.28.518213");
	vrb.bullet("        	         : Nature Genetics 53, 120–126 (2021). DOI: https://doi.org/10.1038/s41588-020-00756-0");
	vrb.bullet("Run date      	 : " + tac.date());

	if (options.count("help")) { std::cout << descriptions << std::endl; std::exit(0); }
}

void caller::check_options() {
	if (!options.count("input-file") && !(options.count("input-region") && options.count("output-region")))
			vrb.error("You must specify input files using one of the following options: --input-file or both --input-region / --output-region");

	if (options.count("input-file") && options.count("input-region") > 1)
			vrb.error("You must specify only input parameter between --input-file, --input-region");

	if (!options.count("output"))
		vrb.error("You must specify a output file with --output");

	if (!options.count("reference"))
		vrb.error("You must use a reference panel using --reference");

	if (options.count("seed") && options["seed"].as < int > () < 0)
		vrb.error("Random number generator needs a positive seed value");

	if (options["thread"].as < int > () < 1)
		vrb.error("Number of threads is a strictly positive number.");

	float s_maf = options["sparse-maf"].as < float > ();
	if (s_maf >= 0.5 || s_maf < 0) vrb.error("The sparse MAF parameter should not be set too high or low. Ideally within the range [1%-0.001%] MAF: 0.1% MAF is the recommended setting]");
}

void caller::verbose_files() {
	vrb.title("Files:");

	if (options.count("input-file"))
		vrb.bullet("Input chunks file    : [" + options["input-file"].as < std::string > () + "]");
	else if(options.count("input-region") && options.count("output-region"))
	{
		vrb.bullet("Input region         : [" + options["input-region"].as < std::string > () + "]");
		vrb.bullet("Output region        : [" + options["output-region"].as < std::string > () + "]");
	}
	else vrb.error("Error with input options");

	vrb.bullet("Reference VCF        : [" + options["reference"].as < std::string > () + "]");
	vrb.bullet("Output binary        : [" + options["output"].as < std::string > () + "]");
	if (options.count("map"))
		vrb.bullet("Genetic Map          : [" + options["map"].as < std::string > () + "]");

	if (options.count("log")) vrb.bullet("Output LOG           : [" + options["log"].as < std::string > () + "]");
}

void caller::verbose_options()
{
	std::array<std::string,2> no_yes = {"NO","YES"};

	vrb.title("GLIMPSE_split_reference parameters:");
	vrb.bullet("Sparse MAF           : [" + stb.str(options["sparse-maf"].as < float > () * 100.0f) + "%]");
	vrb.bullet("Keep monom. ref sites: [" + no_yes[options.count("keep-monomorphic-ref-sites")] + "]");
	vrb.bullet("Seed                 : [" + stb.str(options["seed"].as < int > ()) + "]");
	vrb.bullet("#Threads             : [" + stb.str(options["thread"].as < int > ()) + "]");


}
