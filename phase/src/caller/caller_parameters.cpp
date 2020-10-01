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
#include <caller/caller_header.h>


void caller::declare_options() {
	bpo::options_description opt_base ("Basic parameters");
	opt_base.add_options()
			("help", "Produce help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator")
			("thread", bpo::value<int>()->default_value(1), "Number of threads");

	bpo::options_description opt_input ("Input parameters");
	opt_input.add_options()
			("input,I", bpo::value < string >(), "Genotypes to be phased in VCF/BCF format")
			("input-region", bpo::value < string >(), "Whole genomic region to be phased (including left/right buffers)")
			("reference,R", bpo::value < string >(), "Reference panel of haplotypes in VCF/BCF format")
			("map,M", bpo::value < string >(), "Genetic map")
			("samples-file",  bpo::value < string >(), "File with sample names and ploidy information. One sample per line with a mandatory second column indicating ploidy (1 or 2). Sample names that are not present are assumed to have ploidy 2 (diploids). If the parameter is omitted, all samples are assumed to be diploid. GLIMPSE does NOT handle the use of sex (M/F) instead of ploidy.")
			("impute-reference-only-variants", "Allows imputation at variants only present in the reference panel. The use of this option is intended only to allow imputation at sporadic missing variants. If the number of missing variants is non-sporadic, please re-run the genotype likelihood computation at all reference variants and avoid using this option, since data from the reads should be used. A warning is thrown if reference-only variants are found.");


	bpo::options_description opt_algo ("Model parameters");
	opt_algo.add_options()
			("burnin", bpo::value<int>()->default_value(10), "Number of Burn-in iterations")
			("main", bpo::value<int>()->default_value(10), "Each main iterations contributes to output genotypes. Haplotypes sampled for the last (max 15) iterations are stored in the HS field.")
			("pbwt-depth", bpo::value<int>()->default_value(2), "Number of neighbors to store")
			("pbwt-modulo", bpo::value<int>()->default_value(8), "Frequency of PBWT storage")
			("init-states", bpo::value<int>()->default_value(1000), "Number of states used for initialization")
			("init-pool", bpo::value< string >(), "Pool of samples from which initializing haplotypes should be chosen")
			("ne", bpo::value<float>()->default_value(20000.0), "Effective diploid population size");

	bpo::options_description opt_output ("Output parameters");
	opt_output.add_options()
			("output,O", bpo::value< string >(), "Phased haplotypes in VCF/BCF format")
			("output-region", bpo::value < string >(), "Phased genomic region to output")
			("log", bpo::value< string >(), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_algo).add(opt_output);
}

void caller::parse_command_line(vector < string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { cerr << "Error parsing command line arguments: " << string(e.what()) << endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < string > () +"]");

	vrb.title("[GLIMPSE] Phase and impute low coverage sequencing data");
	vrb.bullet("Author        : Simone RUBINACCI & Olivier DELANEAU, University of Lausanne");
	vrb.bullet("Contact       : simone.rubinacci@unil.ch & olivier.delaneau@unil.ch");
	vrb.bullet("Version       : " + string(PHASE_VERSION));
	vrb.bullet("Run date      : " + tac.date());

	if (options.count("help")) { cout << descriptions << endl; exit(0); }
}

void caller::check_options() {
	if (!options.count("input"))
		vrb.error("You must specify one input file using --input");

	if (!options.count("input-region"))
		vrb.error("You must specify a region to phase using --input-region (this is given by GLIMPSE_chunk)");

	if (!options.count("output-region"))
		vrb.error("You must specify a region to output using --output-region (this is given by GLIMPSE_chunk)");

	if (!options.count("output"))
		vrb.error("You must specify a phased output file with --output");

	if (!options.count("reference"))
		vrb.error("You must use a reference panel using --reference");

	if (options.count("seed") && options["seed"].as < int > () < 0)
		vrb.error("Random number generator needs a positive seed value");

	if (options["main"].as < int > () > 15)
		vrb.error("Maximum value for --main is 15, to run more iteration, increase --burn");
}

void caller::verbose_files() {
	vrb.title("Files:");
	vrb.bullet("Input VCF     : [" + options["input"].as < string > () + "]");
	vrb.bullet("Reference VCF : [" + options["reference"].as < string > () + "]");
	if (options.count("map")) vrb.bullet("Genetic Map   : [" + options["map"].as < string > () + "]");
	vrb.bullet("Output VCF    : [" + options["output"].as < string > () + "]");
	vrb.bullet("Input region  : [" + options["input-region"].as < string > () + "]");
	vrb.bullet("Output region : [" + options["output-region"].as < string > () + "]");
	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < string > () + "]");
}

void caller::verbose_options() {
	vrb.title("Parameters:");
	vrb.bullet("Seed       : " + stb.str(options["seed"].as < int > ()));
	vrb.bullet("#Threads   : " + stb.str(options["thread"].as < int > ()));
	vrb.bullet("#Burnin    : " + stb.str(options["burnin"].as < int > ()));
	vrb.bullet("#Main      : " + stb.str(options["main"].as < int > ()));
	vrb.bullet("PBWT depth : " + stb.str(options["pbwt-depth"].as < int > ()));
	vrb.bullet("PBWT modulo: " + stb.str(options["pbwt-modulo"].as < int > ()));
	if (options.count("map")) vrb.bullet("HMM     : Recombination rates given by genetic map");
	else vrb.bullet("HMM     : Constant recombination rate of 1cM per Mb");
	if (options.count("samples-file")) vrb.bullet("Ploidy     : given by samples file");
	else vrb.bullet("Ploidy     : All samples are diploids in this region");
	vrb.bullet("Init K     : " + stb.str(options["init-states"].as < int > ()));
	if (options.count("impute-reference-variants")) vrb.bullet("Imputation at reference-only variants is performed");
}
