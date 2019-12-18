////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2018 Olivier Delaneau, University of Lausanne
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////

#include <caller/caller_header.h>

void caller::declare_options() {
	bpo::options_description opt_base ("Basic options");
	opt_base.add_options()
			("help", "Produce help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator")
			("thread", bpo::value<int>()->default_value(1), "Number of threads");

	bpo::options_description opt_input ("Input files");
	opt_input.add_options()
			("input,I", bpo::value < std::string >(), "Genotypes to be phased in VCF/BCF format")
			("reference,H", bpo::value < std::string >(), "Reference panel of haplotypes in VCF/BCF format")
			("validation,V", bpo::value < std::string >(), "Validation data in VCF/BCF format")
			("map,M", bpo::value < std::string >(), "Genetic map")
			("region,R", bpo::value < std::string >(), "Target region")
			("buffer", bpo::value < int >()->default_value(200), "Size of the buffer on each side in kilobases (200kb by default)");

	bpo::options_description opt_algo ("Parameters");
	opt_algo.add_options()
			("burnin", bpo::value<int>()->default_value(10), "Burn-in passes")
			("main", bpo::value<int>()->default_value(10), "Main passes")
			("refonly-select", "Selection using only the reference panel")
			("pbwt-depth", bpo::value<int>()->default_value(2), "Number of neighbors to store")
			("pbwt-modulo", bpo::value<int>()->default_value(8), "Number of neighbors to store")
			("init-states", bpo::value<int>()->default_value(1000), "Number of neighbors to store")
			("phasing-switch", "Phasing using switch likelihoods")
			("phasing-flipandswitch", "Phasing using flip and switch likelihoods")
			//("maf-filter", bpo::value<float>()->default_value(0.0f)->implicit_value(0.01f), "MAF filter of variants in the reference panel (requires INFO/AF field)"); for boost >= 1.65
			("maf-filter", bpo::value<float>()->default_value(0.0f), "MAF filter of variants in the reference panel (requires INFO/AF field)");


	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("output,O", bpo::value< std::string >(), "Phased haplotypes in VCF/BCF format")
			("output-HS", "Output haplotype estimates of the last (max 16) main iterations. Output is an integer per genotype stored in INFO/HS field.")
			("log", bpo::value< std::string >(), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_algo).add(opt_output);
}

void caller::parse_command_line(std::vector < std::string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { std::cerr << "Error parsing command line arguments: " << std::string(e.what()) << std::endl; exit(0); }

	if (options.count("help")) { std::cout << descriptions << std::endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < std::string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < std::string > () +"]");

	vrb.title("LCC GWAS");
	vrb.bullet("Author        : Simone RUBINACCI and Olivier DELANEAU, University of Lausanne");
	vrb.bullet("Contact       : olivier.delaneau@gmail.com");
	vrb.bullet("Version       : 1.0.0");
	vrb.bullet("Run date      : " + tac.date());
}

void caller::check_options() {
	if (!options.count("input"))
		vrb.error("You must specify one input file using --input");

	if (!options.count("region"))
		vrb.error("You must specify a region or chromosome to phase using --region");

	if (!options.count("output"))
		vrb.error("You must specify a phased output file with --output");

	if (!options.count("reference"))
		vrb.error("You must use a reference panel using --reference");

	if (options.count("seed") && options["seed"].as < int > () < 0)
		vrb.error("Random number generator needs a positive seed value");

	if (options.count("maf-filter") && (options["maf-filter"].as < float > () > 0.5 || options["maf-filter"].as < float > () < 0.0))
		vrb.error("maf-filter parameter needs a value bounded between 0 and 0.5");



}

void caller::verbose_files() {
	vrb.title("Files:");
	vrb.bullet("Input VCF     : [" + options["input"].as < std::string > () + "]");
	vrb.bullet("Reference VCF : [" + options["reference"].as < std::string > () + "]");
	if (options.count("map")) vrb.bullet("Genetic Map   : [" + options["map"].as < std::string > () + "]");
	vrb.bullet("Output VCF    : [" + options["output"].as < std::string > () + "]");
	vrb.bullet("Genomic region: [" + options["region"].as < std::string > () + "]");
	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < std::string > () + "]");
}

void caller::verbose_options() {
	vrb.title("Parameters:");
	vrb.bullet("Seed       : " + stb.str(options["seed"].as < int > ()));
	vrb.bullet("#Threads   : " + stb.str(options["thread"].as < int > ()));
	vrb.bullet("#Burnin    : " + stb.str(options["burnin"].as < int > ()));
	vrb.bullet("#Main      : " + stb.str(options["main"].as < int > ()));
	vrb.bullet("PBWT depth : " + stb.str(options["pbwt-depth"].as < int > ()));
	vrb.bullet("PBWT modulo: " + stb.str(options["pbwt-modulo"].as < int > ()));
	vrb.bullet("Init K     : " + stb.str(options["init-states"].as < int > ()));
	if (options["maf-filter"].as < float > () > 0.0)
		vrb.bullet("maf-filter : " + stb.str(options["maf-filter"].as < float > ()));
	else
		vrb.bullet("maf-filter : 0.0 (using all sites)");
}
