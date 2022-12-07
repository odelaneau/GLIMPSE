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
#include <checker/checker_header.h>

void checker::declare_options() {
	bpo::options_description opt_base ("Basic options");
	opt_base.add_options()
			("help", "Produces help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator")
			("threads", bpo::value<int>()->default_value(1), "Number of threads");

	bpo::options_description opt_input ("Input parameters");
	opt_input.add_options()
			("input", bpo::value< std::string >(), "File with four columns listing in order: regions frequencies validation and imputed dataset. For genome-wide concordance, add more lines specifying different chromosomes.")
			("samples", bpo::value< std::string >(), "List of samples to process, one sample ID per line")
			("gt-val", "Uses hard called genotypes rather than phread-scaled likelihoods for the validation dataset, reading them from FORMAT/GT field.")
			("gt-tar", "Uses FORMAT/GT field to determine the best-guess genotype rather than the FORMAT/GP (default). FORMAT/DS are FORMAT/GP fields are still required for calibration and rsquared calculations.");


	bpo::options_description opt_algo ("Other parameters");
	opt_algo.add_options()
			("af-tag", bpo::value< std::string >()->default_value("AF"), "Allele frequency INFO tag to use for binning. By default the allele frequency is estimated from the INFO/AF tag.")
			("use-alt-af", "If specified, the metrics work on the ALT allele frequency (range [0,1]), rather than minor allele frequency (range [0,0.5]).")
			("bins", bpo::value< std::vector < double > >()->multitoken(), "Allele frequency bins used for rsquared computations. By default they should as MAF bins [0-0.5], while they should take the full range [0-1] if --use-ref-alt is used.")
			("ac-bins", bpo::value< std::vector < int > >()->multitoken(), "User-defined allele count bins used for rsquared computations.")
			("allele-counts", "Default allele count bins used for rsquared computations. AN field must be defined in the frequency file.")
			("min-val-gl", bpo::value<double>(), "Minimum genotype likelihood probability P(G|R) in validation data [set to zero to have no filter of if using --gt-validation]")
			("min-val-dp", bpo::value<int>(), "Minimum coverage in validation data. If FORMAT/DP is missing and --minDP > 0, the program exits with an error. [set to zero to have no filter of if using --gt-validation]")
			("min-tar-gp", bpo::value< std::vector < float > >()->multitoken(), "Minimum GP probabilities to be used as a filter. By default it looks at the GP field to specify the filter, but will try to use FORMAT/PL if gt-tar option is specified. Leave empty if no filter is used.")
			("out-r2-per-site", "Output r2 at each site.")
			("out-rej-sites", "Output sites where that cannot be used for the concordance.")
			("out-conc-sites", "Output sites where all target genotypes are concordant with the truth.")
			("out-disc-sites", "Output sites where at least one target genotype is diconcordant with the truth.")
			("groups", bpo::value< std::string >(), "Alternative to frequency bins: group bins are user defined, provided in a file.");

	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("output,O", bpo::value< std::string >(), "Prefix of the output files (extensions are automatically added)")
			("log", bpo::value< std::string >(), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_algo).add(opt_output);
}

void checker::parse_command_line(std::vector < std::string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { std::cerr << "Error parsing command line arguments: " << std::string(e.what()) << std::endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < std::string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < std::string > () +"]");

	vrb.title("[GLIMPSE2] Check concordance of imputed data");
	vrb.bullet("Authors              : Simone RUBINACCI & Olivier DELANEAU, University of Lausanne");
	vrb.bullet("Contact              : simone.rubinacci@unil.ch & olivier.delaneau@unil.ch");
	vrb.bullet("Version       	 : GLIMPSE2_concordance v" + std::string(LIGATE_VERSION) + " / commit = " + std::string(__COMMIT_ID__) + " / release = " + std::string (__COMMIT_DATE__));
	vrb.bullet("Citation	         : BiorXiv, (2022). DOI: https://doi.org/10.1101/2022.11.28.518213");
	vrb.bullet("        	         : Nature Genetics 53, 120–126 (2021). DOI: https://doi.org/10.1038/s41588-020-00756-0");
	vrb.bullet("Run date      	 : " + tac.date());

	if (options.count("help")) { std::cout << descriptions << std::endl; exit(0); }
}

void checker::check_options() {
	if (!options.count("input"))
		vrb.error("You must specify --input");

	if (!options.count("output"))
		vrb.error("You must specify --output");

	if ((options.count("bins")+options.count("allele-counts")+options.count("ac-bins"))+options.count("groups") != 1)
		vrb.error("You must specify only one of --bins or --ac-bins or --allele-counts or --groups");

	if (options["threads"].as < int > () < 1)
		vrb.error("Number of threads is a strictly positive number.");

	//if (options.count("min-tar-gp") &&  options.count("gt-tar"))
	//	vrb.error("Options --gt-tar and --min-tar-gp cannot be used together.");

	if ((options.count("min-val-gl") || options.count("min-val-dp")) &&  options.count("gt-val"))
		vrb.error("Options --gt-val cannot be used when --min-val-gl or min-val-dp are used. If gt-val is used all sites are take in the validation (filtering in the validation dataset has to be done prior to GLIMPSE concordance).");

	if (!options.count("gt-val"))
	{
		if (!options.count("min-val-gl"))
			vrb.error("You must specify --min-val-gl");

		if (!options.count("min-val-dp"))
			vrb.error("You must specify --min-val-dp");
	}

}

void checker::verbose_files() {
	vrb.title("Files are listed in [" + options["input"].as < std::string > () + "]");
	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < std::string > () + "]");
	if (options.count("groups")) {
		vrb.bullet("Groups        : [" + options["groups"].as < std::string > () + "]");
	}
}

void checker::verbose_options() {
	const std::array<std::string,2> no_yes = {"NO","YES"};

	vrb.title("Parameters:");
	vrb.bullet("Output            : " + options["output"].as < std::string > ());
	vrb.bullet("Seed              : " + stb.str(options["seed"].as < int > ()));
	vrb.bullet("#Threads          : " + stb.str(options["threads"].as < int > ()));
	if (!options.count("gt-val"))
	{
		vrb.bullet("Min validation GL : " + stb.str(options["min-val-gl"].as < double > ()));
		vrb.bullet("Min validation DP : " + stb.str(options["min-val-dp"].as < int > ()));
	}
	else
	{
		vrb.bullet("Using all validation GT");
	}


	if (options.count("af-tag")) vrb.bullet("Using INFO/" + options["af-tag"].as < std::string > () + " tag for allele frequency bins.");
	else vrb.bullet("Using INFO/AF tag for binning.");

	if (options.count("bins")) {
		std::vector < double > tmp = options["bins"].as < std::vector < double > > ();
		vrb.bullet("#bins             : " + stb.str(tmp.size()));
	}
	else if (options.count("ac-bins"))
	{
		std::vector < int > tmp = options["ac-bins"].as < std::vector < int > > ();
		vrb.bullet("#AC bins             : " + stb.str(tmp.size()));
	}
	else if (options.count("groups"))
	{
		vrb.bullet("#groups : User defined");
	}
	else
	{
		vrb.bullet("#bins : Fixed bins defined based on allele frequency");
	}

	if (options.count("min-tar-gp")) {
		std::vector < float > tmp = options["min-tar-gp"].as < std::vector < float > > ();
		if  (tmp.size() > 0)
			vrb.bullet("#min-tar-gp       : " + stb.str(tmp.size()));
	}
	vrb.bullet("Output r2 per site: " + no_yes[options.count("out-r2-per-site")]);
	vrb.bullet("Output rej. sites : " + no_yes[options.count("out-rej-sites")]);
	vrb.bullet("Output conc sites : " + no_yes[options.count("out-conc-sites")]);
	vrb.bullet("Output disc sites : " + no_yes[options.count("out-disc-sites")]);
}
