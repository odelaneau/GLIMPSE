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

#include "../../../versions/versions.h"
#include <caller/caller_header.h>


void caller::declare_options() {
	bpo::options_description opt_base ("Basic parameters");
	opt_base.add_options()
			("help", "Produces help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator")
			("threads", bpo::value<int>()->default_value(1), "Number of threads");

	bpo::options_description opt_input ("Input parameters");
	opt_input.add_options()
			("bam-file", bpo::value < std::string >(), "Input BAM/CRAM file containing low-coverage sequencing reads. Only one of the following options can be declared: --input-gl, --bam-file, --bam-list.")
			("bam-list", bpo::value < std::string >(), "List (.txt file) of input BAM/CRAM files containing low-coverage sequencing reads. One file per line. A second column (space separated) can be used to specify the sample name, otherwise the name of the file is used. Only one of the following options can be declared: --input-gl, --bam-file, --bam-list.")
			("input-gl", bpo::value < std::string >(), "VCF/BCF file containing the genotype likelihoods. Only one of the following options can be declared: --input-gl, --bam-file, --bam-list.")
			("reference,R", bpo::value < std::string >(), "Haplotype reference in VCF/BCF or binary format")
			("map,M", bpo::value < std::string >(), "Genetic map")
			("input-region", bpo::value < std::string >(), "Imputation region with buffers")
			("output-region", bpo::value < std::string >(), "Imputation region without buffers")
			("sparse-maf", bpo::value<float>()->default_value(0.001f), "(Expert setting) Rare variant threshold")
			("samples-file",  bpo::value < std::string >(), "File with sample names and ploidy information. One sample per line with a mandatory second column indicating ploidy (1 or 2). Sample names that are not present are assumed to have ploidy 2 (diploids). If the parameter is omitted, all samples are assumed to be diploid. GLIMPSE does NOT handle the use of sex (M/F) instead of ploidy.")
			("ind-name", bpo::value < std::string >(), "Only used together with --bam-file. Name of the sample to be processed. If not specified the prefix of the BAM/CRAM file (--bam-file) is used.")
			("keep-monomorphic-ref-sites", "(Expert setting) Keeps monomorphic markers in the reference panel (removed by default)")
			("checkpoint-file-in", bpo::value < std::string >(), "File to read in checkpoint from");

	bpo::options_description opt_vcf_input ("VCF/BCF genotype likelihoods input parameters");
	opt_vcf_input.add_options()
			("impute-reference-only-variants", "Only used together with --input-gl. Allows imputation at variants only present in the reference panel. The use of this option is intended only to allow imputation at sporadic missing variants. If the number of missing variants is non-sporadic, please re-run the genotype likelihood computation at all reference variants and avoid using this option, since data from the reads should be used. A warning is thrown if reference-only variants are found.")
			("input-field-gl", "Only used together with --input-gl. Use FORMAT/GL field instead of FORMAT/PL to read genotype likelihoods")
			("use-gl-indels", "(Expert setting) Only used together with --input-gl. Use genotype likelihoods at indels from the VCF/BCF file. By default GLIMPSE assumes flat likelihoods at non-SNP variants, as genotype likelihoods from low-coverage data are often miscalibrated, potentially affecting neighbouring variants.");

	bpo::options_description opt_algo ("Model parameters");
	opt_algo.add_options()
			("burnin", bpo::value<int>()->default_value(5), "(Expert setting) Number of burn-in iterations of the Gibbs sampler")
			("main", bpo::value<int>()->default_value(15), "(Expert setting) Number of main iterations of the Gibbs sampler")
			("ne", bpo::value<int>()->default_value(100000), "(Expert setting) Effective diploid population size modelling recombination frequency")
			("min-gl", bpo::value<float>()->default_value(1e-10f), "(Expert setting) Minimim haploid likelihood")
			("err-imp", bpo::value<float>()->default_value(1e-12f), "(Expert setting) Imputation HMM error rate")
			("err-phase", bpo::value<float>()->default_value(0.0001f), "(Expert setting) Phasing HMM error rate");

	bpo::options_description opt_selection ("Selection parameters");
	opt_selection.add_options()
			("pbwt-depth", bpo::value<int>()->default_value(12), "(Expert setting) Number of neighbors in the sparse PBWT selection step (positive number)")
			("pbwt-modulo-cm", bpo::value<float>()->default_value(0.1f), "(Expert setting) Frequency of PBWT selection in cM (positive number). This parameter is automatically adjusted in case of small imputation regions")
			("Kinit", bpo::value<int>()->default_value(1000), "(Expert setting) Number of states used for initialization (positive number). Can be set to zero only when –state-list is set, to skip the selection for the initialization step")
			("Kpbwt", bpo::value<int>()->default_value(2000), "(Expert setting)  Maximum number of states selected from the sparse PBWT (positive number). Can be set to zero only when –state-list is set, to skip the selection for during the Gibbs iterations")
			("state-list", bpo::value < std::string >(), "(Expert setting) List (.txt file) of haplotypes always present in the conditioning set, independent from state selection. Not affected by other selection parameters. Each row is a target haplotype (two lines per sample in case of diploid individuals) each column is a space separated list of reference haplotypes (in numerical order 0-(2N-1) ). Useful when prior knowledge of relatedness between the reference and target panel is known a priori." );


	bpo::options_description opt_filters ("BAM/CRAM options and filters");
	opt_filters.add_options()
			("call-model", boost::program_options::value< std::string >()->default_value("standard"), "Model to use to call the data. Only the standard model is available at the moment.")
			("call-indels", "(Expert setting) Use the calling model to produce genotype likelihoods at indels. However the likelihoods from low-coverage data can be miscalibrated, therefore by default GLIMPSE does only imputation into the haplotype scaffold (assuming flat likelihoods)")
			("fasta,F", boost::program_options::value< std::string >(), "Faidx-indexed reference sequence file in the appropriate genome build. Necessary for CRAM files")
			("mapq", boost::program_options::value< int >()->default_value(10), "Minimum mapping quality for a read to be considered")
			("baseq", boost::program_options::value< int >()->default_value(10), "Minimum phred-scalde based quality to be considered")
			("max-depth", boost::program_options::value< int >()->default_value(40), "(Expert setting) Max read depth at a site. If the number of reads exceeds this number, the reads at the sites are downsampled (e.g. to avoid artifactual coverage increases)")
			("keep-orphan-reads", "(Expert setting) Keep paired end reads where one of mates is unmapped")
			("ignore-orientation", "(Expert setting) Ignore the orientation of mate pairs")
			("check-proper-pairing", "(Expert setting) Discard reads that are not properly paired")
			("keep-failed-qc", "(Expert setting) Keep reads that fail sequencing QC (as indicated by the sequencer)")
			("keep-duplicates", "(Expert setting) Keep duplicate sequencing reads in the process")
			("illumina13+", "(Expert setting) Use illimina 1.3 encoding for the base quality (for older sequencing machines)");

	bpo::options_description opt_output ("Output parameters");
	opt_output.add_options()
			("output,O", bpo::value< std::string >(), "Phased and imputed haplotypes in VCF/BCF/BGEN format")
			("contigs-fai", bpo::value< std::string >(), "If specified, header contig names and their lengths are copied from the provided fasta index file (.fai). This allows to create imputed whole-genome files as contigs are present and can be easily merged by bcftools")
			("bgen-bits", boost::program_options::value< int >()->default_value(8), "(Expert setting) Only used toghether when the output is in BGEN file format. Specifies the number of bits to be used for the encoding probabilities of the output BGEN file. If the output is in the .vcf[.gz]/.bcf format, this value is ignored. Accepted values: 1-32")
			("bgen-compr", boost::program_options::value< std::string >()->default_value("zstd"), "(Expert setting) Only used toghether when the output is in BGEN file format. Specifies the compression of the output BGEN file. If the output is in the .vcf[.gz]/.bcf format, this value is ignored. Accepted values: [no,zlib,zstd]")
			("log", bpo::value< std::string >(), "Log file")
			("checkpoint-file-out", bpo::value < std::string >(), "File to save checkpoint info in.");

	descriptions.add(opt_base).add(opt_input).add(opt_vcf_input).add(opt_algo).add(opt_selection).add(opt_filters).add(opt_output);
}

void caller::parse_command_line(std::vector < std::string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { std::cerr << "Error parsing command line arguments: " << std::string(e.what()) << std::endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < std::string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < std::string > () +"]");

	vrb.title("[GLIMPSE2] Phase and impute low coverage sequencing data");
	vrb.bullet("Authors              : Simone RUBINACCI & Olivier DELANEAU, University of Lausanne");
	vrb.bullet("Contact              : simone.rubinacci@unil.ch & olivier.delaneau@unil.ch");
	vrb.bullet("Version       	 : GLIMPSE2_phase v" + std::string(PHASE_VERSION) + " / commit = " + std::string(__COMMIT_ID__) + " / release = " + std::string (__COMMIT_DATE__));
	vrb.bullet("Citation	         : BiorXiv, (2022). DOI: https://doi.org/10.1101/2022.11.28.518213");
	vrb.bullet("        	         : Nature Genetics 53, 120–126 (2021). DOI: https://doi.org/10.1038/s41588-020-00756-0");
	vrb.bullet("Run date      	 : " + tac.date());

	if (options.count("help")) { std::cout << descriptions << std::endl; exit(0); }
}

void caller::check_options() {
	if (!options.count("bam-file") && !options.count("bam-list") && !options.count("input-gl"))
			vrb.error("You must specify input files using one of the following options: --bam, --bam-list or --input-gl");

	if (options.count("bam-file")  + options.count("bam-list") + options.count("input-gl") > 1)
			vrb.error("You must specify only input parameter between --bam, --bam-list and --input-gl");

	if (!options.count("output"))
		vrb.error("You must specify a phased output file with --output");

	if (!options.count("reference"))
		vrb.error("You must use a reference panel using --reference");

	if (options.count("seed") && options["seed"].as < int > () < 0)
		vrb.error("Random number generator needs a positive seed value");

	if (options["main"].as < int > () > 15)
		vrb.error("Maximum value for --main is 15. To run more iteration, increase --burn");

	if (options["threads"].as < int > () < 1)
		vrb.error("Number of threads is a strictly positive number.");

	std::string reference_filename = options["reference"].as<std::string>();
	std::string ext0 = stb.get_extension(reference_filename);
	if (ext0 == "bcf" || ext0 == "vcf") input_fmt = InputFormat::BCF;
	else if (ext0 =="gz")
	{
		auto position = reference_filename.find_last_of ( '.' );
		if ( position != std::string::npos )
		{
			std::string str0 (reference_filename.substr(0, position) );
			std::string ext1 = stb.get_extension(str0);
			if (ext1 == "vcf") input_fmt = InputFormat::BCF;
			else vrb.error("Reference panel file format is not supported.");
		}
		else vrb.error("Reference panel file format is not supported.");
	}
	else input_fmt = InputFormat::GLIMPSE; //any extension except vcf/bcf/vcf.gz

	//output
	std::string output_filename = options["output"].as<std::string>();
	ext0 = stb.get_extension(output_filename);
	if (ext0 == "bcf" || ext0 == "vcf")
	{
		if (ext0 == "vcf")
		{
			output_fmt=OutputFormat::VCF;
			output_compr=OutputCompression::NONE;
		}
		else
		{
			output_fmt=OutputFormat::BCF;
			output_compr=OutputCompression::ZLIB;
		}
	}
	else if (ext0 =="gz")
	{
		auto position = output_filename.find_last_of ( '.' );
		if ( position != std::string::npos )
		{
			std::string str0 (output_filename.substr(0, position) );
			std::string ext1 = stb.get_extension(str0);
			if (ext1 == "vcf")
			{
				output_fmt = OutputFormat::VCF;
				output_compr = OutputCompression::ZLIB;
			}
			else vrb.error("Output file format is not supported.");
		}
		else vrb.error("Output file format is not supported.");
	}
	else if (ext0 =="bgen")
	{
		output_fmt = OutputFormat::BGEN;
		std::string ocompr = options["bgen-compr"].as<std::string>();
		if (ocompr == "zstd") output_compr = OutputCompression::ZSTD;
		else if (ocompr == "zlib") output_compr = OutputCompression::ZLIB;
		else if (ocompr == "no") output_compr = OutputCompression::NONE;
		else vrb.error("Output compression must have a value of the following: [no, zlib, zstd].");

		int nbits = options["bgen-bits"].as<int>();
		if (nbits > 0 && nbits <=32) bgen_bits=nbits;
		else vrb.error("Parameter number of bits for BGEN file format encoding is invalid");
	}
	else vrb.error("Output file format is not supported.");

	#ifdef __BGEN__
	#else
		if (output_fmt == OutputFormat::BGEN) vrb.error("GLIMPSE has not been compiled with BGEN support for output file. Please recompile the program with the BGEN_SUPPORT=YES flag or switch to VCF/BCF.");
	#endif

	if (options.count("contigs-fai"))
	{
		std::string fai_fname = options["contigs-fai"].as <std::string>();
		faidx_t *fai = fai_load3(fai_fname.c_str(),fai_fname.c_str(),NULL,FAI_FASTA);
		if ( !fai ) vrb.error("Could not parse " + fai_fname);
		fai_destroy(fai);
	}

	if (input_fmt != InputFormat::GLIMPSE)
	{
		if (!options.count("input-region"))
			vrb.error("You must specify a region to phase using --input-region (this is given by GLIMPSE_chunk)");

		if (!options.count("output-region"))
			vrb.error("You must specify a region to output using --output-region (this is given by GLIMPSE_chunk)");

		if (!options.count("sparse-maf"))
			vrb.error("You must specify a sparse MAF threshold (e.g. 0.001) for reference panel representation.");

		float s_maf = options["sparse-maf"].as < float > ();
		if (s_maf >= 0.5 || s_maf < 0) vrb.error("The sparse MAF parameter should not be set too high or low. Ideally within the range [1%-0.001%] MAF: 0.1% MAF is the recommended setting]");
	}
	else
	{
		if (options.count("input-region"))
			vrb.error("The --input-region parameter is not required when the binary reference panel is used. Value is already provided by the binary reference panel.");

		if (options.count("output-region"))
			vrb.error("The --output-region parameter is not required when the binary reference panel is used. Value is already provided by the binary reference panel.");

		if (options.count("sparse-maf") && !options["sparse-maf"].defaulted())
			vrb.error("The --sparse-maf parameter is not required when the binary reference panel is used. Value is already provided by the binary reference panel.");

		if (options.count("map"))
			vrb.error("The --map parameter is not required when the binary reference panel is used. Value is already provided by the binary reference panel.");
	}

	if (options.count("state-list"))
	{
		input_file state_list(options["state-list"].as<std::string>());
		if (state_list.fail()) vrb.error("Impossible to open file: " + options["state-list"].as<std::string>());

		if (!options.count("Kinit") || options["Kinit"].as < int > () < 0) vrb.error("Option --Kinit is not present or set to a negative value");
		if (!options.count("Kpbwt") || options["Kpbwt"].as < int > () < 0) vrb.error("Option --Kpbwt is not present or set to a negative value");
		if (options["Kpbwt"].as < int > () > 0 && options["pbwt-depth"].as < int > () <= 0) vrb.error("Option --pbwt-depth is not present or set to a non-positive value");
		if (options["Kpbwt"].as < int > () > 0 && options["pbwt-modulo-cm"].as < float > () <= 0) vrb.error("Option --pbwt-modulo-cm is not present or set to a non-positive value");
	}
	else
	{
		if (!options.count("Kinit") || options["Kinit"].as < int > () <= 0) vrb.error("Option --Kinit is not present or set to a non-positive value");
		if (!options.count("Kpbwt") || options["Kpbwt"].as < int > () <= 0 ) vrb.error("Option --Kpbwt is not present or set to a non-positive value");
		if (!options.count("pbwt-depth") || options["pbwt-depth"].as < int > () <= 0) vrb.error("Option --pbwt-depth is not present or set to a non-positive value");
		if (!options.count("pbwt-modulo-cm") || options["pbwt-modulo-cm"].as < float > () <= 0) vrb.error("Option --pbwt-modulo-cm is not present or set to a non-positive value");
	}

	if (options["max-depth"].as <int> () < 10) vrb.error("Max depth has been set too low [< 10].");
}

void caller::verbose_files() {
	vrb.title("Files:");
	std::array<std::string, 3> fmt2string = { {"VCF","BCF","BGEN"} };
	std::array<std::string, 3> compr2string = { {"NO","ZLIB","ZSTD"} };
	std::string n_bits = "";
	if (output_fmt == OutputFormat::BGEN) n_bits = " - " + stb.str(bgen_bits) + " bits";

	if (options.count("bam-list"))
		vrb.bullet("List BAM/CRAM        : [" + options["bam-list"].as < std::string > () + "]");
	else if(options.count("bam-file"))
		vrb.bullet("Input BAM/CRAM       : [" + options["bam-file"].as < std::string > () + "]");
	else
		vrb.bullet("Input VCF/BCF likel. : [" + options["input-gl"].as < std::string > () + "]");

	if (options.count("fasta"))
		vrb.bullet("Reference seq. fasta : [" + options["fasta"].as < std::string > () + "]");

	if (input_fmt == InputFormat::BCF) vrb.bullet("Reference VCF        : [" + options["reference"].as < std::string > () + "]");
	else vrb.bullet("Reference binary     : [" + options["reference"].as < std::string > () + "]");

	if (options.count("map"))
		vrb.bullet("Genetic Map          : [" + options["map"].as < std::string > () + "]");

	vrb.bullet("Output file          : [" + options["output"].as < std::string > () + "]");
	vrb.bullet("Output format        : [" + fmt2string[static_cast<int>(output_fmt)] + " format" + n_bits + " | " + compr2string[static_cast<int>(output_compr)] + " compression]");


	if (options.count("log")) vrb.bullet("Output LOG           : [" + options["log"].as < std::string > () + "]");
}

void caller::verbose_options()
{
	std::array<std::string,2> no_yes = {"NO","YES"};
	std::string opt_input_region, opt_output_region, opt_sparse_maf, opt_map, opt_k_init, opt_k_pbwt, opt_pbwt_depth, opt_pbwt_modulo_cm, opt_state_list, opt_samples_file, opt_keep_mono;

	opt_input_region = input_fmt == InputFormat::BCF ? options["input-region"].as < std::string > () : "Given by binary reference panel";
	opt_output_region = input_fmt == InputFormat::BCF ? options["output-region"].as < std::string > () : "Given by binary reference panel";
	opt_sparse_maf = input_fmt == InputFormat::BCF ? stb.str(H.sparse_maf) : "Given by binary reference panel";
	if (input_fmt == InputFormat::BCF) opt_map = options.count("map") ? "Given by genetic map" : "Constant rate of 1cM/Mb";
	else opt_map = "Given by binary reference panel";

	opt_k_init = options["Kinit"].as < int > () > 0 ? stb.str(options["Kinit"].as < int > ()) : std::string("Init selection not enabled") ;
	opt_k_pbwt = options["Kpbwt"].as < int > () > 0 ? stb.str(options["Kpbwt"].as < int > ()) : std::string("PBWT selection not enabled");
	opt_pbwt_depth = options["Kpbwt"].as < int > () > 0 ? stb.str(options["pbwt-depth"].as < int > ()) : std::string("PBWT selection not enabled");
	opt_pbwt_modulo_cm = options["Kpbwt"].as < int > () > 0 ? stb.str(options["pbwt-modulo-cm"].as < float > ()) : std::string("PBWT selection not enabled");
	opt_state_list = options.count("state-list") ? stb.str(options["state-list"].as < std::string > ()) : std::string("No list provided");
	opt_samples_file = options.count("samples-file") ? "Given by samples file" : "Only diploid samples in region";
	const float err_imp = std::clamp(options["err-imp"].as < float > (), 1e-12f, 1e-3f);
	const float err_phase = options["err-phase"].as < float > ();

	vrb.title("GLIMPSE_phase parameters:");
	if (options.count("input-gl") && options.count("input-field-gl")) vrb.bullet("Input                : [FORMAT/GL field for genotype likelihoods]");
	else if (options.count("input-gl")) vrb.bullet("Input                : [FORMAT/PL field for genotype likelihoods]");

	vrb.bullet("Imputation model     : [Reference panel imputation]");
	vrb.bullet("Input region         : [" + opt_input_region + "]");
	vrb.bullet("Output region        : [" + opt_output_region + "]");
	vrb.bullet("Sparse MAF           : [" + opt_sparse_maf + "]");
	vrb.bullet("Recombination rates  : [" + opt_map + "]");
	vrb.bullet("Ploidy               : [" + opt_samples_file + "]");
	vrb.bullet("Keep monom. ref sites: [" + no_yes[options.count("keep-monomorphic-ref-sites")] + "]");

	vrb.title("Model parameters:");
	vrb.bullet("#Burnin iterations   : [" + stb.str(options["burnin"].as < int > ()) + "]");
	vrb.bullet("#Main iterations     : [" + stb.str(options["main"].as < int > ()) + "]");
	vrb.bullet("Ne [eff. pop. size]  : [" + stb.str(options["ne"].as < int > ()) + "]");
	vrb.bullet("Phase error rate     : [" + stb.str(err_phase) + "]");
	vrb.bullet("Imputation error rate: [" + stb.str(err_imp) + "]");
	vrb.bullet("Min value for hap GLs: [" + stb.str(options["min-gl"].as < float > ()) + "]");

	vrb.title("Selection parameters:");
	vrb.bullet("K init               : [" + opt_k_init + "]");
	vrb.bullet("K pbwt               : [" + opt_k_pbwt + "]");
	vrb.bullet("PBWT depth           : [" + opt_pbwt_depth + "]");
	vrb.bullet("PBWT modulo (cM)     : [" + opt_pbwt_modulo_cm + "]");
	vrb.bullet("State list           : [" + opt_state_list + "]");

	if (!options.count("input-gl"))
	{
		vrb.title("Genotype calling:");
		vrb.bullet("Calling model        : [" + options["call-model"].as <std::string> () + "]");
		vrb.bullet("Indels model         : [" + stb.str(options.count("call-indels") ? "Perform calling]" : "Haplotype scaffold]"));

		vrb.title("BAM/CRAM filters and options:");
		vrb.bullet("Min mapping quality  : [" + stb.str(options["mapq"].as <int> ()) + "]");
		vrb.bullet("Min base quality     : [" + stb.str(options["baseq"].as <int> ()) + "]");
		vrb.bullet("Max depth            : [" + stb.str(options["max-depth"].as <int> ()) + "]");

		vrb.bullet("Keep failed QC       : [" + no_yes[options.count("keep-failed-qc")] + "]");
		vrb.bullet("Keep orphan reads    : [" + no_yes[options.count("keep-orphan-read")] + "]");
		vrb.bullet("Keep duplicate reads : [" + no_yes[options.count("keep-duplicates")] + "]");
		vrb.bullet("Keep supp. alignment : [" + no_yes[options.count("keep-supp")] + "]");

		vrb.bullet("Check pairing        : [" + no_yes[options.count("check-proper-pairing")] + "]");
		vrb.bullet("Ignore orientation   : [" + no_yes[options.count("ignore-orientation")] + "]");
		vrb.bullet("Illumina-1.3+        : [" + no_yes[options.count("illumina13")] + "]");
	}
	else
	{
		vrb.bullet("use-gl-indels         : [" + stb.str(options.count("call-indels") ? "Use PLs]" : "Haplotype scaffold]"));
	}

	vrb.title("Other parameters");
	vrb.bullet("Seed                 : [" + stb.str(options["seed"].as < int > ()) + "]");
	vrb.bullet("#Threads             : [" + stb.str(options["threads"].as < int > ()) + "]");

}
