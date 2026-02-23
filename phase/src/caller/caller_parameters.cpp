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
			("seed", bpo::value<uint32_t>(&A.mSeed)->default_value(42), "Seed of the random number generator")
			("threads", bpo::value<uint32_t>(&A.mNumThreads)->default_value(1), "Number of threads");

	bpo::options_description opt_input ("Input parameters");
	opt_input.add_options()
			("bam-file", bpo::value < std::string >(&A.mBamFileFilename), "Input BAM/CRAM file containing low-coverage sequencing reads. Only one of the following options can be declared: --input-gl, --bam-file, --bam-list.")
			("bam-list", bpo::value < std::string >(&A.mBamListFilename), "List (.txt file) of input BAM/CRAM files containing low-coverage sequencing reads. One file per line. A second column (space separated) can be used to specify the sample name, otherwise the name of the file is used. Only one of the following options can be declared: --input-gl, --bam-file, --bam-list.")
			("input-gl", bpo::value < std::string >(&A.mInputGLFilename), "VCF/BCF file containing the genotype likelihoods. Only one of the following options can be declared: --input-gl, --bam-file, --bam-list.")
			("reference,R", bpo::value < std::string >(&A.mRefHapFilename), "Haplotype reference in VCF/BCF or binary format")
			("map,M", bpo::value < std::string >(&A.mMapFilename), "Genetic map")
			("input-region", bpo::value < std::string >(&A.mInputRegion), "Imputation region with buffers")
			("output-region", bpo::value < std::string >(&A.mOutputRegion), "Imputation region without buffers")
			("sparse-maf", bpo::value<float>(&A.mSparseMaf)->default_value(0.001f), "(Expert setting) Rare variant threshold")
			("samples-file",  bpo::value < std::string >(&A.mSamplesFilename), "File with sample names and ploidy information. One sample per line with a mandatory second column indicating ploidy (1 or 2). Sample names that are not present are assumed to have ploidy 2 (diploids). If the parameter is omitted, all samples are assumed to be diploid. GLIMPSE does NOT handle the use of sex (M/F) instead of ploidy.")
			("ind-name", bpo::value < std::string >(&A.mIndName), "Only used together with --bam-file. Name of the sample to be processed. If not specified the prefix of the BAM/CRAM file (--bam-file) is used.")
			("keep-monomorphic-ref-sites", "(Expert setting) Keeps monomorphic markers in the reference panel (removed by default)");


	bpo::options_description opt_vcf_input ("VCF/BCF genotype-likelihoods specific parameters");
	opt_vcf_input.add_options()
			("impute-reference-only-variants", "Only used together with --input-gl. Allows imputation at variants only present in the reference panel. The use of this option is intended only to allow imputation at sporadic missing variants. If the number of missing variants is non-sporadic, please re-run the genotype likelihood computation at all reference variants and avoid using this option, since data from the reads should be used. A warning is thrown if reference-only variants are found.")
			("input-field-gl", "Only used together with --input-gl. Use FORMAT/GL field instead of FORMAT/PL to read genotype likelihoods")
			("use-gl-indels", "(Expert setting) Only used together with --input-gl. Use genotype likelihoods at indels from the VCF/BCF file. By default GLIMPSE assumes flat likelihoods at non-SNP variants, as genotype likelihoods from low-coverage data are often miscalibrated, potentially affecting neighbouring variants.");

	bpo::options_description opt_algo ("Model parameters");
	opt_algo.add_options()
			("burnin", bpo::value<uint32_t>(&A.mBurnIn)->default_value(5), "(Expert setting) Number of burn-in iterations of the Gibbs sampler")
			("main", bpo::value<uint32_t>(&A.mMain)->default_value(15), "(Expert setting) Number of main iterations of the Gibbs sampler")
			("ne", bpo::value<uint32_t>(&A.mNe)->default_value(100000), "(Expert setting) Effective diploid population size modelling recombination frequency")
			("min-gl", bpo::value<float>(&A.mMinGL)->default_value(1e-10f), "(Expert setting) Minimim haploid likelihood")
			("err-imp", bpo::value<float>(&A.mErrImp)->default_value(1e-12f), "(Expert setting) Imputation HMM error rate")
			("err-phase", bpo::value<float>(&A.mErrPhase)->default_value(0.0001f), "(Expert setting) Phasing HMM error rate");

	bpo::options_description opt_selection ("Selection parameters");
	opt_selection.add_options()
			("max-pbwt-depth", boost::program_options::value<uint32_t>(&A.mMaxPbwtDepth)->default_value(16), "(Expert setting) Max depth of PBWT indexes to condition on")
			("min-pbwt-depth", boost::program_options::value<uint32_t>(&A.mMinPbwtDepth)->default_value(4), "(Expert setting) Min depth of PBWT indexes to condition on")
			("pbwt-modulo-cm", boost::program_options::value<float>(&A.mPbwtCM)->default_value(0.05f), "(Expert setting) Frequency of PBWT selection in cM (positive number). This parameter is automatically adjusted in case of small imputation regions")
			("Kpbwt", boost::program_options::value<uint32_t>(&A.Kpbwt)->default_value(1500), "(Expert setting)  Maximum number of states selected from the sparse PBWT (positive number). Can be set to zero only when –state-list is set, to skip the selection for during the Gibbs iterations")
			("Kinit", bpo::value<uint32_t>(&A.Kinit)->default_value(1000), "(Expert setting) Number of states used for initialization (positive number). Can be set to zero only when –state-list is set, to skip the selection for the initialization step")
			("state-list", bpo::value < std::string >(&A.mStateListFilename), "(Expert setting) List (.txt file) of haplotypes always present in the conditioning set, independent from state selection. Not affected by other selection parameters. Each row is a target haplotype (two lines per sample in case of diploid individuals) each column is a space separated list of reference haplotypes (in numerical order 0-(2N-1) ). Useful when prior knowledge of relatedness between the reference and target panel is known a priori." );

	bpo::options_description opt_filters ("BAM/CRAM options and filters");
	opt_filters.add_options()
			("call-model", boost::program_options::value< std::string >(&A.mCallModel)->default_value("standard"), "Model to use to call the data. Only the standard model is available at the moment.")
			("call-indels", "(Expert setting) Use the calling model to produce genotype likelihoods at indels. However the likelihoods from low-coverage data can be miscalibrated, therefore by default GLIMPSE does only imputation into the haplotype scaffold (assuming flat likelihoods)")
			("fasta,F", boost::program_options::value< std::string >(&A.mFastaFilename), "Faidx-indexed reference sequence file in the appropriate genome build. Necessary for CRAM files")
			("mapq", boost::program_options::value< uint32_t >(&A.mMapQ)->default_value(10), "Minimum mapping quality for a read to be considered")
			("baseq", boost::program_options::value< uint32_t >(&A.mBaseQ)->default_value(10), "Minimum phred-scalde based quality to be considered")
			("max-depth", boost::program_options::value< uint32_t >(&A.mMaxDepth)->default_value(40), "(Expert setting) Max read depth at a site. If the number of reads exceeds this number, the reads at the sites are downsampled (e.g. to avoid artifactual coverage increases)")
			("keep-orphan-reads", "(Expert setting) Keep paired end reads where one of mates is unmapped")
			("ignore-orientation", "(Expert setting) Ignore the orientation of mate pairs")
			("check-proper-pairing", "(Expert setting) Discard reads that are not properly paired")
			("keep-failed-qc", "(Expert setting) Keep reads that fail sequencing QC (as indicated by the sequencer)")
			("keep-duplicates", "(Expert setting) Keep duplicate sequencing reads in the process")
			("illumina13+", "(Expert setting) Use illimina 1.3 encoding for the base quality (for older sequencing machines)");

	bpo::options_description opt_output ("Output parameters");
	opt_output.add_options()
			("output,O", bpo::value< std::string >(&A.mOutputFilename), "Phased and imputed haplotypes in VCF/BCF/BGEN format")
			("contigs-fai", bpo::value< std::string >(&A.mContigFaiFilename), "If specified, header contig names and their lengths are copied from the provided fasta index file (.fai) instead of being taken from the reference panel (default behavior). This allows to create files with all the contigs in the header (in the case the contigs in the reference panel are limited to a single chromosome) and therefore quickly merge chromosome-level files with bcftools merge --naive")
			("no-out-index", "Skip the calculation of the csi index")
			("no-out-gp-field", "Output FORMAT/GP field (genotype posterior probabilities) if output is in VCF/BCF format.")
			("no-out-ds-field", "Output FORMAT/DS field (genotype dosage) if output is in VCF/BCF format.")
			("out-ap-field", "Output FORMAT/AP field (ALT haplotype probabilities) if output is in VCF/BCF format or outputs phased BGEN file.")
			("bgen-bits", boost::program_options::value< unsigned int >(&A.mBgenNbits)->default_value(8), "(Expert setting) Only used toghether when the output is in BGEN file format. Specifies the number of bits to be used for the encoding probabilities of the output BGEN file. If the output is in the .vcf[.gz]/.bcf format, this value is ignored. Accepted values: 1-32")
			("bgen-compr", boost::program_options::value< std::string >(&A.mBgenCompression)->default_value("zstd"), "(Expert setting) Only used toghether when the output is in BGEN file format. Specifies the compression of the output BGEN file. If the output is in the .vcf[.gz]/.bcf format, this value is ignored. Accepted values: [no,zlib,zstd]")
			("log", bpo::value< std::string >(&A.mLogFilename), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_vcf_input).add(opt_algo).add(opt_selection).add(opt_filters).add(opt_output);
}

void caller::parse_command_line(std::vector < std::string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { std::cerr << "Error parsing command line arguments: " << std::string(e.what()) << std::endl; exit(0); }

	if (!A.mLogFilename.empty() && !vrb.open_log(A.mLogFilename))
		vrb.error("Impossible to create log file [" + A.mLogFilename +"]");

	vrb.title("[GLIMPSE2] Phase and impute low coverage sequencing data");
	vrb.bullet("Authors              : Simone RUBINACCI & Olivier DELANEAU");
	vrb.bullet("Contact              : srubinac@broadinstitute.org & olivier.delaneau@gmail.com");
	vrb.bullet("Version       	 : GLIMPSE2_phase v" + std::string(PHASE_VERSION) + " / commit = " + std::string(__COMMIT_ID__) + " / release = " + std::string (__COMMIT_DATE__));
	vrb.bullet("Citation	         : BiorXiv, (2022). DOI: https://doi.org/10.1101/2022.11.28.518213");
	vrb.bullet("        	         : Nature Genetics 53, 120–126 (2021). DOI: https://doi.org/10.1038/s41588-020-00756-0");
	vrb.bullet("Run date      	 : " + tac.date());

	if (options.count("help")) { std::cout << descriptions << std::endl; exit(0); }

	//set bool options
	A.mKeepMonomorphicRefSites = options.count("keep-monomorphic-ref-sites");
	A.mImputeReferenceOnlyVariants = options.count("impute-reference-only-variants");
	A.mInputFieldGL = options.count("input-field-gl");
	A.mUseGLIndels = options.count("use-gl-indels");

	A.mOutputIndex = !options.count("no-out-index");
	A.mCallIndels = options.count("call-indels");
	A.mKeepFailedQC = options.count("keep-failed-qc");
	A.mKeepOrphanReads = options.count("keep-orphan-reads");
	A.mKeepDuplicates = options.count("keep-duplicates");
	A.mCheckProperPairing = options.count("check-proper-pairing");
	A.mIgnoreOrientation = options.count("ignore-orientation");
	A.mIllumina13 = options.count("illumina13");
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

	if (A.mSeed < 0)
		vrb.error("Random number generator needs a positive seed value");

	if (A.mMain > 15)
		vrb.error("Maximum value for --main is 15. To run more iteration, increase --burn");

	if (A.mNumThreads < 1)
		vrb.error("Number of threads is a strictly positive number.");

	std::string ext0 = stb.get_extension(A.mRefHapFilename);
	if (ext0 == "bcf" || ext0 == "vcf") A.mInputFormat = InputFormat::BCF;
	else if (ext0 =="gz")
	{
		auto position = A.mRefHapFilename.find_last_of ( '.' );
		if ( position != std::string::npos )
		{
			std::string str0 (A.mRefHapFilename.substr(0, position) );
			std::string ext1 = stb.get_extension(str0);
			if (ext1 == "vcf") A.mInputFormat = InputFormat::BCF;
			else vrb.error("Reference panel file format is not supported.");
		}
		else vrb.error("Reference panel file format is not supported.");
	}
	else A.mInputFormat = InputFormat::GLIMPSE; //any extension except vcf/bcf/vcf.gz

	//output
	ext0 = stb.get_extension(A.mOutputFilename);
	if (ext0 == "bcf" || ext0 == "vcf")
	{
		if (ext0 == "vcf")
		{
			A.mOutputFormat=OutputFormat::VCF;
			A.mOutputCompression=OutputCompression::NONE;
		}
		else
		{
			A.mOutputFormat=OutputFormat::BCF;
			A.mOutputCompression=OutputCompression::ZLIB;
		}
	}
	else if (ext0 =="gz")
	{
		auto position = A.mOutputFilename.find_last_of ( '.' );
		if ( position != std::string::npos )
		{
			std::string str0 (A.mOutputFilename.substr(0, position) );
			std::string ext1 = stb.get_extension(str0);
			if (ext1 == "vcf")
			{
				A.mOutputFormat = OutputFormat::VCF;
				A.mOutputCompression = OutputCompression::ZLIB;
			}
			else vrb.error("Output file format is not supported.");
		}
		else vrb.error("Output file format is not supported.");
	}
	else if (ext0 =="bgen")
	{
		A.mOutputFormat = OutputFormat::BGEN;
		if (A.mBgenCompression == "zstd") A.mOutputCompression = OutputCompression::ZSTD;
		else if (A.mBgenCompression == "zlib") A.mOutputCompression = OutputCompression::ZLIB;
		else if (A.mBgenCompression == "no") A.mOutputCompression = OutputCompression::NONE;
		else vrb.error("Output compression must have a value of the following: [no, zlib, zstd].");

		if (A.mBgenNbits >32) vrb.error("Parameter number of bits for BGEN file format encoding is invalid");
	}
	else vrb.error("Output file format is not supported.");

	#ifdef __BGEN__
	#else
		if (A.mOutputFormat == OutputFormat::BGEN) vrb.error("GLIMPSE has not been compiled with BGEN support for output file. Please recompile the program with the BGEN_SUPPORT=YES flag or switch to VCF/BCF.");
	#endif

	A.mPrintGP = !options.count("no-out-gp-field");
	A.mPrintDS = !options.count("no-out-ds-field");
	A.mPrintAP = options.count("out-ap-field");

	if (options.count("contigs-fai"))
	{
		faidx_t *fai = fai_load3(A.mContigFaiFilename.c_str(),A.mContigFaiFilename.c_str(),NULL,FAI_FASTA);
		if ( !fai ) vrb.error("Could not parse " + A.mContigFaiFilename);
		fai_destroy(fai);
	}

	if (A.mInputFormat != InputFormat::GLIMPSE)
	{
		if (!options.count("input-region"))
			vrb.error("You must specify a region to phase using --input-region (this is given by GLIMPSE_chunk)");

		if (!options.count("output-region"))
			vrb.error("You must specify a region to output using --output-region (this is given by GLIMPSE_chunk)");

		if (!options.count("sparse-maf"))
			vrb.error("You must specify a sparse MAF threshold (e.g. 0.001) for reference panel representation.");

		if (A.mSparseMaf >= 0.5 || A.mSparseMaf < 0) vrb.error("The sparse MAF parameter should not be set too high or low. Ideally within the range [5%-0.01%] MAF.");
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
		input_file state_list(A.mStateListFilename);
		if (state_list.fail()) vrb.error("Impossible to open file: " + A.mStateListFilename);
		if (A.Kpbwt>0 && A.mPbwtCM <= 0) vrb.error("Option --pbwt-modulo-cm is not present or set to a non-positive value");
	}
	else if (A.mPbwtCM <= 0) vrb.error("Option --pbwt-modulo-cm is not present or set to a non-positive value");

	if (A.mMinPbwtDepth > A.mMaxPbwtDepth) vrb.error("Min PBWT depth is larger than max PBWT depth parameter.");
	if (A.mErrImp > 1e-3 || A.mErrImp < 1e-13) vrb.error("Error imp parameter must be defined within the range [1e-12, 1e-3]");
	if (A.mErrPhase > 1e-3 || A.mErrPhase < 1e-13) vrb.error("Error phase parameter must be defined within the range [1e-12, 1e-3]");
}

void caller::verbose_files() {
	vrb.title("Files:");
	std::array<std::string, 3> fmt2string = { {"VCF","BCF","BGEN"} };
	std::array<std::string, 3> compr2string = { {"NO","ZLIB","ZSTD"} };
	std::string n_bits = "";
	if (A.mOutputFormat == OutputFormat::BGEN) n_bits = " - " + stb.str(A.mBgenNbits) + " bits";

	if (options.count("bam-list"))
	{
		vrb.bullet("List BAM/CRAM        : [" + A.mBamListFilename + "]");
	}
	else if(options.count("bam-file"))
	{
		vrb.bullet("Input BAM/CRAM       : [" + A.mBamFileFilename + "]");
	}
	else
	{
		vrb.bullet("Input VCF/BCF likel. : [" + A.mInputGLFilename + " | FORMAT/" + (A.mInputFieldGL? "GL":"PL") + " field]");
	}


	if (options.count("fasta")) vrb.bullet("Reference seq. fasta : [" + A.mFastaFilename + "]");
	else vrb.bullet("Reference seq. fasta : [not declared]");

	if (A.mInputFormat == InputFormat::BCF) vrb.bullet("Reference VCF        : [" + A.mRefHapFilename + "]");
	else vrb.bullet("Reference binary     : [" + A.mRefHapFilename + "]");

	if (options.count("map")) vrb.bullet("Genetic Map          : [" + A.mMapFilename + "]");
	else vrb.bullet("Genetic Map          : [not declared]");

	/*
	vrb.bullet("Output file          : [" + A.mOutputFilename + "]");
	vrb.bullet("Output format        : [" + fmt2string[static_cast<int>(A.mOutputFormat)] + " format" + n_bits + " | " + compr2string[static_cast<int>(A.mOutputCompression)] + " compression]");
	*/
	verbose_output_file();

	if (options.count("log")) vrb.bullet("Output LOG           : [" + A.mLogFilename + "]");
	else vrb.bullet("Output LOG           : [not declared]");
}


void caller::verbose_output_file() {
	std::array<std::string, 3> fmt2string = { {"VCF","BCF","BGEN"} };
	std::array<std::string, 3> compr2string = { {"NO","ZLIB","ZSTD"} };
	std::string n_bits = "";
	std::string phased_bgen = "";
	if (A.mOutputFormat == OutputFormat::BGEN)
	{
		n_bits = " - " + stb.str((int)A.mBgenNbits) + " bits";
		phased_bgen = "Unphased";//A.mPrintAP ? "Phased " : "Unphased ";
	}
	std::string index_bcf = "";
	if (A.mOutputFormat != OutputFormat::BGEN)
	{
		index_bcf = A.mOutputIndex ? " | with CSI index" : " | no CSI index";
	}
	vrb.bullet("Output file          : [" + A.mOutputFilename + "]");
	vrb.bullet("Output format        : [" + phased_bgen + fmt2string[static_cast<int>(A.mOutputFormat)] + " format" + n_bits + " | " + compr2string[static_cast<int>(A.mOutputCompression)] + " compression" + index_bcf + "]");

	if (A.mOutputFormat == OutputFormat::VCF && A.mOutputCompression == OutputCompression::NONE)
		vrb.warning("Plain VCF file format with no compression detected as output: this format is strongly discouraged as does not allow the creation of a CSI index, can create very large files. Please use a compressed file format if possible (vcf.gz, bcf).");
}
/*
void caller::verbose_files() {
	vrb.title("Files:");
	vrb.bullet("Input phased array  : [" + options["g"].as < std::string > () + "]");
	vrb.bullet("Reference panel     : [" + options["h"].as < std::string > () + "]");

	if (options.count("map"))
		vrb.bullet("Genetic Map         : [" + options["m"].as < std::string > () + "]");

	verbose_output_file();

	if (options.count("l")) vrb.bullet("Output LOG           : [" + options["l"].as < std::string > () + "]");
}
*/

void caller::verbose_options()
{
	std::array<std::string,2> no_yes = {"NO","YES"};
	std::string opt_input_region, opt_output_region, opt_sparse_maf, opt_map, opt_k_init, opt_k_pbwt, opt_pbwt_min_depth, opt_pbwt_max_depth, opt_pbwt_modulo_cm, opt_state_list, opt_samples_file, opt_call_indels;

	opt_input_region = A.mInputFormat == InputFormat::BCF ? A.mInputRegion : "Given by binary reference panel";
	opt_output_region = A.mInputFormat == InputFormat::BCF ? A.mOutputRegion : "Given by binary reference panel";
	opt_sparse_maf = A.mInputFormat == InputFormat::BCF ? stb.str(A.mSparseMaf) : "Given by binary reference panel";
	if (A.mInputFormat == InputFormat::BCF) opt_map = options.count("map") ? "Given by genetic map" : "Constant rate of 1cM/Mb";
	else opt_map = "Given by binary reference panel";

	opt_k_init = A.Kinit > 0 ? stb.str(A.Kinit) : std::string("Init selection not enabled") ;
	opt_k_pbwt = A.Kpbwt > 0 ? stb.str(A.Kpbwt) : std::string("PBWT selection not enabled");
	opt_pbwt_min_depth = A.Kpbwt > 0 ? stb.str(A.mMinPbwtDepth) : std::string("PBWT selection not enabled");
	opt_pbwt_max_depth = A.Kpbwt > 0 ? stb.str(A.mMaxPbwtDepth) : std::string("PBWT selection not enabled");
	opt_pbwt_modulo_cm = A.Kpbwt > 0 ? stb.str(A.mPbwtCM) : std::string("PBWT selection not enabled");
	opt_state_list = options.count("state-list") ? A.mStateListFilename : std::string("not declared");
	opt_samples_file = options.count("samples-file") ? "Given by samples file" : "Only diploid samples in region";
	opt_call_indels = A.mCallIndels ? "Perform calling" : "Haplotype scaffold";

	vrb.title("GLIMPSE_phase parameters:");
	vrb.bullet("Imputation model     : [Reference panel imputation]");
	vrb.bullet("Input region         : [" + opt_input_region + "]");
	vrb.bullet("Output region        : [" + opt_output_region + "]");
	vrb.bullet("Sparse MAF           : [" + opt_sparse_maf + "]");
	vrb.bullet("Recombination rates  : [" + opt_map + "]");
	vrb.bullet("Ploidy               : [" + opt_samples_file + "]");
	vrb.bullet("Keep monom. ref sites: [" + no_yes[A.mKeepMonomorphicRefSites] + "]");
	vrb.bullet("Use GL indels        : [" + no_yes[A.mUseGLIndels] + "]");

	vrb.title("Model parameters:");
	vrb.bullet("#Burnin iterations   : [" + stb.str(A.mBurnIn) + "]");
	vrb.bullet("#Main iterations     : [" + stb.str(A.mMain) + "]");
	vrb.bullet("Ne [eff. pop. size]  : [" + stb.str(A.mNe) + "]");
	vrb.bullet("Phase error rate     : [" + stb.str(A.mErrPhase) + "]");
	vrb.bullet("Imputation error rate: [" + stb.str(A.mErrImp) + "]");
	vrb.bullet("Min value for hap GLs: [" + stb.str(A.mMinGL) + "]");

	vrb.title("Selection parameters:");
	vrb.bullet("K init               : [" + opt_k_init + "]");
	vrb.bullet("K pbwt               : [" + opt_k_pbwt + "]");
	vrb.bullet("PBWT depth min       : [" + opt_pbwt_min_depth + "]");
	vrb.bullet("PBWT depth max       : [" + opt_pbwt_max_depth + "]");
	vrb.bullet("PBWT modulo (cM)     : [" + opt_pbwt_modulo_cm + "]");
	vrb.bullet("State list           : [" + opt_state_list + "]");

	if (!options.count("input-gl"))
	{

		vrb.title("Genotype calling:");
		vrb.bullet("Calling model        : [" + A.mCallModel + "]");
		vrb.bullet("Indels model         : [" + opt_call_indels + "]");

		vrb.title("BAM/CRAM filters and options:");
		vrb.bullet("Min mapping quality  : [" + stb.str(A.mMapQ) + "]");
		vrb.bullet("Min base quality     : [" + stb.str(A.mBaseQ) + "]");
		vrb.bullet("Max depth            : [" + stb.str(A.mMaxDepth) + "]");
		vrb.bullet("Keep failed QC       : [" + no_yes[A.mKeepFailedQC] + "]");
		vrb.bullet("Keep orphan reads    : [" + no_yes[A.mKeepOrphanReads] + "]");
		vrb.bullet("Keep duplicate reads : [" + no_yes[A.mKeepDuplicates] + "]");
		vrb.bullet("Check pairing        : [" + no_yes[A.mCheckProperPairing] + "]");
		vrb.bullet("Ignore orientation   : [" + no_yes[A.mIgnoreOrientation] + "]");
		vrb.bullet("Illumina-1.3+        : [" + no_yes[A.mIllumina13] + "]");
	}
	else
	{
		vrb.title("Genotype calling: NOT performed");
		vrb.bullet("Calling model        : [-]");
		vrb.bullet("Indels model         : [-]");

		vrb.title("BAM/CRAM filters and options:");
		vrb.bullet("Min mapping quality  : [-]");
		vrb.bullet("Min base quality     : [-]");
		vrb.bullet("Max depth            : [-]");
		vrb.bullet("Keep failed QC       : [-]");
		vrb.bullet("Keep orphan reads    : [-]");
		vrb.bullet("Keep duplicate reads : [-]");
		vrb.bullet("Check pairing        : [-]");
		vrb.bullet("Ignore orientation   : [-]");
		vrb.bullet("Illumina-1.3+        : [-]");
	}

	vrb.title("Other parameters");
	vrb.bullet("Seed                 : [" + stb.str(A.mSeed) + "]");
	vrb.bullet("#Threads             : [" + stb.str(A.mNumThreads) + "]");

}
