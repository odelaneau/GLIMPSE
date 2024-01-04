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

#include <caller/caller_header.h>

#include <io/genotype_reader.h>
#include <io/genotype_writer.h>
#include <io/gmap_reader.h>
#include <io/genotype_bam_caller.h>
#include <boost/archive/binary_iarchive.hpp>
#include <containers/glimpse_mpileup.h>

void caller::print_ref_panel_info(const std::string ref_string)
{
	vrb.bullet(ref_string + " reference panel [Nrh=" + stb.str(H.n_ref_haps) + "] [L=" + stb.str(H.n_tot_sites) + "] [Lrare= " + stb.str(H.n_rar_sites) + " (" + stb.str(((float)H.n_rar_sites/H.n_tot_sites)*100.0, 1) + "%) - Lcommon= " + stb.str(H.n_com_sites)  + " (" + stb.str(((float)H.n_com_sites/H.n_tot_sites)*100.0, 1) + "%)]");
	vrb.bullet("Input region         : [" + V.input_gregion + "]");
	vrb.bullet("Output region        : [" + V.output_gregion + "]");
	vrb.bullet("Sparse MAF           : [" + stb.str(H.sparse_maf) + "]");
}
void caller::read_files_and_initialise() {
	vrb.title("Initialisation:");

	//step0: Initialize seed & other
	rng.setSeed(options["seed"].as < int > ());
	const int nthreads =  options["threads"].as < int > ();
	iterations_per_stage[STAGE_INIT] = 1;
	iterations_per_stage[STAGE_BURN] = options["burnin"].as < int > ();
	iterations_per_stage[STAGE_MAIN] = options["main"].as < int > ();

	if (nthreads < 1) vrb.error("Error defining the number of threads. Only positive values are accepted.");

	if (nthreads > 1) {
		i_workers = 0; i_jobs = 0;
		id_workers = std::vector < pthread_t > (nthreads);
		pthread_mutex_init(&mutex_workers, NULL);
	}
	min_gl = options["min-gl"].as < float > ();

	//step2: Read input files
	std::string reference_filename = options["reference"].as<std::string>();
	if (input_fmt == InputFormat::BCF)
	{
		buildCoordinates();

		genotype_reader readerG(H, G, V, M, options["sparse-maf"].as < float > (), options.count("input-field-gl"),options.count("impute-reference-only-variants"), options.count("keep-monomorphic-ref-sites"), options.count("use-gl-indels"));
		if (options.count("samples-file")) readerG.readSamplesFilePloidy(options["samples-file"].as < std::string > ());

		if (options.count("input-gl"))
		{
			readerG.readGenotypes(options["input-gl"].as < std::string > (), options["reference"].as < std::string > (), nthreads);
			print_ref_panel_info("VCF/BCF");
		}
		else if (options.count("bam-list") || options.count("bam-file"))
		{
			setup_mpileup();
			readerG.readGenotypesAndBAMs(options["reference"].as < std::string > (), nthreads);
			print_ref_panel_info("VCF/BCF");
			readerG.set_ploidy_tar();
			read_BAMs();
		}
		else vrb.error("No valid input options (input-gl / bam-list / bam)");
	}
	else
	{
		vrb.wait("  * Binary reference panel parsing");
		tac.clock();
		{
			std::ifstream ifs(reference_filename, std::ios::binary | std::ios_base::in);
			if (!ifs.good()) vrb.error("Reading binary reference panel file: [" + reference_filename + "]. File not good(): eofbit, failbit or badbit set or file not found.");
			try
			{
				boost::archive::binary_iarchive ia(ifs);
				ia >> H;
				ia >> V;
			} catch (std::exception& e ) {
				std::stringstream err_str;
				err_str <<"problems reading the binary reference panel (exception triggered by boost archive). Please ensure you are using the same GLIMPSE and boost library version";
				err_str << e.what();
				vrb.error(err_str.str());
			}
			if (H.Ypacked.size()==0) vrb.error("Problem reading binary file format. Empty PBWT detected.");

			vrb.bullet("Binary reference panel parsing [done] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
			print_ref_panel_info("Binary");
		}

		genotype_reader readerG(H, G, V, M, H.sparse_maf, options.count("input-field-gl"),options.count("impute-reference-only-variants"), options.count("keep-monomorphic-ref-sites"),options.count("use-gl-indels"));
		if (options.count("samples-file")) readerG.readSamplesFilePloidy(options["samples-file"].as < std::string > ());

		if (options.count("input-gl"))
		{
			readerG.readTarGenotypes(options["input-gl"].as < std::string > (), nthreads);
			H.allocate_hap_only();
		}
		else if (options.count("bam-list") || options.count("bam-file"))
		{
			setup_mpileup();
			readerG.set_ploidy_tar();
			H.allocate_hap_only();
			read_BAMs();
		}
		else vrb.error("No valid input options (input-gl / bam-list / bam)");
	}
	H.transposeRareRef();

	//step3: Read and initialise genetic map
	if (input_fmt == InputFormat::BCF)
	{
		if (options.count("map"))
		{
			gmap_reader readerGM;
			readerGM.readGeneticMapFile(options["map"].as < std::string > ());
			V.setGeneticMap(readerGM);
		} else V.setGeneticMap();
	}
	vrb.bullet("Region spans " + stb.str(V.length()) + " bp and " + stb.str(V.lengthcM(), 2) + " cM");

	//step4
	const int ne = options["ne"].as < int > ();
	HP0 = std::vector < std::vector < float > > (nthreads, std::vector < float > (H.n_tot_sites * 2, 0.0));
	if (H.max_ploidy > 1) HP1 = std::vector < std::vector < float > > (nthreads, std::vector < float > (H.n_tot_sites * 2, 0.0));

	HLC = std::vector < std::vector < float > > (nthreads, std::vector < float > (H.n_tot_sites * 2));
	HMM = std::vector < imputation_hmm * > (nthreads, nullptr);
	if (H.max_ploidy > 1) DMM = std::vector < phasing_hmm * > (nthreads, nullptr);
	COND = std::vector < conditioning_set * > (nthreads, nullptr);
	int kpbwt = std::min(options["Kpbwt"].as < int > (), (int) H.n_ref_haps);
	int kinit = std::min(options["Kinit"].as < int > (), (int) H.n_ref_haps);
	const float err_imp = std::clamp(options["err-imp"].as < float > (), 1e-12f, 1e-3f);
	const float err_phase = options["err-phase"].as < float > ();
	const bool use_list = options.count("state-list");

	for (int t = 0 ; t < HMM.size() ; t ++)
	{
		COND[t] = new conditioning_set(V,H,H.n_ref_haps,ne, kinit, kpbwt, err_imp, err_phase, use_list);
		HMM[t] = new imputation_hmm(COND[t]);
		if (H.max_ploidy > 1) DMM[t] = new phasing_hmm(COND[t]);
	}

	//step5 PBWT
	H.allocatePBWT(options["pbwt-depth"].as < int > (), options["pbwt-modulo-cm"].as < float > (), V, G, kinit, kpbwt);

	//step6 list states
	if (use_list) H.read_list_states(options["state-list"].as < std::string > ());

	//checksum
	if(options.count("checkpoint-file-in") || options.count("checkpoint-file-out")) {
		H.update_checksum(crc);
		G.update_checksum(crc);
		V.update_checksum(crc);
	}

}


void caller::setup_mpileup()
{
	if (!options.count("download-fasta-ref"))
	{
	// without these 2 lines, htslib sometimes tries to download a part of the sequence
	// even though the -f reference was provided.
		setenv("REF_CACHE", "", 0);
		setenv("REF_PATH", "fake_value_so_no_download", 0);
	}

	//if (!options.count("fasta") && !options.count("download-fasta-ref")) vrb.error("Fasta file should be provided or --download-fasta-ref option must be added.");
	if (options.count("fasta"))
	{
		M.fai_fname = options["fasta"].as < std::string > ();
		M.fai = fai_load(M.fai_fname.c_str());
		if (M.fai == NULL) vrb.error("Error loading reference genome. The reference should be a fasta file with a .fai index.");
	}

	M.fflag = (BAM_FUNMAP | BAM_FSECONDARY);
	if (!options.count("keep-failed-qc")) M.fflag |= BAM_FQCFAIL;
	if (!options.count("keep-duplicates")) M.fflag |= BAM_FDUP;
	//if (!options.count("keep-supp")) M.fflag  |= BAM_FSUPPLEMENTARY;

	M.keep_orphan = options.count("keep-orphan-reads");
	M.check_orientation = !options.count("ignore-orientation");
	M.check_proper_pair = options.count("check-proper-pairing");

	M.min_mq = options["mapq"].as <int> ();
	M.min_bq = options["baseq"].as<int>();
	M.max_dp = options["max-depth"].as<int>();

	if (M.min_mq <0) vrb.error("Error setting min map quality. Value is negative: " + std::to_string(M.min_mq) + ".");
	if (M.min_bq <0) vrb.error("Error setting min base quality. Value is negative: " + std::to_string(M.min_bq) + ".");
	if (M.max_dp <0) vrb.error("Error setting max deapth. Value is negative: " + std::to_string(M.max_dp) + ".");

	M.illumina13 = options.count("illumina13+");

	if (options["call-model"].as<std::string>()!="standard")
	{
		vrb.error("Only standard model is supported for now");
	}

    if (options.count("bam-list"))
    {
    	std::string buffer;
		input_file fd (options["bam-list"].as < std::string >());
	    if (!fd.good()) vrb.error("Reading bam list: [" + options["bam-list"].as < std::string >() + "]");

		std::vector<std::string> btokens(2);
		std::set<std::string> sbams;
		std::set<std::string> snames;
		while (getline(fd, buffer))
		{
			if (stb.split(buffer, btokens, " 	") > 2 || btokens.size() < 1) vrb.error("Bam list should contain only one or two columns");
			auto ret0 = sbams.insert(btokens[0]);
			if (ret0.second==false) vrb.error("Repeated filename in bam list: " + M.bam_fnames.back());
			M.bam_fnames.push_back(btokens[0]);

			std::string name = btokens[btokens.size() > 1];
			if (btokens.size() == 1) name = stb.remove_ext(stb.extract_file_name(name));
			auto ret1 = snames.insert(name);
			if (ret1.second==false) vrb.error("Repeated sample name in bam list: " + name);
			M.tar_sample_names.push_back(name);
		}
		fd.close();
		M.n_tar_samples = M.bam_fnames.size();
    }
    else
    {
    	if (!options.count("bam-file")) vrb.error("Error reading --bam-file option");
    	M.bam_fnames  = std::vector<std::string>(1, options["bam-file"].as < std::string > ());
    	M.n_tar_samples = 1;
    	M.tar_sample_names = std::vector<std::string>(1);
    	if (options.count("ind-name")) M.tar_sample_names[0] = options["ind-name"].as<std::string>();
    	else M.tar_sample_names[0] = stb.remove_ext(stb.extract_file_name(M.bam_fnames[0]));
    }
    if (M.n_tar_samples == 0) vrb.error("No input BAM file given.");
}

void * read_BAMs_callback(void * ptr) {
	caller * S = static_cast< caller * >(ptr);
	int id_worker, id_job;
	float prog_step = 1.0/S->M.n_tar_samples;
	const int n_ind = S->M.n_tar_samples;
	pthread_mutex_lock(&S->mutex_workers);
	id_worker = S->i_workers ++;
	pthread_mutex_unlock(&S->mutex_workers);
	for(;;) {
		pthread_mutex_lock(&S->mutex_workers);
		id_job = S->i_jobs ++;
		if (id_job <= n_ind) vrb.progress("  * Reading BAM files", id_job*prog_step);
		pthread_mutex_unlock(&S->mutex_workers);
		if (id_job < n_ind) S->READER_BAM[id_worker]->call_mpileup(id_job);
		else pthread_exit(NULL);
	}
	return nullptr;
}


void caller::read_BAMs()
{
	tac.clock();

	vrb.bullet("Reading BAM files");

	const int nthreads = options["threads"].as < int > ();
	READER_BAM = std::vector < genotype_bam_caller * > (nthreads, nullptr);
	for (int t = 0 ; t < nthreads ; t ++)
		READER_BAM[t] = new genotype_bam_caller(H,G,V,M, V.input_gregion, V.input_start, V.input_stop,  options["call-model"].as <std::string> (), options.count("call-indels"));

    G.stats.cov_ind = std::vector<stats1D>(M.n_tar_samples);
    G.stats.depth_count = std::vector< std::vector<int>> (M.n_tar_samples, std::vector<int> (M.max_dp+1, 0));

	float prog_step = 1.0/M.n_tar_samples;
	float prog_bar = 0.0;

	if (nthreads > 1)
	{
		for (int t = 0 ; t < id_workers.size() ; t++) pthread_create( &id_workers[t] , NULL, read_BAMs_callback, static_cast<void *>(this));
		for (int t = 0 ; t < id_workers.size() ; t++) pthread_join( id_workers[t] , NULL);
	} else for (int i = 0; i < M.n_tar_samples; ++i)
	{
		READER_BAM[0]->call_mpileup(i);
		prog_bar+=prog_step;
		vrb.progress("  * Reading BAM files", prog_bar);
	}

	for (int t = 0 ; t < nthreads ; t ++)
	{
		if(READER_BAM[t]) delete READER_BAM[t];
	}
	READER_BAM.clear();

    vrb.bullet("Reading BAM files done (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");

	//STATS HERE!
    std::string out_file_name = options["output"].as < std::string > ();
    std::string file_name = stb.remove_ext(out_file_name);
    if (stb.get_extension(file_name) == "vcf") file_name = stb.remove_ext(file_name);
    file_name += "_stats_coverage.txt.gz";
    output_file out_file(file_name);
	double avg_coverage = 0.0;
	double n_sites = V.vec_pos.size();
	out_file << "BAM/CRAM coverage per individual [observed / expected from Pois(coverage)]\n";
	out_file << "ID\tSample\t\tCov.\tNo data%\tCount0\t\tCount1\t\tCount2\t\tCount3\t\tCount4\t\tCount5-10\tCount11-30\tCount30+\n";
	for (int i=0; i<M.n_tar_samples;++i)
	{
		double mean_cov = G.stats.cov_ind[i].mean();
		avg_coverage += mean_cov;
		if (mean_cov == 0)
			vrb.warning("No mapped read available for sample " +  G.vecG[i]->name + " in the current chunk. Check your BAM file prior to perform phasing/imputation.");
		if (mean_cov < 0.01)
			vrb.warning("Sample: " + G.vecG[i]->name + " has " + std::to_string((round( G.stats.depth_count[i][0]/n_sites* 1000.0))/1000.0).substr(0,5) + "% of missing data. Imputation might be non meaningful for this sample due to the large amount of missing data. Check your sample if this does not agree with your expectations.");

		std::array<int,3> counts = {0,0,0};
		std::array<double,3> pois = {0,0,0};

		int j=0;
		for (int k=5; k<std::min(11,M.max_dp);++k)
		{
			counts[j] +=  G.stats.depth_count[i][k];
			pois[j] += rng.dpois(k,mean_cov);
		}
		++j;
		for (int k=11; k<std::min(31,M.max_dp);++k)
		{
			counts[j] +=  G.stats.depth_count[i][k];
			pois[j] += rng.dpois(k,mean_cov);
		}
		++j;
		for (int k=31; k<M.max_dp;++k)
		{
			counts[j] +=  G.stats.depth_count[i][k];
			pois[j] += rng.dpois(k,mean_cov);
		}

		const std::string tabs = G.vecG[i]->name.size() < 8? "\t\t" : "\t";
		out_file <<
				std::to_string(i).substr(0,6) + "\t" + G.vecG[i]->name.substr(0,15) + tabs + std::to_string((round(G.stats.cov_ind[i].mean()*1000.0)/1000.0)).substr(0,5) + "\t" + std::to_string((round(G.stats.depth_count[i][0]*100/n_sites*100.0)/100.0)).substr(0,5) + "%\t" +
				+ "\t[" + std::to_string((round(G.stats.depth_count[i][0]/n_sites* 1000.0))/1000.0).substr(0,5) + " / "+ std::to_string((round(rng.dpois(0,mean_cov)* 1000.0))/1000.0).substr(0,5) +"]"
				+ "\t[" + std::to_string((round(G.stats.depth_count[i][1]/n_sites* 1000.0))/1000.0).substr(0,5) + " / "+ std::to_string((round(rng.dpois(1,mean_cov)* 1000.0))/1000.0).substr(0,5) +"]"
				+ "\t[" + std::to_string((round(G.stats.depth_count[i][2]/n_sites* 1000.0))/1000.0).substr(0,5) + " / "+ std::to_string((round(rng.dpois(2,mean_cov)* 1000.0))/1000.0).substr(0,5) +"]"
				+ "\t[" + std::to_string((round(G.stats.depth_count[i][3]/n_sites* 1000.0))/1000.0).substr(0,5) + " / "+ std::to_string((round(rng.dpois(3,mean_cov)* 1000.0))/1000.0).substr(0,5) +"]"
				+ "\t[" + std::to_string((round(G.stats.depth_count[i][4]/n_sites* 1000.0))/1000.0).substr(0,5) + " / "+ std::to_string((round(rng.dpois(4,mean_cov)* 1000.0))/1000.0).substr(0,5) +"]"
				+ "\t[" + std::to_string((round(counts[0]/n_sites* 1000.0))/1000.0).substr(0,5) + " / "+ std::to_string((round(pois[0]* 1000.0))/1000.0).substr(0,5) +"]"
				+ "\t[" + std::to_string((round(counts[1]/n_sites* 1000.0))/1000.0).substr(0,5) + " / "+ std::to_string((round(pois[1]* 1000.0))/1000.0).substr(0,5) +"]"
				+ "\t[" + std::to_string((round(counts[2]/n_sites* 1000.0))/1000.0).substr(0,5) + " / "+ std::to_string((round(pois[2]* 1000.0))/1000.0).substr(0,5) +"]"
				+ "\n";
	}
	out_file.close();
	vrb.bullet("Coverage statistics reported in: [" + file_name + "]");
	vrb.bullet("Avg. sequencing coverage: " + std::to_string(avg_coverage/M.n_tar_samples));

    //let's deallocate now. Not used afterwards
    G.stats.cov_ind.clear();
    G.stats.cov_ind.shrink_to_fit();
    G.stats.depth_count.clear();
    G.stats.depth_count.shrink_to_fit();

}
