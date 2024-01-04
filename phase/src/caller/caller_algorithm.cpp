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
#include "boost/serialization/serialization.hpp"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <filesystem>

const std::string stage_names[3] = {"Init", "Burn-in", "Main"};

void * phase_callback(void * ptr) {
	caller * S = static_cast< caller * >(ptr);
	int id_worker, id_job;
	float prog_step = 1.0/S->G.n_ind;
	const int n_ind = S->G.n_ind;
	pthread_mutex_lock(&S->mutex_workers);
	id_worker = S->i_workers ++;
	pthread_mutex_unlock(&S->mutex_workers);
	for(;;) {
		pthread_mutex_lock(&S->mutex_workers);
		id_job = S->i_jobs ++;
		if (id_job <= n_ind) vrb.progress("  * HMM imputation", id_job*prog_step);
		pthread_mutex_unlock(&S->mutex_workers);
		if (id_job < n_ind) S->phase_individual(id_worker, id_job);
		else pthread_exit(NULL);
	}
	return nullptr;
}

void caller::phase_individual(const int id_worker, const int id_job) {
	const int ploidy = G.vecG[id_job]->ploidy;

	COND[id_worker]->select(id_job, current_stage);

	if (current_stage == STAGE_INIT) G.vecG[id_job]->initHaplotypeLikelihoods(HLC[id_worker],min_gl);
	else
	{
		if (ploidy > 1) G.vecG[id_job]->makeHaplotypeLikelihoods(HLC[id_worker], true,min_gl);
		else G.vecG[id_job]->initHaplotypeLikelihoods(HLC[id_worker],min_gl);
	}
	HMM[id_worker]->computePosteriors(HLC[id_worker], G.vecG[id_job]->flat, HP0[id_worker]);
	G.vecG[id_job]->sampleHaplotypeH0(HP0[id_worker]);
	if (ploidy > 1)
	{
		G.vecG[id_job]->makeHaplotypeLikelihoods(HLC[id_worker], false,min_gl);
		HMM[id_worker]->computePosteriors(HLC[id_worker], G.vecG[id_job]->flat, HP1[id_worker]);
		G.vecG[id_job]->sampleHaplotypeH1(HP1[id_worker]);
		DMM[id_worker]->rephaseHaplotypes(G.vecG[id_job]->H0, G.vecG[id_job]->H1, G.vecG[id_job]->flat);
	}

	if (current_stage == STAGE_MAIN)
	{
		if (ploidy > 1) G.vecG[id_job]->storeGenotypePosteriorsAndHaplotypes(HP0[id_worker], HP1[id_worker]);
		else G.vecG[id_job]->storeGenotypePosteriorsAndHaplotypes(HP0[id_worker]);
	}

	if (options["threads"].as < int > () > 1) pthread_mutex_lock(&mutex_workers);
	statH.push(COND[id_worker]->n_states*1.0);
	statC.push(COND[id_worker]->polymorphic_sites.size() * 100.0/COND[id_worker]->n_tot_sites);
	if (options["threads"].as < int > () > 1) pthread_mutex_unlock(&mutex_workers);
}

void caller::phase_iteration() {
	if (current_stage == STAGE_INIT) {
		vrb.title("Initializing iteration");

		H.initRareTar(G,V);
		H.performSelection_RARE_INIT_GL(V); //not parallel
	} else {
		vrb.title(stb.str(stage_names[current_stage]) + " iteration [" + stb.str(current_iteration+1) + "/" + stb.str(iterations_per_stage[current_stage]) + "]");
		H.updateHaplotypes(G);
		H.transposeRareTar();
		H.matchHapsFromCompressedPBWTSmall(V, current_stage == STAGE_MAIN);
	}

	tac.clock();
	int n_thread = options["threads"].as < int > ();
	i_workers = 0; i_jobs = 0;
	statH.clear();
	statC.clear();

	float prog_step = 1.0/G.n_ind;
	float prog_bar = 0.0;

	if (n_thread > 1) {
		for (int t = 0 ; t < n_thread ; t++) pthread_create( &id_workers[t] , NULL, phase_callback, static_cast<void *>(this));
		for (int t = 0 ; t < n_thread ; t++) pthread_join( id_workers[t] , NULL);
	} else for (int i = 0 ; i < G.n_ind ; i ++)
	{
		phase_individual(0, i);
		prog_bar+=prog_step;
		vrb.progress("  * HMM imputation", prog_bar);
	}
	vrb.bullet("HMM imputation [#states=" + stb.str(statH.mean(), 1) + " / %poly=" + stb.str(statC.mean(), 1) + "%] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");


	write_checkpoint();

	if (current_stage == STAGE_INIT) {
		H.init_states.clear();
		H.init_states.shrink_to_fit();
	}
}

void caller::increment_iteration() {
	current_iteration++;
	while (current_iteration >= iterations_per_stage[current_stage] && current_stage<=STAGE_MAIN) {
		current_stage++;
		current_iteration = 0;
	}
}

void caller::write_checkpoint() {
	if (options.count("checkpoint-file-out")) {
		
		vrb.bullet("writing out checkpoint");
		std::string cp_filename = options["checkpoint-file-out"].as < std::string > ();
		std::string tmp_cp_filename = cp_filename + ".tmp";
		std::ofstream ofs(tmp_cp_filename, std::ios::binary | std::ios_base::out);
		boost::archive::binary_oarchive oa(ofs);
		tac.clock();
		oa << crc.get_value();
		oa << current_stage;
		oa << current_iteration;
		oa << iterations_per_stage[STAGE_BURN];
		oa << options["ne"].as<int>();
		oa << options["min-gl"].as<float>();
		oa << options["err-imp"].as<float>();
		oa << options["err-phase"].as<float>();
		oa << options["pbwt-depth"].as<int>();
		oa << options["pbwt-modulo-cm"].as<float>();
		oa << options["Kinit"].as<int>();
		oa << options["Kpbwt"].as<int>();
		G.serialize_checkpoint_data(oa);
		std::filesystem::rename(tmp_cp_filename.c_str(), cp_filename.c_str());
		vrb.bullet("checkpoint completed (" + stb.str(tac.rel_time(), 2) + "ms)");
	}
}

void caller::read_checkpoint_if_available() {
	if (options.count("checkpoint-file-in")) {
		vrb.bullet("reading checkpoint");
		std::string cp_filename = options["checkpoint-file-in"].as < std::string > ();
		std::ifstream ifs(cp_filename, std::ios::binary | std::ios_base::in);
		boost::archive::binary_iarchive ia(ifs);
		unsigned long long checkpoint_crc;
		ia >> checkpoint_crc;
		if (checkpoint_crc != crc.get_value()) {
			vrb.error("Input data checksum in checkpoint file does not match "
			"checksum of input data for this run.");
		}
		ia >> current_stage;
		ia >> current_iteration;
		int checkpoint_burnin_iterations;
		ia >> checkpoint_burnin_iterations;
		if (current_iteration >= iterations_per_stage[current_stage]) {
			std::stringstream err_str;
			err_str<<"Checkpoint file has already run "<<current_iteration + 1<<" iterations "
			"for stage"<<stage_names[current_stage]<<", and this run only calls for "<<iterations_per_stage[current_stage]<< 
			". This run must call for at least a many iterations as the checkpoint file already ran in order "
			"to use this checkpoint file.";
			vrb.error(err_str.str());
		}
		if (current_stage == STAGE_MAIN && checkpoint_burnin_iterations != iterations_per_stage[STAGE_BURN]) {
			std::stringstream err_str;
			err_str<<"Checkpoint file is in Main stage, and ran "<<checkpoint_burnin_iterations<<" burn-in iterations, while "
			"this run calls for "<<iterations_per_stage[STAGE_BURN]<<" burn-in iterations.  These values must be "
			"equal to use this checkpoint file.";
			vrb.error(err_str.str());
		}
		confirm_checkpoint_param<int>(ia, "ne");
		confirm_checkpoint_param<float>(ia, "min-gl");
		confirm_checkpoint_param<float>(ia, "err-imp");
		confirm_checkpoint_param<float>(ia, "err-phase");
		confirm_checkpoint_param<int>(ia, "pbwt-depth");
		confirm_checkpoint_param<float>(ia, "pbwt-modulo-cm");
		confirm_checkpoint_param<int>(ia, "Kinit");
		confirm_checkpoint_param<int>(ia, "Kpbwt");
		G.serialize_checkpoint_data(ia);
		vrb.bullet("checkpoint read");
	}
}

void caller::phase_loop() {
	//steup iteration counters
	current_stage = STAGE_INIT;
	current_iteration = -1; //increment will make current_iteration 0 if no checkpoint loaded

	read_checkpoint_if_available();

	increment_iteration();

	while (current_stage <= STAGE_MAIN) {
		phase_iteration();
		increment_iteration();
	}

	//Finalization
	//vrb.title("Finalization");
	for (int i = 0 ; i < G.vecG.size() ; i ++)
		 G.vecG[i]->sortAndNormAndInferGenotype();

	//vrb.bullet("done");
}

