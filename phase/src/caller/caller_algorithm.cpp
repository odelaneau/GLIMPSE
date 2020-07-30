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

#include <caller/caller_header.h>

void * phase_callback(void * ptr) {
	caller * S = static_cast< caller * >(ptr);
	int id_worker, id_job;
	pthread_mutex_lock(&S->mutex_workers);
	id_worker = S->i_workers ++;
	pthread_mutex_unlock(&S->mutex_workers);
	for(;;) {
		pthread_mutex_lock(&S->mutex_workers);
		id_job = S->i_jobs ++;
		if (id_job <= S->G.n_ind) vrb.progress("  * HMM imputation", id_job*1.0/S->G.n_ind);
		pthread_mutex_unlock(&S->mutex_workers);
		if (id_job < S->G.n_ind) S->phase_individual(id_worker, id_job);
		else pthread_exit(NULL);
	}
}

void caller::phase_individual(int id_worker, int id_job) {
	if (current_stage == STAGE_INIT) {
		COND[id_worker]->selectRandom(id_job, options["init-states"].as < int > ());
		G.vecG[id_job]->initHaplotypeLikelihoods(HLC[id_worker]);
	} else {
		//DMM[id_worker]->rephaseHaplotypes(G.vecG[id_job]->H0, G.vecG[id_job]->H1);
		COND[id_worker]->selectPBWT(id_job, options["init-states"].as < int > ());
		G.vecG[id_job]->makeHaplotypeLikelihoods(HLC[id_worker], true);
	}
	HMM[id_worker]->computePosteriors(HLC[id_worker], HP0[id_worker]);
	G.vecG[id_job]->sampleHaplotypeH0(HP0[id_worker]);
	G.vecG[id_job]->makeHaplotypeLikelihoods(HLC[id_worker], false);
	HMM[id_worker]->computePosteriors(HLC[id_worker], HP1[id_worker]);
	G.vecG[id_job]->sampleHaplotypeH1(HP1[id_worker]);

	DMM[id_worker]->rephaseHaplotypes(G.vecG[id_job]->H0, G.vecG[id_job]->H1);
	
	if (current_stage == STAGE_MAIN) G.vecG[id_job]->storeGenotypePosteriorsAndHaplotypes(HP0[id_worker], HP1[id_worker]);
	if (options["thread"].as < int > () > 1) pthread_mutex_lock(&mutex_workers);
	statH.push(COND[id_worker]->n_states*1.0);
	statC.push(COND[id_worker]->n_sites * 100.0/COND[id_worker]->Hmono.size());
	if (options["thread"].as < int > () > 1) pthread_mutex_unlock(&mutex_workers);
}

void caller::phase_iteration() {
	tac.clock();
	int n_thread = options["thread"].as < int > ();
	i_workers = 0; i_jobs = 0;
	statH.clear();
	statC.clear();
	if (n_thread > 1) {
		for (int t = 0 ; t < n_thread ; t++) pthread_create( &id_workers[t] , NULL, phase_callback, static_cast<void *>(this));
		for (int t = 0 ; t < n_thread ; t++) pthread_join( id_workers[t] , NULL);
	} else for (int i = 0 ; i < G.n_ind ; i ++) {
		phase_individual(0, i);
		vrb.progress("  * HMM imputation", (i+1)*1.0/G.n_ind);
	}
	vrb.bullet("HMM imputation [#states=" + stb.str(statH.mean(), 1) + " / %poly=" + stb.str(statC.mean(), 1) + "%] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void caller::phase_loop() {
	//First Iteration
	current_stage = STAGE_INIT;
	vrb.title("Initializing iteration");
	phase_iteration();

	//Burn-in 0
	current_stage = STAGE_BURN;
	int nBurnin = options["burnin"].as < int > ();
	for (int iter = 0 ; iter < nBurnin ; iter ++) {
		vrb.title("Burn-in 1 iteration [" + stb.str(iter+1) + "/" + stb.str(nBurnin) + "]");
		H.updateHaplotypes(G);
		H.updatePositionalBurrowWheelerTransform();
		phase_iteration();
	}

	//Main
	current_stage = STAGE_MAIN;
	int nMain = options["main"].as < int > ();
	for (int iter = 0 ; iter < nMain ; iter ++) {
		vrb.title("Main iteration [" + stb.str(iter+1) + "/" + stb.str(nMain) + "]");
		H.updateHaplotypes(G);
		H.updatePositionalBurrowWheelerTransform();
		phase_iteration();
	}

	//Finalization
	vrb.title("Finalization");
	for (int i = 0 ; i < G.vecG.size() ; i ++) G.vecG[i]->sortAndNormAndInferGenotype();
	vrb.bullet("done");
}
