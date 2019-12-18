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

#include <io/genotype_reader.h>
#include <io/genotype_writer.h>
#include <io/gmap_reader.h>

void caller::read_files_and_initialise() {
	vrb.title("Initialization:");

	//step0: Initialize seed & other
	rng.setSeed(options["seed"].as < int > ());
	if (options["thread"].as < int > () > 1) {
		i_workers = 0; i_jobs = 0;
		id_workers = std::vector < pthread_t > (options["thread"].as < int > ());
		pthread_mutex_init(&mutex_workers, NULL);
	}

	//step1: Parsing region std::string
	buildCoordinates();

	//step2: Read input files
	genotype_reader readerG(H, G, V, gregion, options["maf-filter"].as< float > ());
	readerG.scanGenotypes(options["input"].as < std::string > (), options["reference"].as < std::string > ());
	readerG.allocateGenotypes();
	readerG.readGenotypes(options["input"].as < std::string > (), options["reference"].as < std::string > ());

	perform_delayed_imputation = (options["maf-filter"].as< float >() > 0.0f);

	if (perform_delayed_imputation)
	{
		TPROB = std::vector < probability_set * > (readerG.n_main_samples);
		for (int t = 0 ; t < readerG.n_main_samples ; t ++) {
			TPROB[t] = new probability_set (readerG.n_variants);
		}
	}

	//step3: Read and initialise genetic map
	if (options.count("map")) {
		//gmap_reader readerGM;
		readerGM.readGeneticMapFile(options["map"].as < std::string > ());
		V.setGeneticMap(readerGM);
	}

	//step4
	HP0 = std::vector < std::vector < float > > (options["thread"].as < int > (), std::vector < float > (H.n_site * 2, 0.0));
	HP1 = std::vector < std::vector < float > > (options["thread"].as < int > (), std::vector < float > (H.n_site * 2, 0.0));
	HLC = std::vector < std::vector < float > > (options["thread"].as < int > (), std::vector < float > (H.n_site * 2, 0.0));
	HMM = std::vector < haplotype_hmm * > (options["thread"].as < int > (), NULL);
	DMM = std::vector < diplotype_hmm * > (options["thread"].as < int > (), NULL);
	FMM = std::vector < switchandflipphasing * > (options["thread"].as < int > (), NULL);
	SMM = std::vector < switchphasing * > (options["thread"].as < int > (), NULL);
	COND = std::vector < conditioning_set * > (options["thread"].as < int > (), NULL);
	for (int t = 0 ; t < HMM.size() ; t ++) {
		COND[t] = new conditioning_set(V, (readerG.n_ref_samples+readerG.n_main_samples)*2, 20000);
		HMM[t] = new haplotype_hmm(&H,COND[t],perform_delayed_imputation?TPROB[t]:nullptr);
		DMM[t] = new diplotype_hmm(&H,COND[t]);
		FMM[t] = new switchandflipphasing (&H,COND[t]);
		SMM[t] = new switchphasing (&H,COND[t]);
	}

	//step5
	H.refonly_pbwt = options.count("refonly-select");
	H.store_estimate_HS = options.count("output-HS");
	H.initPositionalBurrowWheelerTransform(options["pbwt-depth"].as < int > (), options["pbwt-modulo"].as < int > ());
}
