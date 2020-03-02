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
		id_workers = vector < pthread_t > (options["thread"].as < int > ());
		pthread_mutex_init(&mutex_workers, NULL);
	}

	//step1: Parsing region string
	buildCoordinates();

	//step2: Read input files
	genotype_reader readerG(H, G, V, input_gregion);
	if (options.count("init-pool")) readerG.readInitializingSamples(options["init-pool"].as < string > ());
	readerG.scanGenotypes(options["input"].as < string > (), options["reference"].as < string > ());
	readerG.allocateGenotypes();
	readerG.readGenotypes(options["input"].as < string > (), options["reference"].as < string > ());

	//step3: Read and initialise genetic map
	if (options.count("map")) {
		gmap_reader readerGM;
		readerGM.readGeneticMapFile(options["map"].as < string > ());
		V.setGeneticMap(readerGM);
	} else V.setGeneticMap();
	vrb.bullet("Region spans " + stb.str(V.length()) + " bp and " + stb.str(V.lengthcM(), 2) + " cM");

	//step4
	HP0 = vector < vector < float > > (options["thread"].as < int > (), vector < float > (H.n_site * 2, 0.0));
	HP1 = vector < vector < float > > (options["thread"].as < int > (), vector < float > (H.n_site * 2, 0.0));
	//HP0noPL = vector < vector < float > > (options["thread"].as < int > (), vector < float > (H.n_site * 2, 0.0));
	//HP1noPL = vector < vector < float > > (options["thread"].as < int > (), vector < float > (H.n_site * 2, 0.0));
	//MISMATCH = vector < float > (H.n_site, 0.0f);
	HLC = vector < vector < float > > (options["thread"].as < int > (), vector < float > (H.n_site * 2, 0.0));
	HMM = vector < haplotype_hmm * > (options["thread"].as < int > (), NULL);
	DMM = vector < diplotype_hmm * > (options["thread"].as < int > (), NULL);
	//FMM = vector < switchandflipphasing * > (options["thread"].as < int > (), NULL);
	//SMM = vector < switchphasing * > (options["thread"].as < int > (), NULL);
	COND = vector < conditioning_set * > (options["thread"].as < int > (), NULL);
	for (int t = 0 ; t < HMM.size() ; t ++) {
		COND[t] = new conditioning_set(V,H,(readerG.n_ref_samples+readerG.n_main_samples)*2, 20000);
		HMM[t] = new haplotype_hmm(COND[t]);
		DMM[t] = new diplotype_hmm(COND[t]);
		//FMM[t] = new switchandflipphasing (COND[t]);
		//SMM[t] = new switchphasing (COND[t]);
	}

	//step5
	H.initPositionalBurrowWheelerTransform(options["pbwt-depth"].as < int > (), options["pbwt-modulo"].as < int > ());
}
