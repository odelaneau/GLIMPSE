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
#include <rephaser/rephaser_header.h>

#include <io/genotype_reader.h>
#include <io/gmap_reader.h>

void rephaser::read_files_and_initialise() {
	vrb.title("Initialization:");

	//step0: Initialize seed & other
	rng.setSeed(options["seed"].as < int > ());
	if (options["thread"].as < int > () > 1) {
		i_workers = 0; i_jobs = 0;
		id_workers = vector < pthread_t > (options["thread"].as < int > ());
		pthread_mutex_init(&mutex_workers, NULL);
	}

	//step2: Read input files
	genotype_reader readerG(H, V, options["region"].as < string > ());
	readerG.scanGenotypes(options["input"].as < string > ());
	readerG.allocateData();
	readerG.readGenotypes(options["input"].as < string > ());

	//step3: Read and initialise genetic map
	if (options.count("map")) {
		gmap_reader readerGM;
		readerGM.readGeneticMapFile(options["map"].as < string > ());
		V.setGeneticMap(readerGM);
	} else V.setGeneticMap();
	vrb.bullet("Region spans " + stb.str(V.length()) + " bp and " + stb.str(V.lengthcM(), 2) + " cM");

	//step4
	HMM = vector < rephase_hmm * > (options["thread"].as < int > (), NULL);
	COND = vector < conditioning_set * > (options["thread"].as < int > (), NULL);
	for (int t = 0 ; t < HMM.size() ; t ++) HMM[t] = new rephase_hmm(H, M);

	//step5
	H.mapRareHets(options["maf"].as < double > ());
	M.initialise(V, H, 20000);
	H.initPositionalBurrowWheelerTransform(options["pbwt-depth"].as < int > (), options["pbwt-modulo"].as < int > ());
}
