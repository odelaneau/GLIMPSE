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
#include <io/gmap_reader.h>
#include "boost/serialization/serialization.hpp"
#include <boost/archive/tmpdir.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <sstream>
#include <iostream>

void caller::read_files_and_initialise() {
	vrb.title("Splitting reference panel:");

	//step0: Initialize seed & other
	rng.setSeed(options["seed"].as < int > ());
	if (options["threads"].as < int > () > 1) {
		i_workers = 0; i_jobs = 0;
		id_workers = std::vector < pthread_t > (options["threads"].as < int > ());
		pthread_mutex_init(&mutex_workers, NULL);
	}

	buildCoordinates();

	gmap_reader readerGM;
	if (options.count("map"))
		readerGM.readGeneticMapFile(options["map"].as < std::string > ());

	for (int i=0; i<input_gregion.size(); ++i)
	{
		//MULTI?
		variant_map V;
		ref_haplotype_set H;

		V.chrid = chrid[i];
		V.input_start = input_start[i];
		V.input_stop = input_stop[i];
		V.output_start = output_start[i];
		V.output_stop = output_stop[i];
		V.input_gregion = input_gregion[i];
		V.output_gregion = output_gregion[i];

		std::string reg_out = V.input_gregion;
		std::replace( reg_out.begin(), reg_out.end(), ':', '_');
		std::replace( reg_out.begin(), reg_out.end(), '-', '_');

		ref_genotype_reader readerG(H,V, input_gregion[i], options["sparse-maf"].as < float > (), options.count("keep-monomorphic-ref-sites"));
		readerG.readRefPanel(options["reference"].as<std::string>(),options["threads"].as < int > () );

		if (readerGM.pos_cm.size() > 1) V.setGeneticMap(readerGM);
		else V.setGeneticMap();

		std::string output_prefix = options["output"].as<std::string>();
		vrb.bullet("Region spans " + stb.str(V.length()) + " bp and " + stb.str(V.lengthcM(), 2) + " cM");

		H.build_sparsePBWT(V);

		vrb.print("Writing file...");
		tac.clock();
		{
			std::ofstream ofs(output_prefix + "_" + reg_out + ".bin", std::ios::binary | std::ios_base::out);
			boost::archive::binary_oarchive oa(ofs);
			oa << H;
			oa << V;
		}
		vrb.print("Writing file [done] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	}
}
