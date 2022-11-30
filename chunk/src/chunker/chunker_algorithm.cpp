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

#include <chunker/chunker_header.h>
#include <io/gmap_reader.h>

void chunker::split(output_file & fd, int & cidx, std::string & chr, int start_idx, int stop_idx) {
	// Compute current window properties
	int curr_window_count = stop_idx - start_idx + 1;
	int curr_window_size = positions_all_mb[stop_idx] - positions_all_mb[start_idx]+1;


	// Check if we can further split
	int mid_idx = start_idx + curr_window_count / 2;
	while (mid_idx+1 < positions_all_mb.size() && positions_all_mb[mid_idx+1] == positions_all_mb[mid_idx])	++mid_idx;

	int next0_window_count = mid_idx - start_idx + 1;
	int next1_window_count = stop_idx - mid_idx + 1;

	int next0_window_size = positions_all_mb[mid_idx] - positions_all_mb[start_idx]+1;
	int next1_window_size = positions_all_mb[stop_idx] - positions_all_mb[mid_idx + 1]+1;

	int next0_window_okay = (next0_window_count >= window_count && next0_window_size >= window_mb);
	int next1_window_okay = (next1_window_count >= window_count && next1_window_size >= window_mb);

	// If we can further split, recursion!
	if (next0_window_okay && next1_window_okay) {
		vrb.bullet("Internal window [" + chr + ":" + stb.str(positions_all_mb[start_idx]) + "-" + stb.str(positions_all_mb[stop_idx]) + "] / L=" + stb.str(curr_window_size) + "bp / C=" + stb.str(curr_window_count));
		split (fd, cidx, chr, start_idx, mid_idx);
		split (fd, cidx, chr, mid_idx+1, stop_idx);
	} else {
		//Get left buffer
		int left_idx = -1, left_size = -1, left_count = -1;
		if (start_idx > buffer_count) {
			left_idx = start_idx - buffer_count;
			do {
				left_idx --;
				left_count = start_idx - left_idx + 1;
				left_size = positions_all_mb[start_idx] - positions_all_mb[left_idx] + 1;
			} while (((left_idx > 0) && (left_count < buffer_count)) || (left_size < buffer_mb));
		} else { left_idx = 0; }
		//Get right buffer
		int right_idx = -1, right_size = -1, right_count = -1;
		if (stop_idx < (positions_all_mb.size() - buffer_count)) {
			right_idx = stop_idx + buffer_count - 1;
			do {
				right_idx ++;
				right_count = right_idx - stop_idx + 1;
				right_size = positions_all_mb[right_idx] - positions_all_mb[stop_idx] + 1;
			} while (((right_idx < (positions_all_mb.size() - 1)) && (right_count < buffer_count)) || (right_size < buffer_mb));
		} else { right_idx = positions_all_mb.size() - 1; }
		//Process window
		vrb.bullet("Terminal window [" + stb.str(cidx) + "] -buffer:[" + chr + ":" + stb.str(positions_all_mb[start_idx]) + "-" + stb.str(positions_all_mb[stop_idx]) + "] / +buffer:[" + chr +  ":" + stb.str(positions_all_mb[left_idx]) + "-" + stb.str(positions_all_mb[right_idx]) + "] / L=" + stb.str(curr_window_size) + "bp / C=" + stb.str(curr_window_count));

		//if (report_common_variants)
			//fd << cidx << "\t" << chr << "\t" << chr << ":"<< positions_all_mb[left_idx] << "-" << positions_all_mb[right_idx] << "\t" << chr << ":" << positions_all_mb[start_idx] << "-" << positions_all_mb[stop_idx] << "\t" << curr_window_size << "\t" << curr_window_count << std::endl;
		//else
			fd << cidx << "\t" << chr << "\t" << chr << ":"<< positions_all_mb[left_idx] << "-" << positions_all_mb[right_idx] << "\t" << chr << ":" << positions_all_mb[start_idx] << "-" << positions_all_mb[stop_idx] << "\t" << curr_window_size << "\t" << curr_window_count << std::endl;
		cidx ++;
	}
}



void chunker::chunk() {
	//Initialize
	rng.setSeed(options["seed"].as < int > ());
	sparse_maf = options["sparse-maf"].as < float > ();

	window_cm = options["window-cm"].as < float > ();
	window_mb = (int) (options["window-mb"].as < float > ()*1e6);
	window_count = options["window-count"].as < int > ();
	buffer_cm = options["buffer-cm"].as < float > ();
	buffer_mb = (int) (options["buffer-mb"].as < float > ()*1e6);
	buffer_count = options["buffer-count"].as < int > ();
	report_common_variants = options.count("report-common-variants");

	//Read input data (overlapping coordinates)
	readData(options["input"].as < std::string > (), options["region"].as < std::string > (), options["threads"].as <int> ());

	if (options.count("map"))
	{
		gmap_reader readerGM;
		readerGM.readGeneticMapFile(options["map"].as < std::string > ());
		setGeneticMap(readerGM, map_positions_all, positions_all_mb, positions_all_cm);
	} else
	{
		setGeneticMap(map_positions_all, positions_all_mb, positions_all_cm);
	}
		//V.setGeneticMap();

	if (positions_all_mb.size() == 0) vrb.error("No markers in region: " + options["region"].as < std::string > () + ". Check chromosome name and start/end positions.");
	// Perform chunking!
	int cidx = 0;
	vrb.title("Splitting data into chunks and writting to [" + options["output"].as < std::string > () + "]");
	output_file fd(options["output"].as < std::string > ());

	if (options.count("recursive"))
	{
		split(fd, cidx, chrID, 0, positions_all_mb.size() - 1);
	}
	else if (options.count("sequential"))
	{
		//TODO
	}

	vrb.bullet("#chunks = " + stb.str(cidx));
	fd.close();

	//Finalize
	vrb.title("Total running time = " + stb.str(tac.abs_time()) + " seconds");
}
