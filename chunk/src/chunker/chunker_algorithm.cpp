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

void chunker::split_recursive(output_file & fd, long int & cidx, std::string & chr, long int start_idx, long int stop_idx) {
	cnk_info.reset();

	// Compute current window properties
	assert(stop_idx>start_idx);
	long int curr_window_count = stop_idx - start_idx + 1;
	float curr_window_cm_size = positions_all_cm[stop_idx] - positions_all_cm[start_idx];
	long int curr_window_mb_size = positions_all_mb[stop_idx] - positions_all_mb[start_idx];

	// Check if we can further split
	long int mid_idx = start_idx + curr_window_count / 2;
	while (mid_idx+1 < positions_all_mb.size() && positions_all_mb[mid_idx+1] == positions_all_mb[mid_idx])	++mid_idx;

	//set common variants indices using mapping
	long int curr_start_common_idx = all2common[start_idx];
	long int curr_stop_common_idx = all2common[stop_idx] - 1;
	long int next0_start_common_idx = all2common[start_idx];
	long int next0_stop_common_idx = all2common[mid_idx] - 1;
	long int next1_start_common_idx = all2common[mid_idx+1];
	long int next1_stop_common_idx = all2common[stop_idx] - 1;

	if (curr_start_common_idx <0 || curr_stop_common_idx <0 || curr_start_common_idx > positions_common_mb.size() || curr_stop_common_idx > positions_common_mb.size() || curr_stop_common_idx <= curr_start_common_idx)
		vrb.error("Error setting indices current");
	if (next0_start_common_idx <0 || next0_stop_common_idx <0 || next0_start_common_idx > positions_common_mb.size() || next0_stop_common_idx > positions_common_mb.size() || next0_stop_common_idx <= next0_start_common_idx)
		vrb.error("Error setting indices next0");
	if (next1_start_common_idx <0 || next1_stop_common_idx <0 || next1_start_common_idx > positions_common_mb.size() || next1_stop_common_idx > positions_common_mb.size() || next1_stop_common_idx <= next1_start_common_idx)
		vrb.error("Error setting indices next1");

	long int curr_window_common_count = curr_stop_common_idx - curr_start_common_idx + 1;
	long int next0_window_common_count = next0_stop_common_idx - next0_start_common_idx + 1;
	long int next1_window_common_count = next1_stop_common_idx - next1_start_common_idx + 1;

	long int next0_window_count = mid_idx - start_idx + 1;
	long int next1_window_count = stop_idx - mid_idx + 1;

	float next0_window_cm_size = positions_all_cm[mid_idx] - positions_all_cm[start_idx];
	float next1_window_cm_size = positions_all_cm[stop_idx] - positions_all_cm[mid_idx + 1];

	long int next0_window_mb_size = positions_all_mb[mid_idx] - positions_all_mb[start_idx];
	long int next1_window_mb_size = positions_all_mb[stop_idx] - positions_all_mb[mid_idx + 1];

	//long int next0_window_okay = (next0_window_cm_size >= window_cm && next0_window_mb_size >= window_mb && next0_window_count >= window_count);
	//long int next1_window_okay = (next1_window_cm_size >= window_cm && next1_window_mb_size >= window_mb && next1_window_count >= window_count);
	long int next0_window_okay = (next0_window_cm_size >= window_cm && next0_window_mb_size >= window_mb && next0_window_common_count >= window_count);
	long int next1_window_okay = (next1_window_cm_size >= window_cm && next1_window_mb_size >= window_mb && next1_window_common_count >= window_count);

	// If we can further split, recursion!
	if (next0_window_okay && next1_window_okay) {
		vrb.bullet("long internal window [" + chr + ":" + stb.str(positions_all_mb[start_idx]) + "-" + stb.str(positions_all_mb[stop_idx]) + "] / L=" + stb.str(curr_window_cm_size) + "cM / L=" + stb.str(curr_window_mb_size) + "bp / C=" + stb.str(curr_window_count));
		split_recursive (fd, cidx, chr, start_idx, mid_idx);
		split_recursive (fd, cidx, chr, mid_idx+1, stop_idx);
	} else {
		long int left_idx = -1; long int right_idx=-1;
		add_buffer(start_idx, stop_idx, left_idx, right_idx);
		vrb.bullet("Terminal window [" + stb.str(cidx) + "] -buffer:[" + chr + ":" + stb.str(positions_all_mb[start_idx]) + "-" + stb.str(positions_all_mb[stop_idx]) + "] / +buffer:[" + chr +  ":" + stb.str(positions_all_mb[left_idx]) + "-" + stb.str(positions_all_mb[right_idx]) + "] / L=" + stb.str(curr_window_cm_size) + "cM / L=" + stb.str(curr_window_mb_size) + "bp / C=" + stb.str(curr_window_count));
		cnk_info.add_chunk(positions_all_mb[left_idx], positions_all_mb[right_idx], positions_all_mb[start_idx], positions_all_mb[stop_idx], curr_window_cm_size, curr_window_mb_size, curr_window_count, curr_window_common_count);
		//fd << cidx << "\t" << chr << "\t" << chr << ":"<< positions_all_mb[left_idx] << "-" << positions_all_mb[right_idx] << "\t" << chr << ":" << positions_all_mb[start_idx] << "-" << positions_all_mb[stop_idx] << "\t" << curr_window_cm_size << "\t" << curr_window_mb_size << "\t" << curr_window_count <<"\t" << curr_window_common_count << std::endl;
		cidx ++;
	}
}

void chunker::split_sequential(output_file & fd, long int & cidx, std::string & chr, long int start_idx, long int stop_idx, const bool output_to_file)
{
	// Compute current window properties
	assert(stop_idx>start_idx);
	long int left_idx = -1; long int right_idx=-1;
	cnk_info.reset();

	for (long int i =start_idx; i<stop_idx; ++i)
	{
		long int curr_window_stop_idx=0;
		long int curr_window_count = 0;
		float curr_window_cm_size = 0;
		long int curr_window_mb_size = 0;
		long int curr_window_start_idx = i;
		long int curr_start_common_idx = all2common[curr_window_start_idx];
		left_idx = -1;
		right_idx = -1;

		if (curr_start_common_idx < (positions_common_mb.size() - window_count))
		{
			if (curr_start_common_idx + window_count < common2all.size())
				curr_window_stop_idx = common2all[curr_start_common_idx + window_count - 1];//condition 1 met
			else
				curr_window_stop_idx = positions_all_mb.size()-2;

		}
		else curr_window_stop_idx=positions_all_mb.size()-2;

		do {
			curr_window_stop_idx ++;
			curr_window_count = curr_window_stop_idx - curr_window_start_idx + 1;
			curr_window_mb_size = positions_all_mb[curr_window_stop_idx] - positions_all_mb[curr_window_start_idx];
			curr_window_cm_size = positions_all_cm[curr_window_stop_idx] - positions_all_cm[curr_window_start_idx];
		} while (((curr_window_stop_idx < (positions_all_mb.size() - 1)) && ((curr_window_cm_size < window_cm) || (curr_window_mb_size < window_mb) || curr_window_count < window_count)));
		if (positions_all_mb[curr_window_stop_idx]+buffer_mb >= positions_all_mb.back())
		{
			curr_window_stop_idx=positions_all_mb.size()-1;
			curr_window_count = curr_window_stop_idx - curr_window_start_idx + 1;
			curr_window_mb_size = positions_all_mb[curr_window_stop_idx] - positions_all_mb[curr_window_start_idx];
			curr_window_cm_size = positions_all_cm[curr_window_stop_idx] - positions_all_cm[curr_window_start_idx];
		}
		long int curr_stop_common_idx = all2common[curr_window_stop_idx];
		if (curr_stop_common_idx >= positions_common_mb.size()) curr_stop_common_idx=positions_common_mb.size()-1;
		long int curr_window_common_count = curr_stop_common_idx - curr_start_common_idx + 1;

		chunk_cm_length.push_back(curr_window_cm_size);
		chunk_mb_length.push_back(curr_window_mb_size);
		chunk_common_count.push_back(curr_window_common_count);

		add_buffer(curr_window_start_idx, curr_window_stop_idx, left_idx, right_idx);
		long int buf_start=positions_all_mb[left_idx];
		long int chk_start=positions_all_mb[curr_window_start_idx];
		long int buf_stop=positions_all_mb[right_idx];
		long int chk_stop=positions_all_mb[curr_window_stop_idx];

		if (whole_chr)
		{
			if (i == start_idx)
			{
				buf_start=1;
				chk_start=1;
			}
			if (curr_window_stop_idx == positions_all_mb.size()-1)
			{
				//buf_stop+=10000000;
				//chk_stop+=10000000;
				buf_stop=std::max(contig_len, buf_stop);
				chk_stop=std::max(contig_len, chk_stop);
			}
		}
		cnk_info.add_chunk(buf_start, buf_stop, chk_start, chk_stop, curr_window_cm_size, curr_window_mb_size, curr_window_count, curr_window_common_count);
		i = curr_window_stop_idx;

		cidx++;
	}
}

void chunker::add_buffer(const long int start_idx, const long int stop_idx, long int& left_idx, long int& right_idx)
{
	long int left_mb_size = -1, left_cm_size = -1, left_count = -1;
	if (start_idx > buffer_count) {
		left_idx = start_idx - buffer_count;
		do {
			left_idx --;
			left_count = start_idx - left_idx + 1;
			left_mb_size = positions_all_mb[start_idx] - positions_all_mb[left_idx];
			left_cm_size = positions_all_cm[start_idx] - positions_all_cm[left_idx];
		} while (((left_idx > 0) && ((left_cm_size < buffer_cm) || (left_mb_size < buffer_mb) || left_count < buffer_count)));
	} else { left_idx = 0; }

	long int right_mb_size = -1, right_cm_size = -1, right_count = -1;
	if (stop_idx < (positions_all_mb.size() - buffer_count)) {
		right_idx = stop_idx + buffer_count - 1;
		do {
			right_idx ++;
			right_count = right_idx - stop_idx + 1;
			right_mb_size = positions_all_mb[right_idx] - positions_all_mb[stop_idx];
			right_cm_size = positions_all_cm[right_idx] - positions_all_cm[stop_idx];
		} while (((right_idx < (positions_all_mb.size() - 1)) && ((right_cm_size < buffer_cm) || (right_mb_size < buffer_mb) || right_count < buffer_count)));
	} else { right_idx = positions_all_mb.size() - 1; }
}

void chunker::chunk() {
	//Initialize
	rng.setSeed(options["seed"].as < long int > ());
	sparse_maf = options["sparse-maf"].as < float > ();

	window_cm = options["window-cm"].as < float > ();
	window_mb = (long int) (options["window-mb"].as < float > ()*1e6);
	window_count = options["window-count"].as < long int > ();
	buffer_cm = options["buffer-cm"].as < float > ();
	buffer_mb = (long int) (options["buffer-mb"].as < float > ()*1e6);
	buffer_count = options["buffer-count"].as < long int > ();

	buildCoordinates();

	//Read input data (overlapping coordinates)
	readData(options["input"].as < std::string > (), options["region"].as < std::string > (), options["threads"].as <long int> ());

	if (options.count("map"))
	{
		gmap_reader readerGM;
		readerGM.readGeneticMapFile(options["map"].as < std::string > ());
		setGeneticMap(readerGM, map_positions_all, positions_all_mb, positions_all_cm);
	} else setGeneticMap(map_positions_all, positions_all_mb, positions_all_cm);
	vrb.bullet("Region spans " + stb.str(positions_all_mb.back() - positions_all_mb[0] + 1) + " bp and " + stb.str(positions_all_cm.back() - positions_all_cm[0] + 1, 2) + " cM");
	if (positions_all_mb.size() == 0) vrb.error("No markers in region: " + options["region"].as < std::string > () + ". Check chromosome name and start/end positions.");

	cnk_info.chr = chrID;

	long int cidx = 0;
	long int cidx_uniform=0;

	vrb.title("Splitting data long into chunks and writing to [" + options["output"].as < std::string > () + "]");
	output_file fd(options["output"].as < std::string > ());

	if (options.count("recursive"))
	{
		split_recursive(fd, cidx, chrID, 0, positions_all_mb.size() - 1);
		cnk_info.output_to_file(fd);
	}
	else if (options.count("sequential"))
	{
		split_sequential(fd, cidx, chrID, 0, positions_all_mb.size() - 1, false);
		bool reg_small = chunk_cm_length.size()==1;

		long int i=0;
		while ((chunk_cm_length.back() < window_cm || chunk_mb_length.back() < window_mb) && i<10000 && !reg_small)
		{
			chunk_cm_length.clear();
			chunk_mb_length.clear();
			chunk_common_count.clear();
			window_count+=100;
			cidx=0;
			vrb.progress("Solution failed. Trying again: [" + stb.str(window_cm) + " cM | " + stb.str(window_mb/1e6) + " Mb | " + stb.str(window_count) + " common variants]" + "\tAttempt " + stb.str(i+1) + "/10000.", i/10000);
			split_sequential(fd, cidx, chrID, 0, positions_all_mb.size() - 1, false);
			++i;
		}

		if (reg_small)
		{
			chunk_cm_length.clear();
			chunk_mb_length.clear();
			chunk_common_count.clear();
			cidx=0;
			vrb.bullet("Region appears to small to find a sequential solution (only one chunk detected). Check your parameters if this is not a desired behavior.");
			split_recursive(fd, cidx, chrID, 0, positions_all_mb.size() - 1);
			cnk_info.output_to_file(fd);
		}
		else if (i==10000)
		{
			chunk_cm_length.clear();
			chunk_mb_length.clear();
			chunk_common_count.clear();
			cidx=0;
			vrb.bullet("Could not find a sequential solution to the problem. Running recursive algorithm.");
			split_recursive(fd, cidx, chrID, 0, positions_all_mb.size() - 1);
			cnk_info.output_to_file(fd);
		}
		else
		{
			vrb.bullet("Solution found ["+ stb.str(i) + "/10000]. Writing chunks to file.");
			//split_sequential(fd, cidx, chrID, 0, positions_all_mb.size() - 1, true);
			cnk_info.output_to_file(fd);

			if (options.count("uniform-number-variants"))
			{
				fd.close();
				output_file fd(options["output"].as < std::string > () + "_uniform");
				for (long int j=0; j<chunk_common_count.size(); ++j) if (chunk_common_count[j] > window_count) window_count=chunk_common_count[j];
				cidx_uniform=0;
				i=0;

				split_sequential(fd, cidx_uniform, chrID, 0, positions_all_mb.size() - 1, false);

				while ((chunk_cm_length.back() < window_cm || chunk_mb_length.back() < window_mb) && i<10000)
				{
					chunk_cm_length.clear();
					chunk_mb_length.clear();
					chunk_common_count.clear();
					window_count+=100;
					cidx_uniform=0;
					vrb.progress("First solution failed. Trying again: [" + stb.str(window_cm) + " cM | " + stb.str(window_mb/1e6) + " Mb | " + stb.str(window_count) + " common variants]" + "\tAttempt " + stb.str(i+1) + "/10000.", i/10000);
					split_sequential(fd, cidx_uniform, chrID, 0, positions_all_mb.size() - 1, false);
					++i;
				}
				if (i==10000) vrb.bullet("Could not find a sequential solution to the problem for the uniform pass.");
				else
				{
					vrb.bullet("Uniform solution found ["+ stb.str(i) + "/10000]. Writing chunks to file.");
					cnk_info.output_to_file(fd);
					//split_sequential(fd, cidx_uniform, chrID, 0, positions_all_mb.size() - 1, true);
				}
			}
		}
	}

	vrb.bullet("#chunks = " + stb.str(cidx));
	if (options.count("uniform-number-variants")) vrb.bullet("#chunks uniform solution = " + stb.str(cidx_uniform));

	fd.close();

	//Finalize
	vrb.title("Total running time = " + stb.str(tac.abs_time()) + " seconds");
}
