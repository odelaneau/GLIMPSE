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

#include <ligater/ligater_header.h>

void ligater::read_files_and_initialise() {
	vrb.title("Initialization:");

	//step0: Initialize seed & other
	rng.setSeed(options["seed"].as < int > ());

	//step1:read filenames
	string buffer;
	string filelist = options["input"].as < string > ();
	vrb.title("Read filenames in [" + filelist + "]");
	input_file fd(filelist);
	while (getline(fd, buffer)) filenames.push_back(buffer);
	vrb.bullet("#files = " + stb.str(filenames.size()));
	if (filenames.size() == 0) vrb.error("No filenames in input file.");

	//step2: initilize flags
	vrb.title("Initilialize flags");
	nfiles = filenames.size();
	current_stages = vector < unsigned char > (nfiles, STAGE_NONE);
	vrb.bullet("done");
}
