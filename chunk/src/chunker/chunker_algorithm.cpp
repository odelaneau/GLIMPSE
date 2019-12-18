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
#include <chunker/chunker_header.h>

void chunker::split(output_file & fd, int & cidx, string & chr, vector < int > & vec) {
	int midpoint = vec.size() / 2;
	vector < int > vec1 = vector < int > (vec.begin() , vec.begin() + midpoint);
	vector < int > vec2 = vector < int > (vec.begin() + midpoint, vec.end());

	int size1 = vec1.back() - vec1[0];
	int size2 = vec2.back() - vec2[0];

	if (size1 < chunk_size || size2 < chunk_size) {
		int output_start = vec[0];
		int output_stop = vec.back();
		int input_start = output_start - buffer_size;
		if (input_start < 0) input_start = 0;
		int input_stop = output_stop + buffer_size;
		fd << chr << ":"<< input_start << "-" << input_stop << "\t" << chr << ":"<< output_start << "-" << output_stop << "\t" << input_stop - input_start << "\t" << output_stop - output_start << "\t" << vec.size() << "\t" << cidx << endl;
		cidx ++;
	} else {
		split (fd, cidx, chr, vec1);
		split (fd, cidx, chr, vec2);
	}
}



void chunker::chunk() {
	//Initialize
	rng.setSeed(options["seed"].as < int > ());
	chunk_size = options["window"].as < int > ();
	buffer_size = options["buffer"].as < int > ();

	vector < string > tmp0 = options["input"].as < vector < string > > ();
	vector < string > tmp1 = options["reference"].as < vector < string > > ();
	readData(tmp0, tmp1);


	//
	int cidx = 0;
	output_file fd(options["output"].as < string > ());
	vrb.title("Splitting data into chunks");
	for (int c = 0 ; c < C.size() ; c++) {
		split(fd, cidx, C[c], V[c]);
	}
	vrb.bullet("#chunks = " + stb.str(cidx));
	fd.close();

	//Finalize
	vrb.title("Total running time = " + stb.str(tac.abs_time()) + " seconds");
}
