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

#include <io/gmap_reader.h>

gmap_reader::gmap_reader() {
}

gmap_reader::~gmap_reader() {
	std::vector < int > ().swap(pos_bp);
	std::vector < double > ().swap(pos_cm);
}

void gmap_reader::readGeneticMapFile(const std::string fmap) {
	tac.clock();
	std::string buffer;
	std::vector < std::string > tokens;
	int line = 0;
	int line_used = 0;
	input_file fd_gmap(fmap);
	if (fd_gmap.fail()) vrb.error("Cannot open genetic map file");
	getline(fd_gmap, buffer, '\n');
	int prev_bp = 0;
	double prev_cm = 0;
	while (getline(fd_gmap, buffer, '\n')) {
		if (stb.split(buffer, tokens) == 3) {
			int curr_bp = atoi(tokens[0].c_str());
			double curr_cm = atof(tokens[2].c_str());
			if (curr_bp < prev_bp || curr_cm < prev_cm)
				vrb.error("Wrong order in your genetic map file " + stb.str(prev_bp) + "bp / " + stb.str(prev_cm,5) + "cM > " + stb.str(curr_bp) + "bp / " + stb.str(curr_cm,5) + "cM");
			pos_bp.push_back(curr_bp);
			pos_cm.push_back(curr_cm);
			prev_bp = curr_bp;
			prev_cm = curr_cm;
			++line_used;
		} else vrb.error("Parsing line " + stb.str(line) + " : incorrect number of columns, observed: " + stb.str(tokens.size()) + " expected: 3");
		line++;
	}
	if (line_used < 2) vrb.error("Error parsing genetic file: not enough lines detected. GLIMPSE accepts only files a header line and at least two map entries. Please make sure you have at least three lines in your map file.");
	if (line_used != line) vrb.warning("There are lines in the map file that have been excluded due to wrong formatting (not three columns?). Please check your map file carefully.");
	fd_gmap.close();
	vrb.bullet("Genetic map          : [n=" + stb.str(line_used) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
