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

#ifndef _CHUNKER_H
#define _CHUNKER_H

#include <utils/otools.h>
#include <containers/variant_map.h>

class chunker {
public:
	//COMMAND LINE OPTIONS
	bpo::options_description descriptions;
	bpo::variables_map options;

	//INTERNAL DATA
	//variant_map V;

	float sparse_maf;
	std::string chrID;
	std::vector < float > positions_all_cm;
	std::vector < int > positions_all_mb;

	std::vector < float > positions_common_cm;
	std::vector < int > positions_common_mb;

	std::vector < int > common2all;

	std::multimap < int, int > map_positions_all;	//associative container of variant with position in bp
	std::multimap < int, int > map_positions_common;	//associative container of variant with position in bp

	//std::vector < int > counterR;


	//PARAMETERS
	float window_cm;
	float window_mb;
	int window_count;
	float buffer_cm;
	float buffer_mb;
	int buffer_count;

	bool report_common_variants;

	//CONSTRUCTOR
	chunker();
	~chunker();

	//METHODS
	void readData(std::string fmain, std::string region, int nthreads);

	//PARAMETERS
	void declare_options();
	void parse_command_line(std::vector < std::string > &);
	void check_options();
	void verbose_options();
	void verbose_files();

	//FILE I/O
	void split(output_file &, int &, std::string &, int, int);
	void chunk();
	void chunk(std::vector < std::string > &);

	std::vector < int > getByPos(const int, const std::multimap < int, int >& map_positions);
	void setGeneticMap(const gmap_reader&, const std::multimap < int, int >&, const std::vector<int>& positions_mb, std::vector<float>& positions_cm);
	void setGeneticMap(const std::multimap < int, int >&, const std::vector<int>& positions_mb, std::vector<float>& positions_cm);
	int setCentiMorgan(const std::vector < int > & pos_bp, const std::vector < double > & pos_cM, const std::multimap < int, int >& map_positions, std::vector<float> & positions_cm);
	int interpolateCentiMorgan(const std::vector < int > & pos_bp, const std::vector < double > & pos_cM, const std::multimap < int, int >&, const std::vector<int> & positions_mb, std::vector<float>& positions_cm);
	unsigned int length() const;
	double lengthcM() const;

};


#endif


