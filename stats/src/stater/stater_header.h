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

#ifndef _STATER_H
#define _STATER_H

#include <utils/otools.h>

class stater {
public:
	//COMMAND LINE OPTIONS
	bpo::options_description descriptions;
	bpo::variables_map options;

	//INTERNAL DATA
	vector < string > var_chr;
	vector < int > var_pos;
	vector < string > var_ref;
	vector < string > var_alt;
	vector < string > var_ids;
	vector < double > var_freq;
	vector < double > var_info;
	vector < double > var_avgp;

	vector < string > spl_ids;
	vector < double > spl_info;
	vector < double > spl_avgp;

	//CONSTRUCTOR
	stater();
	~stater();

	//METHODS
	void readDataAndComputeStats(vector < string > files, vector < string > regions);

	//PARAMETERS
	void declare_options();
	void parse_command_line(vector < string > &);
	void check_options();
	void verbose_options();
	void verbose_files();

	//FILE I/O
	void stats();
	void stats(vector < string > &);
};


#endif


