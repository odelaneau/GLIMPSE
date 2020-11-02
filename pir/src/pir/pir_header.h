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

#ifndef _PIR_H
#define _PIR_H

#include <utils/otools.h>

class allele {
public:
	int index;
	int qual;
	int base;

	allele(int _index, int _qual, int _base) {
		index = _index;
		qual = _qual;
		base = _base;
	}
};


class pir {
public:
	//COMMAND LINE OPTIONS
	bpo::options_description descriptions;
	bpo::variables_map options;

	//SITEs
	int snp_index;
	vector < int > positions;
	vector < string > refA;
	vector < string > altA;

	//PIRs
	vector < vector < allele > > dataPIR;		// vector of reads, each being a vector of bases
	unordered_map < string, int  > mapPIR;	// hashmap connecting read ID to their data

	//PARAMs
	int min_baseq;
	int min_mapq;

	//CONSTRUCTOR
	pir();
	~pir();

	//ROUTINES
	char getBase(int code);

	//PARAMETERS
	void declare_options();
	void parse_command_line(vector < string > &);
	void check_options();
	void verbose_options();
	void verbose_files();

	//READ FILES
	void readPositions (string fref, string region);
	void readSequences(string fbam, string region);
	void parseBam(void * d);

	//WRITE FILES
	void writePIRs(string);

	//FILE I/O
	void read_files_and_initialise();
	void extract(vector < string > &);
	void extract();
	void write_files_and_finalise();
};

inline
char pir::getBase(int code) {
	switch (code) {
	case 1: return 'A';
	case 2: return 'C';
	case 4: return 'G';
	case 8: return 'T';
	case 15: return 'N';
	}
	return 0;
}

#endif


