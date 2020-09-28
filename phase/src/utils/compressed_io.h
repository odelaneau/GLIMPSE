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

#ifndef _COMPRESSED_IO_H
#define _COMPRESSED_IO_H

//STL INCLUDES
#include <iostream>
#include <sstream>
#include <fstream>

//BOOST INCLUDES
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

class input_file : public boost::iostreams::filtering_istream {
protected:
	std::ifstream file_descriptor;

public:
	input_file(const std::string filename) {
		if (filename.substr(filename.find_last_of(".") + 1) == "gz") {
			file_descriptor.open(filename.c_str(), std::ios::in | std::ios::binary);
			push(boost::iostreams::gzip_decompressor());
		} else if (filename.substr(filename.find_last_of(".") + 1) == "bz2") {
			file_descriptor.open(filename.c_str(), std::ios::in | std::ios::binary);
			push(boost::iostreams::bzip2_decompressor());
		} else file_descriptor.open(filename.c_str());
		if (!file_descriptor.fail()) push(file_descriptor);
	}

	~input_file() {
		close();
	}

	bool fail() {
		return file_descriptor.fail();
	}

	void close() {
		if (!file_descriptor.fail()) {
			if (!empty()) reset();
			file_descriptor.close();
		}
	}
};

class output_file : public boost::iostreams::filtering_ostream {
protected:
	std::ofstream file_descriptor;

public:
	output_file(std::string filename) {
		if (filename.substr(filename.find_last_of(".") + 1) == "gz") {
			file_descriptor.open(filename.c_str(), std::ios::out | std::ios::binary);
			push(boost::iostreams::gzip_compressor());
		} else if (filename.substr(filename.find_last_of(".") + 1) == "bz2") {
			file_descriptor.open(filename.c_str(), std::ios::out | std::ios::binary);
			push(boost::iostreams::bzip2_compressor());
		} else file_descriptor.open(filename.c_str());
		if (!file_descriptor.fail()) push(file_descriptor);
	}

	~output_file() {
		close();
	}

	bool fail() {
		return file_descriptor.fail();
	}

	void close() {
		if (!file_descriptor.fail()) {
			if (!empty()) reset();
			file_descriptor.close();
		}
	}
};

#endif
