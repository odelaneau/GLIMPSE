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

#ifndef _BITMATRIX_H
#define _BITMATRIX_H

#include <cstdlib>
#include <utils/otools.h>
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/array.hpp"
#include <zlib.h>


inline static unsigned int abracadabra(const unsigned int &i1, const unsigned int &i2) {
	return static_cast<unsigned int>((static_cast<unsigned long int>(i1) * static_cast<unsigned long int>(i2)) >> 32);
}

class bitmatrix
{
public:
	unsigned long int n_bytes, n_cols, n_rows;
	unsigned char * bytes;

	bitmatrix();
	bitmatrix(std::vector<bool> vec); 
	virtual ~bitmatrix();


	void subset(bitmatrix & BM, std::vector < unsigned int > rows);
	void allocate(unsigned int nrow, unsigned int ncol);
	void reallocate(unsigned int nrow, unsigned int ncol);
	void set(unsigned int row, unsigned int col, unsigned char bit);
	void set(unsigned int row, unsigned char bit);
	unsigned char get(unsigned int row, unsigned int col) const;
	unsigned char getByte(unsigned int row, unsigned int col) const;

	void transpose(bitmatrix & BM, unsigned int _min_row, unsigned int _min_col, unsigned int _max_row, unsigned int _max_col);
	void transpose(bitmatrix & BM, unsigned int _max_row, unsigned int _max_col);
	void transpose(bitmatrix & BM);

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & n_bytes;
		ar & n_cols;
		ar & n_rows;

		if (Archive::is_loading::value)
		{
			assert(bytes == nullptr);
			bytes = (unsigned char*)std::malloc(n_bytes*sizeof(unsigned char));
		}
		ar & boost::serialization::make_array<unsigned char>(bytes, n_bytes);
	}

	void update_checksum(checksum &crc)
	{
		crc.process_data(n_bytes);
		crc.process_data(n_cols);
		crc.process_data(n_rows);
		crc.process_data(bytes, n_bytes);
	}
};

inline
void bitmatrix::set(unsigned int row, unsigned int col, unsigned char bit) {
	unsigned int bitcol = col % 8;
	unsigned long targetAddr = ((unsigned long)row) * (n_cols/8) + col/8;
	unsigned char mask = ~(1 << (7 - bitcol));
	this->bytes[targetAddr] &= mask;
	this->bytes[targetAddr] |= (bit << (7 - bitcol));
}

inline
void bitmatrix::set(unsigned int row, unsigned char bit) {
	std::memset(&bytes[(unsigned long)row * (n_cols/8)], bit * 255, n_cols/8);
}

inline
unsigned char bitmatrix::get(unsigned int row, unsigned int col) const {
	unsigned long targetAddr = ((unsigned long)row) * (n_cols>>3) +  (col>>3);
	return (this->bytes[targetAddr] >> (7 - (col%8))) & 1;
}

inline
unsigned char bitmatrix::getByte(unsigned int row, unsigned int col) const {
	return bytes[((unsigned long)row) * (n_cols>>3) +  (col>>3)];
}

#endif
