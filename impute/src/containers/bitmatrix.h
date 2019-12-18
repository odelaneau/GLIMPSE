/*******************************************************************************
 * Copyright (C) 2018 Olivier Delaneau, University of Lausanne
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

#include <utils/otools.h>

inline static unsigned int abracadabra(const unsigned int &i1, const unsigned int &i2) {
	return static_cast<unsigned int>((static_cast<unsigned long int>(i1) * static_cast<unsigned long int>(i2)) >> 32);
}

class bitmatrix	{
public:
	unsigned long int n_bytes, n_cols, n_rows;
	unsigned char * bytes;

	bitmatrix() {
		n_rows = 0;
		n_cols = 0;
		n_bytes = 0;
		bytes = NULL;
	}

	void allocate(unsigned int nrow, unsigned int ncol) {
		n_rows = nrow + ((nrow%8)?(8-(nrow%8)):0);
		n_cols = ncol + ((ncol%8)?(8-(ncol%8)):0);
		n_bytes = (n_cols/8) * (unsigned long)n_rows;
		bytes = (unsigned char*)malloc(n_bytes*sizeof(unsigned char));
		std::memset(bytes, 0, n_bytes);
	}

	~bitmatrix() {
		n_bytes = 0;
		if (bytes != NULL) free(bytes);
	}

	void set(unsigned int row, unsigned int col, unsigned char bit);
	unsigned char get(unsigned int row, unsigned int col);


	/*
	 * This algorithm for transposing bit matrices is adapted from the code of Timur Kristóf
	 * Timur Kristóf: https://github.com/venemo
	 * Original version of the code (MIT license): https://github.com/Venemo/fecmagic/blob/master/src/binarymatrix.h
	 * Of note, function abracadabra is the same than getMultiplyUpperPart function in the original code from Timur Kristóf.
	 */
	void transpose(bitmatrix & BM, unsigned int _max_row, unsigned int _max_col) {
		unsigned int max_row = _max_row + ((_max_row%8)?(8-(_max_row%8)):0);
		unsigned int max_col = _max_col + ((_max_col%8)?(8-(_max_col%8)):0);
		unsigned long targetAddr, sourceAddr;
		union { unsigned int x[2]; unsigned char b[8]; } m4x8d;
		for (unsigned int row = 0; row < max_row; row += 8) {
			for (unsigned int col = 0; col < max_col; col += 8) {
				for (unsigned int i = 0; i < 8; i++) {
					sourceAddr = (row+i) * ((unsigned long)(n_cols/8)) + col/8;
					m4x8d.b[7 - i] = this->bytes[sourceAddr];
				}
				for (unsigned int i = 0; i < 7; i++) {
					targetAddr = ((col+i) * ((unsigned long)(n_rows/8)) + (row) / 8);
					BM.bytes[targetAddr]  = static_cast<unsigned char>(abracadabra(m4x8d.x[1] & (0x80808080 >> i), (0x02040810 << i)) & 0x0f) << 4;
	                BM.bytes[targetAddr] |= static_cast<unsigned char>(abracadabra(m4x8d.x[0] & (0x80808080 >> i), (0x02040810 << i)) & 0x0f) << 0;
				}
				targetAddr = ((col+7) * ((unsigned long)(n_rows/8)) + (row) / 8);
				BM.bytes[targetAddr]  = static_cast<unsigned char>(abracadabra((m4x8d.x[1] << 7) & (0x80808080 >> 0), (0x02040810 << 0)) & 0x0f) << 4;
	            BM.bytes[targetAddr] |= static_cast<unsigned char>(abracadabra((m4x8d.x[0] << 7) & (0x80808080 >> 0), (0x02040810 << 0)) & 0x0f) << 0;
			}
		}
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
unsigned char bitmatrix::get(unsigned int row, unsigned int col) {
	unsigned int bitcol = col  % 8;
	unsigned long targetAddr = ((unsigned long)row) * (n_cols/8) +  col/8;
	unsigned char result = (this->bytes[targetAddr] >> (7 - bitcol)) & 1;
	return result;
}

#endif
