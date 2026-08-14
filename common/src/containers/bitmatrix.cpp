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

#include <containers/bitmatrix.h>
#include <sys/mman.h>

bitmatrix::bitmatrix() {
	n_rows = 0;
	n_cols = 0;
	n_bytes = 0;
	bytes = nullptr;
	owns_mmap = false;
	mmap_length = 0;
	skip_and_discard_on_load = false;
}

bitmatrix::~bitmatrix() {
	release();
}

void bitmatrix::release() {
	if (bytes) {
		if (owns_mmap) ::munmap(bytes, mmap_length);
		else std::free(bytes);
	}
	bytes = nullptr;
	owns_mmap = false;
	mmap_length = 0;
	n_bytes = 0;
	n_cols = 0;
	n_rows = 0;
}

void bitmatrix::adopt_mmap(unsigned char * mapped_bytes, unsigned long mapped_length, unsigned long nrow, unsigned long ncol, unsigned long nbytes_logical) {
	release();		//tear down any previous state first (retry-safety)
	bytes = mapped_bytes;
	owns_mmap = true;
	mmap_length = mapped_length;
	n_rows = nrow;
	n_cols = ncol;
	n_bytes = nbytes_logical;
}

unsigned char * bitmatrix::rowPtr(unsigned int row) {
	return &bytes[((unsigned long)row) * (n_cols/8)];
}

const unsigned char * bitmatrix::rowPtr(unsigned int row) const {
	return &bytes[((unsigned long)row) * (n_cols/8)];
}

void bitmatrix::subset(const bitmatrix & BM, std::vector < unsigned int > rows) {
	assert(!owns_mmap);
	n_rows = rows.size() + ((rows.size()%8)?(8-(rows.size()%8)):0);
	n_cols = BM.n_cols;
	unsigned long new_n_bytes = (n_cols/8) * (unsigned long)n_rows;
	if (new_n_bytes > n_bytes) bytes = (unsigned char*)realloc(bytes, new_n_bytes*sizeof(unsigned char));
	n_bytes = new_n_bytes;
	unsigned long offset_addr = 0;
	for (int r = 0 ; r < rows.size() ; r ++) {
		std::memcpy(&bytes[offset_addr], &BM.bytes[((unsigned long)rows[r]) * (BM.n_cols/8)], n_cols/8);
		offset_addr += n_cols/8;
	}
}

/*
void bitmatrix::subset(bitmatrix & BM, vector < int > rows, unsigned int col_from, unsigned int col_to) {
	n_rows = rows.size() + ((rows.size()%8)?(8-(rows.size()%8)):0);
	unsigned long row_start = col_from/8;
	unsigned long row_end = col_to/8;
	unsigned long n_bytes_per_row = row_end - row_start + 1;
	n_cols = n_bytes_per_row * 8;
	n_bytes = n_bytes_per_row * n_rows;
	cout << row_start << " " << row_end << " " << n_bytes_per_row << endl;
	bytes = (unsigned char*)realloc(bytes, n_bytes*sizeof(unsigned char));
	cout << "OKAY!" << n_bytes_per_row << endl;
	unsigned long offset_addr = 0;
	for (int r = 0 ; r < rows.size() ; r ++) {
		row_start = ((unsigned long)rows[r]) * (BM.n_cols/8) + col_from/8;
		row_end = ((unsigned long)rows[r]) * (BM.n_cols/8) + col_to/8;
		memcpy(&bytes[offset_addr], &BM.bytes[row_start], n_bytes_per_row);
		offset_addr += n_bytes_per_row;
	}
	//return col_from % 8;
}
*/


void bitmatrix::allocate(unsigned int nrow, unsigned int ncol) {
	assert(!owns_mmap && "allocate() called on an mmap-backed bitmatrix");
	n_rows = nrow + ((nrow%8)?(8-(nrow%8)):0);
	n_cols = ncol + ((ncol%8)?(8-(ncol%8)):0);
	n_bytes = (n_cols/8) * (unsigned long)n_rows;
	bytes = (unsigned char*)std::malloc(n_bytes*sizeof(unsigned char));
	std::memset(bytes, 0, n_bytes);
}

void bitmatrix::reallocate(unsigned int nrow, unsigned int ncol) {
	assert(!owns_mmap && "reallocate() called on an mmap-backed bitmatrix");
	n_rows = nrow + ((nrow%8)?(8-(nrow%8)):0);
	n_cols = ncol + ((ncol%8)?(8-(ncol%8)):0);
	unsigned long int new_n_bytes = (n_cols/8) * (unsigned long)n_rows;
	if (new_n_bytes > n_bytes) bytes = (unsigned char*)std::realloc(bytes, new_n_bytes*sizeof(unsigned char));
	n_bytes = new_n_bytes;
}


/*
 * This algorithm for transposing bit matrices is adapted from the code of Timur Kristóf
 * Timur Kristóf: https://github.com/venemo
 * Original version of the code (MIT license): https://github.com/Venemo/fecmagic/blob/master/src/binarymatrix.h
 * Of note, function abracadabra is the same than getMultiplyUpperPart function in the original code from Timur Kristóf.
 */
void bitmatrix::transpose(bitmatrix & BM, unsigned int _max_row, unsigned int _max_col) {
	assert(!BM.owns_mmap && "transpose() destination is an mmap-backed bitmatrix");
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

void bitmatrix::transpose(bitmatrix & BM, unsigned int _min_row, unsigned int _min_col, unsigned int _max_row, unsigned int _max_col) {
	assert(!BM.owns_mmap && "transpose() destination is an mmap-backed bitmatrix");
	unsigned int min_row = _min_row - ((_min_row%8)?(8-(_min_row%8)):0);
	unsigned int min_col = _min_col - ((_min_col%8)?(8-(_min_col%8)):0);
	unsigned int max_row = _max_row + ((_max_row%8)?(8-(_max_row%8)):0);
	unsigned int max_col = _max_col + ((_max_col%8)?(8-(_max_col%8)):0);

	unsigned long targetAddr, sourceAddr;
	union { unsigned int x[2]; unsigned char b[8]; } m4x8d;
	for (unsigned int row = min_row; row < max_row; row += 8) {
		for (unsigned int col = min_col; col < max_col; col += 8) {
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

void bitmatrix::transpose(bitmatrix & BM) {
	transpose(BM, n_rows, n_cols);
}
