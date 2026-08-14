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
#include <utils/checksum_utils.h>
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/array.hpp"


inline static unsigned int abracadabra(const unsigned int &i1, const unsigned int &i2) {
	return static_cast<unsigned int>((static_cast<unsigned long int>(i1) * static_cast<unsigned long int>(i2)) >> 32);
}

class bitmatrix
{
public:
	unsigned long int n_bytes, n_cols, n_rows;
	unsigned char * bytes;

	//MMAP SUPPORT
	//When owns_mmap is true, `bytes` points into an mmap()'d region (adopted via
	//adopt_mmap(), set up by the caller doing the actual open()/mmap()/madvise()
	//syscalls) rather than a malloc'd buffer; release()/the destructor must munmap()
	//it instead of free()'ing it. Such an instance is logically read-only for its
	//whole lifetime -- see the asserts on the mutating methods below.
	bool owns_mmap;
	unsigned long mmap_length;
	//When true, serialize()'s load branch consumes (but never materializes) the
	//bytes for this member -- used when a separate mmap'd cache is going to be
	//adopted afterward instead. See serialize() below.
	bool skip_and_discard_on_load;

	bitmatrix();
	virtual ~bitmatrix();

	//A bitmatrix owns a raw resource (malloc'd buffer or, now, an mmap()'d region);
	//copying it by value would double-free/double-munmap on destruction of both
	//copies. No code path relies on copying a bitmatrix by value.
	bitmatrix(const bitmatrix &) = delete;
	bitmatrix & operator=(const bitmatrix &) = delete;

	void subset(const bitmatrix & BM, std::vector < unsigned int > rows);
	void allocate(unsigned int nrow, unsigned int ncol);
	void reallocate(unsigned int nrow, unsigned int ncol);
	void set(unsigned int row, unsigned int col, unsigned char bit);
	void set(unsigned int row, unsigned char bit);
	unsigned char get(unsigned int row, unsigned int col) const;
	unsigned char getByte(unsigned int row, unsigned int col) const;
	void setByte(unsigned int row, unsigned int col, unsigned char byte);
	unsigned char * rowPtr(unsigned int row);
	const unsigned char * rowPtr(unsigned int row) const;

	//Adopt an externally-created mmap() mapping (the caller performs the actual
	//open()/mmap()/madvise() syscalls; this just records ownership so release()/the
	//destructor tears it down correctly, and sets up the logical row/col view).
	//nbytes_logical may be smaller than mapped_length (mmap rounds length up to a
	//page) -- get()/getByte()/rowPtr() addressing uses n_cols/n_rows, not mmap_length.
	void adopt_mmap(unsigned char * mapped_bytes, unsigned long mapped_length, unsigned long nrow, unsigned long ncol, unsigned long nbytes_logical);
	//Explicit teardown (munmap or free, depending on owns_mmap), callable before a
	//fresh load (e.g. a retry loop reusing the same object). Also invoked by ~bitmatrix.
	void release();

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
			if (skip_and_discard_on_load)
			{
				//Consume exactly n_bytes from the archive stream via a small reusable
				//scratch buffer, so peak RSS never includes this member's full byte
				//array. Used when a separate (e.g. mmap'd) representation of this data
				//is going to be adopted afterward -- see bitmatrix::adopt_mmap() and
				//caller::read_binary_reference_panel().
				constexpr unsigned long CHUNK = 1u << 20; //1 MB
				std::vector<unsigned char> scratch(std::min<unsigned long>(n_bytes, CHUNK));
				for (unsigned long done = 0 ; done < n_bytes ; )
				{
					unsigned long take = std::min<unsigned long>(n_bytes - done, CHUNK);
					ar & boost::serialization::make_array<unsigned char>(scratch.data(), take);
					done += take;
				}
				bytes = nullptr;
				return;
			}
			assert(bytes == nullptr);
			bytes = (unsigned char*)std::malloc(n_bytes*sizeof(unsigned char));
		}
		ar & boost::serialization::make_array<unsigned char>(bytes, n_bytes);
	}

	void update_checksum(checksum &crc) const
	{
		crc.process_data(n_bytes);
		crc.process_data(n_cols);
		crc.process_data(n_rows);
		//NB: if this instance is mmap-backed (owns_mmap), this still works correctly
		//but pages in the entire mapped region to compute the checksum -- combining
		//--checkpoint-file-in/out with the mmap'd hvar-cache path loses the memory
		//benefit for this object. Not addressed here; rare/orthogonal combination.
		crc.process_data(bytes, n_bytes*sizeof(unsigned char));
	}
};

inline
void bitmatrix::set(unsigned int row, unsigned int col, unsigned char bit) {
	assert(!owns_mmap && "set() called on an mmap-backed bitmatrix");
	unsigned int bitcol = col % 8;
	unsigned long targetAddr = ((unsigned long)row) * (n_cols/8) + col/8;
	unsigned char mask = ~(1 << (7 - bitcol));
	this->bytes[targetAddr] &= mask;
	this->bytes[targetAddr] |= (bit << (7 - bitcol));
}

inline
void bitmatrix::set(unsigned int row, unsigned char bit) {
	assert(!owns_mmap && "set() called on an mmap-backed bitmatrix");
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

//Writes 8 packed bits at column col (must be a multiple of 8). Avoids the
//read-modify-write of set() when all 8 bits in the byte are known up front.
inline
void bitmatrix::setByte(unsigned int row, unsigned int col, unsigned char byte) {
	assert((col & 7) == 0);
	assert(!owns_mmap && "setByte() called on an mmap-backed bitmatrix");
	bytes[((unsigned long)row) * (n_cols>>3) +  (col>>3)] = byte;
}

#endif
