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
#include "otools.h"
#include "checksum_utils.h"
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/array.hpp"

/**
 * @brief Multiply two unsigned integers and take the upper 32 bits of the 64-bit result.
 * 
 * This function performs a 32-bit multiplication and returns the high 32 bits,
 * which can be useful in certain hashing or arithmetic contexts.
 * 
 * @param i1 First unsigned integer operand.
 * @param i2 Second unsigned integer operand.
 * @return Upper 32 bits of the 64-bit product of i1 and i2.
 */
inline static unsigned int abracadabra(const unsigned int &i1, const unsigned int &i2) {
	return static_cast<unsigned int>((static_cast<unsigned long int>(i1) * static_cast<unsigned long int>(i2)) >> 32);
}

/**
 * @class bitmatrix
 * @brief Efficient bit matrix container storing bits in a contiguous byte array.
 * 
 * This class manages a 2D bit matrix of size n_rows by n_cols,
 * stored in row-major order, 8 bits per byte.
 * It provides operations for setting and getting bits, transposing,
 * and serializing the matrix.
 */
class bitmatrix
{
public:
	/** Total number of bytes allocated for bit storage */
	unsigned long int n_bytes;

	/** Number of columns (bits per row) */
	unsigned long int n_cols;

	/** Number of rows */
	unsigned long int n_rows;

	/** Pointer to raw byte data storing bits */
	unsigned char * bytes;

	/** Default constructor initializes empty matrix */
	bitmatrix();

	/** Virtual destructor frees allocated memory */
	virtual ~bitmatrix();

	/**
	 * @brief Copy a subset of rows from another bitmatrix into this one.
	 * 
	 * The rows vector specifies the indices of rows in BM to copy.
	 * The resulting matrix is resized to fit the subset with rows padded
	 * to multiple of 8 for alignment.
	 * 
	 * @param BM Source bitmatrix to copy rows from.
	 * @param rows Vector of row indices to include in the subset.
	 */
	void subset(bitmatrix & BM, std::vector < unsigned int > rows);

	/**
	 * @brief Allocate memory for a bitmatrix with given rows and columns.
	 * 
	 * Rows and columns are padded to multiples of 8 to align with byte boundaries.
	 * All bits are initialized to zero.
	 * 
	 * @param nrow Number of rows.
	 * @param ncol Number of columns.
	 */
	void allocate(unsigned int nrow, unsigned int ncol);

	/**
	 * @brief Reallocate memory to resize the bitmatrix to new dimensions.
	 * 
	 * Rows and columns are padded to multiples of 8.
	 * The allocated byte array may be expanded if needed.
	 * 
	 * @param nrow New number of rows.
	 * @param ncol New number of columns.
	 */
	void reallocate(unsigned int nrow, unsigned int ncol);

	/**
	 * @brief Set a single bit in the matrix at (row, col).
	 * 
	 * @param row Row index.
	 * @param col Column index.
	 * @param bit Bit value to set (0 or 1).
	 */
	void set(unsigned int row, unsigned int col, unsigned char bit);

	/**
	 * @brief Set all bits in a given row to the same bit value.
	 * 
	 * @param row Row index.
	 * @param bit Bit value to set (0 or 1).
	 */
	void set(unsigned int row, unsigned char bit);

	/**
	 * @brief Get the bit value at position (row, col).
	 * 
	 * @param row Row index.
	 * @param col Column index.
	 * @return Bit value (0 or 1).
	 */
	unsigned char get(unsigned int row, unsigned int col) const;

	/**
	 * @brief Get the byte containing bits for a given row and column byte index.
	 * 
	 * @param row Row index.
	 * @param col Byte index within the row (column / 8).
	 * @return Byte value containing 8 bits.
	 */
	unsigned char getByte(unsigned int row, unsigned int col) const;

	/**
	 * @brief Transpose a subregion of the bitmatrix from BM into this bitmatrix.
	 * 
	 * The region is defined by min/max row and column indices.
	 * The output matrix will have rows and columns swapped in that region.
	 * 
	 * @param BM Source bitmatrix.
	 * @param _min_row Minimum row index in source.
	 * @param _min_col Minimum column index in source.
	 * @param _max_row Maximum row index in source.
	 * @param _max_col Maximum column index in source.
	 */
	void transpose(bitmatrix & BM, unsigned int _min_row, unsigned int _min_col, unsigned int _max_row, unsigned int _max_col);

	/**
	 * @brief Transpose a subregion [0:_max_row, 0:_max_col] from BM into this bitmatrix.
	 * 
	 * @param BM Source bitmatrix.
	 * @param _max_row Maximum row index.
	 * @param _max_col Maximum column index.
	 */
	void transpose(bitmatrix & BM, unsigned int _max_row, unsigned int _max_col);

	/**
	 * @brief Transpose the entire bitmatrix BM into this bitmatrix.
	 * 
	 * @param BM Source bitmatrix.
	 */
	void transpose(bitmatrix & BM);

	/// Boost Serialization access
	friend class boost::serialization::access;

	/**
	 * @brief Serialize or deserialize the bitmatrix object.
	 * 
	 * During loading, allocates memory for bytes array.
	 * 
	 * @tparam Archive Serialization archive type.
	 * @param ar Archive instance.
	 * @param version Serialization version number.
	 */
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

	/**
	 * @brief Update a checksum with the content of the bitmatrix for data integrity.
	 * 
	 * @param crc Checksum object to update.
	 */
	void update_checksum(checksum &crc) const
	{
		crc.process_data(n_bytes);
		crc.process_data(n_cols);
		crc.process_data(n_rows);
		crc.process_data(bytes, n_bytes*sizeof(unsigned char));
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
