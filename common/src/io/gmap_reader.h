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

#ifndef _GMAP_DATA_H
#define _GMAP_DATA_H

#include "otools.h"

/**
 * @class gmap_reader
 * @brief Class to manage genetic map data with base pair positions and corresponding centiMorgan positions.
 * 
 * This class stores vectors of base pair positions and their corresponding genetic map positions in centiMorgans,
 * and provides functionality to read genetic map files.
 */
class gmap_reader {
public:
	/// Vector of base pair positions (bp)
	std::vector<long int> pos_bp;

	/// Vector of centiMorgan positions (cM) corresponding to pos_bp
	std::vector<double> pos_cm;

	/// Constructor
	gmap_reader();

	/// Destructor
	~gmap_reader();

	/**
	 * @brief Read genetic map data from a file.
	 * 
	 * Parses a genetic map file and fills pos_bp and pos_cm vectors accordingly.
	 * 
	 * @param filename Path to the genetic map file.
	 */
	void readGeneticMapFile(const std::string filename);
};

#endif
