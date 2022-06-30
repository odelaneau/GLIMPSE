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

#ifndef _BASIC_ALGOS_H
#define _BASIC_ALGOS_H

#include <vector>

class basic_algos {
public:
	basic_algos () {};
	~basic_algos () {};

	template < class T >
	unsigned int imax(std::vector < T > & vec) {
		T maxValue = vec[0];
		int maxIndex = 0;
		for (unsigned int i = 1; i < vec.size() ; i ++) if (vec[i] > maxValue) {
			maxValue = vec[i];
			maxIndex = i;
		}
		return maxIndex;
	}

	template<class T>
	inline const T& clamp( const T& v, const T& lo, const T& hi ) const
	{
	    return std::max(lo, std::min(v, hi));
	}
};

#endif

