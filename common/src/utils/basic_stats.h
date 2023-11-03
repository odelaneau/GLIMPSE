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

#ifndef _BASIC_STATS_H
#define _BASIC_STATS_H

#include <vector>
#include <utils/checksum_utils.h>

class stats2D {
protected:
	unsigned long int m_n;
	double m_oldMx, m_newMx;
	double m_oldSx, m_newSx;
	double m_oldMy, m_newMy;
	double m_oldSy, m_newSy;
	double m_oldC, m_newC;

public:

	stats2D() {
		m_n = 0;
		m_oldMx = 0; m_newMx = 0;
		m_oldMy = 0; m_newMy = 0;
		m_oldSx = 0; m_newSx = 0;
		m_oldSy = 0; m_newSy = 0;
		m_oldC = 0; m_newC = 0;
	}

	void clear() {
		m_n = 0;
		m_oldMx = 0; m_newMx = 0;
		m_oldMy = 0; m_newMy = 0;
		m_oldSx = 0; m_newSx = 0;
		m_oldSy = 0; m_newSy = 0;
		m_oldC = 0; m_newC = 0;
	}

	template <class T>
	void push(T x, T y) {
		m_n++;
		if (m_n == 1) {
			m_oldMx = m_newMx = x;
			m_oldMy = m_newMy = y;
			m_oldSx = m_oldSy = 0.0;
			m_oldC = m_newC = 0.0;
		} else {
			m_newMy = m_oldMy + (y - m_oldMy)/m_n;
			m_newC = m_oldC + (x - m_oldMx) * (y - m_newMy);
			m_newMx = m_oldMx + (x - m_oldMx)/m_n;
			m_newSx = m_oldSx + (x - m_oldMx)*(x - m_newMx);
			m_newSy = m_oldSy + (y - m_oldMy)*(y - m_newMy);

			m_oldC = m_newC;
			m_oldMx = m_newMx;
            m_oldSx = m_newSx;
            m_oldMy = m_newMy;
            m_oldSy = m_newSy;
		}
	}

	unsigned long int size() const {
		return m_n;
	}

	double meanX() const {
		return (m_n > 0) ? m_newMx : 0.0;
	}

	double meanY() const {
		return (m_n > 0) ? m_newMy : 0.0;
	}

	double varX() const {
		return ( (m_n > 1) ? m_newSx/(m_n-1) : 0.0 );
	}

	double varY() const {
		return ( (m_n > 1) ? m_newSy/(m_n-1) : 0.0 );
	}

	double sdX() const {
		return sqrt( varX() );
	}

	double sdY() const {
		return sqrt( varY() );
	}

	bool sdNAx() const {
		return !(m_n > 1 && m_newSx >0);
	}

	bool sdNAy() const {
		return !(m_n > 1 && m_newSy >0);
	}

	double corrXY() const {
		const double sdx = (m_n > 1 && m_newSx >0) ? sdX() : 1e-10;
		const double sdy = (m_n > 1 && m_newSy >0) ? sdY() : 1e-10;

		return ( (m_n > 0) ? (m_newC/((m_n-1)*sdx*sdy)) : 0.0 );
	}

	double corrXY_square() const {
		const double sdx = (m_n > 1 && m_newSx >0) ? sdX() : 1e-10;
		const double sdy = (m_n > 1 && m_newSy >0) ? sdY() : 1e-10;
		const double r = (m_n > 0) ? (m_newC/((m_n-1)*sdx*sdy)) : 0.0;
		return ( r*r );
	}

	double getSE95() const {
		const double k = 2.0;
		const double R2 = corrXY_square();
		const double SE = (m_n > 2) ? sqrt((4*R2*((1-R2)*(1-R2))*((m_n-k-1)*(m_n-k-1)))/((m_n*m_n-1)*(m_n + 3))) : 0.0;
		return SE;
		//CI == 0.95:
        //upper = R2 + 2*SE
        //lower = R2 - 2*SE
	}
};

// Code taken from there: https://www.johndcook.com/blog/standard_deviation/ and there: https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
class stats1D {
protected:
	unsigned long int m_n;
	double m_oldM;
	double m_newM;
	double m_oldS;
	double m_newS;

public:
	stats1D() {
		m_n = 0;
		m_oldM = 0;
		m_newM = 0;
		m_oldS = 0;
		m_newS = 0;
	}

	template <class T>
	stats1D(std::vector < T > & X) {
		m_n = 0;
		m_oldM = 0;
		m_newM = 0;
		m_oldS = 0;
		m_newS = 0;
		for (uint32_t e = 0 ; e < X.size() ; e ++) push(X[e]);
	}

	void clear() {
		m_n = 0;
		m_oldM = 0;
		m_newM = 0;
		m_oldS = 0;
		m_newS = 0;
	}

	template <class T>
	void push(T x) {
		m_n++;
		if (m_n == 1) {
			m_oldM = m_newM = x;
			m_oldS = 0.0;
		} else {
			m_newM = m_oldM + (x - m_oldM)/m_n;
            m_newS = m_oldS + (x - m_oldM)*(x - m_newM);
            m_oldM = m_newM;
            m_oldS = m_newS;
		}
	}

	unsigned long int size() const {
		return m_n;
	}

	double mean() const {
		return (m_n > 0) ? m_newM : 0.0;
	}

	double variance() const {
		return ( (m_n > 1) ? m_newS/(m_n - 1) : 0.0 );
	}

	double sd() const {
		return sqrt( variance() );
	}

	void update_checksum(checksum &crc)
	{
		crc.process_data(m_n);
		crc.process_data(m_oldM);
		crc.process_data(m_newM);
		crc.process_data(m_oldS);
		crc.process_data(m_newS);
	}
};

#endif
