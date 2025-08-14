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
#include "checksum_utils.h"

/**
 * @class stats2D
 * @brief Calculates online statistics for two-dimensional data, including means,
 * variances, standard deviations, covariance, and correlation.
 */
class stats2D {
protected:
	unsigned long int m_n;    ///< Number of samples processed

	double m_oldMx, m_newMx; ///< Old and new means for X
	double m_oldMy, m_newMy; ///< Old and new means for Y

	double m_oldSx, m_newSx; ///< Old and new sums of squares for X
	double m_oldSy, m_newSy; ///< Old and new sums of squares for Y

	double m_oldC, m_newC;   ///< Old and new covariance sums

public:

	/**
	 * @brief Default constructor initializes statistics.
	 */
	stats2D() {
		m_n = 0;
		m_oldMx = 0; m_newMx = 0;
		m_oldMy = 0; m_newMy = 0;
		m_oldSx = 0; m_newSx = 0;
		m_oldSy = 0; m_newSy = 0;
		m_oldC = 0; m_newC = 0;
	}

	/**
	 * @brief Reset all statistics to zero.
	 */
	void clear() {
		m_n = 0;
		m_oldMx = 0; m_newMx = 0;
		m_oldMy = 0; m_newMy = 0;
		m_oldSx = 0; m_newSx = 0;
		m_oldSy = 0; m_newSy = 0;
		m_oldC = 0; m_newC = 0;
	}

	/**
	 * @brief Add a new sample (x,y) to the statistics.
	 * @tparam T Numeric type of input data
	 * @param x Value for X variable
	 * @param y Value for Y variable
	 */
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

	/**
	 * @return Number of samples processed
	 */
	unsigned long int size() const {
		return m_n;
	}

	/**
	 * @return Mean of X values
	 */
	double meanX() const {
		return (m_n > 0) ? m_newMx : 0.0;
	}

	/**
	 * @return Mean of Y values
	 */
	double meanY() const {
		return (m_n > 0) ? m_newMy : 0.0;
	}

	/**
	 * @return Sample variance of X
	 */
	double varX() const {
		return ( (m_n > 1) ? m_newSx/(m_n-1) : 0.0 );
	}

	/**
	 * @return Sample variance of Y
	 */
	double varY() const {
		return ( (m_n > 1) ? m_newSy/(m_n-1) : 0.0 );
	}

	/**
	 * @return Sample standard deviation of X
	 */
	double sdX() const {
		return sqrt( varX() );
	}

	/**
	 * @return Sample standard deviation of Y
	 */
	double sdY() const {
		return sqrt( varY() );
	}

	/**
	 * @return Whether standard deviation of X is not available or zero
	 */
	bool sdNAx() const {
		return !(m_n > 1 && m_newSx >0);
	}

	/**
	 * @return Whether standard deviation of Y is not available or zero
	 */
	bool sdNAy() const {
		return !(m_n > 1 && m_newSy >0);
	}

	/**
	 * @return Pearson correlation coefficient between X and Y
	 */
	double corrXY() const {
		const double sdx = (m_n > 1 && m_newSx >0) ? sdX() : 1e-10;
		const double sdy = (m_n > 1 && m_newSy >0) ? sdY() : 1e-10;

		return ( (m_n > 0) ? (m_newC/((m_n-1)*sdx*sdy)) : 0.0 );
	}

	/**
	 * @return Square of Pearson correlation coefficient (R²)
	 */
	double corrXY_square() const {
		const double sdx = (m_n > 1 && m_newSx >0) ? sdX() : 1e-10;
		const double sdy = (m_n > 1 && m_newSy >0) ? sdY() : 1e-10;
		const double r = (m_n > 0) ? (m_newC/((m_n-1)*sdx*sdy)) : 0.0;
		return ( r*r );
	}

	/**
	 * @brief Computes approximate standard error for 95% confidence interval of R².
	 * @return Standard error for 95% confidence interval of R²
	 */
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


/**
 * @class stats1D
 * @brief Calculates online statistics for one-dimensional data, including mean,
 * variance, and standard deviation.
 *
 * Algorithm taken from:
 * https://www.johndcook.com/blog/standard_deviation/
 * and https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
 */
class stats1D {
protected:
	unsigned long int m_n;    ///< Number of samples processed
	double m_oldM;            ///< Previous mean
	double m_newM;            ///< Current mean
	double m_oldS;            ///< Previous sum of squares
	double m_newS;            ///< Current sum of squares

public:
	/**
	 * @brief Default constructor initializes statistics.
	 */
	stats1D() {
		m_n = 0;
		m_oldM = 0;
		m_newM = 0;
		m_oldS = 0;
		m_newS = 0;
	}

	/**
	 * @brief Constructor that computes stats from vector input.
	 * @tparam T Numeric type of input vector elements
	 * @param X Vector of values
	 */
	template <class T>
	stats1D(std::vector < T > & X) {
		m_n = 0;
		m_oldM = 0;
		m_newM = 0;
		m_oldS = 0;
		m_newS = 0;
		for (uint32_t e = 0 ; e < X.size() ; e ++) push(X[e]);
	}

	/**
	 * @brief Reset all statistics to zero.
	 */
	void clear() {
		m_n = 0;
		m_oldM = 0;
		m_newM = 0;
		m_oldS = 0;
		m_newS = 0;
	}

	/**
	 * @brief Add a new sample x to the statistics.
	 * @tparam T Numeric type of input data
	 * @param x Value to add
	 */
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

	/**
	 * @return Number of samples processed
	 */
	unsigned long int size() const {
		return m_n;
	}

	/**
	 * @return Mean of the samples
	 */
	double mean() const {
		return (m_n > 0) ? m_newM : 0.0;
	}

	/**
	 * @return Sample variance
	 */
	double variance() const {
		return ( (m_n > 1) ? m_newS/(m_n - 1) : 0.0 );
	}

	/**
	 * @return Sample standard deviation
	 */
	double sd() const {
		return sqrt( variance() );
	}

	friend class boost::serialization::access;

	/**
	 * @brief Serialization function for checkpointing
	 */
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & m_n;
		ar & m_oldM;
		ar & m_newM;
		ar & m_oldS;
		ar & m_newS;
	}

	/**
	 * @brief Update checksum with internal state data
	 */
	void update_checksum(checksum &crc) const
	{
		crc.process_data(m_n);
		crc.process_data(m_oldM);
		crc.process_data(m_newM);
		crc.process_data(m_oldS);
		crc.process_data(m_newS);
	}
};

#endif
