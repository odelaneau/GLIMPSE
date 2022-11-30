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

#ifndef _RANDOM_NUMBER_H
#define _RANDOM_NUMBER_H

#include <cfloat>
#include <cstdint>
#include <random>
#include <vector>


class random_number_generator {
protected:
	unsigned int seed;
	std::uniform_int_distribution < unsigned int > uniformDistributionInt;
	std::uniform_real_distribution < double > uniformDistributionDouble;
	std::uniform_real_distribution < float > uniformDistributionFloat;

public:

	std::mt19937 randomEngine;

	random_number_generator(unsigned int _seed = 15052011) : seed(_seed), randomEngine(_seed), uniformDistributionInt(0, 32768), uniformDistributionDouble(0, 1.0), uniformDistributionFloat(0.0f,1.0f) {
	}

	~random_number_generator(){
	}

	void setSeed(unsigned int _seed) {
		seed = _seed;
		randomEngine.seed(seed);
	}

	unsigned int getSeed() {
		return seed;
	}

	std::mt19937 & getEngine() {
		return randomEngine;
	}

	unsigned int getInt(unsigned int imin, unsigned int imax) {
		return uniformDistributionInt(randomEngine, std::uniform_int_distribution < unsigned int > {imin, imax}.param());
	}

	unsigned int getInt(unsigned int isize) {
		return getInt(0, isize - 1);
	}

	double getFloat(float fmin=0.0f, float fmax=1.0f) {
		return uniformDistributionFloat(randomEngine, std::uniform_real_distribution < float > {fmin, fmax}.param());
	}

	double getDouble(double fmin, double fmax) {
		return uniformDistributionDouble(randomEngine, std::uniform_real_distribution < double > {fmin, fmax}.param());
	}

	double getDouble() {
		return getDouble(0.0, 1.0);
	}

	bool flipCoin() {
		return (getDouble() < 0.5);
	}

	int sample(std::vector < float > & vec, float sum) {
		float csum = vec[0];
		float u = getFloat() * sum;
		for (int i = 0; i < vec.size() - 1; ++i) {
			if ( u <= csum ) return i;
			csum += vec[i+1];
		}
		return vec.size() - 1;
	}

	int sample(std::vector < double > & vec, double sum) {
		double csum = vec[0];
		double u = getDouble() * sum;
		for (int i = 0; i < vec.size() - 1; ++i) {
			if ( u < csum ) return i;
			csum += vec[i+1];
		}
		return vec.size() - 1;
	}

	int sample4(const double * vec, double sum) {
		double csum = vec[0];
		double u = getDouble() * sum;
		for (int i = 0; i < 3; ++i) {
			if ( u < csum ) return i;
			csum += vec[i+1];
		}
		return 3;
	}

	double dpois(int n, double lambda)
	{
		if (n<0) return 0.0;
		else if (n == 0) return 1/exp(lambda);
		else return exp(n*log(lambda)-lambda-lgamma(n+1));
	}
};

#endif

