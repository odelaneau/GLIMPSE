/*******************************************************************************
 * @file random_number.h
 * @brief Random number generator utility class providing various distributions and sampling methods.
 *
 * This class wraps the C++11 <random> library to offer uniform integer and floating-point
 * random numbers, coin flips, and weighted sampling from vectors.
 *
 * @copyright Copyright (C) 2022-2023 Simone Rubinacci and Olivier Delaneau
 * @license MIT License
 ******************************************************************************/

#ifndef _RANDOM_NUMBER_H
#define _RANDOM_NUMBER_H

#include <cfloat>
#include <cstdint>
#include <random>
#include <vector>

/**
 * @class random_number_generator
 * @brief Utility class for generating random numbers and sampling.
 *
 * Provides methods to:
 * - Generate uniformly distributed integers and floating-point numbers.
 * - Flip a coin (true/false with equal probability).
 * - Sample an index from a weighted distribution given by a vector.
 * - Compute Poisson probability mass function.
 */
class random_number_generator {
protected:
	/** Seed used for random number engine */
	unsigned int seed;

	/** Uniform integer distribution object */
	std::uniform_int_distribution<unsigned int> uniformDistributionInt;

	/** Uniform double precision floating-point distribution */
	std::uniform_real_distribution<double> uniformDistributionDouble;

	/** Uniform single precision floating-point distribution */
	std::uniform_real_distribution<float> uniformDistributionFloat;

public:
	/** Mersenne Twister random engine */
	std::mt19937 randomEngine;

	/**
	 * @brief Constructor initializes random engine with given seed.
	 * @param _seed Initial seed value (default 15052011).
	 */
	random_number_generator(unsigned int _seed = 15052011) :
		seed(_seed),
		randomEngine(_seed),
		uniformDistributionInt(0, 32768),
		uniformDistributionDouble(0, 1.0),
		uniformDistributionFloat(0.0f, 1.0f) {
	}

	/** Destructor */
	~random_number_generator() {}

	/**
	 * @brief Set a new seed for the random engine.
	 * @param _seed New seed value.
	 */
	void setSeed(unsigned int _seed) {
		seed = _seed;
		randomEngine.seed(seed);
	}

	/**
	 * @brief Get current seed.
	 * @return Current seed value.
	 */
	unsigned int getSeed() {
		return seed;
	}

	/**
	 * @brief Access the underlying random engine.
	 * @return Reference to the std::mt19937 engine.
	 */
	std::mt19937& getEngine() {
		return randomEngine;
	}

	/**
	 * @brief Get a random integer in [imin, imax].
	 * @param imin Minimum integer (inclusive).
	 * @param imax Maximum integer (inclusive).
	 * @return Random integer between imin and imax.
	 */
	unsigned int getInt(unsigned int imin, unsigned int imax) {
		return uniformDistributionInt(randomEngine, std::uniform_int_distribution<unsigned int>{imin, imax}.param());
	}

	/**
	 * @brief Get a random integer in [0, isize-1].
	 * @param isize Upper bound exclusive.
	 * @return Random integer between 0 and isize-1.
	 */
	unsigned int getInt(unsigned int isize) {
		return getInt(0, isize - 1);
	}

	/**
	 * @brief Get a random float in [fmin, fmax].
	 * @param fmin Minimum float (default 0.0f).
	 * @param fmax Maximum float (default 1.0f).
	 * @return Random float between fmin and fmax.
	 */
	double getFloat(float fmin = 0.0f, float fmax = 1.0f) {
		return uniformDistributionFloat(randomEngine, std::uniform_real_distribution<float>{fmin, fmax}.param());
	}

	/**
	 * @brief Get a random double in [fmin, fmax].
	 * @param fmin Minimum double.
	 * @param fmax Maximum double.
	 * @return Random double between fmin and fmax.
	 */
	double getDouble(double fmin, double fmax) {
		return uniformDistributionDouble(randomEngine, std::uniform_real_distribution<double>{fmin, fmax}.param());
	}

	/**
	 * @brief Get a random double in [0.0, 1.0].
	 * @return Random double between 0 and 1.
	 */
	double getDouble() {
		return getDouble(0.0, 1.0);
	}

	/**
	 * @brief Simulate a coin flip (true with probability 0.5).
	 * @return true or false with equal probability.
	 */
	bool flipCoin() {
		return (getDouble() < 0.5);
	}

	/**
	 * @brief Sample an index from a weighted vector of floats.
	 * @param vec Vector of weights.
	 * @param sum Sum of all weights in vec.
	 * @return Selected index based on weighted probability.
	 */
	int sample(std::vector<float>& vec, float sum) {
		float csum = vec[0];
		float u = getFloat() * sum;
		for (int i = 0; i < (int)vec.size() - 1; ++i) {
			if (u <= csum) return i;
			csum += vec[i + 1];
		}
		return (int)vec.size() - 1;
	}

	/**
	 * @brief Sample an index from a weighted vector of doubles.
	 * @param vec Vector of weights.
	 * @param sum Sum of all weights in vec.
	 * @return Selected index based on weighted probability.
	 */
	int sample(std::vector<double>& vec, double sum) {
		double csum = vec[0];
		double u = getDouble() * sum;
		for (int i = 0; i < (int)vec.size() - 1; ++i) {
			if (u < csum) return i;
			csum += vec[i + 1];
		}
		return (int)vec.size() - 1;
	}

	/**
	 * @brief Sample an index from a weighted array of 4 doubles.
	 * @param vec Pointer to array of 4 weights.
	 * @param sum Sum of all weights in vec.
	 * @return Selected index in [0..3] based on weighted probability.
	 */
	int sample4(const double* vec, double sum) {
		double csum = vec[0];
		double u = getDouble() * sum;
		for (int i = 0; i < 3; ++i) {
			if (u < csum) return i;
			csum += vec[i + 1];
		}
		return 3;
	}

	/**
	 * @brief Compute the Poisson probability mass function P(X = n) for parameter lambda.
	 * @param n Number of occurrences (non-negative integer).
	 * @param lambda Expected number of occurrences (rate parameter).
	 * @return Probability of observing exactly n occurrences.
	 */
	double dpois(int n, double lambda) {
		if (n < 0) return 0.0;
		else if (n == 0) return 1 / exp(lambda);
		else return exp(n * log(lambda) - lambda - lgamma(n + 1));
	}
};

#endif
