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

#ifndef _GENOTYPE_BAM_CALLER_H
#define _GENOTYPE_BAM_CALLER_H

#include "otools.h"
#include "variant_map.h"
#include "containers/haplotype_set.h"
#include "containers/glimpse_mpileup.h"

/**
 * @namespace calling_utils
 * @brief Utility constants for genotype likelihood and error model calculations.
 *
 * This namespace contains precomputed lookup tables used for genotype calling likelihood 
 * calculations and error modeling in variant calling.
 * 
 * Each array contains 256 elements indexed by quality score or similar metric.
 * These tables avoid recalculating expensive logarithmic or probability values repeatedly.
 *
 * Arrays include:
 * - unphred: Conversion of phred quality scores to error probabilities (log-scale)
 * - log10_one_minus_err: Log10 of (1 - error probability) for each quality score
 * - log10_err_div_three: Log10 of (error probability / 3), assuming errors are equally distributed across 3 alternate alleles
 * - log10_one_div_two_minus_err_div_three: Precomputed values used in diploid genotype likelihoods involving error terms
 * - log10_err: Log10 error probability for each quality score
 * - log10_one_div_two_minus_err: Other precomputed terms related to genotype likelihood calculations
 *
 * These tables improve efficiency by replacing runtime computations with array lookups.
 */
namespace calling_utils
{
	//equivalent to: for (int i=0; i<256; i++) unphred[i] = pow(10., -i/10.);//element zero is set to one.
	const static float unphred[256] = { 1.000000e+00, 7.943282e-01, 6.309573e-01, 5.011872e-01, 3.981072e-01, 3.162278e-01, 2.511886e-01, 1.995262e-01, 1.584893e-01, 1.258925e-01, 1.000000e-01, 7.943282e-02, 6.309573e-02, 5.011872e-02, 3.981072e-02, 3.162278e-02, 2.511886e-02, 1.995262e-02, 1.584893e-02, 1.258925e-02, 1.000000e-02, 7.943282e-03, 6.309573e-03, 5.011872e-03, 3.981072e-03, 3.162278e-03, 2.511886e-03, 1.995262e-03, 1.584893e-03, 1.258925e-03, 1.000000e-03, 7.943282e-04, 6.309573e-04, 5.011872e-04, 3.981072e-04, 3.162278e-04, 2.511886e-04, 1.995262e-04, 1.584893e-04, 1.258925e-04, 1.000000e-04, 7.943282e-05, 6.309573e-05, 5.011872e-05, 3.981072e-05, 3.162278e-05, 2.511886e-05, 1.995262e-05, 1.584893e-05, 1.258925e-05, 1.000000e-05, 7.943282e-06, 6.309573e-06, 5.011872e-06, 3.981072e-06, 3.162278e-06, 2.511886e-06, 1.995262e-06, 1.584893e-06, 1.258925e-06, 1.000000e-06, 7.943282e-07, 6.309573e-07, 5.011872e-07, 3.981072e-07, 3.162278e-07, 2.511886e-07, 1.995262e-07, 1.584893e-07, 1.258925e-07, 1.000000e-07, 7.943282e-08, 6.309573e-08, 5.011872e-08, 3.981072e-08, 3.162278e-08, 2.511886e-08, 1.995262e-08, 1.584893e-08, 1.258925e-08, 1.000000e-08, 7.943282e-09, 6.309573e-09, 5.011872e-09, 3.981072e-09, 3.162278e-09, 2.511886e-09, 1.995262e-09, 1.584893e-09, 1.258925e-09, 1.000000e-09, 7.943282e-10, 6.309573e-10, 5.011872e-10, 3.981072e-10, 3.162278e-10, 2.511886e-10, 1.995262e-10, 1.584893e-10, 1.258925e-10, 1.000000e-10, 7.943282e-11, 6.309573e-11, 5.011872e-11, 3.981072e-11, 3.162278e-11, 2.511886e-11, 1.995262e-11, 1.584893e-11, 1.258925e-11, 1.000000e-11, 7.943282e-12, 6.309573e-12, 5.011872e-12, 3.981072e-12, 3.162278e-12, 2.511886e-12, 1.995262e-12, 1.584893e-12, 1.258925e-12, 1.000000e-12, 7.943282e-13, 6.309573e-13, 5.011872e-13, 3.981072e-13, 3.162278e-13, 2.511886e-13, 1.995262e-13, 1.584893e-13, 1.258925e-13, 1.000000e-13, 7.943282e-14, 6.309573e-14, 5.011872e-14, 3.981072e-14, 3.162278e-14, 2.511886e-14, 1.995262e-14, 1.584893e-14, 1.258925e-14, 1.000000e-14, 7.943282e-15, 6.309573e-15, 5.011872e-15, 3.981072e-15, 3.162278e-15, 2.511886e-15, 1.995262e-15, 1.584893e-15, 1.258925e-15, 1.000000e-15, 7.943282e-16, 6.309573e-16, 5.011872e-16, 3.981072e-16, 3.162278e-16, 2.511886e-16, 1.995262e-16, 1.584893e-16, 1.258925e-16, 1.000000e-16, 7.943282e-17, 6.309573e-17, 5.011872e-17, 3.981072e-17, 3.162278e-17, 2.511886e-17, 1.995262e-17, 1.584893e-17, 1.258925e-17, 1.000000e-17, 7.943282e-18, 6.309573e-18, 5.011872e-18, 3.981072e-18, 3.162278e-18, 2.511886e-18, 1.995262e-18, 1.584893e-18, 1.258925e-18, 1.000000e-18, 7.943282e-19, 6.309573e-19, 5.011872e-19, 3.981072e-19, 3.162278e-19, 2.511886e-19, 1.995262e-19, 1.584893e-19, 1.258925e-19, 1.000000e-19, 7.943282e-20, 6.309573e-20, 5.011872e-20, 3.981072e-20, 3.162278e-20, 2.511886e-20, 1.995262e-20, 1.584893e-20, 1.258925e-20, 1.000000e-20, 7.943282e-21, 6.309573e-21, 5.011872e-21, 3.981072e-21, 3.162278e-21, 2.511886e-21, 1.995262e-21, 1.584893e-21, 1.258925e-21, 1.000000e-21, 7.943282e-22, 6.309573e-22, 5.011872e-22, 3.981072e-22, 3.162278e-22, 2.511886e-22, 1.995262e-22, 1.584893e-22, 1.258925e-22, 1.000000e-22, 7.943282e-23, 6.309573e-23, 5.011872e-23, 3.981072e-23, 3.162278e-23, 2.511886e-23, 1.995262e-23, 1.584893e-23, 1.258925e-23, 1.000000e-23, 7.943282e-24, 6.309573e-24, 5.011872e-24, 3.981072e-24, 3.162278e-24, 2.511886e-24, 1.995262e-24, 1.584893e-24, 1.258925e-24, 1.000000e-24, 7.943282e-25, 6.309573e-25, 5.011872e-25, 3.981072e-25, 3.162278e-25, 2.511886e-25, 1.995262e-25, 1.584893e-25, 1.258925e-25, 1.000000e-25, 7.943282e-26, 6.309573e-26, 5.011872e-26, 3.981072e-26, 3.162278e-26};
	const static float log10_one_minus_err[256] = {0.99999 , -0.686825 , -0.432923 , -0.302062 , -0.220481 , -0.165089 , -0.125628 , -0.0966529 , -0.0749404 , -0.0584352 , -0.0457575 , -0.0359445 , -0.0283048 , -0.0223307 , -0.0176431 , -0.0139554 , -0.0110483 , -0.00875293 , -0.00693823 , -0.00550215 , -0.00436481 , -0.0034635 , -0.00274889 , -0.0021821 , -0.00173241 , -0.00137554 , -0.00109227 , -0.000867397 , -0.000688856 , -0.000547089 , -0.000434512 , -0.000345109 , -0.000274108 , -0.000217717 , -0.00017293 , -0.000137358 , -0.000109104 , -8.66618e-05 , -6.88365e-05 , -5.46779e-05 , -4.34316e-05 , -3.44986e-05 , -2.7403e-05 , -2.17668e-05 , -1.72899e-05 , -1.37338e-05 , -1.09091e-05 , -8.6654e-06 , -6.88316e-06 , -5.46748e-06 , -4.34297e-06 , -3.44974e-06 , -2.74022e-06 , -2.17663e-06 , -1.72896e-06 , -1.37336e-06 , -1.0909e-06 , -8.66532e-07 , -6.88311e-07 , -5.46745e-07 , -4.34295e-07 , -3.44972e-07 , -2.74021e-07 , -2.17663e-07 , -1.72896e-07 , -1.37336e-07 , -1.0909e-07 , -8.66531e-08 , -6.8831e-08 , -5.46744e-08 , -4.34295e-08 , -3.44972e-08 , -2.74021e-08 , -2.17663e-08 , -1.72896e-08 , -1.37336e-08 , -1.0909e-08 , -8.66531e-09 , -6.8831e-09 , -5.46744e-09 , -4.34294e-09 , -3.44972e-09 , -2.74021e-09 , -2.17663e-09 , -1.72896e-09 , -1.37336e-09 , -1.0909e-09 , -8.66531e-10 , -6.8831e-10 , -5.46744e-10 , -4.34294e-10 , -3.44972e-10 , -2.74021e-10 , -2.17663e-10 , -1.72896e-10 , -1.37336e-10 , -1.0909e-10 , -8.66531e-11 , -6.8831e-11 , -5.46744e-11 , -4.34295e-11 , -3.44972e-11 , -2.74021e-11 , -2.17663e-11 , -1.72896e-11 , -1.37336e-11 , -1.0909e-11 , -8.6653e-12 , -6.88308e-12 , -5.46745e-12 , -4.34295e-12 , -3.44974e-12 , -2.74023e-12 , -2.17663e-12 , -1.72894e-12 , -1.37335e-12 , -1.0909e-12 , -8.66545e-13 , -6.88289e-13 , -5.46725e-13 , -4.34285e-13 , -3.44988e-13 , -2.74014e-13 , -2.17649e-13 , -1.72904e-13 , -1.3732e-13 , -1.09114e-13 , -8.66448e-14 , -6.8853e-14 , -5.46774e-14 , -4.3443e-14 , -3.44747e-14 , -2.73869e-14 , -2.17456e-14 , -1.73097e-14 , -1.37417e-14 , -1.08969e-14 , -8.67895e-15 , -6.89494e-15 , -5.44845e-15 , -4.33947e-15 , -3.47158e-15 , -2.74833e-15 , -2.16974e-15 , -1.73579e-15 , -1.35006e-15 , -1.10898e-15 , -8.67895e-16 , -6.75029e-16 , -5.3038e-16 , -4.33947e-16 , -3.37515e-16 , -2.89298e-16 , -2.41082e-16 , -1.92865e-16 , -1.44649e-16 , -9.64327e-17 , -9.64327e-17 , -4.82164e-17 , -4.82164e-17 , -4.82164e-17 , -4.82164e-17 , -4.82164e-17 , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17 , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17  , -4.82164e-17 };
	const static float log10_err_div_three[256] = {-0.477121 , -0.577121 , -0.677121 , -0.777121 , -0.877121 , -0.977121 , -1.07712 , -1.17712 , -1.27712 , -1.37712 , -1.47712 , -1.57712 , -1.67712 , -1.77712 , -1.87712 , -1.97712 , -2.07712 , -2.17712 , -2.27712 , -2.37712 , -2.47712 , -2.57712 , -2.67712 , -2.77712 , -2.87712 , -2.97712 , -3.07712 , -3.17712 , -3.27712 , -3.37712 , -3.47712 , -3.57712 , -3.67712 , -3.77712 , -3.87712 , -3.97712 , -4.07712 , -4.17712 , -4.27712 , -4.37712 , -4.47712 , -4.57712 , -4.67712 , -4.77712 , -4.87712 , -4.97712 , -5.07712 , -5.17712 , -5.27712 , -5.37712 , -5.47712 , -5.57712 , -5.67712 , -5.77712 , -5.87712 , -5.97712 , -6.07712 , -6.17712 , -6.27712 , -6.37712 , -6.47712 , -6.57712 , -6.67712 , -6.77712 , -6.87712 , -6.97712 , -7.07712 , -7.17712 , -7.27712 , -7.37712 , -7.47712 , -7.57712 , -7.67712 , -7.77712 , -7.87712 , -7.97712 , -8.07712 , -8.17712 , -8.27712 , -8.37712 , -8.47712 , -8.57712 , -8.67712 , -8.77712 , -8.87712 , -8.97712 , -9.07712 , -9.17712 , -9.27712 , -9.37712 , -9.47712 , -9.57712 , -9.67712 , -9.77712 , -9.87712 , -9.97712 , -10.0771 , -10.1771 , -10.2771 , -10.3771 , -10.4771 , -10.5771 , -10.6771 , -10.7771 , -10.8771 , -10.9771 , -11.0771 , -11.1771 , -11.2771 , -11.3771 , -11.4771 , -11.5771 , -11.6771 , -11.7771 , -11.8771 , -11.9771 , -12.0771 , -12.1771 , -12.2771 , -12.3771 , -12.4771 , -12.5771 , -12.6771 , -12.7771 , -12.8771 , -12.9771 , -13.0771 , -13.1771 , -13.2771 , -13.3771 , -13.4771 , -13.5771 , -13.6771 , -13.7771 , -13.8771 , -13.9771 , -14.0771 , -14.1771 , -14.2771 , -14.3771 , -14.4771 , -14.5771 , -14.6771 , -14.7771 , -14.8771 , -14.9771 , -15.0771 , -15.1771 , -15.2771 , -15.3771 , -15.4771 , -15.5771 , -15.6771 , -15.7771 , -15.8771 , -15.9771 , -16.0771 , -16.1771 , -16.2771 , -16.3771 , -16.4771 , -16.5771 , -16.6771 , -16.7771 , -16.8771 , -16.9771 , -17.0771 , -17.1771 , -17.2771 , -17.3771 , -17.4771 , -17.5771 , -17.6771 , -17.7771 , -17.8771 , -17.9771 , -18.0771 , -18.1771 , -18.2771 , -18.3771 , -18.4771 , -18.5771 , -18.6771 , -18.7771 , -18.8771 , -18.9771 , -19.0771 , -19.1771 , -19.2771 , -19.3771 , -19.4771 , -19.5771 , -19.6771 , -19.7771 , -19.8771 , -19.9771 , -20.0771 , -20.1771 , -20.2771 , -20.3771 , -20.4771 , -20.5771 , -20.6771 , -20.7771 , -20.8771 , -20.9771 , -21.0771 , -21.1771 , -21.2771 , -21.3771 , -21.4771 , -21.5771 , -21.6771 , -21.7771 , -21.8771 , -21.9771 , -22.0771 , -22.1771 , -22.2771 , -22.3771 , -22.4771 , -22.5771 , -22.6771 , -22.7771 , -22.8771 , -22.9771 , -23.0771 , -23.1771 , -23.2771 , -23.3771 , -23.4771 , -23.5771 , -23.6771 , -23.7771 , -23.8771 , -23.9771 , -24.0771 , -24.1771 , -24.2771 , -24.3771 , -24.4771 , -24.5771 , -24.6771 , -24.7771 , -24.8771 , -24.9771 , -25.0771 , -25.1771 , -25.2771 , -25.3771 , -25.4771 , -25.5771 , -25.6771 , -25.7771 , -25.8771 , -25.9771};
	const static float log10_one_div_two_minus_err_div_three[256] = {-0.778151 , -0.628519 , -0.53808 , -0.477637 , -0.434982 , -0.403853 , -0.380624 , -0.36302 , -0.349527 , -0.339101 , -0.330993 , -0.324659 , -0.319693 , -0.315789 , -0.312712 , -0.310284 , -0.308364 , -0.306846 , -0.305643 , -0.30469 , -0.303935 , -0.303336 , -0.302861 , -0.302483 , -0.302184 , -0.301947 , -0.301758 , -0.301608 , -0.301489 , -0.301395 , -0.30132 , -0.30126 , -0.301213 , -0.301175 , -0.301145 , -0.301122 , -0.301103 , -0.301088 , -0.301076 , -0.301066 , -0.301059 , -0.301053 , -0.301048 , -0.301044 , -0.301042 , -0.301039 , -0.301037 , -0.301036 , -0.301035 , -0.301034 , -0.301033 , -0.301032 , -0.301032 , -0.301031 , -0.301031 , -0.301031 , -0.301031 , -0.301031 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103 , -0.30103};
	const static float log10_err[256]={0, -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1, -1.1, -1.2, -1.3, -1.4, -1.5, -1.6, -1.7, -1.8, -1.9, -2, -2.1, -2.2, -2.3, -2.4, -2.5, -2.6, -2.7, -2.8, -2.9, -3, -3.1, -3.2, -3.3, -3.4, -3.5, -3.6, -3.7, -3.8, -3.9, -4, -4.1, -4.2, -4.3, -4.4, -4.5, -4.6, -4.7, -4.8, -4.9, -5, -5.1, -5.2, -5.3, -5.4, -5.5, -5.6, -5.7, -5.8, -5.9, -6, -6.1, -6.2, -6.3, -6.4, -6.5, -6.6, -6.7, -6.8, -6.9, -7, -7.1, -7.2, -7.3, -7.4, -7.5, -7.6, -7.7, -7.8, -7.9, -8, -8.1, -8.2, -8.3, -8.4, -8.5, -8.6, -8.7, -8.8, -8.9, -9, -9.1, -9.2, -9.3, -9.4, -9.5, -9.6, -9.7, -9.8, -9.9, -10, -10.1, -10.2, -10.3, -10.4, -10.5, -10.6, -10.7, -10.8, -10.9, -11, -11.1, -11.2, -11.3, -11.4, -11.5, -11.6, -11.7, -11.8, -11.9, -12, -12.1, -12.2, -12.3, -12.4, -12.5, -12.6, -12.7, -12.8, -12.9, -13, -13.1, -13.2, -13.3, -13.4, -13.5, -13.6, -13.7, -13.8, -13.9, -14, -14.1, -14.2, -14.3, -14.4, -14.5, -14.6, -14.7, -14.8, -14.9, -15, -15.1, -15.2, -15.3, -15.4, -15.5, -15.6, -15.7, -15.8, -15.9, -16, -16.1, -16.2, -16.3, -16.4, -16.5, -16.6, -16.7, -16.8, -16.9, -17, -17.1, -17.2, -17.3, -17.4, -17.5, -17.6, -17.7, -17.8, -17.9, -18, -18.1, -18.2, -18.3, -18.4, -18.5, -18.6, -18.7, -18.8, -18.9, -19, -19.1, -19.2, -19.3, -19.4, -19.5, -19.6, -19.7, -19.8, -19.9, -20, -20.1, -20.2, -20.3, -20.4, -20.5, -20.6, -20.7, -20.8, -20.9, -21, -21.1, -21.2, -21.3, -21.4, -21.5, -21.6, -21.7, -21.8, -21.9, -22, -22.1, -22.2, -22.3, -22.4, -22.5, -22.6, -22.7, -22.8, -22.9, -23, -23.1, -23.2, -23.3, -23.4, -23.5, -23.6, -23.7, -23.8, -23.9, -24, -24.1, -24.2, -24.3, -24.4, -24.5, -24.6, -24.7, -24.8, -24.9, -25, -25.1, -25.2, -25.3, -25.4, -25.5};
	const static float log10_one_div_two_minus_err[256] = {-1, -1, -1, -1, -0.991857, -0.73572, -0.60413, -0.522193, -0.466596, -0.427004, -0.39794, -0.376165, -0.359614, -0.346902, -0.337063, -0.329404, -0.323415, -0.318716, -0.315019, -0.312105, -0.309804, -0.307985, -0.306545, -0.305405, -0.304502, -0.303785, -0.303217, -0.302767, -0.302409, -0.302125, -0.301899, -0.30172, -0.301578, -0.301466, -0.301376, -0.301305, -0.301248, -0.301203, -0.301168, -0.301139, -0.301117, -0.301099, -0.301085, -0.301074, -0.301065, -0.301057, -0.301052, -0.301047, -0.301044, -0.301041, -0.301039, -0.301037, -0.301035, -0.301034, -0.301033, -0.301033, -0.301032, -0.301032, -0.301031, -0.301031, -0.301031, -0.301031, -0.301031, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103, -0.30103};
}

/**
 * @brief Minimal GLF (Genotype Likelihood Format) generator function for genotype calling.
 * 
 * This function processes the pileup reads at a variant position and extracts
 * high-quality base calls to prepare data for genotype likelihood computations.
 * It filters out low-quality bases, reference skips, and deletions.
 * The base information is encoded into an integer containing base quality, strand,
 * and base identity.
 * 
 * @param[in,out] caller A reference to the call_t object containing pileup data and storage
 *                       for base calls and genotype likelihood statistics.
 *                       - Uses caller.snp_called to skip processing if the variant is already called.
 *                       - Uses caller.n_plp and caller.v_plp for the pileup reads.
 *                       - Clears and fills caller.bca.bases with encoded bases.
 *                       - Updates caller.gls.read_depth to the number of high-quality bases.
 * @param[in] variant Pointer to the variant object being processed (unused in this minimal implementation).
 * 
 * @note The encoded integer stored in caller.bca.bases encodes information as follows:
 *       - Bits 0-3: Base identity as integer (0=A,1=C,2=G,3=T)
 *       - Bit 4: Reverse strand flag (1 if read aligned to reverse strand, 0 otherwise)
 *       - Bits 5 and above: Base quality score (phred-scaled)
 * 
 * @details
 * The function loops over all reads in the pileup:
 * - Skips reads that are reference skips or deletions.
 * - Extracts the base and base quality at the variant position.
 * - Filters out bases not representing A, C, G, or T.
 * - Filters out bases below the minimum base quality threshold.
 * - Encodes and stores the base information in caller.bca.bases.
 * - Updates the read depth metric after processing.
 */
void bcf_call_glfgen_min(call_t& caller, const variant* variant);

/**
 * @brief Minimal GLF generator function specialized for insertion variants.
 * 
 * This function processes pileup reads at an insertion variant site, filtering reads
 * to find those supporting the specific insertion allele. It distinguishes reads
 * matching the reference base from those supporting the insertion event, and encodes
 * base qualities and strand information accordingly.
 * 
 * @param[in,out] caller Reference to the call_t object containing pileup data and storage
 *                       for base calls and genotype likelihood statistics.
 *                       - Uses caller.n_plp and caller.v_plp to access pileup reads.
 *                       - Clears and fills caller.bca.bases with encoded base/indel info.
 *                       - Updates caller.gls.read_depth with the count of informative bases.
 * @param[in] variant Pointer to the variant object representing the insertion to be called.
 *                    The variant's reference and alternate allele sequences are used
 *                    to identify and validate supporting reads.
 * 
 * @details
 * The function operates as follows:
 * - If there are no pileup reads, returns immediately.
 * - Only processes variants where the reference allele is a single base and
 *   the alternate allele length is greater than 1 (i.e., an insertion).
 * - Iterates over all pileup reads:
 *   - Skips reads flagged as reference skips or deletions.
 *   - Checks if the base in the read matches the reference base at the variant position.
 *   - If the read has an indel (insertion/deletion) at this position:
 *     - Validates that the indel length matches the expected insertion length.
 *     - Verifies that the inserted sequence matches the alternate allele bases.
 *     - Aggregates base quality scores across the inserted bases and averages them.
 *   - Filters out reads with average base quality below the minimum threshold.
 *   - Encodes each base/indel call into an integer with the following bits:
 *     - Bits 0-3: Indicator flag (0 for reference, 1 for indel present)
 *     - Bit 4: Strand information (1 if reverse strand, 0 otherwise)
 *     - Bits 5 and above: Base quality score (phred-scaled)
 * - Stores encoded values in caller.bca.bases.
 * - Updates caller.gls.read_depth to the number of processed reads.
 * 
 * @note Encoding differs from SNP calls: the lowest bit indicates indel presence (1)
 *       or absence (0), instead of base identity.
 */
void bcf_call_glfgen_min_ins(call_t& caller, const variant* variant);

/**
 * @brief Minimal GLF generator function specialized for deletion variants.
 * 
 * This function processes pileup reads at a deletion variant site, filtering reads
 * that support either the reference or deletion alleles. It evaluates base matches,
 * indel sizes, and base qualities to encode reads contributing evidence for genotype calling.
 * 
 * @param[in,out] caller Reference to the call_t object containing pileup data and storage
 *                       for base calls and genotype likelihood statistics.
 *                       - Uses caller.n_plp and caller.v_plp to access pileup reads.
 *                       - Clears and fills caller.bca.bases with encoded base/indel info.
 *                       - Updates caller.gls.read_depth with the count of informative bases.
 * @param[in] variant Pointer to the variant object representing the deletion to be called.
 *                    The variant's reference and alternate allele sequences are used
 *                    to identify and validate supporting reads.
 * 
 * @details
 * The function operates as follows:
 * - Returns immediately if no pileup reads.
 * - Only processes variants where the reference allele length is greater than 1
 *   (indicating a deletion) and the alternate allele is a single base.
 * - Iterates over all pileup reads:
 *   - Skips reads flagged as reference skips or deletions.
 *   - Checks if the base in the read matches the reference base at the variant position.
 *   - If the read has no indel:
 *     - Verifies that the bases following the variant position in the read match the
 *       full reference allele sequence.
 *     - Aggregates and averages base quality scores across these bases.
 *   - If the read has an indel:
 *     - Accepts reads with an indel length consistent with the expected deletion length.
 *     - Skips reads with indels inconsistent with the deletion.
 *   - Filters out reads with base quality below the minimum threshold.
 *   - Encodes each base/indel call into an integer with the following bits:
 *     - Bits 0-3: Indicator flag (0 for reference, 1 for indel presence)
 *     - Bit 4: Strand information (1 if reverse strand, 0 otherwise)
 *     - Bits 5 and above: Base quality score (phred-scaled)
 * - Stores encoded values in caller.bca.bases.
 * - Updates caller.gls.read_depth to the number of processed reads.
 * 
 * @note Encoding uses bit 0 as an indel presence flag (1 = indel, 0 = reference).
 */
void bcf_call_glfgen_min_del(call_t& caller, const variant* variant);

/**
 * @brief Clears the current base calls and resets read depth.
 * 
 * This function effectively empties any stored base information in the caller's
 * base call array and resets the genotype likelihood read depth to zero.
 * It can be used as a placeholder or to explicitly indicate that no
 * bases are called for the given variant.
 * 
 * @param[in,out] caller Reference to the call_t object whose base calls and
 *                      read depth will be cleared.
 * @param[in] variant Pointer to the variant object (unused in this function,
 *                    included for interface consistency).
 */
void bcf_call_void(call_t& caller, const variant* variant);

/**
 * @brief Compute genotype likelihoods for a variant using a standard sequencing error model.
 * 
 * This function calculates the log-likelihoods (GLs) of possible genotypes for a specific variant
 * based on the aligned sequencing reads at that variant site. The likelihoods incorporate
 * sequencing base qualities and model errors probabilistically to distinguish true variants
 * from sequencing noise.
 * 
 * The method:
 * 1. Extracts the encoded reference and alternate alleles from the variant.
 * 2. Checks validity of these alleles.
 * 3. Limits the number of reads considered for likelihood calculation if exceeding a threshold,
 *    to avoid excessive computation.
 * 4. Sorts the base calls from aligned reads and groups identical bases together.
 * 5. For each group of identical bases, calculates the contribution to the likelihoods based
 *    on base quality scores and an error model, differentiating between homozygous and 
 *    heterozygous genotype states.
 * 6. Updates the caller's genotype likelihoods for all genotype possibilities.
 * 
 * The likelihoods are represented in log10 scale, where higher values indicate greater support
 * for a genotype given the observed reads and qualities.
 * 
 * @param[in,out] caller A data structure representing the current variant calling state,
 *                      containing base calls, genotype likelihoods, read depths, and 
 *                      intermediate computation buffers. This function updates the genotype 
 *                      likelihoods within this object.
 * @param[in] variant A pointer to the variant being evaluated, which contains the reference
 *                    and alternate allele sequences.
 * 
 * @details
 * - The function uses precomputed lookup tables (`calling_utils::log10_one_minus_err`, 
 *   `calling_utils::log10_err_div_three`, etc.) to translate base quality scores into error 
 *   probabilities in log10 space.
 * - For diploid samples (ploidy > 1), it models all possible genotype combinations between 
 *   the reference and alternate alleles, including homozygous reference, heterozygous, and 
 *   homozygous alternate.
 * - For haploid samples (ploidy = 1), it only considers homozygous genotypes.
 * - The function also keeps track of the number of reads supporting the reference allele 
 *   (`ref_read`) and the total depth (`dp_ind`) after filtering.
 * - If the number of reads exceeds the maximum allowed (`caller.bca.max_dp`), a random 
 *   subset of reads is used by shuffling to avoid bias.
 * 
 * @throws runtime_error If the reference or alternate alleles contain invalid bases.
 * 
 * @note This function assumes base calls and qualities have already been collected and filtered 
 *       into `caller.bca.bases` before being called.
 */
void standard_errmod(call_t& caller, const variant* variant);

/**
 * @brief Compute genotype likelihoods for small indel variants using a simplified error model.
 *
 * This function calculates log-likelihoods for possible genotypes at an indel variant site 
 * based on observed sequencing reads, considering base quality scores and read counts.
 * It supports both haploid and diploid samples.
 *
 * The model:
 * - Treats bases as either reference (0) or alternate (1).
 * - Groups identical bases together to efficiently compute likelihood contributions.
 * - Uses precomputed log10 error probabilities based on base quality scores.
 * - Applies a probabilistic model to compute likelihoods for homozygous reference, 
 *   heterozygous, and homozygous alternate genotypes (diploid), or homozygous genotypes (haploid).
 * - Caps the maximum number of reads considered to avoid excessive computation by randomly 
 *   sampling if necessary.
 *
 * @param[in,out] caller The variant calling context containing read bases, quality info, 
 *                       genotype likelihood buffers, and ploidy information.
 * @param[in] variant Pointer to the variant under consideration (not directly used here 
 *                    but kept for API consistency).
 *
 * @details
 * - Reads are encoded in `caller.bca.bases` with qualities and base calls packed into a uint16_t.
 * - Reads are sorted and processed in groups of identical encoded bases.
 * - For diploid, likelihoods for genotypes are updated as:
 *      - homozygous genotype likelihoods are increased proportionally to error-free base quality.
 *      - heterozygous likelihood incorporates half the error-adjusted probability.
 * - For haploid, only two genotypes are modeled.
 * - The depth of informative reads and number of reference-supporting reads are tracked.
 * - If ploidy is haploid, the likelihood for a third genotype is set to a large negative value 
 *   (effectively ignored).
 *
 * @note This function assumes the bases and qualities have been pre-collected and filtered.
 */
void standard_errmod_indels_v(call_t& caller, const variant* variant);

/**
 * @brief Reset genotype likelihoods and read counts for indel variant calling.
 *
 * This function clears all genotype likelihood values and resets read depth and reference 
 * allele read counts to zero. It effectively initializes the calling context to a neutral 
 * state without incorporating any observed data.
 *
 * This can be used as a placeholder or fallback error model where no data is available or 
 * when an indel site should be treated as missing data.
 *
 * @param[in,out] caller The variant calling context whose genotype likelihoods and read counts 
 *                       will be reset.
 * @param[in] variant Pointer to the variant being evaluated (unused in this function).
 */
void flat_errmod_indels_v(call_t& caller, const variant* varian);


/**
 * @class genotype_bam_caller
 * @brief Class to perform genotype calling from BAM files using haplotype and variant information.
 * 
 * This class encapsulates all necessary data structures and methods
 * to call genotypes in a specified genomic region based on sequencing
 * reads processed by mpileup and haplotype sets.
 */
class genotype_bam_caller
{
public:
    // References to input data
    haplotype_set & H;                 ///< Reference to a set of haplotypes used for calling
    genotype_set & G;                  ///< Reference to the output genotype set to be filled
    const variant_map& V;              ///< Reference to the map of known variants
    const glimpse_mpileup& mpileup_data; ///< Reference to mpileup data for the target region
    const std::string region;          ///< Genomic region string (e.g., "chr1:1000-2000")
    const uint32_t beg;                ///< Begin coordinate of the calling region (1-based)
    const uint32_t end;                ///< End coordinate of the calling region (1-based, inclusive)

    aux_t aux_data;                   ///< Auxiliary data container for intermediate calculations

    // Caller object to hold state of the current genotype call
    call_t caller;                   ///< The call state that tracks the genotype calling progress

    // Model for genotype likelihood calculations
    call_model model;                ///< Model used for calling, contains parameters and methods for likelihoods

    // Vector of functions to group reads according to some criteria before calling
    std::vector<std::function<void(call_t& caller, const variant* variant)>> group_reads;

    // Vector of functions to compute log-likelihoods (LLK) for genotype calls
    std::vector<std::function<void(call_t& caller, const variant* variant)>> compute_llk;

	/**
	 * @brief Constructs a genotype_bam_caller object for genotype calling from BAM mpileup data.
	 * 
	 * Initializes internal state and sets up function pointers for different
	 * read grouping and likelihood computation models based on parameters.
	 * 
	 * @param _H Reference to a haplotype_set object containing haplotype data.
	 * @param _G Reference to a genotype_set object where called genotypes will be stored.
	 * @param _V Reference to a variant_map object describing variant positions and info.
	 * @param _mpileup_data Reference to glimpse_mpileup object containing mpileup read data and parameters.
	 * @param _region Genomic region string (e.g., "chr1:100000-200000") to call variants within.
	 * @param _beg Start coordinate of the region (0-based).
	 * @param _end End coordinate of the region (0-based, exclusive).
	 * @param _callmodel String specifying which error model to use for genotype likelihood computation.
	 *                   Currently only supports "standard".
	 * @param call_indels Boolean flag indicating whether indels (insertions/deletions) should be called.
	 * 
	 * @throws std::runtime_error if the specified call model is not supported.
	 * 
	 * @details
	 * Sets up internal buffers and parameters from mpileup data, reserves space for base qualities,
	 * and assigns function pointers for handling SNPs and indels. Error models for likelihood
	 * computation are set according to the chosen model and whether indels are called.
	 */
	genotype_bam_caller(haplotype_set & H, genotype_set & G, const variant_map& V, const glimpse_mpileup& _mpileup_data, const std::string region, const uint32_t beg, const uint32_t end, std::string callmodel, bool call_indels);
	
	/**
	 * @brief Destructor for genotype_bam_caller class.
	 * 
	 * Cleans up and releases resources associated with BAM file reading:
	 * closes file pointers, destroys BAM header, iterator, and index objects
	 * if they were initialized during the caller's lifetime.
	 */
	virtual ~genotype_bam_caller();

	/**
	 * @brief Perform genotype calling using mpileup for the i-th BAM/CRAM file.
	 * 
	 * This function attempts to open and process the BAM/CRAM file for the specified sample index,
	 * retrying up to 3 times with exponential backoff in case of failures (e.g., streaming issues).
	 * It opens the file, sets options, reads headers, loads indexes, and queries the specified region.
	 * Then, it calls glimpse_mpileup_reg to perform the pileup-based genotype calling.
	 * 
	 * If any step fails, it retries after cleaning resources and waits increasingly longer.
	 * After max retries, the function reports a fatal error.
	 * 
	 * @param i Index of the sample/BAM file to process.
	 */
	void call_mpileup(int i);

	/**
	 * @brief Clean up and release all resources used during BAM/CRAM processing.
	 *
	 * This method releases any open iterators, indexes, headers, and file pointers
	 * associated with the BAM/CRAM file being processed. It also resets internal
	 * caller statistics such as read depth, per-individual depth, likelihoods, and ploidy.
	 * This function is typically called after each BAM file processing attempt or
	 * upon failure to ensure a clean state for retries or termination.
	 */
	void clean();

	/**
	 * @brief Reset genotype calling results for the specified sample index.
	 *
	 * Clears the coverage data and resets the depth counts to zero for sample `i`.
	 * This is typically used to discard partial or failed calling results before retrying
	 * or aborting the processing for a particular sample.
	 *
	 * @param i Index of the sample whose results are to be reset.
	 */
	void reset_results(int i);

	/**
	 * @brief Perform pileup-based genotype likelihood calculation for a single sample over a genomic region.
	 * 
	 * This function iterates over sequencing read pileups in the BAM/CRAM file for sample index `i`, 
	 * matching the pileup positions to known variant sites in the region. For each variant site, it
	 * computes genotype likelihoods (GLs) using a specified model and updates genotype and coverage statistics.
	 * 
	 * The method handles both SNPs and indels and supports diploid or haploid samples depending on `mpileup_data.tar_ploidy[i]`.
	 * 
	 * @param i Index of the sample to process.
	 * 
	 * @return int Returns 0 on success, or -1 on failure (e.g., pileup stream error).
	 * 
	 * @details
	 * - Initializes a BAM pileup iterator over the specified region.
	 * - For each pileup position, it synchronizes with the next variant site from `V.vec_pos`.
	 * - If the pileup position matches a variant site, genotype likelihoods are computed using
	 *   the configured likelihood model for SNPs or indels.
	 * - Updates `G.vecG[i]->GL` with scaled likelihoods.
	 * - Updates coverage statistics in `G.stats`.
	 * - Uses likelihood scaling to store phred-scaled genotype likelihoods as unsigned chars.
	 * - Supports multiallelic genotypes based on sample ploidy.
	 * 
	 * @note
	 * - The function uses htslib pileup APIs for read iteration.
	 * - Sites outside the region or without reads are filled with zero coverage and likelihoods.
	 */
	int glimpse_mpileup_reg(int i);
};

#endif

