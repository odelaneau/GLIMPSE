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

#ifndef _OLIVIER_TOOLS_H
#define _OLIVIER_TOOLS_H

//INCLUDE STANDARD TEMPLATE LIBRARY USEFULL STUFFS (STL)
#include <vector>
#include <list>
#include <queue>
#include <stack>
#include <bitset>
#include <set>
#include <tuple>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <string>
#include <exception>
#include <cassert>
#include <limits>
#include <cstdint>

//INCLUDE BOOST USEFULL STUFFS (BOOST)
#include <boost/program_options.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/align/aligned_allocator.hpp>

//INCLUDE HTS LIBRARY
#include <htslib/hts.h>
#include <htslib/kseq.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/regidx.h>
extern "C" {
	#include <htslib/vcf_sweep.h>
	#include <htslib/synced_bcf_reader.h>
	#include <htslib/vcf.h>
	#include <htslib/vcfutils.h>
}

//INCLUDES BASE STUFFS
#include <utils/compressed_io.h>
#include <utils/random_number.h>
#include <utils/basic_stats.h>
#include <utils/basic_algos.h>
#include <utils/string_utils.h>
#include <utils/timer.h>
#include <utils/verbose.h>

//TYPEDEFS
template <typename T>
using aligned_vector32 = std::vector<T, boost::alignment::aligned_allocator < T, 32 > >;

//CONSTANTS
#define RARE_VARIANT_FREQ	0.001f
#define HAP_NUMBER			8
#define MAX_AMB				22
#define STAGE_INIT	0
#define STAGE_BURN	1
#define STAGE_MAIN	2
#define STAGE_RESTRICT	3

//MACROS
#define DIV2(v)	((v)>>1)
#define MOD2(v)	((v)&1)

//NAMESPACE
namespace bio = boost::iostreams;
namespace bpo = boost::program_options;
namespace bid = boost::uuids;

//equivalent to: for (int i=0; i<256; i++) unphred[i] = pow(10., -i/10.);
const static float unphred[256] = { 1.000000e+00, 7.943282e-01, 6.309573e-01, 5.011872e-01, 3.981072e-01, 3.162278e-01, 2.511886e-01, 1.995262e-01, 1.584893e-01, 1.258925e-01, 1.000000e-01, 7.943282e-02, 6.309573e-02, 5.011872e-02, 3.981072e-02, 3.162278e-02, 2.511886e-02, 1.995262e-02, 1.584893e-02, 1.258925e-02, 1.000000e-02, 7.943282e-03, 6.309573e-03, 5.011872e-03, 3.981072e-03, 3.162278e-03, 2.511886e-03, 1.995262e-03, 1.584893e-03, 1.258925e-03, 1.000000e-03, 7.943282e-04, 6.309573e-04, 5.011872e-04, 3.981072e-04, 3.162278e-04, 2.511886e-04, 1.995262e-04, 1.584893e-04, 1.258925e-04, 1.000000e-04, 7.943282e-05, 6.309573e-05, 5.011872e-05, 3.981072e-05, 3.162278e-05, 2.511886e-05, 1.995262e-05, 1.584893e-05, 1.258925e-05, 1.000000e-05, 7.943282e-06, 6.309573e-06, 5.011872e-06, 3.981072e-06, 3.162278e-06, 2.511886e-06, 1.995262e-06, 1.584893e-06, 1.258925e-06, 1.000000e-06, 7.943282e-07, 6.309573e-07, 5.011872e-07, 3.981072e-07, 3.162278e-07, 2.511886e-07, 1.995262e-07, 1.584893e-07, 1.258925e-07, 1.000000e-07, 7.943282e-08, 6.309573e-08, 5.011872e-08, 3.981072e-08, 3.162278e-08, 2.511886e-08, 1.995262e-08, 1.584893e-08, 1.258925e-08, 1.000000e-08, 7.943282e-09, 6.309573e-09, 5.011872e-09, 3.981072e-09, 3.162278e-09, 2.511886e-09, 1.995262e-09, 1.584893e-09, 1.258925e-09, 1.000000e-09, 7.943282e-10, 6.309573e-10, 5.011872e-10, 3.981072e-10, 3.162278e-10, 2.511886e-10, 1.995262e-10, 1.584893e-10, 1.258925e-10, 1.000000e-10, 7.943282e-11, 6.309573e-11, 5.011872e-11, 3.981072e-11, 3.162278e-11, 2.511886e-11, 1.995262e-11, 1.584893e-11, 1.258925e-11, 1.000000e-11, 7.943282e-12, 6.309573e-12, 5.011872e-12, 3.981072e-12, 3.162278e-12, 2.511886e-12, 1.995262e-12, 1.584893e-12, 1.258925e-12, 1.000000e-12, 7.943282e-13, 6.309573e-13, 5.011872e-13, 3.981072e-13, 3.162278e-13, 2.511886e-13, 1.995262e-13, 1.584893e-13, 1.258925e-13, 1.000000e-13, 7.943282e-14, 6.309573e-14, 5.011872e-14, 3.981072e-14, 3.162278e-14, 2.511886e-14, 1.995262e-14, 1.584893e-14, 1.258925e-14, 1.000000e-14, 7.943282e-15, 6.309573e-15, 5.011872e-15, 3.981072e-15, 3.162278e-15, 2.511886e-15, 1.995262e-15, 1.584893e-15, 1.258925e-15, 1.000000e-15, 7.943282e-16, 6.309573e-16, 5.011872e-16, 3.981072e-16, 3.162278e-16, 2.511886e-16, 1.995262e-16, 1.584893e-16, 1.258925e-16, 1.000000e-16, 7.943282e-17, 6.309573e-17, 5.011872e-17, 3.981072e-17, 3.162278e-17, 2.511886e-17, 1.995262e-17, 1.584893e-17, 1.258925e-17, 1.000000e-17, 7.943282e-18, 6.309573e-18, 5.011872e-18, 3.981072e-18, 3.162278e-18, 2.511886e-18, 1.995262e-18, 1.584893e-18, 1.258925e-18, 1.000000e-18, 7.943282e-19, 6.309573e-19, 5.011872e-19, 3.981072e-19, 3.162278e-19, 2.511886e-19, 1.995262e-19, 1.584893e-19, 1.258925e-19, 1.000000e-19, 7.943282e-20, 6.309573e-20, 5.011872e-20, 3.981072e-20, 3.162278e-20, 2.511886e-20, 1.995262e-20, 1.584893e-20, 1.258925e-20, 1.000000e-20, 7.943282e-21, 6.309573e-21, 5.011872e-21, 3.981072e-21, 3.162278e-21, 2.511886e-21, 1.995262e-21, 1.584893e-21, 1.258925e-21, 1.000000e-21, 7.943282e-22, 6.309573e-22, 5.011872e-22, 3.981072e-22, 3.162278e-22, 2.511886e-22, 1.995262e-22, 1.584893e-22, 1.258925e-22, 1.000000e-22, 7.943282e-23, 6.309573e-23, 5.011872e-23, 3.981072e-23, 3.162278e-23, 2.511886e-23, 1.995262e-23, 1.584893e-23, 1.258925e-23, 1.000000e-23, 7.943282e-24, 6.309573e-24, 5.011872e-24, 3.981072e-24, 3.162278e-24, 2.511886e-24, 1.995262e-24, 1.584893e-24, 1.258925e-24, 1.000000e-24, 7.943282e-25, 6.309573e-25, 5.011872e-25, 3.981072e-25, 3.162278e-25, 2.511886e-25, 1.995262e-25, 1.584893e-25, 1.258925e-25, 1.000000e-25, 7.943282e-26, 6.309573e-26, 5.011872e-26, 3.981072e-26, 3.162278e-26};

enum class InputFormat
{
	BCF = 0,
	GLIMPSE = 1
};

enum class OutputFormat
{
	VCF = 0,
	BCF = 1,
	BGEN = 2
};

enum class OutputCompression
{
	NONE = 0,
	ZLIB = 1,
	ZSTD = 2
};

//MAKE SOME TOOLS FULLY ACCESSIBLE THROUGHOUT THE SOFTWARE
#ifdef _DECLARE_TOOLBOX_HERE
	random_number_generator rng;	//Random number generator
	string_utils stb;				//String manipulation
	basic_algos alg;				//Basic algorithms
	verbose vrb;					//Verbose
	timer tac;						//Timer
#else
	extern random_number_generator rng;
	extern string_utils stb;
	extern basic_algos alg;
	extern verbose vrb;
	extern timer tac;
#endif

#endif
