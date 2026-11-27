/*******************************************************************************
 * Copyright (C) 2026 Simone Rubinacci
 * Copyright (C) 2026 Olivier Delaneau
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

#ifndef _LIBRARY_VERSIONS_H
#define _LIBRARY_VERSIONS_H

#include <string>

#include <htslib/hts.h>
#include <boost/version.hpp>

/**
 * Versions of the libraries that matter when diagnosing problems, for the log banner,
 * e.g. "HTSlib 1.25 / Boost 1.92.0". HTSlib is the version loaded at run time. Boost is the
 * version the tool was built with, which decides the archive version of the binary reference
 * panels (.bin): panels written with Boost 1.84 or newer cannot be read by Boost 1.83 or older.
 */
inline std::string library_versions() {
	return "HTSlib " + std::string(hts_version()) +
		" / Boost " + std::to_string(BOOST_VERSION / 100000) + "." + std::to_string(BOOST_VERSION / 100 % 1000) + "." + std::to_string(BOOST_VERSION % 100);
}

#endif
