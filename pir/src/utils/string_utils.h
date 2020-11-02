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

#ifndef _STRING_UTILS_H
#define _STRING_UTILS_H

#include <sstream>
#include <iomanip>
#include <string>
#include <vector>

//using namespace std;

class string_utils {
public:
	string_utils () {};
	~string_utils () {};

	int split(const std::string & str, std::vector < std::string > & tokens, std::string sep = " 	", unsigned int n_max_tokens = 1000000) {
		tokens.clear();
		if (str == ""){
			tokens.push_back("");
			return tokens.size();
		}
		std::string::size_type p_last = str.find_first_not_of(sep, 0);
		std::string::size_type p_curr = str.find_first_of(sep, p_last);
		while ((std::string::npos != p_curr || std::string::npos != p_last) && tokens.size() < n_max_tokens) {
			tokens.push_back(str.substr(p_last, p_curr - p_last));
			p_last = str.find_first_not_of(sep, p_curr);
			p_curr = str.find_first_of(sep, p_last);
		}
		if (tokens.back()[tokens.back().size()-1] == '\r') tokens.back() = tokens.back().substr(0, tokens.back().size()-1);
		return tokens.size();
	}

	bool numeric(std::string & str) {
		double n;
		std::istringstream in(str);
		if (!(in >> n)) return false;
		return true;
    }

	template < class T >
	std::string str(T n, int prec = -1) {
		std::ostringstream ss( std::stringstream::out );
		if (prec >= 0) { ss << setiosflags( std::ios::fixed ); ss.precision(prec); }
		ss << n;
		return ss.str();
	}

	template < class T >
	std::string str(std::vector < T > & v, int prec = -1) {
		std::ostringstream ss( std::stringstream::out );
		if (prec >= 0) { ss << setiosflags( std::ios::fixed ); ss.precision(prec); }
		for (int e = 0 ; e < v.size() ; e ++) ss << (e>0?" ":"") << v[e] ;
		return ss.str();
	}
};

#endif
