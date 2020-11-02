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

#ifndef _TIMER_H
#define _TIMER_H

#include <chrono>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <string>

class timer {
protected:
	std::chrono::time_point<std::chrono::high_resolution_clock> start_timing_clock, prev_timing_clock;

public:
	timer () {
		start_timing_clock = std::chrono::high_resolution_clock::now();
	}

	~timer() {
	}

	void clock() {
		prev_timing_clock = std::chrono::high_resolution_clock::now();
	}

	unsigned int rel_time() {
		return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - prev_timing_clock).count();
	}

	unsigned int abs_time() {
		return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_timing_clock).count();
	}

	std::string date() {
		auto now = std::chrono::system_clock::now();
		auto in_time_t = std::chrono::system_clock::to_time_t(now);
		std::stringstream ss;
	    ss << std::put_time(std::localtime(&in_time_t), "%d/%m/%Y - %X");
	    return ss.str();
	}
};

#endif
