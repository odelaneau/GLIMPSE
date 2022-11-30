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

	unsigned int rel_time_micro() {
		return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - prev_timing_clock).count();
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

	std::string abs_time_display()
	{
		std::chrono::seconds input_seconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_timing_clock);
	    using namespace std::chrono;
	    typedef duration<int, std::ratio<86400>> days;
	    auto d = duration_cast<days>(input_seconds);
	    input_seconds -= d;
	    auto h = duration_cast<hours>(input_seconds);
	    input_seconds -= h;
	    auto m = duration_cast<minutes>(input_seconds);
	    input_seconds -= m;
	    auto s = duration_cast<seconds>(input_seconds);

	    auto dc = d.count();
	    auto hc = h.count();
	    auto mc = m.count();
	    auto sc = s.count();

	    std::stringstream ss;
	    ss.fill('0');
	    if (dc) {
	        ss << d.count() << "d";
	    }
	    if (dc || hc) {
	        if (dc) { ss << std::setw(2); } //pad if second set of numbers
	        ss << h.count() << " hour";
	        if (h.count() > 1) ss << "s";
	        ss << " ";
	    }
	    if (dc || hc || mc) {
	        if (dc || hc) { ss << std::setw(2); }
	        ss << m.count() << " minute";
	        if (m.count() > 1) ss << "s";
	        ss << " ";
	    }
	    if (dc || hc || mc || sc) {
	        if (dc || hc || mc) { ss << std::setw(2); }
	        ss << s.count() << " second";
	        if (s.count() > 1) ss << "s";
	    }

	    return ss.str();
	}
};

#endif
