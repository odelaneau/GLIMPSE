/*******************************************************************************
 * Copyright (C) 2018 Olivier Delaneau, University of Lausanne
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
/*Copyright (C) 2015 Olivier Delaneau, University of Lausanne, Halit Ongen, Emmanouil T. Dermitzakis
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

#ifndef _VERBOSE_H
#define _VERBOSE_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

class verbose {
protected:
	std::ofstream log;
	bool verbose_on_screen;
	bool verbose_on_log;
	int prev_percent;

public:
	verbose() {
		verbose_on_screen = true;
		verbose_on_log = false;
		prev_percent = -1;
	}

	~verbose() {
		close_log();
	}

	bool open_log(std::string fname) {
		log.open(fname.c_str());
		if (log.fail()) return false;
		else return (verbose_on_log = true);
	}

	void close_log() {
		log.close();
	}

	void set_silent() {
		verbose_on_screen = false;
	}

	void print(std::string s) {
		if (verbose_on_screen) std::cout << s << std::endl;
		if (verbose_on_log) log << s << std::endl;
	}

	void ctitle(std::string s) {
		if (verbose_on_screen) std::cout << std::endl << "\x1B[32m" << s <<  "\033[0m" << std::endl;
		if (verbose_on_log) log << std::endl << s << std::endl;
	}

	void title(std::string s) {
		if (verbose_on_screen) std::cout << std::endl << s << std::endl;
		if (verbose_on_log) log << std::endl << s << std::endl;
	}

	void bullet(std::string s) {
		if (verbose_on_screen) std::cout << "  * " << s << std::endl;
		if (verbose_on_log) log << "  * " << s << std::endl;
	}

	void warning(std::string s) {
		if (verbose_on_screen) std::cout << std::endl << "\x1B[33m" << "WARNING: " <<  "\033[0m" << s << std::endl;
		if (verbose_on_log) log << std::endl << "WARNING: " << s << std::endl;
	}

	void leave(std::string s) {
		if (verbose_on_screen) std::cout << std::endl << "\x1B[33m" << "EXITED: " <<  "\033[0m" << s << std::endl;
		if (verbose_on_log) log << std::endl << "EXITED: " << s << std::endl;
		exit(EXIT_SUCCESS);
	}

	void error(std::string s) {
		if (verbose_on_screen) std::cout << std::endl << "\x1B[31m" << "ERROR: " <<  "\033[0m" << s << std::endl;
		if (verbose_on_log) log << std::endl << "ERROR: " << s << std::endl;
		exit(EXIT_FAILURE);
	}

	void done(std::string s) {
		if (verbose_on_screen) std::cout << std::endl << "\x1B[32m" << "DONE: " <<  "\033[0m" << s << std::endl;
		if (verbose_on_log) log << std::endl << "DONE: " << s << std::endl;
		exit(EXIT_SUCCESS);
	}

	void wait(std::string s) {
		if (verbose_on_screen) {
			std::cout << s << " ...\r";
			std::cout.flush();
		}
	}

	void progress(std::string prefix, float percent) {
		if (verbose_on_screen) {
			int curr_percent = int(percent * 100.0);
			if (prev_percent > curr_percent) prev_percent = -1;
			if (curr_percent > prev_percent) {
				std::cout << prefix << " [" << curr_percent << "%]\r";
				std::cout.flush();
				prev_percent = curr_percent;
			}
		}
	}
};
#endif
