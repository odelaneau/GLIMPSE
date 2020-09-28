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

#ifndef _VERBOSE_H
#define _VERBOSE_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

using namespace std;

class verbose {
protected:
	ofstream log;
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

	bool open_log(string fname) {
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

	void print(string s) {
		if (verbose_on_screen) cout << s << endl;
		if (verbose_on_log) log << s << endl;
	}

	void ctitle(string s) {
		if (verbose_on_screen) cout << endl << "\x1B[32m" << s <<  "\033[0m" << endl;
		if (verbose_on_log) log << endl << s << endl;
	}

	void title(string s) {
		if (verbose_on_screen) cout << endl << s << endl;
		if (verbose_on_log) log << endl << s << endl;
	}

	void bullet(string s) {
		if (verbose_on_screen) cout << "  * " << s << endl;
		if (verbose_on_log) log << "  * " << s << endl;
	}

	void warning(string s) {
		if (verbose_on_screen) cout << endl << "\x1B[33m" << "WARNING: " <<  "\033[0m" << s << endl;
		if (verbose_on_log) log << endl << "WARNING: " << s << endl;
	}

	void leave(string s) {
		if (verbose_on_screen) cout << endl << "\x1B[33m" << "EXITED: " <<  "\033[0m" << s << endl;
		if (verbose_on_log) log << endl << "EXITED: " << s << endl;
		exit(EXIT_SUCCESS);
	}

	void error(string s) {
		if (verbose_on_screen) cout << endl << "\x1B[31m" << "ERROR: " <<  "\033[0m" << s << endl;
		if (verbose_on_log) log << endl << "ERROR: " << s << endl;
		exit(EXIT_FAILURE);
	}

	void done(string s) {
		if (verbose_on_screen) cout << endl << "\x1B[32m" << "DONE: " <<  "\033[0m" << s << endl;
		if (verbose_on_log) log << endl << "DONE: " << s << endl;
		exit(EXIT_SUCCESS);
	}

	void wait(string s) {
		if (verbose_on_screen) {
			cout << s << " ...\r";
			cout.flush();
		}
	}

	void progress(const string prefix, const float percent) {
		if (verbose_on_screen) {
			int curr_percent = int(percent * 100.0);
			if (prev_percent > curr_percent) prev_percent = -1;
			if (curr_percent > prev_percent) {
				cout << prefix << " [" << curr_percent << "%]\r";
				cout.flush();
				prev_percent = curr_percent;
			}
		}
	}
};
#endif
