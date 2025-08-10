/*******************************************************************************
 * @file verbose.h
 * @brief Simple verbosity and logging utility for console and file output.
 *
 * Provides methods to print messages, warnings, errors, progress, and formatted titles,
 * optionally both on screen (console) and into a log file.
 *
 * Supports colored output on terminal for warnings, errors, and done messages.
 *
 * Usage:
 * - Create an instance of verbose.
 * - Enable logging to a file with open_log().
 * - Control verbosity to screen with set_silent().
 * - Use print(), warning(), error(), progress(), etc. to output messages.
 *
 * @copyright Copyright (C) 2022-2023 Simone Rubinacci and Olivier Delaneau
 * @license MIT License
 ******************************************************************************/

#ifndef _VERBOSE_H
#define _VERBOSE_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

/**
 * @class verbose
 * @brief Utility for controlling program output verbosity and logging.
 *
 * Can print messages to standard output and/or to a log file.
 * Supports different message types: normal, warnings, errors, titles, progress.
 * Includes colored output for warnings, errors, done messages on supporting terminals.
 */
class verbose {
protected:
	/// Output log file stream
	std::ofstream log;

	/// Flag to enable verbose output on screen (console)
	bool verbose_on_screen;

	/// Flag to enable verbose output on log file
	bool verbose_on_log;

	/// Stores last printed progress percentage (rounded)
	int prev_percent;

public:
	/// Constructor: enables screen output by default, disables logging
	verbose() {
		verbose_on_screen = true;
		verbose_on_log = false;
		prev_percent = -1;
	}

	/// Destructor: closes log file if open
	~verbose() {
		close_log();
	}

	/**
	 * @brief Open a log file to write verbose output.
	 * @param fname Path to log file
	 * @return true if log file opened successfully, false otherwise
	 */
	bool open_log(std::string fname) {
		log.open(fname.c_str());
		if (log.fail()) return false;
		else return (verbose_on_log = true);
	}

	/// Close the log file stream if open
	void close_log() {
		log.close();
	}

	/// Disable output to screen
	void set_silent() {
		verbose_on_screen = false;
	}

	/**
	 * @brief Print a generic message to screen and/or log.
	 * @param s Message string to print
	 */
	void print(std::string s) {
		if (verbose_on_screen) std::cout << std::setprecision(16)  << s << std::endl;
		if (verbose_on_log) log << s << std::endl;
	}

	/**
	 * @brief Print a green colored title to screen and plain to log.
	 * @param s Title string
	 */
	void ctitle(std::string s) {
		if (verbose_on_screen) std::cout << std::endl << "\x1B[32m" << s <<  "\033[0m" << std::endl;
		if (verbose_on_log) log << std::endl << s << std::endl;
	}

	/**
	 * @brief Print a normal title to screen and log.
	 * @param s Title string
	 */
	void title(std::string s) {
		if (verbose_on_screen) std::cout << std::endl << s << std::endl;
		if (verbose_on_log) log << std::endl << s << std::endl;
	}

	/**
	 * @brief Print a bullet-point message to screen and log.
	 * @param s Message string
	 */
	void bullet(std::string s) {
		if (verbose_on_screen) std::cout << "  * " << s << std::endl;
		if (verbose_on_log) log << "  * " << s << std::endl;
	}

	/**
	 * @brief Print a yellow warning message to screen and log.
	 * @param s Warning message string
	 */
	void warning(std::string s) {
		if (verbose_on_screen) std::cout << std::endl << "\x1B[33m" << "WARNING: " <<  "\033[0m" << s << std::endl;
		if (verbose_on_log) log << std::endl << "WARNING: " << s << std::endl;
	}

	/**
	 * @brief Print a yellow exit message to screen and log, then exit program successfully.
	 * @param s Exit message string
	 */
	void leave(std::string s) {
		if (verbose_on_screen) std::cout << std::endl << "\x1B[33m" << "EXITED: " <<  "\033[0m" << s << std::endl;
		if (verbose_on_log) log << std::endl << "EXITED: " << s << std::endl;
		exit(EXIT_SUCCESS);
	}

	/**
	 * @brief Print a red error message to screen and log, then exit program with failure.
	 * @param s Error message string
	 */
	void error(std::string s) {
		if (verbose_on_screen) std::cout << std::endl << "\x1B[31m" << "ERROR: " <<  "\033[0m" << s << std::endl;
		if (verbose_on_log) log << std::endl << "ERROR: " << s << std::endl;
		exit(EXIT_FAILURE);
	}

	/**
	 * @brief Print a green done message to screen and log, then exit program successfully.
	 * @param s Done message string
	 */
	void done(std::string s) {
		if (verbose_on_screen) std::cout << std::endl << "\x1B[32m" << "DONE: " <<  "\033[0m" << s << std::endl;
		if (verbose_on_log) log << std::endl << "DONE: " << s << std::endl;
		exit(EXIT_SUCCESS);
	}

	/**
	 * @brief Print a waiting message (without newline), useful for showing progress.
	 * @param s Waiting message string
	 */
	void wait(std::string s) {
		if (verbose_on_screen) {
			std::cout << s << " ...\r";
			std::cout.flush();
		}
	}

	/**
	 * @brief Print a progress percentage update, throttled to every 10% increment.
	 * @param prefix Message prefix
	 * @param percent Progress as a float between 0 and 1
	 */
	void progress(const std::string prefix, const float percent) {
		if (verbose_on_screen) {
			int curr_percent = int(std::round(percent*10)*10.0f);
			if (prev_percent > curr_percent)
			{
				std::cout << prefix << "\t[" << 0.0 << "%]\r";
				std::cout.flush();
				prev_percent = 0.0;
				return;
			}
			if (curr_percent > prev_percent + 9.99) {
				std::cout << prefix << "\t[" << curr_percent << "%]\r";
				std::cout.flush();
				prev_percent = curr_percent;
			}
		}
	}

	/**
	 * @brief Print a generic progress message with no percentage.
	 * @param prefix Message prefix
	 */
	void progress(const std::string prefix) {
		if (verbose_on_screen) {
			std::cout << prefix << "\t...\r";
			std::cout.flush();
		}
	}
};

#endif
