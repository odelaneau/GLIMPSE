/*******************************************************************************
 * @file timer.h
 * @brief Timer utility class for measuring elapsed time and formatting timestamps.
 *
 * Provides functionality for:
 * - Recording start and intermediate time points.
 * - Measuring elapsed time in milliseconds, microseconds, and seconds.
 * - Getting current date and time as formatted string.
 * - Formatting elapsed absolute time into human-readable days, hours, minutes, seconds.
 *
 * @copyright Copyright (C) 2022-2023 Simone Rubinacci and Olivier Delaneau
 * @license MIT License
 ******************************************************************************/

#ifndef _TIMER_H
#define _TIMER_H

#include <chrono>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <string>

/**
 * @class timer
 * @brief Measures elapsed time with high-resolution clock and provides formatting utilities.
 *
 * Usage:
 * - Instantiate an object to start the timer.
 * - Call clock() to reset intermediate timer.
 * - Use rel_time() or rel_time_micro() to get milliseconds/microseconds elapsed since last clock().
 * - Use abs_time() to get total seconds elapsed since creation.
 * - Use date() to get current system date/time as string.
 * - Use abs_time_display() for nicely formatted elapsed time string.
 */
class timer {
protected:
	/// Time point for start of timer
	std::chrono::time_point<std::chrono::high_resolution_clock> start_timing_clock;

	/// Time point for last clock() call
	std::chrono::time_point<std::chrono::high_resolution_clock> prev_timing_clock;

public:
	/// Constructor: initializes timer start and previous clock time to now
	timer () {
		start_timing_clock = std::chrono::high_resolution_clock::now();
	}

	/// Destructor
	~timer() {
	}

	/**
	 * @brief Resets the intermediate timer to current time.
	 */
	void clock() {
		prev_timing_clock = std::chrono::high_resolution_clock::now();
	}

	/**
	 * @brief Get elapsed time in milliseconds since last clock() call.
	 * @return Elapsed time in milliseconds.
	 */
	unsigned int rel_time() {
		return std::chrono::duration_cast<std::chrono::milliseconds>(
			std::chrono::high_resolution_clock::now() - prev_timing_clock).count();
	}

	/**
	 * @brief Get elapsed time in microseconds since last clock() call.
	 * @return Elapsed time in microseconds.
	 */
	unsigned int rel_time_micro() {
		return std::chrono::duration_cast<std::chrono::microseconds>(
			std::chrono::high_resolution_clock::now() - prev_timing_clock).count();
	}

	/**
	 * @brief Get total elapsed time in seconds since timer start.
	 * @return Total elapsed time in seconds.
	 */
	unsigned int abs_time() {
		return std::chrono::duration_cast<std::chrono::seconds>(
			std::chrono::high_resolution_clock::now() - start_timing_clock).count();
	}

	/**
	 * @brief Get the current system date and time as a formatted string.
	 * @return Date/time string in "dd/mm/yyyy - HH:MM:SS" format.
	 */
	std::string date() {
		auto now = std::chrono::system_clock::now();
		auto in_time_t = std::chrono::system_clock::to_time_t(now);
		std::stringstream ss;
	    ss << std::put_time(std::localtime(&in_time_t), "%d/%m/%Y - %X");
	    return ss.str();
	}

	/**
	 * @brief Format the total elapsed time since timer start as a human-readable string.
	 *
	 * Format includes days (d), hours, minutes, and seconds with appropriate pluralization.
	 * Examples: "1d02 hours 03 minutes 04 seconds", "03 minutes 12 seconds", "12 seconds"
	 *
	 * @return Formatted elapsed time string.
	 */
	std::string abs_time_display()
	{
		std::chrono::seconds input_seconds = std::chrono::duration_cast<std::chrono::seconds>(
			std::chrono::high_resolution_clock::now() - start_timing_clock);
	    using namespace std::chrono;
	    typedef duration<int, std::ratio<86400>> days;
	    auto d = duration_cast<days>(input_seconds);
	    input_seconds -= d;
	    auto h = duration_cast<hours>(input_seconds);
	    input_seconds -= h;
	    auto m = duration_cast<minutes>(input_seconds);
	    input_seconds -= m;
	    auto s = duration_cast<seconds>(input_seconds);

	    std::stringstream ss;
	    ss.fill('0');
	    if (d.count()) {
	        ss << d.count() << "d";
	    }
	    if (d.count() || h.count()) {
	        if (d.count()) { ss << std::setw(2); } // pad if day was present
	        ss << h.count() << " hour";
	        if (h.count() > 1) ss << "s";
	        ss << " ";
	    }
	    if (d.count() || h.count() || m.count()) {
	        if (d.count() || h.count()) { ss << std::setw(2); }
	        ss << m.count() << " minute";
	        if (m.count() > 1) ss << "s";
	        ss << " ";
	    }
	    if (d.count() || h.count() || m.count() || s.count()) {
	        if (d.count() || h.count() || m.count()) { ss << std::setw(2); }
	        ss << s.count() << " second";
	        if (s.count() > 1) ss << "s";
	    }

	    return ss.str();
	}
};

#endif
