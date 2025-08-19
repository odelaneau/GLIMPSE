/*******************************************************************************
 * @file string_utils.h
 * @brief Utility class providing common string manipulation functions.
 *
 * Includes methods for splitting strings, checking numeric content,
 * converting numbers and vectors to strings, and extracting file names
 * and extensions from paths.
 *
 * @copyright Copyright (C) 2022-2023 Simone Rubinacci and Olivier Delaneau
 * @license MIT License
 ******************************************************************************/

#ifndef _STRING_UTILS_H
#define _STRING_UTILS_H

#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <regex>

/**
 * @class string_utils
 * @brief Collection of utility functions for string processing.
 *
 * This class offers:
 * - Splitting strings by a single character or set of separator characters.
 * - Checking if a string represents a numeric value.
 * - Converting numbers or vectors to strings with optional precision.
 * - Extracting file name, base name, and extension from file paths.
 */
class string_utils {
public:
	/// Default constructor
	string_utils () {};

	/// Default destructor
	~string_utils () {};

	/**
	 * @brief Split a string by a single character separator.
	 * @param str Input string to split.
	 * @param tokens Vector to store resulting tokens.
	 * @param sep Separator character.
	 * @param n_max_tokens Maximum number of tokens to split (default large).
	 * @return Number of tokens parsed.
	 *
	 * Splits the string by the given separator character, skipping
	 * consecutive separators, and trims trailing carriage returns.
	 */
	int split(const std::string & str, std::vector<std::string> & tokens, char sep, unsigned int n_max_tokens = 1000000) {
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
		if (!tokens.empty() && tokens.back().size() > 0 && tokens.back().back() == '\r')
			tokens.back().pop_back();
		return tokens.size();
	}

	/**
	 * @brief Split a string by any character in a set of separators.
	 * @param str Input string to split.
	 * @param tokens Vector to store resulting tokens.
	 * @param sep String of separator characters (default space and tab).
	 * @param n_max_tokens Maximum number of tokens to split (default large).
	 * @return Number of tokens parsed.
	 *
	 * Splits the string by any character contained in the separator string,
	 * skipping consecutive separators, and trims trailing carriage returns.
	 */
	int split(const std::string & str, std::vector<std::string> & tokens, std::string sep = " \t", unsigned int n_max_tokens = 1000000) {
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
		if (!tokens.empty() && tokens.back().size() > 0 && tokens.back().back() == '\r')
			tokens.back().pop_back();
		return tokens.size();
	}

	/**
	 * @brief Check if a string can be parsed as a numeric value.
	 * @param str Input string.
	 * @return True if string represents a number, false otherwise.
	 */
	bool numeric(std::string & str) {
		double n;
		std::istringstream in(str);
		return (in >> n) ? true : false;
	}

	/**
	 * @brief Convert a value to a string with optional fixed decimal precision.
	 * @tparam T Type of the value.
	 * @param n Value to convert.
	 * @param prec Number of decimal places (default no fixed precision).
	 * @return String representation of the value.
	 */
	template <class T>
	std::string str(T n, int prec = -1) {
		std::ostringstream ss(std::stringstream::out);
		if (prec >= 0) {
			ss << std::fixed << std::setprecision(prec);
		}
		ss << n;
		return ss.str();
	}

	/**
	 * @brief Convert a vector of values to a space-separated string with optional precision.
	 * @tparam T Type of elements in the vector.
	 * @param v Vector of values.
	 * @param prec Number of decimal places (default no fixed precision).
	 * @return String of space-separated values.
	 */
	template <class T>
	std::string str(std::vector<T> & v, int prec = -1) {
		std::ostringstream ss(std::stringstream::out);
		if (prec >= 0) {
			ss << std::fixed << std::setprecision(prec);
		}
		for (size_t e = 0; e < v.size(); e++) {
			if (e > 0) ss << " ";
			ss << v[e];
		}
		return ss.str();
	}

	/**
	 * @brief Extract the file name (last component) from a full file path.
	 * @param fullPath Full path string.
	 * @return File name without directory.
	 */
	std::string extract_file_name(const std::string& fullPath) {
		const size_t lastSlashIndex = fullPath.find_last_of("/\\");
		return fullPath.substr(lastSlashIndex + 1);
	}

	/**
	 * @brief Remove the extension from a file name.
	 * @param fileName File name string.
	 * @return File name without extension.
	 */
	std::string remove_ext(const std::string& fileName) {
		const size_t lastSlashIndex = fileName.find_last_of(".");
		return fileName.substr(0, lastSlashIndex);
	}

	/**
	 * @brief Get the base name (file name) from a path.
	 * @param path Full path string.
	 * @return Base name (last component of path).
	 */
	std::string base_name(std::string const & path) {
		return path.substr(path.find_last_of("/\\") + 1);
	}

	/**
	 * @brief Remove extension from a filename.
	 * @param filename File name string.
	 * @return Filename without extension or empty string if none.
	 */
	std::string remove_extension(std::string const & filename) {
		typename std::string::size_type const p(filename.find_last_of('.'));
		return p > 0 && p != std::string::npos ? filename.substr(0, p) : "";
	}

	/**
	 * @brief Get the extension of a file name if valid.
	 * @param filename File name string.
	 * @return Extension string if alphanumeric, else empty string.
	 */
	std::string get_extension(const std::string & filename) {
		auto position = filename.find_last_of('.');
		if (position == std::string::npos)
			return "";
		else {
			std::string extension(filename.substr(position + 1));
			if (std::regex_search(extension, std::regex("[^A-Za-z0-9]")))
				return "";
			else
				return extension;
		}
	}
};

#endif
