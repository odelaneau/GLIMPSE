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

#ifndef _CHECKSUM_UTILS_H
#define _CHECKSUM_UTILS_H

#include <vector>
#include <map>
#include <zlib.h>
#include <boost/archive/text_oarchive.hpp>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <functional>
#include <type_traits>

class checksum {
protected:
    unsigned long value;       // Current CRC32 checksum value
    std::string tmp_filename;  // (Unused) placeholder for filename (optional)
    bool new_data;             // (Unused) flag for new data (optional)

public:
    /**
     * @brief Constructor initializes checksum to initial CRC32 value.
     */
    checksum() : value(crc32(0L, Z_NULL, 0)), new_data(false) {}

    /**
     * @brief Get the current checksum value.
     * @return CRC32 checksum as unsigned long.
     */
    unsigned long get_value() const {
        return value;
    }

    /**
     * @brief Process raw data buffer of fundamental type T.
     * @tparam T Fundamental type (int, float, char, etc.)
     * @param buf Pointer to buffer to process
     * @param nbytes Number of bytes to process
     */
    template<typename T>
    void process_data(const T *buf, unsigned int nbytes)
    {
        static_assert(std::is_fundamental_v<T>, "process_data requires fundamental type");
        value = crc32(value, reinterpret_cast<const Bytef*>(buf), nbytes);
    }

    /**
     * @brief Process a fundamental type object.
     * @tparam T Fundamental type
     * @param obj Object reference
     */
    template<typename T>
    void process_data(const T &obj)
    {
        process_data(&obj, sizeof(T));
    }

    /**
     * @brief Process a std::string by feeding its raw data.
     * @param str Input string
     */
    void process_data(const std::string& str) {
        process_data(str.data(), static_cast<unsigned int>(str.size()));
    }

    /**
     * @brief Process a std::vector of fundamental types as raw contiguous bytes.
     * @tparam T Element type
     * @param vec Vector to process
     */
    template <typename T> 
    void process_data(const std::vector<T> &vec)
    {
        if (!vec.empty()) {
            process_data(vec.data(), static_cast<unsigned int>(sizeof(T) * vec.size()));
        }
    }

    /**
     * @brief Special handling for std::vector<bool> due to its non-contiguous storage.
     * Processes each boolean individually.
     * @param vec Vector of bool
     */
    void process_data(const std::vector<bool> &vec)
    {
        for(bool b : vec) {
            process_data(b);
        }
    }

    /**
     * @brief Process a vector of vectors recursively.
     * @tparam T Element type of inner vectors
     * @param vec Vector of vectors
     */
    template <typename T>
    void process_data(const std::vector< std::vector<T> > &vec)
    {
        for (const auto& inner_vec : vec) {
            process_data(inner_vec);
        }
    }

    /**
     * @brief Process a std::multimap by processing keys and values.
     * @tparam K Key type
     * @tparam V Value type
     * @param mmap Multimap to process
     */
    template<typename K, typename V>
    void process_data(const std::multimap<K, V> &mmap)
    {
        for (const auto& pair : mmap) {
            process_data(pair.first);
            process_data(pair.second);
        }
    }
};

#endif
