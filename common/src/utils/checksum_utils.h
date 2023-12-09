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

class checksum {
    protected:
        unsigned long value;
        std::string tmp_filename;
        bool new_data;
    
    public:
        checksum() {
            value = crc32(0L, Z_NULL, 0);
        }

        unsigned long get_value() {
            return value;
        }

        template<typename T>
        void process_data(const T *buf, unsigned int nbytes)
        {
            static_assert(std::is_fundamental_v<T>);
            value = crc32(value, reinterpret_cast<const Bytef*>(buf), nbytes);
        }

        template<typename T>
        void process_data(const T &obj)
        {
            process_data(&obj, sizeof(T));
        }

        void process_data(const std::string str) {
            process_data(str.data(), str.size());
        }

        template <typename T> 
        void process_data(std::vector<T> const &vec)
        {
            process_data(vec.data(), sizeof(T)*vec.size());
        }

        //we need this because specialization of vector<bool> means that it is not guaranteed to be in contiguous memory, as other vectors are.
        void process_data(std::vector<bool> const &vec)
        {
            for(bool b : vec) {
                process_data(b);
            }
        }

        template <typename T>
        void process_data(std::vector< std::vector<T> > const &vec)
        {
            for (const std::vector<T> inner_vec : vec) {
                process_data(inner_vec);
            }
        }

        template<typename K, typename V>
        void process_data(std::multimap<K, V> const &mmap)
        {
            for (const auto pair : mmap) {
                process_data(pair.first);
                process_data(pair.second);
            }
        }
        
};

#endif