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
#include <containers/bitmatrix.h>
#include <zlib.h>

class checksum {
    protected:
        unsigned long value;
    
    public:
        checksum() {
            value = crc32(0L, Z_NULL, 0);
        }

        unsigned long value() {
            return value;
        }

        template<typename T>
        void process_data(const T *buf, unsigned int nbytes)
        {
            static_assert(std::is_fundamental_v<T>);
            value = crc32(value, static_cast<const uint8_t*>buf, nbytes);
        }

        template<typename T>
        void process_data(const T &obj)
        {
            process_data(&obj, sizeof(T));
        }

        template<>
        void process_data(const std::string str) {
            process_data(str.data(), str.size())
        }

        template <typename T> 
        void process_data(std::vector<T> const &vec)
        {
            process_data(vec.data(), sizeof(T)*vec.size())
        }

        //we need this because specialization of vector<bool> means that it is not guaranteed to be in contiguous memory, as other vectors are.
        template <>
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

// namespace checksum_utils {
//     template <typename T> 
//     unsigned long checksum_vector(unsigned long crc_init, std::vector<T> const &vec)
//     {
//         static_assert(std::is_fundamental_v<T>)
//         return crc32(checksum, &vec[0], sizeof(T)*vec.size())
//     }

//     //we need this because specialization of vector<bool> means that it is not guaranteed to be in contiguous memory, as other vectors are.
//     template <>
//     unsigned long checksum_vector(unsigned long crc_init, std::vector<bool> const &vec)
//     {
//         bitmatrix  vec_bitmatrix;
//         vec_bitmatrix.allocate(1,vec.size())
//         for (size_t i =0; i<vec.size(); i++) {
// 			char c = vec[i]? 1:0;
// 			vec_bitmatrix.set(1, i, c);
// 		}

//         return vec_bitmatrix.checksum(crc_init)
//     }

//     teamplate <typename T>
//     unsigned long checksum_vector(unsigned long crc_init, std::vector< std::vector<T> > const &vec)
//     {
//         unsigned long checksum = crc_init;
//         for (std::vector<T> inner_vec : vec) {
//             checksum = checksum_vector(checksum, inner_vec);
//         }

//         return checksum;
//     }

//     template <typename T>
//     unsigned long checksum_vector(std::vector<T> const &vec)
//     {
//         unsigned long crc_init = crc32(0L, Z_NULL, 0);

// 		return checksum_vector(crc_init, vec);
//     }

    
//     //we need this because specialization of vector<bool> means that it is not guaranteed to be in contiguous memory, as other vectors are.
//     unsigned long checksum_vector_bool(unsigned long crc_init, std::vector<bool> const &vec)
//     {
//         bitmatrix  vec_bitmatrix;
//         vec_bitmatrix.allocate(1,vec.size())
//         for (size_t i =0; i<vec.size(); i++) {
// 			char c = vec[i]? 1:0;
// 			vec_bitmatrix.set(1, i, c);
// 		}

//         return vec_bitmatrix.checksum(crc_init)
//     }

//     unsigned long checksum_vector_bool(std::vector<bool> const &vec)
//     {
//         unsigned long crc_init = crc32(0L, Z_NULL, 0);

// 		return checksum_vector_bool(crc_init, vec);
//     }
// }
#endif