/*******************************************************************************
 * Copyright (C) 2021 Rick Wertenbroek, University of Lausanne (UNIL),
 * University of Applied Sciences and Arts Western Switzerland (HES-SO),
 * School of Management and Engineering Vaud (HEIG-VD).
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

#include "constexpr.hpp"

#ifndef __WAH_HPP__
#define __WAH_HPP__

template <typename T>
void print_vector_(const std::vector<T>& v) {
    for (auto & e : v) {
        std::cout << e << " ";
    }
    std::cout << std::endl;
}

/*
 * Word Aligned Hybrid (WAH)
 * If you have no idea about WAH check out :
 * https://www.youtube.com/watch?v=BsZJ51JoYU0
 *
 * It is very similar to Run Length Encoding (RLE) as it also encodes runs of
 * repeating symbols but is much more flexible when the symbols wiggle (change
 * with high frequency). Also the words are always aligned to processor words
 * e.g., 8, 16, 32 bits etc. This allows for very efficient implementations of
 * operations see : Notes on design and implementation of compressed bit vectors
 * https://escholarship.org/uc/item/9zz3r6k7
 *
 *
 * The Word Aligned Hybrid (WAH) format used here is as follows :
 *
 * NOT USED ANYMORE BECAUSE WAH2 (see below) IS BETTER SUITED
 *
 * - The MSB bit indicates if the other bits represent a pattern or a counter
 * - - If 0 we have a pattern, copy the remaining bits
 * - - If 1 we have a counter, the remaning bits are the counter,
 *                             copy (#bits-1 * counter) zeroes
 *
 * E.g., for 8-bit WAH word :
 * WAH 0b0110'1010 (8 bits) => Decoded 0b110'1010 (7 bits)
 *        ^^^ ^^^^ pattern
 * WAH 0b1000'0011 (8 bits) => Decoded 0b000'0000 000'0000 000'0000 (21 bits)
 *        ^^^ ^^^^ counter
 *
 * MAX counter = 0b_111'1111 (127) which encodes 127*15 = 1905 consecutive zeroes
 *                 ^
 *                 |
 *                 +-- Bit that indicates a counter for repeating value
 *
 *
 *
 *
 *
 * The Word Aligned Hybrid 2 (WAH2) format used here is as follows :
 *
 * - The MSB bit indicates if the others bits represent a pattern or a counter
 * - - If 0 we have a pattern, copy the remaining bits
 * - - If 1 we have a counter, check the next bit, the remaining bits are the counter
 * - - - If it is 0, copy (#bits-1 * counter) zeroes
 * - - - If it is 1, copy (#bits-1 * counter) ones
 *
 * E.g., for 8-bit WAH2 word :
 * WAH 0b0110'1010 (8 bits) => Decoded 0b110'1010 (7 bits)
 *       |^^^ ^^^^ pattern
 * WAH 0b1000'0011 (8 bits) => Decoded 0b000'0000 000'0000 000'0000 (21 bits)
 *       ||^^ ^^^^ counter
 * WAH 0b1100'0011 (8 bits) => Decoded 0b111'1111 111'1111 111'1111 (21 bits)
 *       ||^^ ^^^^ counter
 *       |+-- Bit that indicates the repeating value
 *       |
 *       +--- Bit that indicates a counter
 *
 * MAX counter = 0b__11'1111 = (63) which encodes 63*15 = 945 repeating values
 *                 ^^
 *                 ||
 *                 |+-- Bit that indicates the repeating value
 *                 |
 *                 +--- Bit that indicates a counter
 *
 * Note : Since the values 0 is never used for a counter (you never encode a
 *        stretch of 0 times repeating values) the effective value of the counter
 *        could be taken as +1 (i.e., a 0 counter means 1 and 41 counter means 42)
 *        but it was chosen not to because this doesn't make much of a difference
 *        and makes the encoding more confusing.
 */
namespace wah {

    template<typename T>
    void print_wah2(T wah_word) {
        constexpr size_t WAH_BITS = sizeof(T)*8-1;
        // 0b1000'0000 for 8b
        constexpr T WAH_HIGH_BIT = 1 << WAH_BITS; // Solved at compile time
        // 0b0011'1111 for 8b
        constexpr T WAH_MAX_COUNTER = (WAH_HIGH_BIT>>1)-1;
        constexpr T WAH_COUNT_1_BIT = WAH_HIGH_BIT >> 1;

        if (WAH_HIGH_BIT & wah_word) {
            std::cout << size_t(wah_word & WAH_MAX_COUNTER) * size_t(WAH_BITS) << "*" << (wah_word & WAH_COUNT_1_BIT ? "1" : "0") << " ";
        } else {
            std::cout << std::bitset<WAH_BITS>(wah_word) << " ";
        }
    }

    enum class Wah2State {NONE, IN_WAH_WORD, IN_WAH_0_COUNTER, IN_WAH_1_COUNTER};

    typedef struct Wah2State_t {
        Wah2State state = Wah2State::NONE;
        size_t counter = 0;
    } Wah2State_t;

    template<typename T = uint16_t>
    inline size_t wah2_advance_pointer_count_ones(T*& wah_p, size_t size) {
        constexpr size_t WAH_BITS = sizeof(T)*8-1;
        constexpr T WAH_HIGH_BIT = 1 << WAH_BITS;
        constexpr T WAH_MAX_COUNTER = (WAH_HIGH_BIT>>1)-1;
        constexpr T WAH_COUNT_1_BIT = WAH_HIGH_BIT >> 1;
        T word;
        size_t count = 0;
        size_t bit_position = 0;
            while(bit_position < size) {
            word = *wah_p;
            if (word & WAH_HIGH_BIT) {
                auto counter = (word & WAH_MAX_COUNTER)*WAH_BITS;
                bit_position += counter;
                if (word & WAH_COUNT_1_BIT) {
                    count += counter;
                }
            } else {
                bit_position += WAH_BITS;
                //count += std::popcount(word); // C++20
                count += __builtin_popcount(word);
            }
            wah_p++;
        }
        return count;
    }

    template<typename T = uint16_t>
    inline size_t wah2_count_ones(const T*& wah_p, size_t size) {
        T* p = wah_p;
        return wah2_advance_pointer_count_ones(p, size);
    }

    template<typename T = uint16_t>
    inline void wah2_advance_pointer(T*& wah_p, size_t size) {
        constexpr size_t WAH_BITS = sizeof(T)*8-1;
        constexpr T WAH_HIGH_BIT = 1 << WAH_BITS;
        constexpr T WAH_MAX_COUNTER = (WAH_HIGH_BIT>>1)-1;
        T word;
        size_t bit_position = 0;
        while(bit_position < size) {
            word = *wah_p;
            if (word & WAH_HIGH_BIT) {
                bit_position += (word & WAH_MAX_COUNTER)*WAH_BITS;
            } else {
                bit_position += WAH_BITS;
            }
            wah_p++;
        }
    }

    // Note that bits should be padded to accomodate the last encoding (easiest is to add sizeof(T))
    template<typename T = uint16_t, bool DO_COUNT = false>
    static inline T* wah2_extract_template(T* wah_p, std::vector<bool>& bits, size_t size, size_t& count) {
        constexpr size_t WAH_BITS = sizeof(T)*8-1;
        constexpr T WAH_HIGH_BIT = 1 << WAH_BITS;
        constexpr T WAH_COUNT_1_BIT = WAH_HIGH_BIT >> 1;
        constexpr T WAH_MAX_COUNTER = (WAH_HIGH_BIT>>1)-1;
        if CONSTEXPR_IF (DO_COUNT) {
            count = 0;
        }

        T word;
        size_t bit_position = 0;
        while(bit_position < size) {
            word = *wah_p;
            if (word & WAH_HIGH_BIT) {
                const size_t stop = bit_position + (word & WAH_MAX_COUNTER)*WAH_BITS;
                if (word & WAH_COUNT_1_BIT) {
                    // Expand with ones
                    for (size_t _ = bit_position; _ < stop; ++_) {
                        bits[_] = 1;
                    }
                    if CONSTEXPR_IF (DO_COUNT) {
                        count += (word & WAH_MAX_COUNTER)*WAH_BITS;
                    }
                } else {
                    // Expand with zeroes
                    for (size_t _ = bit_position; _ < stop; ++_) {
                        bits[_] = 0;
                    }
                }
                bit_position = stop;
            } else {
                // Expand with value
                for (size_t _ = bit_position; _ < bit_position+WAH_BITS; ++_) {
                    bits[_] = word & 1; // May not be the most effective way
                    if CONSTEXPR_IF (DO_COUNT) {
                        count += word & 1;
                    }
                    word >>= 1;
                }
                bit_position += WAH_BITS;
            }
            wah_p++;
        }

        return wah_p;
    }

    // Note that bits should be padded to accomodate the last encoding (easiest is to add sizeof(T))
    template<typename T = uint16_t>
    inline T* wah2_extract(T* wah_p, std::vector<bool>& bits, size_t size) {
        size_t _;
        return wah2_extract_template<T>(wah_p, bits, size, _);
    }

    template<typename T = uint16_t>
    inline T* wah2_extract_count_ones(T* wah_p, std::vector<bool>& bits, size_t size, size_t& count) {
        return wah2_extract_template<T, true>(wah_p, bits, size, count);
    }

    // This is a stream based on the current wah pointer and state
    // If finished pulling and state is not none, it means the encoded vector is not
    // a multiple of WAH_BITS, just increment the pointer and set state to NONE to continue
    // with the next vector
    template<typename T = uint16_t>
    inline bool wah2_pull(T*& wah_p, Wah2State_t& state) {
        constexpr size_t WAH_BITS = sizeof(T)*8-1;
        // 0b1000'0000 for 8b
        constexpr T WAH_HIGH_BIT = 1 << WAH_BITS; // Solved at compile time
        // 0b0011'1111 for 8b
        constexpr T WAH_MAX_COUNTER = (WAH_HIGH_BIT>>1)-1;
        constexpr T WAH_COUNT_1_BIT = WAH_HIGH_BIT >> 1;

        const T wah_word = *wah_p;
        const size_t wah_counter_val = size_t(wah_word & (WAH_MAX_COUNTER))*size_t(WAH_BITS);
        bool res = false;

        if (state.state == Wah2State::NONE) {
            if (wah_word & WAH_HIGH_BIT) {
                if (wah_word & WAH_COUNT_1_BIT) {
                    state.state = Wah2State::IN_WAH_1_COUNTER;
                } else {
                    state.state = Wah2State::IN_WAH_0_COUNTER;
                }
            } else {
                state.state = Wah2State::IN_WAH_WORD;
            }
            state.counter = 0;
        }

        switch (state.state) {
            case Wah2State::IN_WAH_WORD :
                res = (wah_word >> state.counter++) & 1;
                if (state.counter == WAH_BITS) {
                    state.state = Wah2State::NONE;
                    wah_p++;
                }
                break;
            case Wah2State::IN_WAH_0_COUNTER :
                res = false;
                state.counter++;
                if (state.counter == wah_counter_val) {
                    state.state = Wah2State::NONE;
                    wah_p++;
                }
                break;
            case Wah2State::IN_WAH_1_COUNTER :
                res = true;
                state.counter++;
                if (state.counter == wah_counter_val) {
                    state.state = Wah2State::NONE;
                    wah_p++;
                }
                break;
            default :
                std::cerr << "Broken WAH 2 State" << std::endl;
                throw "Broken WAH 2 State";
                break;
        }
        return res;
    }

#if 0
    // Allows to get samples, requires initial permutations a, length in max columns, pointer to WAH2 data
    template <typename AET = size_t, typename WAH_T = uint16_t>
    class DecompressPointer {
    public:
        DecompressPointer(const std::vector<AET>& a_i, const size_t len, WAH_T* wah_p) : wah_p(wah_p), wah_origin(wah_p) {
            N = a_i.size();
            a_origin.resize(N);
            a.resize(N);
            b.resize(N);
            samples_sorted = false;
            y.resize(N+sizeof(AET));
            samples.resize(N);
            // This may be made copyless by passing a pointer for the first a
            std::copy(a_i.begin(), a_i.end(), a.begin());
            std::copy(a_i.begin(), a_i.end(), a_origin.begin());
            position = 0;
            length = len;
            // Default should already be set
            state.state = Wah2State::NONE;
            state.counter = 0;
        }

        // Will update a
        bool advance() { // Algorithm 1
            if (position >= length) {
                std::cerr << "Advance called without advancing" << std::endl;
                return false;
            } // Did not advance

            //print_vector_(a);

            AET u = 0;
            AET v = 0;
            wah_p = wah2_extract(wah_p, y, N);
            for (size_t i = 0; i < N; ++i) {
                //bool x = wah2_pull<WAH_T>(wah_p, state);
                if (y[i] == 0) {
                    a[u++] = a[i];
                } else {
                    b[v++] = a[i];
                }
                samples[a[i]] = y[i]; // Samples sorted by id
            }
            std::copy(b.begin(), b.begin()+v, a.begin()+u);
            samples_sorted = true;
            position++;

            // This is an edge case because number of samples may not be a multiple of WAH_BITS
            //if (state.state != Wah2State::NONE) {
            //    wah_p++;
            //    state.counter = 0;
            //    state.state = Wah2State::NONE;
            //}

            return true;
        }

        // Puts the pointer back at the start
        void reset() {
            wah_p = wah_origin;
            position = 0;
            samples_sorted = false;
            std::copy(a_origin.begin(), a_origin.end(), a.begin());
            state.state = Wah2State::NONE;
            state.counter = 0;
        }

        // Check if samples are ready (needs at least one advance)
        bool samples_ready() const { return samples_sorted; }
        // Offset of samples relative to start (samples need to be ready for it to be valid)
        size_t get_position() const { return position-1; }
        // Return the samples
        const std::vector<bool>& get_samples_at_position() const { return samples; }
        // Returns true if this pointer is done (went through all)
        bool done() const { return position >= length; }

    protected:
        AET N;
        size_t position;
        size_t length;
        std::vector<AET> a_origin;
        std::vector<AET> a;
        std::vector<AET> b;
        bool samples_sorted = false;
        std::vector<bool> y;
        std::vector<bool> samples;
        WAH_T* wah_p;
        WAH_T* wah_origin;
        Wah2State_t state;
    };
#endif

    /*
    * Strategy : Allow a maximum block size (this will also reduce the RAM burden)
    * The BCF file is read by variant for all samples (possible to stop early but not optimal)
    *
    * Input parameter : Block size in number of samples and number of variants
    * Let's start by the number of samples
    *
    * */

    /// @brief Word Aligned Hybrid encoding of a bit vector
    template <typename T = uint8_t>
    inline std::vector<T> wah_encode(std::vector<bool>& bits) {
        constexpr size_t WAH_BITS = sizeof(T)*8-1;
        /* Wahbit
         *             ,\
         *             \\\,_
         *              \` ,\
         *         __,.-" =__)
         *       ."        )
         *    ,_/   ,    \/\_
         *    \_|    )_-\ \_-`
         *jgs    `-----` `--`  By Joan Stark */
        // 0b1000'0000 for 8b
        constexpr T WAH_HIGH_BIT = 1 << WAH_BITS; // Solved at compile time
        // 0b0111'1111 for 8b
        constexpr T WAH_MAX_COUNTER = WAH_HIGH_BIT-1;
        constexpr T WAH_MAX_ENCODE = ~((T)(0));

        // Resize to have no problems accessing in the loop below (removes conditionals from loop)
        size_t BITS_WAH_SIZE = bits.size() / WAH_BITS;
        const size_t BITS_REM = bits.size() % WAH_BITS; // Hope compiler combines the divide above and this mod
        const size_t ORIGINAL_SIZE = bits.size();
        if (BITS_REM) {
            BITS_WAH_SIZE += 1;
            bits.resize(bits.size() + (WAH_BITS-BITS_REM), 0);
        }

        std::vector<T> wah;

        T not_set_counter = 0;
        size_t b = 0; // b = i*WAH_BITS+j // but counter reduces number of ops
        for (size_t i = 0; i < BITS_WAH_SIZE; ++i) { // Process loop
            T word = 0;

            // Scan WAH-BITS in bits (e.g., 7 for uint8_t, 31 for uint32_t)
            for (size_t j = 0; j < WAH_BITS; ++j) { // WAH word loop
                // b will always be in the vector because it has been resized above
                if (bits[b++]) { /// @todo Check if if need to/can be optimized
                    word |= WAH_HIGH_BIT;
                }
                word >>= 1;
            }
            // If bits found, encode them on N+1 bits
            if (word) {
                // Encode the number of previous blocks of WAH_BITS null bits
                if (not_set_counter) {
                    wah.push_back(WAH_HIGH_BIT | not_set_counter);
                    not_set_counter = 0;
                }
                wah.push_back(word);
            } else {
                // Handle possible counter overflow (unlikely)
                if (not_set_counter == WAH_MAX_COUNTER) {
                    wah.push_back(WAH_MAX_ENCODE); // = WAH_HIGH_BIT | not_set_counter = WAH_HIGH_BIT | WAH_MAX_COUNTER
                    not_set_counter = 0;
                }
                // Increment counter
                not_set_counter++;
            }
        }

        // Note that this could be deduced from known size, i.e., non encoded bits could be supposed zero by default
        if (not_set_counter) {
            wah.push_back(WAH_HIGH_BIT | not_set_counter);
        }

        if (BITS_REM) {
            bits.resize(ORIGINAL_SIZE);
        }

        return wah;
    }

    template <typename T = uint8_t>
    inline std::vector<bool> wah_decode(const std::vector<T>& wah) {
        constexpr size_t WAH_BITS = sizeof(T)*8-1;
        constexpr T WAH_HIGH_BIT = 1 << WAH_BITS;
        constexpr T WAH_MAX_COUNTER = WAH_HIGH_BIT-1;

        const size_t WAH_SIZE = wah.size();

        std::vector<bool> bits;
        for (size_t i = 0; i < WAH_SIZE; ++i) {
            T word = wah[i];
            if (word & WAH_HIGH_BIT) {
                // Expand with zeroes
                bits.resize(bits.size() + (wah[i] & WAH_MAX_COUNTER)*WAH_BITS, 0);
            } else {
                // Expand with value
                for (size_t j = 0; j < WAH_BITS; ++j) {
                    bits.push_back(word & 1); // May not be the most effective way
                    word >>= 1;
                }
            }
        }

        return bits;
    }

    /// @brief Word Aligned Hybrid encoding of a bit vector
    template <typename T = uint8_t> // This should be one of uint8-16-32-64_t
    inline std::vector<T> wah_encode2(std::vector<bool>& bits) {
        // This is a second version where a counter is used both for 0's and 1's (but a N-2 bit counter)
        constexpr size_t WAH_BITS = sizeof(T)*8-1;
        /* Wahbit
         *             ,\
         *             \\\,_
         *              \` ,\
         *         __,.-" =__)
         *       ."        )
         *    ,_/   ,    \/\_
         *    \_|    )_-\ \_-`
         *jgs    `-----` `--`  By Joan Stark */
        // 0b1000'0000 for 8b
        constexpr T WAH_HIGH_BIT = 1 << WAH_BITS; // Solved at compile time
        // 0b0011'1111 for 8b
        constexpr T WAH_MAX_COUNTER = (WAH_HIGH_BIT>>1)-1;
        constexpr T WAH_ALL_BITS_SET = T(~T(0)) & ~WAH_HIGH_BIT;
        constexpr T WAH_COUNT_1_BIT = WAH_HIGH_BIT >> 1;
        constexpr T WAH_MAX_ENCODE_0 = ~T(0) & ~WAH_COUNT_1_BIT;
        constexpr T WAH_MAX_ENCODE_1 = ~T(0);

        // Resize to have no problems accessing in the loop below (removes conditionals from loop)
        size_t BITS_WAH_SIZE = bits.size() / WAH_BITS;
        const size_t BITS_REM = bits.size() % WAH_BITS; // Hope compiler combines the divide above and this mod
        const size_t ORIGINAL_SIZE = bits.size();
        if (BITS_REM) {
            BITS_WAH_SIZE += 1;
            bits.resize(bits.size() + (WAH_BITS-BITS_REM), 0);
        }

        std::vector<T> wah;

        T not_set_counter = 0;
        T all_set_counter = 0;
        size_t b = 0; // b = i*WAH_BITS+j // but counter reduces number of ops
        for (size_t i = 0; i < BITS_WAH_SIZE; ++i) { // Process loop
            T word = 0;

            // Scan WAH-BITS in bits (e.g., 7 for uint8_t, 31 for uint32_t)
            for (size_t j = 0; j < WAH_BITS; ++j) { // WAH word loop
                // b will always be in the vector because it has been resized above
                if (bits[b++]) { /// @todo Check if if need to/can be optimized
                    word |= WAH_HIGH_BIT;
                }
                word >>= 1;
            }

            /// @note : Inner ifs should be mutually exclusive, counters should not be non zero at the same time ! (This could be reduced to a single counter and a boolean value to indicate what we count)

            // If no bits in word increment block counter
            if (word == (T)(0)) {
                // Encode the number of previous blocks of WAH_BITS set bits
                if (all_set_counter) {
                    wah.push_back(WAH_HIGH_BIT | WAH_COUNT_1_BIT | all_set_counter);
                    all_set_counter = 0;
                }
                // Handle possible counter overflow (unlikely)
                if (not_set_counter == WAH_MAX_COUNTER) {
                    wah.push_back(WAH_MAX_ENCODE_0);
                    not_set_counter = 0;
                }
                // Increment counter
                not_set_counter++;
            // If all bits are set in word increment block counter
            } else if (word == WAH_ALL_BITS_SET) {
                // Encode the number of previous blocks of WAH_BITS null bits
                if (not_set_counter) {
                    wah.push_back(WAH_HIGH_BIT | not_set_counter);
                    not_set_counter = 0;
                }
                // Handle possible counter overflow (unlikely)
                if (all_set_counter == WAH_MAX_COUNTER) {
                    wah.push_back(WAH_MAX_ENCODE_1);
                    all_set_counter = 0;
                }
                // Increment counter
                all_set_counter++;
            } else {
                if (all_set_counter) {
                    wah.push_back(WAH_HIGH_BIT | WAH_COUNT_1_BIT | all_set_counter);
                    all_set_counter = 0;
                }
                // Encode the number of previous blocks of WAH_BITS null bits
                if (not_set_counter) {
                    wah.push_back(WAH_HIGH_BIT | not_set_counter);
                    not_set_counter = 0;
                }
                wah.push_back(word);
            }
        }

        if (not_set_counter) {
            wah.push_back(WAH_HIGH_BIT | not_set_counter);
        }
        if (all_set_counter) {
            wah.push_back(WAH_HIGH_BIT | WAH_COUNT_1_BIT | all_set_counter);
        }

        if (BITS_REM) {
            bits.resize(ORIGINAL_SIZE);
        }

        return wah;
    }

    template <typename T = uint8_t>
    inline std::vector<bool> wah_decode2(const std::vector<T>& wah) {
        constexpr size_t WAH_BITS = sizeof(T)*8-1;
        constexpr T WAH_HIGH_BIT = 1 << WAH_BITS;
        constexpr T WAH_COUNT_1_BIT = WAH_HIGH_BIT >> 1;
        constexpr T WAH_MAX_COUNTER = (WAH_HIGH_BIT>>1)-1;

        const size_t WAH_SIZE = wah.size();

        std::vector<bool> bits;
        for (size_t i = 0; i < WAH_SIZE; ++i) {
            T word = wah[i];
            if (word & WAH_HIGH_BIT) {
                if (word & WAH_COUNT_1_BIT) {
                    // Expand with ones
                    bits.resize(bits.size() + (wah[i] & WAH_MAX_COUNTER)*WAH_BITS, 1);
                } else {
                    // Expand with zeroes
                    bits.resize(bits.size() + (wah[i] & WAH_MAX_COUNTER)*WAH_BITS, 0);
                }
            } else {
                // Expand with value
                for (size_t j = 0; j < WAH_BITS; ++j) {
                    bits.push_back(word & 1); // May not be the most effective way
                    word >>= 1;
                }
            }
        }

        return bits;
    }

    template <typename T = uint16_t>
    inline void process_wah_word(T& word, T& all_set_counter, T& not_set_counter, std::vector<T>& wah) {
        // This is a second version where a counter is used both for 0's and 1's (but a N-2 bit counter)
        constexpr size_t WAH_BITS = sizeof(T)*8-1;
        // 0b1000'0000 for 8b
        constexpr T WAH_HIGH_BIT = 1 << WAH_BITS; // Solved at compile time
        // 0b0011'1111 for 8b
        constexpr T WAH_MAX_COUNTER = (WAH_HIGH_BIT>>1)-1;
        constexpr T WAH_ALL_BITS_SET = T(~T(0)) & ~WAH_HIGH_BIT;
        constexpr T WAH_COUNT_1_BIT = WAH_HIGH_BIT >> 1;
        constexpr T WAH_MAX_ENCODE_0 = ~T(0) & ~WAH_COUNT_1_BIT;
        constexpr T WAH_MAX_ENCODE_1 = ~T(0);

        // If no bits in word increment block counter
        if (word == (T)(0)) {
            // Encode the number of previous blocks of WAH_BITS set bits
            if (all_set_counter) {
                wah.push_back(WAH_HIGH_BIT | WAH_COUNT_1_BIT | all_set_counter);
                all_set_counter = 0;
            }
            // Handle possible counter overflow (unlikely)
            if (not_set_counter == WAH_MAX_COUNTER) {
                wah.push_back(WAH_MAX_ENCODE_0);
                not_set_counter = 0;
            }
            // Increment counter
            not_set_counter++;
        // If all bits are set in word increment block counter
        } else if (word == WAH_ALL_BITS_SET) {
            // Encode the number of previous blocks of WAH_BITS null bits
            if (not_set_counter) {
                wah.push_back(WAH_HIGH_BIT | not_set_counter);
                not_set_counter = 0;
            }
            // Handle possible counter overflow (unlikely)
            if (all_set_counter == WAH_MAX_COUNTER) {
                wah.push_back(WAH_MAX_ENCODE_1);
                all_set_counter = 0;
            }
            // Increment counter
            all_set_counter++;
        } else {
            if (all_set_counter) {
                wah.push_back(WAH_HIGH_BIT | WAH_COUNT_1_BIT | all_set_counter);
                all_set_counter = 0;
            }
            // Encode the number of previous blocks of WAH_BITS null bits
            if (not_set_counter) {
                wah.push_back(WAH_HIGH_BIT | not_set_counter);
                not_set_counter = 0;
            }
            wah.push_back(word);
        }
    }

    struct DefaultPred {
        static inline bool check(const int32_t gt_arr_entry, const int32_t alt_allele) {
            return bcf_gt_allele(gt_arr_entry) == alt_allele;
        }
    };


    /**
     * @brief Copy-less non reordering wah encoder (counts number of satisfied pred)
     * */
    template <typename T = uint16_t, class Pred>
    inline std::vector<T> wah_encode2_with_size(int32_t* gt_array, const int32_t& alt_allele, const size_t size, uint32_t& pred_count) {
    // This is a second version where a counter is used both for 0's and 1's (but a N-2 bit counter)
        constexpr size_t WAH_BITS = sizeof(T)*8-1;
        // 0b1000'0000 for 8b
        constexpr T WAH_HIGH_BIT = 1 << WAH_BITS; // Solved at compile time
        constexpr T WAH_COUNT_1_BIT = WAH_HIGH_BIT >> 1;

        // Resize to have no problems accessing in the loop below (removes conditionals from loop)
        size_t BITS_WAH_SIZE = size / WAH_BITS;
        const size_t BITS_REM = size % WAH_BITS; // Hope compiler combines the divide above and this mod

        std::vector<T> wah; // Output

        T not_set_counter = 0;
        T all_set_counter = 0;
        size_t b = 0; // b = i*WAH_BITS+j // but counter reduces number of ops
        size_t pred_counter = 0;
        for (size_t i = 0; i < BITS_WAH_SIZE; ++i) { // Process loop
            T word = 0;

            // Scan WAH-BITS in bits (e.g., 7 for uint8_t, 31 for uint32_t)
            for (size_t j = 0; j < WAH_BITS; ++j) { // WAH word loop
                if (Pred::check_widx(b, gt_array[b], alt_allele)) {
                    word |= WAH_HIGH_BIT;
                    pred_counter++;
                }
                b++;
                word >>= 1;
            }

            process_wah_word(word, all_set_counter, not_set_counter, wah);
        }

        // Edge case of remaining bits (pad with 0's) // Is separate from above to avoid the check in the main loop
        if (BITS_REM) {
            T word = 0;
            for (size_t j = 0; j < WAH_BITS; ++j) {
                if (j < BITS_REM) {
                    if (Pred::check_widx(b, gt_array[b], alt_allele)) {
                        word |= WAH_HIGH_BIT;
                        pred_counter++;
                    }
                }
                b++;
                word >>= 1;
            }
            process_wah_word(word, all_set_counter, not_set_counter, wah);
        }

        // Push counters (they should be mutually exclusive, they cannot both be non zero)
        if (not_set_counter) {
            wah.push_back(WAH_HIGH_BIT | not_set_counter);
        }
        if (all_set_counter) {
            wah.push_back(WAH_HIGH_BIT | WAH_COUNT_1_BIT | all_set_counter);
        }

        pred_count = pred_counter;
        return wah;
    }

    /**
     * @brief Copy-less reordering wah encoder
     * */
    template <typename T = uint16_t, typename A = uint16_t, class Pred = DefaultPred>
    inline std::vector<T> wah_encode2_with_size(int32_t* gt_array, const int32_t& alt_allele, const std::vector<A>& a, const size_t size, uint32_t& alt_allele_count, bool& has_missing) {
    // This is a second version where a counter is used both for 0's and 1's (but a N-2 bit counter)
        constexpr size_t WAH_BITS = sizeof(T)*8-1;
        // 0b1000'0000 for 8b
        constexpr T WAH_HIGH_BIT = 1 << WAH_BITS; // Solved at compile time
        constexpr T WAH_COUNT_1_BIT = WAH_HIGH_BIT >> 1;

        // Resize to have no problems accessing in the loop below (removes conditionals from loop)
        size_t BITS_WAH_SIZE = size / WAH_BITS;
        const size_t BITS_REM = size % WAH_BITS; // Hope compiler combines the divide above and this mod

        std::vector<T> wah; // Output

        T not_set_counter = 0;
        T all_set_counter = 0;
        size_t b = 0; // b = i*WAH_BITS+j // but counter reduces number of ops
        size_t alt_allele_counter = 0;
        for (size_t i = 0; i < BITS_WAH_SIZE; ++i) { // Process loop
            T word = 0;

            // Scan WAH-BITS in bits (e.g., 7 for uint8_t, 31 for uint32_t)
            for (size_t j = 0; j < WAH_BITS; ++j) { // WAH word loop
                // Check if missing value is found, note that bcf_gt_missing is unphased missing !
                if (bcf_gt_allele(gt_array[a[b]]) == -1 /*bcf_gt_missing*/) {
                    //std::cout << "Missing value found" << std::endl;
                    has_missing = true;
                }
                //if (bcf_gt_allele(gt_array[a[b++]]) == alt_allele) {
                if (Pred::check(gt_array[a[b]], alt_allele)) {
                    word |= WAH_HIGH_BIT;
                    alt_allele_counter++;
                }
                word >>= 1;
                b++;
            }

            process_wah_word(word, all_set_counter, not_set_counter, wah);
        }

        // Edge case of remaining bits (pad with 0's) // Is separate from above to avoid the check in the main loop
        if (BITS_REM) {
            T word = 0;
            for (size_t j = 0; j < WAH_BITS; ++j) {
                if (j < BITS_REM) {
                    if (bcf_gt_allele(gt_array[a[b]]) == -1 /*bcf_gt_missing*/) {
                        //std::cout << "Missing value found" << std::endl;
                        has_missing = true;
                    }
                    //if (bcf_gt_allele(gt_array[a[b++]]) == alt_allele) {
                    if (Pred::check(gt_array[a[b]], alt_allele)) {
                        word |= WAH_HIGH_BIT;
                        alt_allele_counter++;
                    }
                }
                word >>= 1;
                b++;
            }
            process_wah_word(word, all_set_counter, not_set_counter, wah);
        }

        // Push counters (they should be mutually exclusive, they cannot both be non zero)
        if (not_set_counter) {
            wah.push_back(WAH_HIGH_BIT | not_set_counter);
        }
        if (all_set_counter) {
            wah.push_back(WAH_HIGH_BIT | WAH_COUNT_1_BIT | all_set_counter);
        }

        alt_allele_count = alt_allele_counter;
        //minor_allele_count = std::min(uint32_t(size) - alt_allele_counter, alt_allele_counter);
        return wah;
    }

    /**
     * @brief Copy-less reordering wah encoder
     * */
    template <typename T = uint16_t, typename A = uint16_t>
    inline std::vector<T> wah_encode2(int32_t* gt_array, const int32_t& alt_allele, const std::vector<A>& a, uint32_t& alt_allele_count, bool& has_missing) {
        return wah_encode2_with_size<T, A>(gt_array, alt_allele, a, a.size(), alt_allele_count, has_missing);
    }

    /**
     * @brief Encodes a number of values the are all the same as WAH
     * */
    template <typename T = uint16_t>
    inline std::vector<T> wah_encode2_all_same_value(size_t number, bool value) {
        // This is a second version where a counter is used both for 0's and 1's (but a N-2 bit counter)
        constexpr size_t WAH_BITS = sizeof(T)*8-1;
        // 0b1000'0000 for 8b
        constexpr T WAH_HIGH_BIT = 1 << WAH_BITS; // Solved at compile time
        constexpr T WAH_COUNT_1_BIT = WAH_HIGH_BIT >> 1;
        constexpr T WAH_MAX_COUNTER = (WAH_HIGH_BIT>>1)-1;

        std::vector<T> wah; // Output

        // Compute the number of WAH words needed for the output
        size_t WAH_WORDS = number / WAH_BITS;
        const size_t BITS_REM = number % WAH_BITS;
        if (BITS_REM) {
            WAH_WORDS++;
        }

        // Compute the number of full counters required
        const T FULL_COUNTER = WAH_HIGH_BIT | (value ? WAH_COUNT_1_BIT : 0) | WAH_MAX_COUNTER;
        const size_t NUMBER_OF_FULL_COUNTERS = WAH_WORDS / WAH_MAX_COUNTER;
        const size_t LAST_COUNTER = WAH_WORDS % WAH_MAX_COUNTER;

        // Add the full counters
        for (size_t i = 0; i < NUMBER_OF_FULL_COUNTERS; ++i) {
            wah.push_back(FULL_COUNTER);
        }

        // Add the last counter of WAH words
        wah.push_back(WAH_HIGH_BIT | LAST_COUNTER);

        return wah;
    }

} /* namespace wah */

#endif /* __WAH_HPP__ */