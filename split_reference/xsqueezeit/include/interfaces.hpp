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

#ifndef __INTERFACES_HPP__
#define __INTERFACES_HPP__

#include <fstream>
#include <unordered_map>
#include <sys/mman.h>

#include "xcf.hpp"

// This could be changed and is not necessarily useful, it helps recognizing and checking the metadata
const uint32_t DICTIONARY_SIZE_SYMBOL = -1;

#define GEN_WRITE_DICTIONARY(type) \
template<typename T_K, typename T_V> \
inline size_t write_dictionary(std::fstream& ofs, const type<T_K, T_V>& dictionary) { \
    /* This metadata is always 32-bits */ \
    uint32_t key = DICTIONARY_SIZE_SYMBOL; \
    uint32_t val = dictionary.size(); \
    ofs.write(reinterpret_cast<const char*>(&key), sizeof(uint32_t)); \
    ofs.write(reinterpret_cast<const char*>(&val), sizeof(uint32_t)); \
     \
    size_t dictionary_offset = ofs.tellp(); \
     \
    for (const auto& kv : dictionary) { \
        ofs.write(reinterpret_cast<const char*>(&(kv.first)), sizeof(T_K)); \
        ofs.write(reinterpret_cast<const char*>(&(kv.second)), sizeof(T_V)); \
    } \
     \
    return dictionary_offset; \
}

// Generate the above with macros because I don't know how with templates
GEN_WRITE_DICTIONARY(std::map)
GEN_WRITE_DICTIONARY(std::unordered_map)

/// @todo this could check the size of the dictionnary, to make sure it didn't change in between
#define GEN_UPDATE_DICTIONARY_TYPE(type) \
template<typename T_K, typename T_V> \
inline void update_dictionary(std::fstream& ofs, const size_t& dictionary_pos, const type<T_K, T_V>& dictionary) { \
    const size_t old_pos = ofs.tellp(); \
    ofs.seekp(dictionary_pos, std::ios_base::beg); \
    for (const auto& kv : dictionary) { \
        ofs.write(reinterpret_cast<const char*>(&(kv.first)), sizeof(T_K)); \
        ofs.write(reinterpret_cast<const char*>(&(kv.second)), sizeof(T_V)); \
    } \
    ofs.seekp(old_pos, std::ios_base::beg); \
}

// Generate the above with macros because I don't know how with templates
GEN_UPDATE_DICTIONARY_TYPE(std::map)
GEN_UPDATE_DICTIONARY_TYPE(std::unordered_map)

template<typename MAP_T>
inline void read_dictionary(MAP_T& map, uint32_t* p) {
    map.clear();
    // Ignore first word
    (void)*p++;
    // Get dict size
    uint32_t size = *p++;
    for (uint32_t i = 0; i < size; ++i) {
        uint32_t key = *p++;
        uint32_t val = *p++;
        map[key] = val;
    }
    //return p;
}

class IBCFLineEncoder {
public:
    /**
     * @brief encode a BCF line internally
     * @param bcf_fri structure that allows access to BCF header, line record, etc.
     */
    virtual void encode_line(const bcf_file_reader_info_t& bcf_fri) = 0;

    virtual ~IBCFLineEncoder() {}
};

class BCFBlock {
public:
    BCFBlock(const size_t BLOCK_BCF_LINES) :
        effective_bcf_lines_in_block(0),
        BLOCK_BCF_LINES(BLOCK_BCF_LINES)
    {}

    size_t get_effective_bcf_lines_in_block() const {
        return effective_bcf_lines_in_block;
    }

protected:
    size_t effective_bcf_lines_in_block;
    const size_t BLOCK_BCF_LINES;
};

class IWritable {
public:
    /**
     * @brief write a completed writable to stream
     * @param ofs the output stream to write to
     */
    virtual void write_to_stream(std::fstream& ofs) = 0;

    /**
     * @brief gets the unique ID of the writable. If the writable is top level,
     * i.e., in the top level binary blocks of th XSI file. The ID must be
     * requested and officially assigned otherwise compatibility with the XSI
     * format cannot be guaranteed. Request an ID range through an issue here :
     * https://github.com/rwk-unil/xSqueezeIt/issues
     * Range 0-1024 is reserved.
     * Subsequent ranges up to 2^32-1 can be assigned per 1024 to any developer
     * / organization requesting them to extend the XSI format with their own
     * data / formats.
     *
     * This ensures retrocompatibility and compatibility of the extended files
     * with the original implementation even if the extended data / format is
     * not open-source or encrypted.
     * */
    virtual uint32_t get_id() const = 0;

    virtual ~IWritable() {}
};

class IWritableBCFLineEncoder : public IWritable, public IBCFLineEncoder {
public:
    // Abstract

    virtual ~IWritableBCFLineEncoder() {}
};

template <typename T_KEY = uint32_t, typename T_VAL = uint32_t>
class IBinaryBlock {
public:
    enum Dictionary_Keys : T_KEY {
        KEY_DICTIONNARY_SIZE = (T_KEY)-1,
        KEY_BCF_LINES = 0,
        KEY_GT_ENTRY = 256,
    };

    enum Dictionary_Vals : T_VAL {
        VAL_UNDEFINED = (T_VAL)-1,
    };

    virtual ~IBinaryBlock() {}

    void write_to_file(std::fstream& ofs, bool compressed, int compression_level) {
        int fd(0);
        auto ts = get_temporary_file(&fd);
        std::fstream& s = compressed ? ts.stream : ofs;

        size_t block_start_pos = s.tellp();
        size_t block_end_pos(0);
        size_t dictionary_pos(0);

        // Refresh dictionary
        //dictionary[KEY_BCF_LINES] = block_size;
        for (const auto& kv : writable_dictionary) {
            // Make sure the writables are also in the other dictionary
            dictionary[kv.first] = VAL_UNDEFINED;
        }

        // Block starts with dictionary size (key could be removed but helps recognize this)
        dictionary.erase(KEY_DICTIONNARY_SIZE); // Make sure the size is not in the dictionary
#if 0
        s.write(reinterpret_cast<const char*>(&DICTIONNARY_SIZE_KEY), sizeof(T_KEY));
        T_VAL dictionary_size = dictionary.size();
        s.write(reinterpret_cast<const char*>(&dictionary_size), sizeof(T_VAL));

        dictionary_pos = s.tellp();

        // Write dictionary (both combined)
        for (const auto& kv : dictionary) {
            s.write(reinterpret_cast<const char*>(&(kv.first)), sizeof(T_KEY));
            s.write(reinterpret_cast<const char*>(&(kv.second)), sizeof(T_VAL));
        }
#else
        dictionary_pos = write_dictionary(s, dictionary);
#endif

        // Write all the writables (extended entries)
        for (const auto& kv : writable_dictionary) {
            try {
                // Update the current offset in the dictionary
                dictionary.at(kv.first) = (T_VAL)((size_t)s.tellp()-block_start_pos);
            } catch (...) {
                std::cerr << "Entry is missing from dictionnary !" << std::endl;
                throw "Writable error";
            }
            // Write the writable at the current offset
            kv.second->write_to_stream(s);
        }

        block_end_pos = s.tellp();

        // Update the entries in dictionary on file
#if 0
        s.seekp(dictionary_pos, std::ios_base::beg);

        /// @todo here all entries are rewritten, this can be optimized if needed
        for (const auto& kv : dictionary) {
            s.write(reinterpret_cast<const char*>(&(kv.first)), sizeof(T_KEY));
            s.write(reinterpret_cast<const char*>(&(kv.second)), sizeof(T_VAL));
        }

        // Set stream back to end of file
        s.seekp(block_end_pos, std::io_bas::beg);
#else
        update_dictionary(s, dictionary_pos, dictionary);
#endif

        if (compressed) {
            size_t tmp_file_size = block_end_pos - block_start_pos;
            auto file_mmap = mmap(NULL, tmp_file_size, PROT_READ, MAP_SHARED, fd, 0);
            if (file_mmap == NULL) {
                std::cerr << "Failed to memory map file " << ts.filename << std::endl;
                throw "Failed to compress block";
            }

            compress_and_write(ofs, file_mmap, tmp_file_size, compression_level);

            munmap(file_mmap, tmp_file_size);
        }

        // Alignment padding...
        size_t mod_uint32 = size_t(ofs.tellp()) % sizeof(uint32_t);
        //std::cerr << "mod : " << mod_uint32 << std::endl;
        if (mod_uint32) {
            size_t padding = sizeof(uint32_t) - mod_uint32;
            for (size_t i = 0; i < padding; ++i) {
                //std::cerr << "A byte of padding was written" << std::endl;
                ofs.write("", sizeof(char));
            }
        }
        if (fd) {
            close(fd);
        }
        remove(ts.filename.c_str()); // Delete temp file
    }

    //size_t block_size;
    std::unordered_map<T_KEY, T_VAL> dictionary;
    std::unordered_map<T_KEY, std::shared_ptr<IWritable> > writable_dictionary;
    const T_KEY DICTIONNARY_SIZE_KEY = KEY_DICTIONNARY_SIZE;

protected:
    virtual void compress_and_write(std::fstream& ofs, void* data, size_t data_size, int compression_level) = 0;
};

template <typename TK, typename TV>
class BlockEntry {

};

// Todo move this guy
#include <iostream>
#include <zstd.h>
template<typename T_KEY, typename T_VAL> /// @todo maybe not template this
class BlockWithZstdCompressor : public IBinaryBlock<T_KEY, T_VAL> {
    typedef uint32_t T;
    static_assert(std::numeric_limits<T>::is_integer, "");

    void compress_and_write(std::fstream& ofs, void* data, size_t data_size, int compression_level) override {
        size_t output_buffer_size = data_size * 2;
        void *output_buffer = malloc(output_buffer_size);
        if (!output_buffer) {
            std::cerr << "Failed to allocate memory for output" << std::endl;
            throw "Failed to compress block";
        }

        auto result = ZSTD_compress(output_buffer, output_buffer_size, data, data_size, compression_level);
        if (ZSTD_isError(result)) {
            std::cerr << "Failed to compress file" << std::endl;
            std::cerr << "Error : " << ZSTD_getErrorName(result) << std::endl;
            throw "Failed to compress block";
        }

        if (std::numeric_limits<T>::max() < data_size or std::numeric_limits<T>::max() < result) {
            std::cerr << "Size too big to be encoded with " << typeid(T).name() << std::endl;
            throw "Failed to write compressed block";
        }
        T original_size = (T)data_size;
        T compressed_size = (T)result;

        ofs.write(reinterpret_cast<const char*>(&compressed_size), sizeof(T));
        ofs.write(reinterpret_cast<const char*>(&original_size), sizeof(T));
        ofs.write(reinterpret_cast<const char*>(output_buffer), compressed_size);

        free(output_buffer);
    }
};

/// @brief creates haploid arrangement from diploid arrangement
template<typename A_T>
std::vector<A_T> haploid_rearrangement_from_diploid(const std::vector<A_T>& a) {
    std::vector<A_T> a1(a.size() / 2);

    size_t a1_index = 0;
    // Go through diploid arrangement
    for (size_t i = 0; i < a.size(); ++i) {
        // Take the even arrangement (arbitrary)
        if ((a[i] & 1) == 0) {
            a1[a1_index] = a[i] / 2;
            a1_index++;
        }
    }

    return a1;
}

// Decompression related

#endif /* __INTERFACES_HPP__ */