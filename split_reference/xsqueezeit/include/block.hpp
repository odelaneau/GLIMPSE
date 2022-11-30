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

#ifndef __BLOCK_HPP__
#define __BLOCK_HPP__

#include <unordered_map>
#include <iostream>
#include <zstd.h>
#include <string>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

#include "fs.hpp"

#include "wah.hpp"
using namespace wah;

class Writable {
public:
    virtual void write_to_stream(std::fstream& os) = 0;

    virtual ~Writable() {};
};

template<typename T = uint32_t>
class SparseGtLine {
public:
    SparseGtLine() {}

    SparseGtLine(uint32_t index, int32_t* gt_array, int32_t ngt, int32_t sparse_allele) : index(index), sparse_allele(sparse_allele) {
        for (int32_t i = 0; i < ngt; ++i) {
            if (bcf_gt_allele(gt_array[i]) == sparse_allele) {
                sparse_encoding.push_back(i);
            }
        }
    }

    size_t index = 0;
    int32_t sparse_allele = 0;
    std::vector<T> sparse_encoding;
};

class Block {
public:

    Block(size_t BLOCK_SIZE) : BLOCK_SIZE(BLOCK_SIZE) {}

    enum Rearrangement : bool {
        NOT_REARRANGED = false,
        REARRANGED = true,
    };

    enum Dictionnary_Keys : uint32_t {
        KEY_UNUSED = (uint32_t)-1,
        KEY_SORT = 0,
        KEY_SELECT = 1,
        KEY_WAH = 2,
        KEY_SPARSE = 3,
        KEY_SMALL_BLOCK = 4,
        KEY_PLOIDY = 5,
        KEY_MISSING = 6,
        KEY_END_OF_VECTORS = 7
    };

    static void fill_dictionnary(void *block, std::unordered_map<uint32_t, uint32_t>& dict) {
        dict.clear();
        uint32_t* ptr = (uint32_t*)block;
        size_t _ = 0;
        // While there is a dictionnary entry
        while(ptr[_] != KEY_UNUSED and ptr[_+1] != KEY_UNUSED) {
            // Update dictionnary
            dict[ptr[_]] = ptr[_+1];
            _ += 2;
        }
    }

protected:
    const size_t BLOCK_SIZE;
    std::unordered_map<uint32_t, uint32_t> dictionnary;
    std::unordered_map<uint32_t, std::shared_ptr<Writable> > writable_dictionnary;
};

template <typename A_T, typename WAH_T = uint16_t>
class EncodedBlock : public Block {
public:
    EncodedBlock(size_t BLOCK_SIZE) : Block(BLOCK_SIZE) {}

    void reset() {
        rearrangement_track.clear();
        wahs.clear();
        sparse_lines.clear();
        dictionnary.clear();
        writable_dictionnary.clear();
    }

    void set_non_uniform_ploidy() {dictionnary[KEY_PLOIDY] = -1;}
    void set_block_ploidy(int32_t ploidy) {dictionnary[KEY_PLOIDY] = ploidy; }

    std::vector<bool> rearrangement_track;
    std::vector<std::vector<WAH_T> > wahs;
    std::vector<SparseGtLine<A_T> > sparse_lines;

    void write_to_file(std::fstream& os, bool compressed, int compression_level) {
        int fd = 0;
        auto ts = get_temporary_file(&fd);
        std::fstream& s = compressed ? ts.stream : os;

        size_t block_start_pos = 0;
        size_t block_end_pos = 0;

        block_start_pos = s.tellp();

        // Default values
        dictionnary[KEY_SORT] = KEY_UNUSED;
        dictionnary[KEY_SELECT] = KEY_UNUSED;
        dictionnary[KEY_WAH] = KEY_UNUSED;
        dictionnary[KEY_SPARSE] = KEY_UNUSED;
        if (this->rearrangement_track.size() < BLOCK_SIZE) {
            dictionnary[KEY_SMALL_BLOCK] = this->rearrangement_track.size();
        }
        for (const auto& kv : writable_dictionnary) {
            // Make sure the writables are in the other dictionnary
            dictionnary[kv.first] = KEY_UNUSED;
        }

        // Write dictionnary
        for (const auto& kv : dictionnary) {
            s.write(reinterpret_cast<const char*>(&(kv.first)), sizeof(uint32_t));
            s.write(reinterpret_cast<const char*>(&(kv.second)), sizeof(uint32_t));
        }
        // Write end of dictionnary
        uint32_t _ = KEY_UNUSED;
        s.write(reinterpret_cast<const char*>(&_), sizeof(uint32_t));
        s.write(reinterpret_cast<const char*>(&_), sizeof(uint32_t));

        size_t written_bytes = size_t(s.tellp());
        size_t total_bytes = size_t(s.tellp());

        // Write sort
        dictionnary[KEY_SORT] = (uint32_t)((size_t)s.tellp()-block_start_pos);
        auto sort = wah::wah_encode2<uint16_t>(this->rearrangement_track);
        s.write(reinterpret_cast<const char*>(sort.data()), sort.size() * sizeof(decltype(sort.back())));

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << "sort " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        // Write select
        dictionnary[KEY_SELECT] = dictionnary[KEY_SORT]; // Same is used
        // Write wah
        dictionnary[KEY_WAH] = (uint32_t)((size_t)s.tellp()-block_start_pos);
        for (const auto& wah : wahs) {
            s.write(reinterpret_cast<const char*>(wah.data()), wah.size() * sizeof(decltype(wah.back())));
        }

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << "WAH " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        // Write sparse
        dictionnary[KEY_SPARSE] = (uint32_t)((size_t)s.tellp()-block_start_pos);
        for (const auto& sparse_line : sparse_lines) {
            const auto& sparse = sparse_line.sparse_encoding;
            A_T number_of_positions = sparse.size();

            if (sparse_line.sparse_allele == 0) {
                //if (DEBUG_COMPRESSION) std::cerr << "NEGATED ";
                // Set the MSB Bit
                // This will always work as long as MAF is < 0.5
                // Do not set MAF to higher, that makes no sense because if will no longer be a MINOR ALLELE FREQUENCY
                /// @todo Check for this if user can set MAF
                number_of_positions |= (A_T)1 << (sizeof(A_T)*8-1);
            }
            s.write(reinterpret_cast<const char*>(&number_of_positions), sizeof(A_T));
            s.write(reinterpret_cast<const char*>(sparse.data()), sparse.size() * sizeof(decltype(sparse.back())));
        }

        std::cout << "Written " << sparse_lines.size() << " sparse lines" << std::endl;

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << "sparse " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        for (const auto& kv : writable_dictionnary) {
            // Add the current offset to the dictionnary
            dictionnary[kv.first] = (uint32_t)((size_t)s.tellp()-block_start_pos);
            // Write the writable
            kv.second->write_to_stream(s);
        }

        block_end_pos = s.tellp();

        // Rewrite the updated dictionnary
        s.seekp(block_start_pos, std::ios_base::beg);
        for (const auto& kv : dictionnary) {
            s.write(reinterpret_cast<const char*>(&(kv.first)), sizeof(uint32_t));
            s.write(reinterpret_cast<const char*>(&(kv.second)), sizeof(uint32_t));
        }
        // Set stream to end of block
        s.seekp(block_end_pos, std::ios_base::beg);

        // Funky as f...
        if (compressed) {
            size_t file_size = block_end_pos-block_start_pos;
            auto file_mmap = mmap(NULL, file_size, PROT_READ, MAP_SHARED, fd, 0);
            if (file_mmap == NULL) {
                std::cerr << "Failed to memory map file " << ts.filename << std::endl;
                throw "Failed to compress block";
            }

            size_t output_buffer_size = file_size * 2;
            void *output_buffer = malloc(output_buffer_size);
            if (!output_buffer) {
                std::cerr << "Failed to allocate memory for output" << std::endl;
                throw "Failed to compress block";
            }

            auto result = ZSTD_compress(output_buffer, output_buffer_size, file_mmap, file_size, compression_level);
            if (ZSTD_isError(result)) {
                std::cerr << "Failed to compress file" << std::endl;
                std::cerr << "Error : " << ZSTD_getErrorName(result) << std::endl;
                throw "Failed to compress block";
            }

            uint32_t size = (uint32_t)file_size;
            uint32_t compressed_size = (uint32_t)result;
            size_t start = os.tellp();
            (void)start;
            os.write(reinterpret_cast<const char*>(&compressed_size), sizeof(uint32_t));
            os.write(reinterpret_cast<const char*>(&size), sizeof(uint32_t));
            os.write(reinterpret_cast<const char*>(output_buffer), result);

            size_t stop = os.tellp();
            (void)stop;
            //std::cerr << "Start : " << start << " stop : " << stop << std::endl;
            //std::cerr << "compressed block of " << stop-start << " bytes written" << std::endl;
            free(output_buffer);
            munmap(file_mmap, file_size);
        } else {
            //std::cerr << "block of " << block_end_pos-block_start_pos << " bytes written" << std::endl;
        }

        // Alignment padding...
        size_t mod_uint32 = size_t(os.tellp()) % sizeof(uint32_t);
        //std::cerr << "mod : " << mod_uint32 << std::endl;
        if (mod_uint32) {
            size_t padding = sizeof(uint32_t) - mod_uint32;
            for (size_t i = 0; i < padding; ++i) {
                //std::cerr << "A byte of padding was written" << std::endl;
                os.write("", sizeof(char));
            }
        }
        if (fd) {
            close(fd);
        }
        remove(ts.filename.c_str()); // Delete temp file
    }
};

#if 0
class BinaryMatrixBlock : public Block {
public:
    BinaryMatrixBlock(size_t BLOCK_SIZE, const size_t N_HAPS) :
    Block(BLOCK_SIZE), N_HAPS(N_HAPS) {
    }

    template<typename A_T, typename WAH_T>
    load_from_memory(void *addr, bool compressed) {
        void *block;
        if (compressed) {
            size_t compressed_block_size = *(uint32_t*)addr;
            size_t uncompressed_block_size = *(uint32_t*)((uint8_t*)addr + sizeof(uint32_t));
            void *block_ptr = (uint8_t*)addr + sizeof(uint32_t)*2;
            //std::cerr << "Compressed block found, compressed size : " << compressed_block_size << ", uncompressed size : " << uncompressed_block_size << std::endl;
            if (block) {
                free(block);
                block = NULL;
            }
            block = malloc(uncompressed_block_size);
            if (!block) {
                std::cerr << "Failed to allocate memory to decompress block" << std::endl;
                throw "Failed to allocate memory";
            }
            auto result = ZSTD_decompress(block, uncompressed_block_size, block_ptr, compressed_block_size);
            if (ZSTD_isError(result)) {
                std::cerr << "Failed to decompress block" << std::endl;
                std::cerr << "Error : " << ZSTD_getErrorName(result) << std::endl;
                throw "Failed to decompress block";
            }
        } else {
            // Set block pointer
            block = addr;
        }

        // Fill dictionnary
        Block::fill_dictionnary(block, dictionnary);

        // Get pointers to the data
        WAH_T* wah_p = (WAH_T*)((uint8_t*)block + dictionnary[Block::Dictionnary_Keys::KEY_WAH]);
        A_T* sparse_p = (A_T*)((uint8_t*)block + dictionnary[Block::Dictionnary_Keys::KEY_SPARSE]);

        if (dictionnary.find(Block::Dictionnary_Keys::KEY_SMALL_BLOCK) != dictionnary.end()) {
            block_length = dictionnary[Block::Dictionnary_Keys::KEY_SMALL_BLOCK];
        } else {
            block_length = BLOCK_SIZE;
        }

        binary_matrix.clear();

        std::vector<bool> select_track(block_length + sizeof(WAH_T)*8);  // With some extra space because of WAH alignment
        WAH_T* st_p = (WAH_T*)((uint8_t*)block + dictionnary[Block::Dictionnary_Keys::KEY_SELECT]);
        wah2_extract<WAH_T>(st_p, select_track, block_length);
        select_track.resize(block_length);

        std::vector<bool> y(N_HAPS);
        SparseBoolVec sbv();
        for (size_t i = 0; i < block_length; ++i) {
            binary_matrix.push_back(std::vector<bool>(N_HAPS));
            if (select_track[i]) {
                wah_p = wah2_extract(wah_p, y, N_HAPS);
            } else {

            }
        }
        //if (select_track[0]) {
        //    this->wah_p = wah2_extract(this->wah_p, y, N_HAPS); // Extract current values
        //} else {
        //    this->sparse_p = sparse_extract(sparse_p);
        //}
    }
protected:
    const size_t N_HAPS;
    size_t block_length; // Effective block length
    std::vector<std::vector<bool> > binary_matrix;
}
#endif

#endif /* __BLOCK_HPP__ */