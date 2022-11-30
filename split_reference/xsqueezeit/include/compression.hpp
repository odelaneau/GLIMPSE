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

#ifndef __COMPRESSION_HPP__
#define __COMPRESSION_HPP__

#include <cstdint>
#include <vector>
#include <fstream>
#include <iostream>

typedef std::vector<size_t> ppa_t;

const uint32_t ENDIANNESS = 0xaabbccdd;
const uint32_t MAGIC = 0xfeed1767;
const uint32_t VERSION = 1;
const uint8_t  PLOIDY_DEFAULT = 2;

struct header_s {
    // "rsvd" fields are "reserved", unused for the moment and kept for future additions

    // 32 bytes
    uint32_t endianness = ENDIANNESS; // Endianness sanity check
    uint32_t first_magic = MAGIC;     // File format sanity check
    uint32_t version = VERSION;       // File format version
    uint8_t  ploidy = PLOIDY_DEFAULT; // Ploidy of samples encoded
    uint8_t  ind_bytes = 0;           // Number of bytes used to save indices
    uint8_t  aet_bytes = 0;           // Number of bytes used to save positions
    uint8_t  wah_bytes = 0;           // Number of bytes used for WAH
    union {
        uint8_t special_bitset = 0;
        struct {
            bool has_missing : 1;     // The input file had missing data
            bool non_uniform_phasing : 1; // The input file has phased / unphased mixed data
            bool default_phased : 1;
            uint8_t rsvd__1 : 5;
        };
    };
    union {
        uint8_t specific_bitset = 0;
        struct {
            bool iota_ppa : 1;        // Reset sort instead of saving permutation arrays
            bool no_sort : 1;         // Data is not permutated
            bool zstd : 1;
            uint8_t rsvd__2 : 5;
        };
    };
    uint8_t  rsvd_bs[2] = {0,};
    uint32_t rsvd_1[3] = {0,};

    // 64 bytes
    uint64_t hap_samples = 0;         // Number of haplotypes
    uint64_t num_variants = 0;        // Number of variants (total number of ALTs)
    uint32_t block_size = 0;          // DEPRECATED
    uint32_t number_of_blocks = 0;    // DEPRECATED
    uint32_t ss_rate = 0;             // Sub Sample rate of permutation arrays / reset sort rate
    // Offsets, positions of data in the binary file
    uint32_t number_of_ssas = 0;      // Number of sampled loci for random access = ceil(num_variants/ss_rate)
    uint32_t indices_offset = 0;      // Position in the binary file of WAH indices
    uint32_t ssas_offset = 0;         // Position in the binary file of sub sampled permutation arrays (if any)
    uint32_t wahs_offset = 0;         // Position in the binary file of WAH data
    uint32_t samples_offset = 0;      // Position in the binary file of samples (e.g., "NA12878", "HG00101")
    uint32_t indices_sparse_offset = 0; // Position in the binary file of indices for the sparse data
    uint32_t missing_offset = 0;
    uint32_t rearrangement_track_offset = 0; // Position in the binary file of the rearrangement track
    uint32_t sparse_offset = 0;       // Position in the binary file of the sparse data

    // 128 bytes
    uint32_t rare_threshold = 0;      // Threshold for the rearrangement track / sorting / wah vs sparse
    uint64_t xcf_entries = 0;         // Num entries in the BCF file (may be less than num_variants if multi-allelic)
    uint32_t phase_info_offset = 0;
    uint64_t num_samples = 0;
    uint8_t rsvd_3[104] = {0,};

    // 32 bytes
    uint32_t rsvd_4[3] = {0,};
    uint32_t sample_name_chksum = 0;  // Checksum unused for now
    uint32_t bcf_file_chksum = 0;     // Checksum unused for now
    uint32_t data_chksum = 0;         // Checksum unused for now
    uint32_t header_chksum = 0;       // Checksum unused for now
    uint32_t last_magic = MAGIC;      // Sanity check magic
} __attribute__((__packed__));

typedef struct header_s header_t;

static_assert(sizeof(header_t) == 256, "Header is not 256 bytes");

template<typename _ = uint32_t> /// @todo remove template, this is lazyness...
void print_header_info(const header_t& header) {
    std::cerr << "Version : " << header.version << std::endl;
    std::cerr << "Ploidy : " << (size_t)header.ploidy << std::endl;
    std::cerr << "Indice bytes : " << (size_t)header.ind_bytes << std::endl;
    std::cerr << "Sample id bytes : " << (size_t)header.aet_bytes << std::endl;
    std::cerr << "WAH bytes : " << (size_t)header.wah_bytes << std::endl;
    std::cerr << "--" << std::endl;
    //std::cerr << "Has missing : " << (header.has_missing ? "yes" : "no") << std::endl;
    //std::cerr << "Has non uniform phasing : " << (header.non_uniform_phasing ? "yes" : "no") << std::endl;
    //std::cerr << "Uses PPA's : " << (header.iota_ppa ? "no" : "yes" ) << std::endl;
    //std::cerr << "Is not sorted : " << (header.no_sort ? "yes" : "no" ) << std::endl;
    std::cerr << "Has a zstd compression layer : " << (header.zstd ? "yes" : "no") << std::endl;
    std::cerr << "--" << std::endl;
    std::cerr << "Haplotype samples  : " << header.hap_samples << std::endl;
    std::cerr << "Number of samples  : " << header.num_samples << std::endl;
    std::cerr << "Number of variants : " << header.num_variants << std::endl;
    std::cerr << "--" << std::endl;
    std::cerr << "VCF records : " << header.xcf_entries << std::endl;
    //std::cerr << "Permutation arrays  : " << header.wahs_offset - header.ssas_offset << " bytes" << std::endl;
    std::cerr << "GT Data WAH encoded : " << header.samples_offset - header.wahs_offset << " bytes" << std::endl;
}

template<typename _ = uint32_t> /// @todo remove template, this is lazyness...
int fill_header_from_file(const std::string filename, header_t& header) {
    std::fstream s(filename, s.binary | s.in);
    if (!s.is_open()) {
        std::cerr << "Failed to open file " << filename << std::endl;
        throw "Failed to open file";
    }

    // Read the header
    s.read((char *)(&header), sizeof(header_t));
    s.close();

    if ((header.first_magic != MAGIC) or (header.last_magic != MAGIC)) {
        return -1;
    }
    return 0;
}

#endif /* __COMPRESSION_HPP__ */