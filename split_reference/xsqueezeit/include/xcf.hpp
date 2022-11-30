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

#ifndef __XCF_HPP__
#define __XCF_HPP__

#include <cstdint>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <regex>
#include <map>

#include "vcf.h"
#include "hts.h"
#include "synced_bcf_reader.h"

bool has_extension(const std::string& filename, const std::string& extension);

/**
 * @brief Creates csi index for given file
 *
 * @param filename file to index
 * @param n_threads optional parameters, number of threads for indexing, default 1
 * */
void create_index_file(std::string filename, int n_threads = 1);

typedef struct bcf_file_reader_info_t {
    bcf_srs_t* sr = nullptr; /* The BCF Synced reader */
    size_t n_samples = 0; /* The number of samples */
    size_t var_count = 0; /* The number of variant sites extracted */
    int* gt_arr = nullptr; /* Pointer on genotype data array */
    int size_gt_arr = 0; /* Size of above array as given by bcf_get_genotypes() */
    int ngt = 0;
    bcf1_t* line = nullptr; /* Current line pointer */
    size_t line_num = 0; /* Current line number */
    int line_alt_alleles_extracted = 0; /* For multi ALT alleles */

} bcf_file_reader_info_t;

void initialize_bcf_file_reader(bcf_file_reader_info_t& bcf_fri, const std::string& filename);

void destroy_bcf_file_reader(bcf_file_reader_info_t& bcf_fri);

void initialize_bcf_file_reader_with_region(bcf_file_reader_info_t& bcf_fri, const std::string& filename, const std::string& region, bool is_file = 0);

unsigned int bcf_next_line(bcf_file_reader_info_t& bcf_fri);

/// @todo check if better to return a std::vector or simply reuse one passed by reference
void extract_next_variant_and_update_bcf_sr(std::vector<bool>& samples, bcf_file_reader_info_t& bcf_fri);

/**
 * @brief extract_samples Returns a vector with the sample IDs
 * @param bcf_fri BCF File Reader Info (must already by initialized)
 * */
std::vector<std::string> extract_samples(const bcf_file_reader_info_t& bcf_fri);

/**
 * @brief writes a vector of strings to a given file
 *
 * The strings are newline separated in the output file
 *
 * @param v the vector of strings to be written
 * @param ofname the file to write the vector to
 * */
void string_vector_to_file(const std::vector<std::string>& v, const std::string& ofname);

/**
 * @brief reads a file and outputs a vector of strings that corresponds to the lines in the file
 * @param ifname the file to be read
 * @return the vector with the lines of the file
 * */
std::vector<std::string> string_vector_from_file(const std::string& ifname);

/**
 * @brief extract_samples Returns a vector with the sample IDs of file
 * @param fname File name of VCF / BCF file for sample ID extraction
 * */
std::vector<std::string> extract_samples(const std::string& fname);

/**
 * @brief remove_samples Removes all the samples from a VCF / BCF file and saves the result in a VCF / BCF file
 * @param ifname Input file name
 * @param ofname Output file name
 * @return The number of lines (entries)
 * */
size_t remove_samples(const std::string& ifname, const std::string& ofname);

std::string get_entry_from_bcf(const std::string& filename, const char *entry_key);

/**
 * @brief counts the number of entries (lines) in VCF / BCF file
 * @param ifname Input file name
 * @return number of entries
 * */
size_t count_entries(const std::string& ifname);

/**
 * @brief extract a matrix of bits for the genotype data, outer index is variant,
 *        inner index is samples (two bits per diploid individual)
 *
 * @param filename The input VCF/BCF file
 * */
std::vector<std::vector<bool> > extract_matrix(std::string filename);

/**
 * @brief utility function to check if two matrices are the same
 * */
template <typename T>
bool matrices_differ(std::vector<std::vector<T> >& m1, std::vector<std::vector<T> >& m2) {
    if (m1.size() != m2.size()) {
        std::cerr << "Different outer size" << std::endl;
        return true;
    } else {
        //std::cerr << "Outer sizes are the same : " << m1.size() << std::endl;
    }

    for (size_t i = 0; i < m1.size(); ++i) {
        if (m1[i].size() != m2[i].size()) {
            std::cerr << "Different inner size at " << i << std::endl;
            return true;
        }

        for (size_t j = 0; j < m1[i].size(); ++j) {
            if (m1[i][j] != m2[i][j]) {
                std::cerr << "Different value at " << i << ", " << j << std::endl;
                return true;
            }
        }
    }

    return false;
}

/**
 * @brief utility function to check if two VCF/BCF files have the same genotype dat
 *        matrix (only GT data checked, no variant, extra info etc.)
 * */
bool matrices_differ(std::string f1, std::string f2);

template<typename P = uint32_t, typename I = uint32_t>
std::map<P, I> create_map(const std::string& filename, const P threshold = 1000) {
    std::map<P, I> map;

    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, filename);

    int ret = bcf_hdr_set_samples(bcf_fri.sr->readers[0].header, NULL, 0 /* 0 is file 1 is list */); // All samples is "-", NULL is none
    if (ret < 0) {
        std::cerr << "Failed to set no samples in header for file " << filename << std::endl;
        throw "Failed to remove samples";
    }

    P last_mapped_position = 0;
    P position = 0;
    I index = 0;

    while(bcf_next_line(bcf_fri)) {
        position = bcf_fri.line->pos;

        if (last_mapped_position > position) {
            std::cerr << "Failed to map file because it is unordered" << std::endl;
            throw "Unordered input file";
        }

        if ((index == 0) or ((position - last_mapped_position) > threshold)) {
            last_mapped_position = position;
            map.insert({position, index});
        }

        index += bcf_fri.line->n_allele-1;
    }

    destroy_bcf_file_reader(bcf_fri);

    return map;
}

template<typename I = size_t>
struct line_index_t {
    I line = 0; // Line in vcf/bcf file
    I index = 0; // Index in bit matrix
};

/// @todo return line number and not only index (for region extraction)
/// @todo this is not optimal anyways
template<typename I = size_t, typename P = size_t>
struct line_index_t<I> find_index(const std::string& filename, const P position) {
    struct line_index_t<I> line_index = {0,0};

    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, filename);

    int ret = bcf_hdr_set_samples(bcf_fri.sr->readers[0].header, NULL, 0 /* 0 is file 1 is list */); // All samples is "-", NULL is none
    if (ret < 0) {
        std::cerr << "Failed to set no samples in header for file " << filename << std::endl;
        throw "Failed to remove samples";
    }

    while(bcf_next_line(bcf_fri)) {
        if (bcf_fri.line->pos+1 >= position) {
            destroy_bcf_file_reader(bcf_fri);
            return line_index;
        }
        line_index.index += bcf_fri.line->n_allele-1;
        line_index.line++;
    }

    destroy_bcf_file_reader(bcf_fri);
    throw "Position not found";
    return {};
}

std::string unique_id(bcf1_t* bcf_entry_p);

#include <unordered_map>

template <typename T = size_t>
std::unordered_map<std::string, line_index_t<T> > create_variant_map(std::string ifname) {
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, ifname);
    std::unordered_map<std::string, line_index_t<T> > map;
    bcf_fri.sr->max_unpack = BCF_UN_STR; // Minimal, but does not really impact perf because of header without samples below

    int ret = bcf_hdr_set_samples(bcf_fri.sr->readers[0].header, NULL, 0 /* 0 is file 1 is list */); // All samples is "-", NULL is none
    if (ret < 0) {
        std::cerr << "Failed to set no samples in header for file " << ifname << std::endl;
        throw "Failed to set samples";
    }

    T lines = 0;
    T matrix_pos = 0;

    while (bcf_next_line(bcf_fri)) {
        auto id = unique_id(bcf_fri.line);

        //std::cout << bcf_fri.line->d.id << std::endl;
        if (map.find(id) != map.end()) {
            std::cerr << "Duplicate ID : " << id << " at lines : " << map.find(id)->second.line << " and " << lines << std::endl;
        }
        map.insert({id, {lines, matrix_pos}});
        lines++;
        matrix_pos += bcf_fri.line->n_allele-1;
    }

    destroy_bcf_file_reader(bcf_fri);
    return map;
}

void unphase_xcf(const std::string& ifname, const std::string& ofname);

#include <random>
void unphase_xcf_random(const std::string& ifname, const std::string& ofname);

void sprinkle_missing_xcf(const std::string& ifname, const std::string& ofname);

std::vector<std::vector<bool> > extract_common_to_matrix(const std::string& ifname, const double MAF = 0.01);

/**
 * @brief This function replaces the samples by a single sample that represent the
 *        position of the now missing samples data in the binary matrix. This is allows
 *        to make a link between the vcf entry and the binary matrix where the sample
 *        data is actually stored. This allows to tabix/csi index the variant file and
 *        execute region queries on it. The resulting records then point to a position
 *        in the matrix, which can be used to decompress the missing data.
 *
 * */
size_t replace_samples_by_pos_in_binary_matrix(const std::string& ifname, const std::string& ofname, std::string xsi_fname = "", const bool new_version = false, const size_t BLOCK_LENGTH = 8192, const size_t BM_BLOCK_BITS = 15);

std::vector<std::vector<bool> > extract_phase_vectors(const std::string& ifname);

size_t compute_phase_switch_errors(const std::vector<bool>& testseq, const std::vector<bool>& refseq);

void compute_phase_switch_errors(const std::string& testFile, const std::string& refFile);

int32_t seek_default_phased(const std::string& filename, size_t limit = 3);

size_t seek_max_ploidy_from_first_entry(const std::string& filename);

bool file_has_no_samples(const std::string& filename);

bool file_has_no_entries(const std::string& filename);

#endif /* __XCF_HPP__ */