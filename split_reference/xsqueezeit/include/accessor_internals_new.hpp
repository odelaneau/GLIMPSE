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
#ifndef __ACCESSOR_INTERNALS_NEW_HPP__
#define __ACCESSOR_INTERNALS_NEW_HPP__

#include <iostream>

#include <string>
#include <unordered_map>
#include "compression.hpp"
#include "xcf.hpp"
#include "gt_block.hpp"
#include "make_unique.hpp"

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include "fs.hpp"
#include "wah.hpp"

using namespace wah;

#include "accessor_internals.hpp" // For DecompressPointer, AccessorInternals
#include "interfaces.hpp"

/// @todo DecompressPointer NEW version
template <typename A_T = uint32_t, typename WAH_T = uint16_t>
class DecompressPointerGTBlock : /* public DecompressPointer<A_T, WAH_T>, */ public GTBlockDict, protected PBWTSorter {
public:
    DecompressPointerGTBlock(const header_t& header, void* block_p) :
        header(header), block_p(block_p), N_SAMPLES(header.num_samples), N_HAPS(N_SAMPLES ? N_SAMPLES*2 : header.hap_samples), /// @todo fix ploidy
        internal_binary_gt_line_position(0),
        internal_binary_weirdness_position(0),
        internal_binary_phase_position(0),
        y(N_HAPS+sizeof(WAH_T)*8-1, false),
        a(N_HAPS), b(N_HAPS),
        y_missing(N_HAPS+sizeof(WAH_T)*8-1, false),
        y_eovs(N_HAPS+sizeof(WAH_T)*8-1, false),
        y_phase(N_HAPS+sizeof(WAH_T)*8-1, false),
        a_weird(N_HAPS), b_weird(N_HAPS) {
        // Load dictionary
        read_dictionary(dictionary, (uint32_t*)block_p);

        //std::cerr << "[DEBUG] Dictionnary size : " << dictionary.size() << std::endl;

        bcf_lines_in_block = dictionary.at(KEY_BCF_LINES);
        binary_gt_lines_in_block = dictionary.at(KEY_BINARY_LINES);

        MAX_PLOIDY = dictionary.at(KEY_MAX_LINE_PLOIDY);
        if (MAX_PLOIDY == VAL_UNDEFINED) {
            std::cerr << "[DEBUG] Line max ploidy not found, setting to 2..." << std::endl;
            MAX_PLOIDY = 2;
        }

        DEFAULT_PHASING = dictionary.at(KEY_DEFAULT_PHASING);
        // Defalut ploidy should be 0 (unphased) or 1 (phased)
        if ((DEFAULT_PHASING != 1) or (DEFAULT_PHASING != 1)) {
            DEFAULT_PHASING = 0;
        }

        //std::cerr << "Created a new GTB decompress pointer with " << bcf_lines_in_block << " bcf lines and " << binary_gt_lines_in_block << " binary lines" << std::endl;

        // Load dim-1 structures
        fill_bool_vector_from_1d_dict_key(KEY_LINE_SELECT, binary_gt_line_is_wah, binary_gt_lines_in_block);
        if (!fill_bool_vector_from_1d_dict_key(KEY_LINE_SORT, binary_gt_line_is_sorting, binary_gt_lines_in_block)) {
            // By default only the wah lines are sorting
            binary_gt_line_is_sorting = binary_gt_line_is_wah;
        }

        // Check for weirdness
        block_has_weirdness = false;
        block_has_weirdness |= fill_bool_vector_from_1d_dict_key(KEY_LINE_MISSING, line_has_missing, binary_gt_lines_in_block);
        block_has_weirdness |= fill_bool_vector_from_1d_dict_key(KEY_LINE_END_OF_VECTORS, line_has_end_of_vector, binary_gt_lines_in_block);
        //std::cerr << "Block has weirdness : " << (block_has_weirdness ? "yes" : "no") << std::endl;
        //for (auto v : line_has_missing) { std::cerr << (v ? "1" : "0"); }
        //std::cerr << std::endl;

        block_has_non_uniform_phasing = fill_bool_vector_from_1d_dict_key(KEY_LINE_NON_UNIFORM_PHASING, line_has_non_uniform_phasing, binary_gt_lines_in_block);
        //if (block_has_non_uniform_phasing) { std::cerr << "Block has non uniform phasing" << std::endl; }

        // Handle fully haploid lines
        fill_bool_vector_from_1d_dict_key(KEY_LINE_HAPLOID, haploid_binary_gt_line, binary_gt_lines_in_block);
        if (!haploid_binary_gt_line.size()) { haploid_binary_gt_line.resize(binary_gt_lines_in_block, false); }

        /// @todo non default vector lengths (ploidy over 2)

        // Set access to 2D structures (don't decompress them unless needed)
        wah_origin_p = get_pointer_from_dict<WAH_T>(KEY_MATRIX_WAH);
        wah_p = wah_origin_p;

        sparse_origin_p = get_pointer_from_dict<A_T>(KEY_MATRIX_SPARSE);
        sparse_p = sparse_origin_p;

        // Set the weirdness matrices (will be nullptr if not present)
        missing_origin_p = get_pointer_from_dict<WAH_T>(KEY_MATRIX_MISSING);
        //if (missing_origin_p) std::cerr << "Missing 2D array found" << std::endl;
        missing_p = missing_origin_p;

        eovs_origin_p = get_pointer_from_dict<WAH_T>(KEY_MATRIX_END_OF_VECTORS);
        eovs_p = eovs_origin_p;

        non_uniform_phasing_origin_p = get_pointer_from_dict<WAH_T>(KEY_MATRIX_NON_UNIFORM_PHASING);
        non_uniform_phasing_p = non_uniform_phasing_origin_p;
        //if (non_uniform_phasing_origin_p) { std::cerr << "Block has non uniform phasing data" << std::endl; }

        std::iota(a.begin(), a.end(), 0);
        if (block_has_weirdness) {
            std::iota(a_weird.begin(), a_weird.end(), 0);
        }
    }
    virtual ~DecompressPointerGTBlock() {}

    /**
     * @brief Updates all internal structures to point to the requested binary gt entry
     * */
    inline void seek(const size_t position) /*override*/ {
        if (internal_binary_gt_line_position == position) {
            return;
        } else {
            if (internal_binary_gt_line_position > position) {
                std::cerr << "Slow backwards seek !" << std::endl;
                std::cerr << "Current position is : " << internal_binary_gt_line_position << std::endl;
                std::cerr << "Requested position is : " << position << std::endl;
                reset();
            }
            while (internal_binary_gt_line_position < position) {
                const size_t CURRENT_N_HAPS = ((haploid_binary_gt_line[internal_binary_gt_line_position]) ? N_SAMPLES : N_HAPS);
                if (binary_gt_line_is_wah[internal_binary_gt_line_position]) {
                    /// @todo
                    // Resize y based on the vector length
                    if (binary_gt_line_is_sorting[internal_binary_gt_line_position]) {
                        wah_p = wah2_extract(wah_p, y, CURRENT_N_HAPS);
                    } else {
                        /* reference advance */ wah2_advance_pointer(wah_p, CURRENT_N_HAPS);
                    }
                } else {
                    // Is sparse
                    if (binary_gt_line_is_sorting[internal_binary_gt_line_position]) {
                        sparse_p = sparse_extract(sparse_p);
                    } else {
                        sparse_p = sparse_advance_pointer(sparse_p);
                    }
                }
                update_a_if_needed();

                if (block_has_weirdness) {
                    // Advance weirdly
                    weirdness_advance(1, CURRENT_N_HAPS);
                }

                if (block_has_non_uniform_phasing) {
                    phase_advance(1, CURRENT_N_HAPS);
                }

                internal_binary_gt_line_position++;
            }
        }
    }

    inline size_t fill_genotype_array_advance(int32_t* gt_arr, size_t gt_arr_size, size_t n_alleles) {
        allele_counts.resize(n_alleles);
        size_t total_alt = 0;
        size_t n_missing = 0;
        size_t n_eovs = 0;

        const size_t CURRENT_N_HAPS = ((haploid_binary_gt_line[internal_binary_gt_line_position]) ? N_SAMPLES : N_HAPS);
        const size_t START_OFFSET = internal_binary_gt_line_position;

        // Set REF / first ALT
        if (!binary_gt_line_is_wah[internal_binary_gt_line_position]) { /* SPARSE */
            sparse_p = sparse_extract(sparse_p);
            int32_t default_gt = sparse_negated ? 1 : 0;
            int32_t sparse_gt = sparse_negated ? 0 : 1;

            for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
                gt_arr[i] = bcf_gt_unphased(default_gt) | ((i & 1) & DEFAULT_PHASING);
            }
            for (const auto& i : sparse) {
                //if constexpr (DEBUG_DECOMP) std::cerr << "Setting variant at " << i << std::endl;
                gt_arr[i] = bcf_gt_unphased(sparse_gt) | ((i & 1) & DEFAULT_PHASING);
            }
        } else { /* SORTED WAH */
            wah_p = wah2_extract_count_ones(wah_p, y, CURRENT_N_HAPS, ones);
            if (haploid_binary_gt_line[internal_binary_gt_line_position]) {
                auto a1 = haploid_rearrangement_from_diploid(a);
                for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
                    gt_arr[a1[i]] = bcf_gt_unphased(y[i]); // Haploids don't require phase bit
                }
            } else {
                for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
                    gt_arr[a[i]] = bcf_gt_unphased(y[i]) | ((a[i] & 1) & DEFAULT_PHASING);
                }
            }
        }

        allele_counts[1] = ones;
        total_alt = ones;
        update_a_if_needed();
        internal_binary_gt_line_position++;

        // If other ALTs (ALTs are 1 indexed, because 0 is REF)
        for (size_t alt_allele = 2; alt_allele < n_alleles; ++alt_allele) {
            if (!binary_gt_line_is_wah[internal_binary_gt_line_position]) { /* SPARSE */
                sparse_p = sparse_extract(sparse_p);
                if (sparse_negated) { // There can only be one negated because must be more than all others combined
                    // All non set positions are now filled
                    for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
                        // Only overwrite refs
                        if (bcf_gt_allele(gt_arr[i]) == 0) {
                            gt_arr[i] = bcf_gt_unphased(alt_allele) | ((i & 1) & DEFAULT_PHASING);
                        }
                    }
                    for (const auto& i : sparse) {
                        // Restore overwritten refs
                        if (bcf_gt_allele(gt_arr[i]) == (int)alt_allele) {
                            gt_arr[i] = bcf_gt_unphased(0) | ((i & 1) & DEFAULT_PHASING);
                        }
                    }
                } else {
                    // Fill normally
                    for (const auto& i : sparse) {
                        gt_arr[i] = bcf_gt_unphased(alt_allele) | ((i & 1) & DEFAULT_PHASING);
                    }
                }
            } else { /* SORTED WAH */
                wah_p = wah2_extract_count_ones(wah_p, y, CURRENT_N_HAPS, ones);
                if (haploid_binary_gt_line[internal_binary_gt_line_position]) {
                    auto a1 = haploid_rearrangement_from_diploid(a);
                    for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
                        if (y[i]) {
                            gt_arr[a1[i]] = bcf_gt_unphased(y[i]); // Haploids don't require phase bit
                        }
                    }
                } else {
                    for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
                        if (y[i]) {
                            gt_arr[a[i]] = bcf_gt_unphased(alt_allele) | ((a[i] & 1) & DEFAULT_PHASING); /// @todo Phase
                        }
                    }
                }
            }
            allele_counts[alt_allele] = ones;
            total_alt += ones;
            update_a_if_needed();
            internal_binary_gt_line_position++;
        }

        //std::cerr << "[DEBUG] Start offset : " << START_OFFSET << " internal gt line " << internal_binary_gt_line_position << std::endl;
        //for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
        //    std::cerr << gt_arr[i] << " ";
        //}
        //std::cerr << std::endl;

        // Apply missing, eovs
        if (block_has_weirdness) {
            if (START_OFFSET != internal_binary_weirdness_position) {
                std::cerr << "Block decompression corruption on missing or end of vectors" << std::endl;
            }
            // weirdness is either missing or end of vector

            if (line_has_missing.size() and line_has_missing[START_OFFSET]) {
                // Fill missing without advance
                //std::cerr << "Extracting missing from wah p" << std::endl;
                /*missing_p =*/ (void) wah2_extract_count_ones(missing_p, y_missing, CURRENT_N_HAPS, n_missing);
                for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
                    if (y_missing[i]) {
                        const auto index = a_weird[i];
                        //std::cerr << "Filling a missing position" << std::endl;
                        gt_arr[index] = bcf_gt_missing | ((index & 1) & DEFAULT_PHASING);
                    }
                }
            }
            if (line_has_end_of_vector.size() and line_has_end_of_vector[START_OFFSET]) {
                // Fill eovs without advance
                /*eovs_p =*/ (void) wah2_extract_count_ones(eovs_p, y_eovs, CURRENT_N_HAPS, n_eovs);
                for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
                    if (y_eovs[i]) {
                        const auto index = a_weird[i];
                        gt_arr[index] = bcf_int32_vector_end;
                    }
                }
            }

            //std::cerr << "[DEBUG] Weirdness a : ";
            //for (auto& e : a_weird) std::cerr << e << " ";
            //std::cerr << std::endl;

            // PBWT weirdness advance
            // This is not the most optimal because double extraction, but weirdness is an edge case so we don't care
            weirdness_advance(n_alleles-1, CURRENT_N_HAPS);
        }

        // Apply phase info
        if (block_has_non_uniform_phasing) {
            if (START_OFFSET != internal_binary_phase_position) {
                std::cerr << "Block decompression corruption on phase information" << std::endl;
            }

            if (line_has_non_uniform_phasing.size() and line_has_non_uniform_phasing[START_OFFSET]) {
                (void) wah2_extract(non_uniform_phasing_p, y_phase, CURRENT_N_HAPS);
                for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
                    if (y_phase[i]) {
                        //std::cerr << "Toggling phase bit" << std::endl;
                        // Toggle phase bit
                        /// @todo only works for PLOIDY 1 and 2
                        if (gt_arr[i] != bcf_int32_vector_end) {
                            gt_arr[i] ^= (i & 1); // if non default phase toggle bit
                        } // Don't phase end of vector !
                    }
                }
            }

            phase_advance(n_alleles-1, CURRENT_N_HAPS);
        }

        //for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
        //    std::cerr << gt_arr[i] << " ";
        //}
        //std::cerr << std::endl;

        allele_counts[0] = CURRENT_N_HAPS - (total_alt + n_missing + n_eovs);

        return CURRENT_N_HAPS;
    }

    void reset() {
        // Reset internal structures
        std::iota(a.begin(), a.end(), 0);
        internal_binary_gt_line_position = 0;
        wah_p = wah_origin_p;
        sparse_p = sparse_origin_p;

        if (block_has_weirdness) {
            std::iota(a_weird.begin(), a_weird.end(), 0);
            internal_binary_weirdness_position = 0;
            missing_p = missing_origin_p;
            eovs_p = eovs_origin_p;
        }
        if (block_has_non_uniform_phasing) {
            internal_binary_phase_position = 0;
            non_uniform_phasing_p = non_uniform_phasing_origin_p;
        }
    }

    inline void fill_allele_counts_advance(const size_t n_alleles) {
        allele_counts.resize(n_alleles);
        size_t total_alt = 0;

        const size_t CURRENT_N_HAPS = ((haploid_binary_gt_line[internal_binary_gt_line_position]) ? N_SAMPLES : N_HAPS);

        for (size_t alt_allele = 1; alt_allele < n_alleles; ++alt_allele) {
            if (binary_gt_line_is_wah[internal_binary_gt_line_position]) {
                /// @todo
                // Resize y based on the vector length
                if (binary_gt_line_is_sorting[internal_binary_gt_line_position]) {
                    wah_p = wah2_extract_count_ones(wah_p, y, CURRENT_N_HAPS, ones);
                } else {
                    /* reference advance */ ones = wah2_advance_pointer_count_ones(wah_p, CURRENT_N_HAPS);
                }
            } else {
                // Is sparse (both methods count ones)
                if (binary_gt_line_is_sorting[internal_binary_gt_line_position]) {
                    sparse_p = sparse_extract(sparse_p);
                } else {
                    sparse_p = sparse_advance_pointer(sparse_p);
                }
            }
            update_a_if_needed();
            internal_binary_gt_line_position++;

            allele_counts[alt_allele] = ones;
            total_alt += ones;
        }

        allele_counts[0] = CURRENT_N_HAPS - total_alt; // - total missing ?
    }

    inline const std::vector<size_t>& get_allele_count_ref() const {
        return allele_counts;
    }

protected:
    inline void weirdness_advance(const size_t STEPS, const size_t CURRENT_N_HAPS) {
        // Update pointers and PBWT weirdness
        for (size_t i = 0; i < STEPS; ++i) {
            bool current_line_has_missing = false;
            bool current_line_has_eovs = false;

            if (line_has_missing.size() and line_has_missing[internal_binary_weirdness_position]) {
                current_line_has_missing = true;
                // Advance missing pointer
                missing_p = wah2_extract(missing_p, y_missing, CURRENT_N_HAPS);
            }
            if (line_has_end_of_vector.size() and line_has_end_of_vector[internal_binary_weirdness_position]) {
                current_line_has_eovs = true;
                // Fill eovs;
                eovs_p = wah2_extract(eovs_p, y_eovs, CURRENT_N_HAPS);
            }

            // Update PBWT weirdness
            /// @todo refactor all this
            if (current_line_has_missing and current_line_has_eovs) {
                // PBWT on both
                if (haploid_binary_gt_line[internal_binary_weirdness_position]) {
                    //bool_pbwt_sort_two<A_T, 2>(a_weird, b_weird, y_missing, y_eovs, N_HAPS);
                }
                else {
                    bool_pbwt_sort_two<A_T>(a_weird, b_weird, y_missing, y_eovs, N_HAPS);
                }
            } else if (current_line_has_missing) {
                // PBWT on missing
                if (haploid_binary_gt_line[internal_binary_weirdness_position]) {
                    //bool_pbwt_sort<A_T, 2>(a_weird, b_weird, y_missing, N_HAPS);
                }
                else {
                    bool_pbwt_sort<A_T>(a_weird, b_weird, y_missing, N_HAPS);
                }

            } else if (current_line_has_eovs) {
                // PBWT on eovs
                if (haploid_binary_gt_line[internal_binary_weirdness_position]) {
                    //bool_pbwt_sort<A_T, 2>(a_weird, b_weird, y_eovs, N_HAPS);
                }
                else {
                    bool_pbwt_sort<A_T>(a_weird, b_weird, y_eovs, N_HAPS);
                }
            }

            internal_binary_weirdness_position++;
        }
    }

    inline void phase_advance(const size_t STEPS, const size_t CURRENT_N_HAPS) {
        for (size_t i = 0; i < STEPS; ++i) {
            if (line_has_non_uniform_phasing.size() and line_has_non_uniform_phasing[internal_binary_phase_position]) {
                wah2_advance_pointer(non_uniform_phasing_p, CURRENT_N_HAPS);
            }
            internal_binary_phase_position++;
        }
    }

    template<const size_t V_LEN_RATIO = 1>
    inline void private_pbwt_sort() {
        static_assert(V_LEN_RATIO <= 2, "Is not meant to be");
        if CONSTEXPR_IF (V_LEN_RATIO == 1) {
            bool_pbwt_sort<A_T>(a, b, y, N_HAPS);
        } else if CONSTEXPR_IF (V_LEN_RATIO == 2) {
            auto a1 = haploid_rearrangement_from_diploid(a);
            /// @todo find a better solution ?
            std::vector<bool> x(N_SAMPLES);
            for (size_t i = 0; i < N_SAMPLES; ++i) {
                x[a1[i]] = y[i];
            }
            size_t u = 0;
            size_t v = 0;
            for (size_t i = 0; i < N_SAMPLES*V_LEN_RATIO; ++i) {
                if (x[a[i]/V_LEN_RATIO] == 0) {
                    a[u++] = a[i];
                } else {
                    b[v++] = a[i];
                }
            }
            std::copy(b.begin(), b.begin()+v, a.begin()+u);
        }
    }

    inline void update_a_if_needed() {
        // Extracted line is used to sort
        if (binary_gt_line_is_sorting[internal_binary_gt_line_position]) {
            //std::cerr << "[DEBUG] a : ";
            //for (auto& e : a) std::cerr << e << " ";
            //std::cerr << std::endl;
            //std::cerr << "[DEBUG] y : ";
            //for (auto e : y) std::cerr << e << " ";
            //std::cerr << std::endl;
            if (haploid_binary_gt_line[internal_binary_gt_line_position]) {
                //std::cerr << "Sort with VLENRATIO2 for line " << internal_binary_gt_line_position << std::endl;
                private_pbwt_sort<2>();
            } else {
                private_pbwt_sort<1>();
            }
        }
    }

    inline bool fill_bool_vector_from_1d_dict_key(enum Dictionary_Keys key, std::vector<bool>& v, const size_t size) {
        if (dictionary.find(key) != dictionary.end()) {
            if (dictionary[key] != VAL_UNDEFINED) {
                v.resize(size+sizeof(WAH_T)*8-1);
                WAH_T* wah_p = (WAH_T*)(((char*)block_p)+dictionary[key]);
                wah2_extract<WAH_T>(wah_p, v, size);
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    }

    template<typename T>
    inline T* get_pointer_from_dict(enum Dictionary_Keys key) {
        if (dictionary.find(key) != dictionary.end()) {
            if (dictionary[key] != VAL_UNDEFINED) {
                return (T*)(((char*)block_p)+dictionary[key]);
            } else {
                return nullptr;
            }
        } else {
            return nullptr;
        }
    }

    A_T* sparse_extract(A_T* s_p) {
        constexpr A_T MSB_BIT = (A_T)1 << (sizeof(A_T)*8-1);
        A_T num = *s_p;
        s_p++;

        sparse_negated = (num & MSB_BIT);
        num &= ~MSB_BIT; // Remove the bit !

        sparse.clear();
        for (A_T i = 0; i < num; i++) {
            sparse.push_back(*s_p);
            s_p++;
        }

        const size_t CURRENT_N_HAPS = ((haploid_binary_gt_line[internal_binary_gt_line_position]) ? N_SAMPLES : N_HAPS);
        ones = (sparse_negated ? CURRENT_N_HAPS-num : num);

        return s_p;
    }

    A_T* sparse_advance_pointer(A_T* s_p) {
        constexpr A_T MSB_BIT = (A_T)1 << (sizeof(A_T)*8-1);
        A_T num = *s_p;
        s_p++;

        sparse_negated = (num & MSB_BIT);
        num &= ~MSB_BIT; // Remove the bit !

        s_p += num;

        const size_t CURRENT_N_HAPS = ((haploid_binary_gt_line[internal_binary_gt_line_position]) ? N_SAMPLES : N_HAPS);
        ones = (sparse_negated ? CURRENT_N_HAPS-num : num);

        return s_p;
    }

protected:
    const header_t& header;
    const void* block_p;
    const size_t N_SAMPLES;
    const size_t N_HAPS;
    size_t MAX_PLOIDY;
    size_t bcf_lines_in_block;
    size_t binary_gt_lines_in_block;

    std::map<uint32_t, uint32_t> dictionary;

    size_t internal_binary_gt_line_position;

    // WAH
    WAH_T* wah_origin_p;
    WAH_T* wah_p;

    // Sparse
    A_T* sparse_origin_p;
    A_T* sparse_p;
    std::vector<size_t> sparse;
    bool sparse_negated;

    // Weirdness
    bool block_has_weirdness;
    WAH_T* missing_origin_p;
    WAH_T* missing_p;
    WAH_T* eovs_origin_p;
    WAH_T* eovs_p;
    int32_t DEFAULT_PHASING;
    bool block_has_non_uniform_phasing;
    WAH_T* non_uniform_phasing_origin_p;
    WAH_T* non_uniform_phasing_p;

    size_t internal_binary_weirdness_position;
    size_t internal_binary_phase_position;
    std::vector<bool> binary_gt_line_is_wah;
    std::vector<bool> binary_gt_line_is_sorting;
    std::vector<bool> line_has_missing;
    std::vector<bool> line_has_non_uniform_phasing;
    std::vector<bool> line_has_end_of_vector;
    //std::map<size_t, int32_t> non_default_vector_length_positions;
    std::vector<bool> haploid_binary_gt_line;


    std::vector<size_t> allele_counts;
    size_t ones;

    // Internal
    std::vector<bool> y;
    std::vector<A_T> a, b;
    std::vector<bool> y_missing;
    std::vector<bool> y_eovs;
    std::vector<bool> y_phase;
    std::vector<A_T> a_weird, b_weird;
};

template <typename A_T = uint32_t, typename WAH_T = uint16_t>
class AccessorInternalsNewTemplate : public AccessorInternals {
private:
    inline void seek(const size_t& new_position) {
        const size_t OFFSET_MASK = ~((((size_t)-1) >> BM_BLOCK_BITS) << BM_BLOCK_BITS);
        size_t block_id = ((new_position & 0xFFFFFFFF) >> BM_BLOCK_BITS);
        // The offset is relative to the start of the block and is binary gt lines
        uint32_t offset = new_position & OFFSET_MASK;

        // If block ID is not current block
        if (!dp or current_block != block_id) {
            set_gt_block_ptr(block_id);

            // Make DecompressPointer
            dp = make_unique<DecompressPointerGTBlock<A_T, WAH_T> >(header, gt_block_p);
            //std::cerr << "Block ID : " << block_id << " offset : " << offset << std::endl;
        }

        dp->seek(offset);
    }

public:
    size_t fill_genotype_array(int32_t* gt_arr, size_t gt_arr_size, size_t n_alleles, size_t new_position) override {
        seek(new_position);

        return dp->fill_genotype_array_advance(gt_arr, gt_arr_size, n_alleles);
    }

    void fill_allele_counts(size_t n_alleles, size_t new_position) override {
        seek(new_position);

        dp->fill_allele_counts_advance(n_alleles);
    }

    // Directly pass the DecompressPointer Allele counts
    virtual inline const std::vector<size_t>& get_allele_counts() const override {
        return dp->get_allele_count_ref();
    }

    AccessorInternalsNewTemplate(std::string filename) {
        std::fstream s(filename, s.binary | s.in);
        if (!s.is_open()) {
            std::cerr << "Failed to open file " << filename << std::endl;
            throw "Failed to open file";
        }

        // Read the header
        s.read((char *)(&(this->header)), sizeof(header_t));
        s.close();

        // Check magic
        if ((header.first_magic != MAGIC) or (header.last_magic != MAGIC)) {
            std::cerr << "Bad magic" << std::endl;
            std::cerr << "Expected : " << MAGIC << " got : " << header.first_magic << ", " << header.last_magic << std::endl;
            throw "Bad magic";
        }

        // Check version
        if (header.version != 4) {
            std::cerr << "Bad version" << std::endl;
            throw "Bad version";
        }

        file_size = fs::file_size(filename);
        fd = open(filename.c_str(), O_RDONLY, 0);
        if (fd < 0) {
            std::cerr << "Failed to open file " << filename << std::endl;
            throw "Failed to open file";
        }

        // Memory map the file
        file_mmap_p = mmap(NULL, file_size, PROT_READ, MAP_SHARED, fd, 0);
        if (file_mmap_p == NULL) {
            std::cerr << "Failed to memory map file " << filename << std::endl;
            close(fd);
            throw "Failed to mmap file";
        }

        // Test the memory map (first thing is the endianness in the header)
        uint32_t endianness = *(uint32_t*)(file_mmap_p);
        if (endianness != ENDIANNESS) {
            std::cerr << "Bad endianness in memory map" << std::endl;
            throw "Bad endianness";
        }

        if (header.hap_samples == 0) {
            std::cerr << "No samples" << std::endl;
            // Can still be used to "extract" the variant BCF (i.e. loop through the variant BCF and copy it to output... which is useless but ok)
        }

        if (header.ploidy == 0) {
            std::cerr << "Ploidy in header is set to 0 !" << std::endl;
            throw "PLOIDY ERROR";
        }
    }

    virtual ~AccessorInternalsNewTemplate() {
        if (header.zstd and block_p) {
            free(block_p);
            block_p = nullptr;
        }
        munmap(file_mmap_p, file_size);
        close(fd);
    }

protected:
    inline void set_gt_block_ptr(const size_t block_id) {
        set_block_ptr(block_id);
        current_block = block_id;
        char* p = (char*)block_p;

        try {
            p += block_dictionary.at(IBinaryBlock<uint32_t, uint32_t>::KEY_GT_ENTRY);
        } catch (...) {
            std::cerr << "Binary block does not have GT block" << std::endl;
            throw "block error";
        }

        gt_block_p = p;
    }

    inline void set_block_ptr(const size_t block_id) {
        uint32_t* indices_p = (uint32_t*)((uint8_t*)file_mmap_p + header.indices_offset);
        // Find out the block offset
        size_t offset = indices_p[block_id];

        if (header.zstd) {
            size_t compressed_block_size = *(uint32_t*)(((uint8_t*)file_mmap_p) + offset);
            size_t uncompressed_block_size = *(uint32_t*)(((uint8_t*)file_mmap_p) + offset + sizeof(uint32_t));
            void *block_ptr = ((uint8_t*)file_mmap_p) + offset + sizeof(uint32_t)*2;

            if (block_p) {
                free(block_p);
                block_p = nullptr;
            }

            block_p = malloc(uncompressed_block_size);
            if (!block_p) {
                std::cerr << "Failed to allocate memory to decompress block" << std::endl;
                throw "Failed to allocate memory";
            }
            auto result = ZSTD_decompress(block_p, uncompressed_block_size, block_ptr, compressed_block_size);
            if (ZSTD_isError(result)) {
                std::cerr << "Failed to decompress block" << std::endl;
                std::cerr << "Error : " << ZSTD_getErrorName(result) << std::endl;
                throw "Failed to decompress block";
            }
        } else {
            // Set block pointer
            block_p = ((uint8_t*)file_mmap_p) + offset;
        }

        read_dictionary(block_dictionary, (uint32_t*)block_p);
    }

    std::string filename;
    header_t header;
    size_t file_size;
    int fd;
    void* file_mmap_p = nullptr;

    void* block_p = nullptr;
    void* gt_block_p = nullptr;
    std::unique_ptr<DecompressPointerGTBlock<A_T, WAH_T> > dp = nullptr;
    size_t current_block = -1;
    std::map<uint32_t, uint32_t> block_dictionary;
};

#endif /* __ACCESSOR_INTERNALS_NEW_HPP__ */