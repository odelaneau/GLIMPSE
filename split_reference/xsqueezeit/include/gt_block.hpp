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

#ifndef __GT_BLOCK_HPP__
#define __GT_BLOCK_HPP__

#include "interfaces.hpp"
#include "internal_gt_record.hpp"

class GTBlockDict {
public:
    enum Dictionary_Keys : uint32_t {
        // Element (Scalar) keys
        KEY_DICTIONNARY_SIZE = (uint32_t)-1,
        KEY_BCF_LINES = 0,
        KEY_BINARY_LINES = 1,
        KEY_MAX_LINE_PLOIDY = 2,
        KEY_DEFAULT_PHASING = 3,
        // Line (Vector) keys
        KEY_LINE_SORT = 0x10,
        KEY_LINE_SELECT = 0x11,
        KEY_LINE_HAPLOID = 0x12,
        KEY_LINE_VECTOR_LENGTH = 0x15,
        KEY_LINE_MISSING = 0x16,
        KEY_LINE_NON_UNIFORM_PHASING = 0x17,
        KEY_LINE_END_OF_VECTORS = 0x18,
        // Matrix keys
        KEY_MATRIX_WAH = 0x20,
        KEY_MATRIX_SPARSE = 0x21,
        KEY_MATRIX_MISSING = 0x26,
        KEY_MATRIX_NON_UNIFORM_PHASING = 0x27,
        KEY_MATRIX_END_OF_VECTORS = 0x28,
    };

    enum Dictionary_Vals : uint32_t {
        VAL_UNDEFINED = (uint32_t)-1,
    };
};

class PBWTSorter {
public:
    struct MissingPred {
        static inline bool check(int32_t gt_arr_entry, int32_t _) {
            (void)_; // Unused
            return bcf_gt_is_missing(gt_arr_entry);
        }
    };
    struct RawPred {
        static inline bool check(const int32_t gt_arr_entry, const int32_t raw) {
            return gt_arr_entry == raw;
        }
    };
    struct NonDefaultPhasingPred {
        static inline bool check_widx(const size_t& idx, const int32_t gt_arr_entry, const int32_t default_phasing) {
            // This is because there is no phasing information on the first allele
            /// @todo This only works for PLOIDY 1 and 2
            // For ploidy > 2 requires modulo PLOIDY which is slow AF
            return (idx & 1) and (bcf_gt_is_phased(gt_arr_entry) != default_phasing);
        }
    };
    struct WeirdnessPred {
        static inline bool check(const int32_t gt_arr_entry, const int32_t _ /*, const int32_t default_phasing*/) {
            return bcf_gt_is_missing(gt_arr_entry) or (gt_arr_entry == bcf_int32_missing) or gt_arr_entry == bcf_int32_vector_end; // or bcf_gt_is_phased(gt_arr_entry) != default_phasing;
        }
    };
    template<typename T, class Pred, const size_t V_LEN_RATIO = 1>
    inline void pred_pbwt_sort(std::vector<T>& a, std::vector<T>& b, int32_t* gt_arr, const size_t N, const int32_t _) {
        size_t u = 0;
        size_t v = 0;

        for (size_t j = 0; j < N; ++j) {
            if (!Pred::check(gt_arr[a[j]/V_LEN_RATIO], _)) { // If non pred
                a[u] = a[j];
                u++;
            } else { // if pred
                b[v] = a[j];
                v++;
            }
        }
        std::copy(b.begin(), b.begin()+v, a.begin()+u);
    }

    /// @todo V_LEN_RATIO doesn't work on y
    template<typename T>
    inline void bool_pbwt_sort(std::vector<T>& a, std::vector<T>& b, const std::vector<bool>& y, const size_t N) {
        size_t u = 0;
        size_t v = 0;
        for (size_t i = 0; i < N; ++i) {
            if (y[i] == 0) {
                a[u++] = a[i];
            } else {
                b[v++] = a[i];
            }
        }
        std::copy(b.begin(), b.begin()+v, a.begin()+u);
    }

    template<typename T>
    inline void bool_pbwt_sort_two(std::vector<T>& a, std::vector<T>& b, const std::vector<bool>& y1, const std::vector<bool>& y2, const size_t N) {
        size_t u = 0;
        size_t v = 0;
        for (size_t i = 0; i < N; ++i) {
            bool y = y1[i] or y2[i];
            if (y == 0) {
                a[u++] = a[i];
            } else {
                b[v++] = a[i];
            }
        }
        std::copy(b.begin(), b.begin()+v, a.begin()+u);
    }
};

template<typename A_T = uint32_t, typename WAH_T = uint16_t>
class GtBlock : public IWritableBCFLineEncoder, public BCFBlock, public GTBlockDict, protected PBWTSorter {
public:
    const size_t PLOIDY_2 = 2;

    GtBlock(const size_t NUM_SAMPLES, const size_t BLOCK_BCF_LINES, const size_t MAC_THRESHOLD, const int32_t default_phasing = 0) :
        BCFBlock(BLOCK_BCF_LINES),
        MAC_THRESHOLD(MAC_THRESHOLD),
        default_ploidy(PLOIDY_2),
        default_phasing(default_phasing),
        effective_binary_gt_lines_in_block(0),
        line_has_missing(BLOCK_BCF_LINES, false),
        line_has_non_uniform_phasing(BLOCK_BCF_LINES, false),
        line_has_end_of_vector(BLOCK_BCF_LINES, false),
        default_vector_length(PLOIDY_2), max_vector_length(1),
        line_allele_counts(BLOCK_BCF_LINES),
        line_alt_alleles_number(BLOCK_BCF_LINES, 0),
        a(NUM_SAMPLES*PLOIDY_2), b(NUM_SAMPLES*PLOIDY_2),
        a_weirdness(NUM_SAMPLES*PLOIDY_2), b_weirdness(NUM_SAMPLES*PLOIDY_2)
    {
        //std::cerr << "Block MAC Thr : " << MAC_THRESHOLD << std::endl;
        // Reset a
        std::iota(a.begin(), a.end(), 0);
        std::iota(a_weirdness.begin(), a_weirdness.end(), 0);
    }

    inline uint32_t get_id() const override { return IBinaryBlock<uint32_t, uint32_t>::KEY_GT_ENTRY; }

    void write_to_stream(std::fstream& ofs) override {
        //if (effective_bcf_lines_in_block != BLOCK_BCF_LINES) {
        //    std::cerr << "Block with fewer BCF lines written to stream" << std::endl;
        //} else {
        //    std::cerr << "Block written to stream" << std::endl;
        //}
        //std::cerr << "BCF lines : " << effective_bcf_lines_in_block << " binary lines : " << effective_binary_gt_lines_in_block << std::endl;

        size_t block_start_pos = ofs.tellp();
        //size_t block_end_pos(0);
        size_t dictionary_pos(0);

        fill_dictionary();

        dictionary_pos = write_dictionary(ofs, dictionary);

        write_writables(ofs, block_start_pos); // Updates dictionary

        update_dictionary(ofs, dictionary_pos, dictionary);
    }

private:
    inline void scan_genotypes(const bcf_file_reader_info_t& bcf_fri) {
        const auto LINE_MAX_PLOIDY = bcf_fri.ngt / bcf_fri.n_samples;
        if (LINE_MAX_PLOIDY > max_vector_length) {
            max_vector_length = LINE_MAX_PLOIDY;
        }
        //if (LINE_MAX_PLOIDY != default_vector_length) {
            //non_default_vector_length_positions[effective_bcf_lines_in_block] = LINE_MAX_PLOIDY;
        //}

        /// @todo (performance?) we could use reindexing the BCF indices instead of push_back
        if (LINE_MAX_PLOIDY == 1) {
            haploid_line_found = true;
            haploid_binary_gt_line.push_back(true);
        } else {
            haploid_binary_gt_line.push_back(false);
        }

        auto& allele_counts = line_allele_counts[effective_bcf_lines_in_block];
        allele_counts.resize(bcf_fri.line->n_allele, 0);

        line_alt_alleles_number[effective_bcf_lines_in_block] = bcf_fri.line->n_allele-1;

        for (size_t i = 0; i < bcf_fri.n_samples; ++i) {
            // Go over all alleles
            for (size_t j = 0; j < LINE_MAX_PLOIDY; ++j) {
                const size_t index = i*LINE_MAX_PLOIDY+j;
                auto bcf_allele = bcf_fri.gt_arr[index];
                if (j) {
                    // Phasing only applies on the second haplotypes and following for polyploid
                    // This is a quirk of the BCF format, the phase bit of the first genotype is not used...
                    // 0/1 => 0x02 04 and 0|1 => 0x02 05, see VCF / BCF specifications
                    // https://samtools.github.io/hts-specs/
                    // Will be set to non zero if phase changes
                    if (bcf_gt_is_phased(bcf_allele) != default_phasing) {
                        //std::cerr << "default : " << default_phasing << " found : " << bcf_gt_is_phased(bcf_allele) << std::endl;
                        non_uniform_phasing = true;
                        line_has_non_uniform_phasing[effective_bcf_lines_in_block] = true;
                    }
                }
                if (bcf_gt_is_missing(bcf_allele) or (bcf_allele == bcf_int32_missing)) {
                    missing_found = true;
                    line_has_missing[effective_bcf_lines_in_block] = true;
                } else if (bcf_allele == bcf_int32_vector_end) {
                    end_of_vector_found = true;
                    line_has_end_of_vector[effective_bcf_lines_in_block] = true;

                    /// @todo set info for non uniform vector lengths (in other function)
                } else {
                    try {
                        allele_counts.at(bcf_gt_allele(bcf_allele))++;
                    } catch (std::exception& e) {
                        printf("Bad access at : 0x%08x\n", bcf_gt_allele(bcf_allele));
                        printf("Value in array is : 0x%08x\n", bcf_allele);
                        throw "Unknown allele error !";
                    }
                }
            }
        }
    }

public:
    inline void encode_line(const bcf_file_reader_info_t& bcf_fri) override {
        scan_genotypes(bcf_fri);

        auto& allele_counts = line_allele_counts[effective_bcf_lines_in_block];
        const auto LINE_MAX_PLOIDY = bcf_fri.ngt / bcf_fri.n_samples;
        //std::cerr << "[DEBUG] : Line " << effective_bcf_lines_in_block
        //          << " ngt : " << bcf_fri.ngt << std::endl;
        //for (size_t i = 0; i < bcf_fri.ngt; ++i) {
        //    std::cerr << bcf_fri.gt_arr[i] << " ";
        //}
        //std::cerr << std::endl;

        // For all alt alleles (1 if bi-allelic variant site)
        for (size_t alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {

            //std::cerr << "[DEBUG] a : ";
            //for (auto& e : a) std::cerr << e << " ";
            //std::cerr << std::endl;

            const size_t minor_allele_count = std::min(allele_counts[alt_allele], bcf_fri.ngt - allele_counts[alt_allele]);
            if (minor_allele_count > MAC_THRESHOLD) {
                uint32_t _; // Unused
                bool __; // Unused

                binary_gt_line_is_wah.push_back(true);
                if (LINE_MAX_PLOIDY == 1) {
                    // The order is given by a1 instead of a
                    auto a1 = haploid_rearrangement_from_diploid(a);
                    wah_encoded_binary_gt_lines.push_back(wah::wah_encode2_with_size<WAH_T>(bcf_fri.gt_arr, alt_allele, a1, bcf_fri.ngt, _, __));
                    // Here pbwt_sort1 does the logic for the sorting no need to use a1
                    pbwt_sort1(a, b, bcf_fri.gt_arr, bcf_fri.ngt, alt_allele);
                } else if (LINE_MAX_PLOIDY == 2) {
                    wah_encoded_binary_gt_lines.push_back(wah::wah_encode2_with_size<WAH_T>(bcf_fri.gt_arr, alt_allele, a, bcf_fri.ngt, _, __));
                    pbwt_sort(a, b, bcf_fri.gt_arr, bcf_fri.ngt, alt_allele);
                } else {
                    std::cerr << "Cannot handle ploidy of " << LINE_MAX_PLOIDY << " with default ploidy " << default_ploidy << std::endl;
                    throw "PLOIDY ERROR";
                }
            } else {
                int32_t sparse_allele = 0; // If 0 means sparse is negated
                if (allele_counts[alt_allele] == minor_allele_count) {
                    sparse_allele = alt_allele;
                }
                // Sparse does not depend on ploidy because we don't use the rearrangement
                // We directly encode the correct number of alleles
                sparse_encoded_binary_gt_lines.emplace_back(SparseGtLine<A_T>(effective_binary_gt_lines_in_block, bcf_fri.gt_arr, bcf_fri.ngt, sparse_allele));
                binary_gt_line_is_wah.push_back(false);
            }
            effective_binary_gt_lines_in_block++;
        }

        bool weird_line = false;
        uint32_t _(0); // Unused
        bool __(false); // Unused
        if (line_has_missing[effective_bcf_lines_in_block]) {
            weird_line = true;
            if (LINE_MAX_PLOIDY == 1) {
                auto a1 = haploid_rearrangement_from_diploid(a_weirdness);
                wah_encoded_missing_lines.push_back(wah::wah_encode2_with_size<WAH_T, A_T, MissingPred>(bcf_fri.gt_arr, _, a1, bcf_fri.ngt, _, __));
            } else {
                wah_encoded_missing_lines.push_back(wah::wah_encode2_with_size<WAH_T, A_T, MissingPred>(bcf_fri.gt_arr, _, a_weirdness, bcf_fri.ngt, _, __));
            }
        }
        if (line_has_end_of_vector[effective_bcf_lines_in_block]) {
            weird_line = true;
            if (LINE_MAX_PLOIDY == 1) {
                auto a1 = haploid_rearrangement_from_diploid(a_weirdness);
                wah_encoded_end_of_vector_lines.push_back(wah::wah_encode2_with_size<WAH_T, A_T, RawPred>(bcf_fri.gt_arr, bcf_int32_vector_end, a1, bcf_fri.ngt, _, __));
            } else {
                wah_encoded_end_of_vector_lines.push_back(wah::wah_encode2_with_size<WAH_T, A_T, RawPred>(bcf_fri.gt_arr, bcf_int32_vector_end, a_weirdness, bcf_fri.ngt, _, __));
            }
        }

        // Phasing info is not PBWT reordered with weirdness
        /// @todo double compress this structure would save space, e.g., if all the entries are the same
        if (line_has_non_uniform_phasing[effective_bcf_lines_in_block]) {
            wah_encoded_non_uniform_phasing_lines.push_back(wah::wah_encode2_with_size<WAH_T, NonDefaultPhasingPred>(bcf_fri.gt_arr, default_phasing, bcf_fri.ngt, _));
        }

        //std::cerr << "[DEBUG] Weirdness a : ";
        //for (auto& e : a_weirdness) std::cerr << e << " ";
        //std::cerr << std::endl;

        if (weird_line) {
            // PBWT sort on weirdness
            if (LINE_MAX_PLOIDY != default_ploidy) {
                if (LINE_MAX_PLOIDY == 1 and default_ploidy == 2) {
                    // Sort a, b as if they were homozygous with ploidy 2
                    //pred_pbwt_sort<A_T, WeirdnessPred, 2>(a_weirdness, b_weirdness, bcf_fri.gt_arr, a_weirdness.size(), 0 /*unused*/);
                } else {
                    /// @todo add and handle polyploid PBWT sort
                    std::cerr << "Cannot handle ploidy of " << LINE_MAX_PLOIDY << " with default ploidy " << default_ploidy << std::endl;
                    throw "PLOIDY ERROR";
                }
            } else {
                pred_pbwt_sort<A_T, WeirdnessPred>(a_weirdness, b_weirdness, bcf_fri.gt_arr, a_weirdness.size(), 0 /*unused*/);
            }
        }

        effective_bcf_lines_in_block++;
    }

    virtual ~GtBlock() {}

protected:
    const size_t MAC_THRESHOLD;
    size_t default_ploidy;
    int32_t default_phasing;

    size_t effective_binary_gt_lines_in_block;

    // For handling missing
    bool missing_found = false;
    std::vector<bool> line_has_missing;

    // For handling non default phase
    bool non_uniform_phasing = false;
    std::vector<bool> line_has_non_uniform_phasing;

    // For handling end of vector
    bool end_of_vector_found = false;
    std::vector<bool> line_has_end_of_vector;

    // For handling mixed ploidy
    uint32_t default_vector_length;
    uint32_t max_vector_length;
    //std::map<size_t, int32_t> non_default_vector_length_positions;
    bool haploid_line_found = false;
    std::vector<bool> haploid_binary_gt_line;
    //std::set<int32_t> vector_lengths;

    // Set when the binary gt line is WAH encoded (else sparse)
    std::vector<bool> binary_gt_line_is_wah;
    // Set when the binary gt line PBWT sorts the samples
    //std::vector<bool> binary_gt_line_sorts;

    // 2D Structures
    std::vector<std::vector<WAH_T> > wah_encoded_binary_gt_lines;
    std::vector<SparseGtLine<A_T> > sparse_encoded_binary_gt_lines;

    std::vector<std::vector<WAH_T> > wah_encoded_missing_lines;
    std::vector<std::vector<WAH_T> > wah_encoded_non_uniform_phasing_lines;
    std::vector<std::vector<WAH_T> > wah_encoded_end_of_vector_lines;

    // Internal
    std::vector<std::vector<size_t> > line_allele_counts;
    std::vector<size_t> line_alt_alleles_number;

    std::unordered_map<uint32_t, uint32_t> dictionary;

private:
    inline void fill_dictionary() {
        dictionary[KEY_BCF_LINES] = effective_bcf_lines_in_block;
        dictionary[KEY_BINARY_LINES] = effective_binary_gt_lines_in_block;
        dictionary[KEY_MAX_LINE_PLOIDY] = max_vector_length;
        dictionary[KEY_DEFAULT_PHASING] = default_phasing;

        // Those are offsets
        dictionary[KEY_LINE_SORT] = VAL_UNDEFINED;
        dictionary[KEY_LINE_SELECT] = VAL_UNDEFINED;
        dictionary[KEY_MATRIX_WAH] = VAL_UNDEFINED;
        dictionary[KEY_MATRIX_SPARSE] = VAL_UNDEFINED;

        if (missing_found) {
            //std::cerr << "[DEBUG] Missing found" << std::endl;
            // Is an offset
            dictionary[KEY_LINE_MISSING] = VAL_UNDEFINED;
            dictionary[KEY_MATRIX_MISSING] = VAL_UNDEFINED;
        }

        if (end_of_vector_found) {
            //std::cerr << "[DEBUG] EOV found" << std::endl;
            // Is an offset
            dictionary[KEY_LINE_END_OF_VECTORS] = VAL_UNDEFINED;
            dictionary[KEY_MATRIX_END_OF_VECTORS] = VAL_UNDEFINED;
        }

        if (non_uniform_phasing) {
            //std::cerr << "[DEBUG] Non uniform phasing found" << std::endl;
            // Is an offset
            dictionary[KEY_LINE_NON_UNIFORM_PHASING] = VAL_UNDEFINED;
            dictionary[KEY_MATRIX_NON_UNIFORM_PHASING] = VAL_UNDEFINED;
        }

        //if (non_default_vector_length_positions.size()) {
        //    std::cerr << "[DEBUG] Non default vector length found" << std::endl;
        //    // Is an offset
        //    dictionary[KEY_LINE_VECTOR_LENGTH] = VAL_UNDEFINED;
        //}

        if (haploid_line_found) {
            //std::cerr << "[DEBUG] Haploid line found" << std::endl;
            dictionary[KEY_LINE_HAPLOID] = VAL_UNDEFINED;
        }
    }

    inline void write_writables(std::fstream& s, const size_t& block_start_pos) {
        /// @note .at(key) is used to make sure the key is in the dictionary !
        /// @todo handle the exceptions, however there should be none if this class is implemented correctly

        size_t written_bytes = size_t(s.tellp());
        size_t total_bytes = size_t(s.tellp());

        // Write Sort
        dictionary.at(KEY_LINE_SORT) = (uint32_t)((size_t)s.tellp()-block_start_pos);
        write_boolean_vector_as_wah(s, binary_gt_line_is_wah);

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        //std::cout << "sort " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        // Write Select
        dictionary.at(KEY_LINE_SELECT) = dictionary[KEY_LINE_SORT]; // Same is used

        // Write WAH
        dictionary.at(KEY_MATRIX_WAH) = (uint32_t)((size_t)s.tellp()-block_start_pos);
        for (const auto& wah : wah_encoded_binary_gt_lines) {
            //for (auto& wah_word : wah) print_wah2(wah_word);
            write_vector(s, wah);
        }

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        //std::cout << "WAH " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        // Write Sparse
        dictionary.at(KEY_MATRIX_SPARSE) = (uint32_t)((size_t)s.tellp()-block_start_pos);
        for (const auto& sparse_line : sparse_encoded_binary_gt_lines) {
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
            write_vector(s, sparse);
        }
        //std::cout << "Written " << sparse_encoded_binary_gt_lines.size() << " sparse lines" << std::endl;

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        //std::cout << "sparse " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        // Optional write missing
        if (missing_found) {
            //for (auto v : line_has_missing) { std::cerr << (v ? "1" : "0"); }
            //std::cerr << std::endl;
            auto v = reindex_binary_vector_from_bcf_to_binary_lines(line_has_missing);
            //for (auto _ : v) { std::cerr << (_ ? "1" : "0"); }
            //std::cerr << std::endl;
            dictionary.at(KEY_LINE_MISSING) = (uint32_t)((size_t)s.tellp()-block_start_pos);
            write_boolean_vector_as_wah(s, v);
            dictionary.at(KEY_MATRIX_MISSING) = (uint32_t)((size_t)s.tellp()-block_start_pos);
            write_vector_of_vectors(s, wah_encoded_missing_lines);

            written_bytes = size_t(s.tellp()) - total_bytes;
            total_bytes += written_bytes;
            //std::cout << "missing " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;
        }

        // Optional write end of vectors
        if (end_of_vector_found) {
            auto v = reindex_binary_vector_from_bcf_to_binary_lines(line_has_end_of_vector);
            dictionary.at(KEY_LINE_END_OF_VECTORS) = (uint32_t)((size_t)s.tellp()-block_start_pos);
            write_boolean_vector_as_wah(s, v);
            dictionary.at(KEY_MATRIX_END_OF_VECTORS) = (uint32_t)((size_t)s.tellp()-block_start_pos);
            write_vector_of_vectors(s, wah_encoded_end_of_vector_lines);

            written_bytes = size_t(s.tellp()) - total_bytes;
            total_bytes += written_bytes;
            //std::cout << "end of vectors " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;
        }

        // Optional write non uniform phasing
        if (non_uniform_phasing) {
            auto v = reindex_binary_vector_from_bcf_to_binary_lines(line_has_non_uniform_phasing);
            dictionary.at(KEY_LINE_NON_UNIFORM_PHASING) = (uint32_t)((size_t)s.tellp()-block_start_pos);
            write_boolean_vector_as_wah(s, v);
            dictionary.at(KEY_MATRIX_NON_UNIFORM_PHASING) = (uint32_t)((size_t)s.tellp()-block_start_pos);
            write_vector_of_vectors(s, wah_encoded_non_uniform_phasing_lines);
            // for (auto& v : wah_encoded_non_uniform_phasing_lines) {
            //     for (auto& w : v) {
            //         print_wah2(w);
            //     }
            // }

            written_bytes = size_t(s.tellp()) - total_bytes;
            total_bytes += written_bytes;
            //std::cout << "non uniform phasing " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;
        }

        // Optional write non default vector lengths
        //if (non_default_vector_length_positions.size()) {
        //    dictionary.at(KEY_LINE_VECTOR_LENGTH) = (uint32_t)((size_t)s.tellp()-block_start_pos);
        //
        //    //TODO
        //    throw "Not implemented yet !";
        //}

        if (haploid_line_found) {
            dictionary.at(KEY_LINE_HAPLOID) = (uint32_t)((size_t)s.tellp()-block_start_pos);
            write_boolean_vector_as_wah(s, haploid_binary_gt_line);
        }

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        //std::cout << "others " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;
    }

    // Binary lines >= BCF lines
    std::vector<bool> reindex_binary_vector_from_bcf_to_binary_lines(const std::vector<bool>& v) {
        std::vector<bool> result(effective_binary_gt_lines_in_block);

        size_t binary_offset = 0;
        for (size_t i = 0; i < effective_bcf_lines_in_block; ++i) {
            result[binary_offset++] = v[i];
            for (size_t _ = 1; _ < line_alt_alleles_number[i]; ++_) {
                result[binary_offset++] = 0; // Fill
            }
        }

        if (binary_offset != effective_binary_gt_lines_in_block) {
            std::cerr << "Block internal data is corrupted" << std::endl;
        }

        return result;
    }

    template<typename T>
    inline void write_vector(std::fstream& s, const std::vector<T>& v) {
        static_assert(!std::is_same<T, bool>::value, "bool is implementation defined therefore is not portable");
        s.write(reinterpret_cast<const char*>(v.data()), v.size() * sizeof(decltype(v.back())));
    }

    template<typename T>
    inline void write_vector_of_vectors(std::fstream& s, const std::vector<std::vector<T> >& cv) {
        for (const auto& v : cv) {
            write_vector(s, v);
        }
    }

    template<typename _WAH_T = WAH_T>
    inline void write_boolean_vector_as_wah(std::fstream& s, std::vector<bool>& v) {
        auto wah = wah::wah_encode2<_WAH_T>(v);
        write_vector(s, wah);
    }

    // Internal PBWT ordering vectors
    std::vector<A_T> a;
    std::vector<A_T> b;

    std::vector<A_T> a_weirdness;
    std::vector<A_T> b_weirdness;
};

#endif /* __GT_BLOCK_HPP__ */