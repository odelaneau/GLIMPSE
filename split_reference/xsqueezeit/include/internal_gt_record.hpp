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

#ifndef __INTERNAL_GT_RECORD_HPP__
#define __INTERNAL_GT_RECORD_HPP__

#include "xcf.hpp"
#include "block.hpp" // For SparseGTLine ...

/// @todo remove unused variable ngt
template<typename T, const size_t V_LEN_RATIO = 1>
inline void pbwt_sort_(std::vector<T>& a, std::vector<T>& b, int32_t* gt_arr, const size_t ngt, int32_t alt_allele) {
    size_t u = 0;
    size_t v = 0;

    for (size_t i = 0; i < a.size(); ++i) {
        const auto haplotype_id = a[i];
        if (bcf_gt_allele(gt_arr[haplotype_id/V_LEN_RATIO]) != alt_allele) { // If non alt allele
            a[u] = a[i];
            u++;
        } else { // if alt allele
            b[v] = a[i];
            v++;
        }
    }
    std::copy(b.begin(), b.begin()+v, a.begin()+u);
}

template<typename T>
inline void pbwt_sort(std::vector<T>& a, std::vector<T>& b, int32_t* gt_arr, const size_t ngt, int32_t alt_allele) {
    pbwt_sort_<T, 1>(a, b, gt_arr, ngt, alt_allele);
}

// Sort a haplotype pairs given single haplotype
template<typename T>
inline void pbwt_sort1(std::vector<T>&a, std::vector<T>&b, int32_t* gt_arr, const size_t ngt, int32_t alt_allele) {
    pbwt_sort_<T, 2>(a, b, gt_arr, ngt, alt_allele);
}

template<typename T = uint32_t>
class InternalGtRecord {
private:

    // Scans the genotypes for missing data and phasing as well as does the allele counts
    inline void scan_genotypes(const bcf_file_reader_info_t& bcf_fri) {
        for (size_t i = 0; i < bcf_fri.n_samples; ++i) {
            for (size_t j = 0; j < PLOIDY; ++j) {
                const size_t index = i*PLOIDY+j;
                auto bcf_allele = bcf_fri.gt_arr[index];
                if (j) {
                    // Phasing only applies on the second haplotypes and following for polyploid
                    // This is a quirk of the BCF format, the phase bit of the first genotype is not used...
                    // 0/1 => 0x02 04 and 0|1 => 0x02 05, see VCF / BCF specifications
                    // https://samtools.github.io/hts-specs/
                    // Will be set to non zero if phase changes
                    if (bcf_gt_is_phased(bcf_allele) != default_is_phased) {
                        sparse_non_default_phasing.push_back(index);
                    }
                }
                /// @todo check if this works with END_OF_VECTOR for male samples on chrX
                /// the bcf_get_genotypes() call should handle this and return "missing"
                /// see vcf.h in htslib 0x80000000 is missing and 0x80000001 is END_OF_VECTOR in bcf
                /// in htslib 0 is missing, see bcf_gt_is_missing
                if (bcf_gt_is_missing(bcf_allele) or (bcf_allele == bcf_int32_missing)) {
                    sparse_missing.push_back(index);
                } else if (bcf_allele == bcf_int32_vector_end) {
                    std::cerr << "End of vector index : " << index << std::endl;
                    std::cerr << "Mixed ploidy detected" << std::endl;
                    has_mixed_ploidy = true;
                    end_of_vector[index] = true;
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
    InternalGtRecord(const bcf_file_reader_info_t& bcf_fri, std::vector<T>& a, std::vector<T>& b, int32_t default_is_phased, const size_t MAC_THRESHOLD, size_t& variant_counter, const size_t RESET_SORT_BLOCK_LENGTH) :
    PLOIDY(bcf_fri.ngt/bcf_fri.n_samples), n_alleles(bcf_fri.line->n_allele), allele_counts(bcf_fri.line->n_allele, 0), rearrangements(bcf_fri.line->n_allele-1, false), default_is_phased(default_is_phased),
    end_of_vector(PLOIDY, false) {
        scan_genotypes(bcf_fri);

        // For all alt alleles (1 if bi-allelic variant site)
        for (size_t alt_allele = 1; alt_allele < n_alleles; ++alt_allele) {
            if ((variant_counter % RESET_SORT_BLOCK_LENGTH) == 0) {
                // Restart from natural order
                std::iota(a.begin(), a.end(), 0);
            }

            const size_t minor_allele_count = std::min(allele_counts[alt_allele], bcf_fri.ngt - allele_counts[alt_allele]);
            if (minor_allele_count > MAC_THRESHOLD) {
                uint32_t _; // Unused
                bool __; // Unused
                wahs.push_back(wah::wah_encode2_with_size(bcf_fri.gt_arr, alt_allele, a, bcf_fri.ngt, _, __));
                const size_t SORT_THRESHOLD = MAC_THRESHOLD; // For next version
                if (minor_allele_count > SORT_THRESHOLD) {
                    rearrangements[alt_allele-1] = true;
                    pbwt_sort(a, b, bcf_fri.gt_arr, bcf_fri.ngt, alt_allele);
                }
            } else {
                int32_t sparse_allele = 0; // If 0 means sparse is negated
                if (allele_counts[alt_allele] == minor_allele_count) {
                    sparse_allele = alt_allele;
                }
                sparse_lines.emplace_back(SparseGtLine<T>(variant_counter, bcf_fri.gt_arr, bcf_fri.ngt, sparse_allele));
            }

            variant_counter++;
        }
    }

    /// @todo wah/sparse and rearrangements may not be the same in the future

    const size_t PLOIDY = 0;
    const size_t n_alleles = 0;
    std::vector<size_t> allele_counts;
    std::vector<bool> rearrangements;
    int32_t default_is_phased = 0;
    std::vector<std::vector<uint16_t> > wahs;
    std::vector<SparseGtLine<T> > sparse_lines;
    std::vector<T> sparse_missing;
    std::vector<T> sparse_non_default_phasing;
    bool has_mixed_ploidy = false;
    std::vector<bool> end_of_vector;
};

#endif /* __INTERNAL_GT_RECORD_HPP__ */