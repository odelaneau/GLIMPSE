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
#ifndef __ACCESSOR_INTERNALS_HPP__
#define __ACCESSOR_INTERNALS_HPP__

#include <iostream>

#include <string>
#include <unordered_map>
#include "compression.hpp"
#include "xcf.hpp"
#include "block.hpp"
#include "make_unique.hpp"

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include "fs.hpp"
#include "wah.hpp"

#ifndef DEBUGGG
static constexpr bool DEBUG_DECOMP = false;
#else
static constexpr bool DEBUG_DECOMP = true;
#endif
using namespace wah;

template <typename A_T = uint32_t, typename WAH_T = uint16_t>
class DecompressPointer {
public:
    virtual ~DecompressPointer() {}

    virtual void seek(const size_t position) = 0;
    virtual void advance() = 0;
    virtual bool position_is_sparse() const = 0;

    const std::vector<A_T>& get_ref_on_a() const {return a;}
    const std::vector<bool>& get_ref_on_y() const {return y;}
    const std::vector<A_T>& get_sparse_ref() const {return sparse;}
    bool is_negated() const {return sparse_negated;}

    size_t get_current_position() const {return current_position;}
    size_t get_current_position_ones() const {return count_ones;}

protected:
    size_t current_position = 0;
    std::vector<A_T> a;
    std::vector<A_T> b;

    std::vector<bool> y; // Values as arranged by a, can be larger than N_HAPS
    std::vector<A_T> sparse; // Values as sparse
    bool sparse_negated;

    size_t count_ones = 0;
};

template <typename A_T = uint32_t, typename WAH_T = uint16_t>
class DecompressPointerV2 : public DecompressPointer<A_T, WAH_T> {
protected:
    A_T* sparse_extract(A_T* s_p) {
        constexpr A_T MSB_BIT = (A_T)1 << (sizeof(A_T)*8-1);
        A_T num = *s_p;
        s_p++;

        this->sparse_negated = (num & MSB_BIT);
        num &= ~MSB_BIT; // Remove the bit !

        this->sparse.clear();
        if (DEBUG_DECOMP) std::cerr << "DEBUG : Position : " << this->current_position << " Extracted " << num << " Sparse values ";
        for (A_T i = 0; i < num; i++) {
            if (DEBUG_DECOMP) std::cerr << *s_p << " ";
            this->sparse.push_back(*s_p);
            s_p++;
        }
        if (DEBUG_DECOMP) std::cerr << std::endl;

        this->count_ones = (this->sparse_negated ? this->N_HAPS-num : num);

        return s_p;
    }

    A_T* sparse_advance_pointer(A_T* s_p) {
        constexpr A_T MSB_BIT = (A_T)1 << (sizeof(A_T)*8-1);
        A_T num = *s_p;
        s_p++;

        num &= ~MSB_BIT; // Remove the bit !

        s_p += num;
        return s_p;
    }

    virtual void seek_sampled_arrangement(const size_t num = 0) {
        std::iota(this->a.begin(), this->a.end(), 0);
        this->current_position = arrangement_sample_rate * num;

        // Get pointers to the data
        wah_p = wah_origin_p + (indices_p[num]); // Pointer arithmetic handles sizeof(WAH_T)
        sparse_p = sparse_origin_p + (indices_sparse_p[num]);

        if (rearrangement_track[this->current_position]) {
            wah_p = wah2_extract_count_ones(wah_p, this->y, N_HAPS, this->count_ones); // Extract current values
        } else {
            sparse_p = sparse_extract(sparse_p);
        }
    }

public:
    // Decompress Pointer from memory mapped compressed file
    /// @todo pass sparse pointer
    DecompressPointerV2(const size_t N_SITES, const size_t N_HAPS, WAH_T* wah_origin_p, uint32_t* indices_p, A_T* sparse_origin_p, uint32_t* indices_sparse_p, const size_t arrangement_sample_rate, const std::vector<bool>& rearrangement_track):
        N_SITES(N_SITES), N_HAPS(N_HAPS), wah_origin_p(wah_origin_p), indices_p(indices_p), sparse_origin_p(sparse_origin_p), indices_sparse_p(indices_sparse_p), arrangement_sample_rate(arrangement_sample_rate), rearrangement_track(rearrangement_track) {
        this->a.resize(N_HAPS);
        this->b.resize(N_HAPS);
        this->y.resize(N_HAPS + sizeof(WAH_T)*8, 0); // Get some extra space

        wah_p = wah_origin_p;
        sparse_p = sparse_origin_p;
        // Fill with original arrangement
        seek_sampled_arrangement();
    }

    virtual ~DecompressPointerV2() {}

    // Seek out a given position
    void seek(const size_t position) override {
        if (position == this->current_position) { return; }
        size_t advance_steps = 0;
        if ((position > this->current_position) and ((position - this->current_position) < arrangement_sample_rate)) {
            advance_steps = position - this->current_position;
        } else {
            size_t previous_arrangement = position / arrangement_sample_rate;
            advance_steps = position % arrangement_sample_rate;

            //std::cerr << "Seeking position : " << position << std::endl;
            seek_sampled_arrangement(previous_arrangement);
        }

        for (size_t i = 1 /* last step is done below */; i < advance_steps; ++i) {
            private_advance(false);
        }
        advance(); // equiv to private_advance(true)
    }

    // Advance and update inner data structures
    void advance() override {
        private_advance();
    }

    virtual bool position_is_sparse() const override {return !rearrangement_track[this->current_position];}

protected:

    DecompressPointerV2(const size_t N_SITES, const size_t N_HAPS, const size_t arrangement_sample_rate) :
        N_SITES(N_SITES), N_HAPS(N_HAPS), wah_origin_p(nullptr), indices_p(nullptr), sparse_origin_p(nullptr), indices_sparse_p(nullptr), arrangement_sample_rate(arrangement_sample_rate), rearrangement_track(__rt__), __rt__(0) {
        this->a.resize(N_HAPS);
        this->b.resize(N_HAPS);
        this->y.resize(N_HAPS + sizeof(WAH_T)*8, 0); // Get some extra space

        wah_p = nullptr;
        sparse_p = nullptr;
    }

    virtual inline void private_advance(bool extract = true) {
        if (this->current_position >= N_SITES) {
            std::cerr << "Advance called but already at end" << std::endl;
            return;
        }

        A_T u = 0;
        A_T v = 0;

        // Edge case
        if (((this->current_position+1) % arrangement_sample_rate) == 0) {
            std::iota(this->a.begin(), this->a.end(), 0);
        } else {
            if (rearrangement_track[this->current_position]) {
                // PBWT sort
                for (size_t i = 0; i < N_HAPS; ++i) {
                    if (this->y[i] == 0) {
                        this->a[u++] = this->a[i];
                    } else {
                        this->b[v++] = this->a[i];
                    }
                }
                std::copy(this->b.begin(), this->b.begin()+v, this->a.begin()+u);
            }
        }
        if (this->current_position < N_SITES-1) {
            if (rearrangement_track[this->current_position+1]) {
                // Optimisation : Only extract if needed to advance further
                wah_p = wah2_extract_count_ones(wah_p, this->y, N_HAPS, this->count_ones);
                // (in V2 only rearrangement positions are in WAH, everything else is in sparse)
            } else {
                if (extract) {
                    sparse_p = sparse_extract(sparse_p);
                } else {
                    sparse_p = sparse_advance_pointer(sparse_p);
                }
            }
        }
        this->current_position++;
    }

    // Constants, referencing memory mapped file
    const size_t N_SITES;
    const size_t N_HAPS;
    WAH_T* const wah_origin_p;
    uint32_t* const indices_p;
    A_T* const sparse_origin_p;
    uint32_t* const indices_sparse_p;
    const size_t arrangement_sample_rate;

    WAH_T* wah_p; // Gets updated by wah2_extract
    A_T* sparse_p;
private:
    const std::vector<bool>& rearrangement_track;
    std::vector<bool> __rt__;
};

template <typename A_T = uint32_t, typename WAH_T = uint16_t>
class DecompressPointerV3 : public DecompressPointerV2<A_T, WAH_T> {
public:
    DecompressPointerV3(const size_t N_SITES, const size_t N_HAPS, const size_t arrangement_sample_rate, uint32_t *indices, void *file_mmap_p, bool compressed) :
    DecompressPointerV2<A_T, WAH_T>(N_SITES, N_HAPS, arrangement_sample_rate), indices(indices), file_mmap_p(file_mmap_p), compressed(compressed) {
        seek_sampled_arrangement();
    }

    virtual ~DecompressPointerV3() {
        if (compressed and block) {
            free(block);
        }
    }

    bool position_is_sparse() const override {return !select_track[this->current_position % this->arrangement_sample_rate]; }

protected:
    void seek_sampled_arrangement(const size_t num = 0) override {
        std::iota(this->a.begin(), this->a.end(), 0);
        this->current_position = this->arrangement_sample_rate * num;

        // Find out the block offset
        size_t offset = indices[num];

        if (compressed) {
            size_t compressed_block_size = *(uint32_t*)(((uint8_t*)file_mmap_p) + offset);
            size_t uncompressed_block_size = *(uint32_t*)(((uint8_t*)file_mmap_p) + offset + sizeof(uint32_t));
            void *block_ptr = ((uint8_t*)file_mmap_p) + offset + sizeof(uint32_t)*2;
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
            block = ((uint8_t*)file_mmap_p) + offset;
        }

        // Fill dictionnary
        Block::fill_dictionnary(block, dictionnary);

        // Get pointers to the data
        this->wah_p = (WAH_T*)((uint8_t*)block + dictionnary[Block::Dictionnary_Keys::KEY_WAH]);
        this->sparse_p = (A_T*)((uint8_t*)block + dictionnary[Block::Dictionnary_Keys::KEY_SPARSE]);

        if (dictionnary.find(Block::Dictionnary_Keys::KEY_SMALL_BLOCK) != dictionnary.end()) {
            block_length = dictionnary[Block::Dictionnary_Keys::KEY_SMALL_BLOCK];
        } else {
            block_length = this->arrangement_sample_rate;
        }

        select_track.resize(block_length + sizeof(WAH_T)*8);  // With some extra space because of WAH alignment
        WAH_T* st_p = (WAH_T*)((uint8_t*)block + dictionnary[Block::Dictionnary_Keys::KEY_SELECT]);
        wah2_extract<WAH_T>(st_p, select_track, block_length);
        select_track.resize(block_length);

        if (select_track[0]) {
            this->wah_p = wah2_extract_count_ones(this->wah_p, this->y, this->N_HAPS, this->count_ones); // Extract current values
        } else {
            this->sparse_p = this->sparse_extract(this->sparse_p);
        }
    }

    inline void private_advance(bool extract = true) override {
        if (this->current_position >= this->N_SITES) {
            std::cerr << "Advance called but already at end" << std::endl;
            return;
        }

        A_T u = 0;
        A_T v = 0;

        size_t next_position = this->current_position+1;

        // Edge case
        if ((next_position % this->arrangement_sample_rate) == 0) {
            // Because block may need to be decompressed
            seek_sampled_arrangement(next_position / this->arrangement_sample_rate);
        } else {
            if (select_track[this->current_position % this->arrangement_sample_rate]) {
                // PBWT sort
                for (size_t i = 0; i < this->N_HAPS; ++i) {
                    if (this->y[i] == 0) {
                        this->a[u++] = this->a[i];
                    } else {
                        this->b[v++] = this->a[i];
                    }
                }
                std::copy(this->b.begin(), this->b.begin()+v, this->a.begin()+u);
            }

            if (this->current_position < this->N_SITES-1) {
                if (select_track[next_position % this->arrangement_sample_rate]) {
                    // Optimisation : Only extract if needed to advance further
                    this->wah_p = wah2_extract_count_ones(this->wah_p, this->y, this->N_HAPS, this->count_ones);
                    // (in V2 only rearrangement positions are in WAH, everything else is in sparse)
                } else {
                    if (extract) {
                        this->sparse_p = this->sparse_extract(this->sparse_p);
                    } else {
                        this->sparse_p = this->sparse_advance_pointer(this->sparse_p);
                    }
                }
            }
        }
        this->current_position = next_position;
    }

    std::unordered_map<uint32_t, uint32_t> dictionnary;
    // Pointer to current block
    uint32_t *indices;
    void *block = NULL;
    size_t block_length;
    std::vector<bool> select_track;
    void *file_mmap_p;
    const bool compressed;
};

class AccessorInternals {
public:
    virtual ~AccessorInternals() {}
    virtual size_t fill_genotype_array(int32_t* gt_arr, size_t gt_arr_size, size_t n_alleles, size_t position) = 0;
    // Fill genotype array also fills allele counts, so this is only to be used when fill_genotype_array is not called (e.g., to recompute AC only)
    virtual void fill_allele_counts(size_t n_alleles, size_t position) = 0;
    virtual inline const std::vector<size_t>& get_allele_counts() const {return allele_counts;}
    //virtual const std::unordered_map<size_t, std::vector<size_t> >& get_missing_sparse_map() const = 0;
    //virtual const std::unordered_map<size_t, std::vector<size_t> >& get_phase_sparse_map() const = 0;
protected:
    std::vector<size_t> allele_counts;

    const size_t BM_BLOCK_BITS = 15;
};

template <typename A_T = uint32_t, typename WAH_T = uint16_t>
class AccessorInternalsTemplate : public AccessorInternals {
public:
    size_t fill_genotype_array(int32_t* gt_arr, size_t gt_arr_size, size_t n_alleles, size_t new_position) override {
        // Conversion from new to old, for the moment
        const size_t OFFSET_MASK = ~((((size_t)-1) >> BM_BLOCK_BITS) << BM_BLOCK_BITS);
        size_t block_id = ((new_position & 0xFFFFFFFF) >> BM_BLOCK_BITS);
        uint32_t offset = new_position & OFFSET_MASK;
        size_t position = block_id * header.ss_rate + offset;
        // The conversion will be unnecessary in the new version

        this->allele_counts.resize(n_alleles);
        size_t total_alt = 0;

        dp->seek(position);

        // Set REF / first ALT
        if (dp->position_is_sparse()) { /* SPARSE */
            int32_t default_gt = dp->is_negated() ? 1 : 0;
            int32_t sparse_gt = dp->is_negated() ? 0 : 1;
            for (size_t i = 0; i < N_HAPS; ++i) {
                gt_arr[i] = bcf_gt_unphased(default_gt) | ((i & 1) & DEFAULT_PHASED);
            }
            for (const auto& i : dp->get_sparse_ref()) {
                //if constexpr (DEBUG_DECOMP) std::cerr << "Setting variant at " << i << std::endl;
                gt_arr[i] = bcf_gt_unphased(sparse_gt) | ((i & 1) & DEFAULT_PHASED);
            }
        } else { /* SORTED WAH */
            auto& a = dp->get_ref_on_a();
            auto& y = dp->get_ref_on_y();
            for (size_t i = 0; i < N_HAPS; ++i) {
                gt_arr[a[i]] = bcf_gt_unphased(y[i]) | ((a[i] & 1) & DEFAULT_PHASED); /// @todo Phase
            }
        }
        const auto ones = dp->get_current_position_ones();
        allele_counts[1] = ones;
        total_alt = ones;
        dp->advance();

        // If other ALTs (ALTs are 1 indexed, because 0 is REF)
        for (size_t alt_allele = 2; alt_allele < n_alleles; ++alt_allele) {
            if (dp->position_is_sparse()) { /* SPARSE */
                if (dp->is_negated()) { // There can only be one negated because must be more than all others combined
                    // All non set positions are now filled
                    for (size_t i = 0; i < N_HAPS; ++i) {
                        // Only overwrite refs
                        if (bcf_gt_allele(gt_arr[i]) == 0) {
                            gt_arr[i] = bcf_gt_unphased(alt_allele) | ((i & 1) & DEFAULT_PHASED);
                        }
                    }
                    for (const auto& i : dp->get_sparse_ref()) {
                        // Restore overwritten refs
                        if (bcf_gt_allele(gt_arr[i]) == (int)alt_allele) {
                            gt_arr[i] = bcf_gt_unphased(0) | ((i & 1) & DEFAULT_PHASED);
                        }
                    }
                } else {
                    // Fill normally
                    for (const auto& i : dp->get_sparse_ref()) {
                        gt_arr[i] = bcf_gt_unphased(alt_allele) | ((i & 1) & DEFAULT_PHASED);
                    }
                }
            } else { /* SORTED WAH */
                auto& a = dp->get_ref_on_a();
                auto& y = dp->get_ref_on_y();
                for (size_t i = 0; i < N_HAPS; ++i) {
                    if (y[i]) {
                        gt_arr[a[i]] = bcf_gt_unphased(alt_allele) | ((a[i] & 1) & DEFAULT_PHASED); /// @todo Phase
                    }
                }
            }
            const auto ones = dp->get_current_position_ones();
            allele_counts[alt_allele] = ones;
            total_alt += ones;
            dp->advance();
        }

        // Set missing info
        size_t total_missing = 0;
        if (missing_map.find(position) != missing_map.end()) {
            total_missing = missing_map.at(position).size();
            for (auto pos : missing_map.at(position)) {
                gt_arr[pos] = bcf_gt_missing | ((pos & 1) & DEFAULT_PHASED);
            }
        }

        // Set non default phase info
        if (non_default_phase_map.find(position) != non_default_phase_map.end()) {
            for (auto pos : non_default_phase_map.at(position)) {
                gt_arr[pos] ^= 0x1; // Toggle phase bit
            }
        }

        // Set ref allele count (all haps that don't have an alt allele or missing)
        allele_counts[0] = this->N_HAPS - total_alt - total_missing;

        return gt_arr_size; // Old version
    }

    void fill_allele_counts(size_t n_alleles, size_t new_position) override {
        // Conversion from new to old, for the moment
        size_t block_id = new_position >> BM_BLOCK_BITS;
        uint32_t offset = ((uint32_t)new_position << (32-BM_BLOCK_BITS)) >> (32-BM_BLOCK_BITS); // Remove block id bits
        size_t position = block_id * header.ss_rate + offset;
        // The conversion will be unnecessary in the new version

        this->allele_counts.resize(n_alleles);
        size_t total_alt = 0;

        // First allele count
        dp->seek(position);
        const auto ones = dp->get_current_position_ones();
        allele_counts[1] = ones;
        total_alt = ones;
        dp->advance();

        // Subsequent allele counts
        for (size_t alt_allele = 2; alt_allele < n_alleles; ++alt_allele) {
            const auto ones = dp->get_current_position_ones();
            allele_counts[alt_allele] = ones;
            total_alt += ones;
            dp->advance();
        }

        // Get missing info
        size_t total_missing = 0;
        if (missing_map.find(position) != missing_map.end()) {
            total_missing = missing_map.at(position).size();
        }

        // Set ref allele count (all haps that don't have an alt allele or missing)
        allele_counts[0] = this->N_HAPS - total_alt - total_missing;
    }

    AccessorInternalsTemplate(std::string filename) {
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
        if (header.version != 2 and header.version != 3) {
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
        file_mmap = mmap(NULL, file_size, PROT_READ, MAP_SHARED, fd, 0);
        if (file_mmap == NULL) {
            std::cerr << "Failed to memory map file " << filename << std::endl;
            close(fd);
            throw "Failed to mmap file";
        }

        // Test the memory map (first thing is the endianness in the header)
        uint32_t endianness = *(uint32_t*)(file_mmap);
        if (endianness != ENDIANNESS) {
            std::cerr << "Bad endianness in memory map" << std::endl;
            throw "Bad endianness";
        }

        if (header.hap_samples == 0) {
            std::cerr << "No samples" << std::endl;
            // Can still be used to "extract" the variant BCF (i.e. loop through the variant BCF and copy it to output... which is useless but ok)
        }

        // Rearrangement track decompression
        if (header.version == 2) {
            rearrangement_track.resize(header.num_variants + header.wah_bytes*8);
            WAH_T* rt_p = (WAH_T*)((uint8_t*)file_mmap + header.rearrangement_track_offset);
            wah2_extract<WAH_T>(rt_p, rearrangement_track, header.num_variants);
        }

        bool use_ppas = !header.iota_ppa;
        if (use_ppas) {
            std::cerr << "Versions 2-3 do not use ppas" << std::endl;
            throw "Bad file format";
        }

        if (header.has_missing) {
            /// @todo
            // Extract the missing data as a hash table
            // Key is BM index (cannot be entry number because of the -r option)
            // val is vector of sparse entries of missing data
            //std::cerr << "Has missing !" << std::endl;
            uint32_t limit = header.phase_info_offset; // Next offset in file is phase info
            fill_sparse_map(header.missing_offset, limit, missing_map);
            //std::cerr << "missing map size : " << missing_map.size() << std::endl;
        }

        if (header.non_uniform_phasing) {
            /// @todo
            // Extract the non uniform phasing data as a hash table
            // Key is BM index (cannot be entry number because of the -r option)
            // Val is vector of sparse entries of non default phasing
            uint32_t limit = file_size; // This is last in file for the moment
            //std::cerr << "Has non uniform phasing" << std::endl;
            fill_sparse_map(header.phase_info_offset, limit, non_default_phase_map);
        }

        // This may seem silly but is to make sure the cast is ok
        // Because in the header the value is a bool inside a union with a uint8_t
        DEFAULT_PHASED = header.default_phased ? 1 : 0;

        PLOIDY = header.ploidy;
        if (PLOIDY == 0) {
            std::cerr << "Ploidy in header is set to 0 !" << std::endl;
            throw "PLOIDY ERROR";
        }

        N_HAPS = header.hap_samples;

        dp = generate_decompress_pointer();
    }

    virtual ~AccessorInternalsTemplate() {
        munmap(file_mmap, file_size);
        close(fd);
    }

    //inline const std::unordered_map<size_t, std::vector<size_t> >& get_missing_sparse_map() const override {
    //    return this->missing_map;
    //}
    //inline const std::unordered_map<size_t, std::vector<size_t> >& get_phase_sparse_map() const override {
    //    return this->non_default_phase_map;
    //}

protected:
    void fill_sparse_map(uint32_t offset, uint32_t end, std::unordered_map<size_t, std::vector<size_t> >& map) {
        map.clear();
        uint8_t* ptr = ((uint8_t*)file_mmap + offset);
        uint8_t* end_ptr = ((uint8_t*)file_mmap + end);
        if (DEBUG_DECOMP) printf("DEBUG : ptr is : %p end ptr is %p\n", ptr, end_ptr);

        while(ptr < end_ptr) {
            // Decode index
            uint32_t bm_index = *((uint32_t*)ptr);
            ptr += sizeof(uint32_t);
            // Decode quantity
            A_T qty = *((A_T*)ptr);
            ptr += sizeof(A_T);
            // Fill sparse vector
            std::vector<size_t> sparse_pos_vector;
            for (A_T i = 0; i < qty; ++i) {
                A_T pos = *((A_T*)ptr);
                ptr += sizeof(A_T);
                sparse_pos_vector.push_back(pos);
            }
            if (DEBUG_DECOMP) std::cerr << "DEBUG : sparse entry at BM " << bm_index << ", " << qty << " : ";
            if (DEBUG_DECOMP) for (auto s : sparse_pos_vector) {std::cerr << s << " ";}
            if (DEBUG_DECOMP) std::cerr << std::endl;
            // Add to map
            map.insert({bm_index, sparse_pos_vector});
        }
    }

    inline std::unique_ptr<DecompressPointer<A_T, WAH_T> > generate_decompress_pointer(size_t offset = 0) {
        const size_t N_SITES = header.num_variants;
        const size_t N_HAPS = header.hap_samples;
        const size_t arrangements_sample_rate = header.ss_rate;

        std::unique_ptr<DecompressPointer<A_T, WAH_T> > dp = nullptr;

        if (header.version == 3) {
            bool compressed = header.zstd;
            dp = make_unique<DecompressPointerV3<A_T, WAH_T> >(N_SITES, N_HAPS, arrangements_sample_rate, (uint32_t*)((uint8_t*)file_mmap + header.indices_offset), file_mmap, compressed);
        } else {
            WAH_T* wah_origin_p = (WAH_T*)((uint8_t*)file_mmap + header.wahs_offset);
            A_T* sparse_origin_p = (A_T*)((uint8_t*)file_mmap + header.sparse_offset);
            uint32_t* indices_p = (uint32_t*)((uint8_t*)file_mmap + header.indices_offset);
            uint32_t* indices_sparse_p = (uint32_t*)((uint8_t*)file_mmap + header.indices_sparse_offset);

            dp = make_unique<DecompressPointerV2<A_T, WAH_T> >(N_SITES, N_HAPS, wah_origin_p, indices_p, sparse_origin_p, indices_sparse_p, arrangements_sample_rate, rearrangement_track);
        }

        dp->seek(offset);

        return dp;
    }

    std::string filename;
    header_t header;
    size_t file_size;
    int fd;
    void* file_mmap = NULL;

    size_t PLOIDY = 2;
    size_t N_HAPS;

    std::vector<bool> rearrangement_track;

    std::unordered_map<size_t, std::vector<size_t> > missing_map;
    std::unordered_map<size_t, std::vector<size_t> > non_default_phase_map;
    int32_t DEFAULT_PHASED;

    std::unique_ptr<DecompressPointer<A_T, WAH_T> > dp;
};

#endif /* __ACCESSOR_INTERNALS_HPP__ */