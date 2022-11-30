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

#include "xcf.hpp"
#include "fs.hpp"

bool has_extension(const std::string& filename, const std::string& extension) {
    const std::regex ext_regex(std::string(".+\\") + extension);
    return std::regex_match(filename.c_str(), ext_regex);
}

/**
 * @brief Creates csi index for given file
 *
 * @param filename file to index
 * @param n_threads optional parameters, number of threads for indexing, default 1
 * */
void create_index_file(std::string filename, int n_threads) {
    int ret = bcf_index_build3(filename.c_str() /* input */,
                               NULL /* Output filename, or NULL to add .csi/.tbi */,
                               14 /* Positive to generate CSI, or 0 to generate TBI, CSI bin size (CSI default is 14) */,
                               n_threads /* n_threads */);

    if (ret != 0) {
        if (ret == -2) {
            std::cerr << "index: failed to open " << filename << std::endl;
            throw "Failed to open file";
        } else if (ret == -3) {
            std::cerr << "index: " << filename << " is in a format that cannot be usefully indexed" << std::endl;
            throw "Failed to index";
        } else {
            std::cerr << "index: failed to create index for " << filename << std::endl;
            throw "Failed to index";
        }
    }
}

static void initialize_bcf_file_reader_common(bcf_file_reader_info_t& bcf_fri, const std::string& filename) {
    while(!bcf_sr_add_reader(bcf_fri.sr, filename.c_str())) {
        if (bcf_fri.sr->errnum == idx_load_failed) {
            bcf_sr_destroy(bcf_fri.sr);
            bcf_fri.sr = bcf_sr_init();
            bcf_fri.sr->collapse = COLLAPSE_NONE;
            bcf_fri.sr->require_index = 1;
            std::cerr << "Index is missing, indexing " << filename << std::endl;
            try {
                create_index_file(filename);
            } catch (const char* e) {
                std::cerr << "Failed to index" << std::endl << e << std::endl;
                throw "bcf_synced_reader read error";
            }
        } else {
            std::cerr << "Failed to read file " << filename << std::endl;
            std::cerr << "Reason : " << bcf_sr_strerror(bcf_fri.sr->errnum) << std::endl;
            throw "bcf_synced_reader read error";
        }
    }
    bcf_fri.n_samples = bcf_hdr_nsamples(bcf_fri.sr->readers[0].header);

    bcf_fri.var_count = 0; // Already set by default
    bcf_fri.gt_arr = nullptr; // Already set by default
    bcf_fri.size_gt_arr = 0; // Already set by default
    bcf_fri.ngt = 0; // Already set by default
    bcf_fri.line = nullptr; // Already set by default
    bcf_fri.line_num = 0;
    bcf_fri.line_alt_alleles_extracted = 0; // Already set by default
}

void initialize_bcf_file_reader(bcf_file_reader_info_t& bcf_fri, const std::string& filename) {
    bcf_fri.sr = bcf_sr_init();
    bcf_fri.sr->collapse = COLLAPSE_NONE;
    bcf_fri.sr->require_index = 0;

    initialize_bcf_file_reader_common(bcf_fri, filename);
}

void destroy_bcf_file_reader(bcf_file_reader_info_t& bcf_fri) {
    if (bcf_fri.gt_arr != nullptr) {
        free(bcf_fri.gt_arr); // C allocation from realloc() inside htslib
        bcf_fri.gt_arr = nullptr;
    }
    if (bcf_fri.sr != nullptr) {
        bcf_sr_destroy(bcf_fri.sr);
        bcf_fri.sr = nullptr;
    }
    bcf_fri.var_count = 0;
    bcf_fri.size_gt_arr = 0;
    bcf_fri.ngt = 0;
    bcf_fri.line = nullptr; /// @todo check if should be freed, probably not
    bcf_fri.line_num = 0;
    bcf_fri.line_alt_alleles_extracted = 0;
}

void initialize_bcf_file_reader_with_region(bcf_file_reader_info_t& bcf_fri, const std::string& filename, const std::string& region, bool is_file) {
    bcf_fri.sr = bcf_sr_init();
    bcf_fri.sr->collapse = COLLAPSE_NONE;
    bcf_fri.sr->require_index = 1;
    if (bcf_sr_set_regions(bcf_fri.sr, region.c_str(), is_file) < 0) {
        std::cerr << "Failed to read the regions : " << region << std::endl;
        destroy_bcf_file_reader(bcf_fri);
        throw "Failed to read the regions";
    }

    initialize_bcf_file_reader_common(bcf_fri, filename);
}

unsigned int bcf_next_line(bcf_file_reader_info_t& bcf_fri) {
    unsigned int nset = bcf_sr_next_line(bcf_fri.sr);
    if (nset) {
        bcf_fri.line_num++;
        bcf_fri.line = bcf_sr_get_line(bcf_fri.sr, 0 /* First file of possibly more in sr */);
        bcf_fri.line_alt_alleles_extracted = 0;
    }
    return nset;
}

/// @todo check if better to return a std::vector or simply reuse one passed by reference
inline void extract_next_variant_and_update_bcf_sr(std::vector<bool>& samples, bcf_file_reader_info_t& bcf_fri) {
    // No checks are done on bcf_fri in this function because it is supposed to be called in large loops
    // Check the sanity of bcf_fri before calling this function

    samples.resize(bcf_fri.n_samples * 2 /* Diploid */);
    bool do_extract = true;

    // If there is no current line or line has been fully extracted
    if ((!bcf_fri.line) or (bcf_fri.line_alt_alleles_extracted+1 == bcf_fri.line->n_allele)) {
        // Go to next line
        if (bcf_next_line(bcf_fri) == 0) { // End of file
            do_extract = false;
        }
    }

    if (do_extract) {
        const int alt_allele = bcf_fri.line_alt_alleles_extracted+1;

        if (bcf_fri.line_alt_alleles_extracted == 0) {
            bcf_unpack(bcf_fri.line, BCF_UN_STR);
            // Here info about the variant could be extracted
            bcf_fri.ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.size_gt_arr));
            int line_max_ploidy = bcf_fri.ngt / bcf_fri.n_samples;

            if (line_max_ploidy != 2) {
                std::cerr << "[ERROR] Ploidy of samples is different than 2" << std::endl;
                throw "BadPloidy";
            }
        } // Else is already unpacked (when multiple ALT alleles)

        for (size_t i = 0; i < bcf_fri.n_samples; ++i) {
            int32_t* ptr = (bcf_fri.gt_arr) + i * 2; /* line_max_ploidy */
            for (size_t j = 0; j < 2 /* ploidy */; ++j) {
                bool a = (bcf_gt_allele(ptr[j]) == alt_allele);
                samples[i*2+j] = a;
            }
        }

        bcf_fri.line_alt_alleles_extracted++;
        bcf_fri.var_count++;
    } else {
        // No next line, indicate this by returning empty samples.
        samples.clear();
    }
}

/**
 * @brief extract_samples Returns a vector with the sample IDs
 * @param bcf_fri BCF File Reader Info (must already by initialized)
 * */
std::vector<std::string> extract_samples(const bcf_file_reader_info_t& bcf_fri) {
    std::vector<std::string> res;
    for (size_t i = 0; i < bcf_fri.n_samples; ++i) {
        res.push_back(std::string(bcf_fri.sr->readers[0].header->samples[i]));
    }
    return res;
}

/**
 * @brief writes a vector of strings to a given file
 *
 * The strings are newline separated in the output file
 *
 * @param v the vector of strings to be written
 * @param ofname the file to write the vector to
 * */
void string_vector_to_file(const std::vector<std::string>& v, const std::string& ofname) {
    std::fstream s(ofname, s.out);
    if (!s.is_open()) {
        std::cerr << "Failed to open file " << ofname << std::endl;
        throw "Failed to open file";
    }

    for (const auto& str : v) {
        s << str << std::endl;
    }
    s.close();
}

/**
 * @brief reads a file and outputs a vector of strings that corresponds to the lines in the file
 * @param ifname the file to be read
 * @return the vector with the lines of the file
 * */
std::vector<std::string> string_vector_from_file(const std::string& ifname) {
    std::fstream s(ifname, s.in);
    if (!s.is_open()) {
        std::cerr << "Failed to open file " << ifname << std::endl;
        throw "Failed to open file";
    }

    std::vector<std::string> result;

    for (std::string line; std::getline(s, line); ) {
        result.push_back(line);
    }

    s.close();

    return result;
}

/**
 * @brief extract_samples Returns a vector with the sample IDs of file
 * @param fname File name of VCF / BCF file for sample ID extraction
 * */
std::vector<std::string> extract_samples(const std::string& fname) {
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, fname);
    auto res = extract_samples(bcf_fri);
    destroy_bcf_file_reader(bcf_fri);
    return res;
}

/**
 * @brief remove_samples Removes all the samples from a VCF / BCF file and saves the result in a VCF / BCF file
 * @param ifname Input file name
 * @param ofname Output file name
 * @return The number of lines (entries)
 * */
size_t remove_samples(const std::string& ifname, const std::string& ofname) {
    // Input file
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, ifname);

    // Output file
    htsFile *fp = hts_open(ofname.c_str(), "wb"); /// @todo wb wz or other

    /* The bottleneck of VCF reading is parsing of genotype fields. If the reader knows in advance that only subset of samples is needed (possibly no samples at all), the performance of bcf_read() can be significantly improved by calling bcf_hdr_set_samples after bcf_hdr_read(). */
    int ret = bcf_hdr_set_samples(bcf_fri.sr->readers[0].header, NULL, 0 /* 0 is file 1 is list */); // All samples is "-", NULL is none
    if (ret < 0) {
        std::cerr << "Failed to set no samples in header for file " << ifname << std::endl;
        throw "Failed to remove samples";
    }
    bcf_hdr_t *hdr = bcf_hdr_dup(bcf_fri.sr->readers[0].header);

    // Write the header
    ret = bcf_hdr_write(fp, hdr);
    if (ret < 0) {
        std::cerr << "Failed to write header to file " << ofname << std::endl;
        throw "Failed to remove samples";
    }

    // Write the variants
    size_t lines = 0;
    while (bcf_next_line(bcf_fri)) {
        /// @note this could maybe be improved by removing the dup / destroy
        bcf1_t *rec = bcf_dup(bcf_fri.line);
        ret = bcf_write1(fp, hdr, rec);
        bcf_destroy(rec);
        lines++;
    }

    // Close everything
    hts_close(fp);
    bcf_hdr_destroy(hdr);
    destroy_bcf_file_reader(bcf_fri);

    return lines;
}

/**
 * @brief counts the number of entries (lines) in VCF / BCF file
 * @param ifname Input file name
 * @return number of entries
 * */
size_t count_entries(const std::string& ifname) {
    // Input file
    bcf_file_reader_info_t bcf_fri;

    initialize_bcf_file_reader(bcf_fri, ifname);
    bcf_fri.sr->max_unpack = BCF_UN_STR; // Minimal, but does not really impact perf

    /* The bottleneck of VCF reading is parsing of genotype fields. If the reader knows in advance that only subset of samples is needed (possibly no samples at all), the performance of bcf_read() can be significantly improved by calling bcf_hdr_set_samples after bcf_hdr_read(). */
    int ret = bcf_hdr_set_samples(bcf_fri.sr->readers[0].header, NULL, 0 /* 0 is file 1 is list */); // All samples is "-", NULL is none
    if (ret < 0) {
        std::cerr << "Failed to set no samples in header for file " << ifname << std::endl;
        throw "Failed to count entries";
    }

    size_t lines = 0;
    while (bcf_next_line(bcf_fri)) {
        lines++;
    }

    destroy_bcf_file_reader(bcf_fri);

    return lines;
}

/**
 * @brief extract a matrix of bits for the genotype data, outer index is variant,
 *        inner index is samples (two bits per diploid individual)
 *
 * @param filename The input VCF/BCF file
 * */
std::vector<std::vector<bool> > extract_matrix(std::string filename) {
    std::vector<std::vector<bool> > matrix(1);

    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, filename);

    std::vector<std::vector<bool> > line;

    extract_next_variant_and_update_bcf_sr(matrix.back(), bcf_fri);
    for(; !matrix.back().empty(); matrix.push_back(std::vector<bool>(matrix.back().size())), extract_next_variant_and_update_bcf_sr(matrix.back(), bcf_fri));
    matrix.pop_back();

    destroy_bcf_file_reader(bcf_fri);
    return matrix;
}

/**
 * @brief utility function to check if two VCF/BCF files have the same genotype dat
 *        matrix (only GT data checked, no variant, extra info etc.)
 * */
bool matrices_differ(std::string f1, std::string f2) {
    auto m1 = extract_matrix(f1);
    auto m2 = extract_matrix(f2);

    return matrices_differ(m1, m2);
}

std::string unique_id(bcf1_t* bcf_entry_p) {
    std::stringstream ss;
    ss << bcf_entry_p->rid << "_"; // CHROM
    ss << bcf_entry_p->pos << "_"; // POS
    for (size_t i = 0; i < bcf_entry_p->n_allele; ++i) { // REF ALTs
        ss << bcf_entry_p->d.allele[i] << "_";
    }
    return ss.str();
}

void unphase_xcf(const std::string& ifname, const std::string& ofname) {
    bcf_file_reader_info_t bcf_fri;
    htsFile* fp = NULL;
    bcf_hdr_t* hdr = NULL;
    initialize_bcf_file_reader(bcf_fri, ifname);

    fp = hts_open(ofname.c_str(), ofname.compare("-") ? "wb" : "wu"); // "-" for stdout
    if (fp == NULL) {
        std::cerr << "Could not open " << ofname << std::endl;
        throw "File open error";
    }

    // Duplicate the header from the input bcf
    hdr = bcf_hdr_dup(bcf_fri.sr->readers[0].header);

    // Write the header to the new file
    int ret = bcf_hdr_write(fp, hdr);
    if (ret < 0) {
        std::cerr << "Failed to write header !" << std::endl;
        throw "File write error";
    }

    while(bcf_next_line(bcf_fri)) {
        const int32_t PLOIDY = 2;
        bcf1_t *rec = bcf_fri.line;

        // Unpack the line and get genotypes
        bcf_unpack(bcf_fri.line, BCF_UN_STR);
        int ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.size_gt_arr));
        int line_max_ploidy = ngt / bcf_fri.n_samples;

        // Check ploidy, only support diploid for the moment
        if (line_max_ploidy != PLOIDY) {
            std::cerr << "[ERROR] Ploidy of samples is different than 2" << std::endl;
            exit(-1); // Change this
        }

        for (size_t i = 0; i < bcf_fri.n_samples; ++i) {
            int32_t allele_0 = bcf_gt_allele(bcf_fri.gt_arr[i*2]);
            int32_t allele_1 = bcf_gt_allele(bcf_fri.gt_arr[i*2+1]);
            bcf_fri.gt_arr[i*2] = bcf_gt_unphased(std::min(allele_0, allele_1));
            bcf_fri.gt_arr[i*2+1] = bcf_gt_unphased(std::max(allele_0, allele_1));
        }
        bcf_update_genotypes(hdr, rec, bcf_fri.gt_arr, bcf_hdr_nsamples(hdr) * PLOIDY);

        ret = bcf_write1(fp, hdr, rec);
        if (ret < 0) {
            std::cerr << "Failed to write record !" << std::endl;
            throw "File write error";
        }
    }

    // Close / Release ressources
    hts_close(fp);
    bcf_hdr_destroy(hdr);
    destroy_bcf_file_reader(bcf_fri);
}

#include <random>
void unphase_xcf_random(const std::string& ifname, const std::string& ofname) {
    bcf_file_reader_info_t bcf_fri;
    htsFile* fp = NULL;
    bcf_hdr_t* hdr = NULL;
    initialize_bcf_file_reader(bcf_fri, ifname);

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(1, 2);

    fp = hts_open(ofname.c_str(), ofname.compare("-") ? "wb" : "wu"); // "-" for stdout
    if (fp == NULL) {
        std::cerr << "Could not open " << ofname << std::endl;
        throw "File open error";
    }

    // Duplicate the header from the input bcf
    hdr = bcf_hdr_dup(bcf_fri.sr->readers[0].header);

    // Write the header to the new file
    int ret = bcf_hdr_write(fp, hdr);
    if (ret < 0) {
        std::cerr << "Failed to write header !" << std::endl;
        throw "File write error";
    }

    while(bcf_next_line(bcf_fri)) {
        const int32_t PLOIDY = 2;
        bcf1_t *rec = bcf_fri.line;

        // Unpack the line and get genotypes
        bcf_unpack(bcf_fri.line, BCF_UN_STR);
        bcf_fri.ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.size_gt_arr));
        int line_max_ploidy = bcf_fri.ngt / bcf_fri.n_samples;

        // Check ploidy, only support diploid for the moment
        if (line_max_ploidy != PLOIDY) {
            std::cerr << "[ERROR] Ploidy of samples is different than 2" << std::endl;
            exit(-1); // Change this
        }

        for (size_t i = 0; i < bcf_fri.n_samples; ++i) {
            int32_t allele_0 = bcf_gt_allele(bcf_fri.gt_arr[i*2]);
            int32_t allele_1 = bcf_gt_allele(bcf_fri.gt_arr[i*2+1]);
            if (distrib(gen) == 1) {
                bcf_fri.gt_arr[i*2] = bcf_gt_unphased(allele_0);
                bcf_fri.gt_arr[i*2+1] = bcf_gt_unphased(allele_1);
            } else {
                bcf_fri.gt_arr[i*2] = bcf_gt_unphased(allele_1);
                bcf_fri.gt_arr[i*2+1] = bcf_gt_unphased(allele_0);
            }
        }
        bcf_update_genotypes(hdr, rec, bcf_fri.gt_arr, bcf_hdr_nsamples(hdr) * PLOIDY);

        ret = bcf_write1(fp, hdr, rec);
        if (ret < 0) {
            std::cerr << "Failed to write record !" << std::endl;
            throw "File write error";
        }
    }

    // Close / Release ressources
    hts_close(fp);
    bcf_hdr_destroy(hdr);
    destroy_bcf_file_reader(bcf_fri);
}

void sprinkle_missing_xcf(const std::string& ifname, const std::string& ofname) {
    bcf_file_reader_info_t bcf_fri;
    htsFile* fp = NULL;
    bcf_hdr_t* hdr = NULL;
    initialize_bcf_file_reader(bcf_fri, ifname);

    fp = hts_open(ofname.c_str(), ofname.compare("-") ? "wb" : "wu"); // "-" for stdout
    if (fp == NULL) {
        std::cerr << "Could not open " << ofname << std::endl;
        throw "File open error";
    }

    // Duplicate the header from the input bcf
    hdr = bcf_hdr_dup(bcf_fri.sr->readers[0].header);

    // Write the header to the new file
    int ret = bcf_hdr_write(fp, hdr);
    if (ret < 0) {
        std::cerr << "Failed to write header !" << std::endl;
        throw "File write error";
    }

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(1, 100);

    while(bcf_next_line(bcf_fri)) {
        const int32_t PLOIDY = 2;
        bcf1_t *rec = bcf_fri.line;

        // Unpack the line and get genotypes
        bcf_unpack(bcf_fri.line, BCF_UN_STR);
        bcf_fri.ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.size_gt_arr));
        int line_max_ploidy = bcf_fri.ngt / bcf_fri.n_samples;

        // Check ploidy, only support diploid for the moment
        if (line_max_ploidy != PLOIDY) {
            std::cerr << "[ERROR] Ploidy of samples is different than 2" << std::endl;
            exit(-1); // Change this
        }

        for (size_t i = 0; i < bcf_fri.n_samples; ++i) {
            if (distrib(gen) == 1) { // 1% chance to happen
                bcf_fri.gt_arr[i*2] = bcf_gt_is_phased(bcf_fri.gt_arr[i*2]) ?
                                      bcf_gt_phased(-1) :
                                      bcf_gt_unphased(-1);
            }
            if (distrib(gen) == 1) { // 1% chance to happen
                bcf_fri.gt_arr[i*2+1] = bcf_gt_is_phased(bcf_fri.gt_arr[i*2+1]) ?
                                        bcf_gt_phased(-1) :
                                        bcf_gt_unphased(-1);
            }
            // Note : "#define bcf_gt_missing 0" is the unphased version, therefore we use the -1 idx instead
        }
        bcf_update_genotypes(hdr, rec, bcf_fri.gt_arr, bcf_hdr_nsamples(hdr) * PLOIDY);

        ret = bcf_write1(fp, hdr, rec);
        if (ret < 0) {
            std::cerr << "Failed to write record !" << std::endl;
            throw "File write error";
        }
    }

    // Close / Release ressources
    hts_close(fp);
    bcf_hdr_destroy(hdr);
    destroy_bcf_file_reader(bcf_fri);
}

std::vector<std::vector<bool> > extract_common_to_matrix(const std::string& ifname, const double MAF) {
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, ifname);

    const int32_t PLOIDY = 2;
    const size_t N_HAPS = bcf_fri.n_samples*PLOIDY;

    std::vector<std::vector<bool> > output;

    size_t extraction_counter = 0;
    while(bcf_next_line(bcf_fri)) {
        // Unpack the line and get genotypes
        bcf_unpack(bcf_fri.line, BCF_UN_STR);
        bcf_fri.ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.size_gt_arr));
        int line_max_ploidy = bcf_fri.ngt / bcf_fri.n_samples;

        // Check ploidy, only support diploid for the moment
        if (line_max_ploidy != PLOIDY) {
            std::cerr << "[ERROR] Ploidy of samples is different than 2" << std::endl;
            exit(-1); // Change this
        }

        for (int alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {
            uint32_t alt_allele_counter = 0;
            for (size_t i = 0; i < N_HAPS; ++i) {
                if (bcf_gt_allele(bcf_fri.gt_arr[i]) == alt_allele) {
                    alt_allele_counter++;
                }
            }
            uint32_t minor_allele_count = std::min((uint32_t)bcf_fri.n_samples - alt_allele_counter, alt_allele_counter);

            if (minor_allele_count > N_HAPS*MAF) {
                extraction_counter++;
                output.push_back(std::vector<bool>(N_HAPS, false));

                for (size_t i = 0; i < N_HAPS; ++i) {
                    if (bcf_gt_allele(bcf_fri.gt_arr[i]) == alt_allele) {
                        output.back()[i] = true;
                    } else {
                        // default is false
                    }
                }
            }
        }
    }

    // Close / Release ressources
    destroy_bcf_file_reader(bcf_fri);

    return output;
}

/**
 * @brief This function replaces the samples by a single sample that represent the
 *        position of the now missing samples data in the binary matrix. This is allows
 *        to make a link between the vcf entry and the binary matrix where the sample
 *        data is actually stored. This allows to tabix/csi index the variant file and
 *        execute region queries on it. The resulting records then point to a position
 *        in the matrix, which can be used to decompress the missing data.
 *
 * */
size_t replace_samples_by_pos_in_binary_matrix(const std::string& ifname, const std::string& ofname, std::string xsi_fname, const bool new_version, const size_t BLOCK_LENGTH, const size_t BM_BLOCK_BITS) {
    // Input file
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, ifname);

    // Output file
    htsFile *fp = hts_open(ofname.c_str(), "wz"); /// @todo wb wz or other

    /* The bottleneck of VCF reading is parsing of genotype fields. If the reader knows in advance that only subset of samples is needed (possibly no samples at all), the performance of bcf_read() can be significantly improved by calling bcf_hdr_set_samples after bcf_hdr_read(). */
    int ret = bcf_hdr_set_samples(bcf_fri.sr->readers[0].header, NULL, 0 /* 0 is file 1 is list */); // All samples is "-", NULL is none
    if (ret < 0) {
        std::cerr << "Failed to set pseudo sample in header for file " << ifname << std::endl;
        throw "Failed to remove samples";
    }
    bcf_hdr_t *hdr = bcf_hdr_dup(bcf_fri.sr->readers[0].header);
    bcf_hdr_add_sample(hdr, "BIN_MATRIX_POS");
    bcf_hdr_append(hdr, "##FORMAT=<ID=BM,Number=1,Type=Integer,Description=\"Position in GT Binary Matrix\">");

    if (xsi_fname.compare("")) {
        bcf_hdr_append(hdr, std::string("##XSI=").append(std::string(basename((char*)xsi_fname.c_str()))).c_str());
    }

    if (bcf_hdr_sync(hdr) < 0) {
        std::cerr << "bcf_hdr_sync() failed ... oh well" << std::endl;
    }

    // Write the header
    ret = bcf_hdr_write(fp, hdr);
    if (ret < 0) {
        std::cerr << "Failed to write header to file " << ofname << std::endl;
        throw "Failed to remove samples";
    }

    // Write the variants
    size_t pos = 0;
    size_t line = 0;
    int32_t offset = 0;
    int32_t block = 0;
    while (bcf_next_line(bcf_fri)) {
        /// @note this could maybe be improved by removing the dup / destroy
        bcf1_t *rec = bcf_dup(bcf_fri.line);
        bcf_unpack(rec, BCF_UN_STR);
        rec->n_sample = 1;

        size_t idx = new_version ? line : pos;
        //if (line and ((line % BLOCK_LENGTH) == 0)) { // New version
        //if (pos and ((pos % BLOCK_LENGTH) == 0)) { // Old version
        if (idx and ((idx % BLOCK_LENGTH) == 0)) {
            block++;
            offset = 0; // New block
        }
        if (offset >> BM_BLOCK_BITS) {
            std::cerr << "Offset cannot be represented on " << BM_BLOCK_BITS << " bits !" << std::endl;
            throw "Variant BCF generation error, BM bits";
        }
        int32_t _ = block << BM_BLOCK_BITS | offset;

        bcf_update_format_int32(hdr, rec, "BM", &_, 1);
        ret = bcf_write1(fp, hdr, rec);
        bcf_destroy(rec);
        if (bcf_fri.line->n_allele) {
            pos += bcf_fri.line->n_allele-1;
            offset += bcf_fri.line->n_allele-1;
        }
        line++;
    }

    // Close everything
    hts_close(fp);
    bcf_hdr_destroy(hdr);
    destroy_bcf_file_reader(bcf_fri);

    return pos;
}

std::string get_entry_from_bcf(const std::string& filename, const char *entry_key) {
    // Input file
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, filename);

    auto hrec = bcf_hdr_get_hrec(bcf_fri.sr->readers[0].header, BCF_HL_GEN, entry_key, NULL, NULL);
    if (!hrec) {
        destroy_bcf_file_reader(bcf_fri);
        throw "Could not extract entry from header";
    }
    std::string value(hrec->value);

    destroy_bcf_file_reader(bcf_fri);
    return value;
}

std::vector<std::vector<bool> > extract_phase_vectors(const std::string& ifname) {
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, ifname);

    const int32_t PLOIDY = 2;

    std::vector<std::vector<bool> > phase_vectors(bcf_fri.n_samples);

    while(bcf_next_line(bcf_fri)) {
        // Unpack the line and get genotypes
        bcf_unpack(bcf_fri.line, BCF_UN_STR);
        bcf_fri.ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.size_gt_arr));
        int line_max_ploidy = bcf_fri.ngt / bcf_fri.n_samples;

        // Check ploidy, only support diploid for the moment
        if (line_max_ploidy != PLOIDY) {
            std::cerr << "[ERROR] Ploidy of samples is different than 2" << std::endl;
            exit(-1); // Change this
        }

        for (size_t sample_index = 0; sample_index < bcf_fri.n_samples; ++sample_index) {
            int difference = bcf_gt_allele(bcf_fri.gt_arr[sample_index*2+1]) - bcf_gt_allele(bcf_fri.gt_arr[sample_index*2]);
            if (difference > 0) {
                phase_vectors[sample_index].push_back(1);
            } else if (difference < 0) {
                phase_vectors[sample_index].push_back(0);
            }
        }
    }

    // Close / Release ressources
    destroy_bcf_file_reader(bcf_fri);

    return phase_vectors;
}

size_t compute_phase_switch_errors(const std::vector<bool>& testseq, const std::vector<bool>& refseq) {
    if (testseq.size() != refseq.size()) {
        std::cerr << "Seqs are not the same size, cannot compute" << std::endl;
        throw "Size problem";
    }

    size_t phase_switch_error_counter = 0;
    for (size_t i = 1; i < testseq.size(); ++i) {
        if (testseq[i-1] xor testseq[i] xor refseq[i-1] xor refseq[i]) {
            phase_switch_error_counter++;
        }
    }
    return phase_switch_error_counter;
}

void compute_phase_switch_errors(const std::string& testFile, const std::string& refFile) {
    auto test = extract_phase_vectors(testFile);
    auto ref  = extract_phase_vectors(refFile);
    // This array is if we want some more statistics and per sample info
    std::vector<size_t> phase_switch_errors(test.size(), 0);

    if (test.size() != ref.size()) {
        std::cerr << "Results are not the same size, cannot compute" << std::endl;
        throw "Size problem";
    }

    size_t total_errors = 0;
    size_t total_switches = 0;
    for (size_t sample_index = 0; sample_index < test.size(); sample_index++) {
        phase_switch_errors[sample_index] = compute_phase_switch_errors(test[sample_index], ref[sample_index]);
        if (sample_index < 1000) {
            std::cerr << "Sample " << sample_index << " : " << phase_switch_errors[sample_index]
                      << " errors on " << test[sample_index].size() << " sites" << std::endl;
        }
        total_errors += phase_switch_errors[sample_index];
        total_switches += test[sample_index].size() ? test[sample_index].size()-1 : 0;
    }

    double percentage = (double)total_errors / (double)total_switches * 100.0;

    std::cerr << "The error percentage is : " << percentage << " %" << std::endl;
}

int32_t seek_default_phased(const std::string& filename, size_t limit) {
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, filename);
    std::vector<size_t> counts(2, 0);
    while(bcf_next_line(bcf_fri) && limit--) {
        // Unpack the line and get genotypes
        bcf_unpack(bcf_fri.line, BCF_UN_STR);
        bcf_fri.ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.size_gt_arr));
        int line_max_ploidy = bcf_fri.ngt / bcf_fri.n_samples;

        if (line_max_ploidy == 1) {
            destroy_bcf_file_reader(bcf_fri);
            return 0; // Not phased, does not matter for haploid
        }

        for (size_t i = 0; i < bcf_fri.n_samples; ++i) {
            counts[bcf_gt_is_phased(bcf_fri.gt_arr[i*line_max_ploidy+1])]++;
        }
    }
    destroy_bcf_file_reader(bcf_fri);
    if (counts[0] > counts[1]) {
        return 0; // Majority is unphased
    } else {
        return 1; // Majority is phased
    }
}

size_t seek_max_ploidy_from_first_entry(const std::string& filename) {
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, filename);

    size_t max_ploidy = 0;
    if (bcf_fri.n_samples == 0) {
        std::cerr << "The file " << filename << "has no samples" << std::endl;
        destroy_bcf_file_reader(bcf_fri);
        return 0;
    }
    if (bcf_next_line(bcf_fri) == 0) {
        std::cerr << "The file " << filename << "has no entries" << std::endl;
        destroy_bcf_file_reader(bcf_fri);
        return 0;
    } else {
        // Unpack the line and get genotypes
        bcf_unpack(bcf_fri.line, BCF_UN_STR);
        bcf_fri.ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.size_gt_arr));
        max_ploidy = bcf_fri.ngt / bcf_fri.n_samples;
    }

    destroy_bcf_file_reader(bcf_fri);

    return max_ploidy;
}

bool file_has_no_samples(const std::string& filename) {
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, filename);
    if (bcf_fri.n_samples == 0) {
        std::cerr << "The file " << filename << " has no samples" << std::endl;
        destroy_bcf_file_reader(bcf_fri);
        return true;
    }
    destroy_bcf_file_reader(bcf_fri);
    return false;
}

bool file_has_no_entries(const std::string& filename) {
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, filename);
    if (bcf_next_line(bcf_fri) == 0) {
        std::cerr << "The file " << filename << " has no entries" << std::endl;
        destroy_bcf_file_reader(bcf_fri);
        return true;
    }
    destroy_bcf_file_reader(bcf_fri);
    return false;
}

#if 0
uint64_t get_number_of_record_from_indexed_xcf(std::string filename)
{
    const char **seq = NULL;
    hts_idx_t *idx = NULL;
    tbx_t *tbx = NULL;
    htsFile *fp = NULL;

    int tid, nseq = 0, ret = 0;
    bcf_hdr_t *hdr = NULL;

    uint64_t sum = 0;

    fp = hts_open(filename.c_str(), "r");
    if (!fp) {
        std::cerr << "Could not open " << filename << std::endl;
        throw "File error";
    }

    if (hts_get_format(fp)->format == vcf) {
        std::stringstream ss;
        ss << filename << ".tbi";
        tbx = tbx_index_load2(filename.c_str(), ss.str().c_str());
        if (!tbx) {
            std::cerr << "Could not load index for VCF: " << filename << std::endl;
            throw "File error";
        }
    } else if (hts_get_format(fp)->format == bcf) {
        std::stringstream ss;
        ss << filename << ".csi";
        idx = bcf_index_load2(filename.c_str(), ss.str().c_str());
        if (!idx) {
            std::cerr << "Could not load index for VCF: " << filename << std::endl;
            throw "File error";
        }
    } else {
        std::cerr << "Bad file type" << std::endl;
        throw "File error";
    }

    if (tbx) {
        seq = tbx_seqnames(tbx, &nseq);
    } else {
        nseq = hts_idx_nseq(idx);
    }

    for (tid=0; tid<nseq; tid++) {
        uint64_t records, v;
        hts_idx_get_stat(tbx ? tbx->idx : idx, tid, &records, &v);
        sum += records;
    }

    if (fp && hts_close(fp) != 0) {
        std::cerr << "Could not close file..." << std::endl;
    }
    if (tbx) {
        tbx_destroy(tbx);
    }
    if (idx) {
        hts_idx_destroy(idx);
    }
    return sum;
}
#endif