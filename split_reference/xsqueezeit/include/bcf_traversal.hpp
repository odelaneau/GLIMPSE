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

#ifndef __BCF_TRAVERSAL_HPP__
#define __BCF_TRAVERSAL_HPP__

#include "vcf.h"
#include "hts.h"

#include "xcf.hpp"

#include <random>

/// @todo Clean the usage of PLOIDY !

class BcfTraversal {
public:

    BcfTraversal() {}

    // Todo set region / samples

    void traverse(const std::string filename);

    virtual ~BcfTraversal() {}

protected:
    virtual void handle_bcf_file_reader() {}
    virtual void handle_bcf_line() {}

    bcf_file_reader_info_t bcf_fri;
    int ngt = 0;
    size_t line_max_ploidy = 0;
};

class BcfTransformer : public BcfTraversal {
public:
    BcfTransformer() : fp(NULL), hdr(NULL) {}
    virtual ~BcfTransformer() {}

    virtual void transform_header() {}
    virtual void transform_record() {}

    void transform(const std::string& ifname, const std::string& ofname);

protected:
    void handle_bcf_file_reader() override;
    void handle_bcf_line() override;

    bcf1_t *rec;
    htsFile* fp;
    bcf_hdr_t* hdr;
};

#include <random>
class BcfUnphaser : protected BcfTransformer {
public:
    BcfUnphaser() : rd(), gen(rd()), distrib(1, 2) {}
    virtual ~BcfUnphaser() {}

    void unphase_random(const std::string& ifname, const std::string& ofname) {
        transform(ifname, ofname);
    }

protected:
    const size_t PLOIDY = 2;
    void transform_record() override;

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen; // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib;
};

template <typename T>
class BcfFillMatrix : protected BcfTraversal {
public:
    BcfFillMatrix(std::vector<std::vector<T> >& matrix_ref) : matrix_ref(matrix_ref) {}
    virtual ~BcfFillMatrix() {}

    void fill_matrix_from_file(std::string filename) {
        traverse(filename);
    }
protected:
    void handle_bcf_file_reader() override {
        matrix_ref.clear();
        n_samples = bcf_fri.n_samples;
    }

    virtual void line_to_matrix() {
        matrix_ref.push_back(std::vector<T>(n_samples * PLOIDY, false));
        for (size_t i = 0; i < bcf_fri.n_samples * PLOIDY; ++i) {
            matrix_ref.back().at(i) = bcf_gt_allele(bcf_fri.gt_arr[i]);
        }
    }

    void handle_bcf_line() override {
        if (line_max_ploidy != PLOIDY) {throw "PLOIDY ERROR";}
        line_to_matrix();
    }

    const size_t PLOIDY = 2;
    std::vector<std::vector<T> >& matrix_ref;
    size_t n_samples = 0;
};

class BcfFillPhaseMatrix : public BcfFillMatrix<int8_t> {
public:
    BcfFillPhaseMatrix(std::vector<std::vector<int8_t> >& m) : BcfFillMatrix(m) {}
    virtual ~BcfFillPhaseMatrix() {}

protected:
    void line_to_matrix() override {
        matrix_ref.push_back(std::vector<int8_t>(n_samples, 0));
        for (size_t i = 0; i < bcf_fri.n_samples; ++i) {
            auto allele_0 = bcf_gt_allele(bcf_fri.gt_arr[i*2]);
            auto allele_1 = bcf_gt_allele(bcf_fri.gt_arr[i*2+1]);

            if (allele_0 == allele_1) {
                matrix_ref.back().at(i) = allele_0 ? 2 : 0;
            } else {
                matrix_ref.back().at(i) = allele_0 > allele_1 ? -1 : 1;
            }

        }
    }
};

template<typename T = bool>
class BcfMatrix {
public:
    BcfMatrix() {};
    BcfMatrix(std::string filename) : filename(filename) {
        BcfFillMatrix<T> bfm(matrix);
        bfm.fill_matrix_from_file(filename);
    }
    virtual ~BcfMatrix() {}

    template<const bool verbose = false>
    bool compare(const BcfMatrix& other) const {
        if (matrix.size() != other.matrix.size()) {
            if (verbose) {
                std::cerr << "Matrices differ in size" << std::endl;
            }
            return false;
        }

        for (size_t i = 0; i < matrix.size(); ++i) {
            if (matrix.at(i).size() != other.matrix.at(i).size()) {
                if (verbose) {
                    std::cerr << "Matrices differ in size at " << i << std::endl;
                }
                return false;
            }
            for (size_t j = 0; j < matrix.at(i).size(); ++j) {
                if (matrix.at(i).at(j) != other.matrix.at(i).at(j)) {
                    if (verbose) {
                        std::cerr << "Matrices differ at " << i << "," << j << std::endl;
                    }
                    return false;
                }
            }
        }

        return true;
    }

    bool operator==(const BcfMatrix& other) const {
        return compare(other);
    }

    std::string get_original_filename() const {return filename;}
    const std::vector<std::vector<T> >& get_matrix_const_ref() const {return matrix;}
    std::vector<std::vector<T> >& get_matrix_ref() {return matrix;}

    void inject_phase_switch_errors(double phase_switch_probability) {
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        const size_t RANGE_MAX = 1000000;
        std::uniform_int_distribution<> distrib(0, RANGE_MAX);

        std::vector<bool> samples_phase(matrix.at(0).size(), false);
        for (size_t i = 0; i < matrix.size(); ++i) {
            for (size_t j = 0; j < matrix.at(0).size() / 2; ++j) {
                bool allele_1 = matrix.at(i).at(2*j);
                bool allele_2 = matrix.at(i).at(2*j+1);
                if (allele_1 != allele_2) {
                    if (distrib(gen) < phase_switch_probability*RANGE_MAX) {
                        samples_phase.at(j) = !samples_phase.at(j); // Toggle phase
                    }
                }
                if (samples_phase.at(j)) { // Phase switch
                    matrix.at(i).at(2*j) = allele_2;
                    matrix.at(i).at(2*j+1) = allele_1;
                }
            }
        }
    }
protected:
    std::string filename;
    std::vector<std::vector<T> > matrix;
};

class BcfPhaseMatrix : public BcfMatrix<int8_t> {
    BcfPhaseMatrix(std::string filename) {
        this->filename = filename;
        BcfFillPhaseMatrix bfm(matrix);
        bfm.fill_matrix_from_file(filename);
    }
};

class BcfWriteMatrix : protected BcfTransformer {
public:
    BcfWriteMatrix(const BcfMatrix<bool>& bfm) : filename(bfm.get_original_filename()), matrix(bfm.get_matrix_const_ref()) {}
    virtual ~BcfWriteMatrix() {}

    void write(std::string ofname) {
        record_number = 0;
        transform(filename, ofname);
    }
protected:
    std::string filename;
    const std::vector<std::vector<bool> >& matrix;
    void transform_record() {
        for (size_t i = 0; i < bcf_fri.n_samples * line_max_ploidy; ++i) {
            bcf_fri.gt_arr[i] = matrix.at(record_number).at(i) ? bcf_gt_phased(1) : bcf_gt_phased(0);
        }
        record_number++;
    }
    size_t record_number = 0;
};

#endif /* __BCF_TRAVERSAL_HPP__ */