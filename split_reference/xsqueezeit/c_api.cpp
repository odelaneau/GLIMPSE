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

/**
 * @file c_api.cpp
 * @author Rick Wertenbroek
 * @date 23 Sep 2021
 * @brief C API for xSqueezeIt
 *
 */

#include "c_api.h"

#include "xsi_mixed_vcf.hpp"

extern "C" {
    c_xcf *c_xcf_new() {
        return (c_xcf *)new Xcf();
    }

    // void c_xcf_add_reader(c_xcf *x, bcf_sr_t* reader, int reader_id) {
    //     reinterpret_cast<Xcf*>(x)->add_reader(reader, reader_id);
    // }

    void c_xcf_add_readers(c_xcf *x, bcf_srs_t* readers) {
        reinterpret_cast<Xcf*>(x)->add_readers(readers);
    }

    void c_xcf_update_readers(c_xcf *x, bcf_srs_t* readers) {
        reinterpret_cast<Xcf*>(x)->add_readers(readers);
    }

    const char* c_xcf_sample_name(c_xcf *x, int reader_id, const bcf_hdr_t* hdr, int sample_id) {
        return reinterpret_cast<Xcf*>(x)->sample_name(reader_id, hdr, sample_id);
    }

    /// @todo make this more robust / cleaner
    int c_xcf_nsamples(const char* fname) {
        try {
            auto filename = Accessor::get_filename_from_variant_file(fname);
            if (fs::exists(filename)) {
                return Accessor(filename).get_number_of_samples();
            }
        } catch (...) {}

        bcf_srs_t *sr = bcf_sr_init();
        int ret = bcf_sr_add_reader(sr, fname);
        if (!ret) {
            bcf_sr_destroy(sr);
            return 0;
        }
        int nsamples = bcf_hdr_nsamples(sr->readers[0].header);
        bcf_sr_destroy(sr);
        return nsamples;
    }

    int __c__xcf__get__genotypes__void(c_xcf *x, int reader_id, const bcf_hdr_t *hdr, bcf1_t *line, void **dst, int *ndst) {
        return reinterpret_cast<Xcf*>(x)->get_genotypes(reader_id, hdr, line, dst, ndst);
    }

    void c_xcf_delete(c_xcf *x) {
        delete reinterpret_cast<Xcf*>(x);
    }

}