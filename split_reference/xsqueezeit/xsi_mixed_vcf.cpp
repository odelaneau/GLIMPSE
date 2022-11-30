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
 * @file xsi_mixed_vcf.cpp
 * @author Rick Wertenbroek
 * @date 23 Sep 2021
 * @brief Class to handle mixed xSqueezeIt and VCF/BCF files in synced readers
 *        Implementation of methods
 *
 */

#include "xsi_mixed_vcf.hpp"

#include "fs.hpp"
#include "make_unique.hpp"

Xcf::Xcf() {}

/*
 * Simple check if the bcf has an associated xSqueezeIt file
 * This is not very robust but will do for the moment
 * 
 */
bool Xcf::reader_file_is_xsi(bcf_sr_t* reader) {
    try {
        auto fname = Accessor::get_filename_from_variant_file(reader->fname);
        if (fs::exists(fname)) {
            return true;
        } else {
            return false;
        }
    } catch (...) {
        return false;
    }
}

void Xcf::add_reader(bcf_sr_t* reader, int reader_id) {
    if (reader_id >= entries.size()) {
        entries.resize(reader_id+1);
    }

    // Check if reader is for xSqueezeit file
    if (reader_file_is_xsi(reader)) {
        entries[reader_id].is_xsi = true;
        auto fname = Accessor::get_filename_from_variant_file(reader->fname);
        entries[reader_id].accessor = make_unique<Accessor>(fname);
    } else {
        entries[reader_id].is_xsi = false;
        entries[reader_id].accessor = nullptr;
    }
}

void Xcf::add_readers(bcf_srs_t* readers) {
    for (int i = 0; i < readers->nreaders; ++i) {
        add_reader(&(readers->readers[i]), i);
    }
}

void Xcf::update_readers(bcf_srs_t* readers) {
    add_readers(readers);
}

const char *Xcf::sample_name(int reader_id, const bcf_hdr_t* hdr, int sample_id) {
    if (entries[reader_id].is_xsi) {
        return entries[reader_id].accessor->get_sample_list()[sample_id].c_str();
    } else {
        return hdr->samples[sample_id];
    }
}

int Xcf::get_genotypes(int reader_id, const bcf_hdr_t *hdr, bcf1_t *line, void **dst, int *ndst) {
    if (entries[reader_id].is_xsi) {
        return entries[reader_id].accessor->get_genotypes(hdr, line, dst, ndst);
    } else {
        return bcf_get_genotypes(hdr, line, dst, ndst);
    }
}
