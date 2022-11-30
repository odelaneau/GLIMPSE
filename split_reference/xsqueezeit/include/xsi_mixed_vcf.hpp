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
 * @file xsi_mixed_vcf.hpp
 * @author Rick Wertenbroek
 * @date 23 Sep 2021
 * @brief Class to handle mixed xSqueezeIt and VCF/BCF files in synced readers
 *
 */
#ifndef __XSI_MIXED_VCF_HPP__
#define __XSI_MIXED_VCF_HPP__

#include "vcf.h"

#include <vector>
#include <memory>

#include "accessor.hpp"

/**
 * @brief Helper class for mixed xSqueezeIt and VCF/BCF support
 * 
 */
class Xcf {
public:
    Xcf();

    void add_readers(bcf_srs_t* readers);
    void update_readers(bcf_srs_t* readers);

    const char *sample_name(int reader_id, const bcf_hdr_t *hdr, int sample_id);
    int get_genotypes(int reader_id, const bcf_hdr_t *hdr, bcf1_t *line, void **dst, int *ndst);
protected:
    void add_reader(bcf_sr_t* reader, int reader_id);

    bool reader_file_is_xsi(bcf_sr_t* reader);

    class XcfInternalEntry {
    public:
        bool is_xsi = false;
        std::unique_ptr<Accessor> accessor;
    };

    std::vector<XcfInternalEntry> entries;
};

#endif /* __XSI_MIXED_VCF_HPP__ */