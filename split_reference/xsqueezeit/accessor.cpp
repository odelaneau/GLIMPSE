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
#include "accessor.hpp"

Accessor::Accessor(std::string& filename) : filename(filename) {
    std::fstream s(filename, s.binary | s.in);
    if (!s.is_open()) {
        std::cerr << "Failed to open file " << filename << std::endl;
        throw "Failed to open file";
    }

    // Read the header
    s.read((char *)(&(this->header)), sizeof(header_t));

    // Check magic
    if ((header.first_magic != MAGIC) or (header.last_magic != MAGIC)) {
        std::cerr << "Bad magic" << std::endl;
        std::cerr << "Expected : " << MAGIC << " got : " << header.first_magic << ", " << header.last_magic << std::endl;
        throw "Bad magic";
    }

    // Check version
    if (header.version != 2 and header.version != 3) {
        if (header.version == 4) {
            //std::cerr << "Experimental version" << std::endl;
        } else {
            std::cerr << "Bad version" << std::endl;
            throw "Bad version";
        }
    }

    // Extract the sample list
    sample_list.clear();
    s.seekg(header.samples_offset);
    std::string str;
    while(std::getline(s, str, '\0').good() and (sample_list.size() < (header.hap_samples/header.ploidy))) {
        sample_list.push_back(str);
    }
    s.close();

    if (header.aet_bytes == 2) {
        if (header.version == 4) {
            internals = make_unique<AccessorInternalsNewTemplate<uint16_t> >(filename);
        } else {
            internals = make_unique<AccessorInternalsTemplate<uint16_t> >(filename);
        }
    } else if (header.aet_bytes == 4) {
        if (header.version == 4) {
            internals = make_unique<AccessorInternalsNewTemplate<uint32_t> >(filename);
        } else {
            internals = make_unique<AccessorInternalsTemplate<uint32_t> >(filename);
        }
    } else {
        std::cerr << "Unsupported access type" << std::endl;
        throw "Unsupported A_T";
    }

    values = (int*)malloc(sizeof(int));
}

Accessor::~Accessor() {
    if (values) {
        free(values);
    }
}