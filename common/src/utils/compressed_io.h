/*******************************************************************************
 * Copyright (C) 2022-2023 Simone Rubinacci
 * Copyright (C) 2022-2023 Olivier Delaneau
 *
 * MIT Licence
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

#ifndef _COMPRESSED_IO_H
#define _COMPRESSED_IO_H

// STL includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

// Boost includes
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

/**
 * @brief Input file stream wrapper supporting gzip and bzip2 decompression transparently
 */
class input_file : public boost::iostreams::filtering_istream {
protected:
    std::ifstream file_descriptor;

public:
    explicit input_file(const std::string& filename) {
        // Detect compression by file extension
        const std::string ext = filename.substr(filename.find_last_of(".") + 1);
        if (ext == "gz") {
            file_descriptor.open(filename.c_str(), std::ios::in | std::ios::binary);
            push(boost::iostreams::gzip_decompressor());
        } else if (ext == "bz2") {
            file_descriptor.open(filename.c_str(), std::ios::in | std::ios::binary);
            push(boost::iostreams::bzip2_decompressor());
        } else {
            file_descriptor.open(filename.c_str());
        }
        // Push underlying file descriptor to filtering stream if open succeeded
        if (!file_descriptor.fail()) {
            push(file_descriptor);
        }
    }

    ~input_file() {
        close();
    }

    /**
     * @return true if underlying file failed to open
     */
    bool fail() const {
        return file_descriptor.fail();
    }

    /**
     * Close the stream and underlying file descriptor
     */
    void close() {
        if (!file_descriptor.fail()) {
            if (!empty()) reset();
            file_descriptor.close();
        }
    }
};


/**
 * @brief Output file stream wrapper supporting gzip and bzip2 compression transparently
 */
class output_file : public boost::iostreams::filtering_ostream {
protected:
    std::ofstream file_descriptor;

public:
    output_file() = default;

    explicit output_file(const std::string& filename) {
        open(filename);
    }

    ~output_file() {
        close();
    }

    bool fail() const {
        return file_descriptor.fail();
    }

    void close() {
        if (!file_descriptor.fail()) {
            if (!empty()) reset();
            file_descriptor.close();
        }
    }

    /**
     * Open a file with optional compression based on file extension.
     * @param filename Path of the file to open
     */
    void open(const std::string& filename) {
        const std::string ext = filename.substr(filename.find_last_of(".") + 1);
        if (ext == "gz") {
            file_descriptor.open(filename.c_str(), std::ios::out | std::ios::binary);
            push(boost::iostreams::gzip_compressor());
        } else if (ext == "bz2") {
            file_descriptor.open(filename.c_str(), std::ios::out | std::ios::binary);
            push(boost::iostreams::bzip2_compressor());
        } else {
            file_descriptor.open(filename.c_str());
        }
        if (!file_descriptor.fail()) {
            push(file_descriptor);
        }
    }
};

#endif
