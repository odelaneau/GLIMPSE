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

#ifndef _HVAR_CACHE_H
#define _HVAR_CACHE_H

#include <cstdint>
#include <string>
#include <containers/bitmatrix.h>

//A phase-generated, phase-only cache file holding a haplotype-major (transposed)
//copy of a binary reference panel's HvarRef bitmatrix, built by
//`GLIMPSE2_phase --build-hvar-cache` from an already-built .bin file. This is a
//brand-new file format with no relation to (and no backward-compatibility
//constraints from) the .bin reference-panel format -- the source .bin is only ever
//read, never rewritten. Once built, normal `phase` runs mmap() this file instead of
//materializing the .bin's row-major HvarRef in RAM, which is what actually reduces
//memory usage on very large panels (see phase.md / the design notes in
//caller_initialise.cpp for why the row-major layout can't just be mmap'd directly).
struct hvar_cache_header {
	char     magic[16];         //"GLIMPSE2HVCACHE" + NUL, not further NUL-padded beyond that
	uint32_t format_version;
	uint32_t reserved0;
	//Fingerprint of the source .bin this cache was built from, used to detect a
	//stale cache (source panel replaced/rebuilt) cheaply, without hashing the whole
	//multi-GB .bin file.
	uint64_t source_file_size;
	int64_t  source_mtime;       //seconds since epoch, from stat()
	uint64_t n_ref_haps;          //redundant cross-check against the .bin's own scalars
	uint64_t n_com_sites;
	uint64_t blob_bytes;          //exact byte length of the haplotype-major blob
	uint32_t blob_crc32;
	uint32_t reserved1;
};

inline constexpr char HVAR_CACHE_MAGIC[16] = "GLIMPSE2HVCACHE";
inline constexpr uint32_t HVAR_CACHE_FORMAT_VERSION = 1;
//Fixed constant (not the build host's OS page size) so the blob offset is a multiple
//of every realistic page size across the platforms GLIMPSE2 supports (Linux
//x86_64/aarch64, macOS Apple Silicon: 4K/16K/64K observed in practice).
inline constexpr uint64_t HVAR_CACHE_BLOB_OFFSET = 65536;

//Rounds a dimension up to the next multiple of 8, matching the padding rule used
//throughout bitmatrix (allocate/reallocate): row/col counts must address whole bytes.
inline uint64_t hvar_cache_round8(uint64_t x) {
	return x + ((x % 8) ? (8 - (x % 8)) : 0);
}

//Default cache path for a given reference panel path: "<reference>.hvarT" next to it.
std::string default_hvar_cache_path(const std::string & reference_filename);

//stat() a file for its size + mtime. Returns false if the file doesn't exist / can't be stat'd.
bool stat_file_fingerprint(const std::string & path, uint64_t & size, int64_t & mtime);

//Read + sanity-check (magic/version) a cache file's header. Returns false (with err_msg
//set) if the file doesn't exist, is too short, or fails the magic/version check.
bool read_hvar_cache_header(const std::string & cache_path, hvar_cache_header & hdr, std::string & err_msg);

//True if `hdr` (already read from a candidate cache file) matches the CURRENT state of
//`reference_filename` (size+mtime) and the given panel scalars. False (cache is stale
//or mismatched) otherwise -- never throws, just returns false so the caller can fall
//back to a full load.
bool hvar_cache_matches_source(const hvar_cache_header & hdr, const std::string & reference_filename, uint64_t expected_n_ref_haps, uint64_t expected_n_com_sites);

//mmap() the cache file's blob region (PROT_READ, MAP_PRIVATE, MADV_RANDOM). Returns the
//raw mapped pointer and its mapped length (== hdr.blob_bytes) on success, or nullptr
//(with err_msg set) on failure. Deliberately does NOT adopt the result into a bitmatrix
//-- the caller (caller::read_binary_reference_panel) must be able to attempt this
//*before* discarding whatever bytes the destination bitmatrix currently holds, so a
//failure here can fall back to a normal full load instead of being treated as fatal.
//The caller adopts the returned pointer via bitmatrix::adopt_mmap() once it's safe to
//do so (see caller_initialise.cpp).
void * mmap_hvar_cache_blob_raw(const std::string & cache_path, const hvar_cache_header & hdr, std::string & err_msg);

//Tears down a mapping obtained from mmap_hvar_cache_blob_raw() that is being abandoned
//without ever being adopted into a bitmatrix (e.g. because a subsequent step failed).
//No-op if base is nullptr. `length` must be the same value mmap_hvar_cache_blob_raw()
//mapped (i.e. hdr.blob_bytes).
void unmap_hvar_cache_blob_raw(void * base, unsigned long length);

//Build a cache file at `cache_path` from an already-transposed, in-memory
//haplotype-major bitmatrix (`transposed`: n_rows=n_ref_haps, n_cols=n_com_sites),
//fingerprinted against `reference_filename`. Writes atomically (temp file + rename).
//Returns false (with err_msg set) on failure.
bool write_hvar_cache(const std::string & cache_path, const std::string & reference_filename, uint64_t n_ref_haps, uint64_t n_com_sites, const bitmatrix & transposed, std::string & err_msg);

#endif
