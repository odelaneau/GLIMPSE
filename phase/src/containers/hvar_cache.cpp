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

#include <containers/hvar_cache.h>

#include <cstring>
#include <cerrno>
#include <fstream>
#include <vector>
#include <filesystem>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <zlib.h>

std::string default_hvar_cache_path(const std::string & reference_filename) {
	return reference_filename + ".hvarT";
}

bool stat_file_fingerprint(const std::string & path, uint64_t & size, int64_t & mtime) {
	struct stat st{};
	if (::stat(path.c_str(), &st) != 0) return false;
	size = (uint64_t)st.st_size;
	mtime = (int64_t)st.st_mtime;
	return true;
}

bool read_hvar_cache_header(const std::string & cache_path, hvar_cache_header & hdr, std::string & err_msg) {
	std::ifstream ifs(cache_path, std::ios::binary);
	if (!ifs.good()) {
		err_msg = "cache file does not exist or cannot be opened: " + cache_path;
		return false;
	}
	ifs.read(reinterpret_cast<char*>(&hdr), sizeof(hdr));
	if (!ifs.good()) {
		err_msg = "cache file is truncated (shorter than its header): " + cache_path;
		return false;
	}
	if (std::memcmp(hdr.magic, HVAR_CACHE_MAGIC, sizeof(HVAR_CACHE_MAGIC)) != 0) {
		err_msg = "cache file magic mismatch (not a GLIMPSE2 hvar-cache file): " + cache_path;
		return false;
	}
	if (hdr.format_version != HVAR_CACHE_FORMAT_VERSION) {
		err_msg = "cache file format version (" + std::to_string(hdr.format_version) + ") is not supported by this build (expected " + std::to_string(HVAR_CACHE_FORMAT_VERSION) + "): " + cache_path;
		return false;
	}
	return true;
}

bool hvar_cache_matches_source(const hvar_cache_header & hdr, const std::string & reference_filename, uint64_t expected_n_ref_haps, uint64_t expected_n_com_sites) {
	uint64_t size; int64_t mtime;
	if (!stat_file_fingerprint(reference_filename, size, mtime)) return false;
	if (hdr.source_file_size != size || hdr.source_mtime != mtime) return false;
	if (hdr.n_ref_haps != expected_n_ref_haps || hdr.n_com_sites != expected_n_com_sites) return false;
	return true;
}

void * mmap_hvar_cache_blob_raw(const std::string & cache_path, const hvar_cache_header & hdr, std::string & err_msg) {
	if (hdr.blob_bytes == 0) {
		err_msg = "cache blob is empty";
		return nullptr;
	}

	long page_size = sysconf(_SC_PAGESIZE);
	if (page_size <= 0 || (HVAR_CACHE_BLOB_OFFSET % (unsigned long)page_size) != 0) {
		err_msg = "cache blob offset (" + std::to_string(HVAR_CACHE_BLOB_OFFSET) + ") is not a multiple of this system's page size (" + std::to_string(page_size) + ")";
		return nullptr;
	}

	int fd = ::open(cache_path.c_str(), O_RDONLY);
	if (fd < 0) {
		err_msg = "could not open cache file for mmap: " + cache_path + " (" + std::strerror(errno) + ")";
		return nullptr;
	}

	void * base = ::mmap(nullptr, hdr.blob_bytes, PROT_READ, MAP_PRIVATE, fd, (off_t)HVAR_CACHE_BLOB_OFFSET);
	//Safe to close once mmap() has returned successfully; the mapping stays valid
	//independent of the file descriptor's lifetime.
	::close(fd);
	if (base == MAP_FAILED) {
		err_msg = "mmap() of cache blob failed: " + std::string(std::strerror(errno));
		return nullptr;
	}

	//Haplotype access during state selection is scattered across the file (though
	//contiguous per-haplotype), so disable readahead.
	::madvise(base, hdr.blob_bytes, MADV_RANDOM);
	return base;
}

void unmap_hvar_cache_blob_raw(void * base, unsigned long length) {
	if (base) ::munmap(base, length);
}

bool write_hvar_cache(const std::string & cache_path, const std::string & reference_filename, uint64_t n_ref_haps, uint64_t n_com_sites, const bitmatrix & transposed, std::string & err_msg) {
	uint64_t src_size; int64_t src_mtime;
	if (!stat_file_fingerprint(reference_filename, src_size, src_mtime)) {
		err_msg = "could not stat source reference panel file: " + reference_filename;
		return false;
	}

	hvar_cache_header hdr{};
	std::memset(&hdr, 0, sizeof(hdr));
	std::memcpy(hdr.magic, HVAR_CACHE_MAGIC, sizeof(HVAR_CACHE_MAGIC));
	hdr.format_version = HVAR_CACHE_FORMAT_VERSION;
	hdr.source_file_size = src_size;
	hdr.source_mtime = src_mtime;
	hdr.n_ref_haps = n_ref_haps;
	hdr.n_com_sites = n_com_sites;
	hdr.blob_bytes = transposed.n_bytes;
	hdr.blob_crc32 = (uint32_t)crc32(0L, transposed.bytes, (uInt)transposed.n_bytes);

	const std::string tmp_path = cache_path + ".tmp";
	std::ofstream ofs(tmp_path, std::ios::binary | std::ios::trunc);
	if (!ofs.good()) {
		err_msg = "could not open temp cache file for writing: " + tmp_path;
		return false;
	}

	//Zero-padded, fixed-size header region so the blob that follows starts at a
	//page-aligned offset regardless of sizeof(hvar_cache_header).
	std::vector<char> header_block((size_t)HVAR_CACHE_BLOB_OFFSET, 0);
	std::memcpy(header_block.data(), &hdr, sizeof(hdr));
	ofs.write(header_block.data(), (std::streamsize)header_block.size());
	ofs.write(reinterpret_cast<const char*>(transposed.bytes), (std::streamsize)transposed.n_bytes);
	ofs.close();
	if (!ofs) {
		err_msg = "error writing cache blob to: " + tmp_path;
		return false;
	}

	std::error_code ec;
	std::filesystem::rename(tmp_path, cache_path, ec);
	if (ec) {
		err_msg = "could not rename temp cache file [" + tmp_path + "] into place [" + cache_path + "]: " + ec.message();
		return false;
	}
	return true;
}
