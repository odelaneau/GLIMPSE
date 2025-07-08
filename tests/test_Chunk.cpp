#include "chunker_header.h"
#include "TestFile.hpp"

#include <gtest/gtest.h>

#ifndef MAP_FOLDER
#define MAP_FOLDER
#endif

namespace fs = std::filesystem;

TEST(Chunk, test_chunk_ref_file)
{
    // Create posfile
    fs::path pos_file = TestFile().get_tmp_vcf_file("chunking_file");

    // Create args string
    std::vector<std::string> args{
        "--input", pos_file.string(),
        "--region", "chr1",
        "--sequential",
        "--window-count", "2",
        "--output", "chunks.chr1.txt",
        "--map", std::string(MAP_FOLDER) + "/chr1.b38.gmap.gz",
    };
    chunker().chunk(args);

    fs::remove(pos_file);
}