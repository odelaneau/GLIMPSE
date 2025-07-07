#include <gtest/gtest.h>
#include "TestFile.hpp"
#include "caller_header.h"

// Map folder
#ifndef MAP_FOLDER
#define MAP_FOLDER
#endif

namespace fs = std::filesystem;

TEST(Split_reference, test_split_ref_file)
{   
    // Create posfile
    fs::path pos_file = TestFile().get_tmp_vcf_file("simple_pos_file");

    // Create args string
    std::vector<std::string> args{
        "--input", pos_file.string(),
        "--region", "chr1",
        "--output", "chunks.chr1.txt",
        "--map", MAP_FOLDER + "/chr1.b38.gmap.gz"
    };
    caller().phase(args);
}