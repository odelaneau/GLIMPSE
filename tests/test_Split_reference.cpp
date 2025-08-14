#include <gtest/gtest.h>
#include "TestFile.hpp"
#include "spliter_header.h"

// Map folder
#ifndef MAP_FOLDER
#define MAP_FOLDER
#endif

namespace fs = std::filesystem;

TEST(Split_reference, test_split_ref_file)
{   
    // Create posfile
    fs::path simple_file = TestFile().get_tmp_vcf_file("simple_split_file");

    // Create args string
    std::vector<std::string> args{
        "--input-region", "1:1-4611812489108143365",
        "--output-region", "1:1-225061",
        "--output", "binary_reference_panel", 
        "--reference", simple_file.string(), 
        "--map", std::string(MAP_FOLDER) + "/chr1.b38.gmap.gz"
    };
    spliter().phase(args);

    fs::remove(simple_file);
}