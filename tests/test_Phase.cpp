#include "TestFile.hpp"
#include "caller/caller_header.h"
#include <gtest/gtest.h>

#ifndef BIN_DIR
#define BIN_DIR
#endif

TEST(Phase, test_impute_mmini_file)
{
    // Find input from previous tests
    std::string reference = TestFile().globPattern(std::string(BIN_DIR) + "/binary_reference_panel*.bin");
    std::string bam_list = TestFile().globPattern(std::string(BIN_DIR) + "/chunks.*.txt");

    // Create args string
    std::vector<std::string> args{
        "--bam-list", bam_list,
        "--reference", reference,
        "--output", "imputed_glimpse2_mini.bcf", 
        "--threads", "4"
    };

    caller().phase(args);
}