#include "TestFile.hpp"
#include "caller/caller_header.h"
#include <gtest/gtest.h>

TEST(Phase, test_impute_mmini_file)
{
    // Create args string
    std::vector<std::string> args{
        "--bam-list", "bams_1.0x.txt",
        "--reference", "binary_reference_panel_chr20_7702567_12266861.bin",
        "--output", "imputed_glimpse2_rp140k_1.0x_chr20_7702567_12266861.bcf", 
        "--threads", "4"
    };
}