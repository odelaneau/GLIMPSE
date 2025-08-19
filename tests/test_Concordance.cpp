#include "checker_header.h"
#include "TestFile.hpp"

#include <gtest/gtest.h>

namespace fs = std::filesystem;

// Checking tmp file working normally.
TEST(Initial, test_create_tmp_vcf_file) {
  // Expected context
  std::string expected_context =  "##fileformat=VCFv4.2\n"
                                  "##contig=<ID=1,length=1000000>\n"
                                  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                                  "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype Dosage\">\n"
                                  "##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor Allele Frequency\">\n"
                                  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\tsample3\n"
                                  "1\t123456\t.\tA\tG\t100\tPASS\tMAF=0.1\tGT:DS\t0|0:0.0\t0|1:1.0\t1|1:2.0\n"
                                  "1\t123789\t.\tC\tT\t99\tPASS\tMAF=0.2\tGT:DS\t0|1:0.8\t0|0:0.1\t1|1:1.9\n"
                                  "1\t200123\trs123\tG\tA\t300\tPASS\tMAF=0.3\tGT:DS\t1|1:0.2\t0|1:0.7\t1|0:1.1\n"
                                  "1\t250000\t.\tT\tC\t150\tPASS\tMAF=0.4\tGT:DS\t1|1:1.95\t0|0:0.05\t0|1:1.2\n"
                                  "1\t260000\t.\tA\tC\t150\tPASS\tMAF=0.5\tGT:DS\t0|1:1.55\t1|0:0.15\t1|1:1.5\n";

  // Get context of tmp file
  fs::path vcf_path = TestFile().get_tmp_file("simple_pos_file");
  std::cout << "Temp VCF file created at: " << vcf_path << '\n';
  std::ifstream actual_file_flow(vcf_path);
  std::string actual_context((std::istreambuf_iterator<char>(actual_file_flow)), std::istreambuf_iterator<char>());
  actual_file_flow.close();

  // Clear tmp file
  fs::remove(vcf_path);

  // Check Expected context and actual context
  EXPECT_EQ(actual_context, expected_context);
}

TEST(Concordance, test_vcf_vs_itself)
{
  // Call tmp vcf file
  fs::path vcf_path = TestFile().get_tmp_vcf_file("simple_file");

  // Create sample file
  fs::path sample_list = TestFile().create_tmp_file("sample1\nsample2\nsample3", "_sample_list.txt");

  // Create posfile
  fs::path pos_file = TestFile().get_tmp_vcf_file("simple_pos_file");

  // Create lst filr
  std::string lst_context {"1\t" + pos_file.string() + "\t" +  vcf_path.string() + "\t" + vcf_path.string()};
  fs::path lst_file = TestFile().create_tmp_file(lst_context, "_lst_file.lst");

  // Create args string
  std::vector<std::string> args{
      "--input", lst_file.string(),
      "--gt-val",
      "--af-tag", "MAF",
      "--samples", sample_list.string(),
      "--bins", "0.0", "0.1", "0.2", "0.3", "0.4", "0.5",
      "--output", "test_output"
  };

  // Run concordance
  checker().check(args);

  // Clear tmp file
  fs::remove(vcf_path);
  fs::remove(sample_list);
  fs::remove(pos_file);
  fs::remove(lst_file);

}