#pragma once

#include <string>
#include <map>
#include <filesystem>
#include <glob.h>

class TestFile
{
    public:
    std::filesystem::path get_tmp_file(std::string file_name);
    std::filesystem::path create_tmp_file(std::string file_context, std::string prefix);
    std::filesystem::path get_tmp_vcf_file(std::string file_key);
    std::string globPattern(const std::string& pattern);
    std::string createBamFile(const std::string& output_filename);

    private:
    static std::map<std::string, std::string> file_menu;

};