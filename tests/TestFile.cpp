#include "TestFile.hpp"
#include <fstream>
#include <htslib/bgzf.h>
#include <htslib/tbx.h>
#include <iostream>
#include <htslib/sam.h>
#include <cstring> 
#include <unistd.h>

namespace fs = std::filesystem;

std::map<std::string, std::string> TestFile::file_menu
{
    {"simple_file",
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1,length=1000000>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype Dosage\">\n"
        "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype Probabilities\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\tsample3\n"
        "1\t123456\t.\tA\tG\t100\tPASS\t.\tGT:DS:GP\t0|0:0.0:0.33,0.33,0.34\t0|1:1.0:0.33,0.33,0.34\t1|1:2.0:0.33,0.33,0.34\n"
        "1\t123789\t.\tC\tT\t99\tPASS\t.\tGT:DS:GP\t0|1:0.8:0.33,0.33,0.34\t0|0:0.1:0.33,0.33,0.34\t1|1:1.9:0.33,0.33,0.34\n"
        "1\t200123\trs123\tG\tA\t300\tPASS\t.\tGT:DS:GP\t1|1:0.2:0.33,0.33,0.34\t0|1:0.7:0.33,0.33,0.34\t1|0:1.1:0.33,0.33,0.34\n"
        "1\t250000\t.\tT\tC\t150\tPASS\t.\tGT:DS:GP\t1|1:1.95:0.33,0.33,0.34\t0|0:0.05:0.33,0.33,0.34\t0|1:1.2:0.33,0.33,0.34\n"
    },
    {"simple_pos_file",
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1,length=1000000>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype Dosage\">\n"
        "##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor Allele Frequency\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\tsample3\n"
        "1\t123456\t.\tA\tG\t100\tPASS\tMAF=0.1\tGT:DS\t0|0:0.0\t0|1:1.0\t1|1:2.0\n"
        "1\t123789\t.\tC\tT\t99\tPASS\tMAF=0.2\tGT:DS\t0|1:0.8\t0|0:0.1\t1|1:1.9\n"
        "1\t200123\trs123\tG\tA\t300\tPASS\tMAF=0.3\tGT:DS\t1|1:0.2\t0|1:0.7\t1|0:1.1\n"
        "1\t250000\t.\tT\tC\t150\tPASS\tMAF=0.4\tGT:DS\t1|1:1.95\t0|0:0.05\t0|1:1.2\n"
        "1\t260000\t.\tA\tC\t150\tPASS\tMAF=0.5\tGT:DS\t0|1:1.55\t1|0:0.15\t1|1:1.5\n"
    },
    {"chunking_file",
        "##fileformat=VCFv4.2\n"
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">\n"
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\n"
        "chr1\t123\t.\tA\tG\t100\tPASS\tAC=1;AN=4\tGT\t0/0\t0/1\n"
        "chr1\t789\trs123\tC\tT\t99\tPASS\tAC=2;AN=4\tGT\t0/1\t0/1\n"
        "chr1\t200123\t.\tG\tA\t50\tPASS\tAC=0;AN=4\tGT\t0/0\t0/0\n"
        "chr1\t250000\t.\tT\tC\t80\tPASS\tAC=3;AN=4\tGT\t1/1\t0/1\n"
        "chr1\t500123\t.\tG\tC\t50\tPASS\tAC=0;AN=4\tGT\t0/1\t1/0\n"
        "chr1\t550000\t.\tC\tT\t80\tPASS\tAC=3;AN=4\tGT\t0/1\t0/1\n"
    },
    {"simple_split_file",
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1,length=1000000>\n"
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele Count\">\n"
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total Number of Alleles\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype Dosage\">\n"
        "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype Probabilities\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\tsample3\n"
        "1\t123456\t.\tA\tG\t100\tPASS\tAC=3;AN=6\tGT:DS:GP\t0|0:0.0:0.33,0.33,0.34\t0|1:1.0:0.33,0.33,0.34\t1|1:2.0:0.33,0.33,0.34\n"
        "1\t123789\t.\tC\tT\t99\tPASS\tAC=3;AN=6\tGT:DS:GP\t0|1:0.8:0.33,0.33,0.34\t0|0:0.1:0.33,0.33,0.34\t1|1:1.9:0.33,0.33,0.34\n"
        "1\t200123\trs123\tG\tA\t300\tPASS\tAC=4;AN=6\tGT:DS:GP\t1|1:0.2:0.33,0.33,0.34\t0|1:0.7:0.33,0.33,0.34\t1|0:1.1:0.33,0.33,0.34\n"
        "1\t250000\t.\tT\tC\t150\tPASS\tAC=3;AN=6\tGT:DS:GP\t1|1:1.95:0.33,0.33,0.34\t0|0:0.05:0.33,0.33,0.34\t0|1:1.2:0.33,0.33,0.34\n"
    }

};


fs::path TestFile::get_tmp_file(std::string file_key)
{
    // Generate a unique file name
    fs::path temp_dir = fs::temp_directory_path();
    fs::path temp_file;
    do {
        temp_file = temp_dir / ("tmp_vcf_" + std::to_string(std::rand()) + ".vcf");
    } while (fs::exists(temp_file));

    // Write the VCF content
    std::ofstream out(temp_file);
    out << this->file_menu[file_key];
    out.close();

    return temp_file;

}

fs::path TestFile::create_tmp_file(std::string file_context, std::string prefix)
{
    // Generate a unique file name
    fs::path temp_dir = fs::temp_directory_path();
    fs::path temp_file;
    do {
        temp_file = temp_dir / ("tmp_test" + std::to_string(std::rand()) + prefix);
    } while (fs::exists(temp_file));

    // Write the VCF content
    std::ofstream out(temp_file);
    out << file_context;
    out.close();

    return temp_file;
}


fs::path TestFile::get_tmp_vcf_file(std::string file_key)
{
    // Step 1: Generate a unique .vcf.gz path
    fs::path temp_dir = fs::temp_directory_path();
    fs::path temp_file;
    do {
        temp_file = temp_dir / ("tmp_vcf_" + std::to_string(std::rand()) + ".vcf.gz");
    } while (fs::exists(temp_file));

    // Step 2: Write compressed (bgzipped) content using HTSlib's BGZF
    BGZF* out = bgzf_open(temp_file.c_str(), "w");
    if (!out) {
        throw std::runtime_error("Failed to open BGZF file for writing: " + temp_file.string());
    }

    const std::string& content = this->file_menu[file_key];
    if (bgzf_write(out, content.data(), content.size()) < 0) {
        bgzf_close(out);
        throw std::runtime_error("Failed to write BGZF content");
    }

    bgzf_close(out);

    // Step 3: Index the .vcf.gz file
    int ret = tbx_index_build(temp_file.c_str(), 0, &tbx_conf_vcf);
    if (ret != 0) {
        throw std::runtime_error("Failed to index VCF: " + temp_file.string());
    }

    return temp_file;
}


std::string TestFile::globPattern(const std::string& pattern)
{
    glob_t glob_result;
    std::string file_path;

    int return_value = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    if (return_value == 0) {
        file_path = glob_result.gl_pathv[0];
    } else {
        std::cerr << "glob() failed with return code: " << return_value << std::endl;
    }

    globfree(&glob_result);
    return file_path;
}

std::string TestFile::createBamFile(const std::string& basename = "mini") {
    // 1. Create temporary BAM filename in /tmp
    std::string template_path = "/tmp/" + basename + "_XXXXXX.bam";
    std::string tmp_template = template_path;
    char tmpname[tmp_template.size() + 1];
    std::strcpy(tmpname, tmp_template.c_str());

    // Remove ".bam" before calling mkstemp
    tmpname[strlen(tmpname) - 4] = '\0';

    int fd = mkstemp(tmpname);  // mkstemp creates and opens the file
    if (fd == -1) {
        perror("mkstemp");
        return "";
    }
    close(fd); // We'll reopen with sam_open

    std::string bam_path = std::string(tmpname) + ".bam";

    // 2. Create header
    sam_hdr_t* hdr = sam_hdr_init();
    sam_hdr_add_line(hdr, "HD", "VN", "1.6", NULL);
    sam_hdr_add_line(hdr, "SQ", "SN", "1", "LN", "1000000", NULL);

    // 3. Open BAM for writing
    htsFile* out = sam_open(bam_path.c_str(), "wb");
    if (!out) {
        fprintf(stderr, "Error: Cannot open output file %s\n", bam_path.c_str());
        return "";
    }
    if (sam_hdr_write(out, hdr) < 0) {
        fprintf(stderr, "Error: Cannot write BAM header\n");
        sam_close(out);
        return "";
    }

    // 4. Create dummy read
    bam1_t* read = bam_init1();
    read->core.tid = 0;
    read->core.pos = 123456 - 1;
    read->core.qual = 60;
    read->core.flag = 0;
    read->core.l_qname = 7;
    read->core.n_cigar = 1;
    read->core.l_qseq = 10;
    read->core.mtid = -1;
    read->core.mpos = -1;
    read->core.isize = 0;

    int data_len = 7 + 4 + (10 + 1)/2 + 10;
    uint8_t* data = (uint8_t*)malloc(data_len);
    memcpy(data, "read001", 7);
    ((uint32_t*)(data + 7))[0] = bam_cigar_gen(10, BAM_CMATCH);
    uint8_t* seq = data + 7 + 4;
    seq[0] = (1<<4) | 2;
    seq[1] = (4<<4) | 3;
    seq[2] = (1<<4) | 2;
    seq[3] = (4<<4) | 3;
    seq[4] = (1<<4);
    memcpy(seq + 5, "!!!!!!!!!!", 10);

    read->data = data;
    read->l_data = data_len;

    if (sam_write1(out, hdr, read) < 0) {
        fprintf(stderr, "Error: Cannot write read to BAM\n");
        bam_destroy1(read);
        sam_hdr_destroy(hdr);
        sam_close(out);
        return "";
    }

    // 5. Cleanup
    bam_destroy1(read);
    sam_hdr_destroy(hdr);
    sam_close(out);

    // 6. Create index (.bai)
    if (sam_index_build(bam_path.c_str(), 0) != 0) {
        std::cerr << "Failed to create BAM index\n";
        return "";
    }

    return bam_path;
}