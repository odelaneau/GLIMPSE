#include "bcf_traversal.hpp"

void BcfTraversal::traverse(const std::string filename) {
    initialize_bcf_file_reader(bcf_fri, filename);

    handle_bcf_file_reader();
    while(bcf_next_line(bcf_fri)) {
        // Unpack the line and get genotypes
        bcf_unpack(bcf_fri.line, BCF_UN_STR);
        bcf_fri.ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.size_gt_arr));
        line_max_ploidy = bcf_fri.ngt / bcf_fri.n_samples;

        handle_bcf_line();
    }
    destroy_bcf_file_reader(bcf_fri);
}

void BcfTransformer::transform(const std::string& ifname, const std::string& ofname) {
    fp = hts_open(ofname.c_str(), ofname.compare("-") ? "wb" : "wu"); // "-" for stdout
    if (fp == NULL) {
        std::cerr << "Could not open " << ofname << std::endl;
        throw "File open error";
    }

    traverse(ifname);

    // Close / Release ressources
    hts_close(fp);
    bcf_hdr_destroy(hdr);
}

void BcfTransformer::handle_bcf_file_reader() {
    // Duplicate the header from the input bcf
    hdr = bcf_hdr_dup(bcf_fri.sr->readers[0].header);

    transform_header();

    // Write the header to the new file
    int ret = bcf_hdr_write(fp, hdr);
    if (ret) {
        std::cerr << "Failed to write header" << std::endl;
        exit(-1);
    }
}

void BcfTransformer::handle_bcf_line() {
    rec = bcf_fri.line;

    const size_t PLOIDY = 2;
    // Check ploidy, only support diploid for the moment
    if (line_max_ploidy != PLOIDY) {
        std::cerr << "[ERROR] Ploidy of samples is different than 2" << std::endl;
        exit(-1); // Change this
    }

    transform_record();
    bcf_update_genotypes(hdr, rec, bcf_fri.gt_arr, bcf_hdr_nsamples(hdr) * PLOIDY);

    int ret = bcf_write1(fp, hdr, rec);
    if (ret) {
        std::cerr << "Failed to write record" << std::endl;
        exit(-1);
    }
}

void BcfUnphaser::transform_record() {
    for (size_t i = 0; i < bcf_fri.n_samples; ++i) {
        int32_t allele_0 = bcf_gt_allele(bcf_fri.gt_arr[i*2]);
        int32_t allele_1 = bcf_gt_allele(bcf_fri.gt_arr[i*2+1]);
        if (distrib(gen) == 1) {
            bcf_fri.gt_arr[i*2] = bcf_gt_unphased(allele_0);
            bcf_fri.gt_arr[i*2+1] = bcf_gt_unphased(allele_1);
        } else {
            bcf_fri.gt_arr[i*2] = bcf_gt_unphased(allele_1);
            bcf_fri.gt_arr[i*2+1] = bcf_gt_unphased(allele_0);
        }
    }
}