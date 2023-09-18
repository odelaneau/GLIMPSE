
#define _DECLARE_TOOLBOX_HERE

#include <iostream>
#include <containers/ref_haplotype_set.h>
#include <boost/archive/binary_iarchive.hpp>


int main(int argc, char ** argv) {
    if (argc != 2) {
        vrb.error("Invalid command line options. Usage: GLIMPSE2_extract_sites_from_chunk {path_to_reference.bin}");
    }
    std::cout<<"test"<<std::endl;
	ref_haplotype_set H;

    const std::string reference_filename = argv[1];
    std::cout <<"ref ="<<reference_filename<<std::endl;
	std::ifstream ifs(reference_filename, std::ios::binary | std::ios_base::in);
    std::cout<<"ifstream"<<std::endl;
    if (!ifs.good()) vrb.error("Reading binary reference panel file: [" + reference_filename + "]. File not good(): eofbit, failbit or badbit set or file not found.");
    try
    {
        std::cout<<"in try"<<std::endl;
        boost::archive::binary_iarchive ia(ifs);
        
        ia >> H;
    } catch (std::exception& e ) {
        vrb.error(e.what());
        vrb.error("problems reading the binary reference panel (exception triggered by boost archive). Please ensure you are using the same GLIMPSE and boost library version");
    }
    if (H.Ypacked.size()==0) vrb.error("Problem reading binary file format. Empty PBWT detected.");

    std::cout << "Lrare=" << H.n_rar_sites << std::endl;
    std::cout << "Lcommon=" << H.n_com_sites << std::endl;

    return EXIT_SUCCESS;
}

