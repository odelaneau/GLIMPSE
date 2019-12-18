#ifndef _GENOTYPE_STREAM_H_
#define _GENOTYPE_STREAM_H_

#include <utils/otools.h>

#include <containers/variant_map.h>
#include <containers/haplotype_set.h>
#include <objects/variant.h>

class genotype_stream {
public:
	bool is_open;
	//const variant_map V;
	const std::string region;
	const float maf_common;


	//COUNTS
	const unsigned long n_variants;
	const unsigned long n_main_samples;
	const unsigned long n_ref_samples;

	// streaming obj
	variant curr_variant;
	std::vector<bool> ref_haps;
	std::vector<double> curr_GL;

	//reader
	bcf_srs_t * sr;
	unsigned int i_variant;
	int i_common;
	bool is_common_variant;
	unsigned int nset;
	int ngl_main;
	int ngl_arr_main;
	int *gl_arr_main;
	int ngt_ref;
	int *gt_arr_ref;
	int ngt_arr_ref;
	float *af_ptr;///read MAF
	float maf_ref_variant;
	int nval;
	bcf1_t * line_main;
	bcf1_t * line_ref;

	genotype_stream(std::string _region, float _maf_common, unsigned long _n_variants,unsigned long _n_main_samples, unsigned long _n_ref_samples);
	virtual ~genotype_stream();

	void openStream(std::string funphased, std::string freference);

	bool readMarker();

	void closeStream();

};

#endif /* _GENOTYPE_STREAM_H_ */
