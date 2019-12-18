#include "genotype_stream.h"

const float unphred[256] = { 1.000000e+00, 7.943282e-01, 6.309573e-01, 5.011872e-01, 3.981072e-01, 3.162278e-01, 2.511886e-01, 1.995262e-01, 1.584893e-01, 1.258925e-01, 1.000000e-01, 7.943282e-02, 6.309573e-02, 5.011872e-02, 3.981072e-02, 3.162278e-02, 2.511886e-02, 1.995262e-02, 1.584893e-02, 1.258925e-02, 1.000000e-02, 7.943282e-03, 6.309573e-03, 5.011872e-03, 3.981072e-03, 3.162278e-03, 2.511886e-03, 1.995262e-03, 1.584893e-03, 1.258925e-03, 1.000000e-03, 7.943282e-04, 6.309573e-04, 5.011872e-04, 3.981072e-04, 3.162278e-04, 2.511886e-04, 1.995262e-04, 1.584893e-04, 1.258925e-04, 1.000000e-04, 7.943282e-05, 6.309573e-05, 5.011872e-05, 3.981072e-05, 3.162278e-05, 2.511886e-05, 1.995262e-05, 1.584893e-05, 1.258925e-05, 1.000000e-05, 7.943282e-06, 6.309573e-06, 5.011872e-06, 3.981072e-06, 3.162278e-06, 2.511886e-06, 1.995262e-06, 1.584893e-06, 1.258925e-06, 1.000000e-06, 7.943282e-07, 6.309573e-07, 5.011872e-07, 3.981072e-07, 3.162278e-07, 2.511886e-07, 1.995262e-07, 1.584893e-07, 1.258925e-07, 1.000000e-07, 7.943282e-08, 6.309573e-08, 5.011872e-08, 3.981072e-08, 3.162278e-08, 2.511886e-08, 1.995262e-08, 1.584893e-08, 1.258925e-08, 1.000000e-08, 7.943282e-09, 6.309573e-09, 5.011872e-09, 3.981072e-09, 3.162278e-09, 2.511886e-09, 1.995262e-09, 1.584893e-09, 1.258925e-09, 1.000000e-09, 7.943282e-10, 6.309573e-10, 5.011872e-10, 3.981072e-10, 3.162278e-10, 2.511886e-10, 1.995262e-10, 1.584893e-10, 1.258925e-10, 1.000000e-10, 7.943282e-11, 6.309573e-11, 5.011872e-11, 3.981072e-11, 3.162278e-11, 2.511886e-11, 1.995262e-11, 1.584893e-11, 1.258925e-11, 1.000000e-11, 7.943282e-12, 6.309573e-12, 5.011872e-12, 3.981072e-12, 3.162278e-12, 2.511886e-12, 1.995262e-12, 1.584893e-12, 1.258925e-12, 1.000000e-12, 7.943282e-13, 6.309573e-13, 5.011872e-13, 3.981072e-13, 3.162278e-13, 2.511886e-13, 1.995262e-13, 1.584893e-13, 1.258925e-13, 1.000000e-13, 7.943282e-14, 6.309573e-14, 5.011872e-14, 3.981072e-14, 3.162278e-14, 2.511886e-14, 1.995262e-14, 1.584893e-14, 1.258925e-14, 1.000000e-14, 7.943282e-15, 6.309573e-15, 5.011872e-15, 3.981072e-15, 3.162278e-15, 2.511886e-15, 1.995262e-15, 1.584893e-15, 1.258925e-15, 1.000000e-15, 7.943282e-16, 6.309573e-16, 5.011872e-16, 3.981072e-16, 3.162278e-16, 2.511886e-16, 1.995262e-16, 1.584893e-16, 1.258925e-16, 1.000000e-16, 7.943282e-17, 6.309573e-17, 5.011872e-17, 3.981072e-17, 3.162278e-17, 2.511886e-17, 1.995262e-17, 1.584893e-17, 1.258925e-17, 1.000000e-17, 7.943282e-18, 6.309573e-18, 5.011872e-18, 3.981072e-18, 3.162278e-18, 2.511886e-18, 1.995262e-18, 1.584893e-18, 1.258925e-18, 1.000000e-18, 7.943282e-19, 6.309573e-19, 5.011872e-19, 3.981072e-19, 3.162278e-19, 2.511886e-19, 1.995262e-19, 1.584893e-19, 1.258925e-19, 1.000000e-19, 7.943282e-20, 6.309573e-20, 5.011872e-20, 3.981072e-20, 3.162278e-20, 2.511886e-20, 1.995262e-20, 1.584893e-20, 1.258925e-20, 1.000000e-20, 7.943282e-21, 6.309573e-21, 5.011872e-21, 3.981072e-21, 3.162278e-21, 2.511886e-21, 1.995262e-21, 1.584893e-21, 1.258925e-21, 1.000000e-21, 7.943282e-22, 6.309573e-22, 5.011872e-22, 3.981072e-22, 3.162278e-22, 2.511886e-22, 1.995262e-22, 1.584893e-22, 1.258925e-22, 1.000000e-22, 7.943282e-23, 6.309573e-23, 5.011872e-23, 3.981072e-23, 3.162278e-23, 2.511886e-23, 1.995262e-23, 1.584893e-23, 1.258925e-23, 1.000000e-23, 7.943282e-24, 6.309573e-24, 5.011872e-24, 3.981072e-24, 3.162278e-24, 2.511886e-24, 1.995262e-24, 1.584893e-24, 1.258925e-24, 1.000000e-24, 7.943282e-25, 6.309573e-25, 5.011872e-25, 3.981072e-25, 3.162278e-25, 2.511886e-25, 1.995262e-25, 1.584893e-25, 1.258925e-25, 1.000000e-25, 7.943282e-26, 6.309573e-26, 5.011872e-26, 3.981072e-26, 3.162278e-26};

genotype_stream::genotype_stream(std::string _region, float _maf_common, unsigned long _n_variants,unsigned long _n_main_samples, unsigned long _n_ref_samples) :
	region(_region),
	maf_common(_maf_common),
	n_variants(_n_variants),
	n_main_samples(_n_main_samples),
	n_ref_samples(_n_ref_samples)
{
	is_open = false;
	sr = nullptr;
	i_variant = 0;
	i_common = -1;
	is_common_variant = false;
	nset = 0;
	ngl_main = 0;
	ngl_arr_main = 0;
	gl_arr_main = nullptr;
	ngt_ref = 0;
	gt_arr_ref = nullptr;
	ngt_arr_ref = 0;
	af_ptr=nullptr;///read MAF
	maf_ref_variant = 0.0f;
	nval = 0;
	line_main = nullptr;
	line_ref = nullptr;

}

genotype_stream::~genotype_stream()
{
	if (is_open) closeStream();
}

void genotype_stream::openStream(std::string funphased, std::string freference) {
	sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	bcf_sr_set_regions(sr, region.c_str(), 0);
	bcf_sr_add_reader (sr, funphased.c_str());
	bcf_sr_add_reader (sr, freference.c_str());

	ref_haps = std::vector<bool> (2*n_ref_samples);
	curr_GL = std::vector<double> (3*n_main_samples,1.0);

	is_open = true;
}

bool genotype_stream::readMarker()
{
	if (!is_open) return false;
	is_common_variant = false;
	double lg0, lg1, lg2, sum;
	while ((nset = bcf_sr_next_line (sr))) {
		if (nset == 2) {
			line_main =  bcf_sr_get_line(sr, 0);
			line_ref =  bcf_sr_get_line(sr, 1);
			if (line_main->n_allele == 2 && line_ref->n_allele == 2) {
				bcf_get_info_float(sr->readers[1].header, line_ref, "AF", &af_ptr, &nval);
				maf_ref_variant = std::min(af_ptr[0],1.0f -af_ptr[0]);
				if (maf_common <= maf_ref_variant)
				{
					i_common ++;
					is_common_variant = true;
					//curr_variant = *V.vec_pos[i_common];
					i_variant ++;
					break;
				}

				bcf_unpack(line_main, BCF_UN_STR);
				std::string chr = bcf_hdr_id2name(sr->readers[0].header, line_main->rid);
				int pos = line_main->pos + 1;
				std::string id = std::string(line_main->d.id);
				std::string ref = std::string(line_main->d.allele[0]);
				std::string alt = std::string(line_main->d.allele[1]);
				curr_variant = variant (chr, pos, id, ref, alt, 0);
				unsigned int cref = 0, calt = 0, cmis = 0;

				ngt_ref = bcf_get_genotypes(sr->readers[1].header, line_ref, &gt_arr_ref, &ngt_arr_ref);
				assert(ngt_ref == 2 * n_ref_samples);
				for(int i = 0 ; i < 2 * n_ref_samples ; i += 2) {
					bool a0 = (bcf_gt_allele(gt_arr_ref[i+0])==1);
					bool a1 = (bcf_gt_allele(gt_arr_ref[i+1])==1);
					ref_haps[i+0] = a0;
					ref_haps[i+1] = a1;
					a0?calt++:cref++;
					a1?calt++:cref++;
				}

				ngl_main = bcf_get_format_int32(sr->readers[0].header, line_main, "PL", &gl_arr_main, &ngl_arr_main);
				if (ngl_main == 3 * n_main_samples) {
					for(int i = 0 ; i < 3 * n_main_samples ; i += 3) {
						if (gl_arr_main[i+0] != bcf_float_missing && gl_arr_main[i+1] != bcf_float_missing && gl_arr_main[i+2] != bcf_float_missing) {
							lg0 = unphred[(unsigned char)gl_arr_main[i+0]];
							lg1 = unphred[(unsigned char)gl_arr_main[i+1]];
							lg2 = unphred[(unsigned char)gl_arr_main[i+2]];
							sum = lg0+lg1+lg2;
							if (sum > 0.0f)
							{
								curr_GL[i+0] = lg0 / sum;
								curr_GL[i+1] = lg1 / sum;
								curr_GL[i+2] = lg2 / sum;
							}
							else std::fill(std::next(curr_GL.begin(),i+0), std::next(curr_GL.begin(),i+3), 1.0f);
						}
						else std::fill(std::next(curr_GL.begin(),i+0), std::next(curr_GL.begin(),i+3), 1.0f);
					}
				}
				else
				{
					std::fill(curr_GL.begin(), curr_GL.end(), 1.0f);
				}

				curr_variant.cref = cref;curr_variant.calt = calt;curr_variant.cmis = cmis;
				i_variant ++;
			}
			break;
		}
	}
	//if (!nset) closeStream();
	return nset;
}

void genotype_stream::closeStream()
{
	free(gl_arr_main);
	free(gt_arr_ref);
	bcf_sr_destroy(sr);
	is_open=false;
}

