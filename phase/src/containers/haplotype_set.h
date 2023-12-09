/*******************************************************************************
 * Copyright (C) 2022-2023 Simone Rubinacci
 * Copyright (C) 2022-2023 Olivier Delaneau
 *
 * MIT Licence
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#ifndef _HAPLOTYPE_SET_H
#define _HAPLOTYPE_SET_H

#include <utils/otools.h>
#include <utils/checksum_utils.h>

#include <containers/bitmatrix.h>
#include <containers/genotype_set.h>
#include <containers/ref_haplotype_set.h>

class haplotype_set : public ref_haplotype_set {
public:
	//COUNTS
	unsigned int n_tot_haps;					// #haplotypes [target + reference samples]
	unsigned int n_tar_haps;					// #haplotypes [target samples]
	unsigned int n_tar_samples;

	//HAPLOTYPE DATA [plain/sparse bitmatrix representations]
	bitmatrix HvarTar;

	std::vector < std::vector < int > > ShapTar;				// Rare alleles per haplotype
	std::vector < std::vector < int > > SvarTar;				// Rare alleles per variant
	std::vector < std::vector < int > > SindTarGL;				// Rare alleles per ind from GLs
	std::vector < float > cm_pos;

	//PLOIDY
	int fploidy;								//Format ploidy field, to indicate the ploidy in the sample file: 1=only haploids, 2=only diploids, -2=mixed ploidy (haploids and diploids).
	int max_ploidy;
	std::vector < int > tar_ploidy;
	std::vector < int > tar_ind2hapid;
	std::vector < int > tar_hapid2ind;

	//PBWT
	int pbwt_depth;
	float pbwt_modulo_cm;
	std::vector < int > pbwt_array_V;

	//FM
	std::vector < int > pbwt_index;
	std::vector < int > pbwt_small_index;

	std::vector<int> pbwt_small_V;
	std::vector<int> rareTarHaps;

	std::vector < int > f_k;
	std::vector < int > g_k;
	std::vector<int> f_k_small;
	std::vector<int> g_k_small;
	std::vector<int> last_reset;
	std::vector<int> last_rare;


	//Matchings
	std::vector < bool > pbwt_stored;
	std::vector < int > pbwt_grp;

	int Kinit;
	int Kpbwt;
	int K;
	int nstored;

	int counter_gf;
	int counter_sel_gf;
	int counter_rare_restarts;

	std::vector<std::vector<std::vector<int>>> pbwt_states;
	std::vector<std::set<int>> init_states;
	std::vector<std::vector<int>> list_states;

	std::vector<unsigned char> tar_hap;
	std::map<int, int> rare_idx_to_id;

/*
	std::vector< ref_haplotype_map > ref_hashmap;
	std::vector<std::priority_queue < composite_hap > > q;
	std::vector<std::vector<std::vector<int>>> composite_h;
	std::vector<std::vector<std::vector<int>>> composite_m;
	float reject_threshold;
	int rejected_states;
	int switch_states;
	stats1D comp_dist_cm;
*/

	//stats
	stats1D length_sel_mod;
	stats1D length_sel_gf;
	stats1D length_sel_gf_rare;



	//CONSTRUCTOR/DESTRUCTOR/INITIALIZATION
	haplotype_set();

	virtual ~haplotype_set();
	void allocate();
	void allocate_hap_only();


	//ROUTINES
	void initializeAndAllocate(unsigned int, unsigned int, unsigned int, unsigned int);
	void initRareTar(const genotype_set & G, const variant_map& M);
	void updateHaplotypes(const genotype_set &);
	void transposeRareRef();
	void transposeRareTar();

	//Init
	void performSelection_RARE_INIT_GL(const variant_map & M);
	void read_list_states(const std::string file_list);


	//PBWT ROUTINES
	void allocatePBWT(const int _pbwt_depth, const float _pbwt_modulo_cm, const variant_map & V, const genotype_set & G, const int _Kinit, const int _Kpbwt);
	void matchHapsFromCompressedPBWTSmall(const variant_map & V, const bool main_iteration);
	void init_common(const int k, const int l, const int ref_rac_l_com);
	void init_rare(const variant_map & M,const int k, const int l);
	void read_full_pbwt_av(const unsigned char*& pY, const int ref_rac_l);
	void read_small_pbwt_av(const unsigned char*& pY, const int ref_rac_l, const bool update_v);
	void selectK(const int htr, const int k, const int ref_rac_l, const std::vector < int >& pbwt_array, const int k0, const unsigned char a);
	void selectKrare(const int htr, const int k, const int ref_rac_l, const std::vector < int >& pbwt_array, const int k0, const unsigned char a);
	void select_common_pd_fg(const int k, const int l_hq, const int l_all, const int ref_rac_l, const int prev_ref_rac_l);
	void select_rare_pd_fg(const int k, const int ref_rac_l);

	void update_checksum(checksum &crc) const
	{
		ref_haplotype_set::update_checksum(crc);
		crc.process_data(n_tot_haps);
		crc.process_data(n_tar_haps);
		crc.process_data(n_tar_samples);
		HvarTar.update_checksum(crc);
		crc.process_data(ShapTar);
		crc.process_data(SvarTar);
		crc.process_data(SindTarGL);
		crc.process_data(cm_pos);
		crc.process_data(fploidy);
		crc.process_data(max_ploidy);
		crc.process_data(tar_ploidy);
		crc.process_data(tar_ind2hapid);
		crc.process_data(tar_hapid2ind);
	}
};

#endif
