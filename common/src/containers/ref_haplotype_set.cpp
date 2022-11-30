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

#include <containers/ref_haplotype_set.h>
#include <math.h>

ref_haplotype_set::ref_haplotype_set() : n_tot_sites(0),n_rar_sites(0),n_com_sites(0),n_com_sites_hq(0),n_ref_haps(0),sparse_maf(0.001f)
{
}

ref_haplotype_set::~ref_haplotype_set()
{
	n_tot_sites = 0;
	n_rar_sites = 0;
	n_com_sites = 0;
	n_ref_haps = 0;
}

void ref_haplotype_set::allocate()
{
	HvarRef.allocate(n_com_sites, n_ref_haps);
	ShapRef = std::vector < std::vector < int > > (n_ref_haps);
}

void ref_haplotype_set::build_sparsePBWT(const variant_map & M)
{
	tac.clock();

	if (SvarRef.empty())
	{
		SvarRef= std::vector<std::vector<int>>(n_tot_sites);
		for (int h = 0 ; h < n_ref_haps ; h ++)
			for (int r = 0 ; r < ShapRef[h].size() ; r ++)
				SvarRef[ShapRef[h][r]].push_back(h);
	}

	std::vector<bool> ref_small_hap(n_ref_haps, false);
	std::vector<int> map_big_small(n_ref_haps,-1);
	std::vector<int> pbwt_ref_idx(n_ref_haps);

	pbwt_array_A=std::vector<int>(n_ref_haps);
	pbwt_array_B=std::vector<int>(n_ref_haps);

	std::iota(pbwt_array_A.begin(),pbwt_array_A.end(),0);
	std::iota(pbwt_ref_idx.begin(),pbwt_ref_idx.end(),0);

	Ypacked.reserve(std::rintf(sqrtf(n_ref_haps))*n_tot_sites*2);//we allocate more..
	pack3init();
	A_small_idx = std::vector<std::vector <int>>(n_com_sites_hq+1);

	int l_hq=0;
	int last_k=-1; int l_all=0;
	for (int k = 0 ; k < n_tot_sites ; k ++)
	{
		if (M.vec_pos[k]->LQ)
		{
			l_all+=flag_common[k];
			continue;
		}

		if (flag_common[k])
		{
			if (pbwt_small_A.size() > 0 && !(last_k>=0 && flag_common[last_k])) build_init_common(l_hq);
			const int ref_rac_l = M.vec_pos[k]->cref;
			update_full_pbwt_ay(ref_rac_l, l_all,pbwt_ref_idx);
			++l_hq;++l_all;
		}
		else
		{
			if (last_k<0 || flag_common[last_k]) init_small_rare(M,k,l_hq,pbwt_ref_idx,ref_small_hap, map_big_small);
			std::fill(ref_small_hap.begin(), ref_small_hap.end(), major_alleles[k]);
			for (int j=0; j<SvarRef[k].size(); ++j) ref_small_hap[map_big_small[SvarRef[k][j]]]=!major_alleles[k];
			int ref_rac_l = major_alleles[k] ? n_ref_haps-M.vec_pos[k]->calt : pbwt_small_A.size()-M.vec_pos[k]->calt;
			update_small_pbwt_ay(ref_rac_l, ref_small_hap);
		}
		last_k = k;
	}
	Ypacked.shrink_to_fit();
	vrb.bullet("Size sparsePBWT (" + stb.str(Ypacked.size()/(1024*1024)) + " Mb)");
	vrb.bullet("sparsePBWT built and compressed (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void ref_haplotype_set::update_full_pbwt_ay(const int ref_rac_l, const int l, std::vector<int>& pbwt_ref_idx)
{
	unsigned char z;
	int u = 0;
	int v = ref_rac_l;
	int m = 0;
	const std::array<int*,2> occ = {&u,&v};
	int m0;

	while (m < n_ref_haps)
	{
		z = HvarRef.get(l,pbwt_array_A[m]);
		m0 = m;
		while (++m < n_ref_haps && HvarRef.get(l,pbwt_array_A[m]) == z);
		std::copy(pbwt_array_A.begin()+m0, pbwt_array_A.begin()+m, pbwt_array_B.begin() + *occ[z]);
		pack3Add(z, m-m0, Ypacked);
		*occ[z] += m-m0;
	}
	pbwt_array_B.swap(pbwt_array_A);

	for (int h = 0 ; h < n_ref_haps ; h ++) pbwt_ref_idx[pbwt_array_A[h]] = h;
}

void ref_haplotype_set::update_small_pbwt_ay(const int ref_rac_l, std::vector<bool> rare_small_haps)
{
	unsigned char z;
	int u = 0;
	int v = ref_rac_l;
	int m = 0;
	const std::array<int*,2> occ = {&u,&v};
	int m0;
	const int size_a = pbwt_small_A.size();

	while (m < size_a)
	{
		z = rare_small_haps[pbwt_small_A[m]];
		m0 = m;
		while (++m < size_a && rare_small_haps[pbwt_small_A[m]] == z);
		std::copy(pbwt_small_A.begin()+m0, pbwt_small_A.begin()+m, pbwt_small_B.begin() + *occ[z]);
		pack3Add(z, m-m0, Ypacked);
		*occ[z] += m-m0;
	}
	pbwt_small_B.swap(pbwt_small_A);
}

void ref_haplotype_set::build_init_common(const int l)
{
	std::vector<bool> set_big_small(n_ref_haps, true);
	const std::vector<int>& small_idx = A_small_idx[l];
	for (int htr=0; htr<A_small_idx[l].size(); ++htr) set_big_small[small_idx[htr]] = false;

	int n_zeros=0;
	for (int htr=0; htr<n_ref_haps; ++htr)
		if (set_big_small[htr]) pbwt_array_B[n_zeros++] = pbwt_array_A[htr];

	for (int htr=0; htr<pbwt_small_A.size(); ++htr,++n_zeros)
		pbwt_array_B[n_zeros]=pbwt_array_A[small_idx[pbwt_small_A[htr]]];

	pbwt_array_B.swap(pbwt_array_A);
}

void ref_haplotype_set::init_small_rare(const variant_map & M,const int k, const int l, const std::vector<int>& pbwt_ref_idx, std::vector<bool>& ref_small_hap, std::vector<int>& map_big_small)
{
	std::vector<int> rare_ref_haps;
	std::map<int,int> pbwt_small_ref_idx_map;
	if (l>0) for (int htr=0; htr<A_small_idx[l-1].size(); ++htr) map_big_small[A_small_idx[l-1][htr]] = -1; //reset map

	for (int kk=k; kk < n_tot_sites; ++kk)
	{
		if (M.vec_pos[kk]->LQ) continue;
		if (flag_common[kk]) break;
		rare_ref_haps.insert(rare_ref_haps.end(), SvarRef[kk].begin(), SvarRef[kk].end());
	}
	sort(rare_ref_haps.begin(), rare_ref_haps.end());
	rare_ref_haps.erase(unique(rare_ref_haps.begin(), rare_ref_haps.end()), rare_ref_haps.end());

	for (int h = 0 ; h < rare_ref_haps.size() ; h ++) pbwt_small_ref_idx_map.insert(std::pair<int,int>(pbwt_ref_idx[rare_ref_haps[h]],rare_ref_haps[h]));
	int kk2=0;
	A_small_idx[l] = std::vector<int>(pbwt_small_ref_idx_map.size());
	for (auto iter = pbwt_small_ref_idx_map.begin() ; iter != pbwt_small_ref_idx_map.end() ; iter++, kk2++)
	{
		map_big_small[iter->second] = kk2;
		A_small_idx[l][kk2] = iter->first;
	}
	pbwt_small_A.resize(A_small_idx[l].size());
	std::iota(pbwt_small_A.begin(),pbwt_small_A.end(),0);
	pbwt_small_B.resize(A_small_idx[l].size());
	ref_small_hap.resize(A_small_idx[l].size());
}
