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

#include <containers/conditioning_set.h>
#include <algorithm>

conditioning_set::conditioning_set(const argument_set& _A, const variant_map & _mapG,const haplotype_set & _H):
		A(_A), mapG(_mapG), H(_H), n_ref_haps(H.n_ref_haps),
		n_eff_haps(A.mNe),
		n_com_sites(H.common2tot.size()),
		n_tot_sites(mapG.size()),
		nrho(!A.mStateListFilename.empty() ? (-0.04 * n_eff_haps) / std::max(n_ref_haps,n_eff_haps) : (-0.04f * n_eff_haps) / n_ref_haps),
		one_l(1.0 - 1e-7),
		ed_phs(A.mErrPhase),
		ee_phs(1.0 -A.mErrPhase),
		ed_imp(A.mErrImp),
		ee_imp(1.0 - A.mErrImp)
{
	var_type = std::vector < unsigned char > (n_tot_sites);
	major_alleles = H.major_alleles;

	idxHaps_ref.clear();

	n_states = 0;
	Svar = std::vector < std::vector < unsigned int > > (n_tot_sites);

	//Build flags of HQ or LQ sites
	lq_flag = std::vector < bool > (n_tot_sites, true);
	for (int l = 0 ; l < mapG.vec_pos.size() ; l ++) lq_flag[l] = mapG.vec_pos[l]->LQ;
}

conditioning_set::~conditioning_set() {
	n_eff_haps=0;
	n_tot_sites=0;
	var_type.clear();
	idxHaps_ref.clear();
	polymorphic_sites.clear();
	monomorphic_sites.clear();
	major_alleles.clear();
	n_states=0;
	t.clear();
	nt.clear();
}

void conditioning_set::select(const int ind, const int iter) {
	compactSelection(ind, iter);
	updateTransitions();
}

void conditioning_set::compactSelection(const int ind, const int iter)
{
	bool use_list = true;
	const int hapid=H.tar_ind2hapid[ind];
	const int ploidyM1 = H.tar_ploidy[ind]-1;
	idxHaps_ref.clear();

	if (iter==STAGE_INIT && A.Kinit > 0)
	{
		idxHaps_ref.resize(H.init_states[ind].size());
		std::copy(H.init_states[ind].begin(), H.init_states[ind].end(), idxHaps_ref.begin());
	}
	else if (A.Kpbwt > 0 && A.Kpbwt < H.n_ref_haps)
	{
		std::unordered_set<int> states_ind;
		for (int i=0; i<A.mMinPbwtDepth; ++i)
		{
			for (int j=0; j<H.pbwt_states[ind][i].size(); ++j)
			{
				states_ind.insert(H.pbwt_states[ind][i][j]);
			}
		}
		for (int i=A.mMinPbwtDepth; i<H.pbwt_states[ind].size() && states_ind.size() < A.Kpbwt; ++i)
		{
			for (int j=0; j<H.pbwt_states[ind][i].size() && states_ind.size() < A.Kpbwt; ++j)
			{
				states_ind.insert(H.pbwt_states[ind][i][j]);
			}
		}
		idxHaps_ref.resize(states_ind.size());
		std::copy(states_ind.begin(), states_ind.end(), idxHaps_ref.begin());
		std::sort(idxHaps_ref.begin(),idxHaps_ref.end());
	}
	else if (A.Kpbwt >= H.n_ref_haps)
	{
		idxHaps_ref.resize(H.n_ref_haps);
		std::iota(idxHaps_ref.begin(), idxHaps_ref.end(), 0);
		use_list = false;
	}
	//A.Kpbwt == 0 easy: just go here..
	if (use_list && (H.list_states[hapid].size() > 0 || H.list_states[hapid+ploidyM1].size() > 0))
	{
		std::copy(H.list_states[hapid].begin(), H.list_states[hapid].end(), std::back_inserter(idxHaps_ref));
		if (ploidyM1) std::copy(H.list_states[hapid+1].begin(), H.list_states[hapid+1].end(), std::back_inserter(idxHaps_ref));
		sort(idxHaps_ref.begin(), idxHaps_ref.end());
		idxHaps_ref.erase(unique(idxHaps_ref.begin(), idxHaps_ref.end()), idxHaps_ref.end());
	}

	//Check #states
	n_states = idxHaps_ref.size();
	if (n_states == 0) vrb.error("States for individual " + std::to_string(ind) + " are zero. Error during selection.");

	//Update Svar by looking at sparse matrixes of conditioning haplotypes [this is quite costly, very inefficient at the cache level]
	for (int l = 0 ; l < n_tot_sites ; l++) Svar[l].clear();
	for (int k = 0 ; k < idxHaps_ref.size() ; k++) {
		for (int r = 0 ; r < H.ShapRef[idxHaps_ref[k]].size() ; r++) {
			Svar[H.ShapRef[idxHaps_ref[k]][r]].push_back(k);
		}
	}
	//Update var_type
	monomorphic_sites.clear();
	polymorphic_sites.clear();
	//std::vector<unsigned int> col2row_offsets(n_com_sites);
	for (unsigned int l = 0, lcom=0,lrel=0 ; l < n_tot_sites ; l++) {
		if (H.flag_common[l])
		{
			var_type[l] = TYPE_COMMON;
			//col2row_offsets[lcom++] = lrel;
		}
		else if (Svar[l].size())
		{
			var_type[l] = TYPE_RARE;
		}
		else
		{
			var_type[l] = TYPE_MONO;
			monomorphic_sites.push_back(l);
			continue;
		}
		polymorphic_sites.push_back(l);
		lrel++;
	}
	//vrb.print(stb.str(polymorphic_sites.size()));
	//tac.clock();

	/*
	Hvar.reallocate(polymorphic_sites.size(), n_states);
	Hhapsubset.reallocate(n_states,n_com_sites);
	Hvarsubset.reallocate(n_com_sites,n_states);
	Hhapsubset.subset(H.HhapRef, idxHaps_ref);
	Hhapsubset.transpose(Hvarsubset,n_states,n_com_sites);
	for (int labs = 0, lrel = 0; labs < n_tot_sites ; labs ++)
	{
		if (var_type[labs] == TYPE_MONO) continue;
		else if (var_type[labs] == TYPE_RARE)
		{
			Hvar.set(lrel, major_alleles[labs]);
			for (int r = 0 ; r < Svar[labs].size() ; r++) Hvar.set(lrel, Svar[labs][r], !major_alleles[labs]);
		} //else mono: do nothing
		++lrel;
	}
	*/
	Hvar.reallocate(polymorphic_sites.size(), n_states);
	Hhapsubset.reallocate(n_states,n_com_sites);
	Hvarsubset.reallocate(n_com_sites,n_states);
	Hhapsubset.subset(H.HhapRef, idxHaps_ref);
	Hhapsubset.transpose(Hvarsubset,n_states,n_com_sites);
	for (int labs = 0, lrel = 0, lcom = 0 ; labs < n_tot_sites ; labs ++) {
		if (var_type[labs] == TYPE_COMMON) {
			std::memcpy(Hvar.getRowPtr(lrel), Hvarsubset.getRowPtr(lcom), (Hvarsubset.n_cols>>3));
			//for (int k = 0 ; k < idxHaps_ref.size() ; k++) Hvar1.set(lrel, k, Hvarsubset.get(lcom, k));
			lrel++;
			lcom++;
		} else if (var_type[labs] == TYPE_RARE) {
			Hvar.set(lrel, major_alleles[labs]);
			for (int r = 0 ; r < Svar[labs].size() ; r++) Hvar.set(lrel, Svar[labs][r], !major_alleles[labs]);
			lrel++;
		} //else mono: do nothing
	}
/*
*/
	//vrb.bullet("Hapvar2 (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
/*
	//Build bitmatrix Hvar
	bitmatrix Hvar1;
	tac.clock();

	Hvar.reallocate(polymorphic_sites.size(), n_states);
	for (int labs = 0, lrel = 0, lcom = 0 ; labs < n_tot_sites ; labs ++) {
		if (var_type[labs] == TYPE_COMMON) {
			for (int k = 0 ; k < idxHaps_ref.size() ; k++) Hvar.set(lrel, k, H.HvarRef.get(lcom, idxHaps_ref[k]));
			lrel++;
			lcom++;
		} else if (var_type[labs] == TYPE_RARE) {
			Hvar.set(lrel, major_alleles[labs]);
			for (int r = 0 ; r < Svar[labs].size() ; r++) Hvar.set(lrel, Svar[labs][r], !major_alleles[labs]);
			lrel++;
		} //else mono: do nothing
	}

	vrb.bullet("Hapvar0 (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");

	{
		//Hvar1.allocate(polymorphic_sites.size(), n_states);
		tac.clock();
		Hvar1.reallocate(polymorphic_sites.size(), n_states);
		Hhapsubset.reallocate(n_states,n_com_sites);
		Hvarsubset.reallocate(n_com_sites,n_states);
		Hhapsubset.subset(H.HhapRef, idxHaps_ref);
		Hhapsubset.transpose(Hvarsubset,n_states,n_com_sites);
		for (int labs = 0, lrel = 0, lcom = 0 ; labs < n_tot_sites ; labs ++) {
			if (var_type[labs] == TYPE_COMMON) {
				std::memcpy(Hvar1.getRowPtr(lrel), Hvarsubset.getRowPtr(lcom), (Hvarsubset.n_cols>>3));
				//for (int k = 0 ; k < idxHaps_ref.size() ; k++) Hvar1.set(lrel, k, Hvarsubset.get(lcom, k));
				lrel++;
				lcom++;
			} else if (var_type[labs] == TYPE_RARE) {
				Hvar1.set(lrel, major_alleles[labs]);
				for (int r = 0 ; r < Svar[labs].size() ; r++) Hvar1.set(lrel, Svar[labs][r], !major_alleles[labs]);
				lrel++;
			} //else mono: do nothing
		}
		vrb.bullet("Hapvar1 (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	}
*/
/*
	{

		Hvar2.allocate(polymorphic_sites.size(), n_states);
		tac.clock();
		Hhapsubset.reallocate(n_states,n_com_sites);
		Hhapsubset.subset(H.HhapRef, idxHaps_ref);
		Hhapsubset.transpose(Hvar,n_states,n_com_sites,col2row_offsets);
		for (int labs = 0, lrel = 0; labs < n_tot_sites ; labs ++)
		{
			if (var_type[labs] == TYPE_MONO) continue;
			else if (var_type[labs] == TYPE_RARE)
			{
				Hvar.set(lrel, major_alleles[labs]);
				for (int r = 0 ; r < Svar[labs].size() ; r++) Hvar.set(lrel, Svar[labs][r], !major_alleles[labs]);
			} //else mono: do nothing
			++lrel;
		}
		//vrb.bullet("Hapvar2 (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	}
*/
/*
	for (int labs = 0, lrel = 0, lcom = 0 ; labs < n_tot_sites ; labs ++)
	{
		if (var_type[labs] == TYPE_MONO) continue;
		for (int r = 0 ; r <n_states ; r++)
			if (Hvar1.get(lrel, r)!=Hvar.get(lrel,r))
			{
				vrb.error(stb.str(lrel) + " - " + stb.str(r));
			}
		++lrel;
	}
	vrb.print("OK");
*/
}

void conditioning_set::updateTransitions()
{
	if (polymorphic_sites.size() == 0) return;
	t.resize(polymorphic_sites.size() - 1);
	nt.resize(polymorphic_sites.size() - 1);
	for (int l = 1 ; l < polymorphic_sites.size() ; l ++)
	{
		t[l-1] = std::clamp(-expm1(nrho * (mapG.vec_pos[polymorphic_sites[l]]->cm - mapG.vec_pos[polymorphic_sites[l-1]]->cm)), 1e-7, one_l);
		nt[l-1] = 1.0f-t[l-1];
	}
}
