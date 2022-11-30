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

#ifndef _CONDITIONING_SET_H
#define _CONDITIONING_SET_H

#include <utils/otools.h>
#include <containers/variant_map.h>
#include <containers/haplotype_set.h>

//Types included in vector < char > var_type
#define TYPE_COMMON	0
#define TYPE_RARE	1
#define TYPE_MONO	2

class conditioning_set {
public:
	//STATIC DATA
	const variant_map & mapG;
	const haplotype_set & H;

	//COUNTS
	unsigned int n_ref_haps;
	unsigned int n_eff_haps;
	unsigned int n_com_sites;
	unsigned int n_tot_sites;
	unsigned int n_states;

	//CONST
	const double nrho;
	const double one_l;

	//VARIANT TYPES
	std::vector < unsigned char > var_type;				//Type of variant. 3 options: (1) common / (2) rare / (3) monomorphic
	std::vector < bool > major_alleles;					//Major alleles at all variants / used for direct imputation
	std::vector < unsigned int > polymorphic_sites;		//List of common sites [TYPE_COMMON or TYPE_RARE]
	std::vector < unsigned int > monomorphic_sites;		//List of monomorphic sites [TYPE_MONO]
	std::vector < bool > lq_flag;

	//CONDITIONING STATES
	std::vector < unsigned int > idxHaps_ref;				//Indexes of conditioning_states in ref haplotype_set
	std::vector < std::vector < unsigned int > > Svar;		//Sparse bitmatrix / Variant first
	bitmatrix Hvar;									//Plain bitmatrix / Variant first

	//TRANSITION & EMISSION PROBABILITIES
	std::vector < float > t;
	std::vector < float > nt;
	const float ed_phs;
	const float ee_phs;
	const float ed_imp;
	const float ee_imp;

	int Kinit;
	int Kpbwt;

	std::vector<unsigned int> swap_ref;
	std::vector<unsigned int> swap_tar;

	//CONSTRUCTOR/DESTRUCTOR/INITIALIZATION
	conditioning_set(const variant_map &, const haplotype_set &, const unsigned int, const unsigned int, const int,const int, const float, const float, const bool );
	~conditioning_set();
	void clear();

	//SELECTION ROUTINES
	void compactSelection(const int ind, const int iter);
	void select(const int ind, const int iter);

	//UPDATE TRANSITION PROBS
	void updateTransitions();

	inline
	float getTransition(const int prev_abs_idx, const int next_abs_idx)  const
	{
		return std::clamp(-expm1(nrho * (mapG.vec_pos[next_abs_idx]->cm - mapG.vec_pos[prev_abs_idx]->cm)), 1e-7, one_l);
	}
};

#endif
