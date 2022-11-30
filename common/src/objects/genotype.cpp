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

#include <objects/genotype.h>

genotype::genotype(const std::string _name, const int _index, const int _n_variants, const int _ploidy, const int _hapid) :
		name(_name),
		index(_index),
		n_variants(_n_variants),
		ploidy(_ploidy),
		hapid(_hapid),
		stored_cnt(0)
{
}

genotype::~genotype() {
	free();
}

void genotype::free() {
	GL.clear();
	stored_data.clear();
	H0.clear();
	H1.clear();
}

void genotype::allocate() {
	GL = std::vector < unsigned char > ((ploidy+1) * n_variants, (unsigned char) 0);
	flat = std::vector < bool > (n_variants, true);
	H0 = std::vector < bool > (n_variants, false);
	if (ploidy > 1) H1 = std::vector < bool > (n_variants, false);
}

void genotype::initHaplotypeLikelihoods(std::vector < float > & HL, const float min_gl) {
	std::array<float,3> tmp = {0.0f, 0.0f, 0.0f};
	float sum;
	const int ploidyP1 = ploidy+1;

	for (int l = 0 ; l < n_variants ; l ++) {
		if (!flat[l]) {
			tmp[0]=unphred[GL[ploidyP1*l+0]];
			tmp[1]=unphred[GL[ploidyP1*l+1]];
			tmp[2]=ploidy > 1? unphred[GL[ploidyP1*l+2]] : 0.0f;

			sum = tmp[0] + tmp[1] + tmp[2];
			tmp[0] /= sum;
			tmp[1] /= sum;
			tmp[2] /= sum;

			//for ploidy==1 this is it
			HL[2*l+0] = tmp[0];
			HL[2*l+1] = tmp[1];

			//if ploidy ==2 we need to consider tmp[2]
			if (ploidy > 1)
			{
				HL[2*l+0] += 0.5 * tmp[1];
				HL[2*l+1] = 0.5 * tmp[1] + tmp[2];
			}
		}
		else
		{
			HL[2*l+0] = .5;
			HL[2*l+1] = .5;
		}
		if (HL[2*l+0] < min_gl) {HL[2*l+0]=min_gl; HL[2*l+1]=1.0f-min_gl;}
		if (HL[2*l+1] < min_gl) {HL[2*l+0]=1.0f-min_gl; HL[2*l+1]=min_gl;}
	}
}

//Assumes diploid
void genotype::makeHaplotypeLikelihoods(std::vector < float > & HL, const bool first,const float min_gl) const {
	std::array<float,3> tmp;
	float sum;

	for (int l = 0 ; l < n_variants ; l ++) {
		if (!flat[l]) {
			tmp[0]=unphred[GL[3*l+0]];
			tmp[1]=unphred[GL[3*l+1]];
			tmp[2]=unphred[GL[3*l+2]];
			sum = tmp[0] + tmp[1] + tmp[2];
			tmp[0] /= sum;
			tmp[1] /= sum;
			tmp[2] /= sum;

			const int condAllele = first?H1[l]:H0[l];
			HL[2*l+0] = tmp[0+condAllele] / (tmp[0+condAllele] + tmp[1+condAllele]);
			HL[2*l+1] = tmp[1+condAllele] / (tmp[0+condAllele] + tmp[1+condAllele]);

			if (HL[2*l+0] < min_gl) {HL[2*l+0]=min_gl; HL[2*l+1]=1.0f-min_gl;}
			if (HL[2*l+1] < min_gl) {HL[2*l+0]=1.0f-min_gl; HL[2*l+1]=min_gl;}
		}
	}
}

void genotype::sampleHaplotypeH0(const std::vector < float > & HP0) {
	for (int l = 0 ; l < n_variants ; l ++) H0[l] = (rng.getFloat() > HP0[2*l]);
}

//Assumes diploid, HP1 exists.
void genotype::sampleHaplotypeH1(const std::vector < float > & HP1) {
	for (int l = 0 ; l < n_variants ; l ++) H1[l] = (rng.getFloat() > HP1[2*l]);
}

void genotype::storeGenotypePosteriorsAndHaplotypes(const std::vector < float > & HP0)
{
	//haploid case
	std::vector < bool > flag = std::vector < bool > (n_variants, false);

	// Push already variants at which storage already occured in previous iterations
	for (int e = 0 ; e < stored_data.size() ; e ++) {
		int var_idx = stored_data[e].idx;
		float p0 = HP0[2*var_idx+0];
		float p1 = HP0[2*var_idx+1];

		float sc = 1.0f/(p0+p1);
		// Store GP
		stored_data[e].gp0 += p0*sc;
		stored_data[e].gp1 += p1*sc;
		// Store HS
		//H0[var_idx]?_SET32(stored_data[e].hs,2*(stored_cnt%16)+0):_CLR32(stored_data[e].hs,2*(stored_cnt%16)+0);
		//_CLR32(stored_data[e].hs,2*(stored_cnt%16)+1);//just padding
		stored_data[e].hds = false;
		//Flag variant as already stored GP/HS
		flag[var_idx] = true;
	}
	// Push variants that have never been pushed in previous iterations
	for (int l = 0 ; l < n_variants ; l ++) {
		if (!flag[l]) {
			float p0 = HP0[2*l+0];
			float p1 = HP0[2*l+1];
			float sc = 1.0f/(p0+p1);

			// Trigger storage
			if (p0*sc < 0.99999f) {
				int32_t new_hs = 0;
				//H0[l]?_SET32(new_hs,2*(stored_cnt%16)+0):_CLR32(new_hs,2*(stored_cnt%16)+0);
				//_CLR32(new_hs,2*(stored_cnt%16)+1);//just padding
				stored_data.emplace_back(l, p0*sc + (stored_cnt%16)*1.0f, p1*sc, false);
			}
		}
	}
	stored_cnt ++;
}

void genotype::storeGenotypePosteriorsAndHaplotypes(const std::vector < float > & HP0, const std::vector < float > & HP1) {
	//diploid case
	std::vector < bool > flag = std::vector < bool > (n_variants, false);
	// Push already variants at which storage already occured in previous iterations
	for (int e = 0 ; e < stored_data.size() ; e ++)
	{
		int var_idx = stored_data[e].idx;
		float p0 = std::clamp(HP0[2*var_idx+0] * HP1[2*var_idx+0],0.0f,1.0f);
		float p1 = std::clamp(HP0[2*var_idx+0] * HP1[2*var_idx+1] + HP0[2*var_idx+1] * HP1[2*var_idx+0],0.0f,1.0f);
		float p2 = std::clamp(HP0[2*var_idx+1] * HP1[2*var_idx+1],0.0f,1.0f);
		//float sc = 1.0f/(p0+p1+p2);
		// Store GP
		stored_data[e].gp0 += p0/(p0+p1+p2);
		stored_data[e].gp1 += p1/(p0+p1+p2);
		// Store HS
		//H0[var_idx]?_SET32(stored_data[e].hs,2*(stored_cnt%16)+0):_CLR32(stored_data[e].hs,2*(stored_cnt%16)+0);
		//H1[var_idx]?_SET32(stored_data[e].hs,2*(stored_cnt%16)+1):_CLR32(stored_data[e].hs,2*(stored_cnt%16)+1);
		stored_data[e].hds = HP0[2*var_idx+1] < HP1[2*var_idx+1];
		//Flag variant as already stored GP/HS
		flag[var_idx] = true;
	}
	// Push variants that have never been pushed in previous iterations
	for (int l = 0 ; l < n_variants ; l ++) {
		if (!flag[l])
		{
			float p0 = std::clamp(HP0[2*l+0] * HP1[2*l+0],0.0f,1.0f);
			float p1 = std::clamp(HP0[2*l+0] * HP1[2*l+1] + HP0[2*l+1] * HP1[2*l+0],0.0f,1.0f);
			float p2 = std::clamp(HP0[2*l+1] * HP1[2*l+1],0.0f,1.0f);
			//float sc = 1.0f/(p0+p1+p2);

			// Trigger storage
			if (p0/(p0+p1+p2) < 0.99999f) {
				int32_t new_hs = 0;
				//H0[l]?_SET32(new_hs,2*(stored_cnt%16)+0):_CLR32(new_hs,2*(stored_cnt%16)+0);
				//H1[l]?_SET32(new_hs,2*(stored_cnt%16)+1):_CLR32(new_hs,2*(stored_cnt%16)+1);
				stored_data.emplace_back(l, p0/(p0+p1+p2) + stored_cnt*1.0f, p1/(p0+p1+p2), HP0[2*l+1] < HP1[2*l+1]);

			}
		}
	}
	stored_cnt ++;
}

void genotype::sortAndNormAndInferGenotype() {
	//float scale = 1.0f / stored_cnt;

	// Sorting by variant index
	sort(stored_data.begin(), stored_data.end());

	// Normalizing
	for (int e = 0 ; e < stored_data.size() ; e ++) {
		stored_data[e].gp0 /= stored_cnt;
		stored_data[e].gp1 /= stored_cnt;
	}

	// Infer most likely genotypes
	for (int l = 0, e = 0 ; l < n_variants ; ++l) {
		if ((e==stored_data.size()) || (stored_data[e].idx>l)) {		// No storage happened there => GP=(1.0,0.0,0.0)
			H0[l] = false;
			if (ploidy > 1) H1[l] = false;
		} else {
			// Storage happened there => GP !=(1.0, 0.0,0.0)
			if (ploidy > 1)
			{
				switch (stored_data[e].infer())
				{
					case 0: 	H0[l] = false; H1[l] = false; break;
					case 1: 	H0[l] = false; H1[l] = true;  break;
					case 2:		H0[l] = true;  H1[l] = true;  break;
					default: 	H0[l] = false; H1[l] = false; break;
				}
			}
			else H0[l] = stored_data[e].infer_haploid();
			e++;
		}
	}
}
