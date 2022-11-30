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

#include <models/phasing_hmm.h>

phasing_hmm::phasing_hmm(conditioning_set * _C) :
	n_segs(0),n_miss(0),curr_idx_locus(0),curr_abs_locus(0),curr_rel_locus(0),curr_segment_index(0),curr_segment_locus(0),curr_missing_locus(0),probSumT(0), sumHProbs(0), sumDProbs(0),
	nt(0.0f), yt(0.0f)

{
	C = _C;
	EMIT0 = std::vector<aligned_vector32<float> >(3, aligned_vector32<float> (HAP_NUMBER));
	EMIT1 = std::vector<aligned_vector32<float> >(3, aligned_vector32<float> (HAP_NUMBER));
	HProbs = aligned_vector32<float> (HAP_NUMBER*HAP_NUMBER);

	//INIT EMIT0
	EMIT0[0][0] = C->ed_phs; EMIT0[1][0] = C->ed_phs; EMIT0[2][0] = C->ed_phs;
	EMIT0[0][1] = C->ed_phs; EMIT0[1][1] = C->ed_phs; EMIT0[2][1] = C->ee_phs;
	EMIT0[0][2] = C->ed_phs; EMIT0[1][2] = C->ee_phs; EMIT0[2][2] = C->ed_phs;
	EMIT0[0][3] = C->ed_phs; EMIT0[1][3] = C->ee_phs; EMIT0[2][3] = C->ee_phs;
	EMIT0[0][4] = C->ee_phs; EMIT0[1][4] = C->ed_phs; EMIT0[2][4] = C->ed_phs;
	EMIT0[0][5] = C->ee_phs; EMIT0[1][5] = C->ed_phs; EMIT0[2][5] = C->ee_phs;
	EMIT0[0][6] = C->ee_phs; EMIT0[1][6] = C->ee_phs; EMIT0[2][6] = C->ed_phs;
	EMIT0[0][7] = C->ee_phs; EMIT0[1][7] = C->ee_phs; EMIT0[2][7] = C->ee_phs;

	//INIT EMIT1
	EMIT1[0][0] = C->ee_phs; EMIT1[1][0] = C->ee_phs; EMIT1[2][0] = C->ee_phs;
	EMIT1[0][1] = C->ee_phs; EMIT1[1][1] = C->ee_phs; EMIT1[2][1] = C->ed_phs;
	EMIT1[0][2] = C->ee_phs; EMIT1[1][2] = C->ed_phs; EMIT1[2][2] = C->ee_phs;
	EMIT1[0][3] = C->ee_phs; EMIT1[1][3] = C->ed_phs; EMIT1[2][3] = C->ed_phs;
	EMIT1[0][4] = C->ed_phs; EMIT1[1][4] = C->ee_phs; EMIT1[2][4] = C->ee_phs;
	EMIT1[0][5] = C->ed_phs; EMIT1[1][5] = C->ee_phs; EMIT1[2][5] = C->ed_phs;
	EMIT1[0][6] = C->ed_phs; EMIT1[1][6] = C->ed_phs; EMIT1[2][6] = C->ee_phs;
	EMIT1[0][7] = C->ed_phs; EMIT1[1][7] = C->ed_phs; EMIT1[2][7] = C->ed_phs;

	DProbs = std::vector < float >(HAP_NUMBER);
}

phasing_hmm::~phasing_hmm()
{
}

void phasing_hmm::reallocate(const std::vector < bool > & H0, const std::vector < bool > & H1, std::vector < bool > & flat) {
	//SET VARIANT TYPE AND INDEXING
	VAR_TYP.clear();
	VAR_ALT.clear();
	VAR_ABS.clear();
	VAR_REL.clear();
	for (int l = 0, n_het = 0 ; l < C->polymorphic_sites.size() ; l ++) {
		bool a0 = H0[C->polymorphic_sites[l]];
		bool a1 = H1[C->polymorphic_sites[l]];

		if ((!flat[C->polymorphic_sites[l]]) && (!C->lq_flag[C->polymorphic_sites[l]])) {
			if (a0 != a1) {
				VAR_TYP.push_back((char)(n_het % 3));
				VAR_ALT.push_back(a0);
				VAR_ABS.push_back(C->polymorphic_sites[l]);
				VAR_REL.push_back(l);
				n_het ++;
			} else {
				VAR_TYP.push_back(VAR_PEAK_HOM);
				VAR_ALT.push_back(a0);
				VAR_ABS.push_back(C->polymorphic_sites[l]);
				VAR_REL.push_back(l);
			}
		} else {
			if (a0 != a1) {
				VAR_TYP.push_back(VAR_FLAT_HET);
				VAR_ALT.push_back(a0);
				VAR_ABS.push_back(C->polymorphic_sites[l]);
				VAR_REL.push_back(l);
			}
		}
	}

	//COMPUTE SEGMENTATION
	int nv = 0;
	segments.clear();
	n_miss = 0;
	for (int l = 0, n_hets = 0 ; l < VAR_TYP.size() ;) {
		n_hets += (VAR_TYP[l] >= VAR_PEAK_HET);
		n_miss += (VAR_TYP[l] == VAR_FLAT_HET);
		if (n_hets == 4) {
			segments.push_back(nv);
			n_hets = 0;
			nv = 0;
		} else {
			nv ++;
			l++;
		}
	}
	segments.push_back(nv);
	n_segs = segments.size();
	dip_sampled = std::vector < int > (n_segs, -1);

	//REALLOCATE MEMORY
	prob.resize(C->n_states * HAP_NUMBER);
	probSumH.resize(HAP_NUMBER);
	probSumK.resize(C->n_states);

	//phasingProb = std::vector < aligned_vector32 < float >  > (n_segs, aligned_vector32 < float >  (C->n_states * HAP_NUMBER, 0.0f));
	phasingProb.resize(n_segs*C->n_states * HAP_NUMBER);
	//phasingProbSum = std::vector < aligned_vector32 < float >  > (n_segs, aligned_vector32 < float >  (HAP_NUMBER, 0.0f));
	phasingProbSum.resize(n_segs*HAP_NUMBER);
	phasingProbSumSum.resize(n_segs);

	//imputeProb = std::vector < aligned_vector32 < float >  > (n_miss, aligned_vector32 < float >  (C->n_states * HAP_NUMBER, 0.0f));
	imputeProb.resize(n_miss*C->n_states * HAP_NUMBER);
	//imputeProbSum = std::vector < aligned_vector32 < float >  > (n_miss, aligned_vector32 < float >  (HAP_NUMBER, 0.0f));
	imputeProbSum.resize(n_miss * HAP_NUMBER);
	imputeProbSumSum.resize(n_miss);
	imputeProbOf1s.resize(n_miss * HAP_NUMBER);

	std::fill(phasingProb.begin(), phasingProb.end(), 0.0f);
	std::fill(phasingProbSum.begin(), phasingProbSum.end(), 0.0f);
	std::fill(imputeProb.begin(), imputeProb.end(), 0.0f);
	std::fill(imputeProbSum.begin(), imputeProbSum.end(), 0.0f);
}

void phasing_hmm::forward()
{
	curr_segment_index = 0;
	curr_segment_locus = 0;
	curr_missing_locus = 0;

	for (curr_idx_locus = 0 ; curr_idx_locus < VAR_TYP.size() ; curr_idx_locus ++) {
		curr_abs_locus = VAR_ABS[curr_idx_locus];
		curr_rel_locus = VAR_REL[curr_idx_locus];

		yt = 0.0f;
		if (curr_idx_locus) yt = C->getTransition(VAR_ABS[curr_idx_locus-1], VAR_ABS[curr_idx_locus]);
		nt = 1.0f - yt;

		if (VAR_TYP[curr_idx_locus] >= VAR_PEAK_HET) {
			if (curr_idx_locus == 0) INIT_PEAK_HET(VAR_TYP[curr_idx_locus]);
			else if (curr_segment_locus != 0) RUN_PEAK_HET(VAR_TYP[curr_idx_locus]);
			else COLLAPSE_PEAK_HET(VAR_TYP[curr_idx_locus]);
		} else if (VAR_TYP[curr_idx_locus] == VAR_PEAK_HOM) {
			if (curr_idx_locus == 0) INIT_PEAK_HOM(VAR_ALT[curr_idx_locus]);
			else if (curr_segment_locus != 0) RUN_PEAK_HOM(VAR_ALT[curr_idx_locus]);
			else COLLAPSE_PEAK_HOM(VAR_ALT[curr_idx_locus]);
		} else if (VAR_TYP[curr_idx_locus] == VAR_FLAT_HET) {
			if (curr_idx_locus == 0) INIT_FLAT_HET();
			else if (curr_segment_locus != 0) RUN_FLAT_HET();
			else COLLAPSE_FLAT_HET();
		} else vrb.error("Unknown variant type in phasing forward pass at variant " + stb.str(curr_idx_locus));

		if (curr_segment_locus == (segments[curr_segment_index] - 1)) SUMK();

		//PHASE COMMON HETS
		if ((curr_segment_locus == (segments[curr_segment_index]-1)) && (curr_idx_locus != (VAR_TYP.size()-1))) {
			const bool ret1 = TRANS_HAP(); assert(!ret1);
			const bool ret2 = SAMPLE_DIP(); assert(!ret2);
		}

		//PHASE RARE HETS
		if (VAR_TYP[curr_idx_locus] == VAR_FLAT_HET) {
			IMPUTE_FLAT_HET();
			curr_missing_locus++;
		}

		//UPDATE ITERATORS
		curr_segment_locus ++;
		if (curr_segment_locus >= segments[curr_segment_index]) {
			curr_segment_index++;
			curr_segment_locus = 0;
		}
	}
}

void phasing_hmm::backward() {
	curr_segment_index = n_segs - 1;
	curr_segment_locus = segments.back() - 1;
	curr_missing_locus = n_miss - 1;

	for (curr_idx_locus = VAR_TYP.size()-1 ; curr_idx_locus >= 0 ; curr_idx_locus --) {
		curr_abs_locus = VAR_ABS[curr_idx_locus];
		curr_rel_locus = VAR_REL[curr_idx_locus];

		yt = 0.0f;
		if (curr_idx_locus < (VAR_TYP.size()-1)) yt = C->getTransition(VAR_ABS[curr_idx_locus], VAR_ABS[curr_idx_locus+1]);
		nt = 1.0f - yt;

		if (VAR_TYP[curr_idx_locus] >= VAR_PEAK_HET) {
			if (curr_idx_locus == (VAR_TYP.size()-1)) INIT_PEAK_HET(VAR_TYP[curr_idx_locus]);
			else if (curr_segment_locus != (segments[curr_segment_index]-1)) RUN_PEAK_HET(VAR_TYP[curr_idx_locus]);
			else COLLAPSE_PEAK_HET(VAR_TYP[curr_idx_locus]);
		} else if (VAR_TYP[curr_idx_locus] == VAR_PEAK_HOM) {
			if (curr_idx_locus == (VAR_TYP.size()-1)) INIT_PEAK_HOM(VAR_ALT[curr_idx_locus]);
			else if (curr_segment_locus != (segments[curr_segment_index]-1)) RUN_PEAK_HOM(VAR_ALT[curr_idx_locus]);
			else COLLAPSE_PEAK_HOM(VAR_ALT[curr_idx_locus]);
		} else if (VAR_TYP[curr_idx_locus] == VAR_FLAT_HET) {
			if (curr_idx_locus == (VAR_TYP.size()-1)) INIT_FLAT_HET();
			else if (curr_segment_locus != (segments[curr_segment_index]-1)) RUN_FLAT_HET();
			else COLLAPSE_FLAT_HET();
		} else vrb.error("Unknown variant type in phasing backward pass at variant " + stb.str(curr_idx_locus));

		if (curr_segment_locus == 0) {
			SUMK();
			//phasingProb[curr_segment_index] = prob;
			std::copy(prob.begin(), prob.end(), phasingProb.begin()+curr_segment_index*C->n_states*HAP_NUMBER);
			//phasingProbSum[curr_segment_index] = probSumH;
			std::copy(probSumH.begin(), probSumH.end(), phasingProbSum.begin()+curr_segment_index*HAP_NUMBER);
			phasingProbSumSum[curr_segment_index] = probSumT;
		}

		//STORE PROBS FOR PHASING RARE HETS
		if (VAR_TYP[curr_idx_locus] == VAR_FLAT_HET) {
			//imputeProb[curr_missing_locus] = prob;
			std::copy(prob.begin(), prob.end(), imputeProb.begin()+curr_missing_locus*C->n_states*HAP_NUMBER);
			//imputeProbSum[curr_missing_locus] = probSumH;
			std::copy(probSumH.begin(), probSumH.end(), imputeProbSum.begin()+curr_missing_locus*HAP_NUMBER);
			imputeProbSumSum[curr_missing_locus] = probSumT;
			curr_missing_locus--;
		}

		curr_segment_locus--;
		if (curr_segment_locus < 0 && curr_segment_index > 0) {
			curr_segment_index--;
			curr_segment_locus = segments[curr_segment_index] - 1;
		}
	}
}

void phasing_hmm::rephaseHaplotypes(std::vector < bool > & H0, std::vector < bool > & H1, std::vector < bool > & flat) {
	reallocate(H0, H1, flat);
	backward();
	for (int d = 0 ; d < HAP_NUMBER ; d ++) {
		//DProbs[d] = (phasingProbSum[0][d] / phasingProbSumSum[0]) * (phasingProbSum[0][HAP_NUMBER - d - 1] / phasingProbSumSum[0]);
		DProbs[d] = (phasingProbSum[d] / phasingProbSumSum[0]) * (phasingProbSum[HAP_NUMBER - d - 1] / phasingProbSumSum[0]);
		sumDProbs += DProbs[d];
	}
	dip_sampled[0] = rng.sample(DProbs, sumDProbs);
	forward();
	curr_segment_index = 0;
	curr_segment_locus = 0;
	curr_missing_locus = 0;
	for (curr_idx_locus = 0  ; curr_idx_locus < VAR_TYP.size() ; curr_idx_locus ++) {
		if (VAR_TYP[curr_idx_locus] >= VAR_PEAK_HET) {
			const int idx_h0 = dip_sampled[curr_segment_index];
			const int idx_h1 = HAP_NUMBER - dip_sampled[curr_segment_index] - 1;
			const bool a0 = ALLELE(HAP_NUMBER - idx_h0 - 1, 2 - VAR_TYP[curr_idx_locus]);
			const bool a1 = ALLELE(HAP_NUMBER - idx_h1 - 1, 2 - VAR_TYP[curr_idx_locus]);
			H0[VAR_ABS[curr_idx_locus]] = a0;
			H1[VAR_ABS[curr_idx_locus]] = a1;
		} else if (VAR_TYP[curr_idx_locus] == VAR_FLAT_HET) {
			const int idx_h0 = dip_sampled[curr_segment_index];
			const int idx_h1 = HAP_NUMBER - dip_sampled[curr_segment_index] - 1;
			const float h0a0 = (1.0f - imputeProbOf1s[curr_missing_locus*HAP_NUMBER+idx_h0]);
			const float h0a1 = imputeProbOf1s[curr_missing_locus*HAP_NUMBER+idx_h0];
			const float h1a0 = (1.0f - imputeProbOf1s[curr_missing_locus*HAP_NUMBER+idx_h1]);
			const float h1a1 = imputeProbOf1s[curr_missing_locus*HAP_NUMBER+idx_h1];
			float p01 = (h0a0*C->ee_phs + h0a1*C->ed_phs) * (h1a1*C->ee_phs + h1a0*C->ed_phs);
			float p10 = (h0a1*C->ee_phs + h0a0*C->ed_phs) * (h1a0*C->ee_phs + h1a1*C->ed_phs);
			float sum = p01+p10;
			p01 = std::clamp(p01/sum, 0.0f,1.0f);
			p10 = std::clamp(p10/sum, 0.0f,1.0f);
			sum = p01+p10;
			const bool rf = (rng.getFloat()*sum) < p01;
			H0[VAR_ABS[curr_idx_locus]] = rf;
			H0[VAR_ABS[curr_idx_locus]] = !rf;
			curr_missing_locus++;
		}
		curr_segment_locus ++;
		if (curr_segment_locus >= segments[curr_segment_index]) {
			curr_segment_index++;
			curr_segment_locus = 0;
		}
	}
	//Shuffling het at monomorphic
	for (int l = 0 ; l < C->monomorphic_sites.size() ; l ++) {
		if (H0[C->monomorphic_sites[l]] != H1[C->monomorphic_sites[l]])
		{
			const bool rf = rng.getFloat(0.0f,1.0f) < 0.5f;
			H0[C->monomorphic_sites[l]] = rf;
			H1[C->monomorphic_sites[l]] = !rf;
		}
	}
}
