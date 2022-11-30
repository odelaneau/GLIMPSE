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

#include <models/imputation_hmm.h>

inline
float horizontal_add (const __m256& a)
{
    __m128 vlow = _mm256_castps256_ps128(a);
    __m128 vhigh = _mm256_extractf128_ps(a, 1); // high 128
   vlow = _mm_add_ps(vlow, vhigh);     // add the low 128
   __m128 shuf = _mm_movehdup_ps(vlow);        // broadcast elements 3,1 to 2,0
   __m128 sums = _mm_add_ps(vlow, shuf);
   shuf = _mm_movehl_ps(shuf, sums); // high half -> low half
   sums = _mm_add_ss(sums, shuf);    // (no wasted instructions, and all of them are the 4B minimum)
   return _mm_cvtss_f32(sums);
}

imputation_hmm::imputation_hmm(conditioning_set * _C) {
	C = _C;
	modK=0;
	Emissions = aligned_vector32 < float > (2*C->n_tot_sites);
}

imputation_hmm::~imputation_hmm() {
	Alpha.clear();
	AlphaSum.clear();
	Emissions.clear();
}

void imputation_hmm::resize()
{
	modK = ((C->n_states / 8) + (C->n_states % 8 ? 1 : 0))*8;
	AlphaSum.resize(C->polymorphic_sites.size(), 0.0f);
	Alpha.resize(C->polymorphic_sites.size() * modK, 0.0f);
	Beta.resize(modK);
}

void imputation_hmm::init(const std::vector < float > & HL)
{
	for (int l = 0 ; l < C->n_tot_sites ; l ++)
	{
		float p0 = HL[2*l+0] * C->ee_imp + HL[2*l+1] * C->ed_imp;
		float p1 = HL[2*l+0] * C->ed_imp + HL[2*l+1] * C->ee_imp;
		Emissions[2*l+0] = p0 / (p0+p1);
		Emissions[2*l+1] = p1 / (p0+p1);
	}

}

void imputation_hmm::computePosteriors(const std::vector < float > & HL, std::vector < bool > & flat, std::vector < float > & HP) {
	resize();
	init(HL);
	forward(flat);
	backward(HL, flat, HP);
}

void imputation_hmm::forward(std::vector < bool > & flat) {
	const __m256i _vshift_count = _mm256_set_epi32(31,30,29,28,27,26,25,24);
	const unsigned int nstates = C->n_states;
	const unsigned int nstatesMD8 = (nstates / 8) * 8;

	for (int l = 0 ; l < C->polymorphic_sites.size() ; l ++)
	{
		AlphaSum[l] = 0.0f;
		if (flat[C->polymorphic_sites[l]] || C->lq_flag[C->polymorphic_sites[l]])
		{
			if (l == 0)
			{
				fill(Alpha.begin(), Alpha.begin()+modK, 1.0f / nstates);
				AlphaSum[l] = 1.0f;
			}
			else
			{
				const float fact1 = C->t[l-1] / nstates;
				const float fact2 = C->nt[l-1] / AlphaSum[l-1];
				const __m256 _fact1 = _mm256_set1_ps(fact1);
				const __m256 _fact2 = _mm256_set1_ps(fact2);
				__m256 _sum = _mm256_set1_ps(0.0f);
				int k = 0;
				for (; k < nstatesMD8; k += 8)
				{
					const __m256 _prob_prev = _mm256_load_ps(&Alpha[(l-1)*modK+k]);
					const __m256 _prob_curr = _mm256_fmadd_ps(_prob_prev, _fact2, _fact1);
					_sum = _mm256_add_ps(_sum, _prob_curr);
					_mm256_store_ps(&Alpha[l*modK+k], _prob_curr);
				}
				if (k) AlphaSum[l] = horizontal_add(_sum);
				for (int offset = nstatesMD8; offset < nstates ; offset ++) {
					Alpha[l*modK+offset] = (Alpha[(l-1)*modK+offset] * fact2 + fact1);
					AlphaSum[l] += Alpha[l*modK+offset];
				}
			}
		}
		else
		{
			const std::array<float,2> emit = {Emissions[2*C->polymorphic_sites[l]+0], Emissions[2*C->polymorphic_sites[l]+1]};
			const __m256 _emit0 = _mm256_set1_ps(emit[0]);
			const __m256 _emit1 = _mm256_set1_ps(emit[1]);
			if (l == 0)
			{
				const float fact1 = 1.0f / nstates;
				const __m256 _fact1 = _mm256_set1_ps(fact1);
				__m256 _sum = _mm256_set1_ps(0.0f);
				int k = 0;
				for (; k < nstatesMD8; k += 8)
				{
					const __m256i _bcst = _mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k));
					const __m256i _mask = _mm256_sllv_epi32(_bcst, _vshift_count);
					const __m256 _emiss = _mm256_blendv_ps (_emit0, _emit1, _mm256_castsi256_ps(_mask));
					const __m256 _prob_curr = _mm256_mul_ps(_emiss, _fact1);
					_sum = _mm256_add_ps(_sum, _prob_curr);
					_mm256_store_ps(&Alpha[l*modK+k], _prob_curr);
				}
				if (k) AlphaSum[l] = horizontal_add(_sum);
				for (int offset = nstatesMD8; offset < nstates ; offset ++)
				{
					Alpha[l*modK+offset] = emit[C->Hvar.get(l, offset)] * fact1;
					AlphaSum[l] += Alpha[l*modK+offset];
				}
			}
			else
			{
				const float fact1 = C->t[l-1] / nstates;
				const float fact2 = C->nt[l-1] / AlphaSum[l-1];//AlphaSum2[l-1];// AlphaSum[l-1];
				const __m256 _fact1 = _mm256_set1_ps(fact1);
				const __m256 _fact2 = _mm256_set1_ps(fact2);
				__m256 _sum = _mm256_set1_ps(0.0f);

				int k = 0;
				for (; k < nstatesMD8; k += 8)
				{
					const __m256i _mask = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k)), _vshift_count);
					const __m256 _emiss = _mm256_blendv_ps (_emit0, _emit1, _mm256_castsi256_ps(_mask));
					const __m256 _prob_prev = _mm256_load_ps(&Alpha[(l-1)*modK+k]);
					const __m256 _prob_temp = _mm256_fmadd_ps(_prob_prev, _fact2, _fact1);
					const __m256 _prob_curr = _mm256_mul_ps(_prob_temp, _emiss);
					_sum = _mm256_add_ps(_sum, _prob_curr);
					_mm256_store_ps(&Alpha[l*modK+k], _prob_curr);
				}
				if (k) AlphaSum[l] = horizontal_add(_sum);
				for (int offset = nstatesMD8; offset < nstates ; offset ++)
				{
					Alpha[l*modK+offset] = (Alpha[(l-1)*modK+offset]*fact2+fact1)*emit[C->Hvar.get(l, offset)];
					AlphaSum[l] += Alpha[l*modK+offset];
				}
			}
		}
	}
}

void imputation_hmm::backward(const std::vector < float > & HL, std::vector < bool > & flat, std::vector < float > & HP)
{
	float betaSum = 0.0f, betaSumNext = 0.0f;
	std::array<float,2> prob_hid, prob_obs;
	const __m256i _vshift_count = _mm256_set_epi32(31,30,29,28,27,26,25,24);
	const unsigned int nstates = C->n_states;
	const unsigned int nstatesMD8 = (nstates / 8) * 8;
	const __m256 _zero = _mm256_set1_ps(0.0f);
	const __m256 _one = _mm256_set1_ps(1.0f);
	fill(Beta.begin(), Beta.end(), 1.0f);
	__m256 _sum,  _prob0, _prob1;
	for (int l = C->polymorphic_sites.size()-1 ; l >= 0 ; l --)
	{
		betaSum=0.0f;
		prob_hid[0]=0.0f;
		prob_hid[1]=0.0f;

		_sum = _mm256_set1_ps(0.0f);
		_prob0 = _mm256_set1_ps(0.0f);
		_prob1 = _mm256_set1_ps(0.0f);

		if (flat[C->polymorphic_sites[l]] || C->lq_flag[C->polymorphic_sites[l]])
		{
			if (l == (C->polymorphic_sites.size()-1))
			{
				const float fact1 = 1.0f / nstates;
				const __m256 _fact1 = _mm256_set1_ps(fact1);
				int k = 0;
				for (; k < nstatesMD8; k += 8)
				{
					const __m256i _mask_curr = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k)), _vshift_count);
					const __m256 _mask0 = _mm256_blendv_ps (_one, _zero, _mm256_castsi256_ps(_mask_curr));
					const __m256 _mask1 = _mm256_blendv_ps (_zero, _one, _mm256_castsi256_ps(_mask_curr));
					const __m256 _alphas = _mm256_load_ps(&Alpha[l*modK+k]);
					_prob0 = _mm256_add_ps(_prob0, _mm256_mul_ps(_alphas, _mask0));
					_prob1 = _mm256_add_ps(_prob1, _mm256_mul_ps(_alphas, _mask1));
					_mm256_store_ps(&Beta[k], _fact1);
				}
				if (k)
				{
					prob_hid[0] = horizontal_add(_prob0);
					prob_hid[1] = horizontal_add(_prob1);
				}
				for (int offset = nstatesMD8; offset < nstates ; offset ++)
				{
					Beta[offset] = fact1;
					prob_hid[C->Hvar.get(l, offset)] += Alpha[l*modK+offset];
				}
				betaSum = 1.0f;
			}
			else
			{
				const float fact1 = C->t[l] / nstates;
				const float fact2 = C->nt[l] / betaSumNext;
				const __m256 _fact1 = _mm256_set1_ps(fact1);
				const __m256 _fact2 = _mm256_set1_ps(fact2);

				int k = 0;
				for (; k < nstatesMD8; k += 8)
				{
					const __m256i _mask_curr = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k)), _vshift_count);
					const __m256 _mask0 = _mm256_blendv_ps (_one, _zero, _mm256_castsi256_ps(_mask_curr));
					const __m256 _mask1 = _mm256_blendv_ps (_zero, _one, _mm256_castsi256_ps(_mask_curr));
					const __m256 _prob_prev = _mm256_load_ps(&Beta[k]);
					const __m256 _prob_curr = _mm256_fmadd_ps(_prob_prev, _fact2, _fact1);
					const __m256 _alphas = _mm256_load_ps(&Alpha[l*modK+k]);
					const __m256 _dotprod = _mm256_mul_ps(_alphas, _prob_curr);
					_prob0 = _mm256_add_ps(_prob0, _mm256_mul_ps(_dotprod, _mask0));
					_prob1 = _mm256_add_ps(_prob1, _mm256_mul_ps(_dotprod, _mask1));
					_sum = _mm256_add_ps(_sum, _prob_curr);
					_mm256_store_ps(&Beta[k], _prob_curr);
				}
				if (k)
				{
					prob_hid[0] = horizontal_add(_prob0);
					prob_hid[1] = horizontal_add(_prob1);
					betaSum = horizontal_add(_sum);
				}
				for (int offset = nstatesMD8; offset < nstates ; offset ++)
				{
					Beta[offset] = Beta[offset] * fact2 + fact1;
					prob_hid[C->Hvar.get(l, offset)] += Alpha[l*modK+offset] * Beta[offset];
					betaSum += Beta[offset];
				}
			}

			prob_obs[0] = (prob_hid[0]*C->ee_imp + prob_hid[1]*C->ed_imp);
			prob_obs[1] = (prob_hid[0]*C->ed_imp + prob_hid[1]*C->ee_imp);
			if (!flat[C->polymorphic_sites[l]])
			{
				prob_obs[0] *= HL[2*C->polymorphic_sites[l]+0];
				prob_obs[1] *= HL[2*C->polymorphic_sites[l]+1];
			}
		}
		else
		{
			std::array<float,2> emit = {Emissions[2*C->polymorphic_sites[l]+0], Emissions[2*C->polymorphic_sites[l]+1]};
			const __m256 _emit0 = _mm256_set1_ps(emit[0]);
			const __m256 _emit1 = _mm256_set1_ps(emit[1]);

			if (l == (C->polymorphic_sites.size()-1))
			{
				const float fact1 = 1.0f / nstates;
				const __m256 _fact1 = _mm256_set1_ps(fact1);

				int k = 0;
				for (; k < nstatesMD8; k += 8)
				{
					const __m256i _mask_curr = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k)), _vshift_count);
					const __m256 _emiss = _mm256_blendv_ps (_emit0, _emit1, _mm256_castsi256_ps(_mask_curr));
					const __m256 _mask0 = _mm256_blendv_ps (_one, _zero, _mm256_castsi256_ps(_mask_curr));
					const __m256 _mask1 = _mm256_blendv_ps (_zero, _one, _mm256_castsi256_ps(_mask_curr));
					const __m256 _alphas = _mm256_load_ps(&Alpha[l*modK+k]);
					const __m256 _prob_next = _mm256_mul_ps(_emiss, _fact1);
					_prob0 = _mm256_add_ps(_prob0, _mm256_mul_ps(_alphas, _mask0));
					_prob1 = _mm256_add_ps(_prob1, _mm256_mul_ps(_alphas, _mask1));
					_sum= _mm256_add_ps(_sum, _prob_next);
					_mm256_store_ps(&Beta[k], _prob_next);
				}
				if (k)
				{
					prob_hid[0] = horizontal_add(_prob0);
					prob_hid[1] = horizontal_add(_prob1);
					betaSum = horizontal_add(_sum);
				}
				for (int offset = nstatesMD8; offset < nstates ; offset ++)
				{
					prob_hid[C->Hvar.get(l, offset)] += Alpha[l*modK+offset];
					Beta[offset] = emit[C->Hvar.get(l, offset)] * fact1;
					betaSum += Beta[offset];
				}
			}
			else
			{
				const float fact1 = C->t[l] / nstates;
				const float fact2 = C->nt[l] / betaSumNext;
				const __m256 _fact1 = _mm256_set1_ps(fact1);
				const __m256 _fact2 = _mm256_set1_ps(fact2);

				int k = 0;
				for (; k < nstatesMD8; k += 8)
				{
					const __m256i _mask_curr = _mm256_sllv_epi32(_mm256_set1_epi32((unsigned int )C->Hvar.getByte(l, k)), _vshift_count);
					const __m256 _emiss = _mm256_blendv_ps (_emit0, _emit1, _mm256_castsi256_ps(_mask_curr));
					const __m256 _mask0 = _mm256_blendv_ps (_one, _zero, _mm256_castsi256_ps(_mask_curr));
					const __m256 _mask1 = _mm256_blendv_ps (_zero, _one, _mm256_castsi256_ps(_mask_curr));
					const __m256 _prob_prev = _mm256_load_ps(&Beta[k]);
					const __m256 _alphas = _mm256_load_ps(&Alpha[l*modK+k]);
					const __m256 _prob_curr = _mm256_fmadd_ps(_prob_prev, _fact2, _fact1);
					const __m256 _dotprod = _mm256_mul_ps(_alphas, _prob_curr);
					const __m256 _prob_next = _mm256_mul_ps(_prob_curr, _emiss);
					_prob0 = _mm256_add_ps(_prob0, _mm256_mul_ps(_dotprod, _mask0));
					_prob1 = _mm256_add_ps(_prob1, _mm256_mul_ps(_dotprod, _mask1));
					_sum= _mm256_add_ps(_sum, _prob_next);
					_mm256_store_ps(&Beta[k], _prob_next);
				}
				if (k)
				{
					prob_hid[0] = horizontal_add(_prob0);
					prob_hid[1] = horizontal_add(_prob1);
					betaSum = horizontal_add(_sum);
				}
				for (int offset = nstatesMD8; offset < nstates ; offset ++)
				{
					Beta[offset] = Beta[offset] * fact2 + fact1;
					prob_hid[C->Hvar.get(l, offset)] += Alpha[l*modK+offset] * Beta[offset];
					Beta[offset] *= emit[C->Hvar.get(l, offset)];
					betaSum += Beta[offset];
				}
			}
			prob_hid[0] /= emit[0];
			prob_hid[1] /= emit[1];
			prob_obs[0] = (prob_hid[0]*C->ee_imp + prob_hid[1]*C->ed_imp) * HL[2*C->polymorphic_sites[l]+0];
			prob_obs[1] = (prob_hid[0]*C->ed_imp + prob_hid[1]*C->ee_imp) * HL[2*C->polymorphic_sites[l]+1];

		}
		HP[2*C->polymorphic_sites[l]+0] = prob_obs[0] / (prob_obs[0] + prob_obs[1]);
		HP[2*C->polymorphic_sites[l]+1] = prob_obs[1] / (prob_obs[0] + prob_obs[1]);
		betaSumNext = betaSum;
	}
	// Monomorphic sites
	for (int l = 0 ; l < C->monomorphic_sites.size() ; l ++)
	{
		prob_obs[C->major_alleles[C->monomorphic_sites[l]]] = C->ee_imp;
		prob_obs[!C->major_alleles[C->monomorphic_sites[l]]] = C->ed_imp;
		if (!flat[C->monomorphic_sites[l]])
		{
			prob_obs[0] *= HL[2*C->monomorphic_sites[l]+0];
			prob_obs[1] *= HL[2*C->monomorphic_sites[l]+1];
		}
		HP[2*C->monomorphic_sites[l]+0] = prob_obs[0] / (prob_obs[0] + prob_obs[1]);
		HP[2*C->monomorphic_sites[l]+1] = prob_obs[1] / (prob_obs[0] + prob_obs[1]);
	}
}
