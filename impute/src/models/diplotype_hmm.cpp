#include <models/diplotype_hmm.h>

diplotype_hmm::diplotype_hmm(haplotype_set * _H, conditioning_set * _C) {
	H = _H;
	C = _C;

	//INIT EMIT0
	EMIT0[0][0] = C->ed; EMIT0[1][0] = C->ed; EMIT0[2][0] = C->ed;
	EMIT0[0][1] = C->ed; EMIT0[1][1] = C->ed; EMIT0[2][1] = C->ee;
	EMIT0[0][2] = C->ed; EMIT0[1][2] = C->ee; EMIT0[2][2] = C->ed;
	EMIT0[0][3] = C->ed; EMIT0[1][3] = C->ee; EMIT0[2][3] = C->ee;
	EMIT0[0][4] = C->ee; EMIT0[1][4] = C->ed; EMIT0[2][4] = C->ed;
	EMIT0[0][5] = C->ee; EMIT0[1][5] = C->ed; EMIT0[2][5] = C->ee;
	EMIT0[0][6] = C->ee; EMIT0[1][6] = C->ee; EMIT0[2][6] = C->ed;
	EMIT0[0][7] = C->ee; EMIT0[1][7] = C->ee; EMIT0[2][7] = C->ee;

	//INIT EMIT1
	EMIT1[0][0] = C->ee; EMIT1[1][0] = C->ee; EMIT1[2][0] = C->ee;
	EMIT1[0][1] = C->ee; EMIT1[1][1] = C->ee; EMIT1[2][1] = C->ed;
	EMIT1[0][2] = C->ee; EMIT1[1][2] = C->ed; EMIT1[2][2] = C->ee;
	EMIT1[0][3] = C->ee; EMIT1[1][3] = C->ed; EMIT1[2][3] = C->ed;
	EMIT1[0][4] = C->ed; EMIT1[1][4] = C->ee; EMIT1[2][4] = C->ee;
	EMIT1[0][5] = C->ed; EMIT1[1][5] = C->ee; EMIT1[2][5] = C->ed;
	EMIT1[0][6] = C->ed; EMIT1[1][6] = C->ed; EMIT1[2][6] = C->ee;
	EMIT1[0][7] = C->ed; EMIT1[1][7] = C->ed; EMIT1[2][7] = C->ed;

	//
	BetaSum = std::vector < float > (HAP_NUMBER, 0.0);
	probSumH1 = std::vector < float > (HAP_NUMBER, 1.0);
	probSumH2 = std::vector < float > (HAP_NUMBER, 1.0);
	probSumT1 = 1.0;
	probSumT2 = 1.0;
}

diplotype_hmm::~diplotype_hmm() {
}

void diplotype_hmm::reallocate(std::vector < bool > & H0, std::vector < bool > & H1) {
	//SET HET AND ALT
	HET = std::vector < char > (C->n_sites, -1);
	ALT = std::vector < bool > (C->n_sites, false);
	for (int l = 0, n_het = 0 ; l < C->n_sites ; l ++) {
		if (H0[l] != H1[l]) {
			HET[l] = (char)(n_het % 3);
			ALT[l] = false;
			n_het ++;
		} else {
			HET[l] = -1;
			ALT[l] = H0[l];
		}
	}

	//COMPUTE SEGMENTATION
	int nv = 0;
	segments = std::vector < int > ();
	for (int l = 0, n_hets = 0 ; l < C->n_sites ;) {
		n_hets += (HET[l] >= 0);
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

	//REALLOCATE MEMORY
	Alpha.resize(n_segs * C->n_states * HAP_NUMBER);
	Beta.resize(n_segs * C->n_states * HAP_NUMBER);
	AlphaSum.resize(n_segs * HAP_NUMBER);
	prob1.resize(C->n_states * HAP_NUMBER);
	prob2.resize(C->n_states * HAP_NUMBER);
	probSumK1.resize(C->n_states);
	probSumK2.resize(C->n_states);
}

void diplotype_hmm::forward() {
	curr_segment_index = 0;
	curr_segment_locus = 0;
	for (curr_locus = 0 ; curr_locus < C->n_sites ; curr_locus ++) {
		bool paired = (curr_locus % 2 == 0);

		paired?INIT2():INIT1();
		if (curr_locus != 0 && curr_segment_locus == 0) paired?COLLAPSE2(true):COLLAPSE1(true);
		if (curr_locus != 0 && curr_segment_locus != 0) paired?RUN2(true):RUN1(true);
		paired?SUM2():SUM1();

		if (curr_segment_locus == segments[curr_segment_index] - 1) paired?SUMK2():SUMK1();
		if (curr_segment_locus == segments[curr_segment_index] - 1) {
			if (paired) {
				copy(prob2.begin(), prob2.begin() + C->n_states*HAP_NUMBER, Alpha.begin() + curr_segment_index*C->n_states*HAP_NUMBER);
				copy(probSumH2.begin(), probSumH2.end(), AlphaSum.begin() + curr_segment_index*HAP_NUMBER);
			} else {
				copy(prob1.begin(), prob1.begin() + C->n_states*HAP_NUMBER, Alpha.begin() + curr_segment_index*C->n_states*HAP_NUMBER);
				copy(probSumH1.begin(), probSumH1.end(), AlphaSum.begin() + curr_segment_index*HAP_NUMBER);
			}
		}
		curr_segment_locus ++;
		if (curr_segment_locus >= segments[curr_segment_index]) {
			curr_segment_index++;
			curr_segment_locus = 0;
		}
	}
}

void diplotype_hmm::backward() {
	curr_segment_index = n_segs - 1;
	curr_segment_locus = segments.back() - 1;
	for (curr_locus = C->n_sites - 1 ; curr_locus >= 0 ; curr_locus--) {
		bool paired = (curr_locus % 2 == 0);

		paired?INIT2():INIT1();
		if (curr_locus != C->n_sites - 1 && curr_segment_locus == segments[curr_segment_index] - 1) paired?COLLAPSE2(false):COLLAPSE1(false);
		if (curr_locus != C->n_sites - 1 && curr_segment_locus != segments[curr_segment_index] - 1) paired?RUN2(false):RUN1(false);
		paired?SUM2():SUM1();

		if (curr_segment_locus == 0) paired?SUMK2():SUMK1();
		if (curr_segment_locus == 0 && curr_locus != (C->n_sites - 1)) {
			if (paired) copy(prob2.begin(), prob2.begin() + C->n_states*HAP_NUMBER, Beta.begin() + curr_segment_index*C->n_states*HAP_NUMBER);
			else copy(prob1.begin(), prob1.begin() + C->n_states*HAP_NUMBER, Beta.begin() + curr_segment_index*C->n_states*HAP_NUMBER);
		}
		if (curr_locus == 0) BetaSum=(paired?probSumH2:probSumH1);
		curr_segment_locus--;
		if (curr_segment_locus < 0 && curr_segment_index > 0) {
			curr_segment_index--;
			curr_segment_locus = segments[curr_segment_index] - 1;
		}
	}
}

void diplotype_hmm::rephaseHaplotypes(std::vector < bool > & H0, std::vector < bool > & H1) {
	tac.clock();
	reallocate(H0, H1);
	//vrb.bullet("DIP Resize (" + stb.str(tac.rel_time()*1.0, 1) + "ms)");
	forward();
	//vrb.bullet("DIP Forward (" + stb.str(tac.rel_time()*1.0, 1) + "ms)");tac.clock();
	backward();
	//vrb.bullet("DIP Backward (" + stb.str(tac.rel_time()*1.0, 1) + "ms)");tac.clock();

	std::vector < int > dip_sampled = std::vector < int > (n_segs, -1);
	std::vector < double > dip_probs = std::vector < double > (HAP_NUMBER, 0.0);
	for (curr_segment_index = 0, curr_locus = 0 ; curr_segment_index < n_segs ; curr_segment_index ++) {
		float sumHap = 0.0, sumDip = 0.0;
		if (curr_segment_index == 0) {
			for (int h = 0 ; h < HAP_NUMBER ; h ++) sumHap += BetaSum[h];
			for (int d = 0 ; d < HAP_NUMBER ; d ++) {
				dip_probs[d] = (BetaSum[d] / sumHap) * (BetaSum[HAP_NUMBER - d - 1] / sumHap);
				sumDip += dip_probs[d];
			}
		} else {
			TRANSH();
			for (int d = 0 ; d < HAP_NUMBER ; d ++) {
				int prev_h0 = dip_sampled[curr_segment_index - 1];
				int prev_h1 = HAP_NUMBER - dip_sampled[curr_segment_index - 1] - 1;
				dip_probs[d] = (HPROBS[prev_h0 * HAP_NUMBER + d] / sumHProbs) * (HPROBS[prev_h1 * HAP_NUMBER + (HAP_NUMBER - d - 1)] / sumHProbs);
				sumDip += dip_probs[d];
			}
		}
		dip_sampled[curr_segment_index] = rng.sample(dip_probs, sumDip);
		curr_locus += segments[curr_segment_index];
	}
	curr_segment_index = 0;
	curr_segment_locus = 0;
	for (curr_locus = 0  ; curr_locus < C->n_sites ; curr_locus ++) {
		if (HET[curr_locus] >= 0) {
			int idx_h0 = dip_sampled[curr_segment_index];
			int idx_h1 = HAP_NUMBER - dip_sampled[curr_segment_index] - 1;
			bool a0 = ALLELE(HAP_NUMBER - idx_h0 - 1, 2 - HET[curr_locus]);
			bool a1 = ALLELE(HAP_NUMBER - idx_h1 - 1, 2 - HET[curr_locus]);
			H0[curr_locus] = a0;
			H1[curr_locus] = a1;
		}
		curr_segment_locus ++;
		if (curr_segment_locus >= segments[curr_segment_index]) {
			curr_segment_index++;
			curr_segment_locus = 0;
		}
	}
	vrb.bullet("DIP Phasing meth1 (" + stb.str(tac.rel_time()*1.0, 1) + "ms)");tac.clock();
}
