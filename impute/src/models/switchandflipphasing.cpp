#include <models/switchandflipphasing.h>

switchandflipphasing::switchandflipphasing(haplotype_set * _H, conditioning_set * _C) : H(_H), C(_C) {
}

switchandflipphasing::~switchandflipphasing(){}

void switchandflipphasing::sampling(std::vector < bool > & H0, std::vector < bool > & H1) {
	tac.clock();
	UP = std::vector < unsigned char > (H0.size(), 0);
	HAP = std::vector < std::vector < bool > > (2, std::vector < bool > (H0.size(), false));
	std::vector < bool > BUFFER = std::vector < bool > (H0.size(), false);
	HAP[0] = H0;
	HAP[1] = H1;

	//ALLOCATE MEMORY
	backwardProbs = std::vector < std::vector < std::vector < float > > > (2, std::vector < std::vector < float > > (C->n_sites, std::vector < float > (C->n_states, 1.0)));
	forwardProbs = std::vector < std::vector < std::vector < float > > > (2, std::vector < std::vector < float > > (4, std::vector < float > (C->n_states, 0.0)));

	//BACKWARD PASS FOR H0 & H1
	std::vector < float > fact1 = std::vector < float > (2, 0.0);
	std::vector < float > fact2 = std::vector < float > (2, 0.0);
	std::vector < float > betaSumPrevH = std::vector < float > (2, 0.0);
	std::vector < float > betaSumCurrH = std::vector < float > (2, 0.0);
	for (int l = C->n_sites-1 ; l >= 0 ; l --) {
		std::fill(betaSumCurrH.begin(), betaSumCurrH.end(), 0.0);
		if (l == C->n_sites - 1) {
			for (int k = 0 ; k < C->n_states ; k ++) {
				betaSumCurrH[0] += (HAP[0][l]==H->H_opt_var.get(l,C->idxH[k]))?C->ee:C->ed;
				betaSumCurrH[1] += (HAP[1][l]==H->H_opt_var.get(l,C->idxH[k]))?C->ee:C->ed;
			}
		} else {
			fact1[0] = C->nt[l] / betaSumPrevH[0];
			fact1[1] = C->nt[l] / betaSumPrevH[1];
			for (int k = 0 ; k < C->n_states ; k ++) {
				backwardProbs[0][l][k] = C->t[l] + backwardProbs[0][l+1][k] * fact1[0] * ((HAP[0][l+1]==H->H_opt_var.get(l+1,C->idxH[k]))?C->ee:C->ed);
				backwardProbs[1][l][k] = C->t[l] + backwardProbs[1][l+1][k] * fact1[1] * ((HAP[1][l+1]==H->H_opt_var.get(l+1,C->idxH[k]))?C->ee:C->ed);
				betaSumCurrH[0] += ((HAP[0][l]==H->H_opt_var.get(l,C->idxH[k]))?C->ee:C->ed) * backwardProbs[0][l][k];
				betaSumCurrH[1] += ((HAP[1][l]==H->H_opt_var.get(l,C->idxH[k]))?C->ee:C->ed) * backwardProbs[1][l][k];
			}
		}
		betaSumPrevH[0] = betaSumCurrH[0] / C->n_states;
		betaSumPrevH[1] = betaSumCurrH[1] / C->n_states;
	}

	//FORWARD PASS FOR H0 & H1

	int curr_idx = 1, prev_idx = 0;
	int prev_hap0, prev_hap1;
	std::vector < double > phasingD = std::vector < double > (4, 0.0);
	std::vector < double > phasingH = std::vector < double > (8, 0.0);
	std::vector < std::vector < float > > alphaSumCurrH = std::vector < std::vector < float > > (2, std::vector < float > (4, 0.0));
	for (int l = 0 ; l < C->n_sites ; l++) {
		curr_idx = 1 - curr_idx;
		prev_idx = 1 - prev_idx;
		std::fill(alphaSumCurrH[curr_idx].begin(), alphaSumCurrH[curr_idx].end(), 0.0);

		if (H0[l] == H1[l]) {

			//HOM CASE
			if (l == 0) {
				fact1[0] = 1.0 / C->n_states;
				fact1[1] = 1.0 / C->n_states;
				for (int k = 0 ; k < C->n_states ; k ++) {
					forwardProbs[curr_idx][0][k] = fact1[0] * (HAP[0][l]==H->H_opt_var.get(l,C->idxH[k]))?C->ee:C->ed;
					forwardProbs[curr_idx][1][k] = fact1[1] * (HAP[1][l]==H->H_opt_var.get(l,C->idxH[k]))?C->ee:C->ed;
					alphaSumCurrH[curr_idx][0] += forwardProbs[curr_idx][0][k];
					alphaSumCurrH[curr_idx][1] += forwardProbs[curr_idx][1][k];
				}
			} else {
				fact1[0] = C->t[l-1] / C->n_states;
				fact1[1] = C->t[l-1] / C->n_states;
				fact2[0] = C->nt[l-1] / alphaSumCurrH[prev_idx][prev_hap0];
				fact2[1] = C->nt[l-1] / alphaSumCurrH[prev_idx][prev_hap1];
				for (int k = 0 ; k < C->n_states ; k ++) {
					double emit = (HAP[0][l]==H->H_opt_var.get(l,C->idxH[k]))?C->ee:C->ed;
					forwardProbs[curr_idx][0][k] = (forwardProbs[prev_idx][prev_hap0][k]  * fact2[0] + fact1[0]) * emit;
					forwardProbs[curr_idx][1][k] = (forwardProbs[prev_idx][prev_hap1][k]  * fact2[1] + fact1[1]) * emit;
					alphaSumCurrH[curr_idx][0] += forwardProbs[curr_idx][0][k];
					alphaSumCurrH[curr_idx][1] += forwardProbs[curr_idx][1][k];
				}
			}
			prev_hap0 = 0;
			prev_hap1 = 1;
		} else {
			std::fill(phasingD.begin(), phasingD.end(), 0.0);
			std::fill(phasingH.begin(), phasingH.end(), 0.0);

			//HET CASE
			if (l == 0) {
				fact1[0] = 1.0 / C->n_states;
				fact1[1] = 1.0 / C->n_states;
				for (int k = 0 ; k < C->n_states ; k ++) {
					forwardProbs[curr_idx][0][k] = fact1[0] * ((HAP[0][l]==H->H_opt_var.get(l,C->idxH[k]))?C->ee:C->ed);
					forwardProbs[curr_idx][1][k] = fact1[1] * ((HAP[1][l]==H->H_opt_var.get(l,C->idxH[k]))?C->ee:C->ed);
					forwardProbs[curr_idx][2][k] = fact1[0] * ((HAP[1][l]==H->H_opt_var.get(l,C->idxH[k]))?C->ee:C->ed);
					forwardProbs[curr_idx][3][k] = fact1[1] * ((HAP[0][l]==H->H_opt_var.get(l,C->idxH[k]))?C->ee:C->ed);
					for (int c = 0 ; c < 4 ; c++) alphaSumCurrH[curr_idx][c] += forwardProbs[curr_idx][c][k];
				}
			} else {
				fact1[0] = C->t[l-1] / C->n_states;
				fact1[1] = C->t[l-1] / C->n_states;
				fact2[0] = C->nt[l-1] / alphaSumCurrH[prev_idx][prev_hap0];
				fact2[1] = C->nt[l-1] / alphaSumCurrH[prev_idx][prev_hap1];
				for (int k = 0 ; k < C->n_states ; k ++) {
					forwardProbs[curr_idx][0][k] = (forwardProbs[prev_idx][prev_hap0][k] * fact2[0] + fact1[0]) * ((HAP[0][l]==H->H_opt_var.get(l,C->idxH[k]))?C->ee:C->ed);	//f0 + a0
					forwardProbs[curr_idx][1][k] = (forwardProbs[prev_idx][prev_hap1][k] * fact2[1] + fact1[1]) * ((HAP[1][l]==H->H_opt_var.get(l,C->idxH[k]))?C->ee:C->ed);	//f1 + a1
					forwardProbs[curr_idx][2][k] = (forwardProbs[prev_idx][prev_hap0][k] * fact2[0] + fact1[0]) * ((HAP[1][l]==H->H_opt_var.get(l,C->idxH[k]))?C->ee:C->ed);	//f0 + a1
					forwardProbs[curr_idx][3][k] = (forwardProbs[prev_idx][prev_hap1][k] * fact2[1] + fact1[1]) * ((HAP[0][l]==H->H_opt_var.get(l,C->idxH[k]))?C->ee:C->ed);	//f1 + a0
					for (int c = 0 ; c < 4 ; c++) alphaSumCurrH[curr_idx][c] += forwardProbs[curr_idx][c][k];
				}
			}
			//PHASE
			for (int k = 0 ; k < C->n_states ; k ++) {
				phasingH[0] += forwardProbs[curr_idx][0][k] * backwardProbs[0][l][k];	//f0 + a0 + b0
				phasingH[1] += forwardProbs[curr_idx][2][k] * backwardProbs[0][l][k];	//f0 + a1 + b0
				phasingH[2] += forwardProbs[curr_idx][3][k] * backwardProbs[0][l][k];	//f1 + a0 + b0
				phasingH[3] += forwardProbs[curr_idx][1][k] * backwardProbs[0][l][k];	//f1 + a1 + b0
				phasingH[4] += forwardProbs[curr_idx][0][k] * backwardProbs[1][l][k];	//f0 + a0 + b1
				phasingH[5] += forwardProbs[curr_idx][2][k] * backwardProbs[1][l][k];	//f0 + a1 + b1
				phasingH[6] += forwardProbs[curr_idx][3][k] * backwardProbs[1][l][k];	//f1 + a0 + b1
				phasingH[7] += forwardProbs[curr_idx][1][k] * backwardProbs[1][l][k];	//f1 + a1 + b1


			}
			phasingD[0] = phasingH[0] * phasingH[7];	//No flip / No Switch		// f0 + a0 + b0		f1 + a1 + b1	/	f0 + f1
			phasingD[1] = phasingH[1] * phasingH[6];	//Yes flip / No Switch		// f0 + a1 + b0		f1 + a0 + b1	/	f2 + f3
			phasingD[2] = phasingH[2] * phasingH[5];	//No flip / Yes Switch		// f1 + a0 + b0		f0 + a1 + b1	/ 	f3 + f2
			phasingD[3] = phasingH[3] * phasingH[4];	//Yes flip / Yes Switch		// f1 + a1 + b0		f0 + a0 + b1	/	f1 + f0

			double sumP = accumulate(phasingD.begin(), phasingD.end(), 0.0);
			for (int c = 0 ; c < 4 ; c++) phasingD[c] /= sumP;
			//cout << l << " " << stb.str(phasingD[0], 3) << " " << stb.str(phasingD[1], 3) << " " << stb.str(phasingD[2], 3) << " " << stb.str(phasingD[3], 3) << endl;

			//
			UP[l] = rng.sample(phasingD, 1.0);
			/*
			if (UP[l] == 0) { prevForwardProbs[0] = forwardProbs[0]; prevForwardProbs[1] = forwardProbs[1]; prevAlphaSumCurrH[0] =  alphaSumCurrH[0]; prevAlphaSumCurrH[1] =  alphaSumCurrH[1]; }
			if (UP[l] == 1) { prevForwardProbs[0] = forwardProbs[2]; prevForwardProbs[1] = forwardProbs[3]; prevAlphaSumCurrH[0] =  alphaSumCurrH[2]; prevAlphaSumCurrH[1] =  alphaSumCurrH[3]; }
			if (UP[l] == 2) { prevForwardProbs[0] = forwardProbs[3]; prevForwardProbs[1] = forwardProbs[2]; prevAlphaSumCurrH[0] =  alphaSumCurrH[3]; prevAlphaSumCurrH[1] =  alphaSumCurrH[2]; }
			if (UP[l] == 3) { prevForwardProbs[0] = forwardProbs[1]; prevForwardProbs[1] = forwardProbs[0]; prevAlphaSumCurrH[0] =  alphaSumCurrH[1]; prevAlphaSumCurrH[1] =  alphaSumCurrH[0]; }
			*/
			switch (UP[l]) {
			case 0: prev_hap0 = 0; prev_hap1 = 1; break;
			case 1: prev_hap0 = 2; prev_hap1 = 3; break;
			case 2: prev_hap0 = 3; prev_hap1 = 2; break;
			case 3: prev_hap0 = 1; prev_hap1 = 0; break;
			}
		}
	}

	bool switched = false;
	for (int l = C->n_sites - 1 ; l >= 0 ; l --) {
		if (H0[l] != H1[l]) {
			if (UP[l] == 1 || UP[l] == 3) {
				H0[l] = !H0[l];
				H1[l] = !H1[l];
			}
			if (switched) {
				H0[l] = !H0[l] ;
				H1[l] = !H1[l];
			}
			if (UP[l] >= 2) switched = !switched;
		}
	}

	vrb.bullet("DIP Phasing meth2 (" + stb.str(tac.rel_time()*1.0, 1) + "ms)");tac.clock();
}
