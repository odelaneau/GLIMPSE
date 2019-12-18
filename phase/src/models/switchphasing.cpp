#include <models/switchphasing.h>

switchphasing::switchphasing(conditioning_set * _C) : C(_C) {
}

void switchphasing::sampling(vector < bool > & H0, vector < bool > & H1) {
	UP = vector < unsigned char > (H0.size(), 0);
	HAP = vector < vector < bool > > (2, vector < bool > (H0.size(), false));
	HAP[0] = H0;
	HAP[1] = H1;

	//ALLOCATE MEMORY
	backwardProbs = vector < vector < vector < float > > > (2, vector < vector < float > > (C->n_sites, vector < float > (C->n_states, 1.0)));
	forwardProbs = vector < vector < float > > (2, vector < float > (C->n_states, 0.0));
	prevForwardProbs = vector < vector < float > > (2, vector < float > (C->n_states, 0.0));

	//BACKWARD PASS FOR H0 & H1
	vector < float > fact1 = vector < float > (2, 0.0);
	vector < float > fact2 = vector < float > (2, 0.0);
	vector < float > betaSumPrevH = vector < float > (2, 0.0);
	vector < float > betaSumCurrH = vector < float > (2, 0.0);
	for (int l = C->n_sites-1 ; l >= 0 ; l --) {
		fill(betaSumCurrH.begin(), betaSumCurrH.end(), 0.0);
		if (l == C->n_sites - 1) {
			for (int k = 0 ; k < C->n_states ; k ++) {
				betaSumCurrH[0] += (HAP[0][C->Vpoly[l]]==C->Hpoly[l*C->n_states+k])?C->ee:C->ed;
				betaSumCurrH[1] += (HAP[1][C->Vpoly[l]]==C->Hpoly[l*C->n_states+k])?C->ee:C->ed;
			}
		} else {
			fact1[0] = C->nt[l] / betaSumPrevH[0];
			fact1[1] = C->nt[l] / betaSumPrevH[1];
			for (int k = 0 ; k < C->n_states ; k ++) {
				backwardProbs[0][l][k] = C->t[l] + backwardProbs[0][l+1][k] * fact1[0] * ((HAP[0][C->Vpoly[l+1]]==C->Hpoly[(l+1)*C->n_states+k])?C->ee:C->ed);
				backwardProbs[1][l][k] = C->t[l] + backwardProbs[1][l+1][k] * fact1[1] * ((HAP[1][C->Vpoly[l+1]]==C->Hpoly[(l+1)*C->n_states+k])?C->ee:C->ed);
				betaSumCurrH[0] += ((HAP[0][C->Vpoly[l]]==C->Hpoly[l*C->n_states+k])?C->ee:C->ed) * backwardProbs[0][l][k];
				betaSumCurrH[1] += ((HAP[1][C->Vpoly[l]]==C->Hpoly[l*C->n_states+k])?C->ee:C->ed) * backwardProbs[1][l][k];
			}
		}
		betaSumPrevH[0] = betaSumCurrH[0] / C->n_states;
		betaSumPrevH[1] = betaSumCurrH[1] / C->n_states;
		//cout << l << " " << " " << betaSumPrevH[0] << " " << betaSumPrevH[1] << endl;
	}

	//FORWARD PASS FOR H0 & H1
	bool prevSW = false;
	int idx0 = 0, idx1 = 1;
	vector < double > phasingD = vector < double > (4, 0.0);
	vector < double > phasingH = vector < double > (8, 0.0);

	vector < float > alphaSumCurrH = vector < float > (4, 0.0);
	vector < float > prevAlphaSumCurrH = vector < float > (4, 0.0);
	for (int l = 0 ; l < C->n_sites ; l++) {
		fill(alphaSumCurrH.begin(), alphaSumCurrH.end(), 0.0);
		fill(phasingD.begin(), phasingD.end(), 0.0);
		fill(phasingH.begin(), phasingH.end(), 0.0);

		if (H0[C->Vpoly[l]] == H1[C->Vpoly[l]]) {
			//if (prevForwardProbs[0] == prevForwardProbs[1]) cerr<< l << " HOM MATCH" << endl;
			//else  cerr << l << " HOM UNMATCH" << endl;

			//HOM CASE
			if (l == 0) {
				fact1[0] = 1.0 / C->n_states;
				fact1[1] = 1.0 / C->n_states;
				for (int k = 0 ; k < C->n_states ; k ++) {
					forwardProbs[0][k] = fact1[0] * (HAP[0][C->Vpoly[l]]==C->Hpoly[l*C->n_states+k])?C->ee:C->ed;
					forwardProbs[1][k] = fact1[1] * (HAP[1][C->Vpoly[l]]==C->Hpoly[l*C->n_states+k])?C->ee:C->ed;
					alphaSumCurrH[0] += forwardProbs[0][k];
					alphaSumCurrH[1] += forwardProbs[1][k];
				}
			} else {
				fact1[0] = C->t[l-1] / C->n_states;
				fact1[1] = C->t[l-1] / C->n_states;
				fact2[0] = C->nt[l-1] / prevAlphaSumCurrH[0];
				fact2[1] = C->nt[l-1] / prevAlphaSumCurrH[1];
				//cout << fact2[0] << " " << fact2[1] << " " << C->ee << " " << C->ed << endl;
				for (int k = 0 ; k < C->n_states ; k ++) {
					double emit = (HAP[0][C->Vpoly[l]]==C->Hpoly[l*C->n_states+k])?C->ee:C->ed;
					forwardProbs[0][k] = (prevForwardProbs[0][k]  * fact2[0] + fact1[0]) * emit;
					forwardProbs[1][k] = (prevForwardProbs[1][k]  * fact2[1] + fact1[1]) * emit;
					alphaSumCurrH[0] += forwardProbs[0][k];
					alphaSumCurrH[1] += forwardProbs[1][k];
				}
			}

			//if (forwardProbs[0] == forwardProbs[1]) cerr<< l << " HOM2 MATCH" << endl;
			//else  cerr << l << " HOM2 UNMATCH" << endl;

/*
			//PHASE
			double sumH = 0.0;
			for (int k = 0 ; k < C->n_states ; k ++) {
				phasingH[0] += forwardProbs[0][k] /alphaSumCurrH[0] * backwardProbs[0][l][k];	//f0 + a0 + b0
				phasingH[3] += forwardProbs[1][k] /alphaSumCurrH[1] * backwardProbs[0][l][k];	//f1 + a0 + b0
				phasingH[4] += forwardProbs[0][k] /alphaSumCurrH[0] * backwardProbs[1][l][k];	//f0 + a0 + b1
				phasingH[7] += forwardProbs[1][k] /alphaSumCurrH[1] * backwardProbs[1][l][k];	//f1 + a0 + b1
			}
			phasingD[0] = phasingH[0] * phasingH[7];
			phasingD[3] = phasingH[3] * phasingH[4];
			*/
			prevForwardProbs[0] = forwardProbs[0];
			prevForwardProbs[1] = forwardProbs[1];
			prevAlphaSumCurrH[0] =  alphaSumCurrH[0];
			prevAlphaSumCurrH[1] =  alphaSumCurrH[1];

		} else {
			//if (prevForwardProbs[0] == prevForwardProbs[1]) cerr<< l << " HET MATCH" << endl;
			//else  cerr << l << " HET UNMATCH" << endl;

			//HET CASE
			if (l == 0) {
				fact1[0] = 1.0 / C->n_states;
				fact1[1] = 1.0 / C->n_states;
				for (int k = 0 ; k < C->n_states ; k ++) {
					forwardProbs[0][k] = fact1[0] * ((HAP[0][C->Vpoly[l]]==C->Hpoly[l*C->n_states+k])?C->ee:C->ed);
					forwardProbs[1][k] = fact1[1] * ((HAP[1][C->Vpoly[l]]==C->Hpoly[l*C->n_states+k])?C->ee:C->ed);
					for (int c = 0 ; c < 2 ; c++) alphaSumCurrH[c] += forwardProbs[c][k];
					//forwardProbs[2][k] = fact1[0] * ((HAP[1][C->Vpoly[l]]==C->Hpoly[l*C->n_states+k])?C->ee:C->ed);	//t2
					//forwardProbs[3][k] = fact1[1] * ((HAP[0][C->Vpoly[l]]==C->Hpoly[l*C->n_states+k])?C->ee:C->ed);	//t3
					//for (int c = 0 ; c < 4 ; c++) alphaSumCurrH[c] += forwardProbs[c][k];
				}
			} else {
				fact1[0] = C->t[l-1] / C->n_states;
				fact1[1] = C->t[l-1] / C->n_states;
				fact2[0] = C->nt[l-1] / prevAlphaSumCurrH[0];
				fact2[1] = C->nt[l-1] / prevAlphaSumCurrH[1];
				for (int k = 0 ; k < C->n_states ; k ++) {

					forwardProbs[0][k] = (prevForwardProbs[0][k]  * fact2[0] + fact1[0]) * ((HAP[0][C->Vpoly[l]]==C->Hpoly[l*C->n_states+k])?C->ee:C->ed);	//f0 + a0
					forwardProbs[1][k] = (prevForwardProbs[1][k]  * fact2[1] + fact1[1]) * ((HAP[1][C->Vpoly[l]]==C->Hpoly[l*C->n_states+k])?C->ee:C->ed);	//f1 + a1
					for (int c = 0 ; c < 2 ; c++) alphaSumCurrH[c] += forwardProbs[c][k];
					//forwardProbs[2][k] = (prevForwardProbs[0][k]  * fact2[0] + fact1[0]) * ((HAP[1][C->Vpoly[l]]==C->Hpoly[l*C->n_states+k])?C->ee:C->ed);	//f0 + a1
					//forwardProbs[3][k] = (prevForwardProbs[1][k]  * fact2[1] + fact1[1]) * ((HAP[0][C->Vpoly[l]]==C->Hpoly[l*C->n_states+k])?C->ee:C->ed);	//f1 + a0
					//for (int c = 0 ; c < 4 ; c++) alphaSumCurrH[c] += forwardProbs[c][k];
				}
			}
			//PHASE
			for (int k = 0 ; k < C->n_states ; k ++) {
				phasingH[0] += forwardProbs[0][k] /alphaSumCurrH[0] * backwardProbs[0][l][k];	//f0 + a0 + b0
				//phasingH[1] += forwardProbs[2][k] /alphaSumCurrH[2] * backwardProbs[0][l][k];	//f0 + a1 + b0
				//phasingH[2] += forwardProbs[3][k] /alphaSumCurrH[3] * backwardProbs[0][l][k];	//f1 + a0 + b0
				phasingH[3] += forwardProbs[1][k] /alphaSumCurrH[1] * backwardProbs[0][l][k];	//f1 + a1 + b0

				phasingH[4] += forwardProbs[0][k] /alphaSumCurrH[0] * backwardProbs[1][l][k];	//f0 + a0 + b1
				//phasingH[5] += forwardProbs[2][k] /alphaSumCurrH[2] * backwardProbs[1][l][k];	//f0 + a1 + b1
				//phasingH[6] += forwardProbs[3][k] /alphaSumCurrH[3] * backwardProbs[1][l][k];	//f1 + a0 + b1
				phasingH[7] += forwardProbs[1][k] /alphaSumCurrH[1] * backwardProbs[1][l][k];	//f1 + a1 + b1
			}
			phasingD[0] = phasingH[0] * phasingH[7];	//No flip / No Switch
			//phasingD[1] = phasingH[1] * phasingH[6];	//Yes flip / No Switch
			//phasingD[2] = phasingH[2] * phasingH[5];	//No flip / Yes Switch
			phasingD[3] = phasingH[3] * phasingH[4];	//Yes flip / Yes Switch

			double sumP = accumulate(phasingD.begin(), phasingD.end(), 0.0);
			for (int c = 0 ; c < 4 ; c++) phasingD[c] /= sumP;

			//cout << l << " " << (H0[C->Vpoly[l]] == H1[C->Vpoly[l]])  << " " << stb.str(phasingD[0], 3) << " " << stb.str(phasingD[1], 3) << " " << stb.str(phasingD[2], 3) << " " << stb.str(phasingD[3], 3) << endl;
			//cout << l << " " << stb.str(phasingD[0], 3) << " " << stb.str(phasingD[1], 3) << " " << stb.str(phasingD[2], 3) << " " << stb.str(phasingD[3], 3) << endl;

			//
			UP[l] = rng.sample(phasingD, 1.0);

			/*
			switch (UP[l]) {
			case 0:	idx0 = 0; idx1 = 1; break;
			case 1: idx0 = 2; idx1 = 3; break;
			case 2: idx0 = 1; idx1 = 0; break;
			case 3:	idx0 = 3; idx1 = 2; break;
			}
			*/
			if (UP[l] == 0) { prevForwardProbs[0] = forwardProbs[0]; prevForwardProbs[1] = forwardProbs[1]; prevAlphaSumCurrH[0] =  alphaSumCurrH[0]; prevAlphaSumCurrH[1] =  alphaSumCurrH[1]; }
			//if (UP[l] == 1) { prevForwardProbs[0] = forwardProbs[2]; prevForwardProbs[1] = forwardProbs[3]; prevAlphaSumCurrH[0] =  alphaSumCurrH[2]; prevAlphaSumCurrH[1] =  alphaSumCurrH[3]; }
			//if (UP[l] == 2) { prevForwardProbs[0] = forwardProbs[3]; prevForwardProbs[1] = forwardProbs[2]; prevAlphaSumCurrH[0] =  alphaSumCurrH[3]; prevAlphaSumCurrH[1] =  alphaSumCurrH[2]; }
			if (UP[l] == 3) { prevForwardProbs[0] = forwardProbs[1]; prevForwardProbs[1] = forwardProbs[0]; prevAlphaSumCurrH[0] =  alphaSumCurrH[1]; prevAlphaSumCurrH[1] =  alphaSumCurrH[0]; }

			//if (prevForwardProbs[0] == prevForwardProbs[1]) cerr<< l << " BOT MATCH" << endl;
			//else  cerr << l << " BOT UNMATCH" << endl;


		}
		//cout << l << " " << alphaSumCurrH[0] << " " << alphaSumCurrH[1] << " " << alphaSumCurrH[2] << " " << alphaSumCurrH[3] << endl;

		//double sumP = phasing[0] + phasing[1] + phasing[2] + phasing[3];
		//for (int c = 0 ; c < 4 ; c++) phasing[c] /= sumP;
		//cout << l << " " << (H0[C->Vpoly[l]] == H1[C->Vpoly[l]]) ;
		//for (int c  = 0 ; c < 8 ; c++ ) cout << " " << stb.str(phasingD[c], 3);
		//cout << endl;
	}

	//UPDATING HAPLOTYPES
	bool switched = false;
	for (int l = 0 ; l < C->n_sites ; l++) {
		if (H0[C->Vpoly[l]] != H1[C->Vpoly[l]]) {
			if (switched) {
				H0[C->Vpoly[l]] = HAP[1][C->Vpoly[l]];
				H1[C->Vpoly[l]] = HAP[0][C->Vpoly[l]];
			}
			if (UP[l] == 3) switched = !switched;
		}
	}

	/*
	bool switched = false;
	bool flipped = false;
	for (int l = C->n_sites-1 ; l >= 0 ; l --) {
		//cerr << l << " " << (int)UP[l] << endl;
		if (H0[C->Vpoly[l]] != H1[C->Vpoly[l]]) {
			flipped = (switched && (UP[l]==0||UP[l]==2)) || (!switched && (UP[l]==1||UP[l]==3));
			if (UP[l]>=2) switched = !switched;
			if (flipped) {
				H0[C->Vpoly[l]] = HAP[1][C->Vpoly[l]];
				H1[C->Vpoly[l]] = HAP[0][C->Vpoly[l]];
			}
		}
	}
	*/
}
