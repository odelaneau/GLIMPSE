////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2018 Olivier Delaneau, University of Lausanne
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////
#include <containers/haplotype_set.h>

const uint32_t haps_to_int[4] = {0x0,0x1,0x2,0x3};

haplotype_set::haplotype_set() :
	n_site(0),
	n_hap(0),
	n_targ(0),
	n_ref(0),
	n_haps_pbwt(0),
	refonly_pbwt(false),
	pbwt_built(false),
	pbwt_modulo(1),
	pbwt_depth(1),
	n_iterations_main(0)
{
}

haplotype_set::~haplotype_set()
{
}

void haplotype_set::updateHaplotypes(genotype_set & G) {
	tac.clock();
	for (int i = 0 ; i < G.n_ind ; i ++) {
		for (int v = 0 ; v < n_site ; v ++) {
			bool a0 = G.vecG[i]->H0[v];
			bool a1 = G.vecG[i]->H1[v];
			H_opt_var.set(v, n_ref+2*i+0, a0);
			H_opt_var.set(v, n_ref+2*i+1, a1);
		}
	}
	vrb.bullet("HAP update (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::initPositionalBurrowWheelerTransform(int _pbwt_depth, int _pbwt_modulo) {
	n_haps_pbwt = n_ref;

	pbwt_depth = _pbwt_depth;
	pbwt_modulo = _pbwt_modulo;
	pbwt_arrays = std::vector < std::vector < int > > ((n_site/pbwt_modulo) + 1, std::vector < int > (n_haps_pbwt, 0));
	pbwt_indexes = std::vector < std::vector < int > > ((n_site) + 1, std::vector < int > (n_haps_pbwt+1, 0));

}

void haplotype_set::updatePositionalBurrowWheelerTransform() {
	if (!pbwt_built) {
		tac.clock();
		std::vector < int > A = std::vector < int >(n_haps_pbwt, 0);
		std::vector < int > B = std::vector < int >(n_haps_pbwt, 0);
		for (int l = 0, idx_store = 0 ; l < n_site ; l ++) {
			int u = 0, v = 0;
			for (int h = 0 ; h < n_haps_pbwt ; h ++) {
				pbwt_indexes[l][h] = v;
				int alookup = l?A[h]:h;
				if (!H_opt_var.get(l,alookup)) A[u++] = alookup;
				else B[v++] = alookup;
			}
			std::copy(B.begin(), B.begin()+v, A.begin()+u);
			pbwt_indexes[l][n_haps_pbwt] = v;
			if ((l%pbwt_modulo) == 0) {
				std::copy(A.begin(), A.end(), pbwt_arrays[idx_store].begin());
				idx_store++;
			}
		}
		pbwt_built = true;
		vrb.bullet("PBWT building (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	}
}

void haplotype_set::selectRandomRefOnly(const int K, conditioning_set * C) {
	std::vector < int >& idxH = C->idxH;
	idxH.clear();
	for (int h = 0 ; h < n_ref ; h ++) idxH.push_back(h);
	if (n_ref > K) {
		random_shuffle(idxH.begin(), idxH.end());
		idxH.erase(idxH.begin() + K, idxH.end());
		sort(idxH.begin(), idxH.end());
	}

	//C->clear();
	C->n_states = idxH.size();
	/*
	for (int l = 0 ; l < n_site ; l ++) {
		unsigned int ac = 0;
		C->Hmono.push_back(H_opt_var.get(l,idxH[0]));
		for (int k = 0 ; k < idxH.size() ; k ++) ac += H_opt_var.get(l,idxH[k]);
		if (ac > 0 && ac < idxH.size()) {
			C->Vpoly.push_back(l);
			for (int k = 0 ; k < idxH.size() ; k ++) C->Hpoly.push_back(H_opt_var.get(l,idxH[k]));
		} else C->Vmono.push_back(l);
	}
	*/
	//C->n_sites = C->Vpoly.size();
	//C->updateTransitions();
}

void haplotype_set::selectPositionalBurrowWheelerTransform(int ind, conditioning_set * C) {
	std::vector < int >& idxH = C->idxH;
	idxH.clear();

	selectFmIndexRefOnly(ind, idxH);

	//C->clear();
	C->n_states = idxH.size();
	//std::vector < int > freq = std::vector < int > (10, 0);
	/*
	for (int l = 0 ; l < n_site ; l ++) {
		unsigned int ac = H_opt_var.get(l, 2*ind + n_ref + 0) + H_opt_var.get(l, 2*ind + n_ref + 1);
		unsigned int ac2 = H_opt_var.get(l, 2*ind + n_ref + 0) + H_opt_var.get(l, 2*ind + n_ref + 1);
		C->Hmono.push_back(H_opt_var.get(l,idxH[0]));
		for (int k = 0 ; k < idxH.size() ; k ++) ac2 += H_opt_var.get(l,idxH[k]);
		ac += ac2;
		//freq[(ac2>4)?4:ac2]++;
		if (ac > 0 && ac < (idxH.size()+2)) {
			C->Vpoly.push_back(l);
			for (int k = 0 ; k < idxH.size() ; k ++) C->Hpoly.push_back(H_opt_var.get(l,idxH[k]));
		} else C->Vmono.push_back(l);
	}
	//cout << freq[0] << " " << freq[1] << " " << freq[2] << " " << freq[3] << " " << freq[4] << endl;
	C->n_sites = C->Vpoly.size();
	C->updateTransitions();
	*/
}

void haplotype_set::selectFmIndexRefOnly(const int ind, std::vector< int >& idxH) {
	std::vector < bool > hapFlag (n_haps_pbwt, false);
	std::vector < int > loc_h (2);

	const int max_idx = (int) n_haps_pbwt-1;

	for (int l = 0, idx_store = 0 ; l < n_site ; ++l) {
		int c_b = getC_bounded(l); // c_b current marker
		for (int p = 0 ; p < 2 ; p ++) {
			int h = 2*ind + p + n_ref;
            bool a = H_opt_var.get(l,h);
            int& i = loc_h[p];
            //update pos
    		if (l==0)
    			i = a ? (c_b + n_haps_pbwt) / 2 : c_b /2;
    		else
    			i = a ? std::min(c_b + pbwt_indexes[l][i],max_idx) : i - pbwt_indexes[l][i];// fm index prev marker

            //select pos when pbwt modulo
    		if ((l%pbwt_modulo) == 0) {
    			//Backward - up
    			bool ac;
                int o = 0, c = 0, hc = 0;
                for (;;) {
                	if (i-o >= 0) hc = pbwt_arrays[idx_store][i-o];
                	else break;
                	ac = H_opt_var.get(l,hc);
                	if (ac != a) break;
                	hapFlag[hc] = (hc/2 != h/2);
                	c += (hc/2 != h/2);
                	if (c >= pbwt_depth) break;
                	o++;
                }

                //Forward - down
                o = 1; c = 0;
                for (;;) {
                	if (i+o < n_haps_pbwt) hc = pbwt_arrays[idx_store][i+o];
                	else break;
                	ac = H_opt_var.get(l,hc);
                	if (ac != a) break;
                	hapFlag[hc] = (hc/2 != h/2);
                	c += (hc/2 != h/2);
                	if (c >= pbwt_depth) break;
                	o++;
                }
    			if (p) idx_store++;
    		}
		}
	}
	int c0 = 0;
	int c1 = 0;
	int c01 = 0;
	for (int h = 0 ; h < n_haps_pbwt ; h ++) if (hapFlag[h]) idxH.push_back(h); //hash?
}

void haplotype_set::selectStandardFull(const int ind, std::vector< int >& idxH) {
	std::vector < bool > hapFlag = std::vector < bool > (n_haps_pbwt, false);
	for (int l = 0, idx_store = 0 ; l < n_site ; l += pbwt_modulo, idx_store++) {
		for (int p = 0 ; p < 2 ; p ++) {
			int h = 2*ind + p + n_ref;
			int i = pbwt_indexes[idx_store][h];

			//Backward - up
			bool ac;
            int o = 1, c = 0, hc = 0;
            bool a = H_opt_var.get(l,h);
            for (;;) {
            	if (i-o >= 0) hc = pbwt_arrays[idx_store][i-o];
            	else break;
            	ac = H_opt_var.get(l,hc);
            	if (ac != a) break;
            	hapFlag[hc] = (hc/2 != h/2);
            	c += (hc/2 != h/2);
            	if (c >= pbwt_depth) break;
            	o++;
            }

            //Forward - down
            o = 1; c = 0;
            for (;;) {
            	if (i+o < n_haps_pbwt) hc = pbwt_arrays[idx_store][i+o];
            	else break;
            	ac = H_opt_var.get(l,hc);
            	if (ac != a) break;
            	hapFlag[hc] = (hc/2 != h/2);
            	c += (hc/2 != h/2);
            	if (c >= pbwt_depth) break;
            	o++;
            }
		}
	}
	for (int h = 0 ; h < n_haps_pbwt ; h ++) if (hapFlag[h]) idxH.push_back(h);
}

int haplotype_set::getC_bounded(const int l) const {
	int c = n_ref - pbwt_indexes[l][n_ref];
	if (c >= n_ref) --c;
	return c;
}
