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

#include <containers/variant_map.h>

variant_map::variant_map() {
}

variant_map::~variant_map() {
	for (int s = 0 ; s < vec_pos.size() ; s++) delete vec_pos[s];
	vec_pos.clear();
	map_pos.clear();
}

int variant_map::size() {
	return vec_pos.size();
}

variant * variant_map::getByIndex (int i) {
	return vec_pos[i];
}

vector < variant * > variant_map::getByPos (int pos) {
	vector < variant * > vecS;
	pair < multimap < int , variant * >::iterator , multimap < int , variant * >::iterator > ret = map_pos.equal_range(pos);
	for (multimap < int , variant * >::iterator it = ret.first ; it != ret.second ; ++it) vecS.push_back(it->second);
	return vecS;
}

vector < variant * > variant_map::getByRef(int pos, string & ref, string & alt) {
	vector < variant * > vecS = vector < variant * >();
	pair < multimap < int , variant * >::iterator , multimap < int , variant * >::iterator > ret = map_pos.equal_range(pos);
	for (multimap < int , variant * >::iterator it = ret.first ; it != ret.second ; ++it) {
		if (it->second->ref == ref && it->second->alt == alt) vecS.push_back(it->second);
	}
	return vecS;
}

void variant_map::push(variant * v) {
	vec_pos.push_back(v);
	map_pos.insert(pair < int , variant * > (v->bp, v));
}

int variant_map::setCentiMorgan(vector < int > & pos_bp, vector < double > & pos_cM) {
	int cpt = 0;
	for (int l = 0 ; l < pos_cM.size() ; l ++) {
		vector  < variant * > vecS = getByPos(pos_bp[l]);
		for (int si = 0 ; si < vecS.size() ; si ++) {
			vecS[si]->cm = pos_cM[l];
			cpt++;
		}
	}
	return cpt;
}

int variant_map::interpolateCentiMorgan(vector < int > & pos_bp, vector < double > & pos_cM) {
	int cpt = 0;
	double mean_rate = (pos_cM.back() - pos_cM[0]) / (pos_bp.back() - pos_bp[0]);
	for (int s = 0 ; s < vec_pos.size() ; s ++) {
		if (vec_pos[s]->cm < 0) {
			if (vec_pos[s]->bp < pos_bp[0]) vec_pos[s]->cm = pos_cM[0] - mean_rate * (pos_bp[0] - vec_pos[s]->bp);
			else if (vec_pos[s]->bp > pos_bp.back()) vec_pos[s]->cm = pos_cM.back() + mean_rate * (vec_pos[s]->bp - pos_bp.back());
			else {
				int index_from, index_to;
				for ( index_from = 0 ; index_from < pos_bp.size() && pos_bp[index_from] < vec_pos[s]->bp; ) index_from ++ ;
				for ( index_to = pos_cM.size() - 1 ; index_to >= 0 && pos_bp[index_to] > vec_pos[s]->bp; ) index_to --;
				index_from--;
				index_to++;
				vec_pos[s]->cm = pos_cM[index_from] + (vec_pos[s]->bp - pos_bp[index_from]) * (pos_cM[index_to] - pos_cM[index_from]) / (pos_bp[index_to] - pos_bp[index_from]);
				}
			cpt++;
		}
		vrb.progress("  * cM interpolation", (s+1)*1.0/vec_pos.size());

	}
	if (vec_pos[0]->cm < 0) {
		double socle = -1.0 * vec_pos[0]->cm;
		for (int s = 0 ; s < vec_pos.size() ; s ++) vec_pos[s]->cm += socle;
	}
	return cpt;
}

unsigned int variant_map::length() {
	return vec_pos.back()->bp - vec_pos[0]->bp + 1;
}

void variant_map::setGeneticMap(gmap_reader & readerGM) {
	tac.clock();
	int n_set = setCentiMorgan(readerGM.pos_bp, readerGM.pos_cm);
	int n_interpolated = interpolateCentiMorgan(readerGM.pos_bp, readerGM.pos_cm);
	vrb.bullet("cM interpolation [s=" + stb.str(n_set) + " / i=" + stb.str(n_interpolated) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
