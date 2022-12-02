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

#include <containers/variant_map.h>

variant_map::variant_map() :
	input_start(0),
	input_stop(0),
	output_start(0),
	output_stop(0)
{
}

variant_map::~variant_map() {
	for (int s = 0 ; s < vec_pos.size() ; s++) delete vec_pos[s];
	vec_pos.clear();
	map_pos.clear();
}

std::size_t variant_map::size() const {
	return vec_pos.size();
}

variant * variant_map::getByIndex (const int i) {
	return vec_pos[i];
}

std::vector < variant * > variant_map::getByPos (const int pos) {
	std::vector < variant * > vecS;
	std::pair < std::multimap < int , variant * >::iterator , std::multimap < int , variant * >::iterator > ret = map_pos.equal_range(pos);
	for (std::multimap < int , variant * >::iterator it = ret.first ; it != ret.second ; ++it) vecS.push_back(it->second);
	return vecS;
}

std::vector < variant * > variant_map::getByRef(const int pos, const std::string & ref, const std::string & alt) {
	std::vector < variant * > vecS = std::vector < variant * >();
	std::pair < std::multimap < int , variant * >::iterator , std::multimap < int , variant * >::iterator > ret = map_pos.equal_range(pos);
	for (std::multimap < int , variant * >::iterator it = ret.first ; it != ret.second ; ++it) {
		if (it->second->ref == ref && it->second->alt == alt) vecS.push_back(it->second);
	}
	return vecS;
}

void variant_map::push(variant * v) {
	vec_pos.push_back(v);
	map_pos.insert(std::pair < int , variant * > (v->bp, v));
}

int variant_map::setCentiMorgan(const std::vector < int > & pos_bp, const std::vector < double > & pos_cM) {
	int cpt = 0;
	for (int l = 0 ; l < pos_cM.size() ; l ++) {
		std::vector  < variant * > vecS = getByPos(pos_bp[l]);
		for (int si = 0 ; si < vecS.size() ; si ++) {
			vecS[si]->cm = pos_cM[l];
			cpt++;
		}
	}
	return cpt;
}

int variant_map::interpolateCentiMorgan(const std::vector < int > & pos_bp, const std::vector < double > & pos_cM) {
	float prog_step = 1.0/vec_pos.size();
	float prog_bar = 0.0;
	int n_interpolated = 0, i_locus = 0;
	double base, rate, dist;
	double mean_rate = (pos_cM.back() - pos_cM[0]) / (pos_bp.back() - pos_bp[0]);

	//Set up first positions to be mean rate
	while (i_locus<vec_pos.size() && vec_pos[i_locus]->bp < pos_bp[0]) {
		base = pos_cM[0];
		dist = (pos_bp[0] - vec_pos[i_locus]->bp);
		vec_pos[i_locus]->cm = base - mean_rate * dist;
		n_interpolated ++;
		i_locus ++;
	}

	//Set up middle positions using interpolation
	int closest_pos = 1;
	for (; i_locus < vec_pos.size() ; ) {
		if (vec_pos[i_locus]->cm == -1) {

			//Find suitable interpolation interval
			while (vec_pos[i_locus]->bp > pos_bp[closest_pos] && closest_pos < pos_bp.size()) closest_pos++;

			//Interpolate
			if (closest_pos < pos_bp.size()) {
				assert(vec_pos[i_locus]->bp < pos_bp[closest_pos]);
				assert(vec_pos[i_locus]->bp > pos_bp[closest_pos-1]);
				base = pos_cM[closest_pos-1];
				rate = (pos_cM[closest_pos] - pos_cM[closest_pos-1]) / (pos_bp[closest_pos] - pos_bp[closest_pos-1]);
				dist = (vec_pos[i_locus]->bp - pos_bp[closest_pos-1]);
				vec_pos[i_locus]->cm = base + rate * dist;
				n_interpolated ++;
				i_locus ++;
			} else break;
		} else i_locus ++;
	}

	//Set up last positions to be mean rate
	while (i_locus < vec_pos.size()) {
		base = pos_cM.back();
		dist = (vec_pos[i_locus]->bp - pos_bp.back());
		vec_pos[i_locus]->cm = base + mean_rate * dist;
		n_interpolated ++;
		i_locus ++;
	}
	return n_interpolated;
}

unsigned int variant_map::length() const {
	return vec_pos.back()->bp - vec_pos[0]->bp + 1;
}

double variant_map::lengthcM() const {
	return vec_pos.back()->cm - vec_pos[0]->cm;
}


void variant_map::setGeneticMap(const gmap_reader & readerGM) {
	tac.clock();
	if (vec_pos.size() == 0) vrb.error("No variant in common between reference and target panel. This can indicate a problem in the input files or during the parsing.");
	int n_set = setCentiMorgan(readerGM.pos_bp, readerGM.pos_cm);
	int n_interpolated = interpolateCentiMorgan(readerGM.pos_bp, readerGM.pos_cm);
	double baseline = vec_pos[0]->cm;
	for (int l = 0 ; l < vec_pos.size() ; l ++) vec_pos[l]->cm -= baseline;
	vrb.bullet("cM interpolation [s=" + stb.str(n_set) + " / i=" + stb.str(n_interpolated) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void variant_map::setGeneticMap() {
	tac.clock();
	if (vec_pos.size() == 0) vrb.error("No variant in common between reference and target panel. This can indicate a problem in the input files or during the parsing.");
	for (int l = 0 ; l < vec_pos.size() ; l ++) vec_pos[l]->cm = vec_pos[l]->bp * 1.0 / 1e6;
	double baseline = vec_pos[0]->cm;
	for (int l = 0 ; l < vec_pos.size() ; l ++) vec_pos[l]->cm -= baseline;
	vrb.bullet("cM constant [1cM=1Mb] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
