#include <containers/variant_map.h>

variant_map::variant_map() : baseline(0.0) {
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

std::vector < variant * > variant_map::getByPos (int pos) {
	std::vector < variant * > vecS;
	std::pair < std::multimap < int , variant * >::iterator , std::multimap < int , variant * >::iterator > ret = map_pos.equal_range(pos);
	for (std::multimap < int , variant * >::iterator it = ret.first ; it != ret.second ; ++it) vecS.push_back(it->second);
	return vecS;
}

std::vector < variant * > variant_map::getByRef(int pos, std::string & ref, std::string & alt) {
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

int variant_map::setCentiMorgan(std::vector < int > & pos_bp, std::vector < double > & pos_cM) {
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

double variant_map::getCentiMorganPos(gmap_reader & readerGM, int pos_bp)
{
	double base, rate, dist;
	int idx;
	auto lower = std::lower_bound(readerGM.pos_bp.cbegin(), readerGM.pos_bp.cend(), pos_bp);
	if(lower != readerGM.pos_bp.cend()) idx = lower- readerGM.pos_bp.cbegin();
	else return readerGM.pos_cm.back();

	if (idx <= 0) return readerGM.pos_cm[0];

	dist = pos_bp - readerGM.pos_bp[idx-1];
	base = readerGM.pos_cm[idx-1];
	rate = (readerGM.pos_cm[idx] - readerGM.pos_cm[idx-1]) / (readerGM.pos_bp[idx] -  readerGM.pos_bp[idx-1]);

	return base + rate * dist;
}

int variant_map::interpolateCentiMorgan(std::vector < int > & pos_bp, std::vector < double > & pos_cM) {
	int n_interpolated = 0, i_locus = 0;
	double base, rate, dist;
	double mean_rate = (pos_cM.back() - pos_cM[0]) / (pos_bp.back() - pos_bp[0]);

	//Set up first positions to be mean rate
	while (vec_pos[i_locus]->bp < pos_bp[0]) {
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

unsigned int variant_map::length() {
	return vec_pos.back()->bp - vec_pos[0]->bp + 1;
}

void variant_map::setGeneticMap(gmap_reader & readerGM) {
	tac.clock();
	int n_set = setCentiMorgan(readerGM.pos_bp, readerGM.pos_cm);
	int n_interpolated = interpolateCentiMorgan(readerGM.pos_bp, readerGM.pos_cm);
	baseline = vec_pos[0]->cm;
	for (int l = 0 ; l < vec_pos.size() ; l ++) vec_pos[l]->cm -= baseline;
	vrb.bullet("cM interpolation [s=" + stb.str(n_set) + " / i=" + stb.str(n_interpolated) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void variant_map::setGeneticMap() {
	for (int l = 0 ; l < vec_pos.size() ; l ++) vec_pos[l]->cm = vec_pos[l]->bp * 1.0 / 1e6;
	double baseline = vec_pos[0]->cm;
	for (int l = 0 ; l < vec_pos.size() ; l ++) vec_pos[l]->cm -= baseline;
}
