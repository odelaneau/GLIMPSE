#include <chunker/chunker_header.h>

std::vector < int > chunker::getByPos (const int pos, const std::multimap < int, int >& map_positions) {
	std::vector < int > vecS;
	std::pair < std::multimap < int , int >::const_iterator , std::multimap < int , int >::const_iterator > ret = map_positions.equal_range(pos);
	for (std::multimap < int , int >::const_iterator it = ret.first ; it != ret.second ; ++it) vecS.push_back(it->second);
	return vecS;
}

int chunker::setCentiMorgan(const std::vector < int > & pos_bp, const std::vector < double > & pos_cM, const std::multimap < int, int >& map_positions, std::vector<float>& positions_cm) {
	int cpt = 0;
	for (int l = 0 ; l < pos_cM.size() ; l ++) {
		std::vector  < int > vecS = getByPos(pos_bp[l], map_positions);
		for (int si = 0 ; si < vecS.size() ; si ++) {
			positions_cm[vecS[si]] = pos_cM[l];
			cpt++;
		}
	}
	return cpt;
}

int chunker::interpolateCentiMorgan(const std::vector < int > & pos_bp, const std::vector < double > & pos_cM, const std::multimap < int, int >& map_positions, const std::vector<int>& positions_mb, std::vector<float>& positions_cm) {
	float prog_step = 1.0/positions_mb.size();
	float prog_bar = 0.0;
	int n_interpolated = 0, i_locus = 0;
	double base, rate, dist;
	double mean_rate = (pos_cM.back() - pos_cM[0]) / (pos_bp.back() - pos_bp[0]);

	//Set up first positions to be mean rate
	while (i_locus<positions_mb.size() && positions_mb[i_locus] < pos_bp[0]) {
		base = pos_cM[0];
		dist = (pos_bp[0] - positions_mb[i_locus]);
		positions_cm[i_locus] = std::max(0.0,base - mean_rate * dist);
		n_interpolated ++;
		i_locus ++;
	}

	//Set up middle positions using interpolation
	int closest_pos = 1;
	for (; i_locus < positions_mb.size() ; ) {
		if (positions_cm[i_locus] == -1) {

			//Find suitable interpolation interval
			while (positions_mb[i_locus] > pos_bp[closest_pos] && closest_pos < pos_bp.size()) closest_pos++;

			//Interpolate
			if (closest_pos < pos_bp.size()) {
				assert(positions_mb[i_locus] < pos_bp[closest_pos]);
				assert(positions_mb[i_locus] > pos_bp[closest_pos-1]);
				base = pos_cM[closest_pos-1];
				rate = (pos_cM[closest_pos] - pos_cM[closest_pos-1]) / (pos_bp[closest_pos] - pos_bp[closest_pos-1]);
				dist = (positions_mb[i_locus] - pos_bp[closest_pos-1]);
				positions_cm[i_locus] = base + rate * dist;
				n_interpolated ++;
				i_locus ++;
			} else break;
		} else i_locus ++;
	}

	//Set up last positions to be mean rate
	while (i_locus < positions_mb.size()) {
		base = pos_cM.back();
		dist = (positions_mb[i_locus] - pos_bp.back());
		positions_cm[i_locus] = base + mean_rate * dist;
		n_interpolated ++;
		i_locus ++;
	}
	return n_interpolated;
}

void chunker::setGeneticMap(const gmap_reader & readerGM, const std::multimap < int, int >& map_positions, const std::vector<int>& positions_mb, std::vector<float>& positions_cm) {
	tac.clock();
	if (positions_mb.size() == 0) vrb.error("No variant found in region.");
	positions_cm = std::vector<float>(positions_mb.size(), -1);
	int n_set = setCentiMorgan(readerGM.pos_bp, readerGM.pos_cm, map_positions, positions_cm);
	int n_interpolated = interpolateCentiMorgan(readerGM.pos_bp, readerGM.pos_cm, map_positions, positions_mb, positions_cm);
	vrb.bullet("Map cM interpolation : [s=" + stb.str(n_set) + " / i=" + stb.str(n_interpolated) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void chunker::setGeneticMap(const std::multimap < int, int >& map_positions, const std::vector<int>& positions_mb, std::vector<float>& positions_cm) {
	tac.clock();
	if (positions_mb.size() == 0) vrb.error("No variant found in region.");
	positions_cm = std::vector<float>(positions_mb.size(), -1);
	positions_cm[0] = positions_mb[0] * 1.0f / 1e6;
	const float baseline = std::max(positions_cm[0] - 1e-9,1e-9);
	positions_cm[0] = 0;
	for (int l = 1 ; l < positions_mb.size() ; l ++) positions_cm[l] = (positions_mb[l] * 1.0f / 1e6) - baseline;
	vrb.bullet("Genetic map          : [Not provided]");
	vrb.bullet("Map cM constant      : [1cM=1Mb] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
