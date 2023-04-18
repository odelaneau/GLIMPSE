#include <chunker/chunker_header.h>

std::vector < long int > chunker::getByPos (const long int pos, const std::multimap < long int, long int >& map_positions) {
	std::vector < long int > vecS;
	std::pair < std::multimap < long int , long int >::const_iterator , std::multimap < long int , long int >::const_iterator > ret = map_positions.equal_range(pos);
	for (std::multimap < long int , long int >::const_iterator it = ret.first ; it != ret.second ; ++it) vecS.push_back(it->second);
	return vecS;
}

long int chunker::setCentiMorgan(const std::vector < long int > & pos_bp, const std::vector < double > & pos_cM, const std::multimap < long int, long int >& map_positions, std::vector<float>& positions_cm) {
	long int cpt = 0;
	for (long int l = 0 ; l < pos_cM.size() ; l ++) {
		std::vector  < long int > vecS = getByPos(pos_bp[l], map_positions);
		for (long int si = 0 ; si < vecS.size() ; si ++) {
			positions_cm[vecS[si]] = pos_cM[l];
			cpt++;
		}
	}
	return cpt;
}

long int chunker::interpolateCentiMorgan(const std::vector < long int > & pos_bp, const std::vector < double > & pos_cM, const std::multimap < long int, long int >& map_positions, const std::vector<long int>& positions_mb, std::vector<float>& positions_cm) {
	float prog_step = 1.0/positions_mb.size();
	float prog_bar = 0.0;
	long int n_interpolated = 0, i_locus = 0;
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

	//Set up middle positions using long interpolation
	long int closest_pos = 1;
	for (; i_locus < positions_mb.size() ; ) {
		if (positions_cm[i_locus] == -1) {

			//Find suitable long interpolation long interval
			while (positions_mb[i_locus] > pos_bp[closest_pos] && closest_pos < pos_bp.size()) closest_pos++;

			//long interpolate
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

void chunker::setGeneticMap(const gmap_reader & readerGM, const std::multimap < long int, long int >& map_positions, const std::vector<long int>& positions_mb, std::vector<float>& positions_cm) {
	tac.clock();
	if (positions_mb.size() == 0) vrb.error("No variant found in region.");
	positions_cm = std::vector<float>(positions_mb.size(), -1);
	long int n_set = setCentiMorgan(readerGM.pos_bp, readerGM.pos_cm, map_positions, positions_cm);
	long int n_interpolated = interpolateCentiMorgan(readerGM.pos_bp, readerGM.pos_cm, map_positions, positions_mb, positions_cm);
	vrb.bullet("Map cM long interpolation : [s=" + stb.str(n_set) + " / i=" + stb.str(n_interpolated) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void chunker::setGeneticMap(const std::multimap < long int, long int >& map_positions, const std::vector<long int>& positions_mb, std::vector<float>& positions_cm) {
	tac.clock();
	if (positions_mb.size() == 0) vrb.error("No variant found in region.");
	positions_cm = std::vector<float>(positions_mb.size(), -1);
	positions_cm[0] = positions_mb[0] * 1.0f / 1e6;
	const float baseline = std::max(positions_cm[0] - 1e-9,1e-9);
	positions_cm[0] = 0;
	for (long int l = 1 ; l < positions_mb.size() ; l ++) positions_cm[l] = (positions_mb[l] * 1.0f / 1e6) - baseline;
	vrb.bullet("Genetic map          : [Not provided]");
	vrb.bullet("Map cM constant      : [1cM=1Mb] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
