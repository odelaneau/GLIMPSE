#ifndef _PROBABILITY_SET_H_
#define _PROBABILITY_SET_H_

#include <utils/otools.h>

class probability_set {
public:
	std::vector<std::vector<std::vector<int>>> reference_haps;
	std::vector<std::vector<std::vector<float>>> prob_stateL;
	std::vector<std::vector<std::vector<float>>> prob_stateR;

	probability_set(unsigned long n_site);
	~probability_set();

	void addProbability
	(
		int hap,
		int marker,
		int refHap,
		float probL,
		float probR
	);

	void normaliseProb
	(
		int hap,
		int marker,
		float sumProbL,
		float sumProbR
	);

	void normaliseProbNoSum
	(
		int hap,
		int marker
	);
};

#endif /* SRC_CONTAINERS_PROBABILITY_SET_H_ */
