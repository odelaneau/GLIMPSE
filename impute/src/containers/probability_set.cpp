#include <containers/probability_set.h>

	probability_set::probability_set(unsigned long n_site)
	{
		reference_haps = std::vector<std::vector<std::vector<int>>>(2,std::vector<std::vector<int>>(n_site, std::vector<int>()));
		prob_stateL = std::vector<std::vector<std::vector<float>>>(2,std::vector<std::vector<float>>(n_site, std::vector<float>()));
		prob_stateR = std::vector<std::vector<std::vector<float>>>(2,std::vector<std::vector<float>>(n_site, std::vector<float>()));
	}

	probability_set::~probability_set(){};


	void probability_set::addProbability
	(
		int hap,
		int marker,
		int refHap,
		float probL,
		float probR
	)
	{
		reference_haps[hap][marker].push_back(refHap);
		prob_stateL[hap][marker].push_back(probL);
		prob_stateR[hap][marker].push_back(probR);
	}

	void probability_set::normaliseProb
	(
		int hap,
		int marker,
		float sumProbL,
		float sumProbR
	)
	{
		uint32_t size = reference_haps[hap][marker].size();
		if (size == 1)
		{
			prob_stateL[hap][marker][0] = 1.0f;
			prob_stateR[hap][marker][0] = 1.0f;
			return;
		}

		float divProbL = 1.0f/sumProbL;
		float divProbR = 1.0f/sumProbR;

		#pragma GCC ivdep
		for (uint32_t i = 0; i<size; ++i)
		{
			prob_stateL[hap][marker][i] *= divProbL;
			prob_stateR[hap][marker][i] *= divProbR;
		}
	}

	void probability_set::normaliseProbNoSum
	(
		int hap,
		int marker
	)
	{
		const int size = reference_haps[hap][marker].size();
		if (size == 1)
		{
			prob_stateL[hap][marker][0] = 1.0f;
			prob_stateR[hap][marker][0] = 1.0f;
			return;
		}

		float sumProbL = 0.0f;
		float sumProbR = 0.0f;

		for (int i = 0; i<size; ++i)
		{
			sumProbL+=prob_stateL[hap][marker][i];
			sumProbR+=prob_stateR[hap][marker][i];
		}

		float divProbL = 1.0f/sumProbL;
		float divProbR = 1.0f/sumProbR;

		#pragma GCC ivdep
		for (int i = 0; i<size; ++i)
		{
			prob_stateL[hap][marker][i] *= divProbL;
			prob_stateR[hap][marker][i] *= divProbR;
		}
	}
