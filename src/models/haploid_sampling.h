#ifndef _HAPLOID_SAMPLING_H
#define _HAPLOID_SAMPLING_H

#include <utils/otools.h>
#include <containers/conditioning_set.h>

class haploid_sampling {
private:

	conditioning_set * C;

	//EXTERNAL DATA
	vector < vector < bool > > HAP;
	vector < unsigned char > UP;

	//DYNAMIC ARRAYS
	vector < vector < vector < float > > > backwardProbs;
	vector < vector < float > > forwardProbs;
	vector < vector < float > > prevForwardProbs;

public:
	//CONSTRUCTOR/DESTRUCTOR
	haploid_sampling(conditioning_set * C);
	~haploid_sampling();

	void backward();
	void sampling(vector < bool > &, vector < bool > &);
};

#endif
