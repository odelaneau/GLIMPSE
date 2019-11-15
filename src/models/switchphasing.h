#ifndef _S_PHASING_H
#define _S_PHASING_H

#include <utils/otools.h>
#include <containers/conditioning_set.h>

class switchphasing {
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
	switchphasing(conditioning_set * C);
	~switchphasing();

	void sampling(vector < bool > &, vector < bool > &);
};

#endif
