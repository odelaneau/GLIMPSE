#ifndef _SAF_PHASING_H
#define _SAF_PHASING_H

#include <utils/otools.h>
#include <containers/conditioning_set.h>

class switchandflipphasing {
private:

	conditioning_set * C;

	//EXTERNAL DATA
	vector < vector < bool > > HAP;
	vector < unsigned char > UP;

	//DYNAMIC ARRAYS
	vector < vector < vector < float > > > backwardProbs;
	vector < vector < vector < float > > > forwardProbs;
	//vector < vector < float > > prevForwardProbs;

public:
	//CONSTRUCTOR/DESTRUCTOR
	switchandflipphasing(conditioning_set * C);
	~switchandflipphasing();

	void backward();
	void sampling(vector < bool > &, vector < bool > &);
};
#endif
