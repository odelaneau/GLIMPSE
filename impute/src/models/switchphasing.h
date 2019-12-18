#ifndef _S_PHASING_H
#define _S_PHASING_H

#include <utils/otools.h>
#include <containers/conditioning_set.h>
#include <containers/haplotype_set.h>

class switchphasing {
private:
	haplotype_set * H;
	conditioning_set * C;


	//EXTERNAL DATA
	std::vector < std::vector < bool > > HAP;
	std::vector < unsigned char > UP;

	//DYNAMIC ARRAYS
	std::vector < std::vector < std::vector < float > > > backwardProbs;
	std::vector < std::vector < float > > forwardProbs;
	std::vector < std::vector < float > > prevForwardProbs;

public:
	//CONSTRUCTOR/DESTRUCTOR
	switchphasing(haplotype_set * H, conditioning_set * C);
	~switchphasing();

	void backward();
	void sampling(std::vector < bool > &, std::vector < bool > &);
};

#endif
