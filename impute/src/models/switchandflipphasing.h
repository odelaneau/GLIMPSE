#ifndef _SAF_PHASING_H
#define _SAF_PHASING_H

#include <utils/otools.h>
#include <containers/conditioning_set.h>
#include <containers/haplotype_set.h>

class switchandflipphasing {
private:
	haplotype_set * H;
	conditioning_set * C;

	//EXTERNAL DATA
	std::vector < std::vector < bool > > HAP;
	std::vector < unsigned char > UP;

	//DYNAMIC ARRAYS
	std::vector < std::vector < std::vector < float > > > backwardProbs;
	std::vector < std::vector < std::vector < float > > > forwardProbs;
	//vector < vector < float > > prevForwardProbs;

public:
	//CONSTRUCTOR/DESTRUCTOR
	switchandflipphasing(haplotype_set * H, conditioning_set * C);
	~switchandflipphasing();

	void backward();
	void sampling(std::vector < bool > &, std::vector < bool > &);
};
#endif
