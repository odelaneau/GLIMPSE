#include "call_set_header.h"

call_set::call_set () {
}

void call_set::initialize(vector < double > _bins, double _T) {
	vrb.title("Initializing concordance engine");
	T = _T;
	bins = _bins;
	L = bins.size() - 1;
	vrb.bullet("Threshold = " + stb.str(T));
	vrb.bullet("#bins = " + stb.str(L));
	vrb.bullet("bins [" + stb.str(bins) + "]");
}

call_set::~call_set() {

}
