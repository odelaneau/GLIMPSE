#include "call_set_header.h"

call_set::call_set () {
}

void call_set::initialize(vector < double > _bins, double _T) {
	T = _T;
	bins = _bins;
	L = bins.size() - 1;
}

call_set::~call_set() {

}
