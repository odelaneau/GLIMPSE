#include "call_set_header.h"

call_set::call_set () {
}

void call_set::set(int _minDP, vector < double > _BIN, double _minPROB) {
	minDP = _minDP;
	minPROB = _minPROB;
	BIN = _BIN;
}

call_set::~call_set() {

}
