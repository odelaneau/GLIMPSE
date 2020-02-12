#define _DECLARE_TOOLBOX_HERE
#include <snparray/snparray_header.h>

int main(int argc, char ** argv) {
	vector < string > args;
	for (int a = 1 ; a < argc ; a ++) args.push_back(string(argv[a]));
	snparray().process(args);
	return 0;
}

