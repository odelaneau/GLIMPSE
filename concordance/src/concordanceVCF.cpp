#define _DECLARE_TOOLBOX_HERE
#include <checker/checker_header.h>

int main(int argc, char ** argv) {
	vector < string > args;
	for (int a = 1 ; a < argc ; a ++) args.push_back(string(argv[a]));
	checker().check(args);
	return 0;
}

