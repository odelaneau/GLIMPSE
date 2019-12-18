#define _DECLARE_TOOLBOX_HERE
#include <chunker/chunker_header.h>

int main(int argc, char ** argv) {
	vector < string > args;
	for (int a = 1 ; a < argc ; a ++) args.push_back(string(argv[a]));
	chunker().chunk(args);
	return 0;
}

