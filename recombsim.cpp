#include <cstdio>
#include "recombsim.h"
#include "mtrand.h"

#define FREQ_A 0.33
#define FREQ_C 0.17
#define FREQ_G 0.17
#define FREQ_T 0.33

#define MAX_GENOME 500000


Recombsim::Recombsim() {}

char Recombsim::which_base(double rnd) {
	if (rnd <= FREQ_A) {
		return 'A';
	} else if (rnd <= FREQ_A + FREQ_C) {
		return 'C';
	} else if (rnd <= FREQ_A + FREQ_C + FREQ_G) {
		return 'G';
	} else {
		return 'T';
	}
}

int main(int argc, char** argv) {
	
	int i;
	int seed = 42;
	double dblrand = 0.0;
	MTRand53 twister(seed);

	for (i = 0; i < MAX_GENOME; i++) {
		dblrand = twister();
		printf("%c", Recombsim::which_base(dblrand));
	}
	printf("\n");
	return 0;
}
