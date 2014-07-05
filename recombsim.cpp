/*
 * recombsim.cpp  Recombination in bacteria simulation data generator
 * Copyright (C) 2014  Damien Zammit <damien@zamaudio.com>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * For a full copy of the GNU General Public License
 * see <http://www.gnu.org/licenses/>.
 */

#include <cstdio>
#include <time.h>
#include <unistd.h>
#include "recombsim.h"
#include "mtrand.h"

using namespace std;

#include <Bpp/Seq/Alphabet.all>
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Seq/Container.all>
#include <Bpp/Seq/Io/Fasta.h>

using namespace bpp;

#define FREQ_A 0.30266
#define FREQ_C 0.19811
#define FREQ_G 0.19718
#define FREQ_T 0.30205

typedef struct recomb {
	int start;
	int length;
} recomb_t;

string Recombsim::which_base(double rnd) {
	if (rnd <= FREQ_A) {
		return string("A");
	} else if (rnd <= FREQ_A + FREQ_C) {
		return string("C");
	} else if (rnd <= FREQ_A + FREQ_C + FREQ_G) {
		return string("G");
	} else {
		return string("T");
	}
}

int main(int argc, char *argv[])
{
	int opt;
	int tfnd;
	int genomelen = 0;
	float maxmut = 0;
	int maxsteps = 0;

	tfnd = 0;
	while ((opt = getopt(argc, argv, "g:m:s:")) != -1) {
		switch (opt) {
		case 'g':
			genomelen = atoi(optarg);
			tfnd = 1;
			break;
		case 'm':
			maxmut = atof(optarg);
			tfnd = 1;
			break;
		case 's':
			maxsteps = atoi(optarg);
			tfnd = 1;
			break;
		default: /* '?' */
			fprintf(stderr, "Usage: %s [-g genomelen] [-m mutationscale] [-s maxsteps]\n", argv[0]);
			exit(-1);
		}
	}

	if (optind > argc) {
		fprintf(stderr, "Usage: %s [-g genomelen] [-m mutationscale] [-s maxsteps]\n", argv[0]);
		exit(-1);
	}

	int i;
	int seed = time(NULL);
	double dblrand = 0.0;
	MTRand53 twister(seed);
	MTRand irand(seed);

	recomb_t recomb1 = { 4000, 2000 };
	//recomb_t recomb1 = { 4, 5 };

	Sequence *seq0 = new BasicSequence("ANCESTR", "A", &AlphabetTools::DNA_ALPHABET);
	Sequence *seq1 = new BasicSequence("RECOMB1", "A", &AlphabetTools::DNA_ALPHABET);
	
	for (i = 0; i <= genomelen; i++) {
		dblrand = twister();
		seq0->addElement(i, (Recombsim::which_base(dblrand)));
	}
	for (i = 0; i <= recomb1.length; i++) {
		dblrand = twister();
		seq1->addElement(i, (Recombsim::which_base(dblrand)));
	}

	VectorSequenceContainer* reference = new VectorSequenceContainer(&AlphabetTools::DNA_ALPHABET);
	reference->addSequence(*seq0);

	int j;

	Sequence *seq2 = seq0->clone();
	seq2->setName("RECOMB");

	for (j = 0; j < recomb1.length; j++) {
		seq2->setElement(recomb1.start + j, seq1->getChar(j));
	}
	
	VectorSequenceContainer* mfasta = new VectorSequenceContainer(&AlphabetTools::DNA_ALPHABET);
	mfasta->addSequence(*seq0);
	mfasta->addSequence(*seq2);

	int pos;
	string snp;
	char label[64] = {0};
	Sequence *seqtmp = seq2->clone();
	for (j = 1; j <= maxsteps; j++) {
		sprintf(label, "MUT_%d", j);
		seqtmp->setName(label);
		for (i = 0; i < maxmut; i=i+maxmut) {
			pos = (int) (genomelen*twister());
			snp = Recombsim::which_base(twister());
			seqtmp->setElement(pos, snp);
			//printf("pos=%d snp=%s\n",pos,snp.c_str());
		}
		mfasta->addSequence(*seqtmp);
	}
	delete seqtmp;

	Fasta fasWriter;
	char file1[64];
	char file2[64];
	sprintf(file1, "sim.mfasta");
	sprintf(file2, "sim_ref.fasta");

	fasWriter.writeSequences(string(file1), *mfasta, true);
	fasWriter.writeSequences(string(file2), *reference, true);

	delete seq0;
	delete seq1;
	delete seq2;
	printf("Dumped simulated sequences to output files\n");
	return 0;
}
