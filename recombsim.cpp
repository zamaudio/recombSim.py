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
#include "recombsim.h"
#include "mtrand.h"

using namespace std;

#include <Bpp/Seq/Alphabet.all>
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Seq/Container.all>
#include <Bpp/Seq/Io/Fasta.h>

using namespace bpp;

#define FREQ_A 0.33
#define FREQ_C 0.17
#define FREQ_G 0.17
#define FREQ_T 0.33

#define MAX_GENOME 250000

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

int main(int argc, char** argv) {
	
	int i;
	int seed = 42;
	double dblrand = 0.0;
	MTRand53 twister(seed);
	MTRand irand(seed);

	recomb_t recomb1 = { 20000, 1000 };
	recomb_t recomb2 = { 60000, 8000 };
	recomb_t recomb3 = { 100000, 16000 };

	Sequence *seq0 = new BasicSequence("ANCESTR", "A", &AlphabetTools::DNA_ALPHABET);
	Sequence *seq1 = new BasicSequence("RECOMB1", "A", &AlphabetTools::DNA_ALPHABET);
	Sequence *seq2 = new BasicSequence("RECOMB2", "A", &AlphabetTools::DNA_ALPHABET);
	Sequence *seq3 = new BasicSequence("RECOMB3", "A", &AlphabetTools::DNA_ALPHABET);
	
	for (i = 0; i <= MAX_GENOME; i++) {
		dblrand = twister();
		seq0->addElement(i, (Recombsim::which_base(dblrand)));
	}
	for (i = 0; i <= recomb1.length; i++) {
		dblrand = twister();
		seq1->addElement(i, (Recombsim::which_base(dblrand)));
	}
	for (i = 0; i <= recomb2.length; i++) {
		dblrand = twister();
		seq2->addElement(i, (Recombsim::which_base(dblrand)));
	}
	for (i = 0; i <= recomb3.length; i++) {
		dblrand = twister();
		seq3->addElement(i, (Recombsim::which_base(dblrand)));
	}

	VectorSequenceContainer* reference = new VectorSequenceContainer(&AlphabetTools::DNA_ALPHABET);
	reference->addSequence(*seq0);

	int j;

	Sequence *seq = seq0->clone();
	seq->setName("MUT0");

	for (j = 0; j < recomb1.length; j++) {
		seq->setElement(recomb1.start + j, seq1->getChar(j));
	}
	for (j = 0; j < recomb2.length; j++) {
		seq->setElement(recomb2.start + j, seq2->getChar(j));
	}
	for (j = 0; j < recomb3.length; j++) {
		seq->setElement(recomb3.start + j, seq3->getChar(j));
	}
	
	VectorSequenceContainer* mfasta = new VectorSequenceContainer(&AlphabetTools::DNA_ALPHABET);
	mfasta->addSequence(*seq);
	delete seq;

	int pos;
	string snp;
	char label[64] = {0};
	for (j = 1; j <= 100; j++) {
		Sequence *seq = seq0->clone();
		sprintf(label, "MUT%d",j);
		seq->setName(label);
		for (i = 0; i < j*MAX_GENOME/100; i++) {
			pos = (int) (MAX_GENOME*twister());
			snp = Recombsim::which_base(twister());
			seq->setElement(pos, snp);
			//printf("pos=%d snp=%s\n",pos,snp.c_str());
		}
		mfasta->addSequence(*seq);
		delete seq;
	}

	Fasta fasWriter;
	fasWriter.writeSequences(string("simdata.mfasta"), *mfasta, true);
	fasWriter.writeSequences(string("simref.fasta"), *reference, true);

	delete seq0;
	delete seq1;
	delete seq2;
	delete seq3;
	printf("Dumped simulated sequences to output files\n");
	return 0;
}
