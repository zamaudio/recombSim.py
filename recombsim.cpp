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

#define MAX_GENOME 120 //0000

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

	//recomb_t recomb1 = { 20000, 1000 };
	//recomb_t recomb2 = { 60000, 8000 };
	//recomb_t recomb3 = { 100000, 16000 };
	recomb_t recomb1 = { 2, 10 };
	recomb_t recomb2 = { 60, 8 };
	recomb_t recomb3 = { 100, 16 };

	Sequence *seq0 = new BasicSequence("Ancestral Seq", "A", &AlphabetTools::DNA_ALPHABET);
	Sequence *seq1 = new BasicSequence("Recombination block 1", "A", &AlphabetTools::DNA_ALPHABET);
	Sequence *seq2 = new BasicSequence("Recombination block 2", "A", &AlphabetTools::DNA_ALPHABET);
	Sequence *seq3 = new BasicSequence("Recombination block 3", "A", &AlphabetTools::DNA_ALPHABET);
	
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
	seq->setName("Recombinant without mutations");

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
	for (j = 0; j < 100; j++) {
		Sequence *seq = seq0->clone();
		sprintf(label, "Recombinant with %d%% mutations",j);
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

	delete seq0;
	delete seq1;
	delete seq2;
	delete seq3;

	return 0;
}
