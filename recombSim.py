#!/usr/bin/env python
#
# Copyright (C) 2015  Damien Zammit <damien@zamaudio.com>
#
#  This program is libre software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# This program takes a single ancestral sequence in fasta format.
# It generates a tree of 14 new sequences with recombination at known positions
# so that black box testing of recombination detection algorithms can be performed.
#
# Usage: python recombTest.py
# -i fastafile    (Fasta input file)
# -n <integer>    (Number of recombinations  - 1 to 10 per step)
# -l <integer>    (Length of recombinations  - 1 to 1/10 of genome size)
# -d <float>      (Percentage divergence     - default 1% = 0.01)
# -s <integer>    (Pseudo random number seed - default 42)
# -o prefix       (Output prefix prepended   - default "out")

import string, re, collections
import os, sys, subprocess, datetime
from optparse import OptionParser
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import AlignIO
from Bio.Seq import MutableSeq

def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

### OPTIONS
	parser.add_option("-i", "--input", action="store", dest="inputfasta", help="Input alignment (Fasta)", default="")
	parser.add_option("-o", "--output", action="store", dest="outputfasta", help="Output prefix", default="out")
	parser.add_option("-l", "--lenrecomb", action="store", dest="lenrecomb", help="Length of recombinations (1000)", default="1000")	
	parser.add_option("-n", "--nrecomb", action="store", dest="nrecomb", help="Number of recombinations max 10 (1)", default="1")	
	parser.add_option("-d", "--divergence", action="store", dest="divergence", help="Divergence (1%=0.01)", default="0.01")	
	parser.add_option("-s", "--seed", action="store", dest="seed", help="Pseudo random number generator seed (42)", default="42")
	
	return parser.parse_args()


if __name__ == "__main__":

	(options, args) = main()
	
	### FUNCTIONS
	
	def ThrowError(ErrorString):
		print "Error: " + ErrorString + "\n"
		sys.exit()

        def getFastaLength(input_fasta):
                alignment = AlignIO.read(input_fasta, "fasta")
                totalSites = alignment.get_alignment_length()
                return totalSites

        def mergeFastas(inputs, output_fasta):
                os.system("".join(["echo -n \"\" > ", output_fasta]))
                for fasta in inputs:
                    os.system("".join(["cat ", fasta, " >> ", output_fasta]))

        def replaceChunkFasta(input_fasta, start, output_fasta):
                chunklength = int(options.lenrecomb)
                stop = start + chunklength
                chunk = SeqIO.read(input_fasta, "fasta")
                output = SeqIO.read(output_fasta, "fasta")
                change = output.seq.tomutable()
                chunkseq = chunk.seq
                change[start:stop] = chunkseq[0:chunklength]
                output_handle = open(output_fasta, "w")
                SeqIO.write(SeqRecord(change, name="", id="seq1", description=""), output_handle, "fasta")
                output_handle.close()
                cleanFasta(output_fasta)

        def cloneFasta(input_fasta, output_fasta):
                os.system("".join(["cp ", input_fasta, " ", output_fasta]))

        def cleanFastas(files):
                for f in files:
                    cleanFasta(f)

        def cleanFasta(filename):
                os.system("".join(["tail -n+2 ", filename, " > tmp"]))
                os.system("".join(["echo '>", filename, "' > ", filename]))
                os.system("".join(["tr -d '\\n' < tmp >> ", filename]))
                os.system("".join(["echo '' >> ", filename]))
                os.system("".join(["rm -f tmp "]))

        def subsetFasta(input_fasta, start, output_fasta):
                alignment = SeqIO.read(input_fasta, "fasta")
                stop = start + int(options.lenrecomb)
                subset = alignment[start:stop]
                output_handle = open(output_fasta, "w")
                SeqIO.write(subset, output_handle, "fasta")
                output_handle.close()
                cleanFasta(output_fasta)

        def mutateRecomb(input_fasta, globlength):
                #mutrate = float(options.lenrecomb) * float(options.divergence) / (100 * globlength)
                #mutrate = 1000 * float(options.divergence) / (2000 * float(options.lenrecomb))
                mutrate = float(options.divergence) / 2000
                runNetRecodon(mutrate, input_fasta, input_fasta)

	def runNetRecodon(mut, input_fasta, output_fasta):
                genomelength = float(getFastaLength(input_fasta))
                os.system("".join(["cat ", input_fasta, "|tail -n+2 > seqGMRCA2"]))
                os.system("cat seqGMRCA2 | tr -d '\\n' > seqGMRCA")
                os.system("".join(["NetRecodon6.0.0 -n1 -s2 -l", str(genomelength), " -e1000 -_1 -r0 -u", str(mut), " -f4 0.3731 0.1346 0.2110 0.2813 -v1.1588 3.6558 0.6834 0.4268 7.1333 1.0000 -a1.3381 -xseqGMRCA -*2 -y2 -#", options.seed])) #, " > /dev/null 2>&1"]))
                os.system("".join(["cat Results/sequences |tail -n+3 > ", output_fasta]))
                os.system("rm -fr Results")
                os.system("rm seqGMRCA")
                os.system("rm seqGMRCA2")
                cleanFasta(output_fasta)

        def forwardEvolve(input_fasta, output_fasta):
                length = float(getFastaLength(input_fasta))
                mut = 50 / (2000 * length)
                runNetRecodon(mut, input_fasta, output_fasta)

	#### MAIN PROGRAM
	
        output_fasta = options.outputfasta

        seq11 = "".join([output_fasta, "_1_1"])
        seq12 = "".join([output_fasta , "_1_2"])
        seq211 = "".join([output_fasta , "_2_1_1"])
        seq212 = "".join([output_fasta , "_2_1_2"])
        seq221 = "".join([output_fasta , "_2_2_1"])
        seq222 = "".join([output_fasta , "_2_2_2"])
        seq3111 = "".join([output_fasta , "_3_1_1_1"])
        seq3112 = "".join([output_fasta , "_3_1_1_2"])
        seq3121 = "".join([output_fasta , "_3_1_2_1"])
        seq3122 = "".join([output_fasta , "_3_1_2_2"])
        seq3211 = "".join([output_fasta , "_3_2_1_1"])
        seq3212 = "".join([output_fasta , "_3_2_1_2"])
        seq3221 = "".join([output_fasta , "_3_2_2_1"])
        seq3222 = "".join([output_fasta , "_3_2_2_2"])

        chunk = "".join([output_fasta , "_chunk"])

	if options.inputfasta == "":
		ThrowError("Provide input fasta file with -i")

        if int(options.nrecomb) > 10:
                ThrowError("Number of recombinations must be less than 10")

        globlength = getFastaLength(options.inputfasta)
        
        if int(options.lenrecomb) > globlength / 10:
                ThrowError("Length of recombinations must not exceed 1/10 of whole sample")
        
        a = globlength / 20
        pos = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]
        pos = a*pos

        ancseq0 = options.inputfasta
        forwardEvolve(ancseq0, seq11)
        cloneFasta(seq11, seq12)
        forwardEvolve(seq11, seq211)
        forwardEvolve(seq12, seq221)
        for i in range(0, int(options.nrecomb)):
            subsetFasta(seq221, pos[i], chunk)
            mutateRecomb(chunk, globlength)
            replaceChunkFasta(chunk, pos[i], seq211)
        
        cloneFasta(seq211, seq212)
        cloneFasta(seq221, seq222)
        forwardEvolve(seq211, seq3111)
        forwardEvolve(seq212, seq3121)
        forwardEvolve(seq221, seq3211)
        forwardEvolve(seq222, seq3221)
        for i in range(0, int(options.nrecomb)):
            subsetFasta(seq3121, pos[i], chunk)
            mutateRecomb(chunk, globlength)
            replaceChunkFasta(chunk, pos[i], seq3211)

        cloneFasta(seq3111, seq3112)
        cloneFasta(seq3121, seq3122)
        cloneFasta(seq3211, seq3212)
        cloneFasta(seq3221, seq3222)

        cleanFastas([
            ancseq0,
            seq11,
            seq12,
            seq211,
            seq212,
            seq221,
            seq222,
            seq3111,
            seq3112,
            seq3121,
            seq3122,
            seq3211,
            seq3212,
            seq3221,
            seq3222
            ])
        
        mergeFastas([
            ancseq0, 
            seq11, 
            seq12,
            seq211,
            seq212,
            seq221,
            seq222,
            seq3111,
            seq3112,
            seq3121,
            seq3122,
            seq3211,
            seq3212,
            seq3221,
            seq3222
            ], output_fasta)

        os.system("".join(["rm -f ", output_fasta, "_*"]))
        os.system("".join(["mv ",output_fasta, " ", output_fasta, ".fasta"]))

