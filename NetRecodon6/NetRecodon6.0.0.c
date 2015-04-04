/* NetRecodon.c
// Programmer:	Miguel Arenas and David Posada. NetRecodon is a big extending from the program "Recodon" (M. Arenas and D. Posada) and before, from "Recoal" by D. Posada.
// Purpose:	To simulate codon sequences with recombination (intra and inter codon breakpoints), migration, and complex population models.
//			NETRECODON ALLOWS:
//								- Recombination with breakpoints in intra and inter codon positions.
//							 	- Codon models: M0, M1, M7, M8, discrete gamma, continuous gamma, user-probabilities for user-omegas. Output file with the omega simulate for each codon.
//								- Print ancestral sequences: catMRCA (concatenate sequence from MRCAs of all recombinant fragments) and GMRCA.
//								- Output alignments by phylip, fasta and nexus formats.
//								- Tip Dates. The dates of the tip lineages can be different.
//								- Next versions for users: 1) output file with nodes connections of the networks. 2) Print summary settings to a file.
//
//// HISTORY: 
//
///
// Version 2.0 (February 2008)
//		- Improve kappa Qij matrix for codon models.
//		- NetRecodon2.0.0. Recombination in codons.
//		- Option to print in a file the settings for the programmer (doSettingsFile).
//		- Option to print just the MRCA in separate files (doOutMRCAfiles).
//		- Option to print just the Branches of the Net in separate files (doBranchNetfiles). Bug solved in codMatrix.
//		- Version 2.0.8: The Branches of the Net files are in a new format where a recombinant generates just a recombinant node.
//		- Version 2.0.9: The type of substitutions in codon models is improved.
//
// Version 2.2 (October 2008)	
//		- For codon models (not for nuc models where the printed MRCA printed is the concatenate MRCAs) in presence of recombination, NetRecodon print the GMRCA but no the concatenate MRCAs (sumMRCAs). 
//		- NetRecodon prints both GMRCAs and the concatenate MRCAs for codon models. For nucleotide models only concatenate MRCAs.
//
// Version 2.3 (October 2008)	
//		- More codon models introduced: M1 (omega = 0 + omega = 1), M7 and three for M8: M8 (D.beta + omega > 1), M8a (D.beta + omega = 1), M8b (D.beta + omega < 1). Beta distribution have been developed by David Posada, thanks David!. 
//		- Improved a small bug of Version 2.2. (printing GMRCA file, free memory). 
//		- Output file with the simulated omega for each site, this is for M1, M7 and M8 codon models.
//
//			2.3.2 (October 2008)
//		- All the codon models with omega variable have the omega variable per codon, not per branch.
//		- Output file with the simulated omega for each site for every codon model with variable omega.
//
//			2.3.3 (October 2008)
//		- Option for haploid/diploid from user.
//
//			2.3.4 (November 2008)
//		- Minor fixes.
//		- Print sequences for phylip output files with migration in style: "s00001_p1.. and outgrp_p0".
//	    - Print output sequences in fasta and nexus formats.
//
//			2.3.5 (November 2008)
//		- Improve M1 codon model, now P0 with omega <1 (not only neccesary omega = 0) and P1 with omega = 1.
//
// Version 2.4 (January 2009)
//		- Tip Dates. The dates of the tip lineages can be different. Special thanks to David Posada for the code of this improvement. 
//			In this version the tip dates only are working when there is not migration.
//			e.g. 4 samples; 1995:seq1, 2003:seqs 4 and 6, and so on. Where the generation time is also introduced from the user.
//
//			2.4.1 (January 2009)
//		- Tip dates also in presence of migration. 
//			In presence of convergence of demes, it cannot have tip dates with a time higher (older) that the time of the event of convergence of demes. Because at that time the deme of the tip date does not exist (it was converged)..
//		- Minor fixes about print settings of codon models.
//
//			2.4.2 (March 2009)
//		- Minor fixes (Improved memory).
//
//
//			3.0.0 (June 2009)
//		- Nucleotide evolution using the ARG, the network, such as the codon evolution. It prints now both catMRCA and GMRCA as output files.
//
//			3.0.1 (July 2009)
//		- Breakpoints in trapped material are cosiderated as well, for both Gi and selected breakpoints.
//
//			5.0.5 (October 2009)
//		- Intracodon breakpoints taking into account pseudo-ancestral material (non ancestral material that can affect to athe ancestral material in intracodon recombinations). 
//			Print sequence for the last node (global GMRCA) and for the concatenated MRCA sequence for all material, for codon models (intracodon rec).
//
// 			5.0.6 (October 2009)
//		- Print sequence for the last node of the ancestral material (ancestral material GMRCA), for codon models (intracodon rec).
//
//			5.0.7 (September 2010 - February 2011)
//		- Print Network if "doBranchNetfiles = YES". The particular case of, a recombination followed by a coalescent event of the two previous generated recombinant lineages, does not result in loop (it gives a linear lineage)!!
//		- Fixed a bug in the counter of total non-synonymous and synonymous substitutions. Thanks to Suzanne English for her contribution!
//
//			5.0.8 (February 2011)
//		- Several Memory improvements. Including for genomes.
//
//			5.0.9 (March 2011)
//		- Activation of negative growth rates.
//
//			6.0.0 (March 2011)
//		- Improvements in demography by periods: Nbegin / Nend correction, timeCA when N is constant in the period, consideration of (cumDuration[period] - cumDuration[period-1]).
//
*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>

#ifndef macintosh
	#include <sys/types.h>
	#include <sys/stat.h>
#endif
#ifdef MAC
	#include <sioux.h>
	#include <console.h>
	#include <unix.h>
#endif

#ifdef MPI
	#include "mpi.h"
	#define nmaxproc 16
	#define comm MPI_COMM_WORLD
#endif

#define PROGRAM_NAME		"NetRecodon"
#define VERSION_NUMBER		"6.0.0"
#define	NO					0
#define	YES					1
#define NUMCOD				61
#define	pos(i,j,n)			((i)*(n)+(j))
#define post(i,j,n)			((j)*(n)+(i))
#define INCREMENT_NODES		500
#define INCREMENT_SEGMENTS	500
#define USER_INPUT
#undef	USER_INPUT
#define HUDSON_UNITS
#undef	HUDSON_UNITS  
/**/


typedef struct segment
	{
	struct	segment *before1, *before2, *after1, *after2;
	int		sIndex;
	int		sStart;
	int		sEnd;
	int		sIndexNode;
	}
	TreeSegment;
	
typedef struct node
	{
	struct node		*left, *right, *anc1, *anc2, *outgroup, *sib;
	int				index, label, isOutgroup; 
	int			   	numSegNode;
	double			length, time;
	int				indexOldMigPop;
	int				indexCurrentMigPop;
	int 			class, breakp, passNumber, breakCodon, whereBreakCodon;
	int				NetLabelPrint;
	int				MRCAfrom, MRCAto;
	int				GMRCA_ancestral;
	int				*SitesNonAncHere;
	}
	TreeNode;

typedef struct nodex
	{
	struct nodex	*left, *right, *anc1, *outgroup;
	int				index, label, isOutgroup, NetIndex;
	int				indexOldMigPop; 
	double			length, time;
	int				MRCAfrom, MRCAto;
	/*struct Pij*/ /* Future: Because each node can has several Pij matrixes, a Pij per category.. */
	}
	TreeNodex;

typedef struct QijOmegaCat
	{
	double Root_C_cat[NUMCOD];
	double Cijk_C_cat[NUMCOD*NUMCOD*NUMCOD*NUMCOD];
	}
	Qij_OmegaCat;

typedef struct sample
	{	
	float		time;
	int		size;
	int		*member;
	}
	SampleSt;



	/* Read */
static void 	PrintTitle (FILE *file);
static void 	PrintDate (FILE *file);
static void		PrintUsage();
#ifdef USER_INPUT
static void 	UserInput (long int *seed);
#endif
static void		ReadUntil (FILE *fv, char stopChar, char *what);
static void		PrintRunSettings (FILE *filep, long int seed);
static void		ReadParametersFromFile ();
static void		ReadParametersFromCommandLine (int argc, char **argv);
	/* Nucleotide Models */
static void 	HKY (double Pij[4][4], double time, double kappa, double rate, double p_i[4]);
static void 	GTR (double Pij[4][4], double branchLength, double varRate, double p_i[4]);
static void		GTnR (double Pij[4][4], double branchLength, double varRate, double p_i[4]);
static void 	SubstitutionMatrix (double ch_prob[4][4], double time, double kappa, double rate, double p_i[4]);
/*static void		SimulateDataForSite (TreeNodex *p, int siteNum, int numSites, double mutationRate, double kappa, double p_i[4], double rate, long int *seed, char *MRCAsequence);*/
/*static void		EvolveSequenceOnTree (long int *seed, double mutationRate, double kappa, double alpha, double p_i[4], int numSites, int *arrayIndBreakpointsOrd, char *MRCAsequence);*/
static void		EvolveSequenceOnTree_NEW (long int *seed, double mutationRate, double kappa, double alpha, double p_i[4], int numNuc, int indNumRE, int *arrayIndBreakpointsOrd, char *MRCAsequence, int numSites);
static char		WhichNuc (int nucId);
static int		WhichNucNumber (char siteLetter);
static int		TellMeGMRCALabel (TreeNode *p);
static void 	SimulateDataForSite_Nucleotide_RECURSIVE_NET (TreeNode *p, int siteNucleotide, int numSites, double mutationRate, double rate, double kappa, long int *seed);

/*static void 	SetMatrix(double Pij[4][4], double len);*/
/*static void 	SetVector(double *vector, short base, double len);*/
	/* Codon Models */
/*static void		SimulateDataForSite_Codon (TreeNode *p, int siteCodon, int numSites, double mutationRate, int numOmegaCat, double rate, long int *seed);*/
static void		EvolveSequenceOnTree_Codon (long int *seed, double mutationRate, double alpha, int numNuc, int indNumRE, int *arrayIndBreakpointsOrd, char *MRCAsequence, int numOmegaCat, int numSites);

static void 	SimulateDataForSite_Codon_RECURSIVE_NET (TreeNode *p, int siteCodon, int numSites, double mutationRate, int numOmegaCat, double rate, long int *seed, int Cbroke);
static int		CombineTwoCodons (int InCodon1, int InCodon2, int brokePosition);

static void		CodonModel (double Pij[NUMCOD][NUMCOD], double branchLength, double varRate);
static void		CodonModel_Cat (double Pij[NUMCOD][NUMCOD], double branchLength, double varRate);
static int		EnterCodonMRCA_Freq (TreeNode *p, int siteNum, int sitePosition, int numNuc, int out_C[4], int codon[3]);
static int		EnterCodonMRCA_File (TreeNode *p, int siteNum, int numNuc, char *MRCAsequence, int out_C[4]);
static int		makeCodonFromNuc (int x, int y, int z);
static double	codonTable_frequencies_MRCA (int cod);
static double	codonTable_frequencies (int cod);
static void		buildCodonMatrix_Qij_Cijk ();
static int		numdif_codon (int codon1, int codon2);
static void		number_to_codon (int ind, char out[]);
static void		number_to_codon_MRCA(int ind, int codon[]);
static void 	number_to_codon2(int ind, int out[]);
static int		codonTable_DnDs(int cod);
static int		codon_tr_tv (int indi, int indj);
static double	codon_Rmat(int indi, int indj);
static double	codon_NRmat(int indi, int indj);
static int		gammasCalculate (double alpha_d, int numCategories);
	/* Eigen */
static int		EigenREV (double Root[], double Cijk[]);
static int		EigenREV_Codon (double Root_C[], double Cijk_C[]);
/*static double *CalcREprobs (int numSequences, int rateType);*/
	/* Print Sequences */
static void 	PrintSequences (/*int replicate*/);
static void 	PrintAncestralSequences (/*int replicate*/);
static void 	PrintSequences_C (/*int replicate*/);
static void 	PrintAncestralSequences_C (/*int replicate*/);
static void 	PrintOutMRCAFiles_C (/*int replicate*/);
static void 	PrintOutMRCAFiles_C_Conc (/*int replicate*/);
static void 	PrintOutMRCAFiles (/*int replicate*/);
static void 	PrintOutMRCAFiles_Conc (/*int replicate*/);
static void		PrintOutGMRCAFiles_Codon_AncestralMat();
static void 	PrintSequences_FASTA (/*int replicate*/);
static void 	PrintAncestralSequences_FASTA (/*int replicate*/);
static void 	PrintSequences_C_FASTA (/*int replicate*/);
static void 	PrintAncestralSequences_C_FASTA (/*int replicate*/);
static void 	PrintSequences_NEXUS (/*int replicate*/);
static void 	PrintAncestralSequences_NEXUS (/*int replicate*/);
static void 	PrintSequences_C_NEXUS (/*int replicate*/);
static void 	PrintAncestralSequences_C_NEXUS (/*int replicate*/);
static void 	PrintNEXUS_initial ();
static void 	PrintNEXUS_end ();
	/* Random and others */
static double	RndGamma (double s, long int *seed);
static double	RndGamma1 (double s, long int *seed);
static double	RndGamma2 (double s, long int *seed);
static double	RandomGamma (double shape, long int *seed);
static double	RandomGamma1 (double s, long int *seed);
static double	RandomGamma2 (double s, long int *seed);
static double	RandomBeta (double p, double q, long int *seed);
/*static void	RandomizeArray(int array[], int first, int last);*/
/*static int	Rand (int max);*/
/*static double Rndu (void);*/
static double 	RandomUniform (long int *seed);
static double	RandomExponential (double mean, long int *seed);
	/* Total Segments Functions */
static void		MakeCoalescenceTree (int numSequences, int numSites, int numNuc, int N, double recombinationRate, int numPopulations, long int *seed);
static int		bbin (double dat, double *v);
static int		bbinDemes (double dat, double *v, int n);
static int		bbinInOmegaCat (double dat, double *v, int n);
static int		bbin_EnterMRCA (double dat, double *v);
static int		CalcIndividualGi (int Individual, TreeNode *nodes, int *activeGametes, int numNuc, int *S_MRCA, int sizeNode);
static int		IsValidBreakSite (int *activeGametes, TreeNode *nodes, int whichInd, int whichSite, int *S_MRCA);
static int		CountsForExpNumRec (int *activeGametes, int whichInd, int whichSite, TreeNode *nodes, int *S_MRCA, int sizeNode); 
static int		recSegmentsGeneratesLeft (int nodeValue, TreeSegment *s, TreeSegment *n, int numNuc, int whichSite, int *actSegIndex);
static int		recSegmentsGeneratesRight (int nodeValue, TreeSegment *s, TreeSegment *m, int numNuc, int whichSite, int *actSegIndex);
static int		recSegmentsGeneratesLeftBrokenCodon (int nodeValue, TreeSegment *s, TreeSegment *n, int numNuc, int whichSite, int LeftLess, int RightHigh, int *actSegIndex);
static int		recSegmentsGeneratesRightBrokenCodon (int nodeValue, TreeSegment *s, TreeSegment *m, int numNuc, int whichSite, int LeftLess, int RightHigh, int *actSegIndex);
static int		overLapSegmentsCoalMRCA (TreeNode *p, TreeNode *q, int sizeNode_p, int sizeNode_q, int site);
static void		buildTreeCoal (TreeNode *p, TreeNodex *f, int numSequences, int *numActNodex);
static void		buildTreeInit (TreeNode *p, TreeNodex *f, int numNuc, int *arrayIndBreakpointsOrd, int whoBreakp, int numSequences, int *numActNodex);
static void		buildTreeEnd (TreeNode *p, TreeNodex *f, int numNuc, int *arrayIndBreakpointsOrd, int whoBreakp, int numSequences, int *numActNodex);
static void		buildTreeIntern (TreeNode *p, TreeNodex *f, int numNuc, int *arrayIndBreakpointsOrd, int whoBreakp, int numSequences, int *numActNodex);
static int		IndexSeg (TreeNodex *f);
static int		LabelSeg (TreeNodex *f);
static void		ListTimesSeg (TreeNodex *f);
static void		PrintTimesSeg ();
static void		PrintTreesSeg (int replicate, int indNumRE, int numNuc, int *arrayIndBreakpointsOrd);
static void		WriteTreeSeg (TreeNodex *f);
static void		RelabelNodesSeg (TreeNodex *f);
static double	roundit(double d, int dig);
long double		roundl(long double x);

TreeSegment		*segments;
TreeNode		*nodes, **treeRootInit;
TreeNodex		*nodex, **treeRootNodex;
Qij_OmegaCat	*QijOmegas;
long int		userSeed;
int				*matrix, *matrixC, *MRCA, *matrixCnuc,/*tipLabel,*/ intLabel, *outgroup, *breakpoint, *arrayIndBreakpointsOrd;
int				*S_MRCA, actSegIndex, *OnlyAncS_MRCA;
int				numSequences, N, numSites, numNuc, numDataSets, numRE, indNumRE, recNotToCount, numCA, numMU, noisy, numIndividuals, numCodons, numMIG, numCONV, numREbreakCod, numStopCodonREC; 
int				doRateHet, equalBaseFreq, equalBaseFreqCod;
double 			cumNumRE, cumNumREntc, cumNumCA, cumNumMU, meanNumRE, meanNumREtc, meanNumCA, meanNumMU, meanNumMIG, expNumRE, rho, cumNumMIG, cumNumCONV, meanNumCONV, cumNumREbreakCod, cumNumStopCodonREC, meanNumREbreakCod, meanNumStopCodonREC;
double			kappa, titv, alpha, p_i[4], p_i_codon[12], pinv, mutationRate, migrationRate, OmegaRateHet;
double 			recombinationRate;
double 			Rmat[6], Qij[16], Cijk[256], Root[4], mr, tstv, NRmat[12];
double 			Qij_p[16], Cijk_p[256], Root_p[4], mr_p, tstv_p;
double 			Rmat_C[NUMCOD+2], Qij_CC[NUMCOD*NUMCOD], Cijk_C[NUMCOD*NUMCOD*NUMCOD*NUMCOD], Root_C[NUMCOD], Qij_C[NUMCOD][NUMCOD], OmegaInit;
double 			omega;
char			alignmentFile[20], treeFile[20], timesFile[20], model[4], breakpointFile[20], MRCAFile[20], screenFile[20], settingsFile[20], MRCAfilePrint[20], BranchNetFilePrint[20], ConcMRCAfilePrint[20], OmegasPerSiteFile[30], GMRCAancFilePrint[20];
static double	outgroupBranchLength, expTMRCA;
static int		thereisOutgroup, zeroRec, doExponential, doPrintTrees, doPrintTimes, doPrintBreakpoints, doHKY, doGTR, doCodonModel, doPrintAncestralSequences, 
				doSeparatedSequences, doCountsForExpNumRec, doDemographics, doMRCAFile, doCodon_HKY, doCodon_GTR, doMigration, doOmegaCat, doOmegaProb, doOmegaRateHetCont, 
				doOmegaRateHetDisc, doMPI, doGTnR, doCodon_NGTR, doBadReplicate, doConvergDemes, doConvNext, doSettingsFile, doOutMRCAfiles, doBranchNetfiles, doPrintOmegasPerSitefiles, doPrintFASTA, doPrintNEXUS,
				ThisBreakpIsTrapped, doGMRCAsamp;
static int		*Nbegin, *Nend, *cumDuration, numPeriods;
static double	*periodGrowth, growthRate;
int				fixedNumRecEvents, doFixNumRecEvents;
double			counterTime, counterTimeInit, actualTGMRCA, countTMRCA, countTMRCAReps, varianceGMRCArep, varianceTrep, varianceErep, *convDemTimes, *convDemTimes_old;
int				numNodex, numPopulations, numOmegaCat, *numNodesInitPopul, *initPopulation, numConvergDemes, *deme_a, *deme_b, nextConvNumber, *deme_a_old, *deme_b_old, *currentConvDem, totCurrentConv, distance, *NodesMRCAposit;
static double	*varTimeGMRCA, *varTimeT, *omegaVal, *omegaProb, *gammaRates, *omegaValGammaRate;
static int		*varEvent;
int				freqNumber, nextAvailable, doRepitEvol;
static int		nDIGITS;
int				numEqual2, numEqual1, numDifCodSameAA, numDifCodDifAA, numNonSyn0, numNonSyn1, numNonSyn2;
double 			cumNumEqual2, cumNumEqual1, cumNumDifCodSameAA, cumNumDifCodDifAA, meanNumEqual2, meanNumEqual1, meanNumDifCodSameAA, meanNumDifCodDifAA, meanNumNonSyn0, cumNumNonSyn0, meanNumNonSyn1, cumNumNonSyn1, meanNumNonSyn2, cumNumNonSyn2;
int 			numNetLabelPrint;
int 			numMU_S, numMU_NS; 
double			cumNumMU_S, cumNumMU_NS, meanNumMU_S, meanNumMU_NS, expVarTMRCA;
static int		doM1, doM8, doM7, doM0;
double			M1_P0_omeg0, M1_P1_omeg1, M1_omega0;
int				M1_FinalSite_omega, ProbCategory, GammCategory, formatNumber;
double 			M8_P0_beta, M8_P1_omega, M8_p_beta, M8_q_beta, M8_omegaP1, M8_omegaP0, M8_FinalSite_omega, omegaGammaDisc;
double			M7_p_beta, M7_q_beta, M7_FinalSite_omega;
int 			Nscaling;
SampleSt		*datedSample, tempSample;
static int		doDatedTips, numTipDates;
static double	generationTime, latestSamplingTime;
int				*SitesNonAncHere;
int				numNodes;


FILE			*fpAlignment, *fpTrees, *fpTimes, *fpSNP, *fpBreakpoints, *fpMRCA, *fpScreen, *fpmpi, *fpSettings, *fpMRCAprint, *fpBranchNet, *fpConcMRCAprint, *fpOmegasPerSitePrint, *fpGMRCAancPrint;
/* NumRE = Number of total recombination events */
/* NumCA = Number of coalescence events */
/* NumMU = Number of mutational events */
/* CumNumx = Acumulating */
/* numDataSets = Number of replicates */
/* zeroRec = reps with 0 recombination events */
/* fixedNumRecEvents, doFixNumRecEvents = it's to fix a number of Rec events */

// ************************** Variables for MPI *********************

#ifdef MPI

  MPI_Status  status;
  int p, rank, lcola, ierror, root, rep;
  MPI_Request request,requests[nmaxproc];
  int fila, filas[nmaxproc], nodo, indice, resto;
  int zero=0;
  int ii, flag, nada;
  
//    INTEGER tipo,tipoi,corte,count
//    INTEGER ndims,dims(nmaxdims),nlz,iii,kk,nil,ii,jj
//    INTEGER ncapt,counts(nmaxproc),disps(nmaxproc)

#endif     

// ************************************************************************


/******************** MAIN *******************/
int main(int argc, char **argv)
{
	int			dataSetNum, i, k, j, w, sum, pass, sorted, mmm;
	long int	seed, seedFirst, originalSeed;
	double		a, b;
	float 		start, secs;
	char		File[80];
	char		dir[80];
	char		*MRCAsequence;				/* this array will contain the MRCA sequence from inputFile */
	double		*fragmentTMRCAfraction;		/* this array will have the TMRCA of the agments */
	double		specMigPropPopul;
	FILE 		*fp;
	
	/* defaults */
	numMU_S = numMU_NS = mmm = 0;
	cumNumMU_S = cumNumMU_NS = meanNumMU_S = meanNumMU_NS = meanNumREtc = 0.0;
	varianceGMRCArep = varianceTrep = varianceErep = 0.0;
	MRCAsequence = NULL;
	numDataSets = 10; /* the number of samples to simulate */
	numSequences = 6; /* number of gametes in each data set */
	N = 1000; /* effective population size */
	numPeriods = 0;	/* number of distinct demographic periods */						
	recombinationRate = 0.0; /* recombination rate per site per generation */
	numSites = 201;
	numNuc = 201;
	numCodons = 67;
	mutationRate = 1e-7; /* mutation rate per site per generation */
	doGTR = NO;
	doHKY = NO;
	doCodonModel = NO;
	doCodon_HKY = NO;
	doCodon_GTR = NO;
	doCodon_NGTR = NO;
	doOmegaCat = NO;
	numOmegaCat = 0;
	doOmegaRateHetCont = NO;
	doOmegaRateHetDisc = NO;
	titv = -1;
	omega = OmegaInit = -1;
	doMigration = NO;
	migrationRate = -1;
	numPopulations = 1;
	doConvergDemes = NO;
	doConvNext = NO;
	numConvergDemes = 0;
	growthRate = 0;
	doRateHet = NO;
	alpha = 0.0;
	pinv = 0.0;
	equalBaseFreq = YES;
	p_i[0] = 0.25;
	p_i[1] = 0.25;
	p_i[2] = 0.25;
	p_i[3] = 0.25;
	equalBaseFreqCod = YES;
	for (i = 0; i < 12; i++)
		p_i_codon[i] = 0.25;
	for (i = 0; i < 6; i++)
		Rmat[i] = -1;
	for (i = 0; i < 12; i++)
		NRmat[i] = -1;
	noisy = 1;
	thereisOutgroup = NO;
	outgroupBranchLength = 0.1;
	doDemographics = NO;					
	doExponential = NO;						
	doMRCAFile = NO;
	doPrintTrees = NO;
	doPrintTimes = NO;
	doPrintBreakpoints = NO;
	strcpy(alignmentFile, "sequences");		/* alignmentFile = alignment (char type)*/
	strcpy(MRCAfilePrint, "GMRCA");
	strcpy(ConcMRCAfilePrint, "MRCA");
	strcpy(GMRCAancFilePrint, "ancGMRCA");
	strcpy(BranchNetFilePrint, "");
	strcpy(OmegasPerSiteFile, "SimulateOmegasPerSite");
	doPrintAncestralSequences = NO;
	doSeparatedSequences = NO;				/* it writes the sequences in diferents files, a file per replicate */
	doFixNumRecEvents = NO;					/* do a number of desire recombination events */
	fixedNumRecEvents = 0;
	seed = time(NULL);						/* seed for the random number */
	userSeed = 0;							/* seed entered by the user */
	freqNumber = 4;
	doCountsForExpNumRec = YES;		
	i = k = 0;
	b = specMigPropPopul = 0.0;
	Nscaling = 2; /* controls the scaling of time in Nscaling generations 1=haploids 2=diploids */
	numNodes = 0;
	doBadReplicate = NO;
	doMPI = NO;
	doRepitEvol = NO;
	doGMRCAsamp = NO;
	doPrintFASTA = NO;
	doPrintNEXUS = NO;
	doDatedTips = NO;			/* contemporaneous samples */

	/* MORE OPTIONS */
	doSettingsFile = NO;
	doOutMRCAfiles = NO; /* print MRCA and GMRCA */
	doBranchNetfiles = NO;
	doPrintOmegasPerSitefiles = NO; /* print omegas per site in a file */
	




// ***********************   MPI Inicializations    ****************************
#ifdef MPI
 
  MPI_Init(&argc,&argv);

        lcola=1;
        root=0;

// Number of processors
   MPI_Comm_size( MPI_COMM_WORLD, &p);
   if (p<2) {
   	printf("\n***ERROR: At least 2 processors are needed for parallel execution\n\n");
        MPI_Finalize();	 
	exit(1);
   }
   if (p>nmaxproc) {
   	printf("\n***ERROR: Number of processors higher than the maximum (%d) \n\n",nmaxproc);
        MPI_Finalize();	 
	exit(1);
   }
   

// Determine Process number
   MPI_Comm_rank( MPI_COMM_WORLD, &rank);

        printf("Procesor  %d of %d\n",rank,p);	//Depuracion, puede eliminarse
	doMPI = YES; // Deberia no ser necesaria, aunque por si acaso se puede dejar de momento			Ok.
#endif
// ************************************************************************


		/* arguments from external place */
	/* in macintosh the arguments are read from a file. For the other OS
	arguments are given at the command line */
	#ifdef MAC
	/*	_fcreator = 'R*ch';
		_ftype = 'TEXT';*/
		if ((fp = freopen("parameters", "r", stdin)) == NULL)				/* input from parameters file and it reads */
			{
			/* Anadido Mifuel: Para imprimir por pantalla del master.. */
			#ifdef MPI 
					fprintf (stderr, "\n%d: ERROR: Can't read parameters file.",rank);
			#else
				fprintf (stderr, "\nERROR: Can't read parameters file.");
			#endif
			PrintUsage();
			}
		ReadParametersFromFile();
		Sioux();
		fclose(fp);
	#else
	ReadParametersFromCommandLine (argc, argv);								/* input from external arguments */
		if (argc < 2)
			{
			if ((fp = freopen("parameters", "r", stdin)) != NULL) 	
				ReadParametersFromFile();
			else
				{
				/* Anadido Mifuel: Para imprimir por pantalla del master.. */
				#ifdef MPI
					fprintf (stderr, "%d, \nERROR: No parameters specified (use command line or parameter file)",rank);
				#else
					fprintf (stderr, "\nERROR: No parameters specified (use command line or parameter file)");
				#endif
				PrintUsage();
				}
			}
	#endif 
	/* to compile default input files in NO */			
	#ifdef USER_INPUT /* ask the user to input parameter values */
		UserInput(&seed);
	#endif


	/***** Define model *****/
	OmegaInit = omega;
	if (OmegaInit >= 0.0 || numOmegaCat > 0 || doM1 == YES || doM8 == YES || doM7 == YES) /* codon model */
		{
		doCodonModel = YES;
		if (titv >= 0.0)
			doCodon_HKY = YES;
		for (i = 0; i < 6; i++)
			{
			if (Rmat[i] >= 0)
				{
				doCodon_GTR = YES;
				if (doCodon_HKY == YES)
					{
					#ifdef MPI
					if (rank==root)
						fprintf(stderr, "\n Introduce only a titv or only a R-matrix or only a N-Rmatrix ");
					#else
						fprintf(stderr, "\n Introduce only a titv or only a R-matrix or only a N-Rmatrix ");
					#endif
					PrintUsage();
					}
				break;
				}
			}
		for (i = 0; i < 12; i++)
			{
			if (NRmat[i] >= 0)
				{
				doCodon_NGTR = YES;
				if (doCodon_HKY == YES || doCodon_GTR == YES)
					{
					#ifdef MPI
					if (rank==root)
						fprintf(stderr, "\n Introduce only a titv or only a R-matrix or only a N-Rmatrix ");
					#else
						fprintf(stderr, "\n Introduce only a titv or only a R-matrix or only a N-Rmatrix ");
					#endif
					PrintUsage();
					}
				break;
				}
			}
		if (doCodon_HKY == NO && doCodon_GTR == NO && doCodon_NGTR == NO) /* default model */
			{
			doCodon_HKY = YES;
			titv = 0.5;
			}
		}
	else				/* nucleotide model */
		{
		if (titv >= 0.0)
			doHKY = YES;
		if (freqNumber != 4)
			{
			fprintf (stderr, "PARAMETER ERROR: Bad number of enter frequencies if you do not use a codon model (%d)\n\n", freqNumber);
			PrintUsage();
			}
		for (i = 0; i < 6; i++)
			{
			if (Rmat[i] >= 0)
				{
				doGTR = YES;
				if (doHKY == YES)
					{
					#ifdef MPI
					if (rank==root)
						fprintf(stderr, "\n Introduce only a titv or only a R-matrix or only a N-Rmatrix ");
					#else
						fprintf(stderr, "\n Introduce only a titv or only a R-matrix or only a N-Rmatrix ");
					#endif
					PrintUsage();
					}
				break;
				}
			}
		for (i = 0; i < 12; i++)
			{
			if (NRmat[i] >= 0)
				{
				doGTnR = YES;
				if (doHKY == YES || doGTR == YES)
					{
					#ifdef MPI
					if (rank==root)
						fprintf(stderr, "\n Introduce only a titv or only a R-matrix or only a N-Rmatrix ");
					#else
						fprintf(stderr, "\n Introduce only a titv or only a R-matrix or only a N-Rmatrix ");
					#endif
					PrintUsage();
					}
				break;
				}
			}
		if (doHKY == NO && doGTR == NO && doGTnR == NO) /* default model */
			{
			doHKY = YES;
			titv = 0.5;
			}
		}
	if (migrationRate < 0) /* migration model */
		doMigration = NO;
	else
		doMigration = YES;
	if (numConvergDemes > 0)
		{
		doConvergDemes = YES;
		if (doMigration == NO)
			{
			#ifdef MPI
			if (rank==root)				
				fprintf (stderr, "PARAMETER ERROR: Using convergencies demes events it will need a migration model (-q) also\n\n");
			#else
				fprintf (stderr, "PARAMETER ERROR: Using convergencies demes events it will need a migration model (-q) also\n\n");
			#endif
			PrintUsage();
			}
		for (i=1; i<=numConvergDemes; i++)
			{
			if (deme_a_old[i] > (2*numPopulations-1))
				{
				#ifdef MPI
				if (rank==root)				
					fprintf (stderr, "PARAMETER ERROR: Bad number of deme in convergencies demes events (%d)\n\n", deme_a_old[i]);
				#else
					fprintf (stderr, "PARAMETER ERROR: Bad number of deme in convergencies demes events (%d)\n\n", deme_a_old[i]);
				#endif
				PrintUsage();
				}
			if (deme_b_old[i] > (2*numPopulations-1))
				{
				#ifdef MPI
				if (rank==root)	
					fprintf (stderr, "PARAMETER ERROR: Bad number of deme in convergencies demes events (%d)\n\n", deme_b_old[i]);
				#else
					fprintf (stderr, "PARAMETER ERROR: Bad number of deme in convergencies demes events (%d)\n\n", deme_b_old[i]);
				#endif
				PrintUsage();
				}
			}
		if (doMigration == NO)
			{
			#ifdef MPI
			if (rank==root)
				fprintf (stderr, "PARAMETER ERROR: Cannot have convergence of demes events (-%%) without migration model (-q)\n\n");
			#else
				fprintf (stderr, "PARAMETER ERROR: Cannot have convergence of demes events (-%%) without migration model (-q)\n\n");
			#endif
			PrintUsage();
			}
		}

	if (growthRate != 0)
		doExponential = YES;
	else
		doExponential = NO;
	if (numPeriods == 0)
		doDemographics = NO;
	else
		doDemographics = YES;
	if (doExponential == YES && doDemographics == YES)
		{
		#ifdef MPI
		if (rank==root)
			fprintf (stderr, "PARAMETER ERROR: Cannot have both demographics periods (-p) and other demographics (-g)\n\n");
		#else
			fprintf (stderr, "PARAMETER ERROR: Cannot have both demographics periods (-p) and other demographics (-g)\n\n");
		#endif
		PrintUsage();
		}
	
	/**********/
	/* MPI introducing */
#ifdef MPI
	if (rank==root)  // jmourino: Para que solo escriba por pantalla el procesador master			
		fprintf(stderr, "\n In MPI version, by moment, Recodon write an outputfile per replicate. In next versions, it will concatenate those output files\n ");
		
#endif

	/*********/
	
	start = clock();
#ifdef MPI													
      	if (rank==root)  // jmourino: Para que solo escriba por pantalla el procesador master (Tengo la duda de crear Results en el caso que no haya NFS, que lo tendrían que hacer todos)
	{
#endif
	PrintTitle(stderr);		/* write title */
 	PrintDate(stderr);		/* write date */
	mkdir("Results",S_IRWXU);	/* Create "Results" folder (with type S_IRWXU (read, write and execute)) */ 
	/*mkdir("Results",0);*/
#ifdef MPI													
	}
#endif

	#ifdef MAC					
		strcpy (dir,":Results:");	/* Copy the string in char variable dir = Results (char), is different mac vs windows */
	#else
		strcpy (dir,"Results/");
	#endif


	if (doM1 == NO && doM8 == NO && doM7 == NO && doOmegaRateHetCont == NO && doOmegaProb == NO && doOmegaRateHetDisc == NO && doPrintOmegasPerSitefiles == YES) /* For printing a warning */
		{
		doPrintOmegasPerSitefiles = NO;
		#ifdef MPI
		if (rank==root)
			fprintf(stderr, "\n * You are not using an appropriate model for writing the simulate omega per site \n ");
		#else
			fprintf(stderr, "\n * You are not using an appropriate model for writing the simulate omega per site \n ");
		#endif
		}


	/*fprintf(stderr, "\n originalSeed = %lu \n", originalSeed);*/

	/* Anadido Miguel: Aqui se abren los archivos cuando NO se usa MPI */		
#ifndef MPI			
		if (doSeparatedSequences == NO)
			{
			sprintf(File,"%s%s", dir, alignmentFile);			/* File = dir alignmentFile */
			if ((fpAlignment = fopen(File, "w")) == NULL)		/* if fpAlignment can't open */
				{
				fprintf(stderr, "Can't open %s.\n", File);
				exit(-1);
				}
			}
		if (doPrintBreakpoints == YES)	
			{
			sprintf(File,"%s%s", dir, breakpointFile);
			if ((fpBreakpoints = fopen(File, "w")) == NULL) /* if fpBreakpoints can't open */
				{
				fprintf(stderr, "Can't open %s.\n", File);
				exit(-1);
				}
			}
		if (doPrintTrees == YES)							/* if treeFile " " */
			{
			sprintf(File,"%s%s", dir, treeFile);
			if ((fpTrees = fopen(File, "w")) == NULL) 
				{
				fprintf(stderr, "Can't open %s.\n", File);
				exit(-1);
				}
			}
		if (doPrintTimes == YES)							/* if timesFile " " */
			{
			sprintf(File,"%s%s", dir, timesFile);
			if ((fpTimes = fopen(File, "w")) == NULL) 
				{
				fprintf(stderr, "Can't open %s.\n", File);
				exit(-1);
				}
			}
		if (doSettingsFile == YES)							/* For the programmer, if Settings file " " */
			{
			strcpy(settingsFile, "settings");
			sprintf(File,"%s%s", dir, settingsFile);
			if ((fpSettings = fopen(File, "w")) == NULL) 
				{
				fprintf(stderr, "Can't open %s.\n", File);
				exit(-1);
				}
			}
		
#endif
	kappa=(titv*(p_i[0]+p_i[2])*(p_i[1]+p_i[3]))/(p_i[0]*p_i[2] + p_i[1]*p_i[3]);
	/* titv = kappa * (p_i[0]*p_i[2] + p_i[1]*p_i[3]) / ((p_i[0]+p_i[2])*(p_i[1]+p_i[3])); */
	/*fprintf(stderr, "\n kappa = %lf.\n", kappa);*/
	if (doCodonModel == YES) /* numSites is the number of nucleotides or the number of codons */
		{
		numNuc = numSites*3;
		numCodons = numSites;
		}
	else
		numNuc = numSites;

	/* calculate expected number of recombination events per gene */
	a = 0;
	for (i=1; i<numSequences; i++)
		a += 1.0/(double)i;
	rho = 2.0 * Nscaling * N * recombinationRate * /*numSites*/numNuc;	/* rho is the recombination rate per gene per 4*N generations */
	expNumRE = rho * a;
	/*fprintf (stderr, "\n\nExpected number of rec events       =  %3.2f\n", expNumRE);*/

	expTMRCA = 2 * (1 - 1.0/numSequences)* Nscaling * N; /* multiply by 2N to get time in number of generations */
	expVarTMRCA = 0.0;
    for (i=2; i<= numSequences; i++)
       expVarTMRCA += 4.0 / (pow(i,2) * pow(i-1,2)) * pow(Nscaling * N,2) ;


	
	/* make sure dated tips are well specified  find the latest time, which will be time 0.0 going backwards in time 
		and sort samples */
	if (doDatedTips == YES)
		{
		sum = 0;
		for (i=0; i<numTipDates; i++)
			{
			sum += datedSample[i].size;
			if (datedSample[i].time > latestSamplingTime)
				latestSamplingTime = datedSample[i].time;
			}
		if (sum != numSequences)	
			{
			fprintf (stderr,"ERROR: the sum of the dated samples (%d) is not the same as the number of sequences (%d)", sum, numSequences);	
			exit(1);
			}
	
		/* set times in generation for dated tips; the youngest sample will be generation 0 */		
		for (i=0; i<numTipDates; i++)
			datedSample[i].time = (latestSamplingTime - datedSample[i].time) * generationTime;
				/* to get generations back (assumes generation time and sampling times are in the same unit, e.g., years) */

		/* sort the samples from older to younger */
		pass = 1;
		do 
			{
			sorted = YES;
			for (i=0; i<(numTipDates-pass); i++)
				{
				if (datedSample[i].time < datedSample[i+1].time)
					{		
					tempSample = datedSample[i+1];
					datedSample[i+1] = datedSample[i];
					datedSample[i]  = tempSample;
					sorted = NO;
					}
				}
			pass++;
			}	while (sorted == NO);

		i = 0;
		}
	




  	/* Set seed and spin wheels of pseudorandom number generator */
	/* seed = (unsigned int) clock();*/
	if (userSeed > 0)
		seed = userSeed;
	originalSeed = seed;
	for (i=0; i<10; i++)	
		RandomUniform(&seed);		/* function that generates random seed */  //jmourino: Esto tienen sentido realmente?  Por ahora,lo dejamos asi
	seedFirst = seed;
	/*fprintf(stderr, "\n seedFirst = %lu \n", seedFirst);*/
	
	/* default */
	cumNumRE = 0.0;
	cumNumCA = 0.0;
	cumNumMU = 0.0;
	cumNumREntc = 0.0;
	cumNumMIG = 0.0;
	cumNumCONV = 0.0;
	cumNumREbreakCod = 0.0;
	cumNumStopCodonREC = 0.0;
	cumNumEqual2 = 0.0;
	cumNumEqual1 = 0.0;
	cumNumDifCodSameAA = 0.0;
	cumNumDifCodDifAA = 0.0;
	cumNumNonSyn0 = 0.0;
	cumNumNonSyn1 = 0.0;
	cumNumNonSyn2 = 0.0;
	zeroRec = 0;
	countTMRCAReps = 0.0;
	counterTime = 0.0;
	



	/* preparing populations to migrations */
	if (doMigration == YES)
		{
		i = j = w = 0;
		numNodesInitPopul = (int *) calloc((numPopulations+1),(long) sizeof(int));
		if (!numNodesInitPopul)
			{
			#ifdef MPI
				fprintf (stderr, "%d: Could not allocate numNodesInitPopul (%lu)\n", rank, (numPopulations+1) * (long) sizeof(int)); /* casting with loing to avoid warnings in different OS */
			#else
				fprintf (stderr, "Could not allocate numNodesInitPopul (%lu)\n", (numPopulations+1) * (long) sizeof(int)); /* casting with loing to avoid warnings in different OS */
			#endif
			exit (-1);
			}
		
		for (i = 1; i <= numPopulations; i++)
			numNodesInitPopul[i] = initPopulation[i];
		
		for (i = 1; i <= numPopulations; i++)
			j = j + numNodesInitPopul[i];
		
		j = 0;
		for (i = 1; i <= numPopulations; i++)
			j = j + numNodesInitPopul[i];
		if (j != numSequences)
			{
			#ifdef MPI
				fprintf(stderr, "\n%d: Warning in program. Migration in main.",rank);
			#else
				fprintf(stderr, "\nWarning in program. Migration in main.");
			#endif
			exit (-1);
			}
		i = j = w = 0;
		}
	
	/* Variance memories */
	varEvent = (int *) calloc(numDataSets,(long) sizeof(int));
	if (!varEvent)
		{
		#ifdef MPI
			fprintf (stderr, "%d: Could not allocate varEvent (%lu bytes)\n", rank, numDataSets *(long) sizeof(int));
		#else
			fprintf (stderr, "Could not allocate varEvent (%lu bytes)\n", numDataSets *(long) sizeof(int));
		#endif
		exit (1);
		}
	varTimeGMRCA = (double *) calloc(numDataSets,(long) sizeof(double));
	if (!varTimeGMRCA)
		{
		#ifdef MPI
			fprintf (stderr, "%d: Could not allocate varTimeGMRCA (%lu bytes)\n", rank, numDataSets  * (long) sizeof (double));
		#else
			fprintf (stderr, "Could not allocate varTimeGMRCA (%lu bytes)\n", numDataSets  * (long) sizeof (double));
		#endif
		exit (1);
		}
	varTimeT = (double *) calloc(numDataSets,(long) sizeof(double));
	if (!varTimeT)
		{
		#ifdef MPI
			fprintf (stderr, "%d: Could not allocate varTimeT (%lu bytes)\n", rank, numDataSets  * (long) sizeof (double));
		#else
			fprintf (stderr, "Could not allocate varTimeT (%lu bytes)\n", numDataSets  * (long) sizeof (double));
		#endif
		exit (1);
		}
	
	if (doCodonModel == YES) /* making Qij for codon models */
		{
		if (doOmegaCat == NO && doOmegaRateHetCont == NO && doOmegaRateHetDisc == NO)
			buildCodonMatrix_Qij_Cijk ();
		
		if (doOmegaCat == YES)
			{
			QijOmegas = (Qij_OmegaCat *) calloc (numOmegaCat, sizeof(Qij_OmegaCat));
			if (!QijOmegas)
				{
				#ifdef MPI
					fprintf (stderr, "%d: Could not allocate QijOmegas (%lu)\n", rank, numOmegaCat  * (long) sizeof(Qij_OmegaCat));
				#else		
					fprintf (stderr, "Could not allocate QijOmegas (%lu)\n", numOmegaCat  * (long) sizeof(Qij_OmegaCat));
				#endif		
				exit (1);
				}
			
			for (j=1; j<=numOmegaCat; j++)
				{
				omega = omegaVal[j];
				buildCodonMatrix_Qij_Cijk ();
				
				for (w=0;w<NUMCOD;w++)
					QijOmegas[j-1].Root_C_cat[w] = Root_C[w];
				for (w=0;w<NUMCOD*NUMCOD*NUMCOD*NUMCOD;w++)
					QijOmegas[j-1].Cijk_C_cat[w] = Cijk_C[w];
				}			
			}
		
		if (doOmegaRateHetDisc == YES)
			{
			QijOmegas = (Qij_OmegaCat *) calloc (numOmegaCat, sizeof(Qij_OmegaCat));
			if (!QijOmegas)
				{
				#ifdef MPI
					fprintf (stderr, "%d: Could not allocate QijOmegas (%lu)\n", rank, numOmegaCat  * (long) sizeof(Qij_OmegaCat));
				#else
					fprintf (stderr, "Could not allocate QijOmegas (%lu)\n", numOmegaCat  * (long) sizeof(Qij_OmegaCat));
				#endif
				exit (1);
				}
				
			gammaRates =  (double*) calloc (numOmegaCat, sizeof (double));  
			if (gammaRates == NULL)
				{
				#ifdef MPI
						fprintf (stderr, "%d: Could not allocate %s (%lu bytes)", "gammaRates", rank, numOmegaCat  * (long) sizeof (double));
				#else
						fprintf (stderr, "Could not allocate %s (%lu bytes)", "gammaRates", numOmegaCat  * (long) sizeof (double));
				#endif
				exit(0);
				}
			
			omegaValGammaRate = (double *) calloc((numOmegaCat+1),(long) sizeof(double));
			if (!omegaValGammaRate)
				{
				fprintf (stderr, "PARAMETER ERROR: omegaValGammaRate. Could not allocate omega values of categories (%lu bytes)\n", numOmegaCat *(long) sizeof(double));
				exit (1);
				}
			
			/*fprintf (stderr, " categories %d\n", numOmegaCat);*/

			gammasCalculate (OmegaRateHet, numOmegaCat);
		
			for (j=1; j<=numOmegaCat; j++)
				{
				omegaGammaDisc = OmegaInit*gammaRates[j-1];

				omega = roundit(omegaGammaDisc,5);
				omegaValGammaRate[j] = omega;

				buildCodonMatrix_Qij_Cijk ();
				
				for (w=0;w<NUMCOD;w++)
					QijOmegas[j-1].Root_C_cat[w] = Root_C[w];
				for (w=0;w<NUMCOD*NUMCOD*NUMCOD*NUMCOD;w++)
					QijOmegas[j-1].Cijk_C_cat[w] = Cijk_C[w];
				}			
			}
		}
	
	
	if (doMRCAFile == YES) /* From MRCA file */ 
		{
		MRCAsequence = (char *) calloc((numNuc+1), sizeof(char)); 
		if (!MRCAsequence)
			{
			#ifdef MPI
				fprintf (stderr, "%d: Could not allocate MRCAsequence (%lu bytes)\n", rank, (numNuc+1)  * (long) sizeof(char));
			#else
				fprintf (stderr, "Could not allocate MRCAsequence (%lu bytes)\n", (numNuc+1)  * (long) sizeof(char));
			#endif
			exit (-1);
			}
					
		if ((fpMRCA = fopen(MRCAFile, "r+")) == NULL) /* if MRCAFile can't open */
			{
			#ifdef MPI
				fprintf(stderr, "%d: Can't open %s.\n", rank, MRCAFile);
			#else
				fprintf(stderr, "Can't open %s.\n", MRCAFile);
			#endif
			exit(-1);
			}
		else /* It can open MRCAFile */
			{
			i = 0;
			fscanf(fpMRCA, "%s", MRCAsequence); /* copy MRCA file into the MRCAsequence array */
			fclose(fpMRCA); 
			do 
				i++;
			while (MRCAsequence[i] == 'A' || MRCAsequence[i] == 'C' || MRCAsequence[i] == 'G' || MRCAsequence[i] == 'T');
		
			/* Active this whether we want to see the MRCA sequence from MRCA File */
			/*for (k = 0; k < i; k++)
				printf ("%c", MRCAsequence[k]);
			fprintf(stderr, "\nsites = %d\n", i);
			fprintf(stderr, "\n\n");*/
			/*printf ("%s", MRCAsequence);*/
			/*fprintf(stderr, "\n LongMRCA File = %d, Log Seq nuc = %d \n", i, numNuc);*/
			if (i != numNuc)
				{
				#ifdef MPI
					if (rank==root)		// Lo escribe el cero o todos? Si todos hay que quitar esta linea.			No entiendo muy bien que quieres decir. Si entra en este "if (i != numNuc)" el programa se para, communicando el error (hay dos par‡metros de entrada incompatibles). Aqui todavia deberia estar funcionando solo el master.
						{
						fprintf(stderr, "\nWarning in MRCA file (numNuc (from parametersFile) is different to numNuc (from MRCAFile)), verify:");
						fprintf(stderr, "\n1. numNuc (MRCAFile) = %d, numNuc (parametersFile) = %d", i, numNuc);
						fprintf(stderr, "\nor");
						fprintf(stderr, "\n2. A site of MRCA file is different to A, C, G or T");
						}
				#else
					{
					fprintf(stderr, "\nWarning in MRCA file (numNuc (from parametersFile) is different to numNuc (from MRCAFile)), verify:");
					fprintf(stderr, "\n1. numNuc (MRCAFile) = %d, numNuc (parametersFile) = %d", i, numNuc);
					fprintf(stderr, "\nor");
					fprintf(stderr, "\n2. A site of MRCA file is different to A, C, G or T");
					}
				#endif
				exit (-1);
				}
			}
		}
	i = j = w = 0;




// ****************************************************************************
// ************************* Cola de iteraciones ******************************
// ****************************************************************************

#ifdef MPI
        MPI_Barrier(comm);

// ****************************************************************************
        if (rank!=root) {   // CLIENT
// ****************************************************************************

      		MPI_Recv(&fila,1,MPI_INT,root,rank,comm,&status); //recibo numero de replicas
          do {
            MPI_Recv(&rep,1,MPI_INT,root,rank,comm,&status); /// recibo numero de replica

            // recibo datos (en principio no hace falta)
		
		if (p > numDataSets)
			{
			printf("Warning!. There are more processors (%d) than replicates (%d)", p, numDataSets);
			exit (-1);
			}
			
		// CALCULOS (las fila replicas)
		for (dataSetNum=rep; dataSetNum<rep+fila; dataSetNum++){
		printf("%d: calculo fila %d (%d de %d) en procesador %d\n",rank, dataSetNum, dataSetNum-rep+1, fila, rank);
#else

	  for (dataSetNum = 0; dataSetNum < numDataSets; dataSetNum++)		/**** FOR EACH REPLICATE ****/
		{

#endif		
		seed = seedFirst+dataSetNum+10; /* De esta forma la semilla sera la misma para secuencial que para MPI */ /* Variable seed among processors for MPI */
		do /* it will be a good replicate?. It will be running til a nice replicate (with the chosen number of recombinations) */
		{
		/* fprintf (stderr, "\n replicate = %d\n", dataSetNum+1); */
		
		#ifdef MPI
			{
			/*seed = seedFirst+dataSetNum+10;*/ 

			fprintf(stderr, "\n--------- IN MPI ----------\n");

			strcpy(screenFile, "screen");
			sprintf(File,"%s%s%05d", dir, screenFile, dataSetNum+1);					/* To write the screen information of the loop, for this replicate */
			if	((fpScreen = fopen(File, "w")) == NULL) /* if fpScreen can't open */
				{
				fprintf(stderr, "Can't open %s.\n", File);
				exit(-1);
				}



			if (doSeparatedSequences == NO)
				{
				sprintf(File,"%s%s%05d", dir, alignmentFile, dataSetNum+1);			/* File = dir alignmentFile */
				if ((fpAlignment = fopen(File, "w")) == NULL)		/* if fpAlignment can't open */
					{
					fprintf(fpScreen, "Can't open %s.\n", File);
					exit(-1);
					}
				}
			if (doPrintBreakpoints == YES)	
				{
				sprintf(File,"%s%s%05d", dir, breakpointFile, dataSetNum+1);
				if	((fpBreakpoints = fopen(File, "w")) == NULL) /* if fpBreakpoints can't open */
					{
					fprintf(fpScreen, "Can't open %s.\n", File);
					exit(-1);
					}
				}
			if (doPrintTrees == YES)							/* if treeFile " " */
				{
				sprintf(File,"%s%s%05d", dir, treeFile, dataSetNum+1);
				if ((fpTrees = fopen(File, "w")) == NULL) 
					{
					fprintf(fpScreen, "Can't open %s.\n", File);
					exit(-1);
					}
				}
			if (doPrintTimes == YES)							/* if timesFile " " */
				{
				fprintf(stderr, "\n - MPI.opening TIMES -");
				
				sprintf(File,"%s%s%05d", dir, timesFile, dataSetNum+1);
				if ((fpTimes = fopen(File, "w")) == NULL) 
					{
					fprintf(fpScreen, "Can't open %s.\n", File);
					exit(-1);
					}
				}
			}
		#endif

	/* Printing options in the parallelization area */ 
	#ifdef MPI
       fpmpi = fpScreen;
    #else
       fpmpi = stderr;
    #endif

		if (noisy == 0)
			{
			#ifdef MPI
			fprintf (fpmpi, "\nReplicate #%3d/%d", dataSetNum+1, numDataSets);
			#else
			fprintf (fpmpi, "\rReplicate #%3d/%d", dataSetNum+1, numDataSets);
			fflush (stdout);
			#endif
			}
		varEvent[dataSetNum] = 0;			
		varTimeGMRCA[dataSetNum] = varTimeT[dataSetNum] = 0.0;
		counterTimeInit = 0.0;
		actSegIndex = numRE = indNumRE = recNotToCount = numCA = numMU = numMIG = numCONV = 0;
		numMU_S = numMU_NS = 0;
		actualTGMRCA = countTMRCA = 0.0;
		numEqual2 = numEqual1 = numDifCodSameAA = numDifCodDifAA = numNonSyn0 = numNonSyn1 = numNonSyn2 = 0;		


		if (doBranchNetfiles == YES)							/* For the programmer, if just outfiles " " */
			{
			sprintf(File,"%s%s%05d", dir, BranchNetFilePrint, dataSetNum+1);
			if ((fpBranchNet = fopen(File, "w")) == NULL) 
				{
				fprintf(stderr, "Can't open %s.\n", File);
				exit(-1);
				}
			}
		
		NodesMRCAposit = (int *) calloc((numNuc+1),(long) sizeof(int)); /* this array will contain the labels of the MRCA nodes for each position */
		if (!NodesMRCAposit)
			{
			fprintf (fpmpi, "Could not allocate NodesMRCAposit (%lu)\n",(numNuc+1) * (long) sizeof(int)); /* casting with loing to avoid warnings in different OS */
			exit (-1);
			}


		if (noisy > 1)
			fprintf (fpmpi, "\n\n>> Start making coalescent fragments trees .. \n");
			
		MakeCoalescenceTree (numSequences, numSites, numNuc, N, recombinationRate, numPopulations, &seed);			/* Make the coalescence tree, topology */
		/* check whether we want to keep this replicate, It is possible to choose those replicates with a fix recombinations number*/
		
		//exit (10);
		
		/*matrix = (int *)calloc((2*(numSequences+1) * (numNuc+1)),(long) sizeof(int));
		if (!matrix)
			{
			fprintf (fpmpi, "Could not allocate matrix (%lu bytes)\n", (2*(numSequences+1) * (numNuc+1))  * (long) sizeof(int));
			exit (1);
			}*/

		matrix = (int *)calloc(((nextAvailable+1) * (numNuc+1)),(long) sizeof(int));		/* matrix is a vector with long = 2*(numSequences+1) * numSites */
		if (!matrix)
			{
			fprintf (fpmpi, "Could not allocate matrix (%lu bytes)\n", ((nextAvailable+1) * (numNuc+1))  * (long) sizeof(int));
			exit (1);
			}
		for (i = 0; i < (nextAvailable+1) * (numNuc+1); i++)
			matrix[i] = -1;	
		/*fprintf (fpmpi, "\nMatrix tiene num de labelnodes, tamanho = %d \n", nextAvailable+1);*/

		if (doCodonModel == YES)
			{
			matrixC = (int *)calloc(((nextAvailable+1) * (numSites+2)),(long) sizeof(int));		
			if (!matrixC)
				{
				fprintf (fpmpi, "Could not allocate matrixC (%lu bytes)\n", ((nextAvailable+1) * (numSites+2))  * (long) sizeof(int));
				exit (1);
				}
			/*matrixC = (int *)calloc((2*(numSequences+1) * (numSites+2)),(long) sizeof(int));		
			if (!matrixC)
				{
				fprintf (fpmpi, "Could not allocate matrixC (%lu bytes)\n", (2*(numSequences+1) * (numSites+2))  * (long) sizeof(int));
				exit (1);
				}*/
			for (i = 0; i < (nextAvailable+1) * (numSites+2); i++)
				matrixC[i] = -1;
			}


		if (doFixNumRecEvents == NO || numRE == fixedNumRecEvents)  /* This replicate is good */
			{
			if (noisy > 1)
				{
				//fprintf (fpmpi, "\n>> Simulating sequences for fragments trees .. \n\n");
				fprintf (fpmpi, "\n>> Simulating sequences  .. \n\n");
				}

			if (doSeparatedSequences == YES)
				{
				if (doPrintFASTA == NO && doPrintNEXUS == NO) /* phylip */
					{
					sprintf(File,"%s%s%05d", dir, alignmentFile, dataSetNum+1);			/* File = dir alignmentFile */
					if ((fpAlignment = fopen(File, "w")) == NULL)		/* if fpAlignment can't open */
						{
						fprintf(fpmpi, "Can't open %s.\n", File);
						exit(-1);
						}
					}
				if (doPrintFASTA == YES && doPrintNEXUS == NO) /* fasta */
					{
					sprintf(File,"%s%s%05d.fas", dir, alignmentFile, dataSetNum+1);			/* File = dir alignmentFile */
					if ((fpAlignment = fopen(File, "w")) == NULL)		/* if fpAlignment can't open */
						{
						fprintf(fpmpi, "Can't open %s.\n", File);
						exit(-1);
						}
					}
				if (doPrintNEXUS == YES && doPrintFASTA == NO) /* nexus */
					{
					sprintf(File,"%s%s%05d.nex", dir, alignmentFile, dataSetNum+1);			/* File = dir alignmentFile */
					if ((fpAlignment = fopen(File, "w")) == NULL)		/* if fpAlignment can't open */
						{
						fprintf(fpmpi, "Can't open %s.\n", File);
						exit(-1);
						}
					}
				}

			if (doPrintOmegasPerSitefiles == YES)
				{
				sprintf(File,"%s%s%05d", dir, OmegasPerSiteFile, dataSetNum+1);			/* File = dir  */
				if ((fpOmegasPerSitePrint = fopen(File, "w")) == NULL)		/* if it can't open */
					{
					fprintf(fpmpi, "Can't open %s.\n", File);
					exit(-1);
					}
				}

				
			if (doOutMRCAfiles == YES)							/* For the programmer, if just MRCA outfiles " " */
				{
				sprintf(File,"%s%s%05d", dir, MRCAfilePrint, dataSetNum+1); /* GMRCA */
				if ((fpMRCAprint = fopen(File, "w")) == NULL) 
					{
					fprintf(stderr, "Can't open %s.\n", File);
					exit(-1);
					}
				
				sprintf(File,"%s%s%05d", dir, ConcMRCAfilePrint, dataSetNum+1); /* catMRCA */
				if ((fpConcMRCAprint = fopen(File, "w")) == NULL) 
					{
					fprintf(stderr, "Can't open %s.\n", File);
					exit(-1);
					}
				
				if (doCodonModel == YES && doGMRCAsamp == YES)
					{
					sprintf(File,"%s%s%05d", dir, GMRCAancFilePrint, dataSetNum+1); /* ancestral material GMRCA */
					if ((fpGMRCAancPrint = fopen(File, "w")) == NULL) 
						{
						fprintf(stderr, "Can't open %s.\n", File);
						exit(-1);
						}

					}	

				}
			


									
			/* Generating sequences */
			do
				{
				if (doCodonModel == YES)
					EvolveSequenceOnTree_Codon (&seed, mutationRate, alpha, numNuc, indNumRE, arrayIndBreakpointsOrd, MRCAsequence, numOmegaCat, numSites);
				else
					{
					/*EvolveSequenceOnTree (&seed, mutationRate, kappa, alpha, p_i, numSites, arrayIndBreakpointsOrd, MRCAsequence);*/
					EvolveSequenceOnTree_NEW (&seed, mutationRate, kappa, alpha, p_i, numNuc, indNumRE, arrayIndBreakpointsOrd, MRCAsequence, numSites);
					}
				} while (doRepitEvol == YES);	/* it repit the evolve sequence while there are some stop codon by recombination in the codon*/	
			
			
			cumNumMU += numMU;
			cumNumMU_S += numMU_S;
			cumNumMU_NS += numMU_NS;			
			cumNumRE += numRE;
			cumNumREntc += recNotToCount;
			cumNumCA += numCA;
			cumNumREbreakCod += numREbreakCod;
			cumNumStopCodonREC += numStopCodonREC;
			cumNumEqual2 += numEqual2;
			cumNumEqual1 += numEqual1;
			cumNumDifCodSameAA += numDifCodSameAA;
			cumNumDifCodDifAA += numDifCodDifAA;
			cumNumNonSyn0 += numNonSyn0;
			cumNumNonSyn1 += numNonSyn1;
			cumNumNonSyn2 += numNonSyn2;
			if (doMigration == YES)
				{
				cumNumMIG += numMIG;
				if (doConvergDemes == YES)
					cumNumCONV += numCONV;
				}
			if (numRE == 0)
				zeroRec++;				/* zeroRep = number of replicates that has 0 recombination events */
			


			/* Average time to MRCA into this replicate with one or several trees */			
			if (numRE == 0) /* there aren't recombinations */
				{
				fragmentTMRCAfraction = (double *) calloc(1,(long) sizeof(double));
				if (!fragmentTMRCAfraction)
					{
					fprintf (fpmpi, "Could not allocate fragmentTMRCAfraction (%lu bytes)\n", 1 *(long) sizeof(double));
					exit (-1);
					}
				fragmentTMRCAfraction[0] = treeRootNodex[0]->time;
				countTMRCA = fragmentTMRCAfraction[0];
				}
			else
				{
				fragmentTMRCAfraction = (double *) calloc((indNumRE+1),(long) sizeof(double));
				if (!fragmentTMRCAfraction)
					{
					fprintf (fpmpi, "Could not allocate fragmentTMRCAfraction (%lu bytes)\n", (indNumRE+1) *(long) sizeof(double));
					exit (-1);
					}
				for (i = 0; i <= indNumRE; i++) /* it calculates the time to MRCA */
					{
					if (i == 0)
						fragmentTMRCAfraction[i] = treeRootNodex[0]->time*(arrayIndBreakpointsOrd[0]-1)/numNuc;
					if (i == indNumRE)
						fragmentTMRCAfraction[i] = treeRootNodex[indNumRE]->time*(numNuc-arrayIndBreakpointsOrd[indNumRE-1]+1)/numNuc;
					if (i > 0 && i < indNumRE)
						fragmentTMRCAfraction[i] = treeRootNodex[i]->time*(arrayIndBreakpointsOrd[i]-1 - arrayIndBreakpointsOrd[i-1] + 1)/numNuc;
					countTMRCA = countTMRCA + fragmentTMRCAfraction[i];
					}
				}
			free (fragmentTMRCAfraction);
			
			countTMRCAReps = countTMRCAReps + countTMRCA;
			counterTime = counterTime + counterTimeInit;
		
		
										/* output files */
			if (doPrintTrees == YES)
				PrintTreesSeg(dataSetNum, indNumRE, numNuc, arrayIndBreakpointsOrd);
			if (doPrintTimes == YES)
				PrintTimesSeg(dataSetNum, indNumRE, numNuc, arrayIndBreakpointsOrd);
			
			/*fprintf (fpmpi, "\n nextAvailable = %d", nextAvailable);*/
			if (doCodonModel == YES)
				{
				if (doPrintFASTA == NO && doPrintNEXUS == NO) /* phylip */
					{
					if (doPrintAncestralSequences == YES)
						PrintAncestralSequences_C (/*dataSetNum*/);
					else
						PrintSequences_C (/*dataSetNum*/);
					}
				if (doPrintFASTA == YES && doPrintNEXUS == NO) /* fasta */
					{
					if (doPrintAncestralSequences == YES)
						PrintAncestralSequences_C_FASTA (/*dataSetNum*/);
					else
						PrintSequences_C_FASTA (/*dataSetNum*/);
					}
				if (doPrintNEXUS == YES && doPrintFASTA == NO) /* nexus */
					{
					if (doPrintAncestralSequences == YES)
						PrintAncestralSequences_C_NEXUS (/*dataSetNum*/);
					else
						PrintSequences_C_NEXUS (/*dataSetNum*/);
					}

				if (doOutMRCAfiles == YES)
					{
					PrintOutMRCAFiles_C(/*dataSetNum*/);
					PrintOutMRCAFiles_C_Conc(/*dataSetNum*/);
					if (doGMRCAsamp == YES)
						PrintOutGMRCAFiles_Codon_AncestralMat();
					}
				}
			else /* Print from nucleotide models */
				{
				if (doPrintFASTA == NO && doPrintNEXUS == NO) /* phylip */
					{
					if (doPrintAncestralSequences == YES)
						PrintAncestralSequences (/*dataSetNum*/);
					else
						PrintSequences (/*dataSetNum*/);
					}
				if (doPrintFASTA == YES && doPrintNEXUS == NO) /* fasta */
					{
					if (doPrintAncestralSequences == YES)
						PrintAncestralSequences_FASTA (/*dataSetNum*/);
					else
						PrintSequences_FASTA (/*dataSetNum*/);
					}
				if (doPrintNEXUS == YES && doPrintFASTA == NO) /* nexus */
					{
					if (doPrintAncestralSequences == YES)
						PrintAncestralSequences_NEXUS (/*dataSetNum*/);
					else
						PrintSequences_NEXUS (/*dataSetNum*/);
					}				

				if (doOutMRCAfiles == YES)
					{
					PrintOutMRCAFiles(/*dataSetNum*/);
					PrintOutMRCAFiles_Conc();
					}
				}

			

			varEvent[dataSetNum] = numCA+numRE+numMIG+numCONV;
			#ifdef HUDSON_UNITS
				varTimeGMRCA[dataSetNum] = actualTGMRCA/(2*Nscaling*N); /* Hudson units (/2*Nscaling*N) */
			#else
				varTimeGMRCA[dataSetNum] = actualTGMRCA;
			#endif
			varTimeT[dataSetNum] = countTMRCA;
				
			if (noisy > 0)
				{
				fprintf (fpmpi, "\nData set %d",dataSetNum+1);
				fprintf (fpmpi, "\nNumber of total recombination events =  %d", numRE);
				if (doCountsForExpNumRec == YES)
					fprintf (fpmpi, "\nNumber of count recombination events =  %d", numRE - recNotToCount);
				if ((numRE > 0) && (doCodonModel == YES))
					{
					fprintf (fpmpi, "\n Countable number of recombinations breaking codons           =  %d", numREbreakCod);
					fprintf (fpmpi, "\n Countable number of stop codons generated/restarted by InRec =  %d", numStopCodonREC);
					
					/*fprintf (fpmpi, "\n numEqual2           =  %d", numEqual2);
					fprintf (fpmpi, "\n numEqual1           =  %d", numEqual1);
					fprintf (fpmpi, "\n numDifCodSameAA     =  %d", numDifCodSameAA);
					fprintf (fpmpi, "\n numDifCodDifAA      =  %d", numDifCodDifAA);*/
					fprintf (fpmpi, "\n Countable number of recombinations with 0 nonsynonymous changes =  %d (%3.2f%%)", numNonSyn0, 1.00*numNonSyn0*100/numRE);
					fprintf (fpmpi, "\n Countable number of recombinations with 1 nonsynonymous changes =  %d (%3.2f%%)", numNonSyn1, 1.00*numNonSyn1*100/numRE);
					fprintf (fpmpi, "\n Countable number of recombinations with 2 nonsynonymous changes =  %d (%3.2f%%)", numNonSyn2, 1.00*numNonSyn2*100/numRE);
					/*numEqual2, numEqual1, numDifCodSameAA, numDifCodDifAA; numNonSyn*/
					}
				fprintf (fpmpi, "\nNumber of coalescence events  =  %d", numCA);
				if (doMigration == YES)
					{
					fprintf (fpmpi, "\nNumber of migration events    =  %d", numMIG);
					if (doConvergDemes == YES)
						fprintf (fpmpi, "\nNumber of convergence events  =  %d", numCONV);
					}
				fprintf (fpmpi, "\nNumber of mutational events   =  %d", numMU);
				if (doCodonModel == YES)
					{
					fprintf (fpmpi, "\n  synonymous changes =  %d (%3.2f%%)", numMU_S, 1.00*numMU_S*100/numMU);
					fprintf (fpmpi, "\n  nonsynonymous changes =  %d (%3.2f%%)", numMU_NS, 1.00*numMU_NS*100/numMU);
					}

				if (numRE > 0)
					{
					fprintf (fpmpi, "\nTime to GMRCA                 =  %3.2f", actualTGMRCA);
					fprintf (fpmpi, "\nAverage Time to MRCA          =  %3.2f\n\n", countTMRCA);
					}
				else	
					fprintf (fpmpi, "\nTime to MRCA                  =  %3.2f\n\n", countTMRCA);
				/*fprintf (fpmpi, "\n");*/		
				}
			doBadReplicate = NO; /* This is a good replicate */
			
			/* Breakpoints */
			if (doPrintBreakpoints == YES)
				{
				fprintf (fpBreakpoints, "replicate %d:\t", dataSetNum+1);
				for (i=0; i<numRE; i++)
					fprintf (fpBreakpoints, " %d", breakpoint[i]);
				fprintf (fpBreakpoints, "\n");
				}
			/* Only when it does NOT run MPI */
			#ifndef MPI
			if (doSeparatedSequences == YES)
				fclose(fpAlignment);
			
			if (doOutMRCAfiles == YES)
				{
				fclose (fpMRCAprint);
				fclose (fpConcMRCAprint);

				if (doCodonModel == YES && doGMRCAsamp == YES)
					{
					fclose(fpGMRCAancPrint);
					}	

				}
			

			if (doBranchNetfiles == YES)
				fclose (fpBranchNet);
			if (doPrintOmegasPerSitefiles == YES)
				fclose (fpOmegasPerSitePrint);
			#endif
			}
		else			/* we don't want to keep this replicate */
			{
			if (noisy > 0)
				{
				fprintf (fpmpi, "\nData set discarded ");
				fprintf (fpmpi, "(rec events =  %d)\n", numRE);	
				}

			/* Only when it does NOT run MPI */
			/*#ifndef MPI
			dataSetNum--;
			#endif*/
			doBadReplicate = YES; /* This is a bad replicate, it will try again */
			if (doBranchNetfiles == YES)
				fclose (fpBranchNet);
			}

		/* free memory */
		/*for (mmm=0; mmm<nextAvailable; mmm++)*/
		for (mmm=0; mmm<numNodes; mmm++)
			{
			free (nodes[mmm].SitesNonAncHere);
			}
		/* free (SitesNonAncHere); */

		free (nodex);
		free (treeRootNodex);
		free (arrayIndBreakpointsOrd);
		free (treeRootInit);
		free (segments);
		free (nodes);	
		/* free memory */		
		free (breakpoint);
		free (matrix);
		if (doCodonModel == YES)
			free (matrixC);
		free (NodesMRCAposit);

		
#ifdef MPI
			fprintf(stderr, "\n--------- IN MPI ----------\n");
			if (doSeparatedSequences == NO)
				{
				fclose(fpAlignment);
				fprintf(fpmpi, "\n Sequences printed to file \"%s\"", alignmentFile);
				}
			else
				fprintf(fpmpi, "\n Sequences printed to files \"%s\"", alignmentFile);

			if (doPrintTrees == YES)
				{
				fprintf(fpmpi, "\n Trees printed to file \"%s\"", treeFile);
				fclose(fpTrees);
				}
			if (doPrintTimes == YES)
				{
				fprintf(stderr, "\n MPI.closing TIMES ");
				fprintf(fpmpi, "\n Times printed to file \"%s\"", timesFile);
				fclose(fpTimes);
				}
			if (doPrintBreakpoints == YES)
				{
				fprintf(fpmpi, "\n Breakpoints printed to file \"%s\"", breakpointFile);
				fclose(fpBreakpoints);
				}
			fclose(fpScreen);
#endif

		
		} while (doBadReplicate == YES);
		/*fprintf (stderr, "\n FINAL_replicate = %d\n", dataSetNum+1);*/
		
		} /*** END OF REPLICATES ***/  // Vale para ambos (tanto como si es MPI como si no)			Ok.

#ifdef MPI

            MPI_Isend(&nada,1,MPI_INT,root,rank,comm,&request); // digo que termine
            MPI_Wait(&request,&status);

            // envio datos
		printf("%d: * envio a %d %d filas desde %d\n", rank, root, fila, rep);
   	    MPI_Send(&varEvent[rep],fila,MPI_INT,root,rank,comm);
   	    MPI_Send(&varTimeGMRCA[rep],fila,MPI_DOUBLE,root,rank,comm);
   	    MPI_Send(&varTimeT[rep],fila,MPI_DOUBLE,root,rank,comm);

	    MPI_Recv(&fila,1,MPI_INT,root,rank,comm,&status); //recibo numero de replicas

          } while (fila!=0);

	}	// END CLIENT


#endif


#ifdef MPI

// ***********************************************************************************
         else {     // SERVER
// ***********************************************************************************

          for (fila=0; fila<(lcola*(p-2)+1); fila=fila+lcola) {    // jmourino: Hay que controlar que no haya mas procesadores que réplicas			Ok.
	    printf("mando fila=%d \n",fila);
            nodo=(fila)/lcola+1;
   	    MPI_Send(&lcola,1,MPI_INT,nodo,nodo,comm); // mando numero de replicas

            MPI_Send(&fila,1,MPI_INT,nodo,nodo,comm); // mando numero de replica

            // envio datos iniciales

            MPI_Irecv(&nada, 1, MPI_INT, nodo, nodo,comm, &request);  // recepcion asincrona para cuando acabe)

            requests[nodo-1]=request;
            filas[nodo-1]=fila;
          } //for

          for (fila=((p-1)*lcola); fila<(numDataSets-lcola+1); fila=fila+lcola)  { // mientras queden replicas
 	    MPI_Waitany(p-1, requests, &indice, &status); // cuando acabe alguna 
	    printf("acabo en %d la fila %d \n",indice+1,filas[indice]);

            //recibo datos
		printf("%d: * recibo desde %d, %d filas desde %d\n", rank, indice+1, lcola, filas[indice]);
   	    MPI_Recv(&varEvent[filas[indice]],lcola,MPI_INT,indice+1,indice+1,comm,&status);
   	    MPI_Recv(&varTimeGMRCA[filas[indice]],lcola,MPI_DOUBLE,indice+1,indice+1,comm,&status);
   	    MPI_Recv(&varTimeT[filas[indice]],lcola,MPI_DOUBLE,indice+1,indice+1,comm,&status);

   	    MPI_Send(&lcola,1,MPI_INT,indice+1,indice+1,comm); // mando numero de replicas
            MPI_Send(&fila,1,MPI_INT,indice+1,indice+1,comm); // mando numero de replica
	    printf("mando fila=%d \n",fila);

            // mando datos

            MPI_Irecv(&nada, 1, MPI_INT, indice+1, indice+1, comm, &request);  // recepcion asincrona para cuando acabe
            requests[indice]=request;
            filas[indice]=fila;
          } //for

//    		R E S T O 

       resto=numDataSets-(fila);
       printf("resto=%d\n",resto);
       flag=0;
       if (resto > 0) {
	 flag=1;
         printf("resto=%d\n",resto);
	 MPI_Waitany(p-1, requests, &indice, &status); // cuando acabe alguna 
	 printf("resto: acabo en %d la fila %d \n",indice+1,filas[indice]);

            // recibo datos
		printf("%d: * recibo desde %d, %d filas desde %d\n", rank, indice+1, lcola, filas[indice]);
   	    MPI_Recv(&varEvent[filas[indice]],lcola,MPI_INT,indice+1,indice+1,comm,&status);
   	    MPI_Recv(&varTimeGMRCA[filas[indice]],lcola,MPI_DOUBLE,indice+1,indice+1,comm,&status);
   	    MPI_Recv(&varTimeT[filas[indice]],lcola,MPI_DOUBLE,indice+1,indice+1,comm,&status);


  	  MPI_Send(&resto,1,MPI_INT,indice+1,indice+1,comm); // mando resto de replicas
          MPI_Send(&fila,1,MPI_INT,indice+1,indice+1,comm); // mando numero de replica
          printf("resto: mando fila=%d \n",fila);

            // mando datos

          MPI_Irecv(&nada, 1, MPI_INT, indice+1, indice+1,comm, &request);  // recepcion asincrona para cuando acabe
          filas[indice]=fila;
          MPI_Wait(&request,&status);
 	  printf("rresto: acabo en %d la fila %d \n",indice+1,filas[indice]);
        
            // recibo datos
		printf("%d: * recibo desde %d, %d filas desde %d\n", rank, indice+1, resto, fila);
  	    MPI_Recv(&varEvent[fila],resto,MPI_INT,indice+1,indice+1,comm,&status);
   	    MPI_Recv(&varTimeGMRCA[fila],resto,MPI_DOUBLE,indice+1,indice+1,comm,&status);
   	    MPI_Recv(&varTimeT[fila],resto,MPI_DOUBLE,indice+1,indice+1,comm,&status);

  	MPI_Send(&zero,1,MPI_INT,indice+1,indice+1,comm); // mando cero para terminar
       }


//    Recojo las que quedan

          for (fila=0; fila < (p-1-flag); fila ++)  { 	// recojo los que quedan
	    MPI_Waitany(p-1, requests, &indice, &status); // cuando acabe alguna 
	    printf("acabo en %d la fila %d \n",indice+1,filas[indice]);

            // recibo datos
		printf("%d: * recibo desde %d, %d filas desde %d\n", rank, indice+1, lcola, filas[indice]);
   	    MPI_Recv(&varEvent[filas[indice]],lcola,MPI_INT,indice+1,indice+1,comm,&status);
   	    MPI_Recv(&varTimeGMRCA[filas[indice]],lcola,MPI_DOUBLE,indice+1,indice+1,comm,&status);
   	    MPI_Recv(&varTimeT[filas[indice]],lcola,MPI_DOUBLE,indice+1,indice+1,comm,&status);

    	    MPI_Send(&zero,1,MPI_INT,indice+1,indice+1,comm); // mando cero para terminar
          }

        } // END SERVER

#endif			

// ****************************************************************************
// *********************** Fin Cola de iteraciones ****************************
// ****************************************************************************


	
	/* 11 VARIABLES GLOBALES PENDIENTES DE ARREGLAR (8 escalares y 3 vectores) */
		/* data average */
#ifdef MPI

int				tzeroRec;
double 			tcumNumRE, tcumNumREntc, tcumNumCA, tcumNumMU, tcumNumMIG, tcountTMRCAReps, tcounterTime, tcumNumCONV, tcumNumREbreakCod, tcumNumStopCodonREC, tcumNumEqual2, tcumNumEqual1, tcumNumDifCodSameAA, tcumNumDifCodDifAA, tcumNumNonSyn0, tcumNumNonSyn1, tcumNumNonSyn2; 


//		printf("\n%d: Valores antes: %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f \n",rank,cumNumCA,cumNumRE,cumNumREntc,cumNumMU,cumNumMIG,zeroRec,countTMRCAReps,counterTime, cumNumCONV, cumNumREbreakCod, cumNumStopCodonREC, cumNumEqual2, cumNumEqual1, cumNumDifCodSameAA, cumNumDifCodDifAA, cumNumNonSyn0, cumNumNonSyn1, cumNumNonSyn2); 
		MPI_Reduce(&cumNumCA,&tcumNumCA,1, MPI_DOUBLE,MPI_SUM,root,comm);
		MPI_Reduce(&cumNumRE,&tcumNumRE,1, MPI_DOUBLE,MPI_SUM,root,comm);
		MPI_Reduce(&cumNumREntc,&tcumNumREntc,1, MPI_DOUBLE,MPI_SUM,root,comm); 
		MPI_Reduce(&cumNumMU,&tcumNumMU,1, MPI_DOUBLE,MPI_SUM,root,comm);
		MPI_Reduce(&cumNumMIG,&tcumNumMIG,1, MPI_DOUBLE,MPI_SUM,root,comm);
		MPI_Reduce(&cumNumCONV,&tcumNumCONV,1, MPI_DOUBLE,MPI_SUM,root,comm);
		MPI_Reduce(&cumNumREbreakCod,&tcumNumREbreakCod,1, MPI_DOUBLE,MPI_SUM,root,comm);
		MPI_Reduce(&cumNumStopCodonREC,&tcumNumStopCodonREC,1, MPI_DOUBLE,MPI_SUM,root,comm);
		MPI_Reduce(&cumNumEqual2,&tcumNumEqual2,1, MPI_DOUBLE,MPI_SUM,root,comm);
		MPI_Reduce(&cumNumEqual1,&tcumNumEqual1,1, MPI_DOUBLE,MPI_SUM,root,comm);
		MPI_Reduce(&cumNumDifCodSameAA,&tcumNumDifCodSameAA,1, MPI_DOUBLE,MPI_SUM,root,comm);
		MPI_Reduce(&cumNumDifCodDifAA,&tcumNumDifCodDifAA,1, MPI_DOUBLE,MPI_SUM,root,comm);
		MPI_Reduce(&cumNumNonSyn0,&tcumNumNonSyn0,1, MPI_DOUBLE,MPI_SUM,root,comm);
		MPI_Reduce(&cumNumNonSyn1,&tcumNumNonSyn1,1, MPI_DOUBLE,MPI_SUM,root,comm);
		MPI_Reduce(&cumNumNonSyn2,&tcumNumNonSyn2,1, MPI_DOUBLE,MPI_SUM,root,comm);
		
		MPI_Reduce(&zeroRec,&tzeroRec,1, MPI_INT,MPI_SUM,root,comm);
		MPI_Reduce(&countTMRCAReps,&tcountTMRCAReps,1, MPI_DOUBLE,MPI_SUM,root,comm);
		MPI_Reduce(&counterTime,&tcounterTime,1, MPI_DOUBLE,MPI_SUM,root,comm);

		if (rank==root) {
			cumNumCA = tcumNumCA;
			cumNumRE = tcumNumRE;
			cumNumREntc = tcumNumREntc; 
			cumNumMU = tcumNumMU;
			cumNumMIG = tcumNumMIG;
			cumNumCONV = tcumNumCONV;
			cumNumREbreakCod = tcumNumREbreakCod;
			cumNumStopCodonREC = tcumNumStopCodonREC;
			cumNumEqual2 = tcumNumEqual2;
			cumNumEqual1 = tcumNumEqual1;
			cumNumDifCodSameAA = tcumNumDifCodSameAA;
			cumNumDifCodDifAA = tcumNumDifCodDifAA;
			cumNumNonSyn0 = tcumNumNonSyn0;
			cumNumNonSyn1 = tcumNumNonSyn1;
			cumNumNonSyn2 = tcumNumNonSyn2;
			zeroRec = tzeroRec;
			countTMRCAReps = tcountTMRCAReps;
			counterTime = tcounterTime;

		}
	
//	printf("\n%d: Valores despues: %f,%f,%f,%f,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f \n",rank,cumNumCA,cumNumRE,cumNumREntc,cumNumMU,cumNumMIG,zeroRec,countTMRCAReps,counterTime,cumNumCONV, cumNumREbreakCod, cumNumStopCodonREC, cumNumEqual2, cumNumEqual1, cumNumDifCodSameAA, cumNumDifCodDifAA, cumNumNonSyn0, cumNumNonSyn1, cumNumNonSyn2); 
#endif		

		
//		for (i = 0; i < numDataSets; i++)
//			{
//		printf("\n%d: Vectores: %d, %f,%f\n",varEvent[i],varTimeGMRCA[i],varTimeT[i] );
//			varEvent[i] = 0;
//			varTimeGMRCA[i] = 0.0;
//			varTimeT[i] = 0.0;
//			}
		
		
				
		meanNumRE = cumNumRE / numDataSets;				
		meanNumREtc = (cumNumRE - cumNumREntc) / numDataSets;
		meanNumCA = cumNumCA / numDataSets;
		meanNumMU = cumNumMU / numDataSets;
		meanNumMU_S = cumNumMU_S / numDataSets;
		meanNumMU_NS = cumNumMU_NS / numDataSets;

		if (doMigration == YES)
			{
			meanNumMIG = cumNumMIG / numDataSets;
			if (doConvergDemes == YES)
				meanNumCONV = cumNumCONV / numDataSets;
			}
		meanNumREbreakCod = cumNumREbreakCod / numDataSets;
		meanNumStopCodonREC = cumNumStopCodonREC / numDataSets;		
		meanNumEqual2 = cumNumEqual2 / numDataSets;	
		meanNumEqual1 = cumNumEqual1 / numDataSets;
		meanNumDifCodSameAA = cumNumDifCodSameAA / numDataSets;
		meanNumDifCodDifAA = cumNumDifCodDifAA / numDataSets;
		meanNumNonSyn0 = cumNumNonSyn0 /numDataSets;
		meanNumNonSyn1 = cumNumNonSyn1 /numDataSets;
		meanNumNonSyn2 = cumNumNonSyn2 /numDataSets;

#ifdef MPI	 
		fprintf(stderr, "\n--------- IN MPI ----------\n");  
		if (noisy > 0)	/* De esta forma escribira la funcion para MPI por la pantalla del procesador master */
	        if (rank==root)  // jmourino: Para que solo escriba por pantalla el procesador master
				PrintRunSettings (stdout, originalSeed);		/* this function writes the value of all variables */
#else
		if (noisy > 0) /* De esta forma escribira la funcion para Secuencial por la pantalla de ejecucion */
			PrintRunSettings (stdout, originalSeed);		/* this function writes the value of all variables */
		if (doSettingsFile == YES)
			fclose(fpSettings);
#endif
				
		




#ifdef MPI													
      	if (rank==root)  // jmourino: Para que solo escriba por pantalla el procesador master			
		{
#endif
		fprintf(stderr, "\n\n\n*** Simulations finished ***");
		fprintf(stderr, "\n\nOutput files are in folder \"Results\":");
#ifdef MPI
		}	
#endif
			
												
																														
	

	if (doMRCAFile == YES) 
		free (MRCAsequence);
	free(varEvent);
	free(varTimeGMRCA);
	free(varTimeT);
	if (doDemographics == YES)
		{
		free (Nbegin);
		free (Nend);
		free (cumDuration);
		free (periodGrowth);
		}
	if (doMigration == YES) 
		{
		free (initPopulation);
		free (numNodesInitPopul);
		}
	if (doCodonModel == YES)
		{
		if (doOmegaCat == YES)
			{
			free (QijOmegas);
			free (omegaVal); 
			free (omegaProb); 
			}
		if (doOmegaRateHetDisc == YES)
			{
			free (QijOmegas);
			free (gammaRates);
			free (omegaValGammaRate);
			}
		}
	if (doConvergDemes == YES)
		{
		free (convDemTimes_old);
		free (deme_a_old);
		free (deme_b_old);
		}	

	if (doDatedTips == YES)
		{
		for (j=0; j<numTipDates; j++)
			{
			free (datedSample[j].member);
			}
		free (datedSample);
		}



#ifdef MPI													
      	if (rank==root)  // jmourino: Para que solo escriba por pantalla el procesador master			
	{
#endif


	/* closing files */
#ifndef MPI
		
		if (doPrintTimes == YES)
			{			
			fprintf(stderr, "\n Times printed to file \"%s\"", timesFile);
			/*fclose(fpTimes);*/
			}
		

		if (doSeparatedSequences == NO)
			{
			fprintf(stderr, "\n Sequences printed to file \"%s\"", alignmentFile);
			fclose(fpAlignment);
			}
		else
			fprintf(stderr, "\n Sequences printed to files \"%s\"", alignmentFile);

		if (doPrintTrees == YES)
			{
			fprintf(stderr, "\n Trees printed to file \"%s\"", treeFile);
			fclose(fpTrees);
			}
		
		if (doPrintBreakpoints == YES)
			{
			fprintf(stderr, "\n Breakpoints printed to file \"%s\"", breakpointFile);
			fclose(fpBreakpoints);
			}

		if (doSettingsFile == YES)
			fprintf(stderr, "\n Settings printed to file \"%s\"", settingsFile);
		if (doOutMRCAfiles == YES)
			{
			fprintf(stderr, "\n MRCAs printed to files \"%s\"", ConcMRCAfilePrint);
			fprintf(stderr, "\n GMRCAs printed to files \"%s\"", MRCAfilePrint);
			if (doCodonModel == YES)
				{
				fprintf(stderr, "\n GMRCAs of ancestral material were printed to files \"%s\"", GMRCAancFilePrint);
				}	
			}
		if (doBranchNetfiles == YES)
			fprintf(stderr, "\n Branches of the net printed to files \"%s\"", BranchNetFilePrint);
		if (doPrintOmegasPerSitefiles == YES)
			fprintf(stderr, "\n Simulated Omegas per site printed to files \"%s\"", OmegasPerSiteFile);
				
#endif



	/* execution time */
	secs = (double)(clock() - start) / CLOCKS_PER_SEC;
	fprintf(stderr, "\n\n\n___________________________________________________________________");
	fprintf(stderr, "\nTime processing: %G seconds\n", secs);
	fprintf(stderr, "\nIf you need help type '-?' in the command line of the program\n\n");
#ifdef MPI
	}	
#endif


#ifdef MPI

// Finalizo MPI
MPI_Finalize();	 	

#endif

	return (1);
}








/************************* MakeCoalescenceTree ************************/
/* Builds a genealogy for each site under the coalescent with or without
 recombination and migration */ /* this function go by events */ 
 
void MakeCoalescenceTree (int numSequences, int numSites, int numNuc, int N, double recombinationRate, int numPopulations, long int *seed)
	{
	int			e, c, d, i, j, w, k, a, b, aa, bb, aaa, bbb, ss, step, position, overFirst, overEnd, sss, *activeGametes, isCoalescence, whichInd, 
				whichSite, hasPassedBreakPoint,
				firstInd, secondInd, newInd, firstHalf, secondHalf, foundSuperflousNode, eventNum, numActiveGametes, 
				probRecIndividual, legalBreakpoint, cum_gi, period, isRecombination, isMigration, whichDeme, arrivedDeme, currentBigDeme, currentDemesNumber;
	int			minInit_qq, minInit_pp, maxEnd_pp, maxEnd_qq;
	int			actNumSegments, numTotalSegments, labelNodes;
	double		rateCA, rateRE, rate, rateMIG, timeRE, timeCA, timeMIG, currentTime, eventTime, timeCONV, variable1, variable2;
	TreeNode	*p, *q, *r;
	TreeSegment *s, *n, *m, *z;
	TreeNodex	*f, *g, *h;
	long int 	Gi;
	int			*gi; 
	int			*arrayIndBreakpoints, *numSegTrees, *coalVectorCountStarts, *coalVectorCountEnds;
	int			*coalEqualSegInit_p, *coalEqualSegInit_q, *coalEqualSegEnd_p, *coalEqualSegEnd_q;
	int			*initialVector_pp, *initialVector_qq, *endVector_pp, *endVector_qq;
	int			*endsVectorRec, *startsVectorRec;
	int			sizeNode, sizeNode_p, sizeNode_q, out;
	int			numActNodex, countNumTrees;
	double		ran;
	double		*cumPopulPart;
	int			many, memoryBreakp, maxSegNode;
	int			nodeValue;
	int			*cumInitPopul, *numParcialActiveGametes, *GiPartial, *ArrivedDemesOptions, *CurrentDemesState, *stud;
	double		*rateREpartial, *rateCApartial, *rateMIGpartial, *ratePartial, *cumPopulTase;
	int			currentSample, saveThis;	
	int			doBreakpBroken, LeftLess, LeftHigh, RightLess, RightHigh, LeftLess2, RightHigh2, mmm, stateHere_P, stateHere_Q, sigue, Ok_SMRCA_Codon, AncGMRCA_obtained;	

	
	/* defaults */
	isCoalescence = NO;
	isRecombination = NO;
	isMigration = NO;
	minInit_qq = minInit_pp = maxEnd_pp = maxEnd_qq = newInd = whichDeme = labelNodes = ss = saveThis = 0;
	numNodex = 0;
	maxSegNode = 1;
	currentSample = -1;
	m = NULL;
	n = NULL;
	ratePartial = NULL;
	rateREpartial = NULL;
	rateCApartial = NULL;
	rateMIGpartial = NULL;
	cumPopulTase = NULL;
	GiPartial = NULL;
	numParcialActiveGametes = NULL;
	cumInitPopul = NULL;
	ArrivedDemesOptions = NULL;
	CurrentDemesState = NULL;
	numNodes = 5000;
	memoryBreakp = /*numSites*/ numNuc;
	timeCONV = 0;
	currentBigDeme = 0;
	eventTime = 0.0;
	variable1 = 0.0;
	variable2 = 0.0;
	numREbreakCod = 0;
	numStopCodonREC = 0;
	numNetLabelPrint = 0;
	doBreakpBroken = NO;
	LeftLess = LeftHigh = RightLess = RightHigh = LeftLess2 = RightHigh2 = -1;
	stateHere_P = stateHere_Q = -2;
	mmm = sigue = Ok_SMRCA_Codon = 0;
	AncGMRCA_obtained = NO;
	overFirst = overEnd = sizeNode_p = sizeNode_q = sizeNode = out = numActNodex = many = c = 0;

	currentDemesNumber = numPopulations;

	/** Initial memories **/
	/*distance = numSites*2+2; is the max segments number in a node */
	if (recombinationRate == 0.00)
		distance = 1;
	else /* recombinationRate > 0 */
		{
		if (/*numSites*/numNuc < 400)
			distance = /*numSites*/numNuc*4;
		if (/*numSites*/numNuc >= 400 && /*numSites*/numNuc <= 1000)
			distance = /*numSites*/numNuc*3;
		if (/*numSites*/numNuc > 1000)
			distance = /*numSites*/numNuc*1.5;
		if (/*numSites*/numNuc > 10000)
			distance = /*numSites*/numNuc*0.5;
		if (/*numSites*/numNuc > 100000)
			distance = /*numSites*/numNuc*0.005;
		}

	/* initial segments from numSites and rec rate */
	if (recombinationRate == 0.00) 
		maxSegNode = 1;
	else /* recombinationRate > 0 */
		{
		if (/*numSites*/numNuc < 400)
			{
			if (recombinationRate <= 0.00001)
				maxSegNode = /*numSites*/numNuc*2;
			if (recombinationRate <= 0.0001 && recombinationRate > 0.00001)
				maxSegNode = /*numSites*/numNuc*4;
			if (recombinationRate > 0.0001)
				maxSegNode = /*numSites*/numNuc*5;
			}
		else if (/*numSites*/numNuc >= 400 && /*numSites*/numNuc <= 1000)
			{
			maxSegNode = /*numSites*/numNuc*3;
			if (recombinationRate <= 0.00001)
				maxSegNode = /*numSites*/numNuc*2;			
			if (/*numSites*/numNuc > 800)
				{
				if (recombinationRate <= 0.00001)
					{
					maxSegNode = /*numSites*/numNuc;
					}
				else
					{
					maxSegNode = /*numSites*/numNuc*3;
					if (recombinationRate > 0.0001)
						fprintf(fpmpi, "\nrecombination rate is very high");	
					}
				}
			}
		else if (/*numSites*/numNuc > 1000 && /*numSites*/numNuc <= 10000)  /* NEW */
			{			

			if (recombinationRate <= 0.000001)
				{
				maxSegNode = /*numSites*/numNuc*1.5; /* NEW */
				numNodes = 5000;
				}
			if (recombinationRate <= 0.00001 && recombinationRate > 0.000001)
				{
				maxSegNode = /*numSites*/numNuc*2.3; /* NEW */
				numNodes = 5000;
				}
			if (recombinationRate <= 0.0001 && recombinationRate > 0.00001)
				maxSegNode = /*numSites*/numNuc*2.5;
			if (recombinationRate > 0.0001)
				{
				maxSegNode = /*numSites*/numNuc*3;	
				fprintf(fpmpi, "\nrecombination rate is very high");
				}
			}
		else if (/*numSites*/numNuc > 10000 && /*numSites*/numNuc <= 100000) /* large seqs.1 NEW */
			{
			
			if (recombinationRate <= 0.000001)
				{
				maxSegNode = /*numSites*/numNuc*0.5; /* NEW */
				numNodes = 1000;
				}
			if (recombinationRate <= 0.00001 && recombinationRate > 0.000001)
				{
				maxSegNode = /*numSites*/numNuc*0.5; /* NEW */
				numNodes = 6500;
				}
			if (recombinationRate > 0.0001)
				{
				maxSegNode = /*numSites*/numNuc*1.0;	
				fprintf(fpmpi, "\nrecombination rate is very high");
				}
			}
		else if (/*numSites*/numNuc > 100000 && /*numSites*/numNuc <= 1000000) /* large seqs.2 NEW */
			{
			
			if (recombinationRate <= 0.000001)
				{
				maxSegNode = /*numSites*/numNuc*0.005; /* NEW */
				numNodes = 350;
				}
			if (recombinationRate <= 0.00001 && recombinationRate > 0.000001)
				{
				maxSegNode = /*numSites*/numNuc*0.5; /* NEW */
				numNodes = 5000;
				}
			if (recombinationRate > 0.0001)
				{
				maxSegNode = /*numSites*/numNuc*1.0;	
				fprintf(fpmpi, "\nrecombination rate is very high");
				}
			}
		else
			{
			if (recombinationRate <= 0.00001)
				maxSegNode = /*numSites*/numNuc;
			if (recombinationRate <= 0.0001 && recombinationRate > 0.00001)
				maxSegNode = /*numSites*/numNuc*2;
			if (recombinationRate > 0.0001)
				{
				maxSegNode = /*numSites*/numNuc*3;	
				fprintf(fpmpi, "\nrecombination rate is very high");
				}
			}
		if (N > 1000)
			{
			maxSegNode = maxSegNode*1.5;
			if (N > 10000)
				maxSegNode = maxSegNode*1.5;
			if (N > 50000)
				maxSegNode = maxSegNode*1.1;
			}
		}

	if (doMigration == YES) /* migration implies more rec events because lineages can be isolated for long times */
		{
		numNodes = numNodes * 1.5;
		}


	if (numNodes <= numSequences) /* if there are more sequences than nodes in the initial step */
		numNodes = numSequences + 5000;

	numTotalSegments = (numNodes*maxSegNode)+/*numSites*/numNuc;
	/*fprintf(stderr, "\n numTotalSegments = %d \n", numTotalSegments);*/
	
	if (doMigration == YES) /* memory for migrations */
		{
		cumInitPopul = (int *)calloc((numPopulations+numConvergDemes+1),(long) sizeof(int));
		if (!cumInitPopul)
			{
			fprintf (fpmpi, "Could not allocate cumInitPopul (%lu bytes)\n", (numPopulations+numConvergDemes+1)  * (long) sizeof(int));
			exit (1);
			}
		numParcialActiveGametes = (int *) calloc((numPopulations+numConvergDemes+1),(long) sizeof(int));
		if (!numParcialActiveGametes)
			{
			fprintf (fpmpi, "Could not allocate numParcialActiveGametes (%lu bytes)\n", (numPopulations+numConvergDemes+1)  * (long) sizeof(int));
			exit (1);
			}
			
		GiPartial = (int *)calloc((numPopulations+numConvergDemes+1),(long) sizeof(int));
		if (!GiPartial)
			{
			fprintf (fpmpi, "Could not allocate GiPartial (%lu bytes)\n", (numPopulations+numConvergDemes+1)  * (long) sizeof(int));
			exit (1);
			}
		
		cumPopulTase =  (double*) calloc ((numPopulations+numConvergDemes+1), sizeof (double));  
		if (cumPopulTase == NULL)
			{
			fprintf (fpmpi, "Could not allocate cumPopulTase (%lu bytes)", (numPopulations+numConvergDemes+1)  * (long) sizeof(int));
			exit(1);
			}
		rateREpartial =  (double*) calloc ((numPopulations+numConvergDemes+1), sizeof (double));  
		if (rateREpartial == NULL)
			{
			fprintf (fpmpi, "Could not allocate rateREpartial (%lu bytes)", (numPopulations+numConvergDemes+1) * (long) sizeof (double));
			exit(1);
			}
		rateCApartial =  (double*) calloc ((numPopulations+numConvergDemes+1), sizeof (double));  
		if (rateCApartial == NULL)
			{
			fprintf (fpmpi, "Could not allocate rateCApartial (%lu bytes)", (numPopulations+numConvergDemes+1) * (long) sizeof (double));
			exit(1);
			}
		rateMIGpartial =  (double*) calloc ((numPopulations+numConvergDemes+1), sizeof (double));  
		if (rateMIGpartial == NULL)
			{
			fprintf (fpmpi, "Could not allocate rateMIGpartial (%lu bytes)", (numPopulations+numConvergDemes+1) * (long) sizeof (double));
			exit(1);
			}
		ratePartial =  (double*) calloc ((numPopulations+numConvergDemes+1), sizeof (double));  
		if (ratePartial == NULL)
			{
			fprintf (fpmpi, "Could not allocate ratePartial (%lu bytes)", (numPopulations+numConvergDemes+1) * (long) sizeof (double));
			exit(1);
			}
	
		for (i = 0; i < numPopulations+numConvergDemes+1; i++)
			{
			rateREpartial[i] = 0.00;
			rateCApartial[i] = 0.00;
			rateMIGpartial[i] = 0.00;
			ratePartial[i] = 0.00;
			cumPopulTase[i] = 0.00;
			cumInitPopul[i] = 0.00;
			numParcialActiveGametes[i] = 0.00;
			GiPartial[i] = 0.00;
			}
		}
		
	/* allocate space for tree */


	segments = (TreeSegment *) calloc (numTotalSegments, sizeof(TreeSegment)); /* segments */
	if (!segments)
		{
		fprintf (fpmpi, "Could not allocate segments (%lu bytes)\n", numTotalSegments  * (long) sizeof(TreeSegment));
		fprintf (fpmpi, "\n Too much memory is required, try again by less recombinationRate or less numSites or less numSequences.\n");			
		exit (1);
		}
		
	nodes = (TreeNode *) calloc (numNodes, sizeof(TreeNode)); /* nodes */
	if (!nodes)
		{
		fprintf (fpmpi, "Could not allocate nodes (%lu bytes)\n", numNodes  * (long) sizeof(TreeNode));
		exit (1);
		}
	
	
	
	for (mmm=0; mmm<numNodes; mmm++)
		{
		/*printf("\n   Allocating node %d", mmm+1);*/
		nodes[mmm].SitesNonAncHere = (int *) calloc (numNuc+1, sizeof(int));
		if (!nodes[mmm].SitesNonAncHere)
			{
			fprintf (fpmpi, "Could not allocate nodes[mmm].SitesNonAncHere (%lu)\n", (numNuc+1) * (long) sizeof(int));
			exit (1);
			}
		/*printf("\n      allocated %d units for node %d", (numNuc+1), mmm+1);*/
		}
	/*printf("\nFinished allocating nodes...\n");*/



	activeGametes = (int *) calloc (numNodes,(long) sizeof(int)); /* active nodes */
	if (!activeGametes)
		{
		fprintf (fpmpi, "Could not allocate activeGametes (%lu bytes)\n", numNodes *(long) sizeof(int));
		exit (1);
		}
		
	treeRootInit = (TreeNode **) calloc(1, sizeof(TreeNode *)); /* nodes pointers */
	if (!treeRootInit)
		{
		fprintf (fpmpi, "Could not allocate treeRootInit (%lu bytes)\n", 1  * (long) sizeof(TreeNode));
		exit (1);
		}
	stud = (int *)calloc(numSites,(long) sizeof(int));
	if (!stud)
		{			
		fprintf (fpmpi, "Could not allocate stud (%lu bytes)\n", numSites * (long) sizeof(int));
		exit (1);
		}
	for (w = 1; w < numSites; w++) 
		stud[w-1] = w*3+1; 	/* The breakpoints between codons. "stud" is an array with the possible breakpoints beetween codons*/



	
	/* allocate space for S_MRCA vector = [n1, n2, ..., ni,...nk] 
   	where n is the number of sequences and k is the number of sites
   	When an element of MRCS reaches 1, means that the site i found its MRCA */
	/* The MRCA vector has got the number of coalescences n, when 1 coalescence happen n=n-1. In the last coalescent n<=1, we are in the MRCA */
	S_MRCA = (int *) calloc((numNuc+1),(long) sizeof(int));
	if (!S_MRCA)
		{
		fprintf (fpmpi, "Could not allocate S_MRCA (%lu bytes)\n", (numNuc+1) *(long) sizeof(int));
		exit (-1);
		}
	OnlyAncS_MRCA = (int *) calloc((numNuc+1),(long) sizeof(int));
	if (!OnlyAncS_MRCA)
		{
		fprintf (fpmpi, "Could not allocate OnlyAncS_MRCA (%lu bytes)\n", (numNuc+1) *(long) sizeof(int));
		exit (-1);
		}


	/* allocate space for recombination breakpoints */
	breakpoint = (int *) calloc(memoryBreakp,(long) sizeof(int));
	if (!breakpoint)
		{
		fprintf (fpmpi, "Could not allocate breakpoint (%lu bytes)\n", memoryBreakp*(long) sizeof(int));
		exit (-1);
		}

	/* set everything to null */
	for (i=0; i<numNodes; i++)
		{
		nodes[i].left = NULL;
		nodes[i].right = NULL;
		nodes[i].anc1 = NULL;
		nodes[i].anc2 = NULL;
		nodes[i].outgroup = NULL;
		nodes[i].sib = NULL;
		nodes[i].time = 0;
		nodes[i].length = 0;
		nodes[i].index = 0;
		nodes[i].indexOldMigPop = 0;
		nodes[i].indexCurrentMigPop = 0;
		nodes[i].numSegNode = 0;
		nodes[i].label = 0;
		nodes[i].isOutgroup = NO;
		nodes[i].class = 0;
		nodes[i].breakp = NO;
		nodes[i].passNumber = 0;
		nodes[i].breakCodon = NO;
		nodes[i].breakCodon = 0;
		nodes[i].NetLabelPrint = 0;
		nodes[i].MRCAfrom = -1;
		nodes[i].MRCAto = -1;
		nodes[i].GMRCA_ancestral = NO;
		for (mmm=0; mmm<=numNuc; mmm++)
			nodes[i].SitesNonAncHere[mmm] = 0;
		}


	/* set times in generation for dated tips; the latest sample will be generation 0 */		
	if (doDatedTips == YES)
		{
		/* NEW */
		numActiveGametes = 0;
		actNumSegments = 0;

		if (doMigration == YES)
			{
			if (noisy > 1)
				fprintf (fpmpi,"\n Initial relation nodes-demes:");

			for (i = 1; i <= numPopulations; i++)
				{
				cumInitPopul[i] = cumInitPopul[i-1] + numNodesInitPopul[i];
				/*numParcialActiveGametes[i] = numNodesInitPopul[i];*/ /* do not active here */
				numParcialActiveGametes[i] = 0;
				/*fprintf (fpmpi,"\n cumInitPopul[i]: %d", cumInitPopul[i]);*/
				/*fprintf (fpmpi,"\n numParcialActiveGametes[i]: %d", numParcialActiveGametes[i]);*/
				}
			}

		for (j = 0; j < numSequences; j++)		
			{
			p = nodes + j;
			/*fprintf (stderr, "\n >> nodes + %d", j);*/
			/*activeGametes[numActiveGametes] = j;*/ /* do not active here */
			p->index = j;
			p->label = j;
			p->NetLabelPrint = numNetLabelPrint;
			/*fprintf (fpmpi,"\n numNetLabelPrint = %d", numNetLabelPrint);*/
			numNetLabelPrint++;		
			labelNodes = j;
			p->class = 1;
			p->GMRCA_ancestral = NO;
			p->numSegNode = 1;				/* initial, each node contain 1 segment */
			/*fprintf (fpmpi,"\n TIP, The node %d with class %d", p->index, p->class);*/
			
			if (doMigration == YES) /* initial nodes - populations */
				{
				if ((j+1) <= cumInitPopul[1])
					{
					p->indexOldMigPop = 1;
					p->indexCurrentMigPop = 1;
					}
				else
					{
					for (k = 1; k <= numPopulations; k++)
						{
						if ((j+1) > cumInitPopul[k] && (j+1) <= cumInitPopul[k+1])
							{
							p->indexOldMigPop = k+1;
							p->indexCurrentMigPop = k+1;
							break;
							}
						}
					k = 0;
					}
				if (noisy > 1)
					fprintf (fpmpi,"\n > The node %d belongs to deme %d", p->index, p->indexOldMigPop);
				}
			
			for (mmm=1; mmm<=numNuc; mmm++)
				p->SitesNonAncHere[mmm] = 0;

			for (i = 0; i < p->numSegNode; i++) /* segments of the tip */
				{
				s = segments + post(i,j,distance);	
				s->before1 = NULL;
				s->before2 = NULL;
				s->after1 = NULL;
				s->after2 = NULL;			
				s->sIndexNode = j;
				s->sIndex = actSegIndex;
				s->sStart = 1;
				s->sEnd = numNuc;
				
				actSegIndex++;
				}
			/*numActiveGametes++;*/ /* do not active here */
			}
		labelNodes++;
		/*nextAvailable = numActiveGametes;*/ /* do not active here */
		actNumSegments = numSequences;
		/* end NEW */

		for (i=0; i<numTipDates; i++)
			{
			for (j=0; j<datedSample[i].size; j++)
				{
				p = nodes + datedSample[i].member[j]-1;
				p->time =  datedSample[i].time;

				/*fprintf (stderr, "\n > nodes + %d", datedSample[i].member[j]-1);*/
				/*for (k=0; k<numSites; k++)
					{
					p = nodes + pos(datedSample[i].member[j]-1,k,numSites);
					p->time =  datedSample[i].time; 
					}*/ 	/* to get generations back (assumes generation time and sampling times are in the same unit, e.g., years) */
				}
			}
		/*debug */
		if (noisy > 1)
			fprintf (stderr, "\nTime Tips:");
		for (j=0; j<numSequences; j++)
			{
			p = nodes + j;
			if (noisy > 1)
				fprintf (stderr, "\n => node %d, sampling time %f", j, p->time);
			/*for (k=0; k<3; k++)
				{
				p = nodes + pos(j,k,numSites);
				fprintf (stderr, "\n%d =>seq %d, site %d, sampling time %f", p, j+1, k+1, p->time);
				}*/
			
			/* check - Fix only converge demes respect to tip dates (with less time than the time of the convergence of deme) s*/
			if (doConvergDemes == YES) 
				{
				ss=0;
				for (ss=1; ss<=numConvergDemes; ss++)
					{
					if (convDemTimes_old[ss] <= p->time)
						{
						fprintf (fpmpi,"\n\n Warning. Tip dates AND Convergence of demes.\n The time of any tip date (%f) cannot be higher than the time of the convergence of demes (%f)", p->time, convDemTimes_old[ss]);
						PrintUsage();
						}
					
					}
				ss=0;
				}
			}				
		}

			
		


	/* set up MRCA vector */
	for (i = 0; i <= numNuc; i++)
		{
		S_MRCA[i] = numSequences; /* MRCA[i] = numSequences; example initial values: MRCA[0] = 10, MRCA[1] = 10, MRCA[2] = 10 .. */
		OnlyAncS_MRCA[i] = numSequences; /* MRCA[i] = numSequences; example initial values: MRCA[0] = 10, MRCA[1] = 10, MRCA[2] = 10 .. */
		}
		



	/* tip of the trees */ /* set up initial nodes */
	if (doDatedTips == YES) /** active oldest sample for dated tips **/
		{
		currentSample = numTipDates - 1;
		if (noisy > 2)
			fprintf (stderr, "\nActivating initial sample %d (with time = %6.4f):", currentSample,datedSample[currentSample].time);

		numActiveGametes = 0;

		for (j=0; j<datedSample[currentSample].size; j++)
			{
			if (noisy > 2)
				fprintf (stderr, " %d", datedSample[currentSample].member[j]-1);

			/*fprintf (stderr, "\n > Activated nodes + %d \n", datedSample[currentSample].member[j]-1);*/
			p = nodes + datedSample[currentSample].member[j]-1;
			p->index = datedSample[currentSample].member[j]-1;
			activeGametes[numActiveGametes] = datedSample[currentSample].member[j]-1; 
			numActiveGametes++;
			
			if (doMigration == YES)
				{			
				for (i = 1; i <= numPopulations; i++) 
					{
					if (p->indexOldMigPop == i)
						{
						numParcialActiveGametes[i]++;
						}
					}
				}
			}
		nextAvailable = numSequences;  /* to keep the indexes 1 to numSequences for tips */
		if (doMigration == YES && noisy > 2)
			{
			fprintf (fpmpi,"\nInitial demes (with time = %6.4f): \n", datedSample[currentSample].time);			
			for (i = 1; i <= numPopulations; i++) 
				fprintf (fpmpi," Deme %d with %d nodes\n", i, numParcialActiveGametes[i]);
			}
		/*labelNodes++;
		actNumSegments = numActiveGametes-1;*/
		}	
	else /** Here, not tip dates **/
		{
		numActiveGametes = 0;
		actNumSegments = 0;
		
		if (doMigration == YES)
			{
			for (i = 1; i <= numPopulations; i++)
				{
				cumInitPopul[i] = cumInitPopul[i-1] + numNodesInitPopul[i];
				numParcialActiveGametes[i] = numNodesInitPopul[i];
				/*fprintf (fpmpi,"\n numParcialActiveGametes[i]: %d \n", numParcialActiveGametes[i]);*/
				}
			if (noisy > 1)
				fprintf (fpmpi,"\n Initial relation nodes-demes:");
			}
	
		for (j = 0; j < numSequences; j++)		
			{
			p = nodes + j;
			activeGametes[numActiveGametes] = j;
			p->index = j;
			p->label = j;
		
			p->NetLabelPrint = numNetLabelPrint;
			/*fprintf (fpmpi,"\n numNetLabelPrint = %d", numNetLabelPrint);*/
			numNetLabelPrint++;		

			labelNodes = j;
			p->class = 1;
			p->GMRCA_ancestral = NO;
		
			p->numSegNode = 1;				/* initial, each node contain 1 segment */
			if (doMigration == YES) /* initial nodes - populations */
				{
				if ((j+1) <= cumInitPopul[1])
					{
					p->indexOldMigPop = 1;
					p->indexCurrentMigPop = 1;
					}
				else
					{
					for (k = 1; k <= numPopulations; k++)
						{
						if ((j+1) > cumInitPopul[k] && (j+1) <= cumInitPopul[k+1])
							{
							p->indexOldMigPop = k+1;
							p->indexCurrentMigPop = k+1;
							break;
							}
						}
					k = 0;
					}
				if (noisy > 1)
					fprintf (fpmpi,"\n > The node %d belongs to deme %d", p->index, p->indexOldMigPop);
				
				for (mmm=1; mmm<=numNuc; mmm++)
					p->SitesNonAncHere[mmm] = 0;
				}

			/*fprintf (fpmpi,"\n TIP, The node %d with class %d", p->index, p->class);*/

			for (i = 0; i < p->numSegNode; i++) /* segments of the tip */
				{
				s = segments + post(i,j,distance);	
				s->before1 = NULL;
				s->before2 = NULL;
				s->after1 = NULL;
				s->after2 = NULL;			
				s->sIndexNode = j;
				s->sIndex = actSegIndex;
				s->sStart = 1;
				s->sEnd = numNuc;
				
				actSegIndex++;
				}
			numActiveGametes++;
			}
	
		labelNodes++;
		nextAvailable = numActiveGametes;
		actNumSegments = numSequences;
		}
	

	if (doDatedTips == YES)
		{
		if (noisy > 2)
			{
			fprintf (stderr, "\nActive nodes (%d):", numActiveGametes); 
			for (i=0; i<numActiveGametes; i++)
				fprintf (stderr," %d",activeGametes[i]);
			fprintf (stderr,"   Next node available = %d ", nextAvailable);
			}
		if (noisy > 2)
			{			
			fprintf (fpmpi, "\n\n Initial MRCA in the nodes:\n");
			for (i=0; i<numActiveGametes; i++)
				{
				for (j=1; j<=numNuc; j++)
					{
					p = nodes + activeGametes[i];
					/*p = nodes + pos(activeGametes[i],j,numNuc);*/
					if (j == 1)
						fprintf (fpmpi, "%4d -- (MRCA:)", p->index);
					fprintf (fpmpi, "%d", S_MRCA[j]);
					}
				fprintf (fpmpi, "\n");
				}
			}	
		}
	else
		{
		/* print out ancestral (active) status for each site */
		if (noisy > 2)
			{			
			fprintf (fpmpi, "\n\n Initial MRCA in the nodes:\n");
			for (i=0; i<numActiveGametes; i++)
				{
				for (j=1; j<=numNuc; j++)
					{
					p = nodes + activeGametes[i];
					/*p = nodes + pos(activeGametes[i],j,numNuc);*/
					if (j == 1)
						fprintf (fpmpi, "%4d -- (MRCA:)", p->index);
					fprintf (fpmpi, "%d", S_MRCA[j]);
					}
				fprintf (fpmpi, "\n");
				}
			}	
		}




	/*** make coalescence tree ***/	
	if (doMigration == YES) /* Coalescence with migration */
		{		
		eventNum = 0;		
		currentTime = 0.0;
		period = 1;
		
		if (doConvergDemes == YES) /* these arrays are for the convergence demes events */
			{
			deme_a = (int *) calloc(numConvergDemes+1,(long) sizeof(int));
			if (!deme_a)
				{
				fprintf (stderr, "PARAMETER ERROR: Could not allocate deme_a of convergencies demes events (%lu bytes)\n", (numConvergDemes+1) *(long) sizeof(int));
				exit (1);
				}
			deme_b = (int *) calloc(numConvergDemes+1,(long) sizeof(int));
			if (!deme_b)
				{
				fprintf (stderr, "PARAMETER ERROR: Could not allocate deme_b of convergencies demes events (%lu bytes)\n", (numConvergDemes+1) *(long) sizeof(int));
				exit (1);
				}
			convDemTimes = (double*) calloc ((numConvergDemes+1), sizeof (double)); 
			if (convDemTimes == NULL)
				{
				fprintf (stderr, "PARAMETER ERROR: Could not allocate convDemTimes of convergencies demes events (%lu bytes)\n", (numConvergDemes+1) *(long) sizeof(double));
				exit (1);
				}
			currentConvDem = (int *) calloc(numConvergDemes+1,(long) sizeof(int));
			if (!currentConvDem)
				{
				fprintf (stderr, "PARAMETER ERROR: Could not allocate currentConvDem of convergencies demes events (%lu bytes)\n", (numConvergDemes+1) *(long) sizeof(int));
				exit (1);
				}
			CurrentDemesState = (int *)calloc((2*numPopulations),(long) sizeof(int));
			if (!CurrentDemesState)
				{
				fprintf (fpmpi, "Could not allocate CurrentDemesState (%lu bytes)\n", (2*numPopulations)  * (long) sizeof(int));
				exit (1);
				}
			for (k = 1; k <= 2*numPopulations-1; k++)
				CurrentDemesState[k] = 0;
			CurrentDemesState[0] = 0;
			for (k = 1; k <= numPopulations; k++)
				CurrentDemesState[k] = k;
			for (k=1; k<=numConvergDemes; k++)
				{
				deme_a[k] = deme_a_old[k];
				deme_b[k] = deme_b_old[k];
				convDemTimes[k] = convDemTimes_old[k];
				/*fprintf(stderr, "\n\n deme_a[%d] = %d \n", k, deme_a[k]);
				fprintf(stderr, "\n\n deme_b[%d] = %d \n", k, deme_b[k]);
				fprintf(stderr, "\n\n convDemTimes[%d] = %lf \n", k, convDemTimes[k]);*/
				}
			k = 0;
			}

		while (numActiveGametes > 1)
			{
			/* print out MRCA vector */
			/*if (noisy == 8)
				{	
				fprintf(stderr,"S_MRCA [] = ");
				for (j=1; j<=numNuc; j++)
					fprintf (stderr,"%d ", S_MRCA[j]);
				fprintf(stderr,"\n");
				}*/

			/*for (i = 1; i <= 2*numPopulations-1; i++)
				fprintf (fpmpi,"\n CurrentDemesState[%d] = %d\n", i, CurrentDemesState[i]);*/  /* to see the active demes */

			Gi = 0;
			/* allocate memory for each node gi */
			gi = (int *) calloc(numActiveGametes,(long) sizeof(int));
			if (!gi)
				{
				fprintf (fpmpi, "Could not allocate gi (%lu bytes)\n", numActiveGametes *(long) sizeof(int));	
				exit (-1);
				}
			/* calculate gi for each node and total Gi */	
			/* Gi is the total number of ancestral sites, gi is a vector with ancestral and not found MRCA sites */
		
			for (i = 0; i < numActiveGametes; i++)
				{
				p = nodes + activeGametes[i];
				sizeNode = p->numSegNode;
				gi[i] = CalcIndividualGi (i, nodes, activeGametes, numNuc, S_MRCA, sizeNode);
				Gi += gi[i];
				if (noisy == 4)
					fprintf (fpmpi,"\n%d \n", gi[i]);
				}	
			if (noisy == 4)
				fprintf (fpmpi," Gi = %lu ", Gi);
			for (k = 1; k <= numPopulations+currentBigDeme; k++)
				GiPartial[k] = 0;
			
			if (noisy == 4)
				fprintf (fpmpi, "\n");
			for (i = 0; i < numActiveGametes;i++)
				{
				p = nodes + activeGametes[i];
				for (j = 1; j <= numPopulations+currentBigDeme; j++)
					{
					if (p->indexCurrentMigPop == j)
						GiPartial[j] = GiPartial[j] + gi[i];
					}
				if (noisy > 2)
					fprintf (fpmpi, " \nNode %d of deme %d", p->index, p->indexCurrentMigPop);
				}
			if (noisy == 3)
				fprintf (fpmpi, "\n");			

			/* get rates for events */  
			rateRE = rateCA = rateMIG = rate = 0;
			for (i = 1; i <= numPopulations+currentBigDeme; i++)
				{
				rateREpartial[i] = GiPartial[i] * Nscaling * N * recombinationRate;
				rateRE = rateRE + rateREpartial[i];
				/*fprintf (fpmpi, " \nrateREpartial[%d] = %lf", i, rateREpartial[i]);*/

				rateCApartial[i] = numParcialActiveGametes[i] * (numParcialActiveGametes[i] - 1) / 2.0;
				rateCA = rateCA + rateCApartial[i];
				/*fprintf (fpmpi, " \nrateCApartial[%d] = %lf", i, rateCApartial[i]);*/

				rateMIGpartial[i] = numParcialActiveGametes[i] * Nscaling * N * migrationRate;
				if (doConvergDemes == YES && currentDemesNumber <= 1)  /* No more migrations because there is only one deme */
					rateMIGpartial[i] = 0.00;
				rateMIG = rateMIG + rateMIGpartial[i];
				/*fprintf (fpmpi, " \nnumParcialActiveGametes[%d] = %d", i, numParcialActiveGametes[i]);
				fprintf (fpmpi, " \nrateMIGpartial [%d] = %lf", i, rateMIGpartial[i]);*/

				ratePartial[i] = rateREpartial[i] + rateCApartial[i] + rateMIGpartial[i];
				rate = rate + ratePartial[i];
				}
			rate = rateRE + rateCA + rateMIG;
			/* find out time for coalescence */
			if (doDemographics == YES)
				{
				periodGrowth[period] = -log (Nend[period] / (double) Nbegin[period]) / (double) (cumDuration[period] - cumDuration[period-1]);
				/*fprintf (fpmpi, "\n>>>>1 period growth  = %f (period = %d)", periodGrowth[period], period);*/
				if (isnan(periodGrowth[period]) == YES)
					{						
					fprintf (fpmpi, "\nERROR: period growth (%f) is NaN", periodGrowth[period]);
					fprintf (fpmpi, "\n      This might suggest that the growth rate is too negative");
					fprintf (fpmpi, "\n      and the coalescent time is therefore infinite.");
					fprintf (fpmpi, "\n      Try a smaller value");				
					exit (1);
					}

				if (Nend[period] == Nbegin[period])
					{
					timeCA = RandomExponential (rateCA, seed) * Nscaling * (double) Nbegin[period];
					}
				else
					{
					timeCA = log (1 + RandomExponential (rateCA, seed) * periodGrowth[period] * Nscaling * Nbegin[period] * 
			               exp (-periodGrowth[period] * (currentTime - cumDuration[period-1]))) / periodGrowth[period];
					}

				/*	When growth rate is very negative, coalescent time may be infinite
					this results in log (-x) => timCA = NaN. If this not the last period
					just jump to the next. If this is the last period, we have to exit
					the program */
				if (isnan(timeCA) == YES)
					{
					if (period < numPeriods) 
						{
						currentTime = cumDuration[period];
						period++;
						periodGrowth[period] = -log (Nend[period] / (double) Nbegin[period]) / (double) (cumDuration[period] - cumDuration[period-1]);
						continue;
						}
					else
						{
						fprintf (fpmpi, "\nERROR: Coalescent time (%f) is infinite ", timeCA);
						fprintf (fpmpi, "\n      This might suggest that the growth rate is too negative");
						fprintf (fpmpi, "\n      and the coalescent time is therefore infinite.");
						fprintf (fpmpi, "\n      Try a smaller value");
						exit (1);
						}
					}
				}
			else
				{
				timeCA = RandomExponential (rateCA, seed) * Nscaling * N;
				if (doExponential == YES)
					{
					timeCA = log (exp(growthRate*currentTime) + growthRate * timeCA) / growthRate - currentTime;
		
					/*	When growth rate is very negative, coalescent time may be infinite
						this results in log (-x) => timeCA = NaN. We have to exit
						the program */
					if (isnan(timeCA) == YES)
						{				
						fprintf (fpmpi, "\nERROR: Coalescent time (%f) is infinite ", timeCA);
						fprintf (fpmpi, "\n      This might suggest that the growth rate is too negative");
						fprintf (fpmpi, "\n      and the coalescent time is therefore infinite.");
						fprintf (fpmpi, "\n      Try a smaller value");
						exit (1);
						}
					}
				}

			/* find out time for recombination */
			timeRE = RandomExponential (rateRE, seed) * Nscaling * N; /* It doesn't depend of demographics because the node doesn«t depend of another node to make this event */
			/* find out time for migration */
			timeMIG = RandomExponential (rateMIG, seed) * Nscaling * N; /* It doesn't depend of demographics because the node doesn«t depend of another node to make this event */
			/* find out time for convergence of demes */	/* It doesn't depend of demographics because the node doesn«t depend of another node to make this event */		
			w = j = d = totCurrentConv = 0;
			if (doConvergDemes == YES)
				{
				for (k=1; k<=numConvergDemes; k++)
					currentConvDem[k] = 0;
			
				for (k=1; k<=numConvergDemes; k++) /* looking for the event nearer (with less time and future) */
					{ 
					if (convDemTimes[k] >= currentTime && convDemTimes[k] >= 0)
						{
						if (d == 0)
							j = k;
						d++;
						if (convDemTimes[k] >= currentTime && convDemTimes[k] < j)
							j = k;
						}
					}
				if (d > 0) /* there are convergences */
					{
					//fprintf(fpmpi,"\n proximo de menor tiempo = %d, con tiempo %lf", j, convDemTimes[j]);
					nextConvNumber = j;
					timeCONV = convDemTimes[nextConvNumber] - currentTime;
					/*fprintf(fpmpi,"\n timeCONV = %lf", timeCONV);*/
					}
				else
					timeCONV = 0; /* there is no more convergences */
				}
			/*fprintf(fpmpi,"\n timeCA = %lf, timeRE = %lf, timeMIG = %lf \n", timeCA,timeRE,timeMIG);*/





			if (doDatedTips == YES)
				{
				if (doConvergDemes == YES && timeCONV > 0)
					{
					if (timeCA < timeRE && timeCA < timeMIG && timeCA < timeCONV) /* coalescence event */
						{
						isCoalescence = YES;
						isRecombination = NO;
						isMigration = NO;
						doConvNext = NO;
						eventTime = timeCA;
						}
					else if (timeRE < timeCA && timeRE < timeMIG && timeRE < timeCONV) /* recombination event */
						{
						isCoalescence = NO;
						isRecombination = YES;
						isMigration = NO;
						doConvNext = NO;
						eventTime = timeRE;
						}
					else if (timeMIG < timeCA && timeMIG < timeRE && timeMIG < timeCONV) /* migration event */
						{
						isCoalescence = NO;
						isRecombination = NO;
						isMigration = YES;
						doConvNext = NO;

						eventTime = timeMIG;
						}
					else if (timeCONV < timeCA && timeCONV < timeRE && timeCONV < timeMIG) /* convergence between demes event */
						{
						isCoalescence = NO;
						isRecombination = NO;
						isMigration = NO;
						doConvNext = YES;
						eventTime = timeCONV;
						}
					else
						{
						fprintf(fpmpi, "\n\n Warning choosing the type of event");
						fprintf(fpmpi, "\n Check:: \n 1. Do you have a migration rate = 0 starting by more than one demes and without convergence events for the demes? \n 2. Do you have a migration rate = 0 without convergence events for ALL demes IN TIME BEFORE? \n When two demes converg the result deme has a new number that (if there are more demes) it should be converged or migrated. So, e.g. -%%2 1 2 200 3 4 500, in 200 1 more 2 make 5, then in 500 3 and 4 make 6, so finally you have deme 5 and deme 6 as independent demes, it is a problem with migration rate = 0 \n So, if you are working with migration rate = 0 it is recomending convergence all demes, including the new demes (e.g. - %%3 1 2 200 3 4 500 5 6 2000) ");
						fprintf(fpmpi, "\n 3. If migration rate is different than 0: Check if you are using tip dates and convergence of demes, the sample could belongs to a time older than the existence of its user-given deme. To solve this problem try to erase this tip date (this sample time = 0.00) or to modify the time of the tip node to a value younger than the convergence of demes \n");
						exit (-1);
						}

					}
				if (doConvergDemes == NO || timeCONV == 0)
					{
					if (timeCA < timeRE && timeCA < timeMIG) /* coalescence event */
						{
						isCoalescence = YES;
						isRecombination = NO;
						isMigration = NO;
						eventTime = timeCA;
						}
					else if (timeRE < timeCA && timeRE < timeMIG) /* recombination event */
						{
						isCoalescence = NO;
						isRecombination = YES;
						isMigration = NO;
						eventTime = timeRE;
						}
					else if (timeMIG < timeCA && timeMIG < timeRE) /* migration event */
						{
						isCoalescence = NO;
						isRecombination = NO;
						isMigration = YES;
			
						eventTime = timeMIG;
						}
					else
						{
						fprintf(fpmpi, "\n\n Warning choosing the type of event");
						fprintf(fpmpi, "\n Check::: \n 1. Do you have a migration rate = 0 starting by more than one demes and without convergence events for the demes? \n 2. Do you have a migration rate = 0 without convergence events for ALL demes IN TIME BEFORE? \n When two demes converg the result deme has a new number that (if there are more demes) it should be converged or migrated. So, e.g. -%%2 1 2 200 3 4 500, in 200 1 more 2 make 5, then in 500 3 and 4 make 6, so finally you have deme 5 and deme 6 as independent demes, it is a problem with migration rate = 0 \n So, if you are working with migration rate = 0 it is recomending convergence all demes, including the new demes (e.g. - %%3 1 2 200 3 4 500 5 6 2000) ");
						fprintf(fpmpi, "\n 3. If migration rate is different than 0: Check if you are using tip dates and convergence of demes, the sample could belongs to a time older than the existence of its user-given deme. To solve this problem try to erase this tip date (this sample time = 0.00) or to modify the time of the tip node to a value younger than the convergence of demes \n");
						exit (-1);
						}

					}
				}
			
			if (doDatedTips == YES)
				{
				/* if doing dated tips, check whether we need to activate a new sample, update sampling period and start again */
				if ((currentTime + eventTime) > datedSample[currentSample-1].time && currentSample > 0)
					{
					currentSample--;
					/* activate nodes from this sample */
					if (noisy > 2)
						fprintf (stderr, "\nCumulative time = %6.4f  > sample %d time = %6.4f. Activating tips:", currentTime + eventTime, currentSample, datedSample[currentSample].time);
					
					for (i=0; i<datedSample[currentSample].size; i++)
						{
						if (noisy > 2)
							fprintf (stderr, " %d", datedSample[currentSample].member[i]-1);

						p = nodes + datedSample[currentSample].member[i]-1;
						p->index = datedSample[currentSample].member[i]-1;
						activeGametes[numActiveGametes] = datedSample[currentSample].member[i]-1;  

						numActiveGametes++;
											
						for (ss = 1; ss <= numPopulations /*+ numCONV*/; ss++) 
							{
							if (p->indexOldMigPop == ss)
								{
								numParcialActiveGametes[ss]++;
								/*fprintf (stderr, "\n AQUI numParcialActiveGametes[%d] = %d \n", ss, numParcialActiveGametes[ss]);*/ 
								}
							}
						}
					currentTime = datedSample[currentSample].time;

					if (noisy > 2)
						{
						fprintf (stderr, "\nActive nodes (%d):", numActiveGametes); 
						for (i=0; i<numActiveGametes; i++)
							fprintf (stderr," %d",activeGametes[i]);
						fprintf (stderr,"\nNext node available = %d", nextAvailable);
						fprintf (stderr, "\nSample %d activated and going back to currentTime = %6.4f", currentSample, currentTime);
						}
					continue; /* start again*/
					}
				}

			

			if (doConvergDemes == YES && timeCONV > 0)
				{
				/* event is a coalescence, a recombination or a migration? */
				if (timeCA < timeRE && timeCA < timeMIG && timeCA < timeCONV) /* coalescence event */
					{
					isCoalescence = YES;
					isRecombination = NO;
					isMigration = NO;
					doConvNext = NO; /* new */
					eventTime = timeCA;

					/*	if this period is not the last one and if the event time is outside the current interval,
						update period and start again */
					if (doDemographics == YES && period < numPeriods && (currentTime + eventTime) > cumDuration[period])
						{
						currentTime = cumDuration[period];
						period++;
						periodGrowth[period] = -log (Nend[period] / (double) Nbegin[period]) / (double) (cumDuration[period] - cumDuration[period-1]);
						continue;
						}	
			
					numCA++;
					}
				else if (timeRE < timeCA && timeRE < timeMIG && timeRE < timeCONV) /* recombination event */
					{
					isCoalescence = NO;
					isRecombination = YES;
					isMigration = NO;
					doConvNext = NO; /* new */
					eventTime = timeRE;	

					/*	if this period is not the last one and if the event time is outside the current interval,
						update period and start again */
					if (doDemographics == YES && period < numPeriods && (currentTime + eventTime) > cumDuration[period])
						{
						currentTime = cumDuration[period];
						period++;
						periodGrowth[period] = -log (Nend[period] / (double) Nbegin[period]) / (double) (cumDuration[period] - cumDuration[period-1]);
						continue;
						}	
					numRE++;
				
					/* reallocate space for recombination breakpoints */
					if (numRE >= (memoryBreakp-1))
						{
						memoryBreakp += 50;
					 
						breakpoint = (int *) realloc(breakpoint, memoryBreakp *(long) sizeof(int));
						if (!breakpoint)
							{	
							fprintf (fpmpi, "Could not reallocate breakpoint \n");
							exit (-1);
							}
						if (noisy == 4)
							fprintf (fpmpi, "\n...Doing reallocation of breakponts (1)\n");
						}
				
					if (numRE > (numNuc+1) && many == 0)
						{
						if (noisy > 1)
							fprintf (fpmpi, "\n\n Many recombinations %d (more recombinations that sites)!\n", numRE);
						
						many++;
						}
					}
				else if (timeMIG < timeCA && timeMIG < timeRE && timeMIG < timeCONV) /* migration event */
					{
					isCoalescence = NO;
					isRecombination = NO;
					isMigration = YES;
					doConvNext = NO;

					eventTime = timeMIG;	

					/*	if this period is not the last one and if the event time is outside the current interval,
						update period and start again */
					if (doDemographics == YES && period < numPeriods && (currentTime + eventTime) > cumDuration[period])
						{
						currentTime = cumDuration[period];
						period++;
						periodGrowth[period] = -log (Nend[period] / (double) Nbegin[period]) / (double) (cumDuration[period] - cumDuration[period-1]);
						continue;
						}	
					numMIG++;
					}
				else if (timeCONV < timeCA && timeCONV < timeRE && timeCONV < timeMIG) /* convergence between demes event */
					{
					isCoalescence = NO;
					isRecombination = NO;
					isMigration = NO;
					doConvNext = YES;
					eventTime = timeCONV;

					convDemTimes[nextConvNumber] = -1; /* erasing the event */	
					/*	if this period is not the last one and if the event time is outside the current interval,
						update period and start again */
					if (doDemographics == YES && period < numPeriods && (currentTime + eventTime) > cumDuration[period])
						{
						currentTime = cumDuration[period];
						period++;
						periodGrowth[period] = -log (Nend[period] / (double) Nbegin[period]) / (double) (cumDuration[period] - cumDuration[period-1]);
						continue;
						}	
					numCONV++;
					}
				else
					{
					fprintf(fpmpi, "\n\n Warning choosing the type of event");
					fprintf(fpmpi, "\n Check: \n 1. Do you have a migration rate = 0 starting by more than one demes and without convergence events for the demes? \n 2. Do you have a migration rate = 0 without convergence events for ALL demes IN TIME BEFORE? \n When two demes converg the result deme has a new number that (if there are more demes) it should be converged or migrated. So, e.g. -%%2 1 2 200 3 4 500, in 200 1 more 2 make 5, then in 500 3 and 4 make 6, so finally you have deme 5 and deme 6 as independent demes, it is a problem with migration rate = 0 \n So, if you are working with migration rate = 0 it is recomending convergence all demes, including the new demes (e.g. - %%3 1 2 200 3 4 500 5 6 2000) ");
					fprintf(fpmpi, "\n 3. If migration rate is different than 0: Check if you are using tip dates and convergence of demes, the sample could belongs to a time older than the existence of its user-given deme. To solve this problem try to erase this tip date (this sample time = 0.00) or to modify the time of the tip node to a value younger than the convergence of demes \n");
					exit (-1);
					}
				}
			if (doConvergDemes == NO || timeCONV == 0)
				{
				doConvNext = NO;
				/* event is a coalescence, a recombination or a migration? */
				if (timeCA < timeRE && timeCA < timeMIG) /* coalescence event */
					{
					isCoalescence = YES;
					isRecombination = NO;
					isMigration = NO;
					eventTime = timeCA;
					
					/*	if this period is not the last one and if the event time is outside the current interval,
						update period and start again */
					if (doDemographics == YES && period < numPeriods && (currentTime + eventTime) > cumDuration[period])
						{
						currentTime = cumDuration[period];
						period++;
						periodGrowth[period] = -log (Nend[period] / (double) Nbegin[period]) / (double) (cumDuration[period] - cumDuration[period-1]);
						continue;
						}	
			
					numCA++;
					}
				else if (timeRE < timeCA && timeRE < timeMIG) /* recombination event */
					{
					isCoalescence = NO;
					isRecombination = YES;
					isMigration = NO;
					eventTime = timeRE;
										
					/*	if this period is not the last one and if the event time is outside the current interval,
						update period and start again */
					if (doDemographics == YES && period < numPeriods && (currentTime + eventTime) > cumDuration[period])
						{
						currentTime = cumDuration[period];
						period++;
						periodGrowth[period] = -log (Nend[period] / (double) Nbegin[period]) / (double) (cumDuration[period] - cumDuration[period-1]);
						continue;
						}	
					numRE++;
				
					/* reallocate space for recombination breakpoints */
					if (numRE >= (memoryBreakp-1))
						{
						memoryBreakp += 50;
					 
						breakpoint = (int *) realloc(breakpoint, memoryBreakp *(long) sizeof(int));
						if (!breakpoint)
							{	
							fprintf (fpmpi, "Could not reallocate breakpoint \n");
							exit (-1);
							}
						if (noisy == 4)
							fprintf (fpmpi, "\n...Doing reallocation of breakponts (1)\n");
						}
				
					if (numRE > (numNuc+1) && many == 0)
						{
						if (noisy > 1)
							fprintf (fpmpi, "\n\n Many recombinations %d (more recombinations that sites)!\n", numRE);
						
						many++;
						}
					}
				else if (timeMIG < timeCA && timeMIG < timeRE) /* migration event */
					{
					isCoalescence = NO;
					isRecombination = NO;
					isMigration = YES;
			
					eventTime = timeMIG;	

					/*	if this period is not the last one and if the event time is outside the current interval,
						update period and start again */
					if (doDemographics == YES && period < numPeriods && (currentTime + eventTime) > cumDuration[period])
						{
						currentTime = cumDuration[period];
						period++;
						periodGrowth[period] = -log (Nend[period] / (double) Nbegin[period]) / (double) (cumDuration[period] - cumDuration[period-1]);
						continue;
						}	
					numMIG++;
					}
				else
					{
					fprintf(fpmpi, "\n\n Warning choosing the type of event");
					fprintf(fpmpi, "\n Check: \n 1. Do you have a migration rate = 0 starting by more than one demes and without convergence events for the demes? \n 2. Do you have a migration rate = 0 without convergence events for ALL demes IN TIME BEFORE? \n When two demes converg the result deme has a new number that (if there are more demes) it should be converged or migrated. So, e.g. -%%2 1 2 200 3 4 500, in 200 1 more 2 make 5, then in 500 3 and 4 make 6, so finally you have deme 5 and deme 6 as independent demes, it is a problem with migration rate = 0 \n So, if you are working with migration rate = 0 it is recomending convergence all demes, including the new demes (e.g. - %%3 1 2 200 3 4 500 5 6 2000) ");
					fprintf(fpmpi, "\n 3. If migration rate is different than 0: Check if you are using tip dates and convergence of demes, the sample could belongs to a time older than the existence of its user-given deme. To solve this problem try to erase this tip date (this sample time = 0.00) or to modify the time of the tip node to a value younger than the convergence of demes \n");
					exit (-1);
					}
				}




			/** it chooses the deme for this event **/
			if (isCoalescence == YES)
				{
				for (k = 1; k <= numPopulations+currentBigDeme; k++)
					cumPopulTase[k] = 0;
				cumPopulTase[0] = 0;
				w = 0;
						
				for (k = 1; k <= numPopulations+currentBigDeme; k++)
					{
					cumPopulTase[k] = cumPopulTase[k-1] + rateCApartial[k];
					w = w + numParcialActiveGametes[k];
					}
				/* (cumPopulTase[numPopulations+currentBigDeme] must to be similar to rateCA) */
				if (w != numActiveGametes)
					{
					fprintf(fpmpi, "\n 1 The sum of partial active gametes is different than total gametes number, w %d != numActiveGametes %d. In choose the deme", w, numActiveGametes);
					exit (-1);
					}
				for (k = 1; k <= numPopulations+currentBigDeme; k++)
					cumPopulTase[k] = cumPopulTase[k]/cumPopulTase[numPopulations+currentBigDeme];
					
				ran = RandomUniform(seed);
				whichDeme = bbinDemes(ran, cumPopulTase, numPopulations+currentBigDeme);
				w = 0;
				}
		
			if (isMigration == YES) 
				{
				for (k = 1; k <= numPopulations+currentBigDeme; k++)
					cumPopulTase[k] = 0;
				cumPopulTase[0] = 0;
				w = 0;
			
				for (k = 1; k <= numPopulations+currentBigDeme; k++)
					{
					cumPopulTase[k] = cumPopulTase[k-1] + rateMIGpartial[k];
					w = w + numParcialActiveGametes[k];
					}
				/* (cumPopulTase[numPopulations+currentBigDeme] must to be similar to rateMIG) */
				if (w != numActiveGametes)
					{
					fprintf(fpmpi, "\n 2The sum of partial active gametes is different than total gametes number, w %d != numActiveGametes %d. In choose the deme", w, numActiveGametes);
					exit (-1);
					}

				for (k = 1; k <= numPopulations+currentBigDeme; k++)
					cumPopulTase[k] = cumPopulTase[k]/cumPopulTase[numPopulations+currentBigDeme];
				/*fprintf(fpmpi, "\n cumPopulTase[numPopulations+currentBigDeme] = %lf", cumPopulTase[numPopulations+currentBigDeme]);*/
				ran = RandomUniform(seed);
				whichDeme = bbinDemes(ran, cumPopulTase, numPopulations+currentBigDeme);
				w = 0;
				}
			/* for event of convergent demes we do not gave to choose the deme */
			




			/** set time **/
			currentTime += eventTime; /* the time is accumulated */
			if (noisy > 3)
				fprintf(fpmpi, "\n\n");
			eventNum++;
			if (noisy > 1 && doConvNext == NO)
				fprintf (fpmpi, "\n\n*** Event %3d *** rate = %lf, currentTime = %lf\n", eventNum, rate, currentTime);
			if (noisy > 1 && doConvNext == YES)
				fprintf (fpmpi, "\n\n*** Event %3d *** currentTime = %lf\n", eventNum, currentTime);

				
			/*** if RECOMBINATION ***/
			if (isRecombination == YES) 
				{
				if (noisy == 4)
					fprintf(fpmpi, "\n* Recombination *");
					
				w = 0;
				for (k = 1; k <= numPopulations+currentBigDeme; k++)
					w = w + numParcialActiveGametes[k];
				if (w != numActiveGametes)
					{
					fprintf(fpmpi, "\n 3The sum of partial active gametes is different than total gametes number, w %d != numActiveGametes %d. In choose the deme", w, numActiveGametes);
					exit (-1);
					}	
				w = 0;
				
				/* Which node has the recombination? */
					/* assign probability to each node based on their gi's values */
				cum_gi = 0;
				probRecIndividual = Gi * RandomUniform(seed); /* it calculate the individual probability for the breakpoint */
				for (whichInd=0; whichInd<numActiveGametes; whichInd++)
					{
					cum_gi += gi[whichInd];				/* accumulate gi into cum_gi. whichInd is the choose node. */
					if (probRecIndividual < cum_gi)    /* ok whether the node is choosen */
						break;
					}
			
				if (whichInd >= numActiveGametes)
					{
					fprintf (fpmpi, "\n\nERROR: whichInd out of range1!: whichInd = %d\n", whichInd);
					exit (-1);
					}
			
				/* select a valid breakpoint among potential recombining locations */
				/* to be a potential recombining site, a site has to have ancestral material non-MRCA before and after it */
				legalBreakpoint = NO;
				while (legalBreakpoint == NO)
					{
					do /* oct2009 */
						{
						whichSite = (/*numSites*/numNuc) * RandomUniform(seed); /* it choose a site */
						if (whichSite == 1)
							whichSite = numNuc;
						} while (whichSite == 0);
					/* oct2009 */

					if (whichSite > numNuc)
						{
						fprintf (fpmpi, "\n\nERROR: whichSite out of range! : whichSite = %d\n", whichSite);
						exit (-1);
						}
					p = nodes + activeGametes[whichInd];
					sizeNode = p->numSegNode;
					if	(IsValidBreakSite (activeGametes, nodes, whichInd, whichSite, S_MRCA) == YES)
						legalBreakpoint = YES;
					}
				
				/* should this recombination event be counted in the expected number of recombinations E(R)? */
				/* for E(R) count only events with breakpoints as 1|1, 1|0 or 0|1  (i.e., not 0|0)         */
				/* if 1 represent a site that did found already its MRCA count it as a 0 */
				ThisBreakpIsTrapped = NO;
				if (doCountsForExpNumRec == YES)
					{
					p = nodes + activeGametes[whichInd];
					sizeNode = p->numSegNode;	
					if (CountsForExpNumRec (activeGametes, whichInd, whichSite, nodes, S_MRCA, sizeNode) == NO)
						{
						recNotToCount++;
						ThisBreakpIsTrapped = YES;
						/*fprintf(fpmpi,"..not to count.."); */
						}
					}
			
				/* copy whichIndividual to a new space in memory */
				hasPassedBreakPoint = NO;
			
			
				firstHalf = nextAvailable++; /* firstHalf is the first node that was created by the recombination */
				if (nextAvailable >= numNodes) /* if there aren't enough nodes it go into and it addition more */
					{
					/* ReallocNodes(&numNodes, activeGametes); */
					numNodes += INCREMENT_NODES;
					numTotalSegments += (INCREMENT_NODES*maxSegNode)+numNuc;
					
					segments = (TreeSegment *) realloc (segments, numTotalSegments  * (long) sizeof(TreeSegment)); 
					if (!segments)
						{
						fprintf (fpmpi, "Could not reallocate segments (%lu bytes)\n", ((numNodes*distance)+numNuc)  * (long) sizeof(TreeSegment));
						exit (1);
						}
					nodes = (TreeNode *) realloc (nodes, numNodes  * (long) sizeof(TreeNode));
					if (!nodes)
						{
						fprintf (fpmpi, "Could not reallocate nodes (%lu bytes)\n", numNodes  * (long) sizeof(TreeNode));
						exit (-1);
						}
					activeGametes = (int *) realloc (activeGametes, numNodes *(long) sizeof(int));
					if (!activeGametes)
						{
						fprintf (fpmpi, "Could not reallocate activeGametes (%lu bytes)\n", numNodes *(long) sizeof(int));
						exit (-1);
						}
					if (noisy == 4)
						fprintf (fpmpi, "\n\n...Doing reallocation of nodes (1)\n");
					}
									
				secondHalf = nextAvailable++; /* secondhalf is the second node that was created by the recombination */
				if (nextAvailable >= numNodes) /* if there aren't enough nodes it go into and it addition more */
					{
					/* ReallocNodes(&numNodes, activeGametes); */
					numNodes += INCREMENT_NODES;
					numTotalSegments += (INCREMENT_NODES*maxSegNode)+numNuc;
					
					segments = (TreeSegment *) realloc (segments, numTotalSegments  * (long) sizeof(TreeSegment));
					if (!segments)
						{
						fprintf (fpmpi, "Could not reallocate segments (%lu bytes)\n", ((numNodes*distance)+numNuc)  * (long) sizeof(TreeSegment));
						exit (1);
						}
					nodes = (TreeNode *) realloc (nodes, numNodes  * (long) sizeof(TreeNode));
					if (!nodes)
						{
						fprintf (fpmpi, "Could not reallocate nodes (%lu bytes)\n", numNodes  * (long) sizeof(TreeNode));
						exit (-1);
						}
					activeGametes = (int *) realloc (activeGametes, numNodes *(long) sizeof(int));
					if (!activeGametes)
						{
						fprintf (fpmpi, "Could not reallocate activeGametes (%lu bytes)\n", numNodes *(long) sizeof(int));
						exit (-1);
						}
					if (noisy == 4)
						fprintf (fpmpi, "\n\n...Doing reallocation of nodes (1)\n");
					}
								
				
				p = nodes + activeGametes[whichInd];
				q = nodes + firstHalf;	/* parent1 (new) */
				r = nodes + secondHalf; /* parent2 (new) */
				if (p->numSegNode > maxSegNode) /* checking */
					{
					fprintf (fpmpi, "\n\nWarning, too many segments in this node. max = %d", maxSegNode);
					fprintf (fpmpi, "\n p->numSegNode = %d",p->numSegNode);	
					exit(9);
					}
				q->index = firstHalf;
				r->index = secondHalf;

				if (doBranchNetfiles == YES)
					{
					q->NetLabelPrint = numNetLabelPrint;
					/*numNetLabelPrint++;*/
					r->NetLabelPrint = numNetLabelPrint;
					numNetLabelPrint++;
					}

				q->numSegNode = r->numSegNode = p->numSegNode;		/* Good if there are not nill segments.. then, in its case, it will be modify */
				q->indexOldMigPop = q->indexCurrentMigPop = p->indexCurrentMigPop; /* demes evolution */
				r->indexOldMigPop = r->indexCurrentMigPop = p->indexCurrentMigPop;


				q->class = 3;
				r->class = 3;
				q->GMRCA_ancestral = NO;
				r->GMRCA_ancestral = NO;
				q->breakp = whichSite;		
				r->breakp = whichSite;			
				q->time = currentTime;
				r->time = currentTime;
				q->sib = r;
				r->sib = q;
				q->left = p;
				r->left = p;
				p->anc1 = q;
				p->anc2 = r;				

				k = 0;
				for (w = 1; w < numSites; w++)
					{
					//fprintf (fpmpi, "\nstud[%d] = %d", w-1, stud[w-1]);
					if (whichSite == stud[w-1]) /* The breakpoints BETWEEN codons. "stud" is an array with the possible breakpoints beetween codons*/
						k++;

					/* exception for trapped material, any breakpoint here does not break material codons */
					if (ThisBreakpIsTrapped == YES)
						{
						k++;
						/*fprintf (fpmpi, "\n ThisBreakpIsTrapped \n");*/
						}
					}
				if (k == 0) /* broken codon, this indicates the position of rupture of the codon.. */
					{
					//q->breakCodon = YES;
					//r->breakCodon = YES;
					variable1 = whichSite/3.00 + 0.4;
					//fprintf (fpmpi, "\n variable1 = %lf", variable1);
					q->breakCodon = fabs(variable1);
					r->breakCodon = fabs(variable1);
					numREbreakCod++;
					variable2 = fmod(whichSite,3.00);
					if (variable2 == 0)
						q->whereBreakCodon = 2;
					else
						q->whereBreakCodon = 1;
					r->whereBreakCodon = 3;



					

					/* oct2009 */
					/* int			doBreakpBroken, LeftLess, LeftHigh, RightLess, RightHigh; */
					if (doCodonModel == YES)
						doBreakpBroken = YES;
					LeftLess = LeftHigh = RightLess = RightHigh = LeftLess2 = RightHigh2 = -1;

					if (q->whereBreakCodon == 1) /* first codon position breakp */
						{
						RightHigh = q->breakCodon * 3;
						RightLess = RightHigh - 1;
						LeftLess = RightLess - 1; 
						LeftHigh = RightLess - 1;
						
						
						q->SitesNonAncHere[LeftLess+1] = 1;
						q->SitesNonAncHere[LeftLess+2] = 1;
						r->SitesNonAncHere[RightLess-1] = 1;

						/*fprintf (fpmpi, "\n LeftLess+2 = %d; RightLess-1 = %d \n", LeftLess+2, RightLess-1);*/
						for (mmm=1; mmm<=numNuc; mmm++)
							{
							if (mmm <= RightLess-1)
								q->SitesNonAncHere[mmm] = p->SitesNonAncHere[mmm]; 
							if (mmm > RightLess-1 && mmm <= LeftLess+2)
								q->SitesNonAncHere[mmm] = 1;
							if (mmm > LeftLess+2)
								q->SitesNonAncHere[mmm] = -1;

							if (mmm < RightLess-1)
								r->SitesNonAncHere[mmm] = -1; 
							if (mmm == RightLess-1)
								r->SitesNonAncHere[mmm] = 1;
							if (mmm > RightLess-1)
								r->SitesNonAncHere[mmm] =  p->SitesNonAncHere[mmm];
							}	
						}
					else if (q->whereBreakCodon == 2) /* second codon position breakp */
						{
						RightHigh = q->breakCodon * 3;
						RightLess = q->breakCodon * 3;
						LeftHigh = RightLess - 1;
						LeftLess = LeftHigh - 1;


						q->SitesNonAncHere[LeftHigh+1] = 1;
						r->SitesNonAncHere[RightLess-1] = 1;	
						r->SitesNonAncHere[RightLess-2] = 1;
						
						/*fprintf (fpmpi, "\n LeftHigh+1 = %d; RightLess-2 = %d \n", LeftHigh+1, RightLess-2);*/
						for (mmm=1; mmm<=numNuc; mmm++)
							{
							if (mmm < LeftHigh+1)
								q->SitesNonAncHere[mmm] = p->SitesNonAncHere[mmm]; 
							if (mmm == LeftHigh+1)
								q->SitesNonAncHere[mmm] = 1;
							if (mmm > LeftHigh+1)
								q->SitesNonAncHere[mmm] = -1;

							if (mmm >= LeftHigh+1)
								r->SitesNonAncHere[mmm] = p->SitesNonAncHere[mmm]; 
							if (mmm >= RightLess-2 && mmm < LeftHigh+1)
								r->SitesNonAncHere[mmm] = 1;
							if (mmm < RightLess-2)
								r->SitesNonAncHere[mmm] = -1;

							}	
						}
					else
						{
						fprintf (fpmpi, "error at q->whereBreakCodon intra codon Rec _ Main (%d != 1 or 2)\n", q->whereBreakCodon);
						exit (-1);
						}
					/* fprintf (fpmpi, "\n Left (Less-High) %d-%d; Right (Less-High) %d-%d \n", LeftLess, LeftHigh, RightLess, RightHigh);	*/				
					if (LeftLess == -1 || LeftHigh == -1 || RightLess == -1 || RightHigh == -1)
						{
						fprintf (fpmpi, "\n Error (value = -1): Left (Less-High) %d-%d; Right (Less-High) %d-%d \n", LeftLess, LeftHigh, RightLess, RightHigh);					
						exit (-1);
						}


					}
				if (k > 0 || doCodonModel == NO) /* inter codon rec */
					{
					doBreakpBroken = NO;

					for (mmm=1; mmm<=numNuc; mmm++)
						{
						if (mmm >= whichSite)
							{
							q->SitesNonAncHere[mmm] = -1; /* non anc mat */
							}
						if (mmm < whichSite)
							{
							q->SitesNonAncHere[mmm] = p->SitesNonAncHere[mmm]; /* non anc mat */
							}

						if (mmm < whichSite)
							{
							r->SitesNonAncHere[mmm] = -1; /* non anc mat */
							}
						if (mmm >= whichSite)
							{
							r->SitesNonAncHere[mmm] = p->SitesNonAncHere[mmm]; /* non anc mat */
							}
						}	
					}
				/* oct2009 */

				/*for (mmm=1; mmm<=numNuc; mmm++)
					{
					fprintf (fpmpi, "\n Initial. Site %d. Node: %d, p->SitesNonAncHere = %d \n", mmm, p->index, p->SitesNonAncHere[mmm]);
					}
				for (mmm=1; mmm<=numNuc; mmm++)
					{
					fprintf (fpmpi, "\n Site %d. Node: %d, q->SitesNonAncHere = %d;  Node %d, r->SitesNonAncHere = %d \n", mmm, q->index, q->SitesNonAncHere[mmm], r->index, r->SitesNonAncHere[mmm]);
					}*/
				

				/*fprintf (fpmpi, "\n r->index = %d, r->time = %lf, r->class = %d, r->breakp = %d, r->breakCodon = %d, r->whereBreakCodon = %d", r->index, r->time, r->class, r->breakp, r->breakCodon, r->whereBreakCodon);
				fprintf (fpmpi, "\n q->index = %d, q->time = %lf, q->class = %d, q->breakp = %d, q->breakCodon = %d, q->whereBreakCodon = %d\n", q->index, q->time, q->class, q->breakp, q->breakCodon, q->whereBreakCodon);*/
				k = 0;	

			


				if (noisy == 4)
					{
					fprintf (fpmpi, "\nNode index %d with breakpoint on %d site", p->index, whichSite);
					fprintf (fpmpi, "\nThis node contains %d fragment(s):", p->numSegNode);
					}
				
				for (w = 0; w < p->numSegNode; w++)
					{
					s = segments + post(w,p->index,distance);
					if (post(w,p->index,distance) > numTotalSegments) /* checking */
						{
						fprintf (fpmpi, "\n post = %d > numTotalSegments = %d", post(w,p->index,distance), numTotalSegments);
						exit (-7);
						}
					if (noisy == 4)
						{
						fprintf (fpmpi, "\ns->sIndex = %d", s->sIndex);
						fprintf (fpmpi, "\ns->sStart = %d", s->sStart);
						fprintf (fpmpi, "\ns->sEnd = %d",s->sEnd);
						}
					}
				if (noisy == 4)
					fprintf (fpmpi, "\n\n>> Process evolution..");
					
					
				a = b = aa = bb = aaa = bbb = out = 0;
				startsVectorRec = (int *) calloc((p->numSegNode),(long) sizeof(int));
				if (!startsVectorRec)
					{
					fprintf (fpmpi, "Could not allocate startsVectorRec (%lu bytes)\n", (p->numSegNode) *(long) sizeof(int));	
					exit (-1);
					}
				endsVectorRec = (int *) calloc((p->numSegNode),(long) sizeof(int));
				if (!endsVectorRec)
					{
					fprintf (fpmpi, "Could not allocate endsVectorRec (%lu bytes)\n", (p->numSegNode) *(long) sizeof(int));
					exit (-1);
					}
			
				for (i = 0; i < p->numSegNode; i++)					/* for each segment */
					{
					s = segments + post(i,p->index,distance);										
					startsVectorRec[i] = s->sStart;
					endsVectorRec[i] = s->sEnd;
				
					/* first half */
					if (s->sStart >= whichSite)
						{
						if (noisy == 4)
							fprintf (fpmpi, "\nNil segment in Left, don't make it");
							
						q->numSegNode--;
						a++;
						out = 1;
						} 
				
					if (s->sStart == 1 && s->sEnd >= (whichSite-1) && out == 0)
						{
						if (aa > 0)
							{
							if (noisy == 4)
								fprintf (fpmpi, "\nRepit segment in Left, don't make it");
								
							q->numSegNode--;
							out = 1;
							}
						aa++;
						}
					
					for (w = 0; w < i+1; w++)
						{
						if (out == 0 && s->sStart != 1 && s->sStart == startsVectorRec[w] && s->sStart < whichSite && w != i && s->sEnd >= (whichSite-1) && endsVectorRec[w] >= (whichSite-1) && startsVectorRec[w] != 0)
							{
							if (noisy == 4)
								fprintf (fpmpi, "\nRepit segment in Left, don't make it");
								
							q->numSegNode--;
							aaa++;
							out = 1;
							}
						}
					
					if (out == 0 && s->sStart < whichSite)
						{
						if (aa == 0 && aaa == 0)
							{							
							if (post(i-a-aa-aaa,q->index,distance) > numTotalSegments) /* checking */
								{
								fprintf (fpmpi, "\n post = %d > numTotalSegments = %d", post(i-a-aa-aaa,q->index,distance), numTotalSegments);
								exit (-7);
								}
							
							n = segments + post(i-a-aa-aaa,q->index,distance);
							n->sIndexNode = firstHalf;
							}
						if (aa > 0 && aaa == 0)
							{
							if (post(i+1-a-aa-aaa,q->index,distance) > numTotalSegments) /* checking */
								{
								fprintf (fpmpi, "\n post = %d > numTotalSegments = %d", post(i+1-a-aa-aaa,q->index,distance), numTotalSegments);
								exit (-7);
								}
							
							n = segments + post(i+1-a-aa-aaa,q->index,distance);
							n->sIndexNode = firstHalf;
							}
						if (aaa > 0 && aa > 0)
							{
							if (post(i+1-a-aa-aaa,q->index,distance) > numTotalSegments) /* checking */
								{
								fprintf (fpmpi, "\n post = %d > numTotalSegments = %d", post(i+1-a-aa-aaa,q->index,distance), numTotalSegments);
								exit (-7);
								}
							
							n = segments + post(i+1-a-aa-aaa,q->index,distance);
							n->sIndexNode = firstHalf;
							}
						if (aaa > 0 && aa == 0)
							{
							if (post(i+0-a-aa-aaa,q->index,distance) > numTotalSegments) /* checking */
								{
								fprintf (fpmpi, "\n post = %d > numTotalSegments = %d", post(i+0-a-aa-aaa,q->index,distance), numTotalSegments);
								exit (-7);
								}
							
							n = segments + post(i+0-a-aa-aaa,q->index,distance);
							n->sIndexNode = firstHalf;
							}
						/* n = s; initial, the new segments are similar to the old segments */
						nodeValue = firstHalf;


						if (doBreakpBroken == YES && s->sStart < whichSite && s->sEnd >= whichSite)
							{
							w = recSegmentsGeneratesLeftBrokenCodon(nodeValue, s, n, numNuc, whichSite, LeftLess, RightHigh, &actSegIndex); /* it makes the segments of the left node */
							LeftLess2 = LeftLess;
							RightHigh2 = RightHigh;
							}
						else
							{
							w = recSegmentsGeneratesLeft(nodeValue, s, n, numNuc, whichSite, &actSegIndex); /* it makes the segments of the left node */
							}


						if (w != 1)
							{
							fprintf (fpmpi, "Warning in recSegmentsGeneratesLeft");
							exit (-1);
							}
						actNumSegments++;
						if (n->sStart == 0 && n->sEnd == 0) /* unreal segments */
							{
							fprintf (fpmpi, "\nNot to be here. Segment Left start and end = 0");
							
							n->before1 = n->before2 = n->after1 = n->after2 = NULL;
							actNumSegments--;
							exit(-1);
							}
						if (noisy == 4)			
							{
							fprintf (fpmpi, "\nAfter rec. left");
							fprintf (fpmpi, "\nq->seg->sIndex = %d", n->sIndex);
							fprintf (fpmpi, "\nq->seg->sStart = %d", n->sStart);
							fprintf (fpmpi, "\nq->seg->sEnd = %d\n", n->sEnd);
							}
						}
					w = out = 0;
				
					/* second half */
					if (s->sEnd < whichSite)
						{
						if (noisy == 4)
							fprintf (fpmpi, "\nNil segment in Right, don't make it");
							
						r->numSegNode--;
						b++;
						out = 1;
						}
					if (s->sStart <= whichSite && s->sEnd == numNuc && out == 0)
						{
						if (bb > 0)
							{
							if (noisy == 4)
								fprintf (fpmpi, "\nRepit segment in Right, don't make it");
								
							r->numSegNode--;
							out = 1;
							}
						bb++;
						}
					for (w = 0; w < i+1; w++)
						{
						if (out == 0 && s->sEnd != numNuc && s->sEnd == endsVectorRec[w] && s->sEnd > whichSite && w != i && s->sStart <= whichSite && startsVectorRec[w] <= whichSite && endsVectorRec[w] != 0)
							{
							if (noisy == 4)
								fprintf (fpmpi, "\nRepit segment in Right, don't make it");
								
							r->numSegNode--;
								
							bbb++;
							out = 1;
							}
						}
					for (w = 0; w < i+1; w++)
						{
						if (out == 0 && s->sEnd == endsVectorRec[w] && s->sEnd >= whichSite && w != i && s->sStart <= whichSite && startsVectorRec[w] <= whichSite && endsVectorRec[w] != 0)
							{
							if (noisy == 4)
								fprintf (fpmpi, "\nRepit segment in Right, don't make it");
								
							r->numSegNode--;
								
							bbb++;
							out = 1;
							}
						}
					if (out == 0 && s->sEnd >= whichSite)
						{
						if (bb == 0 && bbb == 0)
							{
							if (post(i-b-bb-bbb,r->index,distance) > numTotalSegments) /* checking */
								{
								fprintf (fpmpi, "\n post = %d > numTotalSegments = %d", post(i-b-bb-bbb,r->index,distance), numTotalSegments);
								exit (-7);
								}
							
							m = segments + post(i-b-bb-bbb,r->index,distance);
							m->sIndexNode = secondHalf;
							}
						if (bb > 0 && bbb == 0)
							{
							if (post(i+1-b-bb-bbb,r->index,distance) > numTotalSegments) /* checking */
								{
								fprintf (fpmpi, "\n post = %d > numTotalSegments = %d", post(i+1-b-bb-bbb,r->index,distance), numTotalSegments);
								exit (-7);
								}
							
							m = segments + post(i+1-b-bb-bbb,r->index,distance);
							m->sIndexNode = secondHalf;
							}
						if (bb > 0 && bbb > 0)
							{
							if (post(i+1-b-bb-bbb,r->index,distance) > numTotalSegments) /* checking */
								{
								fprintf (fpmpi, "\n post = %d > numTotalSegments = %d", post(i+1-b-bb-bbb,r->index,distance), numTotalSegments);
								exit (-7);
								}
							
							m = segments + post(i+1-b-bb-bbb,r->index,distance);
							m->sIndexNode = secondHalf;
							}
						if (bb == 0 && bbb > 0)
							{
							if (post(i+0-b-bb-bbb,r->index,distance) > numTotalSegments) /* checking */
								{
								fprintf (fpmpi, "\n post = %d > numTotalSegments = %d", post(i+0-b-bb-bbb,r->index,distance), numTotalSegments);
								exit (-7);
								}
							
							m = segments + post(i+0-b-bb-bbb,r->index,distance);
							m->sIndexNode = secondHalf;
							}
						/* m = s; initial, the new segments are similar to the old segments */
						nodeValue = secondHalf;


						if (doBreakpBroken == YES && s->sStart < whichSite && s->sEnd >= whichSite)
							{
							w = recSegmentsGeneratesRightBrokenCodon(nodeValue, s, m, numNuc, whichSite, LeftLess, RightHigh, &actSegIndex);
							LeftLess2 = LeftLess;
							RightHigh2 = RightHigh;
							}
						else
							{
							w = recSegmentsGeneratesRight(nodeValue, s, m, numNuc, whichSite, &actSegIndex);
							}


						if (w != 1)
							{
							fprintf (fpmpi, "Warning in recSegmentsGeneratesRight");
							exit (-1);
							}
					
						actNumSegments++;
						if (m->sStart == 0 && m->sEnd == 0) /* unreal segments */
							{
							fprintf (fpmpi, "\nNot to be here. Segment Right start and end = 0");
							
							m->before1 = m->before2 = m->after1 = m->after2 = NULL;
							actNumSegments--;
							exit(-1);
							}
						if (noisy == 4)
							{
							fprintf (fpmpi, "\nAfter rec. right");
							fprintf (fpmpi, "\nr->seg->sIndex = %d", m->sIndex);
							fprintf (fpmpi, "\nr->seg->sStart = %d", m->sStart);
							fprintf (fpmpi, "\nr->seg->sEnd = %d\n", m->sEnd);
							} 
						}
					w = out = 0;
					}
			
				free (startsVectorRec);
				free (endsVectorRec);
				a = b = aa = bb = aaa = bbb = 0;

				/* intra codon breakpoints readjust MRCA oct2009 */
				if (doBreakpBroken == YES) /* The codon positions that recombinaed are increased a unit */
					{
					for (w = 1; w <= numNuc; w++)
						{
						if (w >= LeftLess2 && w <= RightHigh2)
							{
							S_MRCA[w]++;
							/*fprintf (fpmpi, "\n Increasing MRCA to %d: now is: %d \n", w, S_MRCA[w]);*/
							}
						}

					}
				/* oct2009 */


				if (noisy == 4)
					{
					fprintf (fpmpi, "\n>> Recombination Results:");
					fprintf (fpmpi, "\nNew left node with %d fragment(s)", q->numSegNode);
					}
				for (w = 0; w < q->numSegNode;w++)
					{
					n = segments + post(w,q->index,distance);
					if (noisy == 4)
						{
						fprintf (fpmpi, "\nq->seg->sIndex = %d", n->sIndex);
						fprintf (fpmpi, "\nq->seg->sStart = %d", n->sStart);
						fprintf (fpmpi, "\nq->seg->sEnd = %d\n", n->sEnd);
						} 
					}
				
				if (noisy == 4)
					fprintf (fpmpi, "\nNew right node with %d fragment(s)", r->numSegNode);
					
				for (w = 0; w < r->numSegNode;w++)
					{
					m = segments + post(w,r->index,distance);
					if (noisy == 4)
						{
						fprintf (fpmpi, "\nr->seg->sIndex = %d", m->sIndex);
						fprintf (fpmpi, "\nr->seg->sStart = %d", m->sStart);
						fprintf (fpmpi, "\nr->seg->sEnd = %d\n", m->sEnd);
						} 
					}
				if (noisy > 3)
					{
					fprintf (fpmpi,"\n");
					/*fprintf (fpmpi, "the node is whichInd = %d, and the site is whichSite = %d", whichInd+1, whichSite+1);*/
					}
				if (noisy > 1)
					{
					fprintf (fpmpi, "Recombination involving %d (copied to %d and %d) in deme %d", p->index, q->index, r->index, p->indexCurrentMigPop);
					fprintf (fpmpi, "\n Breakpoint was at site %d", whichSite);
					}
				if (noisy > 3)
					fprintf (fpmpi,"\n");
				
				if (doBranchNetfiles == YES)
					{
					fprintf(fpBranchNet,"%d %d\n", q->NetLabelPrint, p->NetLabelPrint);
					/*fprintf(fpBranchNet,"%d %d\n", r->NetLabelPrint, p->NetLabelPrint);*/
					fprintf (fpmpi, "\nNET INFORMATION Recombination involving %d (copied to %d and %d) in deme %d", p->NetLabelPrint, q->NetLabelPrint, r->NetLabelPrint, p->indexCurrentMigPop);
					}
					
					
				breakpoint[numRE-1] = whichSite;  /* the breakpoint site is call breakpoint. breakpoint[0] = 7, breakpoint[1] = 96, breakpoint[2] = 187.. */
				
				/* readjust active sites */
				activeGametes[whichInd] = firstHalf;	/* new active nodes firstHalf and secondHalf */
				activeGametes[numActiveGametes++] = secondHalf;		/* there are 1 active node more (in recombination) */
			
				w = 0;	
				for (k = 1; k <= numPopulations+currentBigDeme; k++)
					{
					if (k == p->indexCurrentMigPop)
						numParcialActiveGametes[k] = numParcialActiveGametes[k] +1;
					w = w + numParcialActiveGametes[k];
					}
				if (w != numActiveGametes)
					{
					fprintf (fpmpi, "\n 4The sum of partial active gametes is different than total gametes number, w %d != numActiveGametes %d. In choose the deme", w, numActiveGametes);
					exit (-1);
					}
				}
		
		
			/*** MIGRATION ***/
			if (isMigration == YES)
				{
				if (noisy == 4)
					fprintf (fpmpi, "\n* Migration *");			
					
				/* Which node has the migration? */
				cumPopulPart = (double *) calloc((numParcialActiveGametes[whichDeme]+1),(long) sizeof(double));
				if (!cumPopulPart)
					{
					fprintf (fpmpi, "Could not allocate cumPopulPart (%lu bytes)\n", (numParcialActiveGametes[whichDeme]+1) *(long) sizeof(double));
					exit (-1);
					}
			
				for (k = 1; k <= numParcialActiveGametes[whichDeme]; k++)
					cumPopulPart[k] = 0;
				cumPopulPart[0] = 0;
				w = 0;
			
				for (k = 1; k <= numParcialActiveGametes[whichDeme]; k++)
					cumPopulPart[k] = cumPopulPart[k-1] + 1.0/numParcialActiveGametes[whichDeme];
				/*fprintf(fpmpi, "\n cumPopulPart[numParcialActiveGametes[whichDeme]] = %lf", cumPopulPart[numParcialActiveGametes[whichDeme]]);*/
				for (k = 1; k <= numPopulations+currentBigDeme; k++)
					w = w + numParcialActiveGametes[k];				
					
				if (w != numActiveGametes)
					{
					fprintf (fpmpi, "\n 5The sum of partial active gametes is different than total gametes number, w %d != numActiveGametes %d. In choose the deme", w, numActiveGametes);
					exit (-1);
					}
				
				ran = RandomUniform(seed);
				whichInd = bbinDemes(ran, cumPopulPart, numParcialActiveGametes[whichDeme]);
				w = 0;
				for (i = 0; i < numActiveGametes; i++)
					{
					p = nodes + activeGametes[i];
			
					if (p->indexCurrentMigPop == whichDeme)
						w++;
						
					if (w == whichInd)
						{
						whichInd = i;
						break;
						}
					}
				free (cumPopulPart);
			
				/* Which arrived Deme? */			
				for (k = 1; k <= numPopulations+currentBigDeme; k++)
					cumPopulTase[k] = 0;
				cumPopulTase[0] = 0;
				w = c = e = 0;
				
				if (doConvergDemes == YES)
					{
					ArrivedDemesOptions = (int *)calloc((currentDemesNumber+1),(long) sizeof(int));
					if (!ArrivedDemesOptions)
						{
						fprintf (fpmpi, "Could not allocate ArrivedDemesOptions (%lu bytes)\n", (currentDemesNumber+1)  * (long) sizeof(int));
						exit (1);
						}
					for (k = 1; k <= currentDemesNumber; k++)
						ArrivedDemesOptions[k] = 0;
					ArrivedDemesOptions[0] = 0;

					e = 1;
					w = 1;
					for (k = 1; k <= 2*numPopulations-1; k++)
						{
						if (CurrentDemesState[k] != 0)
							{
							ArrivedDemesOptions[w] = CurrentDemesState[k];
							w++;
							}
						}
					if (w != currentDemesNumber+1) /* Checking */
						{
						fprintf (fpmpi, "\n Warning in currentDemesNumber: w = %d, currentDemesNumber+1 = %d\n", w, currentDemesNumber+1);
						exit (1);
						}
					
					/*for (k = 0; k < numActiveGametes; k++)
						{
						p = nodes + activeGametes[k];
						if (k == 0)
							{
							ArrivedDemesOptions[e] = p->indexCurrentMigPop;
							e++;
							}
						else
							{
							c = 0;
							for (w = 1; w <= currentDemesNumber; w++)
								if (p->indexCurrentMigPop != ArrivedDemesOptions[w])
									c++;
							if (c == currentDemesNumber)
								{
								ArrivedDemesOptions[e] = p->indexCurrentMigPop;
								e++;
								}
							}
							
						fprintf (fpmpi, " \nNode %d of deme %d", p->index, p->indexCurrentMigPop);
						}*/		
					
					/*for (k = 1; k <= currentDemesNumber; k++)
						fprintf (fpmpi, " \n ArrivedDemesOptions[%d] = %d", k, ArrivedDemesOptions[k]);*/ /* To see the array */
				
					/*if (e != currentDemesNumber+1)
						{
						fprintf (fpmpi, "\n Warning in currentDemesNumber: e = %d, currentDemesNumber+1 = %d\n", e, currentDemesNumber+1);
						exit (1);
						}*/ // This is not an error
					w = e = c = 0;
					for (k = 1; k <= numPopulations+currentBigDeme; k++)
						{
						if (k == whichDeme)
							cumPopulTase[k] = cumPopulTase[k-1];
						else
							{
							c = 0;
							for (w = 1; w <= currentDemesNumber; w++)
								if (k == ArrivedDemesOptions[w])
									c++;
							if (c > 0)
								cumPopulTase[k] = cumPopulTase[k-1] + 1.0/(currentDemesNumber-1);			
							}
						}
					}
				else
					{
					for (k = 1; k <= numPopulations+currentBigDeme; k++)
						{
						if (k == whichDeme)
							cumPopulTase[k] = cumPopulTase[k-1];
						 else
							cumPopulTase[k] = cumPopulTase[k-1] + 1.0/(numPopulations+currentBigDeme-1);
						}					
					}
				/*fprintf(fpmpi, "\n cumPopulTase[numPopulations+currentBigDeme] (2) = %lf", cumPopulTase[numPopulations+currentBigDeme]);*/

				ran = RandomUniform(seed);
				arrivedDeme = bbinDemes(ran, cumPopulTase, numPopulations+currentBigDeme);
				/*fprintf(fpmpi,"\n To deme %d  \n", arrivedDeme);*/
				w = 0;
				
				if (whichDeme == arrivedDeme) /* Cheking */
					{
					fprintf (fpmpi, "\n\nERROR in migration, whichDeme == arrivedDeme");
					exit (-1);
					}
				if (doConvergDemes == YES)
					free (ArrivedDemesOptions);
				/* evolution migration node*/
				p = nodes + activeGametes[whichInd];
				if (p->indexCurrentMigPop != whichDeme)
					{
					fprintf (fpmpi, "\n\nERROR in migration, p->indexCurrentMigPop != whichDeme");
					exit (-1);
					}
			
				p->indexCurrentMigPop = arrivedDeme; /* Changing the deme */
			
				if (noisy == 4)
					{
					fprintf (fpmpi, "\nNode %d migrates from deme %d to deme %d", p->index, whichDeme, arrivedDeme);
					fprintf (fpmpi, "\nThis node contains %d fragment(s):", p->numSegNode);	
					}
				for (w = 0; w < p->numSegNode; w++)
					{
					s = segments + post(w,p->index,distance);
					if (noisy == 4)
						{
						fprintf (fpmpi, "\ns->sIndex = %d", s->sIndex);
						fprintf (fpmpi, "\ns->sStart = %d", s->sStart);
						fprintf (fpmpi, "\ns->sEnd = %d",s->sEnd);
						}
					}
				/*numParcialActiveGametes[k] ++ y --*/
				numParcialActiveGametes[whichDeme] = numParcialActiveGametes[whichDeme]-1;
				numParcialActiveGametes[arrivedDeme] = numParcialActiveGametes[arrivedDeme]+1;								
				k = 0;
				
					
				if (noisy == 4)
					fprintf (fpmpi, "\n");
				if (noisy > 1)
					fprintf (fpmpi, "Migration involving node %d from deme %d to deme %d", p->index, whichDeme, arrivedDeme);
				if (noisy == 4)
					fprintf (fpmpi, "\n");
				}
		
		
			/*** COALESCENCE ***/
			if (isCoalescence == YES)
				{
				if (noisy == 4)
					fprintf (fpmpi, "\n* Coalescence *\n");
					
				/* figure out which two nodes are involved */ 			
				/* intial nodes: firstInd and secondInd (they are the descendants). newInd is the ancestral node, is the new node to make */
				cumPopulPart = (double *) calloc((numParcialActiveGametes[whichDeme]+1),(long) sizeof(double));
				if (!cumPopulPart)
					{
					fprintf (fpmpi, "Could not allocate cumPopulPart (%lu bytes)\n", (numParcialActiveGametes[whichDeme]+1) *(long) sizeof(double));
					exit (-1);
					}
				
				for (k = 1; k <= numParcialActiveGametes[whichDeme]; k++)
					cumPopulPart[k] = 0;
				cumPopulPart[0] = 0;
				w = 0;
			
				for (k = 1; k <= numParcialActiveGametes[whichDeme]; k++)
					cumPopulPart[k] = cumPopulPart[k-1] + 1.0/numParcialActiveGametes[whichDeme];
					
				for (k = 1; k <= numPopulations+currentBigDeme; k++)
					w = w + numParcialActiveGametes[k];
					
				if (w != numActiveGametes)
					{
					fprintf (fpmpi, "\n 6The sum of partial active gametes is different than total gametes number, w %d != numActiveGametes %d. In choose the deme", w, numActiveGametes);
					exit (-1);
					}
				
				ran = RandomUniform(seed);
				firstInd = bbinDemes(ran, cumPopulPart, numParcialActiveGametes[whichDeme]);
				w = 0;
				for (i = 0; i < numActiveGametes; i++)
					{
					p = nodes + activeGametes[i];
			
					if (p->indexCurrentMigPop == whichDeme)
						w++;
						
					if (w == firstInd)
						{
						firstInd = i;
						break;
						}
					}
				if (firstInd >= numActiveGametes) /* checking */
					{
					fprintf (fpmpi, "\n\nERROR: firstInd out of range!\n");
					exit (-1);
					}
					
				do
					{
					for (k = 1; k <= numParcialActiveGametes[whichDeme]; k++)
						cumPopulPart[k] = 0;
					cumPopulPart[0] = 0;
					w = 0;
			
					for (k = 1; k <= numParcialActiveGametes[whichDeme]; k++)
						cumPopulPart[k] = cumPopulPart[k-1] + 1.0/numParcialActiveGametes[whichDeme];
						
					for (k = 1; k <= numPopulations+currentBigDeme; k++)
						w = w + numParcialActiveGametes[k];				
						
					if (w != numActiveGametes)
						{
						fprintf (fpmpi, "\n 7The sum of partial active gametes is different than total gametes number, w %d != numActiveGametes %d. In choose the deme", w, numActiveGametes);
						exit (-1);
						}
				
					ran = RandomUniform(seed);
					secondInd = bbinDemes(ran, cumPopulPart, numParcialActiveGametes[whichDeme]);
					w = 0;
					for (i = 0; i < numActiveGametes; i++)
						{
						p = nodes + activeGametes[i];
			
						if (p->indexCurrentMigPop == whichDeme)
							w++;
							
						if (w == secondInd)
							{
							secondInd = i;
							break;
							}
						}
					} while (firstInd == secondInd);
				free (cumPopulPart);			
			
				newInd = nextAvailable;
				if (noisy > 1)
					fprintf (fpmpi, "Coalescence involving %d and %d to create node %d in deme %d", activeGametes[firstInd], activeGametes[secondInd], newInd, whichDeme);
					
				p = nodes + activeGametes[firstInd];
				q = nodes + activeGametes[secondInd];
				
				if (p->numSegNode > maxSegNode) 
					{
					fprintf (fpmpi, "\n\nWarning, too many segments in this node. max = %d", maxSegNode);
					fprintf (fpmpi, "\n p->numSegNode = %d",p->numSegNode);
					exit(9);
					}
				if (q->numSegNode > maxSegNode)
					{
					fprintf (fpmpi, "\n\nWarning, too many segments in this node. max = %d", maxSegNode);
					fprintf (fpmpi, "\n q->numSegNode = %d",q->numSegNode);
					exit(9);
					}

				r = nodes + newInd;		/* new ancester */
				r->index = nextAvailable;
				r->label = labelNodes++;
				r->indexOldMigPop = r->indexCurrentMigPop = whichDeme;
				
				r->breakp = NO;
				r->breakCodon = NO;
				r->class = 4;
				r->GMRCA_ancestral = NO;
				/*fprintf (fpmpi, "\nobtained r->index = %d: r->breakp = %d, r->breakCodon = %d, r->class = %d \n", r->index, r->breakp, r->breakCodon, r->class);*/
				
				for (mmm = 1; mmm <= numNuc; mmm++)	 /* for each segment of p node going to r node*/
					{
					sigue = 0;
					stateHere_P = -2; /* -2, non ancestral; 0 ancestral */
					stateHere_Q = -2; /* -2, non ancestral; 0 ancestral */

					for (i = 0; i < p->numSegNode; i++)	 /* for each segment of p node going to r node*/
						{
						s = segments + post(i,p->index,distance);
						if (mmm >= s->sStart && mmm <= s->sEnd) /* is ancestral material */
							{
							stateHere_P = 0;
							}
						}
					for (i = 0; i < q->numSegNode; i++)	/* for each segment of q node going to r node */
						{
						n = segments + post(i,q->index,distance);
						if (mmm >= n->sStart && mmm <= n->sEnd) /* is ancestral material */
							{
							stateHere_Q = 0;
							}
						}
					
					
					/*fprintf (fpmpi, "\n Here(%d) stateHere_P = %d and stateHere_Q = %d ", mmm, stateHere_P, stateHere_Q);
					fprintf (fpmpi, "\n  Here(%d), node %d: p->SitesNonAncHere[mmm] = %d && node %d: q->SitesNonAncHere[mmm] = %d \n", mmm, p->index, p->SitesNonAncHere[mmm], q->index, q->SitesNonAncHere[mmm]);*/

					if (stateHere_P < 0 && stateHere_Q < 0) /* non ancestral material */
						{
						r->SitesNonAncHere[mmm] = -1;
						sigue++;
						/*fprintf (fpmpi, "1Position %d is NON anc mat \n", mmm);*/
						}
					if (stateHere_P >= 0 || stateHere_Q >= 0) /* ancestral material (inc pseudo) */
						{
						if (p->SitesNonAncHere[mmm] == 1 && q->SitesNonAncHere[mmm] == 1 && sigue == 0) /* pseudo anc mat */
							{
							r->SitesNonAncHere[mmm] = 1;
							sigue++;
							/*fprintf (fpmpi, "2Position %d is PSEUDO anc mat \n", mmm);*/
							}
						if (p->SitesNonAncHere[mmm] == 1 && stateHere_Q < 0 && sigue == 0) /* pseudo anc mat */
							{
							r->SitesNonAncHere[mmm] = 1;
							sigue++;
							/*fprintf (fpmpi, "3Position %d is PSEUDO anc mat \n", mmm);*/
							}
						if (q->SitesNonAncHere[mmm] == 1 && stateHere_P < 0 && sigue == 0) /* pseudo anc mat */
							{
							r->SitesNonAncHere[mmm] = 1;
							sigue++;
							/*fprintf (fpmpi, "4Position %d is PSEUDO anc mat \n", mmm);*/
							}
						if (sigue == 0)
							{
							if (p->SitesNonAncHere[mmm] == 0 || q->SitesNonAncHere[mmm] == 0)  /* anc mat */
								{
								r->SitesNonAncHere[mmm] = 0;
								sigue++;
								/*fprintf (fpmpi, "5Position %d is ANC mat \n", mmm);*/
								}
							}
						}
					}

				/*fprintf (fpmpi, "\nobtained r->index = %d: r->breakp = %d, r->breakCodon = %d, r->class = %d \n", r->index, r->breakp, r->breakCodon, r->class);*/



				if (doBranchNetfiles == YES)
					{
					r->NetLabelPrint = numNetLabelPrint;
					numNetLabelPrint++;

					if (p->NetLabelPrint == q->NetLabelPrint)
						{
						fprintf(fpBranchNet,"%d %d\n", r->NetLabelPrint, p->NetLabelPrint);
						}
					else
						{
						fprintf(fpBranchNet,"%d %d\n", r->NetLabelPrint, p->NetLabelPrint);
						fprintf(fpBranchNet,"%d %d\n", r->NetLabelPrint, q->NetLabelPrint);
						}
					fprintf (fpmpi, "\nNET INFORMATION Coalescence involving %d and %d to create node %d in deme %d", p->NetLabelPrint, q->NetLabelPrint, r->NetLabelPrint, whichDeme);
					}
				

				coalVectorCountStarts = (int *) calloc((p->numSegNode+q->numSegNode),(long) sizeof(int));
				if (!coalVectorCountStarts)
					{
					fprintf (fpmpi, "Could not allocate coalVectorCountStarts (%lu bytes)\n", (p->numSegNode+q->numSegNode) *(long) sizeof(int));
					exit (-1);
					}
				coalVectorCountEnds = (int *) calloc((p->numSegNode+q->numSegNode),(long) sizeof(int));
				if (!coalVectorCountEnds)
					{
					fprintf (fpmpi, "Could not allocate coalVectorCountEnds (%lu bytes)\n", (p->numSegNode+q->numSegNode) *(long) sizeof(int));
					exit (-1);
					}
				
				for (i = 0; i < p->numSegNode; i++)	 /* for each segment of p node going to r node */
					{
					s = segments + post(i,p->index,distance);
					if (post(i,r->index,distance) > numTotalSegments) /* control */
						{
						fprintf (fpmpi, "\n post = %d > numTotalSegments = %d", post(w,r->index,distance), numTotalSegments);
						exit (-7);
						}
					
					m = segments + post(i,r->index,distance);					/* new ancester */
					m->sIndexNode = newInd;
					
					s->before1 = m;
					m->before1 = NULL;
					m->before2 = NULL;
					m->after1 = s;
					m->after2 = NULL;
					m->sIndex = actSegIndex;
					m->sStart = s->sStart;
					m->sEnd = s->sEnd;				

					actSegIndex++;
					actNumSegments++;
					/*r->numSegNode++;*/
						
					if (m->sStart == 0 && m->sEnd == 0) /* unreal segments */
						{
						fprintf (fpmpi, "\nNot to be here. COAL1, segment start and end = 0");
						
						m->before1 = m->before2 = m->after1 = m->after2 = NULL;
						actNumSegments--;
						exit(-1);
						}
					coalVectorCountStarts[i] = s->sStart;
					coalVectorCountEnds[i] = s->sEnd;
					}
				
				j = p->numSegNode;
				r->numSegNode = j;
				a = b = 0;
				
				for (i = 0; i < q->numSegNode; i++)	/* for each segment of q node going to r node */
					{
					if (post(i,q->index,distance) > numTotalSegments) /* Cheking */
						{
						fprintf (fpmpi, "\n post = %d > numTotalSegments = %d", post(w,q->index,distance), numTotalSegments);
						exit (-7);
						}
					n = segments + post(i,q->index,distance);
				
					for (w = 0; w < j; w++)
						if (n->sStart == coalVectorCountStarts[w] && n->sEnd == coalVectorCountEnds[w]) /* Segmento repetido */
							a++; 
					
					if (a == 0)
						{
						m = segments + post(b+j,r->index,distance);					/* new ancester */						
						
						r->numSegNode++;
						n->before1 = m;
						m->before1 = NULL;
						m->before2 = NULL;
						m->after1 = n;
						m->after2 = NULL;
						m->sIndexNode = newInd;
						m->sIndex = actSegIndex;
						m->sStart = n->sStart;
						m->sEnd = n->sEnd;				
					
						actSegIndex++;
						b++;
						actNumSegments++;
						/*r->numSegNode++;*/
						if (m->sStart == 0 && m->sEnd == 0) /* unreal segments */
							{
							fprintf (fpmpi, "\nNot to be here. COAL2, segment start and end = 0");
							
							m->before1 = m->before2 = m->after1 = m->after2 = NULL;
							actNumSegments--;
							}
						}
					a = 0;
					}
				
				a = b = 0;
				free (coalVectorCountStarts);
				free (coalVectorCountEnds);
				

					/* Segment Bonds when this segment goes to 2 descendants segments */
				coalEqualSegInit_p = (int *) calloc((p->numSegNode+q->numSegNode),(long) sizeof(int));
				if (!coalEqualSegInit_p)
					{
					fprintf (fpmpi, "Could not allocate coalEqualSegInit_p (%lu bytes)\n", (p->numSegNode+q->numSegNode) *(long) sizeof(int));
					exit (-1);
					}
				coalEqualSegEnd_p = (int *) calloc((p->numSegNode+q->numSegNode),(long) sizeof(int));
				if (!coalEqualSegEnd_p)
					{
					fprintf (fpmpi, "Could not allocate coalEqualSegEnd_p (%lu bytes)\n", (p->numSegNode+q->numSegNode) *(long) sizeof(int));
					exit (-1);
					}
				coalEqualSegInit_q = (int *) calloc((p->numSegNode+q->numSegNode),(long) sizeof(int));
				if (!coalEqualSegInit_q)
					{
					fprintf (fpmpi, "Could not allocate coalEqualSegInit_q (%lu bytes)\n", (p->numSegNode+q->numSegNode) *(long) sizeof(int));
					exit (-1);
					}
				coalEqualSegEnd_q = (int *) calloc((p->numSegNode+q->numSegNode),(long) sizeof(int));
				if (!coalEqualSegEnd_q)
					{
					fprintf (fpmpi, "Could not allocate coalEqualSegEnd_q (%lu bytes)\n", (p->numSegNode+q->numSegNode) *(long) sizeof(int));
					exit (-1);
					}
				
				for (w = 0; w < p->numSegNode; w++)	
					{
					s = segments + post(w,p->index,distance);
					coalEqualSegInit_p[w] = s->sStart;
					coalEqualSegEnd_p[w] = s->sEnd;
					}
				for (w = 0; w < q->numSegNode; w++)	
					{
					n = segments + post(w,q->index,distance);
					coalEqualSegInit_q[w] = n->sStart;
					coalEqualSegEnd_q[w] = n->sEnd;
					}
				for (w = 0; w < p->numSegNode; w++)
					{
					for (i = 0; i < q->numSegNode; i++)
						{
						if (coalEqualSegInit_p[w] != 0 && coalEqualSegEnd_p[w] != 0 && coalEqualSegInit_p[w] == coalEqualSegInit_q[i] && coalEqualSegEnd_p[w] == coalEqualSegEnd_q[i])
							{
							for (a = 0; a < r->numSegNode; a++)	
								{
								m = segments + post(a,r->index,distance);
									
								if (m->sStart == coalEqualSegInit_p[w] && m->sEnd == coalEqualSegEnd_p[w])
									{
									if (noisy == 4)
										fprintf (fpmpi, "\nFragment %d links to 2 descendants, ", m->sIndex);
										
									for (b = 0; b < p->numSegNode; b++)
										{
										s = segments + post(b,p->index,distance);
												
										if (coalEqualSegInit_p[w] == s->sStart && coalEqualSegEnd_p[w] == s->sEnd)
											{
											m->after1 = s;
											if (noisy == 4)
												fprintf (fpmpi, "fragment %d", s->sIndex);
											}
										}
									for (b = 0; b < q->numSegNode; b++)
										{
										n = segments + post(b,q->index,distance);
										
										if (coalEqualSegInit_p[w] == n->sStart && coalEqualSegEnd_p[w] == n->sEnd)
											{
											m->after2 = n;
											if (noisy == 4)
												fprintf (fpmpi, " and fragment %d", n->sIndex);
											}
										}
									}
								}
							}
						}
					}
						
				
				for (i = 0; i < r->numSegNode; i++)			/* cheking, if there are 2 similar segments, it keeps only 1 */		
					{
					m = segments + post(i,r->index,distance);
				
					for (w = 0; w < r->numSegNode; w++)
						{
						z = segments + post(w,r->index,distance);
					
						if (w != i && m->sStart == z->sStart && m->sEnd == z->sEnd && m->after1 != NULL) /* Cheking */
							{
							fprintf (fpmpi, "\n1COAL. Not to be here. it does not to be 2 equal segments in this node %d. index %d y %d", m->sIndex, z->sIndex, r->index);
							fprintf (fpmpi, "\n m->sStart = %d, m->sEnd = %d, m->sIndexNode = %d", m->sStart, m->sEnd, m->sIndexNode);
							fprintf (fpmpi, "\n z->sStart = %d, z->sEnd = %d, z->sIndexNode = %d", z->sStart, z->sEnd, z->sIndexNode);
							
							fprintf (fpmpi, "\n post(i,r->index,distance) = %d", post(i,r->index,distance));
							fprintf (fpmpi, "\n post(w,r->index,distance) = %d", post(w,r->index,distance));
							fprintf (fpmpi, "\n r->index = %d", r->index);
								
							z->after1 = z->after2 = z->before1 = z->before2 = NULL; /* continue only the z segment */
							r->numSegNode--;
							exit(-1);
							}
						}
					}
				j = w = a = b = 0;
				free (coalEqualSegInit_p);
				free (coalEqualSegEnd_p);
				free (coalEqualSegInit_q);
				free (coalEqualSegEnd_q);
				
				
				if (noisy == 4)
					{
					fprintf (fpmpi, "\n\nCoalescence Result is the new node %d with %d fragment(s)", r->index, r->numSegNode);
					for (i = 0; i < r->numSegNode; i++)	
						{
						m = segments + post(i,r->index,distance);
						/*fprintf (fpmpi, "\npost(i,r->index,distance) = %d",post(i,r->index,distance));*/
						/*fprintf (fpmpi, "\nNode r->index = %d", r->index);*/
						fprintf (fpmpi, "\nFragment %d, ", m->sIndex);	
						fprintf (fpmpi, "m->sStart = %d", m->sStart);
						fprintf (fpmpi, " and m->sEnd = %d", m->sEnd);
						}
					fprintf (fpmpi, "\n");
					}
					
			
				/* MRCA */
				/* Array with the information about ancestral stuff */
				/* S_MRCA. overFirst is the biggest site by left & overEnd is the smallest site by right. The difference is the overLapSites */
				sizeNode_p = p->numSegNode;
				sizeNode_q = q->numSegNode;
				out = 0;
			
				if (r->numSegNode == 1) /* 1.- fast. Only 1 segment */
					{
					for (i = 0; i < r->numSegNode; i++)	
						m = segments + post(i,r->index,distance);
						
					overFirst = m->sStart;
					overEnd = m->sEnd;
				
					for (j = 1; j <= numNuc; j++)
						if (j >= overFirst && j <= overEnd)
							{
							S_MRCA[j]--;
							}
					out = 1;
					}
				a = b = k = 0;
			
				if (r->numSegNode > 1)	
					{
					initialVector_pp = (int *) calloc((p->numSegNode),(long) sizeof(int));
					if (!initialVector_pp)
						{
						fprintf (fpmpi, "Could not allocate initialVector_pp (%lu bytes)\n", (p->numSegNode) *(long) sizeof(int));
						exit (-1);
						}
					endVector_pp = (int *) calloc((p->numSegNode),(long) sizeof(int));
					if (!endVector_pp)
						{
						fprintf (fpmpi, "Could not allocate endVector_pp (%lu bytes)\n", (p->numSegNode) *(long) sizeof(int));
						exit (-1);
						}
					initialVector_qq = (int *) calloc((q->numSegNode),(long) sizeof(int));
					if (!initialVector_qq)
						{
						fprintf (fpmpi, "Could not allocate initialVector_qq (%lu bytes)\n", (q->numSegNode) *(long) sizeof(int));
						exit (-1);
						}
					endVector_qq = (int *) calloc((q->numSegNode),(long) sizeof(int));
					if (!endVector_qq)
						{
						fprintf (fpmpi, "Could not allocate endVector_qq (%lu bytes)\n", (q->numSegNode) *(long) sizeof(int));
						exit (-1);
						}
				
					/* 2.- Fast. There is a segment more big than the other segments */
					for (position = 0; position < p->numSegNode; position++)
						{
						s = segments + post(position,p->index,distance);
						initialVector_pp[position] = s->sStart;
						endVector_pp[position] = s->sEnd;
						}
					for (position = 0; position < q->numSegNode; position++)
						{
						n = segments + post(position,q->index,distance);
						initialVector_qq[position] = n->sStart;
						endVector_qq[position] = n->sEnd;
						}
					
					/* smaller value for the initial */
					for (position = 0; position < p->numSegNode; position++)
						{
						if (position == 0)
							minInit_pp = initialVector_pp[position];
						if (initialVector_pp[position] < minInit_pp && initialVector_pp[position] != 0)
							minInit_pp = initialVector_pp[position];
						}

					/* older value for the end */
					for (position = 0; position < p->numSegNode; position++)
						{
						if (position == 0)
							maxEnd_pp = endVector_pp[position];
						if (endVector_pp[position] > maxEnd_pp && endVector_pp[position] != 0)
							maxEnd_pp = endVector_pp[position];
						}
		
					/* smaller value for the initial */
					for (position = 0; position < q->numSegNode; position++)
						{
						if (position == 0)
							minInit_qq = initialVector_qq[position];
				
						if (initialVector_qq[position] < minInit_qq && initialVector_qq[position] != 0)
							minInit_qq = initialVector_qq[position];
						}

					/* older value for the end */
					for (position = 0; position < q->numSegNode; position++)
						{
						if (position == 0)
							maxEnd_qq = endVector_qq[position];
	
						if (endVector_qq[position] > maxEnd_qq && endVector_qq[position] != 0)
							maxEnd_qq = endVector_qq[position];
						}
			
					/* fast, there are in 1 big segment that contains all the segments */					
					for (position = 0; position < p->numSegNode; position++) 
						{
						if (endVector_pp[position] == maxEnd_pp && initialVector_pp[position] == minInit_pp)	
							a = 1;
						}
					for (position = 0; position < q->numSegNode; position++) 
						{
						if (endVector_qq[position] == maxEnd_qq && initialVector_qq[position] == minInit_qq)	
							b = 1;
						}
				
					free (initialVector_pp);
					free (endVector_pp);
					free (initialVector_qq);
					free (endVector_qq);
				
				
					/* 2.- Fast with a big segment */
					if (minInit_pp <= minInit_qq && maxEnd_pp >= maxEnd_qq && a == 1 && b == 1)
						{
						for (j = 1; j <= numNuc; j++)
							{
							if (j >= minInit_qq && j <= maxEnd_qq)
								{
								S_MRCA[j]--;
								}
							}
						a = out = 1;	
						}
					if (minInit_qq <= minInit_pp && maxEnd_qq >= maxEnd_pp && a == 1 && b == 1 && out == 0)
						{
						for (j = 1; j <= numNuc; j++)
							{
							if (j >= minInit_pp && j <= maxEnd_pp)
								{
								S_MRCA[j]--;
								}
							}
						a = out = 1;	
						}
					
					/* 3.- Fast, noncoincident nodes. Ex: Coalescence from 2 nodes which come from of the same recombination. MRCA variation = 0 */
					if (maxEnd_pp < minInit_qq && out == 0)
						{
						/*fprintf (fpmpi, "happened for MRCA, MRCA does not change");*/;
						
						a = out = 1;
						}
					if (maxEnd_qq < minInit_pp && out == 0)
						{
						/*fprintf (fpmpi, "happened for MRCA, MRCA does not change");*/;
						
						a = out = 1;
						}
					
					if (minInit_pp < minInit_qq && maxEnd_pp < maxEnd_qq && a == 1 && b == 1 && out == 0)
						{
						for (j = 1; j <= numNuc; j++)
							{
							if (j >= minInit_qq && j <= maxEnd_pp)
								{
								S_MRCA[j]--;
								}
							}
						a = out = 1;
						}
						
					if (minInit_qq < minInit_pp && maxEnd_qq < maxEnd_pp && a == 1 && b == 1 && out == 0)
						{
						for (j = 1; j <= numNuc; j++)
							{
							if (j >= minInit_pp && j <= maxEnd_qq)
								{
								S_MRCA[j]--;
								}
							}
						a = out = 1;
						}
					
					a = b = k = w = i = 0;
				
					/* 3.- complex case */
					if (out == 0)
						{
						for (j = 1; j <= numNuc; j++)
							{
							if (overLapSegmentsCoalMRCA(p, q, sizeNode_p, sizeNode_q, j) == YES)
								{
								S_MRCA[j]--;
								}
							}
						}
					}
				j = 0;
		


				/* Folllowing only ancestral material */
				for (mmm = 1; mmm <= numNuc; mmm++)	 /* Only anc material*/
					{
					if (p->SitesNonAncHere[mmm] == 0 && q->SitesNonAncHere[mmm] == 0)
						{
						/*fprintf (fpmpi, "\n This Position %d is ANC mat: p->index = %d, p->SitesNonAncHere = %d, q->index = %d, q->SitesNonAncHere = %d \n", mmm, p->index, p->SitesNonAncHere[mmm], q->index, q->SitesNonAncHere[mmm]);*/
						OnlyAncS_MRCA[mmm]--;
						}				
					}
				/* Is this the GMRCA of the anc material? */
				Ok_SMRCA_Codon = 0;
				for (mmm = 1; mmm <= numNuc; mmm++)	 /* Only anc material*/
					{
					if (OnlyAncS_MRCA[mmm] == 1)
						{
						/*fprintf (fpmpi, "\n-- Yes anc GMRCA position %d -- \n", mmm);*/
						Ok_SMRCA_Codon++;
						}
					}			
				if (Ok_SMRCA_Codon == numNuc && AncGMRCA_obtained == NO)
					{
					r->GMRCA_ancestral = YES;
					AncGMRCA_obtained = YES;
					if (noisy == 4)
						fprintf (fpmpi, "\n--GMRCA of the ancestral material in node %d-- \n", r->index);
					}




				/* in coalescence is possible link the nodes */
				r->left = p;
				r->right = q;
				p->anc1 = r;
				q->anc1 = r;
				r->time = currentTime;
				/*fprintf (fpmpi, " r->index = %d, r->time = %lf\n", r->index, r->time);*/
													
				/* readjust active nodes */
				activeGametes[firstInd] = newInd;
				activeGametes[secondInd] = activeGametes[numActiveGametes-1];
				numActiveGametes--; /* it lose 1 active node */
				nextAvailable++; /* 1 node more to available */
				numParcialActiveGametes[whichDeme] = numParcialActiveGametes[whichDeme]-1;
			
			
				if (nextAvailable >= numNodes)	/* if there aren't enough nodes it go into and it addition more */
					{
					/* ReallocNodes(&numNodes, activeGametes); */
					numNodes += INCREMENT_NODES;
					numTotalSegments += (INCREMENT_NODES*maxSegNode)+numNuc;
					
					/* REALLOC */
					segments = (TreeSegment *) realloc (segments, numTotalSegments  * (long) sizeof(TreeSegment));
					if (!segments)
						{
						fprintf (fpmpi, "Could not reallocate segments (%lu bytes)\n", ((numNodes*distance)+numNuc)  * (long) sizeof(TreeSegment));
						exit (1);
						}
					nodes = (TreeNode *) realloc (nodes, numNodes  * (long) sizeof(TreeNode));
					if (!nodes)
						{
						fprintf (fpmpi, "Could not reallocate nodes (%lu bytes)\n", numNodes  * (long) sizeof(TreeNode));
						exit (-1);
						}
					activeGametes = (int *) realloc (activeGametes, numNodes *(long) sizeof(int));
					if (!activeGametes)
						{
						fprintf (fpmpi, "Could not reallocate activeGametes (%lu bytes)\n", numNodes *(long) sizeof(int));
						exit (-1);
						}
					if (noisy == 4)
						fprintf (fpmpi, "\n\n...Doing reallocation of nodes (1)\n");
					}
				}	/* end of coalescence */
				

			/*** Special event of CONVERGENCIE OF DEMES ***/
			if (doConvNext == YES)
				{
				j = 0;
				doConvNext = NO;		

				currentBigDeme++;
				currentDemesNumber--;

				if (noisy > 1)
					{
					fprintf (fpmpi, "Convergence of demes %d and %d", deme_a[nextConvNumber], deme_b[nextConvNumber]);
					fprintf (fpmpi, " to deme %d", currentBigDeme+numPopulations);
					}

				CurrentDemesState[deme_a[nextConvNumber]] = CurrentDemesState[deme_b[nextConvNumber]] = 0;
				CurrentDemesState[currentBigDeme+numPopulations] = currentBigDeme+numPopulations;
									
				for (d = 0; d < numActiveGametes; d++)
					{
					p = nodes + activeGametes[d];
					
					if (p->indexCurrentMigPop == deme_a[nextConvNumber])
						{
						p->indexCurrentMigPop = currentBigDeme+numPopulations;
						numParcialActiveGametes[currentBigDeme+numPopulations]++;
						numParcialActiveGametes[deme_a[nextConvNumber]]--;
						}
						
					if (p->indexCurrentMigPop == deme_b[nextConvNumber])
						{
						p->indexCurrentMigPop = currentBigDeme+numPopulations;
						numParcialActiveGametes[currentBigDeme+numPopulations]++;
						numParcialActiveGametes[deme_b[nextConvNumber]]--;
						}
					}

				/* when the tip node has a time higher than convergence demes, its initial deme must be the deme of the convergence demes */
				/*if (doDatedTips == YES) 
					{
					for (ss = 0; ss < numSequences; ss++)		
						{
						p = nodes + ss;
						if (p->time >= currentTime)
							{
							if (p->indexOldMigPop == deme_a[nextConvNumber] || p->indexOldMigPop == deme_b[nextConvNumber])
								{
								if (noisy > 2)
									fprintf (fpmpi, "\nInitial node %d that belongs to the deme %d, is belonging now to the deme ", p->index, p->indexOldMigPop);
								p->indexOldMigPop = currentBigDeme+numPopulations;
								p->indexCurrentMigPop = currentBigDeme+numPopulations;
								if (noisy > 2)
									fprintf (fpmpi, "%d by a convergence of demes", p->indexOldMigPop);
								}
							}
						}
					if (noisy > 2)
						fprintf (fpmpi, "\n");
					}*/


				if (noisy > 3)
					fprintf (fpmpi,"\n");
				}	


			/* print out ancestral (active) status for each site and MRCA vector - only mat anc */
			if (noisy > 2)
				{
				fprintf (fpmpi,"\n - Ancestral MRCA in the nodes: -\n");
				for (i=0; i<numActiveGametes; i++)
					{
					for (j=1; j<=numNuc; j++)
						{
						p = nodes + activeGametes[i];
						if (j == 1)
							fprintf (fpmpi, "%4d -- (MRCA:)", p->index);
						fprintf (fpmpi, "%d", OnlyAncS_MRCA[j]);
						}
					fprintf (fpmpi, "\n");
					}
				/*fprintf (fpmpi, "MRCA   ");*/
				for (j=1; j<=numNuc; j++)
					{
					if (OnlyAncS_MRCA[j] <= 1)
						fprintf (fpmpi, "*");
					else
						fprintf (fpmpi, " ");
					}
				/*fprintf (fpmpi, "\n");*/
				fprintf (fpmpi,"\n\n");
				}



			/* print out ancestral (active) status for each site and MRCA vector */
			if (noisy > 2)
				{
				fprintf (fpmpi,"\nMRCA in the nodes:\n");
				for (i=0; i<numActiveGametes; i++)
					{
					for (j=1; j<=numNuc; j++)
						{
						p = nodes + activeGametes[i];
						/*p = nodes + pos(activeGametes[i],j,numNuc);*/
						if (j == 1)
							fprintf (fpmpi, "%4d -- (MRCA:)", p->index);
						fprintf (fpmpi, "%d", S_MRCA[j]);
						}
					fprintf (fpmpi, "\n");
					}
				/*fprintf (fpmpi, "MRCA   ");*/
				for (j=1; j<=numNuc; j++)
					{
					if (S_MRCA[j] /* MRCA[j]*/ <= 1)
						fprintf (fpmpi, "*");
					else
						fprintf (fpmpi, " ");
					}
				/*fprintf (fpmpi, "\n");*/
				}
				
			sizeNode = sizeNode_p = sizeNode_q = 0;
			free (gi);



			/* If doing dated tips, sometimes the number of actives nodes can be already 1 before all have been activated samples. 
			If this happens we need to move towards the next sample, activate it, and start again */
			if (doDatedTips == YES)
			 if (numActiveGametes == 1 && currentSample > 0)
				{
				if (noisy > 2)
					fprintf (stderr, "\n\nOnly 1 lineage active before activating all samples"/*, currentSample*/);
				currentSample--;

				/* activate nodes from this sample*/
				if (noisy > 2)
					fprintf (stderr, "\n Activating sample %d (time = %6.4f). Tips to activate:", currentSample, datedSample[currentSample].time);
				for (i=0; i<datedSample[currentSample].size; i++)
					{
					if (noisy > 2)
						fprintf (stderr, " %d", datedSample[currentSample].member[i]-1);
					
					p = nodes + datedSample[currentSample].member[i]-1;
					p->index = datedSample[currentSample].member[i]-1;
					activeGametes[numActiveGametes] = datedSample[currentSample].member[i]-1; 

					for (ss = 1; ss <= numPopulations /*+ numCONV*/; ss++) 
						{
						if (p->indexOldMigPop == ss)
							{
							numParcialActiveGametes[ss]++;
							/*fprintf (stderr, "\n AQUI numParcialActiveGametes[%d] = %d \n", ss, numParcialActiveGametes[ss]); */
							}
						} 
						
					numActiveGametes++;
					}
				currentTime = datedSample[currentSample].time;
	
				if (noisy > 2)
					{
					fprintf (stderr, "\nActive nodes (%d):", numActiveGametes); 
					for (i=0; i<numActiveGametes; i++)
						fprintf (stderr," %d",activeGametes[i]);
					fprintf (stderr,"   Next node available = %d ", nextAvailable);
					fprintf (stderr, "\nSetting currentTime = %6.4f", currentTime);
					}
				}
			} /***** coalescent tree finished *****/


		for (w = 1; w <= numNuc; w++)
			{
			if (S_MRCA[w] > 1 || S_MRCA[w] < 1)
				{
				fprintf (fpmpi, "\n Warning S_MRCA in the last node is < > 1, S_MRCA[%d] = %d", w, S_MRCA[w]);
				exit (-1);
				}
			}

		if (noisy > 1)
			fprintf (fpmpi, "\n\n\n>> Coalescent tree/s finished\n");
			
			
		counterTimeInit = counterTimeInit + currentTime;
		actualTGMRCA = currentTime;
		
			/* free memory of migrations */			
		free (rateREpartial);
		free (rateCApartial);
		free (rateMIGpartial);
		free (ratePartial);
		free (cumPopulTase);
		free (cumInitPopul);
		free (numParcialActiveGametes);
		free (GiPartial);
		if (doConvergDemes == YES)
			{
			free (currentConvDem);
			free (convDemTimes);
			free (deme_a);
			free (deme_b);
			free (CurrentDemesState);
			}
		/* Use to see the evolution nodes-segments */
		/*#ifdef MPI
			{
			fprintf (fpmpi, "\n\n\nnextAvailable = %d",nextAvailable);
			for (j = 0; j < numNodes; j++)*/	/* looking for parent Node of segments */	
			/*	{
				p = nodes + j;
				fprintf (fpmpi, "\n\n\n\n**The node %d**",p->index);

				for (i=0; i < p->numSegNode; i++)
					{
					s = segments + post(i,j,numNuc);	
					p->seg = s;
			
					fprintf (fpmpi, "\n\ns->sIndex %d con s->parentNode->index = %d", s->sIndex, s->parentNode->index);
					if (s->after1 != NULL)
						fprintf (fpmpi, "\ns->after1->sIndex %d con s->after1->parentNode->index = %d", s->after1->sIndex, s->after1->parentNode->index);
					if (s->after2 != NULL)
						fprintf (fpmpi, "\ns->after2->sIndex %d con s->after2->parentNode->index = %d", s->after2->sIndex, s->after2->parentNode->index);
					}
				}
			}
		*/
		}
	
	if (doMigration == NO) /* Coalescence without migration */
		{		
		eventNum = 0;		
		currentTime = 0.0;
		period = 1;

		while (numActiveGametes > 1)
			{
			/*fprintf (fpmpi,"\nMMM >>>>> NEW EVENT. numActiveGametes = %d", numActiveGametes);*/

			Gi = 0;
			/* allocate memory for each node gi */
			gi = (int *) calloc(numActiveGametes,(long) sizeof(int));
			if (!gi)
				{
				fprintf (fpmpi, "Could not allocate gi (%lu bytes)\n", numActiveGametes *(long) sizeof(int));
				exit (-1);
				}
			/* calculate gi for each node and total Gi */	
			/* Gi is the total number of ancestral sites, gi is a vector with ancestral and not found MRCA sites */
			for (i = 0; i < numActiveGametes; i++)
				{
				p = nodes + activeGametes[i];
				sizeNode = p->numSegNode;
				gi[i] = CalcIndividualGi (i, nodes, activeGametes, numNuc, S_MRCA, sizeNode);
				Gi += gi[i];
				if (noisy == 4)
					fprintf (fpmpi,"\n%d \n", gi[i]);
				}	
			if (noisy == 4)
				fprintf (fpmpi," Gi = %lu ", Gi);
			

			/*fprintf (fpmpi,"\nMMM Gi = %lu, Nscaling = %d, N = %d, recombinationRate = %lf ", Gi, Nscaling, N, recombinationRate);*/	
			
	
			/* get rates for events */  
			rateRE = 1.0 * Gi * Nscaling * N * recombinationRate;        /* recombinationRate is constant */
			rateCA = numActiveGametes * (numActiveGametes - 1) / 2.0;
			rate = rateCA + rateRE;
			
			/*fprintf (fpmpi,"\nMMM numActiveGametes = %d", numActiveGametes);
			fprintf (fpmpi,"\nMMM rateRE = %3.2f, rateCA = %3.2f, rate = %3.2f \n", rateRE, rateCA, rate);*/


			/* find out time for coalescence */
			if (doDemographics == YES)
				{
				periodGrowth[period] = -log (Nend[period] / (double) Nbegin[period]) / (double) (cumDuration[period] - cumDuration[period-1]);
				/*fprintf (fpmpi, "\n>>>>1 period growth  = %f (period = %d)", periodGrowth[period], period);*/
				if (isnan(periodGrowth[period]) == YES)
					{
					fprintf (fpmpi, "\nERROR: period growth (%f) is NaN", periodGrowth[period]);
					fprintf (fpmpi, "\n      This might suggest that the growth rate is too negative");
					fprintf (fpmpi, "\n      and the coalescent time is therefore infinite.");
					fprintf (fpmpi, "\n      Try a smaller value");
					exit (1);
					}

				if (Nend[period] == Nbegin[period])
					{
					timeCA = RandomExponential (rateCA, seed) * Nscaling * (double) Nbegin[period];
					}
				else
					{
					timeCA = log (1 + RandomExponential (rateCA, seed) * periodGrowth[period] * Nscaling * Nbegin[period] * 
							exp (-periodGrowth[period] * (currentTime - cumDuration[period-1]))) / periodGrowth[period];
					}

				/*	When growth rate is very negative, coalescent time may be infinite
					this results in log (-x) => timCA = NaN. If this not the last period
					just jump to the next. If this is the last period, we have to exit
					the program */
				if (isnan(timeCA) == YES)
					{
					if (period < numPeriods) 
						{
						currentTime = cumDuration[period];
						period++;
						periodGrowth[period] = -log (Nend[period] / (double) Nbegin[period]) / (double) (cumDuration[period] - cumDuration[period-1]);
						continue;
						}	
					else
						{				
						fprintf (fpmpi, "\nERROR: Coalescent time (%f) is infinite ", timeCA);
						fprintf (fpmpi, "\n      This might suggest that the growth rate is too negative");
						fprintf (fpmpi, "\n      and the coalescent time is therefore infinite.");
						fprintf (fpmpi, "\n      Try a smaller value");
						exit (1);
						}
					}
				}
			else
				{
				timeCA = RandomExponential (rateCA, seed) * Nscaling * N;
				if (doExponential == YES)
					{
					timeCA = log (exp(growthRate*currentTime) + growthRate * timeCA) / growthRate - currentTime;
		
					/*	When growth rate is very negative, coalescent time may be infinite
						this results in log (-x) => timeCA = NaN. We have to exit
						the program */
					if (isnan(timeCA) == YES)
						{
						fprintf (fpmpi, "\nERROR: Coalescent time (%f) is infinite ", timeCA);
						fprintf (fpmpi, "\n      This might suggest that the growth rate is too negative");
						fprintf (fpmpi, "\n      and the coalescent time is therefore infinite.");
						fprintf (fpmpi, "\n      Try a smaller value");			
						exit (1);
						}
					}
				}

			/* find out time for recombination */
			timeRE = RandomExponential (rateRE, seed) * Nscaling * N;


			if (doDatedTips == YES)
				{
				if (timeCA < timeRE)
					{
					eventTime = timeCA;
					isCoalescence = YES;	
					}
				else
					{
					eventTime = timeRE;
					isCoalescence = NO;	
					}


				/* if doing dated tips, check whether we need to activate a new sample, update sampling period and start again */
				if ((currentTime + eventTime) > datedSample[currentSample-1].time && currentSample > 0)
					{
					currentSample--;
					/* activate nodes from this sample */
					if (noisy > 2)
						fprintf (stderr, "\nCumulative time = %6.4f  > sample %d time = %6.4f. Activating tips:", currentTime + eventTime, currentSample, datedSample[currentSample].time);
					
					for (i=0; i<datedSample[currentSample].size; i++)
						{
						if (noisy > 2)
							fprintf (stderr, " %d", datedSample[currentSample].member[i]-1);

						p = nodes + datedSample[currentSample].member[i]-1;
						p->index = datedSample[currentSample].member[i]-1;
						activeGametes[numActiveGametes] = datedSample[currentSample].member[i]-1;  

						numActiveGametes++;
						}
					currentTime = datedSample[currentSample].time;

					if (noisy > 2)
						{
						fprintf (stderr, "\nActive nodes (%d):", numActiveGametes); 
						for (i=0; i<numActiveGametes; i++)
							fprintf (stderr," %d",activeGametes[i]);
						fprintf (stderr,"\nNext node available = %d", nextAvailable);
						fprintf (stderr, "\nSample %d activated and going back to currentTime = %6.4f", currentSample, currentTime);
						}
					continue; /* start again*/
					}	


				/* event is a coalescence or a recombination? */
				if (isCoalescence == YES)
					{
					/*	if this period is not the last one and if the event time is outside the current interval,
						update period and start again */
					if (doDemographics == YES && period < numPeriods && (currentTime + eventTime) > cumDuration[period])
						{
						currentTime = cumDuration[period];
						period++;
						periodGrowth[period] = -log (Nend[period] / (double) Nbegin[period]) / (double) (cumDuration[period] - cumDuration[period-1]);
						continue;
						}	
			
					numCA++;
					}
				else
					{
					/*	if this period is not the last one and if the event time is outside the current interval,
						update period and start again */
					if (doDemographics == YES && period < numPeriods && (currentTime + eventTime) > cumDuration[period])
						{
						currentTime = cumDuration[period];
						period++;
						periodGrowth[period] = -log (Nend[period] / (double) Nbegin[period]) / (double) (cumDuration[period] - cumDuration[period-1]);
						continue;
						}	
					numRE++;
					/* reallocate for recombination breakpoints */
					if (numRE >= (memoryBreakp-1))
						{
						memoryBreakp += 50;
					 
						breakpoint = (int *) realloc(breakpoint, memoryBreakp *(long) sizeof(int));
						if (!breakpoint)
							{
							fprintf (fpmpi, "Could not reallocate breakpoint \n");
							exit (-1);
							}
						if (noisy == 4)
							fprintf (fpmpi, "\n...Doing reallocation of breakponts (1)\n");
						}
					if (numRE > (numNuc+1) && many == 0)
						{
						if (noisy > 2)
							fprintf (fpmpi, "\n\n Cheking information: Many recombinations %d (more recombinations that sites)\n", numRE);
						
						many++;
						}
					}
				}
			else /* not tip dates */
				{
				/* event is a coalescence or a recombination? */
				if (timeCA < timeRE)
					{
					isCoalescence = YES;
					eventTime = timeCA;
					/*	if this period is not the last one and if the event time is outside the current interval,
						update period and start again */
					if (doDemographics == YES && period < numPeriods && (currentTime + eventTime) > cumDuration[period])
						{
						currentTime = cumDuration[period];
						period++;
						periodGrowth[period] = -log (Nend[period] / (double) Nbegin[period]) / (double) (cumDuration[period] - cumDuration[period-1]);
						continue;
						}	
			
					numCA++;


					}
				else
					{
					isCoalescence = NO;
					eventTime = timeRE;	
					/*	if this period is not the last one and if the event time is outside the current interval,
						update period and start again */
					if (doDemographics == YES && period < numPeriods && (currentTime + eventTime) > cumDuration[period])
						{
						currentTime = cumDuration[period];
						period++;
						periodGrowth[period] = -log (Nend[period] / (double) Nbegin[period]) / (double) (cumDuration[period] - cumDuration[period-1]);
						continue;
						}	
					numRE++;


					/* reallocate for recombination breakpoints */
					if (numRE >= (memoryBreakp-1))
						{
						memoryBreakp += 50;
					 
						breakpoint = (int *) realloc(breakpoint, memoryBreakp *(long) sizeof(int));
						if (!breakpoint)
							{
							fprintf (fpmpi, "Could not reallocate breakpoint \n");
							exit (-1);
							}
						if (noisy == 4)
							fprintf (fpmpi, "\n...Doing reallocation of breakponts (1)\n");
						}
					if (numRE > (numNuc+1) && many == 0)
						{
						if (noisy > 2)
							fprintf (fpmpi, "\n\n Cheking information: Many recombinations %d (more recombinations that sites)\n", numRE);
						
						many++;
						}
					}

				}





				
			/* set time */
			currentTime += eventTime; /* the time is accumulated */

				
			if (noisy > 3)
				fprintf(fpmpi, "\n\n");
				
			eventNum++;
			if (noisy > 1)
				fprintf (fpmpi, "\n\n*** Event %3d *** rate = %lf, currentTime = %lf\n", eventNum, rate, currentTime);
				
				
		
			/*** if RECOMBINATION, readjust active nodes ***/
			if (isCoalescence == NO)
				{
				if (noisy == 4)
					fprintf (fpmpi, "\n* Recombination *");
					
				/* Which node has the recombination? */
					/* assign probability to each node based on their gi's values */
				cum_gi = 0;
				probRecIndividual = Gi * RandomUniform(seed); /* it calculate the individual probability for the breakpoint */
				for (whichInd=0; whichInd<numActiveGametes; whichInd++)
					{
					cum_gi += gi[whichInd];				/* accumulate gi into cum_gi. whichInd is the choose node. */
					if (probRecIndividual < cum_gi)    /* break whether the node is chosen */
						break;
					}
			
				if (whichInd >= numActiveGametes)
					{
					fprintf (fpmpi, "\n\nERROR: whichInd out of range3!: whichInd = %d\n", whichInd);
					exit (-1);
					}
			
				/* select a valid breakpoint among potential recombining locations */
				/* to be a potential recombining site, a site has to have ancestral material non-MRCA before and after it */
				legalBreakpoint = NO;
				while (legalBreakpoint == NO)
					{
					do /* oct2009 */
						{
						whichSite = (/*numSites*/numNuc) * RandomUniform(seed); /* it choose a site */
						if (whichSite == 1)
							whichSite = numNuc;
						} while (whichSite == 0);
					/* oct2009 */

					if (whichSite > numNuc)
						{
						fprintf (fpmpi, "\n\nERROR: whichSite out of range! : whichSite = %d\n", whichSite);
						exit (-1);
						}
					p = nodes + activeGametes[whichInd];
					sizeNode = p->numSegNode;
					if	(IsValidBreakSite (activeGametes, nodes, whichInd, whichSite, S_MRCA) == YES)
						legalBreakpoint = YES;
					}
				
				/* should this recombination event be counted in the expected number of recombinations E(R)? */
				/* for E(R) count only events with breakpoints as 1|1, 1|0 or 0|1  (i.e., not 0|0)         */
				/* if 1 represent a site that did found already its MRCA count it as a 0 */
				ThisBreakpIsTrapped = NO;
				if (doCountsForExpNumRec == YES)
					{
					p = nodes + activeGametes[whichInd];
					sizeNode = p->numSegNode;	
					if (CountsForExpNumRec (activeGametes, whichInd, whichSite, nodes, S_MRCA, sizeNode) == NO)
						{
						recNotToCount++;
						ThisBreakpIsTrapped = YES;
						/*fprintf(stderr,"\n ..not to count.. \n");*/
						}
					}
				/* copy whichIndividual to a new space in memory */
				hasPassedBreakPoint = NO;
			
				firstHalf = nextAvailable++; /* firstHalf is the first node that was created by the recombination */
				if (nextAvailable >= numNodes) /* if there aren't enough nodes it go into and it addition more */
					{
					/* ReallocNodes(&numNodes, activeGametes); */
					numNodes += INCREMENT_NODES;
					numTotalSegments += (INCREMENT_NODES*maxSegNode)+numNuc;
					
					/* REALLOC */
					segments = (TreeSegment *) realloc (segments, numTotalSegments  * (long) sizeof(TreeSegment)); 
					if (!segments)
						{
						fprintf (fpmpi, "Could not reallocate segments (%lu bytes)\n", ((numNodes*distance)+numNuc)  * (long) sizeof(TreeSegment));
						exit (1);
						}
					nodes = (TreeNode *) realloc (nodes, numNodes  * (long) sizeof(TreeNode));
					if (!nodes)
						{
						fprintf (fpmpi, "Could not reallocate nodes (%lu bytes)\n", numNodes  * (long) sizeof(TreeNode));
						exit (-1);
						}
					activeGametes = (int *) realloc (activeGametes, numNodes *(long) sizeof(int));
					if (!activeGametes)
						{
						fprintf (fpmpi, "Could not reallocate activeGametes (%lu bytes)\n",numNodes *(long) sizeof(int));
						exit (-1);
						}
					if (noisy == 4)
						fprintf (fpmpi, "\n\n...Doing reallocation of nodes (1)\n");
					}
									
				secondHalf = nextAvailable++; /* secondhalf is the second node that was created by the recombination */
				if (nextAvailable >= numNodes) /* if there aren't enough nodes it go into and it addition more */
					{
					/* ReallocNodes(&numNodes, activeGametes); */
					numNodes += INCREMENT_NODES;
					numTotalSegments += (INCREMENT_NODES*maxSegNode)+numNuc;
					
					/* REALLOC */
					segments = (TreeSegment *) realloc (segments, numTotalSegments  * (long) sizeof(TreeSegment)); 
					if (!segments)
						{
						fprintf (fpmpi, "Could not reallocate segments (%lu bytes)\n", ((numNodes*distance)+numNuc)  * (long) sizeof(TreeSegment));
						exit (1);
						}
					nodes = (TreeNode *) realloc (nodes, numNodes  * (long) sizeof(TreeNode));
					if (!nodes)
						{
						fprintf (fpmpi, "Could not reallocate nodes (%lu bytes)\n", numNodes  * (long) sizeof(TreeNode));
						exit (-1);
						}
					activeGametes = (int *) realloc (activeGametes, numNodes *(long) sizeof(int));
					if (!activeGametes)
						{
						fprintf (fpmpi, "Could not reallocate activeGametes (%lu bytes)\n", numNodes *(long) sizeof(int));
						exit (-1);
						}
					if (noisy == 4)
						fprintf (fpmpi, "\n\n...Doing reallocation of nodes (1)\n");
					}
				p = nodes + activeGametes[whichInd];
				q = nodes + firstHalf;	/* parent1 (new) */
				r = nodes + secondHalf; /* parent2 (new) */
			
				q->index = firstHalf;
				r->index = secondHalf;
				q->numSegNode = r->numSegNode = p->numSegNode;		/* Good if there are not nill segments.. then, in its case, it will be modify */
			


				q->class = 3;
				r->class = 3;
				q->GMRCA_ancestral = NO;
				r->GMRCA_ancestral = NO;
				q->breakp = whichSite;		
				r->breakp = whichSite;			
				q->time = currentTime;
				r->time = currentTime;
				q->sib = r;
				r->sib = q;
				q->left = p;
				r->left = p;
				p->anc1 = q;
				p->anc2 = r;
				
				if (doBranchNetfiles == YES)
					{
					q->NetLabelPrint = numNetLabelPrint;
					/*numNetLabelPrint++;*/
					r->NetLabelPrint = numNetLabelPrint;
					numNetLabelPrint++;
					}


				k = 0;
				for (w = 1; w < numSites; w++)
					{
					//fprintf (fpmpi, "\nstud[%d] = %d", w-1, stud[w-1]);
					if (whichSite == stud[w-1]) /* The breakpoints BETWEEN codons. "stud" is an array with the possible breakpoints beetween codons*/
						k++;

					/* exception for trapped material, any breakpoint here does not break material codons */
					if (ThisBreakpIsTrapped == YES)
						{
						k++;
						/*fprintf (fpmpi, "\n ThisBreakpIsTrapped, whichSite = %d, p->index = %d \n", whichSite, p->index);*/
						}
					}

				doBreakpBroken = NO;
				if (k == 0)
					{
					if (noisy == 4)
						{
						fprintf (fpmpi, "\n Broken codon, breakpoint at %d \n", whichSite);
						}
					//q->breakCodon = YES;
					//r->breakCodon = YES;
					variable1 = whichSite/3.00 + 0.4;
					//fprintf (fpmpi, "\n variable1 = %lf", variable1);
					q->breakCodon = fabs(variable1);
					r->breakCodon = fabs(variable1);
					numREbreakCod++;
					variable2 = fmod(whichSite,3.00);
					if (variable2 == 0)
						q->whereBreakCodon = 2;
					else
						q->whereBreakCodon = 1;
					r->whereBreakCodon = 3;
					
					

					/* oct2009 */
					/* int			doBreakpBroken, LeftLess, LeftHigh, RightLess, RightHigh; */
					if (doCodonModel == YES)
						doBreakpBroken = YES;
					LeftLess = LeftHigh = RightLess = RightHigh = LeftLess2 = RightHigh2 = -1;

					if (q->whereBreakCodon == 1) /* first codon position breakp */
						{
						RightHigh = q->breakCodon * 3;
						RightLess = RightHigh - 1;
						LeftLess = RightLess - 1; 
						LeftHigh = RightLess - 1;
						
						
						q->SitesNonAncHere[LeftLess+1] = 1;
						q->SitesNonAncHere[LeftLess+2] = 1;
						r->SitesNonAncHere[RightLess-1] = 1;

						/*fprintf (fpmpi, "\n LeftLess+2 = %d; RightLess-1 = %d \n", LeftLess+2, RightLess-1);*/
						for (mmm=1; mmm<=numNuc; mmm++)
							{
							if (mmm <= RightLess-1)
								q->SitesNonAncHere[mmm] = p->SitesNonAncHere[mmm]; 
							if (mmm > RightLess-1 && mmm <= LeftLess+2)
								q->SitesNonAncHere[mmm] = 1;
							if (mmm > LeftLess+2)
								q->SitesNonAncHere[mmm] = -1;

							if (mmm < RightLess-1)
								r->SitesNonAncHere[mmm] = -1; 
							if (mmm == RightLess-1)
								r->SitesNonAncHere[mmm] = 1;
							if (mmm > RightLess-1)
								r->SitesNonAncHere[mmm] =  p->SitesNonAncHere[mmm];
							}	
						}
					else if (q->whereBreakCodon == 2) /* second codon position breakp */
						{
						RightHigh = q->breakCodon * 3;
						RightLess = q->breakCodon * 3;
						LeftHigh = RightLess - 1;
						LeftLess = LeftHigh - 1;


						q->SitesNonAncHere[LeftHigh+1] = 1;
						r->SitesNonAncHere[RightLess-1] = 1;	
						r->SitesNonAncHere[RightLess-2] = 1;
						
						/*fprintf (fpmpi, "\n LeftHigh+1 = %d; RightLess-2 = %d \n", LeftHigh+1, RightLess-2);*/
						for (mmm=1; mmm<=numNuc; mmm++)
							{
							if (mmm < LeftHigh+1)
								q->SitesNonAncHere[mmm] = p->SitesNonAncHere[mmm]; 
							if (mmm == LeftHigh+1)
								q->SitesNonAncHere[mmm] = 1;
							if (mmm > LeftHigh+1)
								q->SitesNonAncHere[mmm] = -1;

							if (mmm >= LeftHigh+1)
								r->SitesNonAncHere[mmm] = p->SitesNonAncHere[mmm]; 
							if (mmm >= RightLess-2 && mmm < LeftHigh+1)
								r->SitesNonAncHere[mmm] = 1;
							if (mmm < RightLess-2)
								r->SitesNonAncHere[mmm] = -1;

							}	
						}
					else
						{
						fprintf (fpmpi, "error at q->whereBreakCodon intra codon Rec _ Main (%d != 1 or 2)\n", q->whereBreakCodon);
						exit (-1);
						}
					/* fprintf (fpmpi, "\n Left (Less-High) %d-%d; Right (Less-High) %d-%d \n", LeftLess, LeftHigh, RightLess, RightHigh);	*/				
					if (LeftLess == -1 || LeftHigh == -1 || RightLess == -1 || RightHigh == -1)
						{
						fprintf (fpmpi, "\n Error (value = -1): Left (Less-High) %d-%d; Right (Less-High) %d-%d \n", LeftLess, LeftHigh, RightLess, RightHigh);					
						exit (-1);
						}

					}
				if (k > 0 || doCodonModel == NO) /* inter codon rec */
					{
					doBreakpBroken = NO;

					for (mmm=1; mmm<=numNuc; mmm++)
						{
						if (mmm >= whichSite)
							{
							q->SitesNonAncHere[mmm] = -1; /* non anc mat */
							}
						if (mmm < whichSite)
							{
							q->SitesNonAncHere[mmm] = p->SitesNonAncHere[mmm]; /* non anc mat */
							}

						if (mmm < whichSite)
							{
							r->SitesNonAncHere[mmm] = -1; /* non anc mat */
							}
						if (mmm >= whichSite)
							{
							r->SitesNonAncHere[mmm] = p->SitesNonAncHere[mmm]; /* non anc mat */
							}
						}	


					}
				/* oct2009 */

				/*for (mmm=1; mmm<=numNuc; mmm++)
					{
					fprintf (fpmpi, "\n Initial. Site %d. Node: %d, p->SitesNonAncHere = %d \n", mmm, p->index, p->SitesNonAncHere[mmm]);
					}
				for (mmm=1; mmm<=numNuc; mmm++)
					{
					fprintf (fpmpi, "\n Site %d. Node: %d, q->SitesNonAncHere = %d;  Node %d, r->SitesNonAncHere = %d \n", mmm, q->index, q->SitesNonAncHere[mmm], r->index, r->SitesNonAncHere[mmm]);
					}*/
				

				/*fprintf (fpmpi, "\n r->index = %d, r->time = %lf, r->class = %d, r->breakp = %d, r->breakCodon = %d, r->whereBreakCodon = %d", r->index, r->time, r->class, r->breakp, r->breakCodon, r->whereBreakCodon);
				fprintf (fpmpi, "\n q->index = %d, q->time = %lf, q->class = %d, q->breakp = %d, q->breakCodon = %d, q->whereBreakCodon = %d\n", q->index, q->time, q->class, q->breakp, q->breakCodon, q->whereBreakCodon);*/
				k = 0;				



				if (noisy == 4)
					{		
					fprintf (fpmpi, "\nNode index %d with breakpoint on %d site", p->index, whichSite);
					fprintf (fpmpi, "\nThis node contains %d fragment(s):", p->numSegNode);
					}
				for (w = 0; w < p->numSegNode; w++)
					{
					s = segments + post(w,p->index,distance);
					/*p->seg = s;*/
					if (noisy == 4)
						{
						fprintf (fpmpi, "\ns->sIndex = %d", s->sIndex);
						fprintf (fpmpi, "\ns->sStart = %d", s->sStart);
						fprintf (fpmpi, "\ns->sEnd = %d",s->sEnd);
						}
					}
				if (noisy == 4)
					fprintf (fpmpi, "\n\n>> Process evolution..");
					
				

				
				a = b = aa = bb = aaa = bbb = out = 0;
				startsVectorRec = (int *) calloc((p->numSegNode),(long) sizeof(int));
				if (!startsVectorRec)
					{
					fprintf (fpmpi, "Could not allocate startsVectorRec (%lu bytes)\n", (p->numSegNode) *(long) sizeof(int));
					exit (-1);
					}
				endsVectorRec = (int *) calloc((p->numSegNode),(long) sizeof(int));
				if (!endsVectorRec)
					{
					fprintf (fpmpi, "Could not allocate endsVectorRec (%lu bytes)\n", (p->numSegNode) *(long) sizeof(int));
					exit (-1);
					}
				

				for (i = 0; i < p->numSegNode; i++)					/* for each segment */
					{
					s = segments + post(i,p->index,distance);
					startsVectorRec[i] = s->sStart;
					endsVectorRec[i] = s->sEnd;
				


					/* first half */
					if (s->sStart >= whichSite)
						{
						if (noisy == 4)
							fprintf (fpmpi, "\nNil segment in Left, don't make it");
							
						q->numSegNode--;
						a++;
						out = 1;
						} 
				
					if (s->sStart == 1 && s->sEnd >= (whichSite-1) && out == 0)
						{
						if (aa > 0)
							{
							if (noisy == 4)
								fprintf (fpmpi, "\nRepit segment in Left, don't make it");
								
							q->numSegNode--;
							out = 1;
							}
						aa++;
						}
					
					for (w = 0; w < i+1; w++)
						{
						if (out == 0 && s->sStart != 1 && s->sStart == startsVectorRec[w] && s->sStart < whichSite && w != i && s->sEnd >= (whichSite-1) && endsVectorRec[w] >= (whichSite-1) && startsVectorRec[w] != 0)
							{
							if (noisy == 4)
								fprintf (fpmpi, "\nRepit segment in Left, don't make it");
								
							q->numSegNode--;
							aaa++;
							out = 1;
							}
						}				
					
					
					if (out == 0 && s->sStart < whichSite)
						{
						if (aa == 0 && aaa == 0)
							{
							n = segments + post(i-a-aa-aaa,q->index,distance);
							n->sIndexNode = firstHalf;
							}
						if (aa > 0 && aaa == 0)
							{
							n = segments + post(i+1-a-aa-aaa,q->index,distance);
							n->sIndexNode = firstHalf;
							}
						if (aaa > 0 && aa > 0)
							{
							n = segments + post(i+1-a-aa-aaa,q->index,distance);
							n->sIndexNode = firstHalf;
							}
						if (aaa > 0 && aa == 0)
							{
							n = segments + post(i+0-a-aa-aaa,q->index,distance);
							n->sIndexNode = firstHalf;
							}
						nodeValue = firstHalf;


						if (doBreakpBroken == YES && s->sStart < whichSite && s->sEnd >= whichSite)
							{
							w = recSegmentsGeneratesLeftBrokenCodon(nodeValue, s, n, numNuc, whichSite, LeftLess, RightHigh, &actSegIndex); /* it makes the segments of the left node */
							LeftLess2 = LeftLess;
							RightHigh2 = RightHigh;							
							}
						else
							{
							w = recSegmentsGeneratesLeft(nodeValue, s, n, numNuc, whichSite, &actSegIndex); /* it makes the segments of the left node */
							}
						

						if (w != 1) /* Cheking */
							{
							fprintf (fpmpi, "Warning in recSegmentsGeneratesLeft");
							exit (-1);
							}
						actNumSegments++;
						if (n->sStart == 0 && n->sEnd == 0) /* Cheking, unreal segments */
							{
							fprintf (fpmpi, "\nNot to be here. Segment Left start and end = 0");
							
							n->before1 = n->before2 = n->after1 = n->after2 = NULL;
							actNumSegments--;
							exit(-1);
							}
						if (noisy == 4)
							{
							fprintf (fpmpi, "\nAfter rec. left");
							fprintf (fpmpi, "\nq->seg->sIndex = %d", n->sIndex);
							fprintf (fpmpi, "\nq->seg->sStart = %d", n->sStart);
							fprintf (fpmpi, "\nq->seg->sEnd = %d\n", n->sEnd);		
							}		

						}
					w = out = 0;
				


					/* second half */

					if (s->sEnd < whichSite)
						{
						if (noisy == 4)							
							fprintf (fpmpi, "\nNil segment in Right, don't make it");
							
						r->numSegNode--;
						b++;
						out = 1;
						}
					if (s->sStart <= whichSite && s->sEnd == numNuc && out == 0)
						{
						if (bb > 0)
							{
							if (noisy == 4)
								fprintf (fpmpi, "\nRepit segment in Right, don't make it");
								
							r->numSegNode--;
							out = 1;
							}
						bb++;
						}
					for (w = 0; w < i+1; w++)
						{
						if (out == 0 && s->sEnd != numNuc && s->sEnd == endsVectorRec[w] && s->sEnd > whichSite && w != i && s->sStart <= whichSite && startsVectorRec[w] <= whichSite && endsVectorRec[w] != 0)
							{
							if (noisy == 4)
								fprintf (fpmpi, "\nRepit segment in Right, don't make it");								
								
							r->numSegNode--;
								
							bbb++;
							out = 1;
							}
						}
					for (w = 0; w < i+1; w++)
						{
						if (out == 0 && s->sEnd == endsVectorRec[w] && s->sEnd >= whichSite && w != i && s->sStart <= whichSite && startsVectorRec[w] <= whichSite && endsVectorRec[w] != 0)
							{
							if (noisy == 4)								
								fprintf (fpmpi, "\nRepit segment in Right, don't make it");
								
							r->numSegNode--;
								
							bbb++;
							out = 1;
							}
						}
					

					if (out == 0 && s->sEnd >= whichSite)
						{
						if (bb == 0 && bbb == 0)
							{
							m = segments + post(i-b-bb-bbb,r->index,distance);
							m->sIndexNode = secondHalf;
							}
						if (bb > 0 && bbb == 0)
							{
							m = segments + post(i+1-b-bb-bbb,r->index,distance);
							m->sIndexNode = secondHalf;
							}
						if (bb > 0 && bbb > 0)
							{
							m = segments + post(i+1-b-bb-bbb,r->index,distance);
							m->sIndexNode = secondHalf;
							}
						if (bb == 0 && bbb > 0)
							{
							m = segments + post(i+0-b-bb-bbb,r->index,distance);
							m->sIndexNode = secondHalf;
							}
						/* m = s; initial, the new segments are similar at the old segments */
						nodeValue = secondHalf;


						if (doBreakpBroken == YES && s->sStart < whichSite && s->sEnd >= whichSite)
							{
							w = recSegmentsGeneratesRightBrokenCodon(nodeValue, s, m, numNuc, whichSite, LeftLess, RightHigh, &actSegIndex);
							LeftLess2 = LeftLess;
							RightHigh2 = RightHigh;
							}
						else
							{
							w = recSegmentsGeneratesRight(nodeValue, s, m, numNuc, whichSite, &actSegIndex);
							}
						
						if (w != 1)
							{
							fprintf (fpmpi, "Warning in recSegmentsGeneratesRight");
							exit (-1);
							}
						actNumSegments++;
						if (m->sStart == 0 && m->sEnd == 0) /* unreal segments */
							{
							fprintf (fpmpi, "\nNot to be here. segment Right start and end = 0");
							
							m->before1 = m->before2 = m->after1 = m->after2 = NULL;
							actNumSegments--;
							exit(-1);
							}
						if (noisy == 4)
							{
							fprintf (fpmpi, "\nAfter rec. right");
							fprintf (fpmpi, "\nr->seg->sIndex = %d", m->sIndex);
							fprintf (fpmpi, "\nr->seg->sStart = %d", m->sStart);
							fprintf (fpmpi, "\nr->seg->sEnd = %d\n", m->sEnd);	
							}
						
						}
					w = out = 0;


					}
			

				free (startsVectorRec);
				free (endsVectorRec);
				a = b = aa = bb = aaa = bbb = 0;


				/* intra codon breakpoints readjust MRCA oct2009 */
				if (doBreakpBroken == YES) /* The codon positions that recombinaed are increased a unit */
					{
					
					for (w = 1; w <= numNuc; w++)
						{
						if (w >= LeftLess2 && w <= RightHigh2)
							{
							S_MRCA[w]++;
							/*fprintf (fpmpi, "\n Increasing MRCA to %d: now is: %d \n", w, S_MRCA[w]);*/
							}
						}

					}
				/* oct2009 */



				if (noisy == 4)
					{
					fprintf (fpmpi, "\n>> Recombination Results:");
					fprintf (fpmpi, "\nNew left node with %d fragment(s)", q->numSegNode);				
					}
				for (w = 0; w < q->numSegNode;w++)
					{
					n = segments + post(w,q->index,distance);
					if (noisy == 4)
						{
						fprintf (fpmpi, "\nq->seg->sIndex = %d", n->sIndex);
						fprintf (fpmpi, "\nq->seg->sStart = %d", n->sStart);
						fprintf (fpmpi, "\nq->seg->sEnd = %d\n", n->sEnd);				
						}
					}
				
				if (noisy == 4)
					fprintf (fpmpi, "\nNew right node with %d fragment(s)", r->numSegNode);
					
					
				for (w = 0; w < r->numSegNode;w++)
					{
					m = segments + post(w,r->index,distance);
					if (noisy == 4)
						{
						fprintf (fpmpi, "\nr->seg->sIndex = %d", m->sIndex);
						fprintf (fpmpi, "\nr->seg->sStart = %d", m->sStart);
						fprintf (fpmpi, "\nr->seg->sEnd = %d\n", m->sEnd);						
						}
					}
				if (noisy > 3)
					fprintf (fpmpi, "\n");
					
				
				/*fprintf (stderr, "the node is whichInd = %d, and the site is whichSite = %d", whichInd+1, whichSite+1);*/
				if (noisy > 1)
					{
					fprintf (fpmpi, "Recombination involving %d (copied to %d and %d)", p->index, q->index, r->index );
					fprintf (fpmpi, "\n Breakpoint was at site %d", whichSite);
					}
				if (noisy > 3)
					fprintf (fpmpi, "\n");
				
				if (doBranchNetfiles == YES)
					{
					fprintf(fpBranchNet,"%d %d\n", q->NetLabelPrint, p->NetLabelPrint);
					/*fprintf(fpBranchNet,"%d %d\n", r->NetLabelPrint, p->NetLabelPrint);*/
					fprintf (fpmpi, "\nNET INFORMATION Recombination involving %d (copied to %d and %d)", p->NetLabelPrint, q->NetLabelPrint, r->NetLabelPrint);
					}
				breakpoint[numRE-1] = whichSite;  /* the breakpoint site is call breakpoint. breakpoint[0] = 7, breakpoint[1] = 96, breakpoint[2] = 187.. */
				
				/* readjust active sites */
				activeGametes[whichInd] = firstHalf;	/* new active nodes firstHalf and secondHalf */
				activeGametes[numActiveGametes++] = secondHalf;		/* there are 1 active node more (in recombination) */
				}
		
		
			/*** if COALESCENCE, readjust nodes and pointers ***/
			if (isCoalescence == YES)
				{
				if (noisy == 4)
					fprintf (fpmpi, "\n* Coalescence *\n");
					
				/* figure out which two nodes are involved */ 
			
				/* intial nodes: firstInd and secondInd (they are the descendants). newInd is the ancestral node, is the new node to make */
				firstInd = numActiveGametes * RandomUniform(seed);
				if (firstInd >= numActiveGametes)
					{
					fprintf (fpmpi, "\n\nERROR: firstInd out of range!\n");
					exit (-1);
					}
				do
					{
					secondInd = numActiveGametes * RandomUniform(seed);
					} while (firstInd == secondInd);	/* the new nodes must to be diferents */
				
				/*debug*/
				/*if (doDatedTips == YES)
					fprintf (stderr, "\nsecondInd =%d, numActiveGametes=%d\n", secondInd, numActiveGametes);*/

				newInd = nextAvailable;
				if (noisy > 1)					
					fprintf (fpmpi, "Coalescence involving %d and %d to create node %d", activeGametes[firstInd], activeGametes[secondInd], newInd);
					
					
				p = nodes + activeGametes[firstInd];
				q = nodes + activeGametes[secondInd];
				r = nodes + newInd;		/* new ancester */
				r->index = nextAvailable;
				r->label = labelNodes++;
				
				/*fprintf (fpmpi, "\n\nCoalescence node label %d, y index %d\n\n", r->label, r->index);*/
				
				r->breakp = NO;
				r->breakCodon = NO;
				r->class = 4;
				r->GMRCA_ancestral = NO;
				

				for (mmm = 1; mmm <= numNuc; mmm++)	 /* for each segment of p node going to r node*/
					{
					sigue = 0;
					stateHere_P = -2; /* -2, non ancestral; 0 ancestral */
					stateHere_Q = -2; /* -2, non ancestral; 0 ancestral */

					for (i = 0; i < p->numSegNode; i++)	 /* for each segment of p node going to r node*/
						{
						s = segments + post(i,p->index,distance);
						if (mmm >= s->sStart && mmm <= s->sEnd) /* is ancestral material */
							{
							stateHere_P = 0;
							}
						}
					for (i = 0; i < q->numSegNode; i++)	/* for each segment of q node going to r node */
						{
						n = segments + post(i,q->index,distance);
						if (mmm >= n->sStart && mmm <= n->sEnd) /* is ancestral material */
							{
							stateHere_Q = 0;
							}
						}
					
					
					/*fprintf (fpmpi, "\n Here(%d) stateHere_P = %d and stateHere_Q = %d ", mmm, stateHere_P, stateHere_Q);
					fprintf (fpmpi, "\n  Here(%d), node %d: p->SitesNonAncHere[mmm] = %d && node %d: q->SitesNonAncHere[mmm] = %d \n", mmm, p->index, p->SitesNonAncHere[mmm], q->index, q->SitesNonAncHere[mmm]);*/

					if (stateHere_P < 0 && stateHere_Q < 0) /* non ancestral material */
						{
						r->SitesNonAncHere[mmm] = -1;
						sigue++;
						/*fprintf (fpmpi, "1Position %d is NON anc mat \n", mmm);*/
						}
					if (stateHere_P >= 0 || stateHere_Q >= 0) /* ancestral material (inc pseudo) */
						{
						if (p->SitesNonAncHere[mmm] == 1 && q->SitesNonAncHere[mmm] == 1 && sigue == 0) /* pseudo anc mat */
							{
							r->SitesNonAncHere[mmm] = 1;
							sigue++;
							/*fprintf (fpmpi, "2Position %d is PSEUDO anc mat \n", mmm);*/
							}
						if (p->SitesNonAncHere[mmm] == 1 && stateHere_Q < 0 && sigue == 0) /* pseudo anc mat */
							{
							r->SitesNonAncHere[mmm] = 1;
							sigue++;
							/*fprintf (fpmpi, "3Position %d is PSEUDO anc mat \n", mmm);*/
							}
						if (q->SitesNonAncHere[mmm] == 1 && stateHere_P < 0 && sigue == 0) /* pseudo anc mat */
							{
							r->SitesNonAncHere[mmm] = 1;
							sigue++;
							/*fprintf (fpmpi, "4Position %d is PSEUDO anc mat \n", mmm);*/
							}
						if (sigue == 0)
							{
							if (p->SitesNonAncHere[mmm] == 0 || q->SitesNonAncHere[mmm] == 0)  /* anc mat */
								{
								r->SitesNonAncHere[mmm] = 0;
								sigue++;
								/*fprintf (fpmpi, "5Position %d is ANC mat \n", mmm);*/
								}
							}
						}
					}




				/*fprintf (fpmpi, "\nobtained r->index = %d: r->breakp = %d, r->breakCodon = %d, r->class = %d \n", r->index, r->breakp, r->breakCodon, r->class);*/
				
				if (doBranchNetfiles == YES)
					{
					r->NetLabelPrint = numNetLabelPrint;
					numNetLabelPrint++;

					if (p->NetLabelPrint == q->NetLabelPrint)
						{
						fprintf(fpBranchNet,"%d %d\n", r->NetLabelPrint, p->NetLabelPrint);
						}
					else
						{
						fprintf(fpBranchNet,"%d %d\n", r->NetLabelPrint, p->NetLabelPrint);
						fprintf(fpBranchNet,"%d %d\n", r->NetLabelPrint, q->NetLabelPrint);
						}
					fprintf (fpmpi, "\nNET INFORMATION Coalescence involving %d and %d to create node %d", p->NetLabelPrint, q->NetLabelPrint, r->NetLabelPrint);
					}


				coalVectorCountStarts = (int *) calloc(((p->numSegNode)+(q->numSegNode)),(long) sizeof(int));
				if (!coalVectorCountStarts)
					{
					fprintf (fpmpi, "Could not allocate coalVectorCountStarts (%lu bytes)\n", ((p->numSegNode)+(q->numSegNode)) *(long) sizeof(int));
					exit (-1);
					}
				coalVectorCountEnds = (int *) calloc(((p->numSegNode)+(q->numSegNode)),(long) sizeof(int));
				if (!coalVectorCountEnds)
					{
					fprintf (fpmpi, "Could not allocate coalVectorCountEnds (%lu bytes)\n", ((p->numSegNode)+(q->numSegNode)) *(long) sizeof(int));
					exit (-1);
					}
				
				
				for (i = 0; i < p->numSegNode; i++)	 /* for each segment of p node going to r node*/
					{
					s = segments + post(i,p->index,distance);
					m = segments + post(i,r->index,distance);					/* new ancester */
					m->sIndexNode = newInd;
					
					s->before1 = m;
					m->before1 = NULL;
					m->before2 = NULL;
					m->after1 = s;
					m->after2 = NULL;
					m->sIndex = actSegIndex;
					m->sStart = s->sStart;
					m->sEnd = s->sEnd;				

					actSegIndex++;
					actNumSegments++;
						
					if (m->sStart == 0 && m->sEnd == 0) /* unreal segments */
						{
						fprintf (fpmpi, "\nNot to be here. COAL1, segment start and end = 0");
						
						m->before1 = m->before2 = m->after1 = m->after2 = NULL;
						actNumSegments--;
						exit(-1);
						}
						
					coalVectorCountStarts[i] = s->sStart;
					coalVectorCountEnds[i] = s->sEnd;
					}
					
				j = p->numSegNode;
				r->numSegNode = j;
				a = b = 0;
				
				for (i = 0; i < q->numSegNode; i++)	/* for each segment of q node going to r node */
					{
					n = segments + post(i,q->index,distance);
					
					for (w = 0; w < j; w++)
						if (n->sStart == coalVectorCountStarts[w] && n->sEnd == coalVectorCountEnds[w]) /* Repeated segment */
							a++;
						
					if (a == 0)
						{
						m = segments + post(j+b,r->index,distance);					/* new ancester */
						r->numSegNode = r->numSegNode+1;
						
						n->before1 = m;
						m->before1 = NULL;
						m->before2 = NULL;
						m->after1 = n;
						m->after2 = NULL;
						m->sIndex = actSegIndex;
						m->sStart = n->sStart;
						m->sEnd = n->sEnd;				
						m->sIndexNode = newInd;
						
						actSegIndex++;
						b++;
						actNumSegments++;
					
						if (m->sStart == 0 && m->sEnd == 0) /* Cheking */ /* unreal segments */
							{
							fprintf (fpmpi, "\nNot to be here. COAL2, segment start and end = 0");
							
							m->before1 = m->before2 = m->after1 = m->after2 = NULL;
							actNumSegments--;
							}
						}
					a = 0;
					}
				a = b = 0;
				free (coalVectorCountStarts);
				free (coalVectorCountEnds);
								

					/* Segment Bonds when this segment goes to 2 descendants segments */
				coalEqualSegInit_p = (int *) calloc((p->numSegNode+q->numSegNode),(long) sizeof(int));
				if (!coalEqualSegInit_p)
					{
					fprintf (fpmpi, "Could not allocate coalEqualSegInit_p (%lu bytes)\n", (p->numSegNode+q->numSegNode) *(long) sizeof(int));
					exit (-1);
					}
				coalEqualSegEnd_p = (int *) calloc((p->numSegNode+q->numSegNode),(long) sizeof(int));
				if (!coalEqualSegEnd_p)
					{
					fprintf (fpmpi, "Could not allocate coalEqualSegEnd_p (%lu bytes)\n", (p->numSegNode+q->numSegNode) *(long) sizeof(int));
					exit (-1);
					}
				coalEqualSegInit_q = (int *) calloc((p->numSegNode+q->numSegNode),(long) sizeof(int));
				if (!coalEqualSegInit_q)
					{
					fprintf (fpmpi, "Could not allocate coalEqualSegInit_q (%lu bytes)\n", (p->numSegNode+q->numSegNode) *(long) sizeof(int));
					exit (-1);
					}
				coalEqualSegEnd_q = (int *) calloc((p->numSegNode+q->numSegNode),(long) sizeof(int));
				if (!coalEqualSegEnd_q)
					{
					fprintf (fpmpi, "Could not allocate coalEqualSegEnd_q (%lu bytes)\n", (p->numSegNode+q->numSegNode) *(long) sizeof(int));
					exit (-1);
					}
				
				for (w = 0; w < p->numSegNode; w++)	
					{
					s = segments + post(w,p->index,distance);
						
					coalEqualSegInit_p[w] = s->sStart;
					coalEqualSegEnd_p[w] = s->sEnd;
					}
				for (w = 0; w < q->numSegNode; w++)	
					{
					n = segments + post(w,q->index,distance);
						
					coalEqualSegInit_q[w] = n->sStart;
					coalEqualSegEnd_q[w] = n->sEnd;
					}
					
				for (w = 0; w < p->numSegNode; w++)
					{
					for (i = 0; i < q->numSegNode; i++)
						{
						if (coalEqualSegInit_p[w] != 0 && coalEqualSegEnd_p[w] != 0 && coalEqualSegInit_p[w] == coalEqualSegInit_q[i] && coalEqualSegEnd_p[w] == coalEqualSegEnd_q[i])
							{
							for (a = 0; a < r->numSegNode; a++)	
								{
								m = segments + post(a,r->index,distance);									
								if (m->sStart == coalEqualSegInit_p[w] && m->sEnd == coalEqualSegEnd_p[w])
									{
									if (noisy == 4)
										fprintf (fpmpi, "\nFragment %d links to 2 descendants, ", m->sIndex);
										
									for (b = 0; b < p->numSegNode; b++)
										{
										s = segments + post(b,p->index,distance);												
										if (coalEqualSegInit_p[w] == s->sStart && coalEqualSegEnd_p[w] == s->sEnd)
											{
											m->after1 = s;
											if (noisy == 4)
												fprintf (fpmpi, "fragment %d", s->sIndex);
											}
										}
									for (b = 0; b < q->numSegNode; b++)
										{
										n = segments + post(b,q->index,distance);											
										if (coalEqualSegInit_p[w] == n->sStart && coalEqualSegEnd_p[w] == n->sEnd)
											{
											m->after2 = n;
											if (noisy == 4)
												fprintf (fpmpi, " and fragment %d", n->sIndex);
												
											}
										}
									}
								}
							}
						}
					}
						
				
				for (i = 0; i < r->numSegNode; i++)			/* cheking, if there are 2 similar segments, it keeps only 1 */		
					{
					m = segments + post(i,r->index,distance);					
					for (w = 0; w < r->numSegNode; w++)
						{
						z = segments + post(w,r->index,distance);						
						if (w != i && m->sStart == z->sStart && m->sEnd == z->sEnd && m->after1 != NULL)
							{
							fprintf (fpmpi, "\n2COAL. Not to be here. It cannot have 2 equal segments en this node.");
							fprintf (fpmpi, "\n m->sIndex = %d, z->sIndex = %d \n", m->sIndex, z->sIndex);
							fprintf (fpmpi, "\n m->sIndex = %d, z->sIndex = %d \n", m->sIndex, z->sIndex);
							
							z->after1 = z->after2 = z->before1 = z->before2 = NULL; /* continue only the z segment */
							r->numSegNode--;
							exit(-1);
							}
						}
					}
				j = w = a = b = 0;
				free (coalEqualSegInit_p);
				free (coalEqualSegEnd_p);
				free (coalEqualSegInit_q);
				free (coalEqualSegEnd_q);
				
					
				
				if (noisy == 4)
					{
					fprintf (fpmpi, "\n\nCoalescence Result is the new node %d with %d fragment(s)", r->index, r->numSegNode);
					for (i = 0; i < r->numSegNode; i++)	
						{
						m = segments + post(i,r->index,distance);
					
						fprintf (fpmpi, "\nFragment %d, ", m->sIndex);
						fprintf (fpmpi, "m->sStart = %d", m->sStart);
						fprintf (fpmpi, " and m->sEnd = %d", m->sEnd);
						}
					fprintf (fpmpi, "\n");
					}
				
				
		
				/** MRCA **/
				/* Array with the information about ancestral stuff */
				/* S_MRCA. overFirst is the biggest site by left & overEnd is the smallest site by right. The difference is the overLapSites */
				sizeNode_p = p->numSegNode;
				sizeNode_q = q->numSegNode;
				out = 0;
			
				if (r->numSegNode == 1) /* 1.- fast. Only 1 segment */
					{
					for (i = 0; i < r->numSegNode; i++)	
						m = segments + post(i,r->index,distance);
				
					overFirst = m->sStart;
					overEnd = m->sEnd;
				
					for (j = 1; j <= numNuc; j++)
						{
						if (j >= overFirst && j <= overEnd)
							{
							S_MRCA[j]--;
							}
						}
					out = 1;
					}
				a = b = k = 0;
			
				if (r->numSegNode > 1)	
					{
					initialVector_pp = (int *) calloc((p->numSegNode),(long) sizeof(int));
					if (!initialVector_pp)
						{
						fprintf (fpmpi, "Could not allocate initialVector_pp (%lu bytes)\n", (p->numSegNode) *(long) sizeof(int));
						exit (-1);
						}
					endVector_pp = (int *) calloc((p->numSegNode),(long) sizeof(int));
					if (!endVector_pp)
						{
						fprintf (fpmpi, "Could not allocate endVector_pp (%lu bytes)\n", (p->numSegNode) *(long) sizeof(int));
						exit (-1);
						}
					initialVector_qq = (int *) calloc((q->numSegNode),(long) sizeof(int));
					if (!initialVector_qq)
						{
						fprintf (fpmpi, "Could not allocate initialVector_qq (%lu bytes)\n", (q->numSegNode) *(long) sizeof(int));
						exit (-1);
						}
					endVector_qq = (int *) calloc((q->numSegNode),(long) sizeof(int));
					if (!endVector_qq)
						{
						fprintf (fpmpi, "Could not allocate endVector_qq (%lu bytes)\n", (q->numSegNode) *(long) sizeof(int));						
						exit (-1);
						}
				
				
					/* 2.- Fast. There is a segment more big than the other segments */
					for (position = 0; position < p->numSegNode; position++)
						{
						s = segments + post(position,p->index,distance);
						initialVector_pp[position] = s->sStart;
						endVector_pp[position] = s->sEnd;
						}
					for (position = 0; position < q->numSegNode; position++)
						{
						n = segments + post(position,q->index,distance);		
						initialVector_qq[position] = n->sStart;
						endVector_qq[position] = n->sEnd;
						}
					
					/* smaller value for the initial */
					for (position = 0; position < p->numSegNode; position++)
						{
						if (position == 0)
							minInit_pp = initialVector_pp[position];
		
						if (initialVector_pp[position] < minInit_pp && initialVector_pp[position] != 0)
							minInit_pp = initialVector_pp[position];
						}

					/* older value for the end */
					for (position = 0; position < p->numSegNode; position++)
						{
						if (position == 0)
							maxEnd_pp = endVector_pp[position];
	
						if (endVector_pp[position] > maxEnd_pp && endVector_pp[position] != 0)
							maxEnd_pp = endVector_pp[position];
						}
		
						/* smaller value for the initial */
					for (position = 0; position < q->numSegNode; position++)
						{
						if (position == 0)
							minInit_qq = initialVector_qq[position];
				
						if (initialVector_qq[position] < minInit_qq && initialVector_qq[position] != 0)
							minInit_qq = initialVector_qq[position];
						}

					/* older value for the end */
					for (position=0; position < q->numSegNode; position++)
						{
						if (position == 0)
							maxEnd_qq = endVector_qq[position];
	
						if (endVector_qq[position] > maxEnd_qq && endVector_qq[position] != 0)
							maxEnd_qq = endVector_qq[position];
						}
			
	
					/* fast, there are in 1 big segment that contains all the segments */					
					for (position = 0; position < p->numSegNode; position++) 
						{
						if (endVector_pp[position] == maxEnd_pp && initialVector_pp[position] == minInit_pp)	
							a = 1;
						}
					for (position = 0; position < q->numSegNode; position++) 
						{
						if (endVector_qq[position] == maxEnd_qq && initialVector_qq[position] == minInit_qq)	
							b = 1;
						}
				
					free (initialVector_pp);
					free (endVector_pp);
					free (initialVector_qq);
					free (endVector_qq);
				
				
					/* 2.- Fast with a big segment */
					if (minInit_pp <= minInit_qq && maxEnd_pp >= maxEnd_qq && a == 1 && b == 1)
						{
						for (j = 1; j <= numNuc; j++)
							{
							if (j >= minInit_qq && j <= maxEnd_qq)
								{
								S_MRCA[j]--;
								}
							}
						a = out = 1;	
						}
					if (minInit_qq <= minInit_pp && maxEnd_qq >= maxEnd_pp && a == 1 && b == 1 && out == 0)
						{
						for (j = 1; j <= numNuc; j++)
							{
							if (j >= minInit_pp && j <= maxEnd_pp)
								{
								S_MRCA[j]--;
								}
							}
						a = out = 1;	
						}
					
					if (maxEnd_pp < minInit_qq && out == 0)
						a = out = 1;
						
					if (maxEnd_qq < minInit_pp && out == 0)
						a = out = 1;
						
					
					if (minInit_pp < minInit_qq && maxEnd_pp < maxEnd_qq && a == 1 && b == 1 && out == 0)
						{
						for (j = 1; j <= numNuc; j++)
							{
							if (j >= minInit_qq && j <= maxEnd_pp)
								{
								S_MRCA[j]--;
								}
							}
						a = out = 1;
						}
						
					if (minInit_qq < minInit_pp && maxEnd_qq < maxEnd_pp && a == 1 && b == 1 && out == 0)
						{
						for (j = 1; j <= numNuc; j++)
							{
							if (j >= minInit_pp && j <= maxEnd_qq)
								{
								S_MRCA[j]--;
								}
							}
						a = out = 1;
						}
					
					a = b = k = w = i = 0;
				
					/* 3.- complex case */
					if (out == 0)
						{
						for (j = 1; j <= numNuc; j++)
							{
							if (overLapSegmentsCoalMRCA(p, q, sizeNode_p, sizeNode_q, j) == YES)
								{
								S_MRCA[j]--;
								}
							}
						}
					}
				j = 0;




				/* Folllowing only ancestral material */
				for (mmm = 1; mmm <= numNuc; mmm++)	 /* Only anc material*/
					{
					if (p->SitesNonAncHere[mmm] == 0 && q->SitesNonAncHere[mmm] == 0)
						{
						/*fprintf (fpmpi, "\n This Position %d is ANC mat: p->index = %d, p->SitesNonAncHere = %d, q->index = %d, q->SitesNonAncHere = %d \n", mmm, p->index, p->SitesNonAncHere[mmm], q->index, q->SitesNonAncHere[mmm]);*/
						OnlyAncS_MRCA[mmm]--;
						}				
					}
				/* Is this the GMRCA of the anc material? */
				Ok_SMRCA_Codon = 0;
				for (mmm = 1; mmm <= numNuc; mmm++)	 /* Only anc material*/
					{
					if (OnlyAncS_MRCA[mmm] == 1)
						{
						/*fprintf (fpmpi, "\n-- Yes anc GMRCA position %d -- \n", mmm);*/
						Ok_SMRCA_Codon++;
						}
					}			
				if (Ok_SMRCA_Codon == numNuc && AncGMRCA_obtained == NO)
					{
					r->GMRCA_ancestral = YES;
					AncGMRCA_obtained = YES;
					if (noisy == 4)
						fprintf (fpmpi, "\n--GMRCA of the ancestral material in node %d-- \n", r->index);
					}



				/* in coalescence is possible link the nodes */
				r->left = p;
				r->right = q;
				p->anc1 = r;
				q->anc1 = r;
				r->time = currentTime;
				/*fprintf (fpmpi, " r->index = %d, r->time = %lf\n", r->index, r->time);*/

				/* readjust active nodes */
				activeGametes[firstInd] = newInd;
				activeGametes[secondInd] = activeGametes[numActiveGametes-1];
				numActiveGametes--; /* it lose 1 active node */
				nextAvailable++; /* 1 node more to available */
			
				if (nextAvailable >= numNodes)	/* if there aren't enough nodes it go into and it addition more */
					{
					/* ReallocNodes(&numNodes, activeGametes); */
					numNodes += INCREMENT_NODES;
					numTotalSegments += (INCREMENT_NODES*maxSegNode)+numNuc;
					
					/* REALLOC */
					segments = (TreeSegment *) realloc (segments, numTotalSegments  * (long) sizeof(TreeSegment)); 
					if (!segments)
						{
						fprintf (fpmpi, "Could not reallocate segments (%lu bytes)\n", ((numNodes*distance)+numNuc)  * (long) sizeof(TreeSegment));
						exit (1);
						}
					nodes = (TreeNode *) realloc (nodes, numNodes  * (long) sizeof(TreeNode));
					if (!nodes)
						{
						fprintf (fpmpi, "Could not reallocate nodes (%lu bytes)\n", numNodes  * (long) sizeof(TreeNode));
						exit (-1);
						}
					activeGametes = (int *) realloc (activeGametes, numNodes *(long) sizeof(int));
					if (!activeGametes)
						{
						fprintf (fpmpi, "Could not reallocate activeGametes (%lu bytes)\n", numNodes *(long) sizeof(int));
						exit (-1);
						}
					if (noisy == 4)
						fprintf (fpmpi, "\n\n...Doing reallocation of nodes (1)\n");
						
					}
				}	/* end of coalescence */
		
						
			/* print out ancestral (active) status for each site and MRCA vector - only mat anc */
			if (noisy > 2)
				{
				fprintf (fpmpi,"\n - Ancestral MRCA in the nodes: -\n");
				for (i=0; i<numActiveGametes; i++)
					{
					for (j=1; j<=numNuc; j++)
						{
						p = nodes + activeGametes[i];
						if (j == 1)
							fprintf (fpmpi, "%4d -- (MRCA:)", p->index);
						fprintf (fpmpi, "%d", OnlyAncS_MRCA[j]);
						}
					fprintf (fpmpi, "\n");
					}
				/*fprintf (fpmpi, "MRCA   ");*/
				for (j=1; j<=numNuc; j++)
					{
					if (OnlyAncS_MRCA[j] <= 1)
						fprintf (fpmpi, "*");
					else
						fprintf (fpmpi, " ");
					}
				/*fprintf (fpmpi, "\n");*/
				fprintf (fpmpi,"\n\n");
				}

			/* print out ancestral (active) status for each site and MRCA vector */
			if (noisy > 2)
				{
				fprintf (fpmpi,"\nMRCA in the nodes:\n");
				for (i=0; i<numActiveGametes; i++)
					{
					for (j=1; j<=numNuc; j++)
						{
						p = nodes + activeGametes[i];
						if (j == 1)
							fprintf (fpmpi, "%4d -- (MRCA:)", p->index);
						fprintf (fpmpi, "%d", S_MRCA[j]);
						}
					fprintf (fpmpi, "\n");
					}
				/*fprintf (fpmpi, "MRCA   ");*/
				for (j=1; j<=numNuc; j++)
					{
					if (S_MRCA[j] <= 1)
						fprintf (fpmpi, "*");
					else
						fprintf (fpmpi, " ");
					}
				/*fprintf (fpmpi, "\n");*/
				}
			
			
			sizeNode = sizeNode_p = sizeNode_q = 0;
			free (gi);

			

			/* If doing dated tips, sometimes the number of actives nodes can be already 1 before all have been activated samples. 
			If this happens we need to move towards the next sample, activate it, and start again */
			if (doDatedTips == YES)
			 if (numActiveGametes == 1 && currentSample > 0)
				{
				if (noisy > 2)
					fprintf (stderr, "\n\nOnly 1 lineage active before activating all samples"/*, currentSample*/);
				currentSample--;

				/* activate nodes from this sample*/
				if (noisy > 2)
					fprintf (stderr, "\n Activating sample %d (time = %6.4f). Tips to activate:", currentSample, datedSample[currentSample].time);
				for (i=0; i<datedSample[currentSample].size; i++)
					{
					if (noisy > 2)
						fprintf (stderr, " %d", datedSample[currentSample].member[i]-1);
					
					p = nodes + datedSample[currentSample].member[i]-1;
					p->index = datedSample[currentSample].member[i]-1;
					activeGametes[numActiveGametes] = datedSample[currentSample].member[i]-1;  
						
					numActiveGametes++;
					}
				currentTime = datedSample[currentSample].time;
	
				if (noisy > 2)
					{
					fprintf (stderr, "\nActive nodes (%d):", numActiveGametes); 
					for (i=0; i<numActiveGametes; i++)
						fprintf (stderr," %d",activeGametes[i]);
					fprintf (stderr,"   Next node available = %d ", nextAvailable);
					fprintf (stderr, "\nSetting currentTime = %6.4f", currentTime);
					}
				}




			} /* coalescent tree finished */

		
		 /* Cheking */
		/*for (w = 1; w <= numNuc; w++)
			{
			if (S_MRCA[w] > 1 || S_MRCA[w] < 1)
				{*/
				/*
				fprintf (fpmpi, "\n Warning S_MRCA in the last node is > 1, S_MRCA[%d] = %d", w, S_MRCA[w]);
				exit (-1);
				}
			}*/


		/*if (noisy > 1)			
			fprintf (fpmpi, "\n\n\n>> Coalescent tree/s finished\n");*/
		if (noisy > 1)			
			fprintf (fpmpi, "\n\n\n>> Coalescent finished\n");
			
				
		counterTimeInit = counterTimeInit + currentTime;
		actualTGMRCA = currentTime;


		/* Use it to see the evolution nodes-segments */
		
		/*
		fprintf (fpmpi, "\n\n\nnextAvailable = %d",nextAvailable);
		for (j = 0; j < numNodes; j++)*/	/* Looking for parent Node of segments */	
		/*	{
			p = nodes + j;
			fprintf (fpmpi, "\n\n\n\n**The node %d**",p->index);

			for (i=0; i < p->numSegNode; i++)
				{
				s = segments + post(i,j,numNuc);	
				p->seg = s;
			
				fprintf (fpmpi, "\n\ns->sIndex %d con s->parentNode->index = %d", s->sIndex, s->parentNode->index);
				if (s->after1 != NULL)
					fprintf (fpmpi, "\ns->after1->sIndex %d con s->after1->parentNode->index = %d", s->after1->sIndex, s->after1->parentNode->index);
				if (s->after2 != NULL)
					fprintf (fpmpi, "\ns->after2->sIndex %d con s->after2->parentNode->index = %d", s->after2->sIndex, s->after2->parentNode->index);
				}
			}*/
			
		}
	
		



	
	/********** BUILDING TREES ***********/
	/* first, we build an array with order breakpoints from - to + */
	free (activeGametes);
	
	if (doFixNumRecEvents == NO || numRE == fixedNumRecEvents) /* Good replicate */
		{
		if (noisy > 1)
			{
			fprintf (fpmpi, "\n\n>> Net finished\n");
			fprintf (fpmpi, "\n\n>> Building trees ..");
			}
		if (noisy == 4 && numRE > 0)
			fprintf (fpmpi, "\n\n");
		
		
		arrayIndBreakpoints = (int *) calloc((numRE+1),(long) sizeof(int));
		if (!arrayIndBreakpoints)
			{
			fprintf (fpmpi, "Could not allocate arrayIndBreakpoints (%lu bytes)\n", (numRE+1) *(long) sizeof(int));
			exit (-1);
			}
		arrayIndBreakpointsOrd = (int *) calloc((numRE+1),(long) sizeof(int));
		if (!arrayIndBreakpointsOrd)
			{
			fprintf (fpmpi, "Could not allocate arrayIndBreakpointsOrd (%lu bytes)\n", (numRE+1) *(long) sizeof(int));	
			exit (-1);
			}
		for (w = 0; w < numRE; w++)
			{
			arrayIndBreakpoints[w] = breakpoint[w];
			arrayIndBreakpointsOrd[w] = 0;
			
			/*
			fprintf (fpmpi, "arrayIndBreakpoints[%d] = %d\n", w, arrayIndBreakpoints[w]);
			*/
			}
		
		/*
		fprintf (fpmpi, "arrayIndBreakpoints[%d+1] = %d\n", w+1, arrayIndBreakpoints[w+1]);	
		*/	
	
		j = arrayIndBreakpoints[0];
		for (k = 0; k < numRE; k++)
			{
			for (w = 0; w < numRE; w++)		/* order breakpoints, from - to + */
				if (arrayIndBreakpoints[w] < j)
					j = arrayIndBreakpoints[w];
				
			if (j > numNuc) /* end */
				{
				arrayIndBreakpointsOrd[k] = 0;
				break;
				}
					
			for (w = 0; w < numRE; w++)
				if (arrayIndBreakpoints[w] == j)
					arrayIndBreakpoints[w] = numNuc+1;
		
			arrayIndBreakpointsOrd[k] = j;
			j = numNuc+1;
			}
		free (arrayIndBreakpoints);
	
		indNumRE = 0;
		if (noisy == 4 && numRE > 0)
			fprintf (fpmpi, "Breakpoints list (from the smallest to the biggest value):\n");
			
			
		for (k = 0; k < numRE; k++)
			{
			if (arrayIndBreakpointsOrd[k] != 0)
				indNumRE++;
			if (noisy == 4 && numRE > 0 && arrayIndBreakpointsOrd[k] != 0)
				fprintf (fpmpi, "%d\n", arrayIndBreakpointsOrd[k]);
			}
		
		for (k = 0; k < indNumRE; k++)
			{
			if (k != 0)
				if (arrayIndBreakpointsOrd[k] < arrayIndBreakpointsOrd[k-1])
					{
					fprintf (fpmpi, " \n\n Warning in arrayIndBreakpointsOrd[k], 2 similar breakpoints in this array. See the indepBrekp vector");
					/*exit (-1);*/ /* but no problem */
					}
			}
	
		/* about the GMRCA */
		p = nodes + newInd; /* because the last one process is the last one coalescence */
		p->class = 5;
		treeRootInit[0] = p;
		

		/* The trees */
		if (noisy > 3)
			fprintf (fpmpi, "\nIt generates %d tree(s):", indNumRE+1);
			
		j = k = w = countNumTrees = 0;
	
		numSegTrees = (int *) calloc((indNumRE+2),(long) sizeof(int));
		if (!numSegTrees)
			{
			fprintf (fpmpi, "Could not allocate numSegTrees (%lu bytes)\n", (indNumRE+2) *(long) sizeof(int));
			exit (-1);
			}
				/* determination of segments number of each tree, from the GMRCA node */
		for (i = 0; i < p->numSegNode; i++)
			{
			s = segments + post(i,p->index,distance);			
			j = 0;
			
			if (s->sStart == 1 && s->sEnd >= (arrayIndBreakpointsOrd[j]-1) && s->after1 != NULL) /* first tree */
				numSegTrees[j]++;
			if (s->sStart <= arrayIndBreakpointsOrd[indNumRE-1] && s->sEnd == numNuc && s->after1 != NULL) /* last tree */
				numSegTrees[indNumRE]++;
			if (s->after1 == NULL && s->after2 == NULL) /* no conect segments, it must to be 0 */
				numSegTrees[indNumRE+1]++;
			if (indNumRE > 1)
				{
				for (j = 0; j < (indNumRE-1); j++)
					if (s->sStart <= (arrayIndBreakpointsOrd[j]) && s->sEnd >= (arrayIndBreakpointsOrd[j+1]-1) && s->after1 != NULL) /* internal trees */
						numSegTrees[j+1]++;
				}
			}
		
	
		if (noisy == 4)
			for (i = 0; i < (indNumRE+1); i++)
				fprintf (fpmpi, "\nTree %d starts by %d fragment(s)", i+1, numSegTrees[i]);
			
		
		
			/*** building trees ***/
		j = k = w = a = b = 0;
		/* default */
		numNodex = (indNumRE+1)*numNodes;
		if (indNumRE > 1000)
			numNodex = numNodex/5;
		if (indNumRE > 10000)
			numNodex = numNodex/50;
		if (indNumRE > 100000)
			numNodex = numNodex/50000;		
		/*fprintf (stderr, "\nindNumRE = %d, numNodex= %d", indNumRE, numNodex);*/
		
		nodex = (TreeNodex *) calloc (numNodex, sizeof(TreeNodex)); /* nodes */
		if (!nodex)
			{
			fprintf (fpmpi, "\nCould not allocate nodex1 (%lu bytes)\n", numNodex  * (long) sizeof(TreeNodex));
			
			fprintf (fpmpi, "\nSecond try \n");
			free (nodex);
			numNodex = (indNumRE+1)*numNodes/2;
			if (indNumRE > 1000)
				numNodex = numNodex/15;
			if (indNumRE > 10000)
				numNodex = numNodex/100;
			if (indNumRE > 100000)
				numNodex = numNodex/500000;	

			nodex = (TreeNodex *) calloc (numNodex, sizeof(TreeNodex)); /* nodes */
			if (!nodex)
				{
				fprintf (fpmpi, "\nCould not allocate nodex2 (%lu bytes)\n", numNodex  * (long) sizeof(TreeNodex));

				fprintf (fpmpi, "\nThird try \n");
				free (nodex);
				numNodex = (indNumRE+1)*numNodes/5;
				if (indNumRE > 1000)
					numNodex = numNodex/30;
				if (indNumRE > 10000)
					numNodex = numNodex/200;
				if (indNumRE > 100000)
					numNodex = numNodex/700000;	

				nodex = (TreeNodex *) calloc (numNodex, sizeof(TreeNodex)); /* nodes */
				if (!nodex)
					{
					fprintf (fpmpi, "\nCould not allocate nodex3 (%lu bytes)\n", numNodex  * (long) sizeof(TreeNodex));

					fprintf (fpmpi, "\nFourth try \n");
					free (nodex);
					numNodex = (indNumRE+1)*numNodes/50;
					if (indNumRE > 1000)
						numNodex = numNodex/100;
					if (indNumRE > 10000)
						numNodex = numNodex/400;
					if (indNumRE > 100000)
						numNodex = numNodex/900000;	

					nodex = (TreeNodex *) calloc (numNodex, sizeof(TreeNodex)); /* nodes */
					if (!nodex)
						{
						fprintf (fpmpi, "\nCould not allocate nodex4 (%lu bytes)\n", numNodex  * (long) sizeof(TreeNodex));
						exit (1);
						}
					}
				}
			}

		treeRootNodex = (TreeNodex **) calloc(indNumRE+1, sizeof(TreeNodex *)); /* node pointers */
		if (!treeRootNodex)
			{
			fprintf (fpmpi, "Could not allocate treeRootNodex (%lu bytes)\n", indNumRE+1  * (long) sizeof(TreeNodex));			
			exit (1);
			}

		for (i = 0; i < numNodex; i++)
			{
			nodex[i].left = NULL;
			nodex[i].right = NULL;
			nodex[i].anc1 = NULL;
			nodex[i].outgroup = NULL;
			nodex[i].index = 0;
			nodex[i].NetIndex = -1;
			nodex[i].indexOldMigPop = 0;
			nodex[i].label = 0;
			nodex[i].isOutgroup = NO;
			nodex[i].length = 0.0;
			nodex[i].time = 0.0;
			nodex[i].MRCAfrom = -1;
			nodex[i].MRCAto = -1;
			}

		/* making trees */
		/*fprintf (stderr, "\n\n\n *********** MAKING TREES **********\n\n");*/
		f = nodex + a;
		f->index = 0;
		f->NetIndex = p->index;
		f->length = p->length;
		f->time = p->time;
		f->indexOldMigPop = p->indexOldMigPop;
		i = 0;
	
		if (indNumRE == 0) /* there aren't recombinations */
			{
			/*f = p;*/
			f->MRCAfrom = 1;
			f->MRCAto = numNuc;
			treeRootNodex[a] = f;
			/*treeRootNodex[a].MRCAfrom = 1;
			treeRootNodex[a].MRCAto = numNuc; */

			/*if (noisy >= 3)
				fprintf (fpmpi, "\nThere aren't recombinations. treeRootNodex[%d] = %d", a, treeRootNodex[a]->index);*/
			if (noisy > 2)
				fprintf (fpmpi,"\n\n>> Making the tree..");
				
			countNumTrees++;
			buildTreeCoal(p, f, numSequences, &numActNodex);
			numActNodex++;
			i = 0;
			}
	
		/* Several trees */
		if (indNumRE > 0) /* first tree */
			{
			/*f->MRCAfrom = 1;
			f->MRCAto = arrayIndBreakpointsOrd[j]-1;
			fprintf (fpmpi,"\n a = %d, arrayIndBreakpointsOrd[j]-1 = %d, f->index = %d, f->NetIndex = %d, f->time = %lf", a, arrayIndBreakpointsOrd[j]-1, f->index, f->NetIndex, f->time);
			fprintf (fpmpi,"\n in first tree numActNodex = %d, p->index = %d", numActNodex, p->index);*/

			treeRootNodex[a] = f;
			/*fprintf (fpmpi,"\n in first tree treeRootNodex[a]->index = %d", treeRootNodex[a]->index);*/
			treeRootNodex[a]->MRCAfrom = 1;
			treeRootNodex[a]->MRCAto = arrayIndBreakpointsOrd[j]-1; 
			i = a;
			if (noisy > 2)
				fprintf (fpmpi,"\n\n\n>> Making the first tree..");
				
		
			for (step = 0; step < p->numSegNode; step++)				
				{
				s = segments + post(step,p->index,distance);				
				if (s->sStart == 1 && s->sEnd >= (arrayIndBreakpointsOrd[j]-1) && s->after1 != NULL)
					{
					k++;
					if (k == numSegTrees[0])
						{
						countNumTrees++;
						buildTreeInit(p, f, numNuc, arrayIndBreakpointsOrd, j, numSequences, &numActNodex);
						}
					}
				}
			k = 0;
			numActNodex++;
		
			if (indNumRE > 1)		/* Internal trees */
				{
				for (j = 0; j < (indNumRE-1); j++)
					{
					/*fprintf (fpmpi,"\n Desde internal primero trees numActNodex = %d", numActNodex);*/
					f = nodex + numActNodex;
					f->index = numActNodex;
					f->NetIndex = p->index;
					f->length = p->length;
					f->time = p->time;
					f->indexOldMigPop = p->indexOldMigPop;
					numActNodex++;
					
					f->MRCAfrom = arrayIndBreakpointsOrd[j];
					f->MRCAto = arrayIndBreakpointsOrd[j+1]-1;
					treeRootNodex[j+1] = f;
					/*fprintf (fpmpi,"\n in intern tree treeRootNodex[j+1]->index = %d", treeRootNodex[j+1]->index);*/
					/*treeRootNodex[j+1].MRCAfrom = arrayIndBreakpointsOrd[j];
					treeRootNodex[j+1].MRCAto = arrayIndBreakpointsOrd[j+1]-1;*/
					
					if (noisy > 2)
						fprintf (fpmpi,"\n>> Making internal tree %d..", j+1);
						
						
					for (step = 0; step < p->numSegNode; step++)		
						{
						s = segments + post(step,p->index,distance);			 
						if (s->sStart <= (arrayIndBreakpointsOrd[j]) && s->sEnd >= (arrayIndBreakpointsOrd[j+1]-1) && s->after1 != NULL) /* internal trees, it can to be several trees */
							{
							k++;
							if (k == numSegTrees[j+1])	/* the last segment of internal tree */
								{
								countNumTrees++;
								buildTreeIntern(p, f, numNuc, arrayIndBreakpointsOrd, j, numSequences, &numActNodex);
								}
							}
						}
					k = 0;
					numActNodex++;
					}
				j = 0;
				}

		
				/* last tree */
			/*fprintf (fpmpi,"\n in Last tree primero numActNodex = %d", numActNodex);*/
			f = nodex + numActNodex;
			f->index = numActNodex;
			f->NetIndex = p->index;
			f->length = p->length;
			f->time = p->time;
			f->indexOldMigPop = p->indexOldMigPop;
		
			numActNodex++;

			f->MRCAfrom = arrayIndBreakpointsOrd[indNumRE-1];
			f->MRCAto = numNuc;			

			treeRootNodex[indNumRE] = f;
			/*fprintf (fpmpi,"\n in last tree treeRootNodex[indNumRE]->index = %d", treeRootNodex[indNumRE]->index);*/
			/*treeRootNodex[indNumRE].MRCAfrom = arrayIndBreakpointsOrd[indNumRE-1];
			treeRootNodex[indNumRE].MRCAto = numNuc;*/

		
			if (noisy > 2)
				fprintf (fpmpi,"\n>> Making the last tree..");
				
			for (step = 0; step < p->numSegNode; step++)		
				{
				s = segments + post(step,p->index,distance);
				if (s->sStart <= arrayIndBreakpointsOrd[indNumRE-1] && s->sEnd == numNuc && s->after1 != NULL)
					{
					k++;
					if (k == numSegTrees[indNumRE])	/* the last segment of internal tree */
						{
						countNumTrees++;
						j = indNumRE-1;
					
						buildTreeEnd(p, f, numNuc, arrayIndBreakpointsOrd, j, numSequences, &numActNodex);
						}
					}
				}
			free (numSegTrees);
			k = j = 0;
			numActNodex++;
			

			/* get rid of superflous nodes */
			a = 0;
			
			if (noisy > 3)
				fprintf (fpmpi, "\n\n");
			if (noisy > 3)
				fprintf (fpmpi, "\n>> Removing superfluos nodes .."); /* nodex */
				
			
			foundSuperflousNode = YES;
			while (foundSuperflousNode == YES)
				{
				foundSuperflousNode = NO;
		
				for (i = 0; i < numActNodex; i++) /* available all nodes */
					{
					f = nodex + i;
			
					if (f->left == NULL && f->right == NULL && f->anc1 == NULL)
						{
						/* nothing to do with this node because it is not connected to anything */
						}
					else if (f->left == NULL && f->right == NULL && f->anc1 != NULL)
						{
						/* do not do anything with this node because it is a tip */
						}
					else if (f->left != NULL && f->right == NULL && f->anc1 != NULL)
						{
						/* this is a superflous node and can be removed */
						foundSuperflousNode = YES;
					
						g = f->left;
						h = f->anc1;
						if (f->anc1->left == f) /* p->anc up, p->left down, total: up and down for left */
							{
							h->left = g;
							g->anc1 = h;
							f->left = NULL;
							f->anc1 = NULL;
							}
						else
							{
							h->right = g;
							g->anc1 = h;
							f->left = NULL;
							f->anc1 = NULL;
							}
						}
					else if (f->left == NULL && f->right != NULL && f->anc1 != NULL)
						{
						/* this is a superflous node and can be removed */
						foundSuperflousNode = YES;
			
						g = f->right;
						h = f->anc1;
					
						if (f->anc1->left == f)
							{
							h->left = g;
							g->anc1 = h;
							f->right = NULL;
							f->anc1 = NULL;
							}
						else
							{
							h->right = g;
							g->anc1 = h;
							f->right = NULL;
							f->anc1 = NULL;
							}
						}
					else if (f->left != NULL && f->right != NULL && f->anc1 != NULL)
						{
						/* this is an internal node formed by a coalescence event, do not touch */
						}
					else if (f->left != NULL && f->right != NULL && f->anc1 == NULL)
						{
						/*fprintf (fpmpi, "\nEn superfluos. f->index = %d\n",  f->index);*/
						/*fprintf (stderr, "\nGMRCA node %d", f->index);*/
						if (a <= indNumRE)					/* THIS IS VERY IMPORTANT, "HERE" WE GIVE ALL THE MRCAs OF NODEX */
							treeRootNodex[a] = f;
						a++;
					
						/* this is the last (coalescence event) in the tree, GMRCA */
						}
					else if (f->left != NULL && f->right == NULL && f->anc1 == NULL)
						{
						/* Seems to be the last coalescent event among sequences with non-ancestral material */	
						/* it is not superfluous, we just remove it */
						f->left->anc1 = NULL;
						}
					else if (f->left == NULL && f->right != NULL && f->anc1 == NULL)
						{
						/* not clear what this node could be doing, but we will remove it anyway */
						fprintf (fpmpi, "strange\n");
						
						f->left = NULL;
						f->right->anc1 = NULL;
						}
					else
						{
						fprintf (fpmpi, "You should not be here, I think\n");
						fprintf (fpmpi, "%d %d-- %d %d %d\n", IndexSeg(f), a, IndexSeg(f->left), IndexSeg(f->right), IndexSeg(f->anc1));
						exit (-1);
						}
					}
				}
			}
		
		if (noisy == 4)
			{
			fprintf (fpmpi, "\n\nProgram Information. List about the root nodes index: ");
			for (j = 0; j <= indNumRE; j++)		
				fprintf (fpmpi, "\ntreeRootNodex[%d]->index = %d", j, treeRootNodex[j]->index); /* careful with the superfluos nodes*/
			}
		
		/*fprintf (fpmpi, "\n\nProgram Information. List about the root nodes index: ");
		for (j = 0; j <= indNumRE; j++)		
			fprintf (fpmpi, "\ntreeRootNodex[%d]->index = %d, treeRootNodex[%d]->time = %lf, NetIndex = %d, from %d to %d", j, treeRootNodex[j]->index, j, treeRootNodex[j]->time, treeRootNodex[j]->NetIndex, treeRootNodex[j]->MRCAfrom, treeRootNodex[j]->MRCAto); */
		



		/* MRCAconcatenate part */
		NodesMRCAposit[0]= -1; /* en este vector se guardara: posicion 1 es el label (=index) del nodo MRCA del nuc 1, posicion 2... */
		/*fprintf (fpmpi, "\n\n");*/
		for (sss = 1; sss <= numNuc; sss++)
			{

			for (j = 0; j <= indNumRE; j++)
				{
				if (sss >= treeRootNodex[j]->MRCAfrom && sss <= treeRootNodex[j]->MRCAto)
					{
					saveThis = 0;
					saveThis = treeRootNodex[j]->NetIndex;
					NodesMRCAposit[sss] = saveThis;
					break;

					}
				}
			/*fprintf (fpmpi, "\n Site %d for node label %d ", sss, NodesMRCAposit[sss]);*/
			}




		
		
		if (thereisOutgroup == YES)   /* special node -Outgroup- */
			{
			if (noisy > 1)
				fprintf (fpmpi, "\n\n\n>> Attaching outgroup (everytree) .. \n");
				
			if (indNumRE == 0)
				{
				/* f is the outgroup node */
				numActNodex++;
				f = nodex + numActNodex;
		
				f->left = NULL;
				f->right = NULL;
				f->anc1 = treeRootNodex[0];
				f->time = 0;
				f->length = outgroupBranchLength/mutationRate;
				f->indexOldMigPop = 0;
				f->isOutgroup = YES;
				f->index = numActNodex;
				f->NetIndex = nextAvailable++/* -1 */;	
				/* p->label = p->index;*/
				f->label = numSequences;
				f->outgroup = NULL;		/* the outgroup can't has outgroup.. */
				treeRootNodex[0]->outgroup = f;
				}
			if (indNumRE > 0)
				{
				for (step = 0; step <= indNumRE; step++)
					{
					/* p is the outgroup node */
					numActNodex++;
					f = nodex + numActNodex;
			
					f->left = NULL;
					f->right = NULL;
					f->anc1 = treeRootNodex[step];
					f->time = 0;
					f->length = outgroupBranchLength/mutationRate;
					f->indexOldMigPop = 0;
					f->isOutgroup = YES;
					f->index = numActNodex;	
					f->NetIndex = nextAvailable++/* -1 */;
					/* p->label = p->index; */
					f->label = numSequences;
					f->outgroup = NULL;	/* the outgroup can't has outgroup.. */
					treeRootNodex[step]->outgroup = f;
					}
				}
			}
		else
			{
			if (noisy > 1)
				fprintf (fpmpi, "\n");
			}

		/* relabel nodes on tree */
		if (noisy > 1)
			fprintf (fpmpi, "\n\n>> Relabeling nodes (everytree) .. \n\n");

		/* tipLabel = 0; */
		if (indNumRE == 0)
			{
			/*tipLabel = 0;*/
			if (thereisOutgroup == YES)
				intLabel = numSequences+1;
			else
				intLabel = numSequences;
			RelabelNodesSeg(treeRootNodex[0]);
			}
		if (indNumRE != 0)
			{
			for (step = 0; step <= indNumRE; step++)
				{
				/*tipLabel = 0;*/
				if (thereisOutgroup == YES)
					intLabel = numSequences+1;
				else
					intLabel = numSequences;

				RelabelNodesSeg(treeRootNodex[step]);		
				}
			}
		/** NET RECODON **/
		if (thereisOutgroup == YES)   /* NET RECODON - OUTGROUP *//* special "node" -Outgroup- */
			{
			/* r is the outgroup node */
			r = nodes + nextAvailable;
		
			r->left = NULL;
			r->right = NULL;
			r->sib = NULL;
			r->class = 2;
			r->GMRCA_ancestral = NO;
			r->breakp = NO;
			r->breakCodon = NO;
			r->anc1 = treeRootInit[0];
			r->time = 0;
			r->length = outgroupBranchLength/mutationRate;
			r->indexOldMigPop = 0;
			r->isOutgroup = YES;
			r->index = nextAvailable++;
			//r->NetIndex = nextAvailable++/* -1 */;	
			/* p->label = p->index;*/
			//r->label = numSequences;
			r->outgroup = NULL;		/* the outgroup can't has outgroup.. */
			treeRootInit[0]->outgroup = r;				
			/*nextAvailable++;*/
			}
		/*fprintf (fpmpi, "\n\n nextAvailable = %d \n\n", nextAvailable);*/
		for (i=0; i<nextAvailable; i++) /* NET RECODON - labels*/
			{
			nodes[i].label = nodes[i].index;
			/*fprintf (fpmpi, "\n\n Node index = %d, Node label = %d, class = %d, deme = %d \n\n", nodes[i].index, nodes[i].label, nodes[i].class, nodes[i].indexOldMigPop);*/
			}


		/* Cheking */
		/*for (i = 0; i < p->numSegNode; i++)	
			{
			s = segments + post(i,p->index,distance);
			if (s->after1 == NULL && s->after2 == NULL)
				{*/
				/*
				fprintf (fpmpi, "\n\nThe fragments %d isn't in the trees\n\n", s->sIndex);
				exit (-1);
				}
			}*/
		
		/*if (indNumRE+1 != countNumTrees)
			{
			fprintf (fpmpi, "\n\nWarning, indNumRE+1 %d != countNumTrees %d \n\n", indNumRE+1, countNumTrees);
			exit (-1);
			}*/
		/*for (i=0; i<nextAvailable; i++)
			{
			r = nodes+i;
			
			fprintf(fpmpi, "\n r->index = %d, r->label = %d, r->time = %lf", r->index, r->label, r->time);
			}*/
			
			
			
			
		}

	free (stud);
	/*free (treeRootInit);*/

	free (S_MRCA);
	free (OnlyAncS_MRCA);
	
	/*free (activeGametes);	*/
	/*free (segments);
	free (nodes);*/
	
	}





/********************* CalcREprobs ********************/
/* 
	Returns an array (k = 2 to numSequences) with the
	probability that the recombination occurs while there
	are k lineages 
*/

/*double *CalcREprobs (int numSequences, int rateType)
{
	int			j, k;
	double		*probs;
	double		probs_sum;
	double		rate, jrate, sum;*/
	
	
	/*probs = (double *) calloc((numSequences-1),(long) sizeof(double));
		if (!probs)
		{
		fprintf (fpmpi, "Could not allocate probs (%lu bytes)\n", (numSequences-1) *(long) sizeof(double));
		exit (-1);
		}
	
	probs_sum = 0;
	for (k=numSequences; k>1; k--)
		{*/
		/* which rate type */
		/*if (rateType == 1)*/	/* short external branches and long internal ones */
			/*rate = k * (k - 1) / 2.0;
		else if (rateType == 2)*/	/* in between */
			/*rate = 1;
		else if (rateType == 3)*/ /* long external branches and short internal ones */
			/*rate = numSequences - k + 2.0;*/

		/* calculate the sum term */
		/*sum = 0;
		for (j=2; j<=numSequences; j++)
			{
			if (rateType == 1)*/	/* short external branches and long internal ones */
				/*jrate = j * (j - 1) / 2.0;
			else if (rateType == 2)*/	/* in between */
				/*jrate = 1;
			else if (rateType == 3)*/ /* long external branches and short internal ones */
				/*jrate = numSequences - j + 2.0;

			sum += j / jrate;
			}*/

		/* probability that the recombination occurs while there are k lineages */
		/*probs[k] = (k / rate) / (sum);
		probs_sum += probs[k];*/
		
		
		/*
		fprintf (fpmpi, "\nprobs with %-2d lineages = %6.4f", k, probs[k]);
		fprintf (fpmpi, " (%6.4f)", probs_sum);
		*/
		/*}

	return probs;
}*/






/*********************************************** Codon Model ********************************************/
/* Functions to codon model. Some of them were taken from code provided by R. Nielsen and Z. Yang. 
  It will be showed */
/********************************************************************************************************/

/********************************** EvolveSequenceOnTree_Codon ***********************************/
/* This function evolves sequences on the coalescent trees by codon model */

void EvolveSequenceOnTree_Codon (long int *seed, double m, 
							double alpha, int numNuc, int indNumRE, int *arrayIndBreakpointsOrd, char *MRCAsequence, int numOmegaCat, int numSites)
	{
	int			i, w, n, a, j, sitePosition, siteCodon, numCo, controlStopCodon, Cbroke, M8option;
	double		varRate, ran, cumFreq[64];
	int			out_C[4], codon[3];
	int			*codonMRCASeq, *arrayIndBreakpointsOrd_C;
	double 		*cumProbCatM1, *cumProbCatM8, *cumProbCat, *hetProb;
	double		betaVar, GammaVarRateOmega, VarOmega1;
	TreeNode	*q;
	/*double		Qij_C[NUMCOD][NUMCOD];*/ /* global variable */

	w = sitePosition = controlStopCodon = a = j = M8option = 0;
	n = siteCodon = 1;
	out_C[0] = out_C[1] = out_C[2] = out_C[3] = -1;
	Cbroke = NO;
	doRepitEvol = NO;	


	codonMRCASeq = (int *)calloc((numSites+1),(long) sizeof(int));
	if (!codonMRCASeq)
		{
		fprintf (fpmpi, "Could not allocate codonMRCASeq (%lu bytes)\n", (numSites+1)  * (long) sizeof(int));		
		exit (1);
		}
	arrayIndBreakpointsOrd_C = (int *)calloc((indNumRE+1),(long) sizeof(int));
	if (!arrayIndBreakpointsOrd_C)
		{
		fprintf (fpmpi, "Could not allocate arrayIndBreakpointsOrd_C (%lu bytes)\n", (indNumRE+1)  * (long) sizeof(int));
		exit (1);
		}
	

	/* making codon breakpoints */
	for (i = 0; i < indNumRE+1; i++)
		arrayIndBreakpointsOrd_C[i] = fabs(arrayIndBreakpointsOrd[i]/3);
	
	if (doMRCAFile == YES) /* MRCA from File */
		{
		for (i = 1; i <= numNuc; i++)		/* taking MRCA sequence */
			{
			sitePosition = n;	/* site position */
			n++;
			if (n == 4)
				n = 1;
		
			if (i == arrayIndBreakpointsOrd[w]) /* recombinations, breakpoints. CAMBIA LA RAIZ DEL ARBOL */
				w++;
		

			w = 0; /* NETRECODON Ahora solo hay un nodo raiz, "GMRCA" base de toda la red */
			

			if (sitePosition == 1)
				out_C[1] = EnterCodonMRCA_File (treeRootInit[w], i, numNuc, MRCAsequence, out_C);
			else if (sitePosition == 2)
				out_C[2] = EnterCodonMRCA_File (treeRootInit[w], i, numNuc, MRCAsequence, out_C);
			else
				{
				out_C[3] = EnterCodonMRCA_File (treeRootInit[w], i, numNuc, MRCAsequence, out_C);
			
				if ((out_C[1] == 3 && out_C[2] == 0 && out_C[3] == 0) || (out_C[1] == 3 && out_C[2] == 0 && out_C[3] == 2) || (out_C[1] == 3 && out_C[2] == 2 && out_C[3] == 0)) /* Checking stop codons */
					{
					fprintf (fpmpi, "\nWarning: There are some STOP CODON in your MRCA sequence:");
					if (out_C[1] == 3 && out_C[2] == 0 && out_C[3] == 0)
						fprintf (fpmpi, " TAA");
					if (out_C[1] == 3 && out_C[2] == 0 && out_C[3] == 2)
						fprintf (fpmpi, " TAG");
					if (out_C[1] == 3 && out_C[2] == 2 && out_C[3] == 0)
						fprintf (fpmpi, " TGA");
					exit(-1);
					}
				if (out_C[0] != -1 || out_C[1] != -1 || out_C[2] != -1 || out_C[3] != -1)	
					{
					matrixC[pos(out_C[0] /*p->label*/,siteCodon,numSites)] = makeCodonFromNuc(out_C[1], out_C[2], out_C[3]);
					
					if (matrixC[pos(out_C[0] /*p->label*/,siteCodon,numSites)] > 60)
						{
						fprintf (fpmpi, "\n stop codon20 \n");
						exit(-1);
						}
					siteCodon++;
					out_C[0] = out_C[1] = out_C[2] = out_C[3] = -1;
					}
				else
					{
					fprintf (fpmpi, "\nWarning: Error from EnterCodonMRCA function");
					exit (-1);
					}
				}
			}
		}
	else /* MRCA from nucleotide frequencies */
		{

		cumFreq[0] = codonTable_frequencies_MRCA(0);
		for (i = 1; i <= 63; i++)
			cumFreq[i] = cumFreq[i-1] + codonTable_frequencies_MRCA(i);

		/*fprintf (fpmpi, "\n cumFreq[0] = %lf", cumFreq[0]);*/
		for (i = 1; i <= 63; i++)
			{
			cumFreq[i] = cumFreq[i]/cumFreq[63];
			/*fprintf (fpmpi, "\n cumFreq[%d] = %lf", i, cumFreq[i]);*/
			}
		
		
		for (i=1; i <= numSites; i++)
			{
			ran = RandomUniform(seed);
			codonMRCASeq[i] = bbin_EnterMRCA(ran, cumFreq); /* codon MRCA sequence */
			/*if (i == 1)
				{
				fprintf (fpmpi, "\n ran = %lf", ran);
				fprintf (fpmpi, "\n codonMRCASeq[%d] = %d", i, codonMRCASeq[i]);
				}*/
			
			/*fprintf (fpmpi, "\n codon[%d] = %d", i, codonMRCASeq[i]);*/  /* active to see the GMRCA codon sequence */
			
			if (codonMRCASeq[i] == 48 || codonMRCASeq[i] == 50 || codonMRCASeq[i] == 56) /* Cheking stop codons */
				{
				fprintf (fpmpi, "\n Warning by stop codons in EvolveSequenceOnTree_Codon");
				exit (-2);
				}
			}
		

		numCo = 1;
		for (i = 1; i <= numNuc; i++)
			{
			/*fprintf (fpmpi, "\n nuc %d, numCo %d", i, numCo);*/			
			
			sitePosition = n;	/* site position */
			n++;
			
			if (n == 4)
				n = 1;
				
			if (sitePosition == 1 && i != 1)
				numCo++;
				
			if (i == arrayIndBreakpointsOrd[w]) /* recombinations, breakpoints */
				w++;					
			w = 0; /* NETRECODON Ahora solo hay un nodo raiz, "GMRCA" base de toda la red */


			number_to_codon_MRCA(codonMRCASeq[numCo], codon); /* de un codon salen 3 nucleotidos para el vector "codon" */


			if (sitePosition == 1)
				{
				out_C[1] = EnterCodonMRCA_Freq (treeRootInit[w], i, sitePosition, numNuc, out_C, codon);
				/*fprintf (fpmpi, "\n out_C[%d] = %d", sitePosition, out_C[sitePosition]);*/
				}
			else if (sitePosition == 2)
				{
				out_C[2] = EnterCodonMRCA_Freq (treeRootInit[w], i, sitePosition, numNuc, out_C, codon);
				/*fprintf (fpmpi, "\n out_C[%d] = %d", sitePosition, out_C[sitePosition]);*/
				}
			else
				{
				out_C[3] = EnterCodonMRCA_Freq (treeRootInit[w], i, sitePosition, numNuc, out_C, codon);
				/*fprintf (fpmpi, "\n out_C[%d] = %d", sitePosition, out_C[sitePosition]);*/	

				if ((out_C[1] == 3 && out_C[2] == 0 && out_C[3] == 0) || (out_C[1] == 3 && out_C[2] == 0 && out_C[3] == 2) || (out_C[1] == 3 && out_C[2] == 2 && out_C[3] == 0))
					{
					fprintf (fpmpi, "\nWarning: There are some STOP CODONS in the MRCA sequence:");
					if (out_C[1] == 3 && out_C[2] == 0 && out_C[3] == 0)
						fprintf (fpmpi, " TAA");
					if (out_C[1] == 3 && out_C[2] == 0 && out_C[3] == 2)
						fprintf (fpmpi, " TAG");
					if (out_C[1] == 3 && out_C[2] == 2 && out_C[3] == 0)
						fprintf (fpmpi, " TGA");						
					exit(-1);
					}
				
				if (out_C[0] != -1 || out_C[1] != -1 || out_C[2] != -1 || out_C[3] != -1)	
					{
					matrixC[pos(out_C[0] /*p->label*/,siteCodon,numSites)] = makeCodonFromNuc(out_C[1], out_C[2], out_C[3]);
					/*fprintf (fpmpi, "\n from nuc Freq, GMRCA, out_C_0 = %d\n", out_C[0]);*/
					
					if (matrixC[pos(out_C[0] /*p->label*/,siteCodon,numSites)] > 60)
						{
						fprintf (fpmpi, "\n stop codon21 \n");
						exit(-1);
						}
					siteCodon++;
					out_C[0] = out_C[1] = out_C[2] = out_C[3] = -1;
					}
				else
					{
					fprintf (fpmpi, "\nWarning: Error from EnterCodonMRCA function");
					exit (-1);
					}
				}
			}
		

		}
	

	free (codonMRCASeq);
	


	/*fprintf (fpmpi, "\n Entra en la recursion, EvolveSequenceOnTree_Codon, numRE = %d, numREbreakCod = %d", numRE, numREbreakCod);*/
	w = 0;
	if (doPrintOmegasPerSitefiles == YES) /* print omegas to file */
		fprintf(fpOmegasPerSitePrint,"Site      OmegaValue\n");
		
	for (i = 1; i <= numSites; i++)	/************* Other sequences */ /* Hace esto para cada codon FOR EACH SITE (CODON)**********/
		{
		/*fprintf (fpmpi, "\n\n +++++++ CODON = %d ++++++++ \n\n", i);*/
		Cbroke = NO;
		





		/******* ALL VARIABLE CODON MODELS PER SITE. PROB, GAMMA, M1, M7 and M8 codon models *******/
		if (doOmegaProb == YES)
			{
			cumProbCat =  (double*) calloc ((numOmegaCat+1), sizeof (double));  
			if (cumProbCat == NULL)
				{
				fprintf (fpmpi, "Could not allocate cumProbCat (%lu bytes)", (numOmegaCat+1) * (long) sizeof (double));
				exit(1);
				}
			
			for (j=1; j<=numOmegaCat; j++)
				cumProbCat[j] = 0;
			cumProbCat[0] = 0;
			for (j=1; j<=numOmegaCat; j++)
				cumProbCat[j] = cumProbCat[j-1] + omegaProb[j];
			ran = RandomUniform(seed);
			ProbCategory = bbinInOmegaCat (ran, cumProbCat, numOmegaCat);
			j = 0;
			free (cumProbCat);
			/*fprintf (stderr, "\n ProbCategory = %d \n", ProbCategory);*/
			if (doPrintOmegasPerSitefiles == YES) /* print omegas to file */
				{
				fprintf(fpOmegasPerSitePrint,"Site %d - Omega %3.5f\n", i, omegaVal[ProbCategory]);
				}
			}

		if (doOmegaRateHetDisc == YES)
			{
			hetProb =  (double*) calloc ((numOmegaCat+1), sizeof (double));  
			if (hetProb == NULL)
				{
				fprintf (fpmpi, "Could not allocate hetProb (%lu bytes)", (numOmegaCat+1) * (long) sizeof (double));
				exit(1);
				}
			cumProbCat =  (double*) calloc ((numOmegaCat+1), sizeof (double));  
			if (cumProbCat == NULL)
				{
				fprintf (fpmpi, "Could not allocate cumProbCat (%lu bytes)", (numOmegaCat+1) * (long) sizeof (double));
				exit(1);
				}			


			for (j=1; j<=numOmegaCat; j++)
				{
				cumProbCat[j] = 0;
				hetProb[j] = 1.0/numOmegaCat;
				}
			cumProbCat[0] = 0;
			for (j=1; j<=numOmegaCat; j++)
				cumProbCat[j] = cumProbCat[j-1] + hetProb[j];
			ran = RandomUniform(seed);
			GammCategory = bbinInOmegaCat (ran, cumProbCat, numOmegaCat);
			j = 0;
			free (hetProb);
			free (cumProbCat);

			/*fprintf (stderr, "\n GammCategory = %d \n", GammCategory);*/
			if (doPrintOmegasPerSitefiles == YES) /* print omegas to file */
				{
				fprintf(fpOmegasPerSitePrint,"Site %d - Omega %3.5f\n", i, omegaValGammaRate[GammCategory]);
				}
			}



		if (doM1 == YES) /* for each site */ /* Important Note: M1_FinalSite_omega = 1, prop0, omega will be 0; if M1_FinalSite_omega = 2, prop1, omega will be 1  */
			{
			cumProbCatM1 =  (double*) calloc ((numOmegaCat+1), sizeof (double));  
			if (cumProbCatM1 == NULL)
				{
				fprintf (fpmpi, "Could not allocate cumProbCatM1 (%lu bytes)", (numOmegaCat+1) * (long) sizeof (double));
				exit(1);
				}
			
			/*fprintf(fpmpi,"\n Categ number = %d \n", numOmegaCat);*/ 


			for (j=1; j<=numOmegaCat; j++)
				cumProbCatM1[j] = 0;
			cumProbCatM1[0] = 0;
			for (j=1; j<=numOmegaCat; j++)
				{
				cumProbCatM1[j] = cumProbCatM1[j-1] + omegaProb[j];
				/*fprintf (stderr, "\n cumProbCatM1[%d] = %lf \n", j, cumProbCatM1[j]);*/
				}
			ran = RandomUniform(seed);
			M1_FinalSite_omega = bbinInOmegaCat (ran, cumProbCatM1, numOmegaCat);
			j = 0;
			/*fprintf (stderr, "\n M1_FinalSite_omega (or category) = %d, ran = %lf \n", M1_FinalSite_omega, ran);*/
			if (doPrintOmegasPerSitefiles == YES) /* print omegas to file */
				{
				if (M1_FinalSite_omega == 1)
					{
					fprintf(fpOmegasPerSitePrint,"Site %d - Omega %3.5f\n", i, M1_omega0);
					}
				if (M1_FinalSite_omega == 2)
					{
					fprintf(fpOmegasPerSitePrint,"Site %d - Omega 1.00000\n", i);
					}
				}
			
			free (cumProbCatM1);
			}
		

		if (doOmegaRateHetCont == YES) /* for each site and each branch */
			{
			GammaVarRateOmega = RndGamma (OmegaRateHet, seed) / OmegaRateHet; 
			VarOmega1 = OmegaInit*GammaVarRateOmega;
			omega = roundit(VarOmega1,5);
			if (doPrintOmegasPerSitefiles == YES) /* print omegas to file */
				{
				fprintf(fpOmegasPerSitePrint,"Site %d - Omega %3.5f\n", i, omega);
				}
			/*fprintf(stderr,"doOmegaRateHetCont - Site %d - Omega %3.5f\n", i, omega);*/

			/*fprintf(stderr, "\nBEFORE Cijk_C[10*NUMCOD*NUMCOD+10*NUMCOD+10] = %d\n", Cijk_C[10*NUMCOD*NUMCOD+10*NUMCOD+10]);
			fprintf(stderr, "\nBEFORE Cijk_C[5*NUMCOD*NUMCOD+5*NUMCOD+5] = %d\n", Cijk_C[5*NUMCOD*NUMCOD+5*NUMCOD+5]);*/

			buildCodonMatrix_Qij_Cijk ();
			/*fprintf(stderr, "\nAFTER Cijk_C[10*NUMCOD*NUMCOD+10*NUMCOD+10] = %d\n", Cijk_C[10*NUMCOD*NUMCOD+10*NUMCOD+10]);
			fprintf(stderr, "\nAFTER Cijk_C[5*NUMCOD*NUMCOD+5*NUMCOD+5] = %d\n", Cijk_C[5*NUMCOD*NUMCOD+5*NUMCOD+5]);*/	
			}


		if (doM7 == YES) /* for each site */  
			{
			/*fprintf (stderr, "\n M7_p_beta = %lf, M7_q_beta = %lf \n", M7_p_beta, M7_q_beta);*/	
			betaVar = RandomBeta (M7_p_beta, M7_q_beta, seed); /* Thanks David Posada! */
			M7_FinalSite_omega = betaVar;

			omega = roundit(M7_FinalSite_omega,5);

			/*omega = M7_FinalSite_omega;
			if (omega < 0.001)
				omega = 0.00;*/	
			/*fprintf (stderr, "\n M7: Site %d, omega = %lf, M7_FinalSite_omega = %lf  \n", i, omega, M7_FinalSite_omega);*/	
			buildCodonMatrix_Qij_Cijk ();

			if (doPrintOmegasPerSitefiles == YES) /* print omegas to file */
				{
				fprintf(fpOmegasPerSitePrint,"Site %d - Omega %3.5f\n", i, omega);
				}
			}


		if (doM8 == YES) /* for each site */  
			{
			cumProbCatM8 =  (double*) calloc ((3), sizeof (double));  
			if (cumProbCatM8 == NULL)
				{
				fprintf (fpmpi, "Could not allocate cumProbCatM8 (%lu bytes)", (3) * (long) sizeof (double));
				exit(1);
				}

			for (j=0; j<=2; j++)
				cumProbCatM8[j] = 0;
			/*fprintf (stderr, "\n M8_P0_beta = %lf, M8_P1_omega = %lf \n", M8_P0_beta, M8_P1_omega);*/			
			for (j=1; j<=2; j++)
				{
				if (j == 1)
					{
					cumProbCatM8[j] = cumProbCatM8[j-1] + M8_P0_beta;
					}
				if (j == 2)
					{
					cumProbCatM8[j] = cumProbCatM8[j-1] + M8_P1_omega;
					}
				/*fprintf (stderr, "\n cumProbCatM8[%d] = %lf \n", j, cumProbCatM8[j]);*/
				}
			
			ran = RandomUniform(seed);
			M8option = bbinInOmegaCat (ran, cumProbCatM8, 2);
			j = 0;			
			free (cumProbCatM8);
			M8_FinalSite_omega = -1.0;
			/*fprintf (stderr, "\n M8option (or category) = %d, ran = %lf \n", M8option, ran);*/

			if (M8option == 1) /* beta distribution */ /* M8_p_beta, M8_q_beta -> M8_omegaP0 */
				{
				/* From David Posada. Thanks David!! */
				/*int		i, numReps;
				long int	seed;
				double	betaVar, p, q, sum, sumsq, variance;
				seed = time(NULL); 			
				p = 0.3;
				q = 0.5;
				sum = 0;
				numReps = 1000;
				fprintf (stderr, "beta variable = ");
				for (i=0; i<numReps; i++)
					{
					betaVar = RandomBeta (p, q, &seed);
					fprintf (stderr, " %6.4f", betaVar);
					sum += betaVar;
					sumsq += pow(betaVar,2);
					}
				variance = (1.0 / (double) (numReps-1)) * (sumsq - pow(sum,2) / (double) numReps);
				fprintf (stderr, "\nmean = %6.4f  var = %6.4f", sum/=numReps, variance);
				return 0;*/

				betaVar = RandomBeta (M8_p_beta, M8_q_beta, seed);
				M8_omegaP0 = betaVar;
				M8_FinalSite_omega = M8_omegaP0;
				/*fprintf (stderr, "\n option 1, omega = %lf   D.beta: p = %lf, q = %lf \n", M8_FinalSite_omega, M8_p_beta, M8_q_beta);*/
				}
			else if (M8option == 2) /* direct omega */
				{
				M8_FinalSite_omega = M8_omegaP1;
				/*fprintf (stderr, "\n option 2, omega = %lf \n", M8_FinalSite_omega);*/
				}
			else
				{
				fprintf (fpmpi, "\n Error in M8option (%d): EvolveSequenceOnTree_Codon \n", M8option);
				exit(1);
				}
			
			if (M8_FinalSite_omega == -1.0) /* Check point */
				{
				fprintf (fpmpi, "\n Error in M8_FinalSite_omega (-1): EvolveSequenceOnTree_Codon \n");
				exit(1);
				}

			/*fprintf (stderr, "\n M8: Site %d, M8_FinalSite_omega = %lf  \n", i, M8_FinalSite_omega);*/
			/*omega = M8_FinalSite_omega;*/
			
			omega = roundit(M8_FinalSite_omega,5);

			/*fprintf (stderr, "\n M8: Site %d, omega = %lf  \n", i, omega);*/
			/*if (omega < 0.001)
				omega = 0.00;*/

			buildCodonMatrix_Qij_Cijk ();

			if (doPrintOmegasPerSitefiles == YES) /* print omegas to file */
				{
				fprintf(fpOmegasPerSitePrint,"Site %d - Omega %3.5f\n", i, omega);
				}
			}
		/*** End ALL VARIABLE CODON MODELS ***/




		/*fprintf (fpmpi, "\n--- COMIENZA RED SOBRE CODON %d ---\n", i);*/
		for (a=0; a<nextAvailable; a++) /*Para ver si este codon se rompe en su evolucion*/
			{
			q = nodes + a;
			//fprintf (fpmpi, "\n\n q->breakCodon = %d \n\n", q->breakCodon);
			if (i == q->breakCodon)
				{
				Cbroke = YES; /* el codon i se rompe */
				/*fprintf (fpmpi, "Este codon se rompe de este nodo, q->index = %d,  q->label = %d \n", q->index, q->label);*/
				}
			}

		if (RandomUniform(seed) < pinv)		
			varRate = 0.0;
		else
			{
			if (doRateHet == YES)
				varRate = RndGamma (alpha, seed) / alpha; 
			else
				varRate = 1; 
			}
		
		if (i == (arrayIndBreakpointsOrd_C[w]+1) && indNumRE != 0)
			w++;
		
		w = 0; /* NETRECODON Ahora solo hay un nodo raiz, "GMRCA" base de toda la red */
		if (doRepitEvol == NO)
			SimulateDataForSite_Codon_RECURSIVE_NET (treeRootInit[w], i, numSites, m, numOmegaCat, varRate, seed, Cbroke);
		}
	/*fprintf (fpmpi, "\n\n SALE DE LA RECURSION \n\n");*/
	


	if (doRepitEvol == NO) 
		free (arrayIndBreakpointsOrd_C);
	else /* Reorganiza si la evolucion no fue valida */
		{
		/*fprintf (fpmpi, "\n\n SALE CON CODONES STOP \n\n");*/
		/*exit (2);*/
		for (a=0; a<nextAvailable; a++) /*reorganiza los nodos como al inicio*/
			{
			q = nodes + a;
			q->passNumber = 0;
			}
		numMU = 0;
		numMU_S = 0;
		numMU_NS = 0;
		numEqual2 = numEqual1 = numDifCodSameAA = numDifCodDifAA = numNonSyn0 = numNonSyn1 = numNonSyn2 = 0;
		}



	}



/* To round decimals in a number */
double roundit(double d,int dig) /* d is a double number, dig is the number of decimals to cut */
	{
	double c;
	long double m = powl(10,dig);
	c = roundl(d*m) / m;
	return c;
	}








/********************************** SimulateDataForSite_Nucleotide_RECURSIVE_NET ***********************************/
/* Simulates the nucleotide substitution process for a given nucleotide */

void SimulateDataForSite_Nucleotide_RECURSIVE_NET (TreeNode *p, int siteNucleotide, int numSites, double m, double varRate, double kappa, long int *seed)
	{
	int			i, k, step, control;
	double		ran, cumProb[4], Pij[4][4];
	double 		a, b;
	TreeNode *q;
	TreeSegment *s;
	
//	TreeNode *q, *r;


	/*
	fprintf (fpmpi, "\n codon = %d", siteCodon);
	*/
	q = NULL;
	a = 0;
	b = 0;
	step = 0;
	k = 0;
	control = 0;
	

	
	if (numRE == 0) /*** NO HAY NINGUNA REC ***/
		{
		if (p != NULL)
			{
			/*fprintf (fpmpi, "\n\n--Pasa por p->label = %d--\n", p->label);*/
			if (p->anc1 != NULL) // todos aquellos nodos que no son el GMRCA
				{
				/*fprintf (fpmpi, "\n COAL p->anc1 != NULL.. p->label = %d\n", p->label);*/		
				if (p->isOutgroup == YES) // llegada al nodo outgroup
					{
					SubstitutionMatrix (Pij, p->length * m, kappa, varRate, p_i);
					}
				else //llegada a cualquier nodo que no es outgroup
					{
					SubstitutionMatrix (Pij, (p->anc1->time - p->time) * m, kappa, varRate, p_i);
					}
			
				/* Introduce Mutacion */
				cumProb[0] = Pij[matrix[pos(p->anc1->label,siteNucleotide,numSites)]][0];
				for (i=1; i<4; i++)
					cumProb[i] = cumProb[i-1] + Pij[matrix[pos(p->anc1->label,siteNucleotide,numSites)]][i];
				ran = RandomUniform(seed);

				/*if (p->label == 1)
					{
					fprintf(stderr,"\n ***** p->label = %d, p->index = %d,  p->time = %lf, (p->anc1->time - p->time) = %lf, sitioNum = %d, deme = %d, node_original = %d \n", p->label, p->index,  p->time, (p->anc1->time - p->time), siteNum, p->indexOldMigPop, p->NetIndex);
					}*/
				if (ran >= 0.0 && ran <= cumProb[0])
					matrix[pos(p->label,siteNucleotide,numSites)] = 0; 
				else if (ran > cumProb[0] && ran <= cumProb[1])
					matrix[pos(p->label,siteNucleotide,numSites)] = 1; 
				else if (ran > cumProb[1] && ran <= cumProb[2])
					matrix[pos(p->label,siteNucleotide,numSites)] = 2; 
				else
					matrix[pos(p->label,siteNucleotide,numSites)] = 3; 
			

				if (matrix[pos(p->label,siteNucleotide,numSites)] != matrix[pos(p->anc1->label,siteNucleotide,numSites)])
					{
					numMU++;	
					/*fprintf(stderr,"\n < MUT (no rec), numMU = %d \n", numMU);*/
					}
			

				/*fprintf (stderr, " matrix[pos(p->label,siteNucleotide,numSites)] = %d, matrix[pos(p->anc1->label,siteNucleotide,numSites)] = %d, numMU = %d \n", matrix[pos(p->label,siteNucleotide,numSites)], matrix[pos(p->anc1->label,siteNucleotide,numSites)], numMU);
				fprintf (stderr, " (p->anc1->time - p->time) = %lf \n", (p->anc1->time - p->time));*/
				/* fin Mutacion */
				}


			/* It crosses the tree */
			SimulateDataForSite_Nucleotide_RECURSIVE_NET (p->left, siteNucleotide, numSites, m, varRate, kappa, seed);
			SimulateDataForSite_Nucleotide_RECURSIVE_NET (p->right, siteNucleotide, numSites, m, varRate, kappa, seed);	
			if (thereisOutgroup == YES)
				SimulateDataForSite_Nucleotide_RECURSIVE_NET (p->outgroup, siteNucleotide, numSites, m, varRate, kappa, seed);	
			}
		}
	else  /*** HAY RECOMBINACIONES ***/ /* .... EL CODON NO SE ROMPE CON RECOMBINACIONES .... */
		{
			
		if (p != NULL)
			{
			/*fprintf (fpmpi, "\n\n--Pasa por p->label = %d, p->index = %d--\n", p->label, p->index);*/

			k = 0;
			for (step = 0; step < p->numSegNode; step++)	/* pasan aquellos nodos contenidos por el NUCLEOTIDE, que es llamado siteNucleotide */			
				{
				s = segments + post(step,p->index,distance);				
				if (s->sStart <= siteNucleotide && s->sEnd >= siteNucleotide)
					k++;
				}

			if (p->anc1 != NULL && p->isOutgroup == YES) /* pasa el outgroup si lo hay */
				k++;



			if (p->anc1 != NULL && k > 0) // todos aquellos nodos que no son el GMRCA
				{
				/*fprintf (fpmpi, "\n Entra en BETWEEN p->anc1 != NULL.. p->label = %d\n", p->label);*/			
				if (p->isOutgroup == YES) // llegada al nodo outgroup
					{
					SubstitutionMatrix (Pij, p->length * m, kappa, varRate, p_i);
					}
				else // llegada a cualquier nodo que no es outgroup
					{
					SubstitutionMatrix (Pij, (p->anc1->time - p->time) * m, kappa, varRate, p_i);
					}
			


				/* Introduce Mutacion */
				if (matrix[pos(p->anc1->label,siteNucleotide,numSites)] > -1) /* va por anc1 */
					{

					cumProb[0] = Pij[matrix[pos(p->anc1->label,siteNucleotide,numSites)]][0];
					for (i=1; i<4; i++)
						cumProb[i] = cumProb[i-1] + Pij[matrix[pos(p->anc1->label,siteNucleotide,numSites)]][i];
					ran = RandomUniform(seed);

					/*if (p->label == 1)
						{
						fprintf(stderr,"\n ***** p->label = %d, p->index = %d,  p->time = %lf, (p->anc1->time - p->time) = %lf, sitioNum = %d, deme = %d, node_original = %d \n", p->label, p->index,  p->time, (p->anc1->time - p->time), siteNum, p->indexOldMigPop, p->NetIndex);
						}*/
					if (ran >= 0.0 && ran <= cumProb[0])
						matrix[pos(p->label,siteNucleotide,numSites)] = 0; 
					else if (ran > cumProb[0] && ran <= cumProb[1])
						matrix[pos(p->label,siteNucleotide,numSites)] = 1; 
					else if (ran > cumProb[1] && ran <= cumProb[2])
						matrix[pos(p->label,siteNucleotide,numSites)] = 2; 
					else
						matrix[pos(p->label,siteNucleotide,numSites)] = 3; 


					if (matrix[pos(p->label,siteNucleotide,numSites)] != matrix[pos(p->anc1->label,siteNucleotide,numSites)])
						{
						numMU++;
						/*fprintf(stderr,"\n < MUT (no rompe codon, por anc1), numMU = %d \n", numMU);*/
						}
			
					/*fprintf (stderr, " matrix[pos(p->label,siteNucleotide,numSites)] = %d, matrix[pos(p->anc1->label,siteNucleotide,numSites)] = %d, numMU = %d \n", matrix[pos(p->label,siteNucleotide,numSites)], matrix[pos(p->anc1->label,siteNucleotide,numSites)], numMU);
					fprintf (stderr, " (p->anc1->time - p->time) = %lf \n", (p->anc1->time - p->time));*/
						


					/* check  section */
					/*if (matrix[pos(p->anc1->label,siteNucleotide,numSites)] < 0) 
						{
						fprintf (stderr, "\n\n Warning in RECURSIVE_NET function (1): p->anc1->label = %d\n", p->anc1->label);
						exit (-1);
						}
					for (step = 0; step < p->anc1->numSegNode; step++)			
						{
						s = segments + post(step,p->anc1->index,distance);				
						if (s->sStart <= siteNucleotide && s->sEnd >= siteNucleotide)
						control++;
						}
					if (control == 0)
						{
						fprintf (stderr, "\n\n Warning in RECURSIVE_NET function (1): control = %d\n", control);
						exit (-1);
						}*/
					}
				else /* va por anc2 */
					{

					cumProb[0] = Pij[matrix[pos(p->anc2->label,siteNucleotide,numSites)]][0];
					for (i=1; i<4; i++)
						cumProb[i] = cumProb[i-1] + Pij[matrix[pos(p->anc2->label,siteNucleotide,numSites)]][i];
					ran = RandomUniform(seed);								
					/*fprintf (fpmpi, "\n ran = %lf ",ran);*/


					if (ran >= 0.0 && ran <= cumProb[0])
						matrix[pos(p->label,siteNucleotide,numSites)] = 0; 
					else if (ran > cumProb[0] && ran <= cumProb[1])
						matrix[pos(p->label,siteNucleotide,numSites)] = 1; 
					else if (ran > cumProb[1] && ran <= cumProb[2])
						matrix[pos(p->label,siteNucleotide,numSites)] = 2; 
					else
						matrix[pos(p->label,siteNucleotide,numSites)] = 3; 


					if (matrix[pos(p->label,siteNucleotide,numSites)] != matrix[pos(p->anc2->label,siteNucleotide,numSites)])
						{
						numMU++;
						/*fprintf(stderr,"\n < MUT (no rompe codon, por anc2), numMU = %d \n", numMU);*/
						}
					


					/*fprintf (stderr, " matrix[pos(p->label,siteNucleotide,numSites)] = %d, matrix[pos(p->anc2->label,siteNucleotide,numSites)] = %d, numMU = %d \n", matrix[pos(p->label,siteNucleotide,numSites)], matrix[pos(p->anc2->label,siteNucleotide,numSites)], numMU);
					fprintf (stderr, " (p->anc2->time - p->time) = %lf \n", (p->anc2->time - p->time));*/
					
					/* check section */
					/*if (matrixC[pos(p->anc2->label,siteNucleotide,numSites)] < 0) 
						{
						fprintf (stderr, "\n\n Warning in RECURSIVE_NET function (2): p->anc2->label = %d\n", p->anc2->label);
						exit (-1);
						}
					for (step = 0; step < p->anc2->numSegNode; step++)			
						{
						s = segments + post(step,p->anc2->index,distance);				
						if (s->sStart <= siteNucleotide && s->sEnd >= siteNucleotide)
							control++;
						}
					if (control == 0) 
						{
						fprintf (stderr, "\n\n Warning in RECURSIVE_NET function (2): control = %d\n", control);
						exit (-1);
						}*/
					}


				/* fin Mutacion */
				}
		
			/* It crosses the tree */
			if (k > 0)
				{
				SimulateDataForSite_Nucleotide_RECURSIVE_NET (p->left, siteNucleotide, numSites, m, varRate, kappa, seed);
				SimulateDataForSite_Nucleotide_RECURSIVE_NET (p->right, siteNucleotide, numSites, m, varRate, kappa, seed);
				}
			if (thereisOutgroup == YES)
				SimulateDataForSite_Nucleotide_RECURSIVE_NET (p->outgroup, siteNucleotide, numSites, m, varRate, kappa, seed);	
			}
				
				
		}
		
	}














/********************************** SimulateDataForSite_Codon_RECURSIVE_NET ***********************************/
/* Simulates the codon substitution process for a given codon */

void SimulateDataForSite_Codon_RECURSIVE_NET (TreeNode *p, int siteCodon, int numSites, double m, int numOmegaCat, double varRate, long int *seed, int Cbroke)
	{
	int			i, j, k, step, NOcontinueNode, doCombineCodons, control;
	double		ran, cumProb[NUMCOD], Pij[NUMCOD][NUMCOD];
	double a,b;
	int nuc1, nuc2, nuc3, InCodon1, InCodon2, brokePosition, outCodon;
	int aminoacid1, aminoacid2;
	TreeNode *q;
	TreeSegment *s;
	
//	TreeNode *q, *r;


	/*
	fprintf (fpmpi, "\n codon = %d", siteCodon);
	*/
	q = NULL;
	a = 0;
	b = 0;
	step = 0;
	k = 0;
	control = 0;
	nuc3 = siteCodon*3;
	nuc2 = siteCodon*3-1;
	nuc1 = siteCodon*3-2;
	NOcontinueNode = 0;
	doCombineCodons = NO;
	InCodon1 = 0;
	InCodon2 = 0;
	brokePosition = 0;

	if (doRepitEvol == NO) /* control by codon stop formed by recombination inside of codons */
		{
		if (numRE == 0) /*** NO HAY NINGUNA REC ***/
			{
			if (p != NULL)
				{
				/*fprintf (fpmpi, "\n\n--Pasa por p->label = %d--\n", p->label);*/
				if (p->anc1 != NULL) // todos aquellos nodos que no son el GMRCA
					{
					/*fprintf (fpmpi, "\n COAL p->anc1 != NULL.. p->label = %d\n", p->label);*/		
					if (p->isOutgroup == YES) // llegada al nodo outgroup
						{
						if (doOmegaCat == YES || doOmegaRateHetDisc == YES)
							CodonModel_Cat (Pij, p->length * m, varRate);
						else /* omega cte o doOmegaRateHetCont == YES */
							CodonModel (Pij, p->length * m, varRate);
						}
					else //llegada a cualquier nodo que no es outgroup
						{
						if (doOmegaCat == YES || doOmegaRateHetDisc == YES)
							CodonModel_Cat (Pij, (p->anc1->time - p->time) * m, varRate);
						else /* omega cte o doOmegaRateHetCont == YES */
							CodonModel (Pij, (p->anc1->time - p->time) * m, varRate);
						}
			
					/* Introduce Mutacion */
					cumProb[0] = Pij[matrixC[pos(p->anc1->label,siteCodon,numSites)]][0];
					j = matrixC[pos(p->anc1->label,siteCodon,numSites)];
			
					/*fprintf (fpmpi, "\n Pij[%d][x]", j);		
					fprintf (fpmpi, "\n matrixC[pos(p->anc1->label,siteCodon,numSites) = %d", matrixC[pos(p->anc1->label,siteCodon,numSites)]);
					fprintf (fpmpi, "\n cumProb[0] = %lf", cumProb[0]);*/
			
					for (i=1; i<NUMCOD; i++)
						{
						cumProb[i] = cumProb[i-1] + Pij[matrixC[pos(p->anc1->label,siteCodon,numSites)]][i];
						/*fprintf (fpmpi, "\n cumProb[%d] = %lf", i, cumProb[i]);*/
						}
								
					ran = RandomUniform(seed);
					matrixC[pos(p->label,siteCodon,numSites)] = bbin(ran, cumProb); /* binary search in the probabilities */
			
					if (matrixC[pos(p->label,siteCodon,numSites)] > 60)  /* check */
						{
						for (i=0; i<NUMCOD; i++)
							{
							fprintf (stderr, "\n");
							a = 0;
							for (j=0; j<NUMCOD; j++)
								{
								fprintf (stderr, "P[%d][%d] = %3.2f ",i,j,Pij[i][j]);
								a = a + Pij[i][j];
								}
							fprintf (stderr, "\n a = %lf \n\n", a);
							}
				
						fprintf (stderr, "\n SimulateDataForSite_Codon: 1. stop codon22 %d\n", matrixC[pos(p->label,siteCodon,numSites)]);
						exit(-1);
						}
			
					if (matrixC[pos(p->label,siteCodon,numSites)] != matrixC[pos(p->anc1->label,siteCodon,numSites)])
						{
						numMU++;
						/*fprintf(stderr,"\n < MUT (no rec b), numMU = %d \n", numMU);*/

						aminoacid1 = codonTable_DnDs(matrixC[pos(p->label,siteCodon,numSites)]);
						aminoacid2 = codonTable_DnDs(matrixC[pos(p->anc1->label,siteCodon,numSites)]);
						if (aminoacid1 == aminoacid2)
							{
							numMU_S++;
							/*fprintf (stderr, "\n Syn mut: %d to %d, numMU_S = %d \n", aminoacid1, aminoacid2, numMU_S);*/
							}
						else
							{
							numMU_NS++;
							/*fprintf (stderr, "\n NonSyn mut: %d to %d, numMU_NS = %d \n", aminoacid1, aminoacid2, numMU_NS);*/
							}
						if ((aminoacid1 == -1) || (aminoacid2 == -1))
							{
							fprintf (stderr, "\n error in type of mutations %d %d\n", aminoacid1, aminoacid2);
							exit(-1);
							}
						aminoacid1 = aminoacid2 = -1;
						}

				/*	fprintf (stderr, " Node %d to %d - matrixC[pos(p->label,siteCodon,numSites)] = %d, matrixC[pos(p->anc1->label,siteCodon,numSites)] = %d, numMU = %d \n", p->index, p->anc1->index, matrixC[pos(p->label,siteCodon,numSites)], matrixC[pos(p->anc1->label,siteCodon,numSites)], numMU);
					fprintf (stderr, " (p->anc1->time - p->time) = %lf \n", (p->anc1->time - p->time)); */
					/* fin Mutacion */
					}
		
				/* It crosses the tree */
				SimulateDataForSite_Codon_RECURSIVE_NET (p->left, siteCodon, numSites, m, numOmegaCat, varRate, seed, Cbroke);
				SimulateDataForSite_Codon_RECURSIVE_NET (p->right, siteCodon, numSites, m, numOmegaCat, varRate, seed, Cbroke);	
				if (thereisOutgroup == YES)
					SimulateDataForSite_Codon_RECURSIVE_NET (p->outgroup, siteCodon, numSites, m, numOmegaCat, varRate, seed, Cbroke);	
				}
			}
		else  /*** HAY RECOMBINACIONES ***/
			{
			if (Cbroke == NO) /* ESTE CODON NO SE ROMPE CON RECOMBINACIONES */
				{
				if (p != NULL)
					{
					/*fprintf (fpmpi, "\n\n--Pasa por p->label = %d, p->index = %d--\n", p->label, p->index);*/

					k = 0;
					for (step = 0; step < p->numSegNode; step++)	/* pasan aquellos nodos contenidos por el codon */			
						{
						s = segments + post(step,p->index,distance);				
						if (s->sStart <= nuc2 && s->sEnd >= nuc2)
							k++;
						}

					if (p->anc1 != NULL && p->isOutgroup == YES) /* pasa el outgroup si lo hay */
						k++;


					if (p->anc1 != NULL && k > 0) // todos aquellos nodos que no son el GMRCA
						{
						/*fprintf (fpmpi, "\n Entra en BETWEEN p->anc1 != NULL.. p->label = %d\n", p->label);*/			

						if (p->isOutgroup == YES) // llegada al nodo outgroup
							{
							if (doOmegaCat == YES || doOmegaRateHetDisc == YES)
								CodonModel_Cat (Pij, p->length * m, varRate);
							else /* omega cte o doOmegaRateHetCont == YES */
								CodonModel (Pij, p->length * m, varRate);
							}
						else // llegada a cualquier nodo que no es outgroup
							{
							if (doOmegaCat == YES || doOmegaRateHetDisc == YES)
								CodonModel_Cat (Pij, (p->anc1->time - p->time) * m, varRate);
							else /* omega cte o doOmegaRateHetCont == YES */
								CodonModel (Pij, (p->anc1->time - p->time) * m, varRate);
							}
			
						/* Introduce Mutacion */
						if (matrixC[pos(p->anc1->label,siteCodon,numSites)] > -1) /* va por anc1 */
							{
							cumProb[0] = Pij[matrixC[pos(p->anc1->label,siteCodon,numSites)]][0];
							j = matrixC[pos(p->anc1->label,siteCodon,numSites)];
			
							/*fprintf (fpmpi, "\n Pij[%d][x]", j);		
							fprintf (fpmpi, "\n matrixC[pos(p->anc1->label,siteCodon,numSites) = %d", matrixC[pos(p->anc1->label,siteCodon,numSites)]);
							fprintf (fpmpi, "\n cumProb[0] = %lf", cumProb[0]);*/
			
							for (i=1; i<NUMCOD; i++)
								{
								cumProb[i] = cumProb[i-1] + Pij[matrixC[pos(p->anc1->label,siteCodon,numSites)]][i];
								/*fprintf (fpmpi, "\n cumProb[%d] = %lf", i, cumProb[i]);*/
								}
								
							ran = RandomUniform(seed);
							/*fprintf (fpmpi, "\n ran = %lf ",ran);*/
							matrixC[pos(p->label,siteCodon,numSites)] = bbin(ran, cumProb); /* binary search in the probabilities */
						
							/*if (matrixC[pos(p->label,siteCodon,numSites)] > 60)*/  /* check */
							/*	{
								for (i=0; i<NUMCOD; i++)
									{
									fprintf (stderr, "\n");
									a = 0;
									for (j=0; j<NUMCOD; j++)
										{
										fprintf (stderr, "P[%d][%d] = %3.2f ",i,j,Pij[i][j]);
										a = a + Pij[i][j];
										}
									fprintf (stderr, "\n a = %lf \n\n", a);
									}
				
								fprintf (stderr, "\n stop codon22 %d\n", matrixC[pos(p->label,siteCodon,numSites)]);
								exit(-1);
								}*/
			
							if (matrixC[pos(p->label,siteCodon,numSites)] != matrixC[pos(p->anc1->label,siteCodon,numSites)])
								{
								numMU++;
								/*fprintf(stderr,"\n < MUT (rec por anc1 b), numMU = %d \n", numMU);*/

								aminoacid1 = codonTable_DnDs(matrixC[pos(p->label,siteCodon,numSites)]);
								aminoacid2 = codonTable_DnDs(matrixC[pos(p->anc1->label,siteCodon,numSites)]);
								if (aminoacid1 == aminoacid2)
									{
									numMU_S++;
									/*fprintf (stderr, "\n Syn mut: %d to %d, numMU_S = %d \n", aminoacid1, aminoacid2, numMU_S);*/
									}
								else
									{
									numMU_NS++;
									/*fprintf (stderr, "\n NonSyn mut: %d to %d, numMU_NS = %d \n", aminoacid1, aminoacid2, numMU_NS);*/
									}
								if ((aminoacid1 == -1) || (aminoacid2 == -1))
									{
									fprintf (stderr, "\n error in type of mutations %d %d\n", aminoacid1, aminoacid2);
									exit(-1);
									}
								aminoacid1 = aminoacid2 = -1;
								}
						/*	fprintf (stderr, " Node %d to %d - matrixC[pos(p->label,siteCodon,numSites)] = %d, matrixC[pos(p->anc1->label,siteCodon,numSites)] = %d, numMU = %d \n",  p->index, p->anc1->index, matrixC[pos(p->label,siteCodon,numSites)], matrixC[pos(p->anc1->label,siteCodon,numSites)], numMU);
							fprintf (stderr, " (p->anc1->time - p->time) = %lf \n", (p->anc1->time - p->time)); */
						
							/* check  section */
							/*if (matrixC[pos(p->anc1->label,siteCodon,numSites)] < 0) 
								{
								fprintf (stderr, "\n\n Warning in RECURSIVE_NET function (1): p->anc1->label = %d\n", p->anc1->label);
								exit (-1);
								}
							for (step = 0; step < p->anc1->numSegNode; step++)			
								{
								s = segments + post(step,p->anc1->index,distance);				
								if (s->sStart <= nuc2 && s->sEnd >= nuc2)
									control++;
								}
							if (control == 0)
								{
								fprintf (stderr, "\n\n Warning in RECURSIVE_NET function (1): control = %d\n", control);
								exit (-1);
								}*/
							}
						else /* va por anc2 */
							{
							cumProb[0] = Pij[matrixC[pos(p->anc2->label,siteCodon,numSites)]][0];
							j = matrixC[pos(p->anc2->label,siteCodon,numSites)];
			
							/*fprintf (fpmpi, "\n Pij[%d][x]", j);		
							fprintf (fpmpi, "\n matrixC[pos(p->anc1->label,siteCodon,numSites) = %d", matrixC[pos(p->anc1->label,siteCodon,numSites)]);
							fprintf (fpmpi, "\n cumProb[0] = %lf", cumProb[0]);*/
			
							for (i=1; i<NUMCOD; i++)
								{
								cumProb[i] = cumProb[i-1] + Pij[matrixC[pos(p->anc2->label,siteCodon,numSites)]][i];
								/*fprintf (fpmpi, "\n cumProb[%d] = %lf", i, cumProb[i]);*/
								}
								
							ran = RandomUniform(seed);
							/*fprintf (fpmpi, "\n ran = %lf ",ran);*/
							matrixC[pos(p->label,siteCodon,numSites)] = bbin(ran, cumProb); /* binary search in the probabilities */
						
							/*if (matrixC[pos(p->label,siteCodon,numSites)] > 60)*/  /* check */
							/*	{
								for (i=0; i<NUMCOD; i++)
									{
									fprintf (stderr, "\n");
									a = 0;
									for (j=0; j<NUMCOD; j++)
										{
										fprintf (stderr, "P[%d][%d] = %3.2f ",i,j,Pij[i][j]);
										a = a + Pij[i][j];
										}
									fprintf (stderr, "\n a = %lf \n\n", a);
									}
				
								fprintf (stderr, "\n stop codon22 %d\n", matrixC[pos(p->label,siteCodon,numSites)]);
								exit(-1);
								}*/
			
							if (matrixC[pos(p->label,siteCodon,numSites)] != matrixC[pos(p->anc2->label,siteCodon,numSites)])
								{
								numMU++;
								/*fprintf(stderr,"\n < MUT (rec por anc2b), numMU = %d \n", numMU);*/


								/* NEW MA */
								aminoacid1 = codonTable_DnDs(matrixC[pos(p->label,siteCodon,numSites)]);
								aminoacid2 = codonTable_DnDs(matrixC[pos(p->anc2->label,siteCodon,numSites)]);
								if (aminoacid1 == aminoacid2)
									{
									numMU_S++;
									/*fprintf (stderr, "\n Syn mut: %d to %d, numMU_S = %d \n", aminoacid1, aminoacid2, numMU_S);*/
									}
								else
									{
									numMU_NS++;
									/*fprintf (stderr, "\n NonSyn mut: %d to %d, numMU_NS = %d \n", aminoacid1, aminoacid2, numMU_NS);*/
									}
								if ((aminoacid1 == -1) || (aminoacid2 == -1))
									{
									fprintf (stderr, "\n error in type of mutations %d %d\n", aminoacid1, aminoacid2);
									exit(-1);
									}
								aminoacid1 = aminoacid2 = -1;
								}
						/*	fprintf (stderr, " Node %d to %d - matrixC[pos(p->label,siteCodon,numSites)] = %d, matrixC[pos(p->anc2->label,siteCodon,numSites)] = %d, numMU = %d \n",  p->index, p->anc1->index, matrixC[pos(p->label,siteCodon,numSites)], matrixC[pos(p->anc2->label,siteCodon,numSites)], numMU);
							fprintf (stderr, " (p->anc2->time - p->time) = %lf \n", (p->anc2->time - p->time)); */
							
							/* check section */
							/*if (matrixC[pos(p->anc2->label,siteCodon,numSites)] < 0) 
								{
								fprintf (stderr, "\n\n Warning in RECURSIVE_NET function (2): p->anc2->label = %d\n", p->anc2->label);
								exit (-1);
								}
							for (step = 0; step < p->anc2->numSegNode; step++)			
								{
								s = segments + post(step,p->anc2->index,distance);				
								if (s->sStart <= nuc2 && s->sEnd >= nuc2)
									control++;
								}
							if (control == 0) 
								{
								fprintf (stderr, "\n\n Warning in RECURSIVE_NET function (2): control = %d\n", control);
								exit (-1);
								}*/
							}


						/* fin Mutacion */
						}
		
					/* It crosses the tree */
					if (k > 0)
						{
						SimulateDataForSite_Codon_RECURSIVE_NET (p->left, siteCodon, numSites, m, numOmegaCat, varRate, seed, Cbroke);
						SimulateDataForSite_Codon_RECURSIVE_NET (p->right, siteCodon, numSites, m, numOmegaCat, varRate, seed, Cbroke);
						}
					if (thereisOutgroup == YES)
						SimulateDataForSite_Codon_RECURSIVE_NET (p->outgroup, siteCodon, numSites, m, numOmegaCat, varRate, seed, Cbroke);	
					}
				}
			else  /* ESTE CODON SI ROMPE CON RECOMBINACIONES */
				{
				NOcontinueNode = 0;
			
				if (p != NULL)
					{
					k = 0;
					/*fprintf (fpmpi, "\n\n-- DENTRO DE CODON ROTO POR REC: Pasa por p->label = %d, p->index = %d--\n", p->label, p->index);*/

					if (p->anc1 != NULL && p->isOutgroup == YES) /* pasa el outgroup si lo hay */
						k++;
					for (step = 0; step < p->numSegNode; step++)	/* pasan aquellos nodos contenidos por el codon */			
						{
						s = segments + post(step,p->index,distance);				
						if ((s->sStart <= nuc1 && s->sEnd >= nuc1) || (s->sStart <= nuc2 && s->sEnd >= nuc2) || (s->sStart <= nuc3 && s->sEnd >= nuc3)) /* se continua por la evolucion de los 3 nuc q contiene */
							{
							k++;
							}
						}
					/*fprintf (fpmpi, "-- nuc1: %d, nuc2: %d, nuc3: %d --\n", nuc1, nuc2, nuc3);*/

					if (p->breakCodon == siteCodon) /* Special recombinant node */
						{
						
						q = p->sib;

						if (p->passNumber == 0)
							NOcontinueNode = 1;
						if (p->passNumber > 0) /* esta llegando por segunda vez, hay que hacer superposici—n de codones mediante corte nuc */
							{
							NOcontinueNode = 0;
							doCombineCodons = YES;
							/*fprintf (fpmpi, " p->passNumber = %d", p->passNumber);*/
							}
					
						/*p->passNumber++;*/ /* MA add: hoy borre esto. july 2009 */
						q->passNumber++;
						
						/*fprintf (fpmpi, "\n** El nodo p->label = %d (p->index = %d) es recombinante y rompe al codon", p->label, p->index);
						fprintf (fpmpi, " Su Sib es p->passNumber = %d, q->label = %d (q->index = %d), p->whereBreakCodon = %d **", p->passNumber, q->label, q->index, p->whereBreakCodon);*/
						}


					if (p->anc1 != NULL && k > 0) // todos aquellos nodos que no son el GMRCA
						{
						/*fprintf (fpmpi, "\n Entra en INSIDE p->anc1 != NULL.. p->label = %d \n", p->label);*/					
						
						if (p->isOutgroup == YES) // llegada al nodo outgroup
							{
							if (doOmegaCat == YES || doOmegaRateHetDisc == YES)
								CodonModel_Cat (Pij, p->length * m, varRate);
							else /* omega cte o doOmegaRateHetCont == YES */
								CodonModel (Pij, p->length * m, varRate);
							}
						else // llegada a cualquier nodo que no es outgroup
							{
							if (doOmegaCat == YES || doOmegaRateHetDisc == YES)
								CodonModel_Cat (Pij, (p->anc1->time - p->time) * m, varRate);
							else /* omega cte o doOmegaRateHetCont == YES */
								CodonModel (Pij, (p->anc1->time - p->time) * m, varRate);
							}
			
						/* Introduce Mutacion */
						if (matrixC[pos(p->anc1->label,siteCodon,numSites)] > -1) /* va por anc1 */
							{
							cumProb[0] = Pij[matrixC[pos(p->anc1->label,siteCodon,numSites)]][0];
							j = matrixC[pos(p->anc1->label,siteCodon,numSites)];
							/*fprintf (fpmpi, "\n Pij[%d][x]", j);		
							fprintf (fpmpi, "\n matrixC[pos(p->anc1->label,siteCodon,numSites) = %d", matrixC[pos(p->anc1->label,siteCodon,numSites)]);
							fprintf (fpmpi, "\n cumProb[0] = %lf", cumProb[0]);*/
							for (i=1; i<NUMCOD; i++)
								{
								cumProb[i] = cumProb[i-1] + Pij[matrixC[pos(p->anc1->label,siteCodon,numSites)]][i];
								/*fprintf (fpmpi, "\n cumProb[%d] = %lf", i, cumProb[i]);*/
								}
								
							ran = RandomUniform(seed);
							matrixC[pos(p->label,siteCodon,numSites)] = bbin(ran, cumProb); /* binary search in the probabilities */
							/*fprintf (fpmpi, "\n ran = %lf ",ran);*/
					
							/*if (matrixC[pos(p->label,siteCodon,numSites)] > 60)*/  /* check */
							/*	{
								for (i=0; i<NUMCOD; i++)
									{
									fprintf (stderr, "\n");
									a = 0;
									for (j=0; j<NUMCOD; j++)
										{
										fprintf (stderr, "P[%d][%d] = %3.2f ",i,j,Pij[i][j]);
										a = a + Pij[i][j];
										}
									fprintf (stderr, "\n a = %lf \n\n", a);
									}
								fprintf (stderr, "\n stop codon22 %d\n", matrixC[pos(p->label,siteCodon,numSites)]);
								exit(-1);
								}*/
			
							if (matrixC[pos(p->label,siteCodon,numSites)] != matrixC[pos(p->anc1->label,siteCodon,numSites)])
								{
								numMU++;
								/*fprintf(stderr,"\n < MUT (rec SI por anc1), numMU = %d \n", numMU);*/

								aminoacid1 = codonTable_DnDs(matrixC[pos(p->label,siteCodon,numSites)]);
								aminoacid2 = codonTable_DnDs(matrixC[pos(p->anc1->label,siteCodon,numSites)]);
								if (aminoacid1 == aminoacid2)
									{
									numMU_S++;
									/*fprintf (stderr, "\n Syn mut: %d to %d, numMU_S = %d \n", aminoacid1, aminoacid2, numMU_S);*/
									}
								else
									{
									numMU_NS++;
									/*fprintf (stderr, "\n NonSyn mut: %d to %d, numMU_NS = %d \n", aminoacid1, aminoacid2, numMU_NS);*/
									}
								if ((aminoacid1 == -1) || (aminoacid2 == -1))
									{
									fprintf (stderr, "\n error in type of mutations %d %d\n", aminoacid1, aminoacid2);
									exit(-1);
									}
								aminoacid1 = aminoacid2 = -1;
								}
						/*	fprintf (stderr, " Node %d to %d - matrixC[pos(p->label,siteCodon,numSites)] = %d, matrixC[pos(p->anc1->label,siteCodon,numSites)] = %d, numMU = %d \n",  p->index, p->anc1->index, matrixC[pos(p->label,siteCodon,numSites)], matrixC[pos(p->anc1->label,siteCodon,numSites)], numMU);
							fprintf (stderr, " (p->anc1->time - p->time) = %lf \n", (p->anc1->time - p->time)); */

							/* check section */
							if (matrixC[pos(p->anc1->label,siteCodon,numSites)] < 0) 
								{
								fprintf (stderr, "\n\n Warning in RECURSIVE_NET function (3): p->anc1->label = %d\n", p->anc1->label);
								exit (-1);
								}
							for (step = 0; step < p->anc1->numSegNode; step++)			
								{
								s = segments + post(step,p->anc1->index,distance);				
								if ((s->sStart <= nuc1 && s->sEnd >= nuc1) || (s->sStart <= nuc2 && s->sEnd >= nuc2) || (s->sStart <= nuc3 && s->sEnd >= nuc3))
									control++;
								}
							if (control == 0)
								{
								fprintf (stderr, "\n\n Warning in RECURSIVE_NET function (3): control = %d\n", control);
								exit (-1);
								}
							}
						else /* va por anc2 */
							{
							cumProb[0] = Pij[matrixC[pos(p->anc2->label,siteCodon,numSites)]][0];
							j = matrixC[pos(p->anc2->label,siteCodon,numSites)];
							/*fprintf (fpmpi, "\n Pij[%d][x]", j);		
							fprintf (fpmpi, "\n matrixC[pos(p->anc1->label,siteCodon,numSites) = %d", matrixC[pos(p->anc1->label,siteCodon,numSites)]);
							fprintf (fpmpi, "\n cumProb[0] = %lf", cumProb[0]);*/
							for (i=1; i<NUMCOD; i++)
								{
								cumProb[i] = cumProb[i-1] + Pij[matrixC[pos(p->anc2->label,siteCodon,numSites)]][i];
								/*fprintf (fpmpi, "\n cumProb[%d] = %lf", i, cumProb[i]);*/
								}
								
							ran = RandomUniform(seed);
							matrixC[pos(p->label,siteCodon,numSites)] = bbin(ran, cumProb); /* binary search in the probabilities */
							/*fprintf (fpmpi, "\n ran = %lf ",ran);*/
					
							/*if (matrixC[pos(p->label,siteCodon,numSites)] > 60)*/  /* check */
								/*{
								for (i=0; i<NUMCOD; i++)
									{
									fprintf (stderr, "\n");
									a = 0;
									for (j=0; j<NUMCOD; j++)
										{
										fprintf (stderr, "P[%d][%d] = %3.2f ",i,j,Pij[i][j]);
										a = a + Pij[i][j];
										}
									fprintf (stderr, "\n a = %lf \n\n", a);
									}
								fprintf (stderr, "\n stop codon22 %d\n", matrixC[pos(p->label,siteCodon,numSites)]);
								exit(-1);
								}*/
			
							if (matrixC[pos(p->label,siteCodon,numSites)] != matrixC[pos(p->anc2->label,siteCodon,numSites)])
								{
								numMU++;
								/*fprintf(stderr,"\n < MUT (rec SI por anc2), numMU = %d \n", numMU);*/

								/* NEW MA */
								aminoacid1 = codonTable_DnDs(matrixC[pos(p->label,siteCodon,numSites)]);
								aminoacid2 = codonTable_DnDs(matrixC[pos(p->anc2->label,siteCodon,numSites)]);
								if (aminoacid1 == aminoacid2)
									{
									numMU_S++;
									/*fprintf (stderr, "\n Syn mut: %d to %d, numMU_S = %d \n", aminoacid1, aminoacid2, numMU_S);*/
									}
								else
									{
									numMU_NS++;
									/*fprintf (stderr, "\n NonSyn mut: %d to %d, numMU_NS = %d \n", aminoacid1, aminoacid2, numMU_NS);*/
									}
								if ((aminoacid1 == -1) || (aminoacid2 == -1))
									{
									fprintf (stderr, "\n error in type of mutations %d %d\n", aminoacid1, aminoacid2);
									exit(-1);
									}
								aminoacid1 = aminoacid2 = -1;
								}
						/*	fprintf (stderr, " Node %d to %d - matrixC[pos(p->label,siteCodon,numSites)] = %d, matrixC[pos(p->anc2->label,siteCodon,numSites)] = %d, numMU = %d \n",  p->index, p->anc1->index, matrixC[pos(p->label,siteCodon,numSites)], matrixC[pos(p->anc2->label,siteCodon,numSites)], numMU);
							fprintf (stderr, " (p->anc2->time - p->time) = %lf \n", (p->anc2->time - p->time)); */
							
							/* check section */
							if (matrixC[pos(p->anc2->label,siteCodon,numSites)] < 0) 
								{
								fprintf (stderr, "\n\n Warning in RECURSIVE_NET function (4): p->anc2->label = %d\n", p->anc2->label);
								exit (-1);
								}
							for (step = 0; step < p->anc2->numSegNode; step++)			
								{
								s = segments + post(step,p->anc2->index,distance);				
								if ((s->sStart <= nuc1 && s->sEnd >= nuc1) || (s->sStart <= nuc2 && s->sEnd >= nuc2) || (s->sStart <= nuc3 && s->sEnd >= nuc3))
									control++;
								}
							if (control == 0)
								{
								fprintf (stderr, "\n\n Warning in RECURSIVE_NET function (4): control = %d\n", control);
								exit (-1);
								}
							}
						/* fin Mutacion */

						/* Combining the codons of these recombinant nodes */
						if (doCombineCodons == YES) 
							{
							/*fprintf (fpmpi, "\n In doCombineCodons --\n");*/

							if (p->whereBreakCodon == 3)
								{
								InCodon2 = matrixC[pos(p->label,siteCodon,numSites)];
								InCodon1 = matrixC[pos(q->label,siteCodon,numSites)];
								brokePosition = q->whereBreakCodon;
								/*fprintf (stderr, "\n Case P: brokePosition = %d, InCodon1 = %d (q->label = %d, q->index = %d), InCodon2 = %d\n", brokePosition, InCodon1, q->label, q->index, InCodon2);*/
								}
							else if (q->whereBreakCodon == 3)
								{
								InCodon2 = matrixC[pos(q->label,siteCodon,numSites)];
								InCodon1 = matrixC[pos(p->label,siteCodon,numSites)];
								brokePosition = p->whereBreakCodon;
								/*fprintf (stderr, "\n Case q: brokePosition = %d, InCodon1 = %d, InCodon2 = %d\n", brokePosition, InCodon1, InCodon2);*/
								}
							else
								{
								fprintf (stderr, "\n\n Warning in SimulateDataForSite_Codon_RECURSIVE_NET \n");
								exit(-1);
								}
							outCodon = CombineTwoCodons (InCodon1, InCodon2, brokePosition); /*function combine y paste*/
							if (noisy > 3)
								{
								/*fprintf (stderr, " InCodon1 = %d, InCodon2 = %d, brokePosition = %d", InCodon1, InCodon2, brokePosition);*/
								fprintf (stderr, " Combining two codon by recombination: Codon1 = %d, Codon2 = %d", InCodon1, InCodon2);
								fprintf (stderr, " to the new Codon = %d \n", outCodon);
								}

							matrixC[pos(p->label,siteCodon,numSites)] = outCodon;
							matrixC[pos(q->label,siteCodon,numSites)] = outCodon;
							/*fprintf (stderr, " NEW FOR!!: p->label = %d, q->label = %d, matrixC[pos(p->label,siteCodon,numSites)] = %d, numMU = %d \n\n", p->label, q->label, matrixC[pos(p->label,siteCodon,numSites)], numMU);*/
							}
						/* fin Mutacion */
						}
		
					/* It crosses the tree */
					if (doRepitEvol == NO)
						{
						if (NOcontinueNode == 0 && k > 0)
							{
							SimulateDataForSite_Codon_RECURSIVE_NET (p->left, siteCodon, numSites, m, numOmegaCat, varRate, seed, Cbroke);
							SimulateDataForSite_Codon_RECURSIVE_NET (p->right, siteCodon, numSites, m, numOmegaCat, varRate, seed, Cbroke);
							}
						else
							NOcontinueNode = 0;	
						if (thereisOutgroup == YES)
							SimulateDataForSite_Codon_RECURSIVE_NET (p->outgroup, siteCodon, numSites, m, numOmegaCat, varRate, seed, Cbroke);	
						}
					}
				}
			}
		}
	}




/********************************** CombineTwoCodons ***********************************/
/* It combines two codons to make a new codon based in a broke position  */

int		CombineTwoCodons (int InCodon1, int InCodon2, int brokePosition)
	{
	int outCodon, nucNumber1[3], nucNumber2[3], nucOutNumber[3];
	char  nuc1[3], nuc2[3]/*, out[3]*/;
	int i, j, k, a, b, c;

	outCodon = j = k = a = b = c = 0;
	for (i=0; i<=2; i++)
		{
		nucNumber1[i] = -1;
		nucNumber2[i] = -1;
		nucOutNumber[i] = -1;
		}

	if (InCodon1 == InCodon2)
		outCodon = InCodon1;
	else
		{
		number_to_codon(InCodon1, nuc1);
		number_to_codon(InCodon2, nuc2);

		for (i=0; i<=2; i++)
			{
			if (nuc1[i] == 'A')
				nucNumber1[i] = 0;
			else if (nuc1[i] == 'C')
				nucNumber1[i] = 1;
			else if (nuc1[i] == 'G')
				nucNumber1[i] = 2;
			else if (nuc1[i] == 'T')
				nucNumber1[i] = 3;
			else {			
				fprintf (fpmpi, "error in codon of CombineTwoCodons (1)\n"); 			
				exit(-1);
				}
			}
		for (i=0; i<=2; i++)
			{
			if (nuc2[i] == 'A')
				nucNumber2[i] = 0;
			else if (nuc2[i] == 'C')
				nucNumber2[i] = 1;
			else if (nuc2[i] == 'G')
				nucNumber2[i] = 2;
			else if (nuc2[i] == 'T')
				nucNumber2[i] = 3;
			else {			
				fprintf (fpmpi, "error in codon of CombineTwoCodons (2)\n"); 			
				exit(-1);
				}
			}

		/*for (i = 0; i<=2; i++)
			{
			fprintf (stderr, " nucNumber1[%d] = %d, nucNumber2[%d] = %d \n", i, nucNumber1[i], i, nucNumber2[i]);			
			}*/		

		/*j = makeCodonFromNuc (nucNumber1[0], nucNumber1[1], nucNumber1[2]);
		k = makeCodonFromNuc (nucNumber2[0], nucNumber2[1], nucNumber2[2]);
		fprintf (stderr, "\n codon1 = %d, codon2 = %d \n", j, k);*/
		/* Codon crossing */
		if (brokePosition == 1)
			{
			nucOutNumber[0] = nucNumber1[0];
			nucOutNumber[1] = nucNumber2[1];
			nucOutNumber[2] = nucNumber2[2];
			}
		else if (brokePosition == 2)
			{
			nucOutNumber[0] = nucNumber1[0];
			nucOutNumber[1] = nucNumber1[1];
			nucOutNumber[2] = nucNumber2[2];
			}
		else 
			{			
			fprintf (fpmpi, "error in codon of CombineTwoCodons - Across codons\n"); 			
			exit(-1);
			}
		/*for (i = 0; i<=2; i++)
			{
			fprintf (stderr, " nucOutNumber[%d] = %d\n", i, nucOutNumber[i]);			
			}*/
		/*nucOutNumber[0] = 3;
		nucOutNumber[1] = 0;
		nucOutNumber[1] = 0;*/
		
		/* STOP CODONS TAA TAG TGA */
		if (nucOutNumber[0] == 3)
			if ((nucOutNumber[1] == 0 && nucOutNumber[2] == 0) || (nucOutNumber[1] == 0 && nucOutNumber[2] == 2) || (nucOutNumber[1] == 2 && nucOutNumber[2] == 0))
				{
				if (noisy > 1)
					fprintf (fpmpi, "\nStop codon created from the recombination inside of the codon \nRestarting the evolution model from the GMRCA .. \n");
				/*exit(-1);*/
				doRepitEvol = YES; /* Repit the evolution of the sequences */
				numStopCodonREC++; 			
				}


		if (doRepitEvol == NO)
			outCodon = makeCodonFromNuc (nucOutNumber[0], nucOutNumber[1], nucOutNumber[2]); /* Warning and END if the new codon is STOP */
		else
			outCodon = -1;
		
		
		/*number_to_codon(outCodon, out);
		for (i=0; i<=2; i++)
			{
			if (out[i] == 'A')
				{
				nucOutNumber[i] = 0;
				fprintf (fpmpi, " nucOutNumber[%d] = %d\n", i, nucOutNumber[i]); 
				}
			else if (out[i] == 'C')
				{
				nucOutNumber[i] = 1;
				fprintf (fpmpi, " nucOutNumber[%d] = %d\n", i, nucOutNumber[i]);
				}
			else if (out[i] == 'G')
				{
				nucOutNumber[i] = 2;
				fprintf (fpmpi, " nucOutNumber[%d] = %d\n", i, nucOutNumber[i]);
				}
			else if (out[i] == 'T')
				{
				nucOutNumber[i] = 3;
				fprintf (fpmpi, " nucOutNumber[%d] = %d\n", i, nucOutNumber[i]);
				}
			else {			
				fprintf (fpmpi, "error in codon of CombineTwoCodons (2)\n"); 			
				exit(-1);
				}
			}*/

		}

	
	/* Results */
	/*if (InCodon1 == InCodon2) 
		{
		numEqual2++; 
		}
	if ((InCodon1 != InCodon2) && (outCodon > -1))
		{
		if ((InCodon1 == outCodon) || (InCodon2 == outCodon))
			numEqual1++;
		if ((InCodon1 != outCodon) && (InCodon2 != outCodon))
			{
			a = codonTable_DnDs(InCodon1);
			b = codonTable_DnDs(InCodon2);
			c = codonTable_DnDs(outCodon);

			if ((c == a) || (c == b))
				numDifCodSameAA++;
			if ((c != a) && (c != b))
				numDifCodDifAA++;
			}
		}*/
	/*numEqual2, numEqual1, numDifCodSameAA, numDifCodDifAA;*/
	if (outCodon > -1)
		{
		a = codonTable_DnDs(InCodon1);
		b = codonTable_DnDs(InCodon2);
		c = codonTable_DnDs(outCodon);
		if ((a == c) && (b == c))
			numNonSyn0++;
		if ((c != a) && (c != b))
			numNonSyn2++;
		if ((c != a) && (c == b))
			numNonSyn1++;
		if ((c == a) && (c != b))
			numNonSyn1++;
		}


	return (outCodon);
	}












/********************************** SimulateDataForSite_Codon ***********************************/
/* Simulates the nucleotide substitution process for a given codon */

/*void SimulateDataForSite_Codon (TreeNode *p, int siteCodon, int numSites, double m, int numOmegaCat, double varRate, long int *seed)
	{
	int			i, j;
	double		ran, cumProb[NUMCOD], Pij[NUMCOD][NUMCOD];
	double a,b;
	int aminoacid1, aminoacid2;
	
	
	a = 0;
	b = 0;
	
	if (p != NULL)
		{
		if (p->anc1 != NULL) 
			{			
			if (p->isOutgroup == YES)
				{
				if (doOmegaCat == YES || doOmegaRateHetDisc == YES)
					CodonModel_Cat (Pij, p->length * m, varRate, numOmegaCat, seed);
				else 																*/			/* omega cte o doOmegaRateHetCont == YES */
	/*				CodonModel (Pij, p->length * m, varRate, seed);
				}
			else
				{
				if (doOmegaCat == YES || doOmegaRateHetDisc == YES)
					CodonModel_Cat (Pij, (p->anc1->time - p->time) * m, varRate, numOmegaCat, seed);
				else 																*/			/* omega cte o doOmegaRateHetCont == YES */
	/*				CodonModel (Pij, (p->anc1->time - p->time) * m, varRate, seed);
				}
			
			cumProb[0] = Pij[matrixC[pos(p->anc1->label,siteCodon,numSites)]][0];
			j = matrixC[pos(p->anc1->label,siteCodon,numSites)]; */
			
			/*fprintf (fpmpi, "\n Pij[%d][x]", j);		
			fprintf (fpmpi, "\n matrixC[pos(p->anc1->label,siteCodon,numSites) = %d", matrixC[pos(p->anc1->label,siteCodon,numSites)]);
			fprintf (fpmpi, "\n cumProb[0] = %lf", cumProb[0]);
			*/
			
	/*		for (i=1; i<NUMCOD; i++)
				{
				cumProb[i] = cumProb[i-1] + Pij[matrixC[pos(p->anc1->label,siteCodon,numSites)]][i];	*/
				/*
				fprintf (fpmpi, "\n cumProb[%d] = %lf", i, cumProb[i]);
				*/
	/*			}
								
			ran = RandomUniform(seed);
			matrixC[pos(p->label,siteCodon,numSites)] = bbin(ran, cumProb); 	*/								/* binary search in the probabilities */
			
	/*		if (matrixC[pos(p->label,siteCodon,numSites)] > 60)  				*/								/* check */
	/*			{
				for (i=0; i<NUMCOD; i++)
					{
					fprintf (stderr, "\n");
					a = 0;
					for (j=0; j<NUMCOD; j++)
						{
						fprintf (stderr, "P[%d][%d] = %3.2f ",i,j,Pij[i][j]);
						a = a + Pij[i][j];
						}
					fprintf (stderr, "\n a = %lf \n\n", a);
					}
				
				fprintf (stderr, "\n SimulateDataForSite_Codon: 2. stop codon22 %d\n", matrixC[pos(p->label,siteCodon,numSites)]);
				exit(-1);
				}
			
			if (matrixC[pos(p->label,siteCodon,numSites)] != matrixC[pos(p->anc1->label,siteCodon,numSites)])
				{
				numMU++;
				aminoacid1 = codonTable_DnDs(matrixC[pos(p->label,siteCodon,numSites)]);
				aminoacid2 = codonTable_DnDs(matrixC[pos(p->anc1->label,siteCodon,numSites)]);
				if (aminoacid1 == aminoacid2)
					numMU_S++;
				else
					numMU_NS++;
				if ((aminoacid1 == -1) || (aminoacid2 == -1))
					{
					fprintf (stderr, "\n error in type of mutations %d %d\n", aminoacid1, aminoacid2);
					exit(-1);
					}
				aminoacid1 = aminoacid2 = -1;
				}
				
			}	*/
		
		/* It crosses the tree */
	/*	SimulateDataForSite_Codon (p->left, siteCodon, numSites, m, numOmegaCat, varRate, seed);
		SimulateDataForSite_Codon (p->right, siteCodon, numSites, m, numOmegaCat, varRate, seed);	
		if (thereisOutgroup == YES)
			SimulateDataForSite_Codon (p->outgroup, siteCodon, numSites, m, numOmegaCat, varRate, seed);	
		}
	}	*/




/***************** bbin_EnterMRCA *****************/
/* binary search of the codon MRCA sequence */
int bbin_EnterMRCA (double dat, double *v)
	{
     int init,end,middle;
     init = 0;
     end = 63;
    
	 if (dat >= 0 && dat <= v[0])
		return 0;
	 
	 while (init <= end) 
		{
		middle = (init+end)/2;
         
		if (dat > v[middle-1] && dat <= v[middle])
			return (middle);
		else if (dat > v[middle])
			init = middle+1;
		else 
			end = middle-1;
		}
	
	fprintf (fpmpi, "\n Warning in bbin_EnterMRCA function");
     exit (-1);
	 return -1;
	}




/***************** bbin *****************/
/* binary search of the probabilities */
int bbin (double dat, double *v)
	{
     int init,end,middle;
     init = 0;
     end = NUMCOD;
    
	 if (dat >= 0 && dat <= v[0])
		return 0;
	 
	 while (init <= end) 
		{
		middle = (init+end)/2;
         
		if (dat > v[middle-1] && dat <= v[middle])
			return (middle);
		else if (dat > v[middle])
			init = middle+1;
		else 
			end = middle-1;
		}
	
	fprintf (fpmpi, "\n Warning in bbin function");	
     exit (-1);
	 return -1;
	}


/***************** bbinDemes *****************/
/* binary search in the probabilities with demes*/
int bbinDemes (double dat, double *v, int n)
	{
	int init,end,middle;
	 
	if (dat >= 0 && dat <= v[1])
		return (1); /* first population */
	 
     init = 1;
     end = n;
    
	 while (init <= end) 
		{
		middle = (init+end)/2;
         
		if (dat > v[middle-1] && dat <= v[middle])
			return (middle);
		else if (dat > v[middle])
			init = middle+1;
		else 
			end = middle-1;
		}
	
	fprintf (fpmpi, "\n Warning in bbinDemes function");
     exit (-1);
	 return -1;
	}




/***************** bbinInOmegaCat *****************/
/* binary search in the probabilities with demes*/
int bbinInOmegaCat (double dat, double *v, int n)
	{
	int init,end,middle;
	 
	if (dat >= 0 && dat <= v[1])
		return (1); /* first population */
	 
     init = 1;
     end = n;
    
	 while (init <= end) 
		{
		middle = (init+end)/2;
         
		if (dat > v[middle-1] && dat <= v[middle])
			return (middle);
		else if (dat > v[middle])
			init = middle+1;
		else 
			end = middle-1;
		}
	
	fprintf (fpmpi, "\n Warning in bbinInOmegaCat function");
     exit (-1);
	 return -1;
	}







/********************************** EnterCodonMRCA_File ***********************************/
/* Enter the MRCA sequence (from input file or by the nucleotide frequencies) for codon model */

int		EnterCodonMRCA_File (TreeNode *p, int siteNum, int numNuc, char *MRCAsequence, int out_C[4])
	{
	int			j;
	j = -1;
	
	if (p != NULL)
		{
		if (p->anc1 == NULL) /* root */		
			{
			if (doMRCAFile == YES)
				matrix[pos(p->label,siteNum,numNuc)] = WhichNucNumber(MRCAsequence[siteNum-1]);
			else
				{
				fprintf (fpmpi, "\n Warning in EnterCodonMRCA_File");				
				exit (-1);
				}
			out_C[0] = p->label;
			j = matrix[pos(p->label,siteNum,numNuc)];
			}
		}
	return (j);
	}


/********************************** EnterCodonMRCA_Freq ***********************************/
/* Enter the MRCA sequence (from input file or by the nucleotide frequencies) for codon model */

int		EnterCodonMRCA_Freq (TreeNode *p, int siteNum, int sitePosition, int numNuc, int out_C[4], int codon[3])
	{
	int			j;
	j = -1;


	if (p != NULL)
		{
		if (p->anc1 == NULL) /* root */		
			{
			if (doMRCAFile == YES)
				{
				fprintf (fpmpi, "\n Warning in EnterCodonMRCA_Freq");				
				exit (-1);
				}
			else
				{

				/*fprintf (fpmpi, "\nMRCA p->index = %d", p->index);	
				if (p->label == 1)
					{
					fprintf(stderr,"\n ***** p->label = %d, p->index = %d, p->time = %lf, sitioNum = %d \n", p->label, p->index,  p->time, siteNum);
					}*/


				if (sitePosition == 1)
					matrix[pos(p->label,siteNum,numNuc)] = codon[0];
				else if (sitePosition == 2)
					matrix[pos(p->label,siteNum,numNuc)] = codon[1];
				else if (sitePosition == 3)
					matrix[pos(p->label,siteNum,numNuc)] = codon[2];
				else
					{
					fprintf (fpmpi, "\nWarning in generates de MRCA sequence from nucleotide frequencies");					
					exit(-1);
					}
				}
			out_C[0] = p->label;
			j = matrix[pos(p->label,siteNum,numNuc)];
			}
		}
	return (j);
	}





/********************* buildCodonMatrix_Qij_Cijk **********************/
/* It builds the Qij and Cijk matrix to codon Model */
/* This fuction was verified perfectly by Nielsen codon model Code */
static void	buildCodonMatrix_Qij_Cijk ()
	{
	int i, j, m, w;
	double k[NUMCOD];
	
	/*double Qij_C[NUMCOD][NUMCOD], double omega, double titv, double p_i[4], double p_i_codon[12]*/ /* Global Variables */	
	/* rows = leave codons, columns = arrive codons */
		
	/* initiallity */
	m = 0;
	for (i = 0; i < NUMCOD; i++)
		{
		k[i] = 0.0;
		for (j = 0; j < NUMCOD; j++)
			Qij_C[i][j] = 0.0;
		}
	for (w=0;w<NUMCOD;w++)
		Root_C[w] = 0.0;
	for (w=0;w<NUMCOD*NUMCOD*NUMCOD*NUMCOD;w++)
		Cijk_C[w] = 0.0;

	for (i = 0; i < NUMCOD; i++)
		{
		for (j = 0; j < NUMCOD; j++)
			{
			Qij_CC[m] = 0;
			/*fprintf (stderr, "\nQij_CC[%d] = %3.2f ", m, Qij_CC[m]);*/
			m++;
			}
		}
	m = i = j = 0;
	/*fprintf(stderr,"\n \n buildCodonMatrix_Qij_Cijk: omega = %lf\n\n\n", omega);*/	


	/* Qij matrix building */
	for (i = 0; i < NUMCOD; i++)
		{
		for (j = 0; j < NUMCOD; j++)
			{
			if (numdif_codon(i,j) == 1)	/* only when there is 1 change */ 
				{
				Qij_C[i][j] = codonTable_frequencies(j) /*0.25*/; 	/* frequencies of arrive codon. Main diagonal (0 substitutions) and more than 1 substitution = 0 */
																	/* from here single those of 1 single change are affected, more than 1 change = 0 (something*0 = 0) */
				if (codonTable_DnDs(i) != codonTable_DnDs(j))				
					Qij_C[i][j] = Qij_C[i][j]*omega;							/* nonsynonymous substitutions */
				
				/*fprintf(stderr,"\n \n omega = %lf", omega);*/
				
				if (doCodon_HKY == YES)										/* Codon Model HKY */
					{
					if (codon_tr_tv(i, j) == 0) /* only transitions */
						Qij_C[i][j] = Qij_C[i][j]*/*titv*/kappa;
					}
				
				if (doCodon_GTR == YES)										/* Codon Model GTR */
					Qij_C[i][j] = Qij_C[i][j]*codon_Rmat(i,j); 
				
				if (doCodon_NGTR == YES)									/* Codon Model non variable GTR */
					{
					Qij_C[i][j] = Qij_C[i][j]*codon_NRmat(i,j);
					/*fprintf (stderr,"\n ESPECIAL, codon_NRmat = %lf, Qij_C[%d][%d] = %lf", codon_NRmat(i,j), i, j, Qij_C[i][j]);*/
					}
				}
			}
		}	
	
	
	/* MAIN diagonal (row sum = 0) */
	mr = 0;
	for (i = 0; i < NUMCOD; i++)
		for (j = 0; j < NUMCOD; j++)
			k[i] = k[i] + Qij_C[i][j];
	for (i = 0; i < NUMCOD; i++)
		for (j = 0; j < NUMCOD; j++)
			if (i == j)
				{
				Qij_C[i][j] = 0.0;
				Qij_C[i][j] = 0.0-k[i];
				
				mr -= Qij_C[i][j]*codonTable_frequencies(j); /* mr, will be Eigen input */
				}
	/* scala factor is good (Felsenstein book, pag 205) */

	/* Qij_CC[4096], matrix Qij_C in just a vector */
	for (i = 0; i < NUMCOD; i++)
		{
		for (j = 0; j < NUMCOD; j++)
			{
			Qij_CC[m] = Qij_C[i][j];
			/*fprintf (stderr, "\nQij_CC[%d] = %3.2f ", m, Qij_CC[m]);*/
			m++;
			}
		}
	
	/* Active to see the Qij matrix */
	/*fprintf (stderr, "\nIn the matrix Qij_C, End:");
	for (i = 0; i < NUMCOD; i++)
		for (j = 0; j < NUMCOD; j++)
			fprintf (stderr, "\nQij_C[%d][%d] = %lf ", i, j, Qij_C[i][j]);*/
			
	/* Doing eigen */
	EigenREV_Codon(Root_C, Cijk_C); /* Root_C y Cijk_C bilds in eigen. Input Qij_CC and mr. EYE!, Qij_CC change by Eigen */
	
	
	/*for (i = 0; i< NUMCOD*NUMCOD*NUMCOD*NUMCOD;i++)
		fprintf (stderr, " \n Cijk_C [%d] = %3.2f ", i, Cijk_C[i]);*/
	/*for (i = 0; i< NUMCOD;i++)
		fprintf (stderr, " \n Root_C [%d] = %3.2f ", i, Root_C[i]);*/
	}








/************************** makeCodonFromNuc *************************/
/* It Builds the codon from the 3 nucleotides (numbers). It returns the codon number from the three nucleotides number */
int makeCodonFromNuc (int x, int y, int z)
	{

	if ( x < 0 || x > 3 || y < 0 || y > 3 || z < 0 || z > 3)
		{
		fprintf (fpmpi, "\nWarning in makeCodonFromNuc function.");		
		exit(-1);
		}
		
	if (x == 0)
		{
		if (y == 0)
			{
			if (z == 0)
				return (0); /* AAA */
			else if (z == 1)
				return (1); /* AAC */
			else if (z == 2)
				return (2); /* AAG */
			else
				return (3); /* AAT */
			}
		else if (y == 1)
			{
			if (z == 0)
				return (4); /* ACA */
			else if (z == 1)
				return (5); /* ACC */
			else if (z == 2)
				return (6); /* ACG */
			else
				return (7); /* ACT */
			}
		else if (y == 2)
			{
			if (z == 0)
				return (8); /* AGA */
			else if (z == 1)
				return (9); /* AGC */
			else if (z == 2)
				return (10); /* AGG */
			else
				return (11); /* AGT */
			}
		else
			{
			if (z == 0)
				return (12); /* ATA */
			else if (z == 1)
				return (13); /* ATC */
			else if (z == 2)
				return (14); /* ATG */
			else
				return (15); /* ATT */
			}
		}
	else if (x == 1)
		{
		if (y == 0)
			{
			if (z == 0)
				return (16); /* CAA */
			else if (z == 1)
				return (17); /* CAC */
			else if (z == 2)
				return (18); /* CAG */
			else
				return (19); /* CAT */
			}
		else if (y == 1)
			{
			if (z == 0)
				return (20); /* CCA */
			else if (z == 1)
				return (21); /* CCC */
			else if (z == 2)
				return (22); /* CCG */
			else
				return (23); /* CCT */
			}
		else if (y == 2)
			{
			if (z == 0)
				return (24); /* CGA */
			else if (z == 1)
				return (25); /* CGC */
			else if (z == 2)
				return (26); /* CGG */
			else
				return (27); /* CGT */
			}
		else
			{
			if (z == 0)
				return (28); /* CTA */
			else if (z == 1)
				return (29); /* CTC */
			else if (z == 2)
				return (30); /* CTG */
			else
				return (31); /* CTT */
			}
		}
	else if (x == 2)
		{
		if (y == 0)
			{
			if (z == 0)
				return (32); /* GAA */
			else if (z == 1)
				return (33); /* GAC */
			else if (z == 2)
				return (34); /* GAG */
			else
				return (35); /* GAT */
			}
		else if (y == 1)
			{
			if (z == 0)
				return (36); /* GCA */
			else if (z == 1)
				return (37); /* GCC */
			else if (z == 2)
				return (38); /* GCG */
			else
				return (39); /* GCT */
			}
		else if (y == 2)
			{
			if (z == 0)
				return (40); /* GGA */
			else if (z == 1)
				return (41); /* GGC */
			else if (z == 2)
				return (42); /* GGG */
			else
				return (43); /* GGT */
			}
		else
			{
			if (z == 0)
				return (44); /* GTA */
			else if (z == 1)
				return (45); /* GTC */
			else if (z == 2)
				return (46); /* GTG */
			else
				return (47); /* GTT */
			}
		}
	else								/* changing to 61 codons erasing the stop codons */
		{
		if (y == 0)
			{
			if (z == 0) /* TAA */	/* STOP CODON!! */
				{				
				fprintf (fpmpi, "\n warning in makeCodonFromNuc by an TAA stop codon");
				exit (-1);
				} 
			else if (z == 1)
				return (48); /* TAC */
			else if (z == 2) /* TAG */	/* STOP CODON!! */
				{
				fprintf (fpmpi, "\n warning in makeCodonFromNuc by an TAG stop codon");				
				exit (-1);
				}
			else
				return (49); /* TAT */
			}
		else if (y == 1)
			{
			if (z == 0)
				return (50); /* TCA */
			else if (z == 1)
				return (51); /* TCC */
			else if (z == 2)
				return (52); /* TCG */
			else
				return (53); /* TCT */
			}
		else if (y == 2)
			{
			if (z == 0) /* TGA */	/* STOP CODON!! */
				{				
				fprintf (fpmpi, "\n warning in makeCodonFromNuc by an TGA stop codon");
				exit (-1);
				}
			else if (z == 1)
				return (54); /* TGC */
			else if (z == 2)
				return (55); /* TGG */
			else
				return (56); /* TGT */
			}
		else
			{
			if (z == 0)
				return (57); /* TTA */
			else if (z == 1)
				return (58); /* TTC */
			else if (z == 2)
				return (59); /* TTG */
			else
				return (60); /* TTT */
			}
		}
	return (61);
	}




/************************* codonTable_frequencies_MRCA ***************************/
/******** adapted from codontable function written by Rasmus Nielsen **********/
/* This function generates the frequencies for each codon. */

double codonTable_frequencies_MRCA(int cod)
	{
	/*A = 0, C = 1, G = 2, T = 3*/
	int i, j;
	double b;
	double A[3], C[3], G[3], T[3];
	
	j = 0;
	
	
	for (i = 0; i < 3; i++)
		{
		A[i] = p_i_codon[i+j];
		C[i] = p_i_codon[i+1+j];
		G[i] = p_i_codon[i+2+j];
		T[i] = p_i_codon[i+3+j];
		j = j+3;		
		/*fprintf (fpmpi, "\n A[%d] = %lf, C[%d] = %lf, G[%d] = %lf, T[%d] = %lf", i, A[i], i, C[i], i, G[i], i, T[i]);*/
		}
	
	if (cod > 63)
		{		
		fprintf (fpmpi, "\n Warning in codonTable_frequencies_MRCA function. cod >= 64.");
		return -1;
		}
	if (cod < 16){									/*** A ***/
		if (cod < 4){								
			if (cod == 0) b = A[0]*A[1]*A[2];		/* AAA */
			else if (cod == 1) b = A[0]*A[1]*C[2];	/* AAC */
			else if (cod == 2) b = A[0]*A[1]*G[2];	/* AAG */
			else b = A[0]*A[1]*T[2];				/* AAT */
			}
		else if (cod > 11){								
			cod = cod - 12;
			if (cod == 0) b = A[0]*T[1]*A[2];		/* ATA */
			else if (cod == 1) b = A[0]*T[1]*C[2];	/* ATC */
			else if (cod == 2) b = A[0]*T[1]*G[2];	/* ATG */
			else b = A[0]*T[1]*T[2];				/* ATT */
			}
		else if (cod > 7){								
			cod = cod - 8;
			if (cod == 0) b = A[0]*G[1]*A[2];		/* AGA */
			else if (cod == 1) b = A[0]*G[1]*C[2];	/* AGC */
			else if (cod == 2) b = A[0]*G[1]*G[2];	/* AGG */
			else b = A[0]*G[1]*T[2];				/* AGT */
			}
		else {
			if (cod == 4) b = A[0]*C[1]*A[2];		/* ACA */
			else if (cod == 5) b = A[0]*C[1]*C[2];	/* ACC */
			else if (cod == 6) b = A[0]*C[1]*G[2];	/* ACG */
			else b = A[0]*C[1]*T[2];				/* ACT */
			}
		}
	else if (cod > 47){								/*** T ***/
		cod = cod - 48;	
		if (cod < 4){										
			if (cod == 0) 
				{
				b = T[0]*A[1]*A[2];		/* TAA */		/* STOP CODON!! */
				b = 0.0;
				}
			else if (cod == 1) b = T[0]*A[1]*C[2];	/* TAC */
			else if (cod == 2)
				{
				b = T[0]*A[1]*G[2];	/* TAG */		/* STOP CODON!! */
				b = 0.0;
				}
			else b = T[0]*A[1]*T[2];				/* TAT */
			}
		else if (cod > 11){								
			cod = cod - 12;
			if (cod == 0) b = T[0]*T[1]*A[2];		/* TTA */
			else if (cod == 1) b = T[0]*T[1]*C[2];	/* TTC */
			else if (cod == 2) b = T[0]*T[1]*G[2];	/* TTG */
			else b = T[0]*T[1]*T[2];				/* TTT */
			}
		else if (cod > 7){							
			cod = cod - 8;
			if (cod == 0) 
				{
				b = T[0]*G[1]*A[2];		/* TGA */		/* STOP CODON!! */
				b = 0.0;
				}
			else if (cod == 1) b = T[0]*G[1]*C[2];	/* TGC */
			else if (cod == 2) b = T[0]*G[1]*G[2];	/* TGG */
			else b = T[0]*G[1]*T[2];				/* TGT */
			}
		else {
			if (cod == 52) b = T[0]*C[1]*A[2];		/* TCA */
			else if (cod == 53) b = T[0]*C[1]*C[2];	/* TCC */
			else if (cod == 54) b = T[0]*C[1]*G[2];	/* TCG */
			else b = T[0]*C[1]*T[2];				/* TCT */
			}									   	
		}
	else if (cod < 32){								/*** C ***/
		cod = cod - 16;
		if (cod < 4){								
			if (cod == 0) b = C[0]*A[1]*A[2];		/* CAA */
			else if (cod == 1) b = C[0]*A[1]*C[2];	/* CAC */
			else if (cod == 2) b = C[0]*A[1]*G[2];	/* CAG */
			else b = C[0]*A[1]*T[2];				/* CAT */
			}
		else if (cod > 11){ 
			if (cod == 28) b = C[0]*T[1]*A[2];		/* CTA */
			else if (cod == 29) b = C[0]*T[1]*C[2];	/* CTC */
			else if (cod == 30) b = C[0]*T[1]*G[2];	/* CTG */
			else b = C[0]*T[1]*T[2];				/* CTT */
			}
		else if (cod > 7){ 
			if (cod == 24) b = C[0]*G[1]*A[2];		/* CGA */
			else if (cod == 25) b = C[0]*G[1]*C[2];	/* CGC */
			else if (cod == 26) b = C[0]*G[1]*G[2];	/* CGG */
			else b = C[0]*G[1]*T[2];				/* CGT */
			}							
		else {
			if (cod == 20) b = C[0]*C[1]*A[2];		/* CCA */
			else if (cod == 21) b = C[0]*C[1]*C[2];	/* CCC */
			else if (cod == 22) b = C[0]*C[1]*G[2];	/* CCG */
			else b = C[0]*C[1]*T[2];				/* CCT */
			}											
		}
	else if (cod > 31 && cod < 48){					/*** G ***/
		cod = cod - 32;
		if (cod < 4){									
			if (cod == 0) b = G[0]*A[1]*A[2];		/* GAA */
			else if (cod == 1) b = G[0]*A[1]*C[2];	/* GAC */
			else if (cod == 2) b = G[0]*A[1]*G[2];	/* GAG */
			else b = G[0]*A[1]*T[2];				/* GAT */
			}
		else if (cod > 11){
			if (cod == 44) b = G[0]*T[1]*A[2];		/* GTA */
			else if (cod == 45) b = G[0]*T[1]*C[2];	/* GTC */
			else if (cod == 46) b = G[0]*T[1]*G[2];	/* GTG */
			else b = G[0]*T[1]*T[2];				/* GTT */
			}							
		else if (cod > 7){
			if (cod == 40) b = G[0]*G[1]*A[2];		/* GGA */
			else if (cod == 41) b = G[0]*G[1]*C[2];	/* GGC */
			else if (cod == 42) b = G[0]*G[1]*G[2];	/* GGG */
			else b = G[0]*G[1]*T[2];				/* GGT */
			}											
		else{
			if (cod == 36) b = G[0]*C[1]*A[2];		/* GCA */
			else if (cod == 37) b = G[0]*C[1]*C[2];	/* GCC */
			else if (cod == 38) b = G[0]*C[1]*G[2];	/* GCG */
			else b = G[0]*C[1]*T[2];				/* GCT */
			}				 										
		}
	else {
			fprintf (fpmpi, "error in codonTable_frequencies_MRCA function\n"); 
			exit(-1);
		}
	
	/*fprintf (fpmpi, "\n In codonTable_frequencies_MRCA, cod = %d, b = %lf", cod, b);*/
	
	return b;
	}


/************************* codonTable_frequencies ***************************/
/******** adapted from codontable function written by Rasmus Nielsen **********/
/* This function generates the frequencies for each codon. */

double codonTable_frequencies(int cod)
	{
	/*A = 0, C = 1, G = 2, T = 3*/
	int i, j;
	double b;
	double A[3], C[3], G[3], T[3];
	
	j = 0;
	
	/*
	fprintf (fpmpi, "\n In codonTable_frequencies, cod = %d,", cod);
	*/
	
	
	for (i = 0; i < 3; i++)
		{
		A[i] = p_i_codon[i+j];
		C[i] = p_i_codon[i+1+j];
		G[i] = p_i_codon[i+2+j];
		T[i] = p_i_codon[i+3+j];
		j = j+3;
		
		/*
		fprintf (fpmpi, "\n A[%d] = %lf, C[%d] = %lf, G[%d] = %lf, T[%d] = %lf", i, A[i], i, C[i], i, G[i], i, T[i]);		
		*/
		}

	
	if (cod > 60)
		return -1;
	if (cod>=48)
		cod++;
	if (cod>=50)
		cod++;
	if (cod>=56)
		cod++;
	if (cod < 16){									/*** A ***/
		if (cod < 4){								
			if (cod == 0) b = A[0]*A[1]*A[2];		/* AAA */
			else if (cod == 1) b = A[0]*A[1]*C[2];	/* AAC */
			else if (cod == 2) b = A[0]*A[1]*G[2];	/* AAG */
			else b = A[0]*A[1]*T[2];				/* AAT */
			}
		else if (cod > 11){								
			cod = cod - 12;
			if (cod == 0) b = A[0]*T[1]*A[2];		/* ATA */
			else if (cod == 1) b = A[0]*T[1]*C[2];	/* ATC */
			else if (cod == 2) b = A[0]*T[1]*G[2];	/* ATG */
			else b = A[0]*T[1]*T[2];				/* ATT */
			}
		else if (cod > 7){								
			cod = cod - 8;
			if (cod == 0) b = A[0]*G[1]*A[2];		/* AGA */
			else if (cod == 1) b = A[0]*G[1]*C[2];	/* AGC */
			else if (cod == 2) b = A[0]*G[1]*G[2];	/* AGG */
			else b = A[0]*G[1]*T[2];				/* AGT */
			}
		else {
			if (cod == 4) b = A[0]*C[1]*A[2];		/* ACA */
			else if (cod == 5) b = A[0]*C[1]*C[2];	/* ACC */
			else if (cod == 6) b = A[0]*C[1]*G[2];	/* ACG */
			else b = A[0]*C[1]*T[2];				/* ACT */
			}
		}
	else if (cod > 47){								/*** T ***/
		cod = cod - 48;	
		if (cod < 4){										
			if (cod == 0) 
				{
				b = T[0]*A[1]*A[2];		/* TAA */		/* STOP CODON!! */
				b = 0.0;
				
				fprintf (fpmpi, "\n warning in codonTable_frequencies by an TAA stop codon");
				exit (-1);
				}
			else if (cod == 1) b = T[0]*A[1]*C[2];	/* TAC */
			else if (cod == 2)
				{
				b = T[0]*A[1]*G[2];	/* TAG */		/* STOP CODON!! */
				b = 0.0;
				
				fprintf (fpmpi, "\n warning in codonTable_frequencies by an TAG stop codon");
				exit (-1);
				}
			else b = T[0]*A[1]*T[2];				/* TAT */
			}
		else if (cod > 11){								
			cod = cod - 12;
			if (cod == 0) b = T[0]*T[1]*A[2];		/* TTA */
			else if (cod == 1) b = T[0]*T[1]*C[2];	/* TTC */
			else if (cod == 2) b = T[0]*T[1]*G[2];	/* TTG */
			else b = T[0]*T[1]*T[2];				/* TTT */
			}
		else if (cod > 7){							
			cod = cod - 8;
			if (cod == 0) 
				{
				b = T[0]*G[1]*A[2];		/* TGA */		/* STOP CODON!! */
				b = 0.0;
				
				fprintf (fpmpi, "\n warning in codonTable_frequencies by an TGA stop codon");				
				exit (-1);
				}
			else if (cod == 1) b = T[0]*G[1]*C[2];	/* TGC */
			else if (cod == 2) b = T[0]*G[1]*G[2];	/* TGG */
			else b = T[0]*G[1]*T[2];				/* TGT */
			}
		else {
			if (cod == 52) b = T[0]*C[1]*A[2];		/* TCA */
			else if (cod == 53) b = T[0]*C[1]*C[2];	/* TCC */
			else if (cod == 54) b = T[0]*C[1]*G[2];	/* TCG */
			else b = T[0]*C[1]*T[2];				/* TCT */
			}									   	
		}
	else if (cod < 32){								/*** C ***/
		cod = cod - 16;
		if (cod < 4){								
			if (cod == 0) b = C[0]*A[1]*A[2];		/* CAA */
			else if (cod == 1) b = C[0]*A[1]*C[2];	/* CAC */
			else if (cod == 2) b = C[0]*A[1]*G[2];	/* CAG */
			else b = C[0]*A[1]*T[2];				/* CAT */
			}
		else if (cod > 11){ 
			if (cod == 28) b = C[0]*T[1]*A[2];		/* CTA */
			else if (cod == 29) b = C[0]*T[1]*C[2];	/* CTC */
			else if (cod == 30) b = C[0]*T[1]*G[2];	/* CTG */
			else b = C[0]*T[1]*T[2];				/* CTT */
			}
		else if (cod > 7){ 
			if (cod == 24) b = C[0]*G[1]*A[2];		/* CGA */
			else if (cod == 25) b = C[0]*G[1]*C[2];	/* CGC */
			else if (cod == 26) b = C[0]*G[1]*G[2];	/* CGG */
			else b = C[0]*G[1]*T[2];				/* CGT */
			}							
		else {
			if (cod == 20) b = C[0]*C[1]*A[2];		/* CCA */
			else if (cod == 21) b = C[0]*C[1]*C[2];	/* CCC */
			else if (cod == 22) b = C[0]*C[1]*G[2];	/* CCG */
			else b = C[0]*C[1]*T[2];				/* CCT */
			}											
		}
	else if (cod > 31 && cod < 48){					/*** G ***/
		cod = cod - 32;
		if (cod < 4){									
			if (cod == 0) b = G[0]*A[1]*A[2];		/* GAA */
			else if (cod == 1) b = G[0]*A[1]*C[2];	/* GAC */
			else if (cod == 2) b = G[0]*A[1]*G[2];	/* GAG */
			else b = G[0]*A[1]*T[2];				/* GAT */
			}
		else if (cod > 11){
			if (cod == 44) b = G[0]*T[1]*A[2];		/* GTA */
			else if (cod == 45) b = G[0]*T[1]*C[2];	/* GTC */
			else if (cod == 46) b = G[0]*T[1]*G[2];	/* GTG */
			else b = G[0]*T[1]*T[2];				/* GTT */
			}							
		else if (cod > 7){
			if (cod == 40) b = G[0]*G[1]*A[2];		/* GGA */
			else if (cod == 41) b = G[0]*G[1]*C[2];	/* GGC */
			else if (cod == 42) b = G[0]*G[1]*G[2];	/* GGG */
			else b = G[0]*G[1]*T[2];				/* GGT */
			}											
		else{
			if (cod == 36) b = G[0]*C[1]*A[2];		/* GCA */
			else if (cod == 37) b = G[0]*C[1]*C[2];	/* GCC */
			else if (cod == 38) b = G[0]*C[1]*G[2];	/* GCG */
			else b = G[0]*C[1]*T[2];				/* GCT */
			}				 										
		}
	else {
			
		fprintf (fpmpi, "error in codonTable_frequencies function\n"); 			
		exit(-1);
		}
		
	return b;
	}






/****************************** codonTable_DnDs *******************************/
/******** adapted from codontable function written by Rasmus Nielsen **********/
/* This function generates the aminoacid for each codon. It can use to compare two codons and say if they are synonimous or nonsynonymous
For stop codons we will use b = 25 like the number of aa */
int codonTable_DnDs(int cod)
	{
	/*A = 0, C = 1, G = 2, T = 3*/
	int b;

	if (cod > 60)
		return -1;
	if (cod>=48)
		cod++;
	if (cod>=50)
		cod++;
	if (cod>=56)
		cod++;
	if (cod < 16){									/*A*/
		if (cod < 4){										/*A*/
			if (cod==3||cod==1) b=13;
			else b=14;
			}
		else if (cod > 11){								/*T*/
			cod = cod - 12;
			if (cod==2) b=4;
			else b=3;
			}
		else if (cod > 7){								/*G*/
			cod = cod - 8;
			if (cod==3||cod==1) b=6;
			else b=19;
			}
		else b=8;											/*C*/
		}
	else if (cod > 47){							/*T*/
		cod = cod - 48;	
		if (cod < 4){										/*A*/
			if (cod==3||cod==1) b=10;
			else 
				{
				b=25;	/* two stop codons */
				fprintf (stderr, "\n warning in codonTable_DnDs by a stop codon");
				exit (-1);
				}
			}
		else if (cod > 11){								/*T*/
			cod = cod - 12;
			if (cod==3||cod==1) b=1;
			else b=2;
			}
		else if (cod > 7){								/*G*/
			cod = cod - 8;
			if (cod==3||cod==1) b=17;
			else if (cod==0)
				{
				b=25; /* one stop codon */
				fprintf (stderr, "\n warning in codonTable_DnDs by a stop codon");
				exit (-1);
				}
			else b=18;
			}
		else b=6;									   	/*C*/
		}
	else if (cod < 32){							/*C*/
		cod = cod - 16;
		if (cod < 4){										/*A*/
			if (cod==3||cod==1) b=11;
			else b=12;
			}
		else if (cod > 11) b=2;							/*T*/
		else if (cod > 7) b=19;							/*G*/
		else b=7;											/*C*/
		}
	else if (cod > 31 && cod < 48){			/*G*/
		cod = cod - 32;
		if (cod < 4){										/*A*/
			if (cod==3||cod==1) b=15;
			else b=16;
			}
		else if (cod > 11) b=5;							/*T*/
		else if (cod > 7) b=0;							/*G*/
		else b=9;											/*C*/
		}
	else {fprintf (stderr, "error in codon table DnDs\n"); exit(-1);}
	return b;
	}



/****************************** numdif_codon *******************************/
/* adapted from numdif function written by Rasmus Nielsen */
/* Compare 2 codons by nucletotide diferences (return 0-3 diferences) */
int numdif_codon(int codon1, int codon2)
	{
	int i, returner = 0;
	double cod1[3], cod2[3];
	
	
	if (codon1>=48)
		codon1++;
	if (codon1>=50)
		codon1++;
	if (codon1>=56)
		codon1++;	
	if (codon2>=48)
		codon2++;
	if (codon2>=50)
		codon2++;
	if (codon2>=56)
		codon2++;
	cod1[0] = floor((double)codon1/16.0);
	cod1[1] = floor(((double)codon1 - 16.0*cod1[0])/4.0);
	cod1[2] = (double)codon1 - 4.0*cod1[1] - 16.0*cod1[0];
	cod2[0] = floor((double)codon2/16.0);
	cod2[1] = floor(((double)codon2 - 16.0*cod2[0])/4.0);
	cod2[2] = (double)codon2 - 4.0*cod2[1] - 16.0*cod2[0];
	for (i=0; i<3; i++)
		{
		if (cod1[i] != cod2[i])
			returner++;
		}
	return returner; /* 0, 1, 2 or 3 */
	}


/****************************** number_to_codon *******************************/
/* adapted from number_to_codon function written by Rasmus Nielsen */
/* from one codon number to the three nucleotides letters */
 void number_to_codon(int ind, char out[])
	{
	int i, codon[3];
	
	if (ind < 0) /* important check, not eraser */
		{
		fprintf (fpmpi, "\n Warning in number_to_codon: ind = %d", ind);
		exit (-1);
		}


	if (ind > 60)
		{
		fprintf (fpmpi, "\n stop codon in number_to_codon %d\n", ind); 
		exit(-1);
		}
	if (ind>=48)
		ind++;
	if (ind>=50)
		ind++;
	if (ind>=56)
		ind++;
	if (ind == 48 || ind == 50 || ind == 56) /* Cheking */
		{
		fprintf (fpmpi, "\n Warning in number_to_codon2");
		exit (-1);
		}
	codon[0] = floor((double)ind/16.0);
	codon[1] = floor((ind - 16.0*(double)codon[0])/4.0);
	codon[2] = ind - 4*codon[1] - 16*codon[0];
	for (i=0; i<3; i++)
		{
		if (codon[i] == 0)
			out[i] = 'A';
		else if (codon[i] == 1)
			out[i] = 'C';
		else if (codon[i] == 2)
			out[i] = 'G';
		else if (codon[i] == 3)
			out[i] = 'T';
		else {			
			fprintf (fpmpi, "error in codon (%i)\n",ind); 			
			exit(-1);
			}
		}
	}


/****************************** number_to_codon2 *******************************/
/* adapted from number_to_codon function written by Rasmus Nielsen */
/* from one codon number to the three nucleotides letters */
 void number_to_codon2(int ind, int out[])
	{
	/*int i, codon[3];*/
	
	if (ind < 0) /* important check, not eraser */
		{
		fprintf (fpmpi, "\n Warning in number_to_codon2: ind = %d", ind);
		exit (-1);
		}


	if (ind > 60)
		{
		fprintf (fpmpi, "\n stop codon in number_to_codon2 %d\n", ind); 
		exit(-1);
		}
	if (ind>=48)
		ind++;
	if (ind>=50)
		ind++;
	if (ind>=56)
		ind++;
	if (ind == 48 || ind == 50 || ind == 56) /* Cheking */
		{
		fprintf (fpmpi, "\n Warning in number_to_codon2");
		exit (-1);
		}
	out[0] = floor((double)ind/16.0);
	out[1] = floor((ind - 16.0*(double)out[0])/4.0);
	out[2] = ind - 4*out[1] - 16*out[0];
	}









/****************************** number_to_codon_MRCA *******************************/
/* adapted from number_to_codon function written by Rasmus Nielsen */
/* from one codon number to the three nucleotides numbers for MRCA*/
 void number_to_codon_MRCA(int ind, int codon[])
	{
	codon[0] = floor((double)ind/16.0);
	codon[1] = floor((ind - 16.0*(double)codon[0])/4.0);
	codon[2] = ind - 4*codon[1] - 16*codon[0];
	}






/****************************** codon_tr_tv *******************************/
/* adapted from tr_tv function written by Rasmus Nielsen */
/* It determines, between two codons, if 1 mutation is a transversion or a transition 
(return 0 or 1. return 2 if the value for this codon change is null (several substitutions or main diagonal)).
 It's for the codon model version HKY */
 int codon_tr_tv(int indi, int indj)
	{
	char codon1[3], codon2[3];
	int i, cod1, cod2;
	
	if (indi > 60 || indj > 60)
		{
		fprintf (stderr, "\n stop codon1 \n");
		exit(-1);
		}
		
	number_to_codon(indi, codon1);
	number_to_codon(indj, codon2);
	
	for (i=0; i<3; i++)
		{
		if (codon1[i] != codon2[i])
			{
			if (codon1[i] == 'A' || codon1[i] == 'G')
				cod1 = 0;
			else cod1 = 1;
			if (codon2[i] == 'A' || codon2[i] == 'G')
				cod2 = 0;
			else cod2 = 1;
			return abs(cod1 - cod2);
			}
		}
	
	fprintf (stderr, "\nERROR IN codon_tr_tv function!!");
	exit(-1);
	return (0);
	}



/****************************** codon_Rmat *******************************/
/*determines, between two codons, the corresponding rate for the codon model version GTR. 
It returns the rate of the change between these codons*/
 double codon_Rmat(int indi, int indj)
	{
	char codon1[3], codon2[3];
	int i/*, k*/; /* Rmat[6] from parameters input file*/
	double rate;
	
	rate = 1.0;
	/*k = 0;*/
	
	number_to_codon(indi, codon1);
	number_to_codon(indj, codon2);
	if (indi > 60 || indj > 60)
		{
		fprintf (stderr, "\n stop codon2 \n");
		exit(-1);
		}
	
	for (i=0; i<3; i++)
		{
		if (codon1[i] != codon2[i]) /* Important for codon with just 1 difference (more than 1 difference will be null..) */
			{
			/*k++;*/
			if (codon1[i] == 'A' && codon2[i] == 'C')
				rate = Rmat[0];
			else if (codon1[i] == 'A' && codon2[i] == 'G')
				rate = Rmat[1];
			else if (codon1[i] == 'A' && codon2[i] == 'T')
				rate = Rmat[2];
			else if (codon1[i] == 'C' && codon2[i] == 'G')
				rate = Rmat[3];
			else if (codon1[i] == 'C' && codon2[i] == 'T')
				rate = Rmat[4];
			else if (codon1[i] == 'G' && codon2[i] == 'T')
				rate = Rmat[5];
		
			else if (codon2[i] == 'A' && codon1[i] == 'C')
				rate = Rmat[0];
			else if (codon2[i] == 'A' && codon1[i] == 'G')
				rate = Rmat[1];
			else if (codon2[i] == 'A' && codon1[i] == 'T')
				rate = Rmat[2];
			else if (codon2[i] == 'C' && codon1[i] == 'G')
				rate = Rmat[3];
			else if (codon2[i] == 'C' && codon1[i] == 'T')
				rate = Rmat[4];
			else if (codon2[i] == 'G' && codon1[i] == 'T')
				rate = Rmat[5];
			else
				{
				fprintf (stderr, "\nERROR IN codon_Rmat function!!");
				exit(-1);
				}
			}
		}
	/*if (k == 0 || k > 1)*/ /* main diagonal principal and more than one substitution = 0 */
	/*	return 0;*/
	
	return rate;
	}

/****************************** codon_NRmat *******************************/
/*determines, between two codons, the corresponding rate for the codon model version NOT REVERSIBLE GTR. 
It returns the rate of the change between these codons*/
 double codon_NRmat(int indi, int indj)
	{
	char codon1[3], codon2[3];
	int i/*, k*/; /* Rmat[12] from parameters input file*/
	double rate;
	
	rate = 1.0;
	/*k = 0;*/
	
	number_to_codon(indi, codon1);
	number_to_codon(indj, codon2);
	if (indi > 60 || indj > 60)
		{
		fprintf (stderr, "\n stop codon3 \n");
		exit(-1);
		}
		
	for (i=0; i<3; i++) /*	NRmat: AC CA AG GA AT TA CG GC CT TC GT=1 TG */
		{
		if (codon1[i] != codon2[i]) /* Important for codon with just 1 difference (more than 1 difference will be null..) */
			{
			/*k++;*/
			if (codon1[i] == 'A' && codon2[i] == 'C')
				rate = NRmat[0];
			else if (codon1[i] == 'C' && codon2[i] == 'A')
				rate = NRmat[1];
			else if (codon1[i] == 'A' && codon2[i] == 'G')
				rate = NRmat[2];
			else if (codon1[i] == 'G' && codon2[i] == 'A')
				rate = NRmat[3];
			else if (codon1[i] == 'A' && codon2[i] == 'T')
				rate = NRmat[4];
			else if (codon1[i] == 'T' && codon2[i] == 'A')
				rate = NRmat[5];
			else if (codon1[i] == 'C' && codon2[i] == 'G')
				rate = NRmat[6];
			else if (codon1[i] == 'G' && codon2[i] == 'C')
				rate = NRmat[7];
			else if (codon1[i] == 'C' && codon2[i] == 'T')
				rate = NRmat[8];
			else if (codon1[i] == 'T' && codon2[i] == 'C')
				rate = NRmat[9];
			else if (codon1[i] == 'G' && codon2[i] == 'T')
				rate = NRmat[10];
			else if (codon1[i] == 'T' && codon2[i] == 'G')
				rate = NRmat[11];
		
			/*else if (codon2[i] == 'A' && codon1[i] == 'C')
				rate = NRmat[0];
			else if (codon2[i] == 'C' && codon1[i] == 'A')
				rate = NRmat[1];
			else if (codon2[i] == 'A' && codon1[i] == 'G')
				rate = NRmat[2];
			else if (codon2[i] == 'G' && codon1[i] == 'A')
				rate = NRmat[3];
			else if (codon2[i] == 'A' && codon1[i] == 'T')
				rate = NRmat[4];
			else if (codon2[i] == 'T' && codon1[i] == 'A')
				rate = NRmat[5];
			else if (codon2[i] == 'C' && codon1[i] == 'G')
				rate = NRmat[6];
			else if (codon2[i] == 'G' && codon1[i] == 'C')
				rate = NRmat[7];
			else if (codon2[i] == 'C' && codon1[i] == 'T')
				rate = NRmat[8];
			else if (codon2[i] == 'T' && codon1[i] == 'C')
				rate = NRmat[9];
			else if (codon2[i] == 'G' && codon1[i] == 'T')
				rate = NRmat[10];
			else if (codon2[i] == 'T' && codon1[i] == 'G')
				rate = NRmat[11];*/
			else
				{
				fprintf (stderr, "\nERROR IN codon_NRmat function!!");
				exit(-1);
				}
			}
		}
	
	return rate;
	}


/*************** CodonModel **********************/
/* It builds the Pij matrix from Cijk (variable global) and time */
void CodonModel (double Pij[NUMCOD][NUMCOD], double branchLength, double varRate) /*, long int *seed)*/
	{	
	int 	i, j, k;
	double	t, expt[NUMCOD]/*, GammaVarRateOmega*/; 
	/* double 			Rmat_C[66], Qij_CC[4096], Cijk_C[16777216], mr, tstv, Root_C[NUMCOD], Qij_C[NUMCOD][NUMCOD];*/ /* global variables */
	double a;
	
	a = 0;
	k = 0;
	
	/*fprintf (fpmpi, "\nIN CODON MODEL FUNCTION\n");*/
	/*fprintf (fpmpi, "\n*** NORMAL omega = %lf ***\n", omega);*/
	/* for each site and each branch */
	/*if (doOmegaRateHetCont == YES) 
		{
		GammaVarRateOmega = RndGamma (OmegaRateHet, seed) / OmegaRateHet; 
		omega = OmegaInit*GammaVarRateOmega;
		buildCodonMatrix_Qij_Cijk ();
		}*/
	
	/*fprintf (fpmpi, "\n varRate = %lf", varRate);*/
	if (varRate > 0)
		varRate = varRate / (1.0 - pinv);
	t = branchLength * varRate;
	
	if (t<1e-6) 
		{ 
		/*fprintf (fpmpi, "\n TIEMPO BAJO.   branchLength = %lf, varRate = %lf", branchLength, varRate);*/
		for (i=0; i<NUMCOD; i++) 
			for (j=0; j<NUMCOD; j++) 
				{
				if (i==j)
					Pij[i][j] = 1.0;
				else 	
					Pij[i][j] = 0.0;
				}
		}
	else
		{
		for (k=1; k<NUMCOD; k++) 
			expt[k]=exp(t*Root_C[k]);
		
		for (i=0; i<NUMCOD; i++) 
			for (j=0; j<NUMCOD; j++) 
				{
				Pij[i][j] = Cijk_C[i*NUMCOD*NUMCOD+j*NUMCOD+0]; 
				for (k=1; k<NUMCOD; k++)
					Pij[i][j]+=Cijk_C[i*NUMCOD*NUMCOD+j*NUMCOD+k]*expt[k];			/* THIS SLOWS DOWN THE PROGRAM!!! */ /* 15,5 seg, without 0,7 seg*/
				}
		}
	
	
	/* Active this to see the Pij matrix */
	
	/*fprintf (fpmpi, "\n");
	for (k=0; k<NUMCOD; k++)
		{
		for (i=0; i<NUMCOD; i++) 
			fprintf (fpmpi, "\nPij[%d][%d] = %lf ", k, i, Pij[k][i]);
		}*/
	
	
	
	/*for (i=0; i<NUMCOD; i++)
		{
	//	fprintf (stderr, "\n");
		a = 0;
		for (j=0; j<NUMCOD; j++)
			{
	//		fprintf (stderr, "PCodonM[%d][%d] = %3.2f ",i,j,Pij[i][j]);
			a = a + Pij[i][j];
			}
	
		if (a>1.1 || a <0.99)
			{
			fprintf (stderr, "\n a = %lf \n\n", a);
			exit (-5);
			}
		}
	*/
	}



/*************** CodonModel_Cat **********************/
/* It builds the Pij matrix from Cijk (variable global) and time */
void CodonModel_Cat (double Pij[NUMCOD][NUMCOD], double branchLength, double varRate) /* long int *seed, int numOmegaCat */ 
	{	
	int 	i, j, k, category;
	double	t, expt[NUMCOD]/*, ran*/; 
	/*double	*cumProbCat, *hetProb;*/
	/* double 			Rmat_C[66], Qij_CC[4096], Cijk_C[16777216], mr, tstv, Root_C[NUMCOD], Qij_C[NUMCOD][NUMCOD];*/ /* global variables */
	
	k = 0;
	category = 0;
	/*cumProbCat =  (double*) calloc ((numOmegaCat+1), sizeof (double));  
	if (cumProbCat == NULL)
		{
		fprintf (fpmpi, "Could not allocate cumProbCat (%lu bytes)", (numOmegaCat+1) * (long) sizeof (double));
		exit(1);
		}
	hetProb =  (double*) calloc ((numOmegaCat+1), sizeof (double));  
	if (hetProb == NULL)
		{
		fprintf (fpmpi, "Could not allocate hetProb (%lu bytes)", (numOmegaCat+1) * (long) sizeof (double));
		exit(1);
		}*/

	
		
	if (doOmegaProb == YES)
		{

		/*for (j=1; j<=numOmegaCat; j++)
			cumProbCat[j] = 0;
		cumProbCat[0] = 0;
		for (j=1; j<=numOmegaCat; j++)
			cumProbCat[j] = cumProbCat[j-1] + omegaProb[j];
		ran = RandomUniform(seed);
		category = bbinInOmegaCat (ran, cumProbCat, numOmegaCat);
		j = 0;*/
		category = ProbCategory;
		/*fprintf (stderr, "\n category = %d \n", category);*/
		}

	if (doOmegaRateHetDisc == YES)
		{
		/*for (j=1; j<=numOmegaCat; j++)
			{
			cumProbCat[j] = 0;
			hetProb[j] = 1.0/numOmegaCat;
			}
		cumProbCat[0] = 0;
		for (j=1; j<=numOmegaCat; j++)
			cumProbCat[j] = cumProbCat[j-1] + hetProb[j];
		ran = RandomUniform(seed);
		category = bbinInOmegaCat (ran, cumProbCat, numOmegaCat);
		j = 0;*/
		category = GammCategory;
		}
		
	if (doM1 == YES)
		{
		category = M1_FinalSite_omega;
		/*fprintf (stderr, "\n category = %d \n", category);*/
		}
	
	
	/*fprintf (fpmpi, "\n*** CAT category = %d ***\n", category);*/
	
	/*free (cumProbCat);
	free (hetProb);*/
	
	if (varRate > 0)
		varRate = varRate / (1.0 - pinv);
	t = branchLength * varRate;
	
	if (t<1e-6) 
		{ 
		for (i=0; i<NUMCOD; i++) 
			for (j=0; j<NUMCOD; j++) 
				{
				if (i==j)
					Pij[i][j] = 1.0;
				else 	
					Pij[i][j] = 0.0;
				}
		}
	else
		{
		for (k=1; k<NUMCOD; k++)
			expt[k]=exp(t*QijOmegas[category-1].Root_C_cat[k]);
			
		for (i=0; i<NUMCOD; i++) 
			for (j=0; j<NUMCOD; j++) 
				{
				Pij[i][j] = QijOmegas[category-1].Cijk_C_cat[i*NUMCOD*NUMCOD+j*NUMCOD+0]; 
				for (k=1; k<NUMCOD; k++)
					Pij[i][j]+=QijOmegas[category-1].Cijk_C_cat[i*NUMCOD*NUMCOD+j*NUMCOD+k]*expt[k];			/* THIS SLOWS DOWN THE PROGRAM!!! */ /* example, with 15,5 seg, without 0,7 seg*/
				}
		}
		
		
	/* Active this to see the Pij matrix */
	/*		
	fprintf (fpmpi, "\n");
	for (k=0; k<NUMCOD; k++)
		{
		for (i=0; i<NUMCOD; i++) 
			fprintf (fpmpi, "\nPij[%d][%d] = %lf ", k, i, Pij[k][i]);
		}
	*/

	}





/*********************** Nucleotidic Substitution Models ***************************************/
/***********************************************************************************************/

/********************************** EvolveSequenceOnTree_NEW ***********************************/
/* This function evolves sequences on the coalescent trees by nucleotide model */

void EvolveSequenceOnTree_NEW (long int *seed, double m, double kappa, 
							double alpha, double p_i[4], int numNuc, int indNumRE, int *arrayIndBreakpointsOrd, char *MRCAsequence, int numSites)
	{
	int			i, w, j, sitePosition;
	double		varRate, ran, cumPi[4];
	int			*arrayIndBreakpointsOrd_C;
	int 		GMRCA_label;
	/*TreeNode	*q;*/

	w = sitePosition = j = 0;
	GMRCA_label = -1;
	/*fprintf (fpmpi, "En EvolveSequenceOnTree_NEW function \n");*/	

	arrayIndBreakpointsOrd_C = (int *)calloc((indNumRE+1),(long) sizeof(int));
	if (!arrayIndBreakpointsOrd_C)
		{
		fprintf (fpmpi, "Could not allocate arrayIndBreakpointsOrd_C (%lu bytes)\n", (indNumRE+1)  * (long) sizeof(int));
		exit (1);
		}
	

	/* making codon breakpoints */
	for (i = 0; i < indNumRE+1; i++)
		arrayIndBreakpointsOrd_C[i] = fabs(arrayIndBreakpointsOrd[i]/3);
	

	/** Label of GMRCA **/
	w = 0; /* NETRECODON Ahora solo hay un nodo raiz, "GMRCA" base de toda la red */
	GMRCA_label =  TellMeGMRCALabel (treeRootInit[w]); 
	if (GMRCA_label < 0) /* MRCA from File */
		{
		fprintf (fpmpi, "Error in GMRCA_label in the EvolveSequenceOnTree_NEW function \n");
		exit (1);
		}
	/*fprintf (fpmpi, "GMRCA LAbel %d\n", GMRCA_label);*/	

	if (doMRCAFile == YES) /* MRCA from File */
		{
		for (i = 1; i <= numNuc; i++)		/* taking MRCA sequence */
			{	
			/* form OLD . La secuencia de entrada es GMRCA */
			if (i == arrayIndBreakpointsOrd[w]) /* recombinations, breakpoints. CAMBIA LA RAIZ DEL ARBOL */
				w++;
			w = 0; /* NETRECODON Ahora solo hay un nodo raiz, "GMRCA" base de toda la red */
			
			matrix[pos(GMRCA_label,i,numNuc)] = WhichNucNumber(MRCAsequence[i-1]); /* from MRCA input file */
			}
		}
	else /* MRCA from nucleotide frequencies */
		{		
		for (i = 1; i <= numNuc; i++)
			{
			cumPi[0] = p_i[0];
			for (j = 1; j < 4; j++)
				cumPi[j] = cumPi[j-1] + p_i[j];
			
			ran = RandomUniform(seed);

			if (ran >= 0.0 && ran <= cumPi[0])
				matrix[pos(GMRCA_label,i,numSites)] = 0;
			else if (ran > cumPi[0] && ran <= cumPi[1])
				matrix[pos(GMRCA_label,i,numSites)] = 1;
			else if (ran > cumPi[1] && ran <= cumPi[2])
				matrix[pos(GMRCA_label,i,numSites)] = 2;
			else
				matrix[pos(GMRCA_label,i,numSites)] = 3;

			}
		}
	



	/*fprintf (fpmpi, "\n Entra en la recursion nucleotidica, EvolveSequenceOnTree_NEW, numRE = %d", numRE);*/
	w = 0;
	
	for (i = 1; i <= numSites; i++)	/************* Other sequences */ /* Hace esto para cada nucleotido FOR EACH SITE (NUCLEOTIDE)**********/
		{
		/*fprintf (fpmpi, "\n\n +++++++ NUCLEOTIDE = %d ++++++++ \n\n", i);*/

		if (RandomUniform(seed) < pinv)		
			varRate = 0.0;
		else
			{
			if (doRateHet == YES)
				varRate = RndGamma (alpha, seed) / alpha; 
			else
				varRate = 1; 
			}

		/* form OLD . La secuencia de entrada es GMRCA */
		if (i == (arrayIndBreakpointsOrd_C[w]+1) && indNumRE != 0)
			w++;
		w = 0; /* NETRECODON Ahora solo hay un nodo raiz, "GMRCA" base de toda la red */



		SimulateDataForSite_Nucleotide_RECURSIVE_NET (treeRootInit[w], i, numSites, m, varRate, kappa, seed);
		}
	/*fprintf (fpmpi, "\n\n SALE DE LA RECURSION \n\n");*/
	


	free (arrayIndBreakpointsOrd_C);

	}








/********************************** TellMeGMRCALabel ***********************************/
/* It returns the GMRCA label */

int		TellMeGMRCALabel (TreeNode *p)
	{
	int			j;
	j = -1;
	
	if (p != NULL)
		{
		if (p->anc1 == NULL) /* root */		
			{
			j = p->label;
			}
		}
	return (j);
	}





/*void EvolveSequenceOnTree (long int *seed, double m, double kappa, 
							double alpha, double p_i[4], int numSites, int *arrayIndBreakpointsOrd, char *MRCAsequence)
{
	int			i, w;
	double		varRate;
	
	w = 0;
	
	for (i = 1; i <= numSites; i++)
		{
	  	if (RandomUniform(seed) < pinv)		
	  		varRate = 0.0;
	  	else
	  		{
		  	if (doRateHet == YES)
		  		varRate = RndGamma (alpha, seed) / alpha; 
			else
				varRate = 1; 
			}
		
		if (i == arrayIndBreakpointsOrd[w])
			w++;
*/
		/*fprintf (fpmpi, "\n w = %d, treeRootNodex[w] label = %d, index = %d  \n", w, treeRootNodex[w]->label, treeRootNodex[w]->index);*/
/*
		SimulateDataForSite (treeRootNodex[w], i, numSites, m, kappa, p_i, varRate, seed, MRCAsequence);
		}
}*/



/********************************** SimulateDataForSite ***********************************/
/* Simulates the nucleotide substitution process for a given site */
/*
void SimulateDataForSite (TreeNodex *p, int siteNum, int numSites, double m, 
							double kappa, double p_i[4], double varRate, long int *seed, char *MRCAsequence)
{
	int			i;
	double		ran, cumPi[4], cumProb[4], Pij[4][4];
	
		
	if (p != NULL)
		{
		if (p->anc1 == NULL) *//* root */	
/*			{
			if (doMRCAFile == YES)
				{
				matrix[pos(p->label,siteNum,numSites)] = WhichNucNumber(MRCAsequence[siteNum-1]); *//* from MRCA input file */
/*				}
			else
				{
				cumPi[0] = p_i[0];
				for (i = 1; i < 4; i++)
					cumPi[i] = cumPi[i-1] + p_i[i];
				ran = RandomUniform(seed);
*/
				/*if (p->label == 1)
					{
					fprintf(stderr,"\n ***** p->label = %d, p->index = %d,  p->time = %lf, sitioNum = %d, deme = %d, node_original = %d \n", p->label, p->index,  p->time, siteNum, p->indexOldMigPop, p->NetIndex);
					}*/
				/*fprintf (fpmpi, "\n p->label = %d \n", p->label);*/
/*
				if (ran >= 0.0 && ran <= cumPi[0])
					matrix[pos(p->label,siteNum,numSites)] = 0;
				else if (ran > cumPi[0] && ran <= cumPi[1])
					matrix[pos(p->label,siteNum,numSites)] = 1;
				else if (ran > cumPi[1] && ran <= cumPi[2])
					matrix[pos(p->label,siteNum,numSites)] = 2;
				else
					matrix[pos(p->label,siteNum,numSites)] = 3;
				}
			}
		else			
			{
			if (p->isOutgroup == YES)
				SubstitutionMatrix (Pij, p->length * m, kappa, varRate, p_i);
			else 				
				SubstitutionMatrix (Pij, (p->anc1->time - p->time) * m, kappa, varRate, p_i);
			
			cumProb[0] = Pij[matrix[pos(p->anc1->label,siteNum,numSites)]][0];
			for (i=1; i<4; i++)
				cumProb[i] = cumProb[i-1] + Pij[matrix[pos(p->anc1->label,siteNum,numSites)]][i];
			ran = RandomUniform(seed);
*/
			/*if (p->label == 1)
				{
				fprintf(stderr,"\n ***** p->label = %d, p->index = %d,  p->time = %lf, (p->anc1->time - p->time) = %lf, sitioNum = %d, deme = %d, node_original = %d \n", p->label, p->index,  p->time, (p->anc1->time - p->time), siteNum, p->indexOldMigPop, p->NetIndex);
				}*/
/*			if (ran >= 0.0 && ran <= cumProb[0])
				matrix[pos(p->label,siteNum,numSites)] = 0; 
			else if (ran > cumProb[0] && ran <= cumProb[1])
				matrix[pos(p->label,siteNum,numSites)] = 1; 
			else if (ran > cumProb[1] && ran <= cumProb[2])
				matrix[pos(p->label,siteNum,numSites)] = 2; 
			else
				matrix[pos(p->label,siteNum,numSites)] = 3; 
			
			if (matrix[pos(p->label,siteNum,numSites)] != matrix[pos(p->anc1->label,siteNum,numSites)])
				numMU++;	
			}
*/		
		/* It crosses the tree */
/*		SimulateDataForSite (p->left, siteNum, numSites, m, kappa, p_i, varRate, seed, MRCAsequence);
		SimulateDataForSite (p->right, siteNum, numSites, m, kappa, p_i, varRate, seed, MRCAsequence);	
		if (thereisOutgroup == YES)
			SimulateDataForSite (p->outgroup, siteNum, numSites, m, kappa, p_i, varRate, seed, MRCAsequence);	
		}
}
*/


/********************* Substitution matrix **********************/
/* Sets the apropriate model of nucleotide substitution  */
static void SubstitutionMatrix (double ch_prob[4][4], double branchLength, double kappa, double varRate, double p_i[4])
{
	if (doHKY == YES)
		HKY (ch_prob, branchLength, kappa, varRate, p_i);
	if (doGTnR == YES)
		GTnR (ch_prob, branchLength, varRate, p_i);
	if (doGTR == YES)
		GTR (ch_prob, branchLength, varRate, p_i);
}





/*********************************** HKY **************************************/
/*	HKY performs Hasegawa-Kishino-Yano 85 correction */ 
 
void HKY (double Pij[4][4], double branchLength, double kappa, double varRate, double p_i[4])
{
	int			i, j;
	double		A, t, PIj, beta; 

	beta = 0.5 / ((p_i[0] + p_i[2])*(p_i[1] + p_i[3]) + kappa*((p_i[0]*p_i[2]) + (p_i[1]*p_i[3])));

	if (varRate > 0)
		varRate = varRate / (1.0 - pinv);

	t = branchLength * varRate;
	/*fprintf(stderr,"\n t = %lf \n", t);*/	

	if (t == 0.0)		/* no mutations */
		{
		for (i=0; i<4; i++)
			{
			for (j=0; j<4; j++)
				{
				if (i == j)
					Pij[i][j] = 1.0;
				else 
					Pij[i][j] = 0.0;
				}
			}
		}
	else				/*there are mutations */
		for (i=0; i<4; i++)
		  	{
		  	for (j=0; j<4; j++)
				{
				if (j == 0 || j == 2)	/* purine */
					PIj = p_i[0] + p_i[2];
				else
					PIj = p_i[1] + p_i[3]; /* pyrimidine */
					
				A = 1 + PIj*(kappa-1);
				
				if (i==j) /* diagonal principal */
					Pij[i][j] = p_i[j] + p_i[j]*(1/PIj - 1)*exp(-beta*t) + ((PIj-p_i[j])/PIj)*exp(-beta*t*A);
				else if ((i==0 && j==2) || (i==1 && j==3) || (i==2 && j==0) || (i==3 && j==1)) /* transition */
					Pij[i][j] = p_i[j] + p_i[j]*(1/PIj - 1)*exp(-beta*t) - (p_i[j]/PIj)*exp(-beta*t*A);
				else /* transversion */
					Pij[i][j] = p_i[j]*(1-exp(-beta*t));
				}
			}
}






/*************** GTR **********************/
void GTR (double Pij[4][4], double branchLength, double varRate, double p_i[4])
{	
	int 	i, j, k;
	double	t, expt[4];
	/* double Rmat[6], Qij[16], Cijk[256], Root[4], mr, tstv; Global Variables*/
	
	if (varRate > 0)
		varRate = varRate / (1.0 - pinv);

	t = branchLength * varRate;

	k=0;
	for (i=0; i<3; i++) 
		for (j=i+1; j<4; j++)
      		if (i*4+j != 11)
				Qij[i*4+j]=Qij[j*4+i]=Rmat[k++];
				
	Qij[3*4+2]=Qij[2*4+3]=1.0;
	
	for (i=0; i<4; i++) 
		for (j=0; j<4; j++)
			Qij[i*4+j] *= p_i[j];
				
	mr=0;		
	for (i=0; i<4; i++) {
		Qij[i*4+i]=0; 
		Qij[i*4+i]=-(Qij[i*4]+Qij[i*4+1]+Qij[i*4+2]+Qij[i*4+3]); 
		
		mr-=p_i[i]*Qij[i*4+i];
		}

	EigenREV(Root, Cijk);
	
	/* calculate mean ts/tv ratio */
	mr=2*(p_i[3]*Qij[3*4+1]+p_i[0]*Qij[0*4+2]);
	tstv=mr/(1-mr); 
	
	
	/* P(t)ij = SUM Cijk * exp{Root*t}*/
	if (t<1e-6) 
		{ 
		for (i=0; i<4; i++) 
			for (j=0; j<4; j++) 
				{
				if (i==j)
					Pij[i][j] = 1.0;
				else 	
					Pij[i][j] = 0.0;
				}
		}
	else
		{
		for (k=1; k<4; k++) 
			expt[k]=exp(t*Root[k]);
		for (i=0; i<4; i++) 
			for (j=0; j<4; j++) 
				{
				Pij[i][j]=Cijk[i*4*4+j*4+0];
				for (k=1; k<4; k++)
					Pij[i][j]+=Cijk[i*4*4+j*4+k]*expt[k];
				}
		}
		
	
	/*
	for (i=0; i<4; i++)
		{
		fprintf(fpmpi,"\n");
		for (j=0; j<4; j++)
			{
			fprintf(fpmpi,"%5.4f ",Pij[i][j]);
			}
		}
		fprintf(fpmpi,"\n\n");
	*/
}





/*************** GTnR **********************/
/* GTR non reversible */
void GTnR (double Pij[4][4], double branchLength, double varRate, double p_i[4])
{	
	int 	i, j, k;
	double	t, expt[4];
		
	if (varRate > 0)
		varRate = varRate / (1.0 - pinv);

	t = branchLength * varRate;
	
	/*	
		A	C	G	T
	A	0	1	2	3
	C	4	5	6	7
	G	8	9	10	11
	T	12	13	14	15
	*/

	/* fills no symetrical matrix */
	k=0;
	for (i=0; i<3; i++) 
		for (j=i+1; j<4; j++)
			{
			Qij[i*4+j]=NRmat[k++];
			Qij[j*4+i]=NRmat[k++];
			}

	
/*	AC CA AG GA AT TA CG GC CT TC GT=1 TG */
	
	/* all rates relative to GT */
	Qij[11] = 1.0;
	
	for (i=0; i<4; i++) 
		for (j=0; j<4; j++) 
			Qij[i*4+j] *= p_i[j];
			
	mr=0;		
	for (i=0; i<4; i++) 
		{
		Qij[i*4+i]=0; 
		Qij[i*4+i]=-(Qij[i*4]+Qij[i*4+1]+Qij[i*4+2]+Qij[i*4+3]); 
		mr-=p_i[i]*Qij[i*4+i];
		}
	
	EigenREV(Root, Cijk);
	
/* calculate mean ts/tv ratio */ /*double check*/
/*	mr=2*(p_i[3]*Qij[3*4+1]+p_i[0]*Qij[0*4+2]);*/
	mr = p_i[3]*Qij[3*4+1] + p_i[0]*Qij[0*4+2] + p_i[1]*Qij[1*4+3] + p_i[2]*Qij[2*4+0] ;
	tstv=mr/(1-mr);
	
/* P(t)ij = SUM Cijk * exp{Root*t}
*/
	if (t<1e-6) /* too small branch */
		{ 
		for (i=0; i<4; i++) 
			for (j=0; j<4; j++) 
				{
				if (i==j)
					Pij[i][j] = 1.0;
				else 	
					Pij[i][j] = 0.0;
				}
		}
	else
		{
		for (k=1; k<4; k++) 
			expt[k]=exp(t*Root[k]);
		for (i=0; i<4; i++) 
			for (j=0; j<4; j++) 
				{
				Pij[i][j]=Cijk[i*4*4+j*4+0];
				for (k=1; k<4; k++)
					Pij[i][j]+=Cijk[i*4*4+j*4+k]*expt[k];
				}
		}
	
	/*
	for (i=0; i<4; i++)
		{
		fprintf(fpmpi,"\n");
		for (j=0; j<4; j++)
			{
			fprintf(fpmpi,"%5.4f ",Pij[i][j]);
			}
		}
		fprintf(fpmpi,"\n\n");
		
	*/
}







/********************* RandomExponential ********************/
/* Generates a random number from a Poisson distibution with
  mean lambda. 
*/

double RandomExponential (double lambda, long int *seed)
{

	double 	exponentialNumber, U;

	do
		U = RandomUniform (seed);
	while (U == 0);

	exponentialNumber = -log (U) / lambda;

	return exponentialNumber;
}






/***************************** RandomUniform **********************************/
/* It returns a random uniform variate in range 0..1. It is described in
  	Park, S. K. and K. W. Miller. 1988. Random number generators: good
   ones are hard to find. Communications of the ACM, 31(10):1192-1201.
*/

double RandomUniform (long int *seed)
{

	long int	lo, hi, test;
	
	hi = (*seed) / 127773;
	lo = (*seed) % 127773;
	test = 16807 * lo - 2836 * hi;
	if (test > 0)
		*seed = test;
	else
		*seed = test + 2147483647;
	return (double)(*seed) / (double)2147483647;

}






/* Gamma functions written by Ziheng Yang */
double RndGamma (double s, long int *seed)
{
	double r=0.0;
	
	if (s <= 0.0)     
		puts ("jgl gamma..");
	else if (s < 1.0) 
		r = RndGamma1 (s, seed);
	else if (s > 1.0) 
		r = RndGamma2 (s, seed);
	else          
		r =- log(RandomUniform(seed));
	return (r);
}

double RndGamma1 (double s, long int *seed)
{
	double			r, x=0.0, smal=1e-37, w;
	static double  a, p, uf, ss=10.0, d;
	
	if (s!=ss) 
		{
		a = 1.0-s;
		p = a/(a+s*exp(-a));
		uf = p*pow(smal/a,s);
		d = a*log(a);
		ss = s;
		}
	for (;;) 
		{
		r = RandomUniform(seed);
		if (r > p)       
			x = a-log((1.0-r)/(1.0-p)), w=a*log(x)-d;
		else if (r>uf) 
			x = a*pow(r/p,1/s), w=x;
		else           
			return (0.0);
		r = RandomUniform(seed);
		if (1.0-r <= w && r > 0.0)
		if (r*(w+1.0) >= 1.0 || -log(r) <= w) 
			continue;
		break;
		}
	return (x);
}

double RndGamma2 (double s, long int *seed)
{
	double			r ,d, f, g, x;
	static double	b, h, ss=0;
	
	if (s!=ss) 
		{
		b = s-1.0;
		h = sqrt(3.0*s-0.75);
		ss = s;
		}
	for (;;) 
		{
		r = RandomUniform(seed);
		g = r-r*r;
		f = (r-0.5)*h/sqrt(g);
		x = b+f;
		if (x <= 0.0) 
			continue;
		r = RandomUniform(seed);
		d = 64*r*r*g*g*g;
		if (d*x < x-2.0*f*f || log(d) < 2*(b*log(x/b)-f)) 
			break;
		}
	return (x);
}




/*
	Program: SimulateBeta
	Purpose: To simulate random beta variables 
 	Author: David Posada 
 	Started: October 08
*/
/*** Thanks David! ***/

/**************************** RandomBeta*************************/
/*	Generates a beta number 
	Beta (p,q) = gamma (1,p) / [ gamma (1,p) + gamma (1,q)]
	Here  gamma (scale, shape)


	Proof-> In R, type: mean(rbeta(1000,0.3,0.5)) to get the mean of 1000 beta variables
	var(rbeta(1000,0.3,0.5)) to get the variance.
	repeat with 1000000 to get more accuracy.


*/

double	RandomBeta (double p, double q, long int *seed)
{
	double betaNumber = 0;
	double gammaNumber1, gammaNumber2;
	
	/* Note RandomGamm function already assumes a scale=1 */
	gammaNumber1 = RandomGamma (p, seed); /* I think we should not divide by p here */
	gammaNumber2 = RandomGamma (q, seed); /* I think we should not divide by q here */
	betaNumber = gammaNumber1 / (gammaNumber1 + gammaNumber2);
	
	return (betaNumber);
}




/**************************** RandomGamma *************************/
/*	Generates a gamma number using routines in Ziheng's
	Yang tools.h in PAML

	Random standard gamma (Mean=Var=s,  with shape par=s, scale par=1)
	r^(s-1)*exp(-r)
	
	J. Dagpunar (1988) Principles of random variate generation,
	Clarendon Press, Oxford
	
	Calling rndgamma1() if s<1 or rndgamma2() if s>1 or exponential if s=1
*/	

double	RandomGamma (double shape, long int *seed)
{
	double gammaNumber = 0;
	
	if (shape <= 0)
		fprintf (stderr, "ERROR: problems with gamma variable generation, shape < 0");
	else if (shape < 1)
		gammaNumber = RandomGamma1 (shape, seed);
	else if (shape > 1)
		gammaNumber = RandomGamma2 (shape, seed);
	else
		gammaNumber = -log (RandomUniform(seed));
	
	return (gammaNumber);


}

/*************** RandomGamma1 ***************/
double RandomGamma1 (double s, long int *seed)
{
/* Random standard gamma for s<1
   switching method
*/
	double		r, x = 0.0, small2 = 1e-37, w;
	static double	a, p, uf, ss2 = 10.0, d;
	
	if (s != ss2) 
		{
		a  = 1.0 - s;
		p  = a/(a+s*exp(-a));
		uf = p*pow(small2/a,s);
		d  = a*log(a);
		ss2 = s;
		}
	for (;;) 
		{
		r = RandomUniform(seed);
		if (r > p)        
			x = a-log((1.0-r)/(1.0-p)), w=a*log(x)-d;
		else if (r>uf)  
			x = a*pow(r/p,1/s), w=x;
		else            
			return (0.0);
		r = RandomUniform(seed);
		if (1.0-r <= w && r > 0.0)
		if (r*(w+1.0) >= 1.0 || -log(r) <= w)  
			continue;
		break;
		}
	return (x);
}


/*************** RandomGamma2 ***************/
double RandomGamma2 (double s, long int *seed)
{
/* Random standard gamma for s>1
   Best's (1978) t distribution method
*/
	double		r ,d, f, g, x;
	static double	b, h, ss=0;
	
	if (s!=ss) 
		{
		b  = s-1.0;
		h  = sqrt(3.0*s-0.75);
		ss = s;
		}
	for (;;) 
		{
		r = RandomUniform(seed);
		g = r-r*r;
		f = (r-0.5)*h/sqrt(g);
		x = b+f;
		if (x <= 0.0) 
			continue;
		r = RandomUniform(seed);
		d = 64*r*r*g*g*g;
		if (d*x < x-2.0*f*f || log(d) < 2*(b*log(x/b)-f))  
			break;
		}
	return (x);
}







/*void SetMatrix(double Pij[4][4], double len)
{	
	int i,j,k;
	double expt[4];*/
	
/* P(t)ij = SUM Cijk * exp{Root*t}
*/
/*	if (len<1e-6) 
		{ 
		for (i=0; i<4; i++) 
			for (j=0; j<4; j++) 
				{
				if (i==j)
					Pij[i][j] = 1.0;
				else 	
					Pij[i][j] = 0.0;
				}
		return; 
		}
	
	for (k=1; k<4; k++) 
		expt[k]=exp(len*Root[k]);
	for (i=0; i<4; i++) 
		for (j=0; j<4; j++) 
			{
			Pij[i][j]=Cijk[i*4*4+j*4+0];
			for (k=1; k<4; k++)
				Pij[i][j]+=Cijk[i*4*4+j*4+k]*expt[k];
			}*/
	 
	
/* the rows are cumulative to help with picking one using
  a random number */ /*matrix = Pij */
/*	matrix[1]+=matrix[0];
	matrix[2]+=matrix[1];
	matrix[3]+=matrix[2];*/ /* This should always be 1.0... 

	matrix[5]+=matrix[4];
	matrix[6]+=matrix[5];
	matrix[7]+=matrix[6];*/ /* ...but it is easier to spot bugs... 
	
	matrix[9]+=matrix[8];
	matrix[10]+=matrix[9];
	matrix[11]+=matrix[10];*/ /* ...though less efficient... 
	
	matrix[13]+=matrix[12];
	matrix[14]+=matrix[13];
	matrix[15]+=matrix[14];*/ /* ...but probably not much. 
*/
/*}*/


/*void SetVector(double *vector, short base, double len)
{	
	int i,j,k;
	double expt[4];
	double *P;

	P=vector;
	if (len<1e-6) 
		{ 
		for (i=0; i<4; i++) 
			{
			if (i==base)
				*P=1.0;
			else 	
				*P=0.0;
			P++;
			}
		return; 
		}
	for (k=1; k<4; k++) 
		expt[k]=exp(len*Root[k]);

	for (j=0; j<4; j++)
		{
		(*P)=Cijk[base*4*4+j*4+0];
		for (k=1; k<4; k++)
			(*P)+=Cijk[base*4*4+j*4+k]*expt[k];
		P++;
		}
		vector[1]+=vector[0];
		vector[2]+=vector[1];
		vector[3]+=vector[2];
}*/




/* Everything below is shamelessly taken from Yang's Paml package */

int abyx (double a, double x[], int n);
int xtoy (double x[], double y[], int n);
int matinv( double x[], int n, int m, double space[]);
int eigen(int job, double A[], int n, double rr[], double ri[],
         double vr[], double vi[], double w[]);
void balance(double mat[], int n, int *low, int *hi, double scale[]);
void unbalance(int n, double vr[], double vi[], int low, int hi,
              double scale[]);
int realeig(int job, double mat[], int n,int low, int hi, double valr[],
           double vali[], double vr[], double vi[]);
void elemhess(int job, double mat[], int n, int low, int hi, 
           double vr[], double vi[], int work[]);

typedef struct { double re, im; } complex;
#define csize(a) (fabs(a.re)+fabs(a.im))

complex compl (double re,double im);
complex cconj (complex a);
complex cplus (complex a, complex b);
complex cminus (complex a, complex b);
complex cby (complex a, complex b);
complex cdiv (complex a,complex b);
complex ccexp (complex a);
complex cfactor (complex x, double a);
int cxtoy (complex x[], complex y[], int n);
int cmatby (complex a[], complex b[], complex c[], int n,int m,int k);
int cmatout (FILE * fout, complex x[], int n, int m);
int cmatinv( complex x[], int n, int m, double space[]);


/* Eigen function for codon models */
int EigenREV_Codon (double Root_C[], double Cijk_C[])
{

	int i,j,k;
	double U[NUMCOD*NUMCOD], V[NUMCOD*NUMCOD], T1[NUMCOD*NUMCOD], T2[NUMCOD*NUMCOD];
	
	
	nDIGITS = 40; /* no. of digits to the base BASE in the fraction */
	abyx (1/mr, Qij_CC, NUMCOD*NUMCOD);
	
	if ((k=eigen (1, Qij_CC, NUMCOD, Root_C, T1, U, V, T2))!=0) {
		fprintf(stderr, "\ncomplex roots in EigenREV_Codon (k = %d)", k);
		exit(0);
	}
	xtoy (U, V, NUMCOD*NUMCOD);
	matinv (V, NUMCOD, NUMCOD, T1);
	
	for (i=0; i<NUMCOD; i++) 
   		for (j=0; j<NUMCOD; j++) 
   			for (k=0; k<NUMCOD; k++)
   				Cijk_C[i*NUMCOD*NUMCOD+j*NUMCOD+k] = U[i*NUMCOD+k]*V[k*NUMCOD+j];
	
	/*fprintf(stderr, "\nCijk_C[10*NUMCOD*NUMCOD+10*NUMCOD+10] = %d\n", Cijk_C[10*NUMCOD*NUMCOD+10*NUMCOD+10]);
	fprintf(stderr, "\nCijk_C[5*NUMCOD*NUMCOD+5*NUMCOD+5] = %d\n", Cijk_C[5*NUMCOD*NUMCOD+5*NUMCOD+5]);*/			
			
	return (0);
}




/* Eigen function for nucleotide models */
int EigenREV (double Root[], double Cijk[])
{

	int i,j,k;
	double U[16], V[16], T1[16], T2[16];
	
	nDIGITS = 53; /* no. of digits to the base BASE in the fraction */
	abyx (1/mr, Qij, 16);

	if ((k=eigen (1, Qij, 4, Root, T1, U, V, T2))!=0) 
		{
		fprintf(fpmpi, "\ncomplex roots in EigenREV");
		exit(0);
		}
	xtoy (U, V, 16);
	matinv (V, 4, 4, T1);
	for (i=0; i<4; i++) 
   		for (j=0; j<4; j++) 
   			for (k=0; k<4; k++)
   				Cijk[i*4*4+j*4+k] = U[i*4+k]*V[k*4+j];
				
	return (0);
}





int abyx (double a, double x[], int n)
{ int i; for (i=0; i<n; x[i]*=a,i++) ; return(0); }
int xtoy (double x[], double y[], int n)
{ int i; for (i=0; i<n; y[i]=x[i],i++) ; return(0); }
int matinv( double x[], int n, int m, double space[])
{
/* x[n*m] ... m>=n
*/
  register int i,j,k;
  int *irow=(int*) space;
  double ee=1.0e-20, t,t1,xmax;
  double det=1.0;

  for (i=0; i<n; i++) {
     xmax = 0.;
     for (j=i; j<n; j++) {
	 if (xmax < fabs(x[j*m+i])) {
	   xmax = fabs( x[j*m+i] );
	   irow[i] = j;
	 }
     }
     det *= xmax;
     if (xmax < ee)  {
	 
	 #ifdef MPI
		fprintf (stderr, "\nDet becomes zero at %3d!\t\n", i+1);
	 #else
		fprintf (stderr, "\nDet becomes zero at %3d!\t\n", i+1);
	 #endif
	 return(-1);
     }
     if (irow[i] != i) {
	 for (j=0; j<m; j++) {
	   t = x[i*m+j];
	   x[i*m+j] = x[irow[i] * m + j];
	   x[ irow[i] * m + j] = t;
	 }
     }
     t = 1./x[i*m+i];
     for (j=0; j<n; j++) {
	 if (j == i) continue;
	 t1 = t*x[j*m+i];
	 for (k=0; k<n; k++) x[j*m+k] -= t1*x[i*m+k];
	 x[j*m+i] = -t1;
     }
     for (j=0; j<m; j++)  x[i*m+j] *= t;
     x[i*m+i] = t;
  }                           /* i */
  for (i=n-1; i>=0; i--) {
     if (irow[i] == i) continue;
     for (j=0; j<n; j++) {
	 t = x[j*m+i];
	 x[j*m+i] = x[ j*m + irow[i] ];
	 x[ j*m + irow[i] ] = t;
     }
  }
  return (0);
}





/***********************************************************
* This eigen() works for eigenvalue/vector analysis
*        for real general square matrix A
*        A will be destroyed
*        rr,ri are vectors containing eigenvalues
*        vr,vi are matrices containing (right) eigenvectors
*
*             A*[vr+vi*i] = [vr+vi*i] * diag{rr+ri*i}
*
* Algorithm: Handbook for Automatic Computation, vol 2
*            by Wilkinson and Reinsch, 1971
*            most of source codes were taken from a public domain
*            solftware called MATCALC.
* Credits:  to the authors of MATCALC
*
* return    -1 not converged
*             0 no complex eigenvalues/vectors
*             1 complex eigenvalues/vectors
* Tianlin Wang at University of Illinois
* Thu May 6 15:22:31 CDT 1993
***************************************************************/

#define FOR(i,n) for(i=0; i<n; i++)
#define FPN(file) fputc('\n', file)
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

#define pos(i,j,n)     ((i)*(n)+(j))

#define BASE       2   /* base of floating point arithmetic */
#define MAXITER   30   /* max. no. of iterations to converge */	
	
int eigen(int job, double A[], int n, double rr[], double ri[], 
         double vr[], double vi[], double work[])
{   
/* double work[n*2]: working space
*/
   int low,hi,i,j,k, it, istate=0;
   double tiny=sqrt(pow((double)BASE,(double)(1-nDIGITS))), t; 
	
   balance(A,n,&low,&hi,work);
   elemhess(job,A,n,low,hi,vr,vi, (int*)(work+n));

	if (-1 == realeig(job,A,n,low,hi,rr,ri,vr,vi)) return (-1);
	if (job) unbalance(n,vr,vi,low,hi,work);

/* sort, added by Z. Yang */
  for (i=0; i<n; i++) {
      for (j=i+1,it=i,t=rr[i]; j<n; j++)
          if (t<rr[j]) { t=rr[j]; it=j; }
      rr[it]=rr[i];  rr[i]=t;
      t=ri[it];      ri[it]=ri[i]; ri[i]=t;
      for (k=0; k<n; k++) {
         t=vr[k*n+it]; vr[k*n+it]=vr[k*n+i]; vr[k*n+i]=t;
         t=vi[k*n+it]; vi[k*n+it]=vi[k*n+i]; vi[k*n+i]=t;
      }
      if (fabs(ri[i])>tiny) istate=1;
  }

   return (istate) ;
}

/* complex funcctions
*/

complex compl (double re,double im)
{
   complex r;

   r.re = re;
   r.im = im;
   return(r);
}

complex cconj (complex a)
{
   a.im = -a.im;
   return(a);
}

#define csize(a) (fabs(a.re)+fabs(a.im))

complex cplus (complex a, complex b)
{
  complex c;
  c.re = a.re+b.re; 
  c.im = a.im+b.im;  
  return (c);
}

complex cminus (complex a, complex b)
{
  complex c;
  c.re = a.re-b.re; 
  c.im = a.im-b.im;  
  return (c);
}

complex cby (complex a, complex b)
{
  complex c;
  c.re = a.re*b.re-a.im*b.im ;
  c.im = a.re*b.im+a.im*b.re ;
  return (c);
}

complex cdiv (complex a,complex b)
{
   double ratio, den;
   complex c;

   if (fabs(b.re) <= fabs(b.im)) {
       ratio = b.re / b.im;
       den = b.im * (1 + ratio * ratio);
       c.re = (a.re * ratio + a.im) / den;
       c.im = (a.im * ratio - a.re) / den;
   }
   else {
       ratio = b.im / b.re;
       den = b.re * (1 + ratio * ratio);
       c.re = (a.re + a.im * ratio) / den;
       c.im = (a.im - a.re * ratio) / den;
   }
   return(c);
}

complex ccexp (complex a)
{
  complex c;
  c.re = exp(a.re);
  if (fabs(a.im)==0) c.im = 0; 
  else { c.im = c.re*sin(a.im); c.re*=cos(a.im); }
  return (c);
}

complex cfactor (complex x, double a)
{
  complex c;
  c.re = a*x.re; 
  c.im = a*x.im;
  return (c);
}

int cxtoy (complex x[], complex y[], int n)
{
  int i;
  FOR (i,n) y[i]=x[i];
  return (0);
}

int cmatby (complex a[], complex b[], complex c[], int n,int m,int k)
/* a[n*m], b[m*k], c[n*k] ...... c = a*b
*/
{
  int i,j,i1;
  complex t;

  FOR (i,n) FOR(j,k) {
      for (i1=0,t=compl(0,0); i1<m; i1++) 
          t = cplus (t, cby(a[i*m+i1],b[i1*k+j]));
      c[i*k+j] = t;
  }
  return (0);
}

int cmatout (FILE * fout, complex x[], int n, int m)
{
  int i,j;
  for (i=0,FPN(fout); i<n; i++,FPN(fout)) 
       FOR(j,m) fprintf(fout, "%7.3f%7.3f ", x[i*m+j].re, x[i*m+j].im);
  return (0);
}

int cmatinv( complex x[], int n, int m, double space[])
{
/* x[n*m] ... m>=n
*/
  int i,j,k, *irow=(int*) space;
  double xmaxsize, ee=1e-20;
  complex xmax, t,t1;

  FOR(i,n) {
      xmaxsize = 0.;
      for (j=i; j<n; j++) {
         if ( xmaxsize < csize (x[j*m+i])) {
              xmaxsize = csize (x[j*m+i]);
              xmax = x[j*m+i];
              irow[i] = j;
         }
      }
      if (xmaxsize < ee)  {
			
			#ifdef MPI
				fprintf (stderr, "\nDet goes to zero at %8d!\t\n", i+1);
			#else
				fprintf (stderr, "\nDet goes to zero at %8d!\t\n", i+1);
			#endif
          return(-1);
      }
      if (irow[i] != i) {
          FOR(j,m) {
               t = x[i*m+j];
               x[i*m+j] = x[irow[i]*m+j];
               x[ irow[i]*m+j] = t;
          }
      }
      t = cdiv (compl(1,0), x[i*m+i]);
      FOR(j,n) {
          if (j == i) continue;
          t1 = cby (t,x[j*m+i]);
          FOR(k,m) x[j*m+k] = cminus (x[j*m+k], cby(t1,x[i*m+k]));
          x[j*m+i] = cfactor (t1, -1);
      }
      FOR(j,m)  x[i*m+j] = cby (x[i*m+j], t);
      x[i*m+i] = t;
  }                        
  for (i=n-1; i>=0; i--) {
       if (irow[i] == i) continue;
       FOR(j,n) {
          t = x[j*m+i];
          x[j*m+i] = x[j*m+irow[i]];
          x[ j*m+irow[i]] = t;
       }
  }
  return (0);
}


void balance(double mat[], int n,int *low, int *hi, double scale[])
{
/* Balance a matrix for calculation of eigenvalues and eigenvectors
*/
   double c,f,g,r,s;
   int i,j,k,l,done;
       /* search for rows isolating an eigenvalue and push them down */
   for (k = n - 1; k >= 0; k--) {
       for (j = k; j >= 0; j--) {
           for (i = 0; i <= k; i++) {
               if (i != j && fabs(mat[pos(j,i,n)]) != 0) break;
           }

           if (i > k) {
               scale[k] = j;

               if (j != k) {
                   for (i = 0; i <= k; i++) {
                      c = mat[pos(i,j,n)];
                      mat[pos(i,j,n)] = mat[pos(i,k,n)];
                      mat[pos(i,k,n)] = c;
                   }

                   for (i = 0; i < n; i++) {
                      c = mat[pos(j,i,n)];
                      mat[pos(j,i,n)] = mat[pos(k,i,n)];
                      mat[pos(k,i,n)] = c;
                   }
               }
               break;
           }
       }
       if (j < 0) break;
   }

   /* search for columns isolating an eigenvalue and push them left */

   for (l = 0; l <= k; l++) {
       for (j = l; j <= k; j++) {
           for (i = l; i <= k; i++) {
               if (i != j && fabs(mat[pos(i,j,n)]) != 0) break;
           }
           if (i > k) {
               scale[l] = j;
               if (j != l) {
                   for (i = 0; i <= k; i++) {
                      c = mat[pos(i,j,n)];
                      mat[pos(i,j,n)] = mat[pos(i,l,n)];
                      mat[pos(i,l,n)] = c;
                   }

                   for (i = l; i < n; i++) {
                      c = mat[pos(j,i,n)];
                      mat[pos(j,i,n)] = mat[pos(l,i,n)];
                      mat[pos(l,i,n)] = c;
                   }
               }

               break;
           }
       }

       if (j > k) break;
   }

   *hi = k;
   *low = l;

   /* balance the submatrix in rows l through k */

   for (i = l; i <= k; i++) {
       scale[i] = 1;
   }

   do {
       for (done = 1,i = l; i <= k; i++) {
           for (c = 0,r = 0,j = l; j <= k; j++) {
               if (j != i) {
                   c += fabs(mat[pos(j,i,n)]);
                   r += fabs(mat[pos(i,j,n)]);
               }
           }

           if (c != 0 && r != 0) {
               g = r / BASE;
               f = 1;
               s = c + r;

               while (c < g) {
                   f *= BASE;
                   c *= BASE * BASE;
               }

               g = r * BASE;

               while (c >= g) {
                   f /= BASE;
                   c /= BASE * BASE;
               }

               if ((c + r) / f < 0.95 * s) {
                   done = 0;
                   g = 1 / f;
                   scale[i] *= f;

                   for (j = l; j < n; j++) {
                       mat[pos(i,j,n)] *= g;
                   }

                   for (j = 0; j <= k; j++) {
                       mat[pos(j,i,n)] *= f;
                   }
               }
           }
       }
   } while (!done);
}


/*
 * Transform back eigenvectors of a balanced matrix
 * into the eigenvectors of the original matrix
 */
void unbalance(int n,double vr[],double vi[], int low, int hi, double scale[])
{
   int i,j,k;
   double tmp;

   for (i = low; i <= hi; i++) {
       for (j = 0; j < n; j++) {
           vr[pos(i,j,n)] *= scale[i];
           vi[pos(i,j,n)] *= scale[i];
       }
   }

   for (i = low - 1; i >= 0; i--) {
       if ((k = (int)scale[i]) != i) {
           for (j = 0; j < n; j++) {
               tmp = vr[pos(i,j,n)];
               vr[pos(i,j,n)] = vr[pos(k,j,n)];
               vr[pos(k,j,n)] = tmp;

               tmp = vi[pos(i,j,n)];
               vi[pos(i,j,n)] = vi[pos(k,j,n)];
               vi[pos(k,j,n)] = tmp;       
           }
       }
   }

   for (i = hi + 1; i < n; i++) {
       if ((k = (int)scale[i]) != i) {
           for (j = 0; j < n; j++) {
               tmp = vr[pos(i,j,n)];
               vr[pos(i,j,n)] = vr[pos(k,j,n)];
               vr[pos(k,j,n)] = tmp;

               tmp = vi[pos(i,j,n)];
               vi[pos(i,j,n)] = vi[pos(k,j,n)];
               vi[pos(k,j,n)] = tmp;       
           }
       }
   }
}

/*
 * Reduce the submatrix in rows and columns low through hi of real matrix mat to
 * Hessenberg form by elementary similarity transformations
 */
void elemhess(int job,double mat[],int n,int low,int hi, double vr[],
             double vi[], int work[])
{
/* work[n] */
   int i,j,m;
   double x,y;

   for (m = low + 1; m < hi; m++) {
       for (x = 0,i = m,j = m; j <= hi; j++) {
           if (fabs(mat[pos(j,m-1,n)]) > fabs(x)) {
               x = mat[pos(j,m-1,n)];
               i = j;
           }
       }

       if ((work[m] = i) != m) {
           for (j = m - 1; j < n; j++) {
              y = mat[pos(i,j,n)];
              mat[pos(i,j,n)] = mat[pos(m,j,n)];
              mat[pos(m,j,n)] = y;
           }

           for (j = 0; j <= hi; j++) {
              y = mat[pos(j,i,n)];
              mat[pos(j,i,n)] = mat[pos(j,m,n)];
              mat[pos(j,m,n)] = y;
           }
       }

       if (x != 0) {
           for (i = m + 1; i <= hi; i++) {
               if ((y = mat[pos(i,m-1,n)]) != 0) {
                   y = mat[pos(i,m-1,n)] = y / x;

                   for (j = m; j < n; j++) {
                       mat[pos(i,j,n)] -= y * mat[pos(m,j,n)];
                   }

                   for (j = 0; j <= hi; j++) {
                       mat[pos(j,m,n)] += y * mat[pos(j,i,n)];
                   }
               }
           }
       }
   }
   if (job) {
      for (i=0; i<n; i++) {
         for (j=0; j<n; j++) {
            vr[pos(i,j,n)] = 0.0; vi[pos(i,j,n)] = 0.0;
         }
         vr[pos(i,i,n)] = 1.0;
      }

      for (m = hi - 1; m > low; m--) {
         for (i = m + 1; i <= hi; i++) {
            vr[pos(i,m,n)] = mat[pos(i,m-1,n)];
         }

        if ((i = work[m]) != m) {
           for (j = m; j <= hi; j++) {
              vr[pos(m,j,n)] = vr[pos(i,j,n)];
              vr[pos(i,j,n)] = 0.0;
           }
           vr[pos(i,m,n)] = 1.0;
        }
     }
  }
}

/*
 * Calculate eigenvalues and eigenvectors of a real upper Hessenberg matrix
 * Return 1 if converges successfully and 0 otherwise
 */
 
int realeig(int job,double mat[],int n,int low, int hi, double valr[],
     double vali[], double vr[],double vi[])
{
  complex v;
  double p=0,q=0,r=0,s=0,t,w,x,y,z=0,ra,sa,norm,eps;
  int niter,en,i,j,k,l,m;
  double precision = pow((double)BASE,(double)(1-nDIGITS));

  eps = precision;
  for (i=0; i<n; i++) {
     valr[i]=0.0;
     vali[i]=0.0;
  }
     /* store isolated roots and calculate norm */
  for (norm = 0,i = 0; i < n; i++) {
     for (j = MAX(0,i-1); j < n; j++) {
        norm += fabs(mat[pos(i,j,n)]);
     }
     if (i < low || i > hi) valr[i] = mat[pos(i,i,n)];
  }
  t = 0;
  en = hi;

  while (en >= low) {
     niter = 0;
     for (;;) {

      /* look for single small subdiagonal element */

        for (l = en; l > low; l--) {
           s = fabs(mat[pos(l-1,l-1,n)]) + fabs(mat[pos(l,l,n)]);
           if (s == 0) s = norm;
           if (fabs(mat[pos(l,l-1,n)]) <= eps * s) break;
        }

        /* form shift */

        x = mat[pos(en,en,n)];

        if (l == en) {            /* one root found */
           valr[en] = x + t;
           if (job) mat[pos(en,en,n)] = x + t;
           en--;
           break;
        }

        y = mat[pos(en-1,en-1,n)];
        w = mat[pos(en,en-1,n)] * mat[pos(en-1,en,n)];

        if (l == en - 1) {               /* two roots found */
           p = (y - x) / 2;
           q = p * p + w;
           z = sqrt(fabs(q));
           x += t;
           if (job) {
              mat[pos(en,en,n)] = x;
              mat[pos(en-1,en-1,n)] = y + t;
           }
           if (q < 0) {               /* complex pair */
              valr[en-1] = x+p;
              vali[en-1] = z;
              valr[en] = x+p;
              vali[en] = -z;
           }
           else {                     /* real pair */
              z = (p < 0) ? p - z : p + z;
              valr[en-1] = x + z;
              valr[en] = (z == 0) ? x + z : x - w / z;
              if (job) {
                 x = mat[pos(en,en-1,n)];
                 s = fabs(x) + fabs(z);
                 p = x / s;
                 q = z / s;
                 r = sqrt(p*p+q*q);
                 p /= r;
                 q /= r;
                 for (j = en - 1; j < n; j++) {
                    z = mat[pos(en-1,j,n)];
                    mat[pos(en-1,j,n)] = q * z + p *
                    mat[pos(en,j,n)];
                    mat[pos(en,j,n)] = q * mat[pos(en,j,n)] - p*z;
                 }
                 for (i = 0; i <= en; i++) {
                    z = mat[pos(i,en-1,n)];
                    mat[pos(i,en-1,n)] = q * z + p * mat[pos(i,en,n)];
                    mat[pos(i,en,n)] = q * mat[pos(i,en,n)] - p*z;
                 }
                 for (i = low; i <= hi; i++) {
                    z = vr[pos(i,en-1,n)];
                    vr[pos(i,en-1,n)] = q*z + p*vr[pos(i,en,n)];
                    vr[pos(i,en,n)] = q*vr[pos(i,en,n)] - p*z;
                 }
              }
           }
           en -= 2;
           break;
        }
        if (niter == MAXITER) return(-1);
        if (niter != 0 && niter % 10 == 0) {
           t += x;
           for (i = low; i <= en; i++) mat[pos(i,i,n)] -= x;
           s = fabs(mat[pos(en,en-1,n)]) + fabs(mat[pos(en-1,en-2,n)]);
           x = y = 0.75 * s;
           w = -0.4375 * s * s;
        }
        niter++;
          /* look for two consecutive small subdiagonal elements */
        for (m = en - 2; m >= l; m--) {
           z = mat[pos(m,m,n)];
           r = x - z;
           s = y - z;
           p = (r * s - w) / mat[pos(m+1,m,n)] + mat[pos(m,m+1,n)];
           q = mat[pos(m+1,m+1,n)] - z - r - s;
           r = mat[pos(m+2,m+1,n)];
           s = fabs(p) + fabs(q) + fabs(r);
           p /= s;
           q /= s;
           r /= s;
           if (m == l || fabs(mat[pos(m,m-1,n)]) * (fabs(q)+fabs(r)) <=
               eps * (fabs(mat[pos(m-1,m-1,n)]) + fabs(z) +
               fabs(mat[pos(m+1,m+1,n)])) * fabs(p)) break;
        }
        for (i = m + 2; i <= en; i++) mat[pos(i,i-2,n)] = 0;
        for (i = m + 3; i <= en; i++) mat[pos(i,i-3,n)] = 0;
            /* double QR step involving rows l to en and columns m to en */
        for (k = m; k < en; k++) {
           if (k != m) {
              p = mat[pos(k,k-1,n)];
              q = mat[pos(k+1,k-1,n)];
              r = (k == en - 1) ? 0 : mat[pos(k+2,k-1,n)];
              if ((x = fabs(p) + fabs(q) + fabs(r)) == 0) continue;
              p /= x;
              q /= x;
              r /= x;
           }
           s = sqrt(p*p+q*q+r*r);
           if (p < 0) s = -s;
           if (k != m) {
              mat[pos(k,k-1,n)] = -s * x;
           }
           else if (l != m) {
              mat[pos(k,k-1,n)] = -mat[pos(k,k-1,n)];
           }
           p += s;
           x = p / s;
           y = q / s;
           z = r / s;
           q /= p;
           r /= p;
               /* row modification */
           for (j = k; j <= (!job ? en : n-1); j++){
              p = mat[pos(k,j,n)] + q * mat[pos(k+1,j,n)];
              if (k != en - 1) {
                 p += r * mat[pos(k+2,j,n)];
                 mat[pos(k+2,j,n)] -= p * z;
              }
              mat[pos(k+1,j,n)] -= p * y;
              mat[pos(k,j,n)] -= p * x;
           }
           j = MIN(en,k+3);
             /* column modification */
           for (i = (!job ? l : 0); i <= j; i++) {
              p = x * mat[pos(i,k,n)] + y * mat[pos(i,k+1,n)];
              if (k != en - 1) {
                 p += z * mat[pos(i,k+2,n)];
                 mat[pos(i,k+2,n)] -= p*r;
              }
              mat[pos(i,k+1,n)] -= p*q;
              mat[pos(i,k,n)] -= p;
           }
           if (job) {            /* accumulate transformations */
              for (i = low; i <= hi; i++) {
                 p = x * vr[pos(i,k,n)] + y * vr[pos(i,k+1,n)];
                 if (k != en - 1) {
                    p += z * vr[pos(i,k+2,n)];
                    vr[pos(i,k+2,n)] -= p*r;
                 }
                 vr[pos(i,k+1,n)] -= p*q;
                 vr[pos(i,k,n)] -= p;
              }
           }
        }
     }
  }

  if (!job) return(0);
  if (norm != 0) {
      /* back substitute to find vectors of upper triangular form */
     for (en = n-1; en >= 0; en--) {
        p = valr[en];
        if ((q = vali[en]) < 0) {           /* complex vector */
           m = en - 1;
           if (fabs(mat[pos(en,en-1,n)]) > fabs(mat[pos(en-1,en,n)])) {
              mat[pos(en-1,en-1,n)] = q / mat[pos(en,en-1,n)];
              mat[pos(en-1,en,n)] = (p - mat[pos(en,en,n)]) /
                    mat[pos(en,en-1,n)];
           }
           else {
              v = cdiv(compl(0.0,-mat[pos(en-1,en,n)]),
                   compl(mat[pos(en-1,en-1,n)]-p,q));
              mat[pos(en-1,en-1,n)] = v.re;
              mat[pos(en-1,en,n)] = v.im;
           }
           mat[pos(en,en-1,n)] = 0;
           mat[pos(en,en,n)] = 1;
           for (i = en - 2; i >= 0; i--) {
              w = mat[pos(i,i,n)] - p;
              ra = 0;
              sa = mat[pos(i,en,n)];
              for (j = m; j < en; j++) {
                 ra += mat[pos(i,j,n)] * mat[pos(j,en-1,n)];
                 sa += mat[pos(i,j,n)] * mat[pos(j,en,n)];
              }
              if (vali[i] < 0) {
                 z = w;
                 r = ra;
                 s = sa;
              }
              else {
                 m = i;
                 if (vali[i] == 0) {
                    v = cdiv(compl(-ra,-sa),compl(w,q));
                    mat[pos(i,en-1,n)] = v.re;
                    mat[pos(i,en,n)] = v.im;
                 }
                 else {                     /* solve complex equations */
                    x = mat[pos(i,i+1,n)];
                    y = mat[pos(i+1,i,n)];
                    v.re = (valr[i]- p)*(valr[i]-p) + vali[i]*vali[i] - q*q;
                    v.im = (valr[i] - p)*2*q;
                    if ((fabs(v.re) + fabs(v.im)) == 0) {
                       v.re = eps * norm * (fabs(w) +
                               fabs(q) + fabs(x) + fabs(y) + fabs(z));
                    }
                    v = cdiv(compl(x*r-z*ra+q*sa,x*s-z*sa-q*ra),v);
                    mat[pos(i,en-1,n)] = v.re;
                    mat[pos(i,en,n)] = v.im;
                    if (fabs(x) > fabs(z) + fabs(q)) {
                       mat[pos(i+1,en-1,n)] = 
                            (-ra - w * mat[pos(i,en-1,n)] +
                            q * mat[pos(i,en,n)]) / x;
                       mat[pos(i+1,en,n)] = (-sa - w * mat[pos(i,en,n)] -
                            q * mat[pos(i,en-1,n)]) / x;
                    }
                    else {
                       v = cdiv(compl(-r-y*mat[pos(i,en-1,n)],
                            -s-y*mat[pos(i,en,n)]),compl(z,q));
                       mat[pos(i+1,en-1,n)] = v.re;
                       mat[pos(i+1,en,n)] = v.im;
                    }
                 }
              }
           }
        }
        else if (q == 0) {                            /* real vector */
           m = en;
           mat[pos(en,en,n)] = 1;
           for (i = en - 1; i >= 0; i--) {
              w = mat[pos(i,i,n)] - p;
              r = mat[pos(i,en,n)];
              for (j = m; j < en; j++) {
                 r += mat[pos(i,j,n)] * mat[pos(j,en,n)];
              }
              if (vali[i] < 0) {
                 z = w;
                 s = r;
              }
              else {
                 m = i;
                 if (vali[i] == 0) {
                    if ((t = w) == 0) t = eps * norm;
                    mat[pos(i,en,n)] = -r / t;
                 }
                 else {           /* solve real equations */
                    x = mat[pos(i,i+1,n)];
                    y = mat[pos(i+1,i,n)];
                    q = (valr[i] - p) * (valr[i] - p) + vali[i]*vali[i];
                    t = (x * s - z * r) / q;
                    mat[pos(i,en,n)] = t;
                    if (fabs(x) <= fabs(z)) {
                       mat[pos(i+1,en,n)] = (-s - y * t) / z;
                    }
                    else {
                       mat[pos(i+1,en,n)] = (-r - w * t) / x;
                    }
                 }
              }
           }
        }
     }
            /* vectors of isolated roots */
     for (i = 0; i < n; i++) {
        if (i < low || i > hi) {
           for (j = i; j < n; j++) {
              vr[pos(i,j,n)] = mat[pos(i,j,n)];
           }
        }
     }
      /* multiply by transformation matrix */

     for (j = n-1; j >= low; j--) {
        m = MIN(j,hi);
        for (i = low; i <= hi; i++) {
           for (z = 0,k = low; k <= m; k++) {
              z += vr[pos(i,k,n)] * mat[pos(k,j,n)];
           }
           vr[pos(i,j,n)] = z;
        }
     }
  }
   /* rearrange complex eigenvectors */
  for (j = 0; j < n; j++) {
     if (vali[j] != 0) {
        for (i = 0; i < n; i++) {
           vi[pos(i,j,n)] = vr[pos(i,j+1,n)];
           vr[pos(i,j+1,n)] = vr[pos(i,j,n)];
           vi[pos(i,j+1,n)] = -vi[pos(i,j,n)];
        }
        j++;
     }
  }
  return(0);
}
							/* End of Eigen */
/***********************************************************************/






/************************** gammasCalculate **********************************************/
/* Discrete Gamma's from an alpha (hetereogeneous rate) 
Purpose: sample mean (or median) values from a discrete gamma distribution
		 given an alpha and a number of equally-probable categories
*/

#define POINTGAMMA(prob_d,alpha_d,beta_d) 		PointChi2(prob_d,2.0*(alpha_d))/(2.0*(beta_d))

int		DiscreteGamma (double *rK, double alfa_d, double beta_d, int K_d, int median);
double 	IncompleteGamma (double x_d, double alpha_d, double LnGamma_alpha);
double 	PointChi2 (double prob_d, double v_d);
double	LnGamma (double alp);
double 	PointNormal (double prob_d);


int gammasCalculate (double alpha_d, int numCategories)
{
	/*int			i;*/

	DiscreteGamma (gammaRates, alpha_d, alpha_d, numCategories, 0); /* 0=means 1=medians */

	/*fprintf (stderr, "alpha = %.2f\nnumCategories = %d\nrates =",alpha_d, numCategories);
	for (i=0; i<numCategories; i++)
		 fprintf (stderr, " %.4f", gammaRates[i]);*/
	
	return(0);
}


// code below from Ziheng Yang's PAML (if I remember well)

/*---------------------------------------------------------------------------------
|
|  DiscreteGamma
|
|  Discretization of gamma distribution with equal proportions in each
|  category.
|
---------------------------------------------------------------------------------*/
int DiscreteGamma (double *rK, double alfa_d, double beta_d, int K_d, int median)

{

	int 			i;
	double 		gap05 = 1.0/(2.0*K_d), t, factor = alfa_d/beta_d*K_d, lnga1;

	if (median) 
		{
		for (i=0; i<K_d; i++) 
			rK[i] = POINTGAMMA((i*2.0+1.0)*gap05, alfa_d, beta_d);
		for (i=0,t=0; i<K_d; i++) 
			t += rK[i];
		for (i=0; i<K_d; i++)    
			rK[i] *= factor / t;
		}
	else 
		{
		lnga1 = LnGamma(alfa_d+1);
		/* calculate the points in the gamma distribution */
		for (i=0; i<K_d-1; i++) 
			rK[i] = POINTGAMMA((i+1.0)/K_d, alfa_d, beta_d);
		/* calculate the cumulative values */
		for (i=0; i<K_d-1; i++) 
			rK[i] = IncompleteGamma(rK[i] * beta_d, alfa_d + 1.0, lnga1);
		rK[K_d-1] = 1.0;
		/* calculate the relative values and rescale */
		for (i=K_d-1; i>0; i--)
			{
			rK[i] -= rK[i-1];
			rK[i] *= factor;
			}
		rK[0] *= factor;
		}

	return (0);
	
}


/*---------------------------------------------------------------------------------
|
|  IncompleteGamma
|
|  Returns the incomplete gamma ratio I(x_d,alpha) where x_d is the upper
|  limit of the integration and alpha is the shape parameter. Returns (-1)
|  if in error.  
|
|  Bhattacharjee, G. P. 1970. The incomplete gamma integral. Applied
|     Statistics, 19:285-287 (AS32)
|
---------------------------------------------------------------------------------*/
double IncompleteGamma (double x_d, double alpha_d, double LnGamma_alpha)

{

	int 			i;
	double 		p = alpha_d, g = LnGamma_alpha,
					accurate = 1e-8, overflow = 1e30,
					factor, gin = 0.0, rn = 0.0, a = 0.0, b = 0.0, an = 0.0, 
					dif = 0.0, term = 0.0, pn[6];

	if (x_d == 0.0) 
		return (0.0);
	if (x_d < 0 || p <= 0) 
		return (-1.0);

	factor = exp(p*log(x_d)-x_d-g);  
	if (x_d>1 && x_d>=p) 
		goto l30;
	gin = 1.0; 
	term = 1.0; 
	rn = p;
	l20:
		rn++;
		term *= x_d/rn;  
		gin += term;
		if (term > accurate) 
			goto l20;
		gin *= factor/p;
		goto l50;
	l30:
		a = 1.0-p;  
		b = a+x_d+1.0; 
		term = 0.0;
		pn[0] = 1.0; 
		pn[1] = x_d; 
		pn[2] = x_d+1; 
		pn[3] = x_d*b;
		gin = pn[2]/pn[3];
	l32:
		a++; 
		b += 2.0; 
		term++;  
		an = a*term;
		for (i=0; i<2; i++) 
			pn[i+4] = b*pn[i+2]-an*pn[i];
		if (pn[5] == 0) 
			goto l35;
		rn = pn[4]/pn[5];  
		dif = fabs(gin-rn);
		if (dif>accurate) 
			goto l34;
		if (dif<=accurate*rn) 
			goto l42;
	l34:
		gin = rn;
	l35:
		for (i=0; i<4; i++) 
			pn[i] = pn[i+2];
		if (fabs(pn[4]) < overflow) 
			goto l32;
		for (i=0; i<4; i++) 
			pn[i] /= overflow;
		goto l32;
	l42:
		gin = 1.0-factor*gin;
	l50:
		return (gin);

}



/*---------------------------------------------------------------------------------
|
|  PointChi2
|
|  Returns z so that Prob(x < z) = prob where x is Chi2 distributed with df=v. 
|  Returns -1 if in error.  0.000002 < prob < 0.999998.
|
---------------------------------------------------------------------------------*/
double PointChi2 (double prob_d, double v_d)

{

	double 		e = 0.5e-6, aa = 0.6931471805, p = prob_d, g,
					xx, c, ch, a = 0.0, q = 0.0, p1 = 0.0, p2 = 0.0, t = 0.0, 
					x_d = 0.0, b = 0.0, s1, s2, s3, s4, s5, s6;

	if (p < 0.000002 || p > 0.999998 || v_d <= 0.0) 
		return (-1.0);
	g = LnGamma (v_d/2.0);
	xx = v_d/2.0;  
	c = xx - 1.0;
	if (v_d >= -1.24*log(p)) 
		goto l1;
	ch = pow((p*xx*exp(g+xx*aa)), 1.0/xx);
	if (ch-e<0) 
		return (ch);
	goto l4;
	l1:
		if (v_d > 0.32) 
			goto l3;
		ch = 0.4;  
		a = log(1.0-p);
	l2:
		q = ch; 
		p1 = 1.0+ch*(4.67+ch); 
		p2 = ch*(6.73+ch*(6.66+ch));
		t = -0.5+(4.67+2.0*ch)/p1 - (6.73+ch*(13.32+3.0*ch))/p2;
		ch -= (1.0-exp(a+g+0.5*ch+c*aa)*p2/p1)/t;
		if (fabs(q/ch-1.0)-0.01 <= 0.0) 
			goto l4;
		else                      
			goto l2;
	l3: 
		x_d = PointNormal (p);
		p1 = 0.222222/v_d;  
		ch = v_d*pow((x_d*sqrt(p1)+1.0-p1), 3.0);
		if (ch > 2.2*v_d+6.0) 
			ch = -2.0*(log(1.0-p)-c*log(0.5*ch)+g);
	l4:
		q = ch;  
		p1 = 0.5*ch;
		if ((t = IncompleteGamma (p1, xx, g)) < 0.0) 
			{
			#ifdef MPI
			if (rank==root)
				fprintf (stderr, "%s  Error: Problem in PointChi2", " ");
			#else
			fprintf (stderr, "%s  Error: Problem in PointChi2", " ");
			#endif
			return (-1.0);
			}
		p2 = p-t;
		t = p2*exp(xx*aa+g+p1-c*log(ch));  
		b = t/ch; 
		a = 0.5*t-b*c;
		s1 = (210.0+a*(140.0+a*(105.0+a*(84.0+a*(70.0+60.0*a))))) / 420.0;
		s2 = (420.0+a*(735.0+a*(966.0+a*(1141.0+1278.0*a))))/2520.0;
		s3 = (210.0+a*(462.0+a*(707.0+932.0*a)))/2520.0;
		s4 = (252.0+a*(672.0+1182.0*a)+c*(294.0+a*(889.0+1740.0*a)))/5040.0;
		s5 = (84.0+264.0*a+c*(175.0+606.0*a)) / 2520.0;
		s6 = (120.0+c*(346.0+127.0*c)) / 5040.0;
		ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
		if (fabs(q/ch-1.0) > e) 
			goto l4;
		return (ch);

}



/*---------------------------------------------------------------------------------
|
|  LnGamma
|
|  Calculates the log of the gamma function. The Gamma function is equal
|  to:
|
|     Gamma(alp) = {integral from 0 to infinity} t^{alp-1} e^-t dt
|
|  The result is accurate to 10 decimal places. Stirling's formula is used
|  for the central polynomial part of the procedure.
|
|  Pike, M. C. and I. D. Hill. 1966. Algorithm 291: Logarithm of the gamma
|     function. Communications of the Association for Computing
|     Machinery, 9:684.
|     
---------------------------------------------------------------------------------*/
double LnGamma (double alp)

{

	double 		x_d = alp, f = 0.0, z;
	
	if (x_d < 7) 
		{
		f = 1.0; 
		z = x_d-1.0;
		while (++z < 7.0) 
			f *= z;
		x_d = z;  
		f = -log(f);
		}
	z = 1.0 / (x_d*x_d);
	return (f + (x_d-0.5)*log(x_d) - x_d + 0.918938533204673 + 
			(((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z +
			0.083333333333333)/x_d); 

}


/*---------------------------------------------------------------------------------
|
|  PointNormal
|
|  Returns z so That Prob{x<z} = prob where x ~ N(0,1) and
|  (1e-12) < prob < 1-(1e-12). Returns (-9999) if in error. 
|
|  Odeh, R. E. and J. O. Evans. 1974. The percentage points of the normal
|    distribution. Applied Statistics, 22:96-97 (AS70).
|
|  Newer methods:
|
|  Wichura, M. J. 1988. Algorithm AS 241: The percentage points of the
|     normal distribution. 37:477-484.
|  Beasley, JD & S. G. Springer. 1977. Algorithm AS 111: The percentage
|     points of the normal distribution. 26:118-121.
|
---------------------------------------------------------------------------------*/
double PointNormal (double prob_d)

{

	double 		a0 = -0.322232431088, a1 = -1.0, a2 = -0.342242088547, a3 = -0.0204231210245,
 					a4 = -0.453642210148e-4, b0 = 0.0993484626060, b1 = 0.588581570495,
 					b2 = 0.531103462366, b3 = 0.103537752850, b4 = 0.0038560700634,
 					y, z = 0, p = prob_d, p1;

	p1 = (p<0.5 ? p : 1-p);
	if (p1<1e-20) 
	  return (-9999);
	y = sqrt (log(1/(p1*p1)));  
	z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
	
	return (p<0.5 ? -z : z);

}
						/* end of gamma's rate */
/**************************************************************************/









/******************************** CalcIndividualGi *********************************/
/* Calculates G, the number of potential locations in an active ancestral gamete 
  that a recombination event can take place */

static int CalcIndividualGi (int who, TreeNode *nodes, int *activeGametes, int numNuc, int *S_MRCA, int sizeNode)
{
	int			Gi, i, v, vv, position;
	int 		lessValue, MaxValue;
	/*int			*stud;*/
	TreeNode	*p;
	TreeSegment *s;
	lessValue = -1;
	MaxValue = -1;
	Gi = 0;
	position = 0;
	i = v = vv = 0;
		
	p = nodes + activeGametes[who];
	
	if (sizeNode != p->numSegNode)
		fprintf (fpmpi, "\n\nWarning in CalcIndividualGi, sizeNode != p->numSegNode");
	if (noisy == 4) 
		fprintf (fpmpi, "\nComputing Gi, node %d with %d fragment(s):",p->index,p->numSegNode);
	for (position = 0; position < p->numSegNode; position++)
		{
		s = segments + post(position,p->index,distance);
		if (noisy == 4)
			fprintf (fpmpi, "\nFragment %d, s->sStart = %d and s->sEnd = %d",s->sIndex, s->sStart, s->sEnd);
		}
		
		
	/*if (doCodonModel == YES)*/ /* Special to codon model. No breakpoint in codons. */
	/*	{
		stud = (int *)calloc(numSites,(long) sizeof(int));
		if (!stud)
			{			
			fprintf (fpmpi, "Could not allocate stud (%lu bytes)\n", numSites  * (long) sizeof(int));
			exit (1);
			}
		
		for (v = 1; v < numSites; v++) 
			stud[v-1] = v*3+1;*/ 	/* The only ones possible breakpoints (only between codons). "stud" is an array with the possible breakpoints */
				
	/*	for (i = 1; i < numNuc; i++)
			{
			vv = 0;
			for (position = 0; position < p->numSegNode; position++)
				{
				s = segments + post(position,p->index,distance);

				if (i > s->sStart && i <= s->sEnd && S_MRCA[i] > 1)*/	/* a segment with just one site cannot recombine */
	/*				{
					for (v = 1; v < numSites; v++)
						if (i == stud[v-1])
							{
							Gi++;
							vv++;
							break;
							}
					if (vv > 0)
						break;
					}
				}
			}
		free (stud);
		return Gi; 
		}*/
	

	/* breakp between two fragments - trapped material*/
	for (position = 0; position < p->numSegNode; position++)
		{
		s = segments + post(position,p->index,distance);
			
		if (position == 0)
			{
			lessValue = s->sStart;
			MaxValue = s->sEnd;
			}			


		if (s->sStart < lessValue)
			{
			lessValue = s->sStart;
			}

		if (s->sEnd >= MaxValue)
			{
			MaxValue = s->sEnd;
			}
		}



	for (i = 1; i <= /*numSites*/numNuc; i++)
		{
		/*fprintf (fpmpi, "\n ******************* i = %d, MaxValue = %d, lessValue = %d, S_MRCA[i] = %d \n ", i, MaxValue, lessValue, S_MRCA[i]); */
		if (i > lessValue && i <= MaxValue && S_MRCA[i] > 1)	/* a segment with just one site cannot recombine */
			{
			Gi++;
			}
		}
	return Gi; 

}







/******************************** IsValidBreakSite *********************************/
/* returns whether a site is potentially recombining: it has 
	  ancestral material non-MRCA before and after it */

static int 	IsValidBreakSite (int *activeGametes, TreeNode *nodes, int whichInd, int whichSite, int *S_MRCA) /* This function returns YES or NO about if the site is a potencial breakpoint (if it's in the ancestral material) or no */
{
	int			k, a, i, position;
	int 		lessValue, MaxValue;
	TreeNode	*p;
	TreeSegment *s;
	
	position = k = a = i = 0;
	lessValue = -1;
	MaxValue = -1;
	
	p = nodes + activeGametes[whichInd];
	for (position = 0; position < p->numSegNode; position++)
		{
		s = segments + post(position,p->index,distance);
		
		/* the breakpoint must to be in one or more segments & no S_MRCA (= 1) */
		if (s->sStart < whichSite && s->sEnd >= whichSite)	/* Ok, the breakpoint is in segments of this node */
			{
			k++;
			break;
			}
		}
	
	/* breakp between two fragments - trapped material*/
	if (k == 0)
		{
		for (position = 0; position < p->numSegNode; position++)
			{
			s = segments + post(position,p->index,distance);
			
			if (position == 0)
				{
				lessValue = s->sStart;
				MaxValue = s->sEnd;
				}			


			if (s->sStart < lessValue)
				{
				lessValue = s->sStart;
				}

			if (s->sEnd >= MaxValue)
				{
				MaxValue = s->sEnd;
				}
			}

		if (lessValue < whichSite && MaxValue >= whichSite)
			{
			k++;
			}	
		}

	

	if (whichSite < 0) /* Cheking */
		{
		k = 0;		
		fprintf (fpmpi, "\n Warning in IsValidBreakSite function, breakp = %d ", whichSite);
		exit (-1);
		}
	if (k != 0)
		if (S_MRCA[whichSite] > 1)	/* Ok, the breakpoint is not in the MRCA */
			a++;
			
								
	if (k > 0 && a > 0)
		return YES;
	else
		return NO;
}





/*********************************** CountsForExpNumRec **********************************/
/* Checks whether a recombination event should be counted in the expected number of 
  recombinations E(R). For E(R) count only events with breakpoints as X|X, X|0 or 0|X, 
  where X is ancestral material (i.e, 1), but that did not find yet its MRCA  */                                                                  

static int CountsForExpNumRec (int *activeGametes, int whichInd, int whichSite, TreeNode *nodes, int *S_MRCA, int sizeNode)
{
	int			k, j, a, i, position, FirstFragmentPos;
	TreeNode	*p;
	TreeSegment *s;
	int *initialVector, *endVector;

	position = 0;
	k = j = a = i = FirstFragmentPos = 0;
	
	initialVector = (int *)calloc(sizeNode,(long) sizeof(int));
	if (!initialVector)
		{
		fprintf (fpmpi, "Could not allocate initialVector (%lu bytes)\n", sizeNode *(long) sizeof(int));
		exit (1);
		}
	endVector = (int *)calloc(sizeNode,(long) sizeof(int));
	if (!endVector)
		{
		fprintf (fpmpi, "Could not allocate endVector (%lu bytes)\n", sizeNode *(long) sizeof(int));
		exit (1);
		}
	


	p = nodes + activeGametes[whichInd];
	if (p->numSegNode != sizeNode)
		fprintf (fpmpi, "Warning in CountsForExpNumRec function, p->numSegNode != sizeNode");		
		
	for (position = 0; position < p->numSegNode; position++)
		{
		s = segments + post(position,p->index,distance);
		initialVector[position] = s->sStart;
		endVector[position] = s->sEnd;
		}
		


	/* the breakpoint must to be in one or more segments (at first or into the segment) & no MRCA (> 1) */
	for (position = 0; position < p->numSegNode; position++)
		{
		if ((endVector[position] + 1) == whichSite)
			FirstFragmentPos++;
		}

	for (position = 0; position < p->numSegNode; position++)
		{		
		if (initialVector[position] < whichSite && endVector[position] >= whichSite)	/* Ok, the breakpoint is in these segments */
			k++;
		if (initialVector[position] == whichSite && FirstFragmentPos > 0) /* Ok, the breakpoint is in these segments */
			k++;
		}


	free (initialVector);
	free (endVector);
	if (whichSite < 0) /* Cheking */
		{
		k = 0;
		fprintf (fpmpi, "\n Warning in CountsForExpNumRec function ");
		exit (-1);
		}	
	if (k > 0)
		if (S_MRCA[whichSite] > 1 || S_MRCA[whichSite+1] > 1)	/* Ok, the breakpoint isn't S_MRCA */
			a++;
			
	if (k > 0 && a > 0)
		return YES;
	else
		return NO;
}




/********************* WhichNuc ************************/
/* Returns character representation for nucleotides */

char WhichNuc (int nucIeotide)
{
	if (nucIeotide == 0)
		return ('A');
	else if (nucIeotide == 1)
		return ('C');
	else if (nucIeotide == 2)
		return ('G');
	else if (nucIeotide == 3)
		return ('T');
	else
		return ('X');
}



/********************* WhichNucNumber ************************/
/* Returns character representation for nucleotides */

int WhichNucNumber (char siteLetter)
{
	if (siteLetter == 'A')
		return (0);
	else if (siteLetter == 'C')
		return (1);
	else if (siteLetter == 'G')
		return (2);
	else if (siteLetter == 'T')
		return (3);
	else
		return (4);
}







/***************************** PrintSequences *******************************/
/* Prints sequences to alignment file */

static void PrintSequences (/*int replicate*/) 
{
	int		 i, j, m, dem, outgroupLabel;
	TreeNode	*f;
	
	dem = 0;
	outgroupLabel = 0;


	if (thereisOutgroup == YES)
		{
		fprintf(fpAlignment,"%d %d \n", numSequences+1, numNuc);
		/*fprintf(fpAlignment,"Dataset_%d %d %d \n", replicate+1, numSequences+1, numNuc);*/
		}
	else
		{
		/*fprintf(fpAlignment,"Dataset_%d %d %d \n", replicate+1, numSequences, numNuc);*/
		fprintf(fpAlignment,"%d %d \n", numSequences, numNuc);
		}
	
	if (doMigration == NO)
		{
		if(thereisOutgroup == YES)
			{
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 2)
					{
					outgroupLabel = f->label;
					break;
					}
				}
	
			for (i=0; i<numSequences+1; i++)
				{
				if (i == numSequences)
					fprintf(fpAlignment, "outgroup  ");
				else 
					fprintf (fpAlignment,"seq%05d  ", i+1);
				
				if (i == numSequences) /* outgroup */
					{
					for (j=1; j<=numNuc; j++)
						{
						fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(outgroupLabel,j,numNuc)]));
						}
					}
				else /* tip nodes */
					{
					for (j=1; j<=numNuc; j++)
						{
						fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
						}
					}
				fprintf (fpAlignment,"\n");
				}	
			fprintf (fpAlignment,"\n");
			/*fprintf(stderr,"\n sequence %d go with nodo of label %d \n", i+1, i+1);*/
			}
		else
			{
			for (i=0; i<numSequences; i++)
				{
				fprintf (fpAlignment,"seq%05d  ", i+1);

				for (j=1; j<=numNuc; j++)
					fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
				fprintf (fpAlignment,"\n");
				}	
			fprintf (fpAlignment,"\n");
			}
		}
	else /* migration */
		{
		if(thereisOutgroup == YES)
			{	

			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 2)
					{
					outgroupLabel = f->label;
					break;
					}
				}
			/*fprintf (stderr, "\n outgroupLabel = %d \n", outgroupLabel);*/


			for (i=0; i<numSequences+1; i++)
				{

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					if ((f->label == i) && (f->index < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}
				
				if (i == numSequences)
					fprintf(fpAlignment, "outgrp_p0 ");
				else 
					fprintf (fpAlignment,"s%05d_p%d ", i+1, dem);

				
				if (i == numSequences) /* outgroup */
					{
					for (j=1; j<=numNuc; j++)
						{
						fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(outgroupLabel,j,numNuc)]));
						}
					}
				else /* tip nodes */
					{
					for (j=1; j<=numNuc; j++)
						{
						fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
						}
					}

				fprintf (fpAlignment,"\n");
				}	
			fprintf (fpAlignment,"\n");
			}
		else
			{
			for (i=0; i<numSequences; i++)
				{
				/*for (m = 0; m < numNodex; m++)
					{
					f = nodex + m;
					if ((f->label == i) && (f->NetIndex < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}*/

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					if ((f->label == i) && (f->index < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}

				fprintf (fpAlignment,"s%05d_p%d ", i+1, dem);
				
				for (j=1; j<=numNuc; j++)
					{
					fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
					}
				
				fprintf (fpAlignment,"\n");
				}	
			fprintf (fpAlignment,"\n");
			}
		}
}
		


/***************************** PrintAncestralSequences *******************************/
/* Prints ancestral sequences to alignment file */

static void PrintAncestralSequences (/*int replicate*/)   
{
	int		 i, j, a, m, dem, outgroupLabel, rootLabel;
	TreeNode	*f;
	dem = 0;
	outgroupLabel = rootLabel = 0;

	if (numRE == 0) /* There are NOT recombinations */
		{
		if (thereisOutgroup == YES)
			{
			fprintf(fpAlignment,"%d %d\n", 2*numSequences, numNuc);
			/*fprintf(fpAlignment,"Dataset_%d %d %d\n", replicate+1, 2*numSequences, numNuc);*/
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 2)
					outgroupLabel = f->label;
				if (f->class == 5)
					rootLabel = f->label;
				}
			}
		else
			{
			fprintf(fpAlignment,"%d %d\n", 2*numSequences-1, numNuc);
			/*fprintf(fpAlignment,"Dataset_%d %d %d\n",replicate+1, 2*numSequences-1, numNuc);*/
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 5)
					rootLabel = f->label;
				}
			}


		
		if (doMigration == NO)
			{
			if(thereisOutgroup == YES)
				{	
				for (i=0; i<2*numSequences; i++)
					{
					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,"seq%05d  ", i+1);
					else if (i == numSequences)
						fprintf(fpAlignment, "outgroup  ");
					else if (i == 2*numSequences-1)
						fprintf (fpAlignment,"root      ");
					else
						fprintf (fpAlignment,"anc%05d  ", i+1);
		

					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
							}
						}
					else if (i == numSequences) /* is outgroup */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(outgroupLabel,j,numNuc)]));
							}
						}
					else if (i == 2*numSequences-1) /* is root */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
							}
						}
					else  /* is ancestral */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i-1,j,numNuc)]));
							}
						}

					fprintf (fpAlignment,"\n");
					}
				fprintf (fpAlignment,"\n");
				}
			else
				{
				for (i=0; i<2*numSequences-1; i++)
					{
					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,"seq%05d  ", i+1);
					else if (i == 2*numSequences-2)
						fprintf (fpAlignment,"root      ");
					else
						fprintf (fpAlignment,"anc%05d  ", i+1);


					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
							}
						}
					else if (i == 2*numSequences-2) /* is root */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
							}
						}
					else /* is ancestral */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
							}
						}
					fprintf (fpAlignment,"\n");
					}	
				fprintf (fpAlignment,"\n");
				}
			}
		else /* migration */
			{
			if(thereisOutgroup == YES)
				{	
				for (i=0; i<2*numSequences; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodes + m;
						if ((f->label == i) && (f->index <= numSequences*2-2))
							{
							dem = f->indexOldMigPop;
							break;
							}
						}*/

					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;

						if (i < numSequences) /* tip */
							{
							if ((f->label == i) && (f->index <= numSequences*2-2)) 
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == 2*numSequences-1) /* root */
							{
							if (f->class == 5)
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == numSequences) /* outgroup */
							{
							dem = 0;
							break;
							}
						else /* ancestral */
							{
							if (f->label == i-1)
								{
								dem = f->indexOldMigPop;
								if (f->class != 4)
									{
									fprintf (stderr, "\n Warning in PrintAncestralSequences_C. f->label = %d, f->class = %d \n", f->label, f->class);
									exit(-1);
									}
								break;
								}
							}
						}
					
					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,"s%05d_p%d ", i+1, dem);
					else if (i == numSequences)
						fprintf(fpAlignment, "outgrp_p0 ");
					else if (i == 2*numSequences-1)
						fprintf (fpAlignment,"root_p%d   ", dem);
					else
						fprintf (fpAlignment,"a%05d_p%d ", i+1, dem);
		

					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
							}
						}
					else if (i == numSequences) /* is outgroup */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(outgroupLabel,j,numNuc)]));
							}
						}
					else if (i == 2*numSequences-1) /* is root */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
							}
						}
					else  /* is ancestral */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i-1,j,numNuc)]));
							}
						}
					fprintf (fpAlignment,"\n");
					/*fprintf(stderr," \n\n sequence %d with label(in times file) %d", i+1, i+1);*/
					}	
				fprintf (fpAlignment,"\n");
				}
			else
				{
				for (i=0; i<2*numSequences-1; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodes + m;
						if ((f->label == i) && (f->index <= numSequences*2-2))
							{
							dem = f->indexOldMigPop;
							break;
							}
						}*/

					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;

						if (i < numSequences) /* tip */
							{
							if ((f->label == i) && (f->index <= numSequences*2-2)) 
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == 2*numSequences-2) /* root */
							{
							if (f->class == 5)
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else /* ancestral */
							{
							if (f->label == i)
								{
								dem = f->indexOldMigPop;
								if (f->class != 4)
									{
									fprintf (stderr, "\n Warning in PrintAncestralSequences_C. f->label = %d, f->class = %d \n", f->label, f->class);
									exit(-1);
									}
								break;
								}
							}
						}
						
					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,"s%05d_p%d ", i+1, dem);
					else if (i == 2*numSequences-2)
						fprintf (fpAlignment,"root_p%d   ", dem);
					else
						fprintf (fpAlignment,"a%05d_p%d ", i+1, dem);

					
					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
							}
						}
					else if (i == 2*numSequences-2) /* is root */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
							}
						}
					else /* is ancestral */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
							}
						}

					fprintf (fpAlignment,"\n");
					}	
				fprintf (fpAlignment,"\n");
				}
			}
		}
	else /* There are recombinations */
		{
		if (thereisOutgroup == YES)
			fprintf(fpAlignment,"%d %d \n", numSequences+2, numNuc);
			/*fprintf(fpAlignment,"Dataset_%d %d %d \n", replicate+1, numSequences+1, numNuc);*/
		else
			/*fprintf(fpAlignment,"Dataset_%d %d %d \n", replicate+1, numSequences, numNuc);*/
			fprintf(fpAlignment,"%d %d \n", numSequences+1, numNuc);
		
		if (doMigration == NO)
			{
			if(thereisOutgroup == YES)
				{	

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 2)
						outgroupLabel = f->label;
					if (f->class == 5)
						rootLabel = f->label;
					}		


				for (i=0; i<2*numSequences; i++)
					{
					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,"seq%05d  ", i+1);
						a++;
						}
					if (i == numSequences)
						{
						fprintf(fpAlignment, "outgroup  ");
						a++;
						}
					if (i == 2*numSequences-1)
						{
						fprintf (fpAlignment,"root      ");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == numSequences) /* is outgroup */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(outgroupLabel,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-1) /* is root */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}

						/*fprintf(stderr," \n\n sequence %d with label(in times file) %d", i+1, i+1);*/
						}
					}	
				fprintf (fpAlignment,"\n");
				}
			else
				{
				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 5)
						rootLabel = f->label;
					}


				for (i=0; i<2*numSequences-1; i++)
					{
					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,"seq%05d  ", i+1);
						a++;
						}
					if (i == 2*numSequences-2)
						{
						fprintf (fpAlignment,"root      ");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{

						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-2) /* is root */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}


						}
					}	
				fprintf (fpAlignment,"\n");
				}
			}
		else /* migration */
			{
			if(thereisOutgroup == YES)
				{	
				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 2)
						outgroupLabel = f->label;
					if (f->class == 5)
						rootLabel = f->label;
					}	





				for (i=0; i<2*numSequences; i++)
					{					
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodes + m;
						if (f->index < numSequences*/ /*|| i == 2*numSequences-1*//*)
							if ( f->label == i)
								{
								dem = f->indexOldMigPop;*/
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
							/*	break;
								}
						}*/
					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;
						if (f->index < numSequences /*|| i == 2*numSequences-1*/)
							if (f->label == i)
								{
								dem = f->indexOldMigPop;
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								break;
								}
						}

					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,"s%05d_p%d ", i+1, dem);
						a++;
						}
					if (i == numSequences)
						{
						fprintf(fpAlignment, "outgrp_p0 ");
						a++;
						}
					if (i == 2*numSequences-1)
						{
						fprintf (fpAlignment,"root      ");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{

						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == numSequences) /* is outgroup */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(outgroupLabel,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-1) /* is root */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}

						}
					}	
				fprintf (fpAlignment,"\n");
				}
			else
				{

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 5)
						rootLabel = f->label;
					}


				for (i=0; i<2*numSequences-1; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodes + m;
						if (f->index < numSequences*/ /*|| i == 2*numSequences-1*//*)
							if ( f->label == i)
								{
								dem = f->indexOldMigPop;
								break;
								}
						}*/
					
					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;
						if (f->index < numSequences /*|| i == 2*numSequences-1*/)
							if ( f->label == i)
								{
								dem = f->indexOldMigPop;
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								break;
								}
						}


					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,"s%05d_p%d ", i+1, dem);
						a++;
						}
					if (i == 2*numSequences-2)
						{
						fprintf (fpAlignment,"root      ");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-2) /* is root */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}

						}
					}	
				fprintf (fpAlignment,"\n");
				}

			}
		}
}
		







/***************************** PrintSequences for codon Model *******************************/
/* Prints sequences to alignment file in phylip sequential format */

static void PrintSequences_C (/*int replicate*/)
{
	int		 i, j, k, m, dem, outgroupLabel;
	char codon[3];
	TreeNode	*f;
	dem = 0;
	outgroupLabel = 0;
	
	if (thereisOutgroup == YES)
		{
		fprintf(fpAlignment,"%d %d \n", numSequences+1, numNuc);
		/*fprintf(fpAlignment,"Dataset_%d %d %d \n", replicate+1, numSequences+1, numNuc);*/
		}
	else
		{
		/*fprintf(fpAlignment,"Dataset_%d %d %d \n", replicate+1, numSequences, numNuc);*/
		fprintf(fpAlignment,"%d %d \n", numSequences, numNuc);
		}
	
	if (doMigration == NO)
		{
		if(thereisOutgroup == YES)
			{
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 2)
					{
					outgroupLabel = f->label;
					break;
					}
				}
	
			for (i=0; i<numSequences+1; i++)
				{
				if (i == numSequences)
					fprintf(fpAlignment, "outgroup  ");
				else 
					fprintf (fpAlignment,"seq%05d  ", i+1);

				if (i == numSequences) /* outgroup */
					{
					for (j=1; j<=numSites; j++)
						{
						if (matrixC[pos(outgroupLabel,j,numSites)] > 60 || matrixC[pos(outgroupLabel,j,numSites)] > 60)
							{
							fprintf (stderr, "\n stop codon5 \n");
							exit(-1);
							}
						number_to_codon(matrixC[pos(outgroupLabel,j,numSites)], codon);
						for (k = 0; k < 3; k++)
							fprintf (fpAlignment, "%c", codon[k]);
						}
					}
				else /* tip nodes */
					{
					for (j=1; j<=numSites; j++)
						{
						if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
							{
							fprintf (stderr, "\n stop codon5 \n");
							exit(-1);
							}
						number_to_codon(matrixC[pos(i,j,numSites)], codon);
						for (k = 0; k < 3; k++)
							fprintf (fpAlignment, "%c", codon[k]);
						}
					}
				fprintf (fpAlignment,"\n");
				/*fprintf(stderr,"\n sequence %d go with nodo of label %d \n", i+1, i+1);*/
				}	
			fprintf (fpAlignment,"\n");
			}
		else
			{
			for (i=0; i<numSequences; i++)
				{
				fprintf (fpAlignment,"seq%05d  ", i+1);

				for (j=1; j<=numSites; j++)
					{
					if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
						{
						fprintf (stderr, "\n stop codon6 \n");
						exit(-1);
						}
					number_to_codon(matrixC[pos(i,j,numSites)], codon);
					for (k = 0; k < 3; k++)
						fprintf (fpAlignment, "%c", codon[k]);
					}
				fprintf (fpAlignment,"\n");
				}	
			fprintf (fpAlignment,"\n");
			}
		}
	else /* migration */
		{
		if(thereisOutgroup == YES)
			{
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 2)
					{
					outgroupLabel = f->label;
					break;
					}
				}
			/*fprintf (stderr, "\n outgroupLabel = %d \n", outgroupLabel);*/
			

			for (i=0; i<numSequences+1; i++)
				{
				/*for (m = 0; m < numNodex; m++)
					{
					f = nodex + m;
					if ((f->label == i) && (f->NetIndex < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}*/
				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					if ((f->label == i) && (f->index < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}		
		
				if (i == numSequences)
					fprintf(fpAlignment, "outgrp_p0 ");
				else 
					fprintf (fpAlignment,"s%05d_p%d ", i+1, dem);
				
				if (i == numSequences) /* outgroup */
					{
					for (j=1; j<=numSites; j++)
						{
						if (matrixC[pos(outgroupLabel,j,numSites)] > 60 || matrixC[pos(outgroupLabel,j,numSites)] > 60)
							{
							fprintf (stderr, "\n stop codon7 \n");
							exit(-1);
							}
						number_to_codon(matrixC[pos(outgroupLabel,j,numSites)], codon);
						for (k = 0; k < 3; k++)
							fprintf (fpAlignment, "%c", codon[k]);
						}
					}
				else /* tip nodes */
					{
					for (j=1; j<=numSites; j++)
						{
						if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
							{
							fprintf (stderr, "\n stop codon7 \n");
							exit(-1);
							}
						number_to_codon(matrixC[pos(i,j,numSites)], codon);
						for (k = 0; k < 3; k++)
							fprintf (fpAlignment, "%c", codon[k]);
						}
					}

				fprintf (fpAlignment,"\n");
				/*fprintf(stderr,"\n sequence %d go with node label(times file label) %d \n", i+1, i+1);*/
				}	
			fprintf (fpAlignment,"\n");
			}
		else
			{
			for (i=0; i<numSequences; i++)
				{
				/*for (m = 0; m < numNodex; m++)
					{
					f = nodex + m;
					if ((f->label == i) && (f->NetIndex < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}*/
				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					if ((f->label == i) && (f->index < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}
					
				fprintf (fpAlignment,"s%05d_p%d ", i+1, dem);

				for (j=1; j<=numSites; j++)
					{
					if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
						{
						fprintf (stderr, "\n stop codon8 \n");
						exit(-1);
						}
					number_to_codon(matrixC[pos(i,j,numSites)], codon);
					for (k = 0; k < 3; k++)
						fprintf (fpAlignment, "%c", codon[k]);
					}
				fprintf (fpAlignment,"\n");
				}	
			fprintf (fpAlignment,"\n");
			}
		}
}
		


/***************************** PrintAncestralSequences for Codon Model *******************************/
/* Prints ancestral sequences to alignment file for codon Model */

static void PrintAncestralSequences_C (/*int replicate*/)
{
	int		 i, j, k, a, m, dem, outgroupLabel, rootLabel;
	char codon[3];
	TreeNode	*f;
	dem = 0;
	outgroupLabel = rootLabel = 0;	

	if (numRE == 0) /* There are NOT recombinations */
		{
		if (thereisOutgroup == YES)
			{
			fprintf(fpAlignment,"%d %d\n", 2*numSequences, numNuc);
			/*fprintf(fpAlignment,"Dataset_%d %d %d\n", replicate+1, 2*numSequences, numNuc);*/
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 2)
					outgroupLabel = f->label;
				if (f->class == 5)
					rootLabel = f->label;
				}
			}
		else
			{
			fprintf(fpAlignment,"%d %d\n", 2*numSequences-1, numNuc);
			/*fprintf(fpAlignment,"Dataset_%d %d %d\n",replicate+1, 2*numSequences-1, numNuc);*/
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 5)
					rootLabel = f->label;
				}
			}
		

		if (doMigration == NO)
			{
			if (thereisOutgroup == YES)
				{	
				for (i=0; i<2*numSequences; i++)
					{
					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,"seq%05d  ", i+1);
					else if (i == numSequences)
						fprintf(fpAlignment, "outgroup  ");
					else if (i == 2*numSequences-1)
						fprintf (fpAlignment,"root      ");
					else
						fprintf (fpAlignment,"anc%05d  ", i+1);


					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else if (i == numSequences) /* is outgroup */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(outgroupLabel,j,numSites)] > 60 || matrixC[pos(outgroupLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(outgroupLabel,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else if (i == 2*numSequences-1) /* is root */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else  /* is ancestral */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i-1,j,numSites)] > 60 || matrixC[pos(i-1,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i-1,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}

					fprintf (fpAlignment,"\n");
					}	
				fprintf (fpAlignment,"\n");
				}
			else
				{
				for (i=0; i<2*numSequences-1; i++)
					{
					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,"seq%05d  ", i+1);
					else if (i == 2*numSequences-2)
						fprintf (fpAlignment,"root      ");
					else
						fprintf (fpAlignment,"anc%05d  ", i+1);

					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon11 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i,j,numSites)], codon);
							for	(k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else if (i == 2*numSequences-2) /* is root */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon11 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
							for	(k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else /* is ancestral */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon11 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i,j,numSites)], codon);
							for	(k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}

					fprintf (fpAlignment,"\n");
					}	
				fprintf (fpAlignment,"\n");
				}
			}
		else /* migration */
			{
			if (thereisOutgroup == YES)
				{	
				for (i=0; i<2*numSequences; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodex + m;
						if ((f->label == i) && (f->NetIndex <= numSequences*2-2))
							{
							dem = f->indexOldMigPop;
							break;
							}
						}*/
					
					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;

						if (i < numSequences) /* tip */
							{
							if ((f->label == i) && (f->index <= numSequences*2-2)) 
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == 2*numSequences-1) /* root */
							{
							if (f->class == 5)
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == numSequences) /* outgroup */
							{
							dem = 0;
							break;
							}
						else /* ancestral */
							{
							if (f->label == i-1)
								{
								dem = f->indexOldMigPop;
								if (f->class != 4)
									{
									fprintf (stderr, "\n Warning in PrintAncestralSequences_C. f->label = %d, f->class = %d \n", f->label, f->class);
									exit(-1);
									}
								break;
								}
							}
						}



					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,"s%05d_p%d ", i+1, dem);
					else if (i == numSequences)
						fprintf(fpAlignment, "outgrp_p0 ");
					else if (i == 2*numSequences-1)
						fprintf (fpAlignment,"root_p%d   ",dem);
					else
						fprintf (fpAlignment,"a%05d_p%d ", i+1, dem);
		


					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else if (i == numSequences) /* is outgroup */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(outgroupLabel,j,numSites)] > 60 || matrixC[pos(outgroupLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(outgroupLabel,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else if (i == 2*numSequences-1) /* is root */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else  /* is ancestral */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i-1,j,numSites)] > 60 || matrixC[pos(i-1,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i-1,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}


					fprintf (fpAlignment,"\n");
					}	
				fprintf (fpAlignment,"\n");
				}
			else
				{
				for (i=0; i<2*numSequences-1; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodex + m;
						if ((f->label == i) && (f->NetIndex <= numSequences*2-2))
							{
							dem = f->indexOldMigPop;
							break;
							}
						}*/
					
					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;

						if (i < numSequences) /* tip */
							{
							if ((f->label == i) && (f->index <= numSequences*2-2)) 
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == 2*numSequences-2) /* root */
							{
							if (f->class == 5)
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else /* ancestral */
							{
							if (f->label == i)
								{
								dem = f->indexOldMigPop;
								if (f->class != 4)
									{
									fprintf (stderr, "\n Warning in PrintAncestralSequences_C. f->label = %d, f->class = %d \n", f->label, f->class);
									exit(-1);
									}
								break;
								}
							}
						}



					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,"s%05d_p%d ", i+1, dem);
					else if (i == 2*numSequences-2)
						fprintf (fpAlignment,"root_p%d   ", dem);
					else
						fprintf (fpAlignment,"a%05d_p%d ", i+1, dem);


					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon11 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i,j,numSites)], codon);
							for	(k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else if (i == 2*numSequences-2) /* is root */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon11 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
							for	(k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else /* is ancestral */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon11 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i,j,numSites)], codon);
							for	(k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}

					fprintf (fpAlignment,"\n");

					}	
				}
			fprintf (fpAlignment,"\n");
			}
		}
	else /* There are recombinations */
		{
		if (thereisOutgroup == YES)
			fprintf(fpAlignment,"%d %d\n", numSequences+2, numNuc);
		else
			fprintf(fpAlignment,"%d %d\n", numSequences+1, numNuc);
		
		if (doMigration == NO)
			{
			if (thereisOutgroup == YES)
				{
				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 2)
						outgroupLabel = f->label;
					if (f->class == 5)
						rootLabel = f->label;
					}				
	
				for (i=0; i<2*numSequences; i++)
					{
					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,"seq%05d  ", i+1);
						a++;
						}
					if (i == numSequences)
						{
						fprintf(fpAlignment, "outgroup  ");
						a++;
						}
					if (i == 2*numSequences-1)
						{
						fprintf (fpAlignment,"root      ");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon14 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(i,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == numSequences) /* is outgroup */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(outgroupLabel,j,numSites)] > 60 || matrixC[pos(outgroupLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon14 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(outgroupLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-1) /* is root */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon14 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}

						}
					}	
				fprintf (fpAlignment,"\n");
				}
			else
				{

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 5)
						rootLabel = f->label;
					}

				for (i=0; i<2*numSequences-1; i++)
					{
					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,"seq%05d  ", i+1);
						a++;
						}
					if (i == 2*numSequences-2)
						{
						fprintf (fpAlignment,"root      ");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon15 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(i,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-2) /* is root */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon15 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}

						}
					}	
				fprintf (fpAlignment,"\n");
				}
			}
		else /* migration */
			{
			if (thereisOutgroup == YES)
				{	

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 2)
						outgroupLabel = f->label;
					if (f->class == 5)
						rootLabel = f->label;
					}	


				for (i=0; i<2*numSequences; i++)
					{

					/*for (m = 0; m < numNodex; m++)
						{
						f = nodex + m;
						if (f->NetIndex < numSequences*/ /*|| i == 2*numSequences-1*//*
							if ( f->label == i)
								{
								dem = f->indexOldMigPop;*/
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								/*break;
								}
						}*/

					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;
						if (f->index < numSequences /*|| i == 2*numSequences-1*/)
							if (f->label == i)
								{
								dem = f->indexOldMigPop;
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								break;
								}
						}



					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,"s%05d_p%d ", i+1, dem);
						a++;
						}
					if (i == numSequences)
						{
						fprintf(fpAlignment, "outgrp_p0 ");
						a++;
						}
					if (i == 2*numSequences-1)
						{
						fprintf (fpAlignment,"root      ");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{

						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon16 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(i,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == numSequences) /* is outgroup */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(outgroupLabel,j,numSites)] > 60 || matrixC[pos(outgroupLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon16 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(outgroupLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-1) /* is root */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon16 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}

						}
					}	
				fprintf (fpAlignment,"\n");
				}
			else
				{


				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 5)
						rootLabel = f->label;
					}

				for (i=0; i<2*numSequences-1; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodex + m;
						if (f->NetIndex < numSequences*/ /*|| i == 2*numSequences-1*//*)
							if (f->label == i)
								{
								dem = f->indexOldMigPop;*/
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								/*break;
								}
						}*/
					
					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;
						if (f->index < numSequences /*|| i == 2*numSequences-1*/)
							if ( f->label == i)
								{
								dem = f->indexOldMigPop;
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								break;
								}
						}



					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,"s%05d_p%d ", i+1, dem);
						a++;
						}
					if (i == 2*numSequences-2)
						{
						fprintf (fpAlignment,"root      ");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon17 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(i,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-2) /* is root */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon17 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}


						}
					}	
				fprintf (fpAlignment,"\n");
				}
			}
		}
}






/***************************** PrintOutMRCAFiles (GMRCA) for nucleotide Model *******************************/

/* Prints just GMRCA files for nucleotide Model */
static void PrintOutMRCAFiles (/*int replicate*/)
{
int		 i, j, a, m, dem, outgroupLabel, rootLabel;
	TreeNode	*f;
	dem = 0;
	outgroupLabel = rootLabel = 0;

	if (numRE == 0) /* There are NOT recombinations */
		{
		if (thereisOutgroup == YES)
			{
			fprintf(fpMRCAprint,"%d %d\n", 1, numNuc);
			/*fprintf(fpMRCAprint,"Dataset_%d %d %d\n", replicate+1, 2*numSequences, numNuc);*/
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 2)
					outgroupLabel = f->label;
				if (f->class == 5)
					rootLabel = f->label;
				}
			}
		else
			{
			fprintf(fpMRCAprint,"%d %d\n", 1, numNuc);
			/*fprintf(fpMRCAprint,"Dataset_%d %d %d\n",replicate+1, 2*numSequences-1, numNuc);*/
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 5)
					rootLabel = f->label;
				}
			}


		
		if (doMigration == NO)
			{
			if(thereisOutgroup == YES)
				{	
				for (i=0; i<2*numSequences; i++)
					{
					if (i == 2*numSequences-1)
						fprintf (fpMRCAprint,"root      ");
							

					if (i == 2*numSequences-1) /* is root */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpMRCAprint, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
							}
						}
					
					//fprintf (fpMRCAprint,"\n");
					}
			//	fprintf (fpMRCAprint,"\n");
				}
			else
				{
				for (i=0; i<2*numSequences-1; i++)
					{
					if (i == 2*numSequences-2)
						fprintf (fpMRCAprint,"root      ");
					

					if (i == 2*numSequences-2) /* is root */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpMRCAprint, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
							}
						}
					
				//	fprintf (fpMRCAprint,"\n");
					}	
			//	fprintf (fpMRCAprint,"\n");
				}
			}
		else /* migration */
			{
			if(thereisOutgroup == YES)
				{	
				for (i=0; i<2*numSequences; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodes + m;
						if ((f->label == i) && (f->index <= numSequences*2-2))
							{
							dem = f->indexOldMigPop;
							break;
							}
						}*/

					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;

						if (i < numSequences) /* tip */
							{
							if ((f->label == i) && (f->index <= numSequences*2-2)) 
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == 2*numSequences-1) /* root */
							{
							if (f->class == 5)
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == numSequences) /* outgroup */
							{
							dem = 0;
							break;
							}
						else /* ancestral */
							{
							if (f->label == i-1)
								{
								dem = f->indexOldMigPop;
								if (f->class != 4)
									{
									fprintf (stderr, "\n Warning in PrintAncestralSequences_C. f->label = %d, f->class = %d \n", f->label, f->class);
									exit(-1);
									}
								break;
								}
							}
						}
					
				
					if (i == 2*numSequences-1)
						fprintf (fpMRCAprint,"root      ");
					if (i == 2*numSequences-1) /* is root */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpMRCAprint, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
							}
						}
					
				//	fprintf (fpMRCAprint,"\n");
					/*fprintf(stderr," \n\n sequence %d with label(in times file) %d", i+1, i+1);*/
					}	
			//	fprintf (fpMRCAprint,"\n");
				}
			else
				{
				for (i=0; i<2*numSequences-1; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodes + m;
						if ((f->label == i) && (f->index <= numSequences*2-2))
							{
							dem = f->indexOldMigPop;
							break;
							}
						}*/

					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;

						if (i < numSequences) /* tip */
							{
							if ((f->label == i) && (f->index <= numSequences*2-2)) 
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == 2*numSequences-2) /* root */
							{
							if (f->class == 5)
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else /* ancestral */
							{
							if (f->label == i)
								{
								dem = f->indexOldMigPop;
								if (f->class != 4)
									{
									fprintf (stderr, "\n Warning in PrintAncestralSequences_C. f->label = %d, f->class = %d \n", f->label, f->class);
									exit(-1);
									}
								break;
								}
							}
						}
						
					
					if (i == 2*numSequences-2)
						fprintf (fpMRCAprint,"root      ");
										
					else if (i == 2*numSequences-2) /* is root */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpMRCAprint, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
							}
						}
				
				//	fprintf (fpMRCAprint,"\n");
					}	
			//	fprintf (fpMRCAprint,"\n");
				}
			}
		}
	else /* There are recombinations */
		{
		if (thereisOutgroup == YES)
			fprintf(fpMRCAprint,"%d %d \n", 1, numNuc);
			/*fprintf(fpMRCAprint,"Dataset_%d %d %d \n", replicate+1, numSequences+1, numNuc);*/
		else
			/*fprintf(fpMRCAprint,"Dataset_%d %d %d \n", replicate+1, numSequences, numNuc);*/
			fprintf(fpMRCAprint,"%d %d \n", 1, numNuc);
		
		if (doMigration == NO)
			{
			if(thereisOutgroup == YES)
				{	

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 2)
						outgroupLabel = f->label;
					if (f->class == 5)
						rootLabel = f->label;
					}		


				for (i=0; i<2*numSequences; i++)
					{
					a = 0;
					if (i < numSequences) /* is tip */
						{
				//		fprintf (fpMRCAprint,"seq%05d  ", i+1);
						a++;
						}
					if (i == numSequences)
						{
				//		fprintf(fpMRCAprint, "outgroup  ");
						a++;
						}
					if (i == 2*numSequences-1)
						{
						fprintf (fpMRCAprint,"root      ");
						a++;
						}
					/*else
						fprintf (fpMRCAprint,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						if (i == 2*numSequences-1) /* is root */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpMRCAprint, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
								}
					//		fprintf (fpMRCAprint,"\n");
							}

						/*fprintf(stderr," \n\n sequence %d with label(in times file) %d", i+1, i+1);*/
						}
					}	
			//	fprintf (fpMRCAprint,"\n");
				}
			else
				{
				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 5)
						rootLabel = f->label;
					}


				for (i=0; i<2*numSequences-1; i++)
					{
					a = 0;
					if (i < numSequences) /* is tip */
						{
					//	fprintf (fpMRCAprint,"seq%05d  ", i+1);
						a++;
						}
					if (i == 2*numSequences-2)
						{
						fprintf (fpMRCAprint,"root      ");
						a++;
						}
					/*else
						fprintf (fpMRCAprint,"anc%05d  ", i+1);*/
					if (a != 0)
						{

						if (i == 2*numSequences-2) /* is root */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpMRCAprint, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
								}
						//	fprintf (fpMRCAprint,"\n");
							}

						}
					}	
			//	fprintf (fpMRCAprint,"\n");
				}
			}
		else /* migration */
			{
			if(thereisOutgroup == YES)
				{	
				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 2)
						outgroupLabel = f->label;
					if (f->class == 5)
						rootLabel = f->label;
					}	

				for (i=0; i<2*numSequences; i++)
					{					
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodes + m;
						if (f->index < numSequences*/ /*|| i == 2*numSequences-1*//*)
							if ( f->label == i)
								{
								dem = f->indexOldMigPop;*/
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
							/*	break;
								}
						}*/
					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;
						if (f->index < numSequences /*|| i == 2*numSequences-1*/)
							if (f->label == i)
								{
								dem = f->indexOldMigPop;
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								break;
								}
						}

					a = 0;
					if (i < numSequences) /* is tip */
						{
					//	fprintf (fpMRCAprint,"s%05d_p%d ", i+1, dem);
						a++;
						}
					if (i == numSequences)
						{
					//	fprintf(fpMRCAprint, "outgrp_p0 ");
						a++;
						}
					if (i == 2*numSequences-1)
						{
						fprintf (fpMRCAprint,"root      ");
						a++;
						}
					/*else
						fprintf (fpMRCAprint,"anc%05d  ", i+1);*/
					if (a != 0)
						{

						if (i == 2*numSequences-1) /* is root */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpMRCAprint, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
								}
						//	fprintf (fpMRCAprint,"\n");
							}

						}
					}	
		//		fprintf (fpMRCAprint,"\n");
				}
			else
				{

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 5)
						rootLabel = f->label;
					}


				for (i=0; i<2*numSequences-1; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodes + m;
						if (f->index < numSequences*/ /*|| i == 2*numSequences-1*//*)
							if ( f->label == i)
								{
								dem = f->indexOldMigPop;
								break;
								}
						}*/
					
					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;
						if (f->index < numSequences /*|| i == 2*numSequences-1*/)
							if ( f->label == i)
								{
								dem = f->indexOldMigPop;
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								break;
								}
						}


					a = 0;
					if (i < numSequences) /* is tip */
						{
					//	fprintf (fpMRCAprint,"s%05d_p%d ", i+1, dem);
						a++;
						}
					if (i == 2*numSequences-2)
						{
						fprintf (fpMRCAprint,"root      ");
						a++;
						}
					/*else
						fprintf (fpMRCAprint,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						
						if (i == 2*numSequences-2) /* is root */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpMRCAprint, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
								}
					//		fprintf (fpMRCAprint,"\n");
							}

						}
					}	
		//		fprintf (fpMRCAprint,"\n");
				}

			}
		}
}







/***************************** PrintOutMRCAFiles_Conc, it prints MRCAs for nucleotide Model *******************************/
/* Prints just MRCA files for nucleotide Model */

static void PrintOutMRCAFiles_Conc ()
{
	char	*MRCAseq_Array;		/* this array will contain the MRCA sequence from inputFile */

	/*char 	MRCAseq_Array[numNuc];*/
	/*int  	codon[3];*/
	int nucPosition, thisLabel, valueNuc, sss;

	nucPosition = 0;	

	MRCAseq_Array = (char *) calloc((numNuc+1), sizeof(char)); 
	if (!MRCAseq_Array)
		{
		#ifdef MPI
			fprintf (stderr, "%d: Could not allocate MRCAseq_Array (%lu bytes)\n", rank, (numNuc+1)  * (long) sizeof(char));
		#else
			fprintf (stderr, "Could not allocate MRCAseq_Array (%lu bytes)\n", (numNuc+1)  * (long) sizeof(char));
		#endif
		exit (-1);
		}


/*	matrix = (int *)calloc(((nextAvailable+1) * (numNuc+1)),(long) sizeof(int));		
		if (!matrix)
			{
			fprintf (fpmpi, "Could not allocate matrix (%lu bytes)\n", ((nextAvailable+1) * (numNuc+1))  * (long) sizeof(int));
			exit (1);
			}
		for (i = 0; i < (nextAvailable+1) * (numNuc+1); i++)
			matrix[i] = -1;	*/







	/* Writing the MRCA sequence in the array: MRCAseq_Array*/
	/*NodesMRCAposit[0]= -1;*/ /* en este vector se guardara: posicion 1 es el label (=index) del nodo MRCA del nuc 1, posicion 2... */
	/*fprintf (fpmpi, "\n\n Vector of MRCA: ");*/
	for (sss = 1; sss <= numNuc; sss++)
		{
		thisLabel = NodesMRCAposit[sss];
		valueNuc = matrix[pos(thisLabel,sss,numNuc)];
			
		/*fprintf (fpmpi, " %d", valueNuc);*/

		MRCAseq_Array[sss] = WhichNuc(valueNuc);
		}
	/*fprintf (fpmpi, "\n\n Nuc del Vector of MRCA: ");*/
		

	/* print in file */
	fprintf(fpConcMRCAprint,"%d %d\n", 1, numNuc);
	fprintf (fpConcMRCAprint,"root      ");
	
	for (sss = 1; sss <= numNuc; sss++)
		{
		/*fprintf (fpmpi, " %c", MRCAseq_Array[sss]);*/
		fprintf (fpConcMRCAprint, "%c", MRCAseq_Array[sss]);
		}

	free (MRCAseq_Array);
}










/***************************** PrintMRCAfiles_C, it prints the GMRCAs for Codon Model *******************************/
/* Prints just GMRCA files for codon Model */

static void PrintOutMRCAFiles_C (/*int replicate*/)
{
	int		 i, j, k, a, m, dem, outgroupLabel, rootLabel;
	char codon[3];
	TreeNode	*f;
	dem = 0;
	outgroupLabel = rootLabel = 0;	

	if (numRE == 0) /* There are NOT recombinations */
		{
		if (thereisOutgroup == YES)
			{
			fprintf(fpMRCAprint,"%d %d\n", 1, numNuc);
			/*fprintf(fpMRCAprint,"Dataset_%d %d %d\n", replicate+1, 2*numSequences, numNuc);*/
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 2)
					outgroupLabel = f->label;
				if (f->class == 5)
					rootLabel = f->label;
				}
			}
		else
			{
			fprintf(fpMRCAprint,"%d %d\n", 1, numNuc);
			/*fprintf(fpMRCAprint,"Dataset_%d %d %d\n",replicate+1, 2*numSequences-1, numNuc);*/
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 5)
					{
					rootLabel = f->label;
					}
				}
			}
		

		if (doMigration == NO)
			{
			if (thereisOutgroup == YES)
				{	
				for (i=0; i<2*numSequences; i++)
					{
					if (i == 2*numSequences-1)
						fprintf (fpMRCAprint,"root      ");
					

					if (i == 2*numSequences-1) /* is root */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpMRCAprint, "%c", codon[k]);
							}
						}

					//fprintf (fpMRCAprint,"\n");
					}	
			//	fprintf (fpMRCAprint,"\n");
				}
			else
				{
				for (i=0; i<2*numSequences-1; i++)
					{
					if (i == 2*numSequences-2)
						fprintf (fpMRCAprint,"root      ");
					
					
					if (i == 2*numSequences-2) /* is root */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon11 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
							for	(k = 0; k < 3; k++)
								fprintf (fpMRCAprint, "%c", codon[k]);
							}
						}

				//	fprintf (fpMRCAprint,"\n");
					}	
			//	fprintf (fpMRCAprint,"\n");
				}
			}
		else /* migration */
			{
			if (thereisOutgroup == YES)
				{	
				for (i=0; i<2*numSequences; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodex + m;
						if ((f->label == i) && (f->NetIndex <= numSequences*2-2))
							{
							dem = f->indexOldMigPop;
							break;
							}
						}*/
					
					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;

						if (i < numSequences) /* tip */
							{
							if ((f->label == i) && (f->index <= numSequences*2-2)) 
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == 2*numSequences-1) /* root */
							{
							if (f->class == 5)
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == numSequences) /* outgroup */
							{
							dem = 0;
							break;
							}
						else /* ancestral */
							{
							if (f->label == i-1)
								{
								dem = f->indexOldMigPop;
								if (f->class != 4)
									{
									fprintf (stderr, "\n Warning in PrintAncestralSequences_C. f->label = %d, f->class = %d \n", f->label, f->class);
									exit(-1);
									}
								break;
								}
							}
						}



					if (i == 2*numSequences-1)
						{
						fprintf (fpMRCAprint,"root      ");
						//fprintf (fpMRCAprint,"root_pop%d      ",dem);
						}


					if (i == 2*numSequences-1) /* is root */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpMRCAprint, "%c", codon[k]);
							}
						}
					
				//	fprintf (fpMRCAprint,"\n");
					}	
			//	fprintf (fpMRCAprint,"\n");
				}
			else
				{
				for (i=0; i<2*numSequences-1; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodex + m;
						if ((f->label == i) && (f->NetIndex <= numSequences*2-2))
							{
							dem = f->indexOldMigPop;
							break;
							}
						}*/
					
					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;

						if (i < numSequences) /* tip */
							{
							if ((f->label == i) && (f->index <= numSequences*2-2)) 
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == 2*numSequences-2) /* root */
							{
							if (f->class == 5)
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else /* ancestral */
							{
							if (f->label == i)
								{
								dem = f->indexOldMigPop;
								if (f->class != 4)
									{
									fprintf (stderr, "\n Warning in PrintAncestralSequences_C. f->label = %d, f->class = %d \n", f->label, f->class);
									exit(-1);
									}
								break;
								}
							}
						}




					
					if (i == 2*numSequences-2)
						{
						fprintf (fpMRCAprint,"root      ");
						//fprintf (fpMRCAprint,"root_pop%d      ",dem);
						}

					
					if (i == 2*numSequences-2) /* is root */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon11 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
							for	(k = 0; k < 3; k++)
								fprintf (fpMRCAprint, "%c", codon[k]);
							}
						}
					
				//	fprintf (fpMRCAprint,"\n");

					}	
				}
			//fprintf (fpMRCAprint,"\n");
			}
		}
	else /* There are recombinations */
		{
		if (thereisOutgroup == YES)
			fprintf(fpMRCAprint,"%d %d\n", 1, numNuc);
		else
			fprintf(fpMRCAprint,"%d %d\n", 1, numNuc);
		
		if (doMigration == NO)
			{
			if (thereisOutgroup == YES)
				{
				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 2)
						outgroupLabel = f->label;
					if (f->class == 5)
						rootLabel = f->label;
					}				
	
				for (i=0; i<2*numSequences; i++)
					{
					a = 0;
					if (i < numSequences) /* is tip */
						{
						//fprintf (fpMRCAprint,"seq%05d  ", i+1);
						a++;
						}
					if (i == numSequences)
						{
						//fprintf(fpMRCAprint, "outgroup  ");
						a++;
						}
					if (i == 2*numSequences-1)
						{
						fprintf (fpMRCAprint,"root      ");
						a++;
						}
					/*else
						fprintf (fpMRCAprint,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						if (i == 2*numSequences-1) /* is root */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon14 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpMRCAprint, "%c", codon[k]);
								}
						//	fprintf (fpMRCAprint,"\n");
							}

						}
					}	
			//	fprintf (fpMRCAprint,"\n");
				}
			else
				{

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 5)
						{
						rootLabel = f->label;
						}
					}

				for (i=0; i<2*numSequences-1; i++)
					{
					a = 0;
					if (i < numSequences) /* is tip */
						{
						//fprintf (fpMRCAprint,"seq%05d  ", i+1);
						a++;
						}
					if (i == 2*numSequences-2)
						{
						fprintf (fpMRCAprint,"root      ");
						a++;
						}
					/*else
						fprintf (fpMRCAprint,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						if (i == 2*numSequences-2) /* is root */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon15 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpMRCAprint, "%c", codon[k]);
								}
						//	fprintf (fpMRCAprint,"\n");
							}

						}
					}	
			//	fprintf (fpMRCAprint,"\n");
				}
			}
		else /* migration */
			{
			if (thereisOutgroup == YES)
				{	

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 2)
						outgroupLabel = f->label;
					if (f->class == 5)
						rootLabel = f->label;
					}	


				for (i=0; i<2*numSequences; i++)
					{

					/*for (m = 0; m < numNodex; m++)
						{
						f = nodex + m;
						if (f->NetIndex < numSequences*/ /*|| i == 2*numSequences-1*//*
							if ( f->label == i)
								{
								dem = f->indexOldMigPop;*/
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								/*break;
								}
						}*/

					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;
						if (f->index < numSequences /*|| i == 2*numSequences-1*/)
							if (f->label == i)
								{
								dem = f->indexOldMigPop;
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								break;
								}
						}



					a = 0;
					if (i < numSequences) /* is tip */
						{
						//fprintf (fpMRCAprint,"seq%05d_pop%d  ", i+1, dem);
						a++;
						}
					if (i == numSequences)
						{
						//fprintf(fpMRCAprint, "outgroup_pop0  ");
						a++;
						}
					if (i == 2*numSequences-1)
						{
						fprintf (fpMRCAprint,"root      ");
						a++;
						}
					/*else
						fprintf (fpMRCAprint,"anc%05d  ", i+1);*/
					if (a != 0)
						{

						if (i == 2*numSequences-1) /* is root */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon16 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpMRCAprint, "%c", codon[k]);
								}
						//	fprintf (fpMRCAprint,"\n");
							}

						}
					}	
			//	fprintf (fpMRCAprint,"\n");
				}
			else
				{


				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 5)
						rootLabel = f->label;
					}

				for (i=0; i<2*numSequences-1; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodex + m;
						if (f->NetIndex < numSequences*/ /*|| i == 2*numSequences-1*//*)
							if (f->label == i)
								{
								dem = f->indexOldMigPop;*/
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								/*break;
								}
						}*/
					
					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;
						if (f->index < numSequences /*|| i == 2*numSequences-1*/)
							if ( f->label == i)
								{
								dem = f->indexOldMigPop;
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								break;
								}
						}


					a = 0;
					if (i < numSequences) /* is tip */
						{
						//fprintf (fpMRCAprint,"seq%05d_pop%d  ", i+1, dem);
						a++;
						}
					if (i == 2*numSequences-2)
						{
						fprintf (fpMRCAprint,"root      ");
						a++;
						}
					/*else
						fprintf (fpMRCAprint,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						
						if (i == 2*numSequences-2) /* is root */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon17 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpMRCAprint, "%c", codon[k]);
								}
						//	fprintf (fpMRCAprint,"\n");
							}


						}
					}	
				//fprintf (fpMRCAprint,"\n");
				}
			}
		}
}



/***************************** PrintMRCAfiles_C, it prints MRCAs for Codon Model *******************************/
/* Prints just MRCA files for codon Model */

static void PrintOutMRCAFiles_C_Conc (/*int replicate*/)
{
	int		 i, k;
	char	*MRCAseq_Array;		/* this array will contain the MRCA sequence from inputFile */

	/*char 	MRCAseq_Array[numNuc];*/
	int  	codon[3];
	int xx, yy, zz, nucPosition, thisLabel, valueNuc, sss;

	nucPosition = 0;	

	MRCAseq_Array = (char *) calloc((numNuc+1), sizeof(char)); 
	if (!MRCAseq_Array)
		{
		#ifdef MPI
			fprintf (stderr, "%d: Could not allocate MRCAseq_Array (%lu bytes)\n", rank, (numNuc+1)  * (long) sizeof(char));
		#else
			fprintf (stderr, "Could not allocate MRCAseq_Array (%lu bytes)\n", (numNuc+1)  * (long) sizeof(char));
		#endif
		exit (-1);
		}

	matrixCnuc = (int *)calloc(((nextAvailable+1) * (numNuc+2)),(long) sizeof(int));		/* New nuc matrix from codon matrix */
	if (!matrixCnuc)
		{
		fprintf (fpmpi, "Could not allocate matrixC (%lu bytes)\n", ((nextAvailable+1) * (numNuc+2))  * (long) sizeof(int));
		exit (1);
		}
	for (i = 0; i < (nextAvailable+1) * (numNuc+2); i++)
		matrixCnuc[i] = -1;



	/*fprintf(fpmpi,"Tamanho max de matrix: nextAvailable = %d \n", nextAvailable);*/ /* From Matrix Codon to Matrix nuc */ /* matrixC[pos(p->label,siteCodon,numSites)]; */
	for (xx=0; xx<=nextAvailable; xx++)
		{
		/*fprintf(fpmpi,"p->Label %d \n", xx);*/
		for (yy=1; yy<=numSites; yy++)
			{
			zz = matrixC[pos(xx,yy,numSites)];
			/*fprintf(fpmpi,"%d ", zz);*/ 				
			

			/* write the nuc */ /* void number_to_codon2(int ind, int out[])*/
			codon[0] = codon[1] = codon[2] = -1;
			if (zz > -1) /* this node do not have that codon */
				{
				number_to_codon2(zz, codon);
				}
			else
				{
				codon[0] = codon[1] = codon[2] = -1;
				}

			for (k = 0; k < 3; k++)
				{
				if (k == 0)
					{
					nucPosition = (yy*3)-2;
					}
				else if (k == 1)
					{
					nucPosition = (yy*3)-1;
					}
				else if (k == 2)
					{
					nucPosition = yy*3;
					}
				else
					{
					fprintf(fpmpi,"\n never here PrintOutMRCAFiles_C_Conc");
					exit(-1);
					}

				matrixCnuc[pos(xx,nucPosition,numNuc)] = codon[k];
				}
			}
		/*fprintf(fpmpi,"\n");*/
		}

	/* For see the new matrix */
	/*fprintf(fpmpi,"Tamanho max de matrix_NUC: nextAvailable = %d \n", nextAvailable);*/
	/*for (xx=0; xx<=nextAvailable; xx++)
		{
		fprintf(fpmpi,"\np->Label %d \n", xx);
		for (yy=1; yy<=numNuc; yy++)
			{
			zz = matrixCnuc[pos(xx,yy,numNuc)];
			fprintf(fpmpi,"%d ", zz); 				
			}
		}*/




	/* Writing the MRCA sequence in the array: MRCAseq_Array*/
	/*NodesMRCAposit[0]= -1;*/ /* en este vector se guardara: posicion 1 es el label (=index) del nodo MRCA del nuc 1, posicion 2... */
	/*fprintf (fpmpi, "\n\n Vector of MRCA: ");*/
	for (sss = 1; sss <= numNuc; sss++)
		{
		thisLabel = NodesMRCAposit[sss];
		valueNuc = matrixCnuc[pos(thisLabel,sss,numNuc)];
			
		/*fprintf (fpmpi, " %d", valueNuc);*/

		MRCAseq_Array[sss] = WhichNuc(valueNuc);
		}
	/*fprintf (fpmpi, "\n\n Nuc del Vector of MRCA: ");*/
		

	/* print in file */
	fprintf(fpConcMRCAprint,"%d %d\n", 1, numNuc);
	fprintf (fpConcMRCAprint,"root      ");
	
	for (sss = 1; sss <= numNuc; sss++)
		{
		/*fprintf (fpmpi, " %c", MRCAseq_Array[sss]);*/
		fprintf (fpConcMRCAprint, "%c", MRCAseq_Array[sss]);
		}

	free (MRCAseq_Array);
	free (matrixCnuc);
}





/***************************** PrintOutGMRCAFiles_Codon_AncestralMat, it prints the GMRCAs for the ancestral material in Codon Models *******************************/
/* Prints just GMRCA files for codon Model */

static void PrintOutGMRCAFiles_Codon_AncestralMat (/*int replicate*/)
{
	int		 j, k, m, rootLabel, Ok;
	char codon[3];
	TreeNode	*f;
	rootLabel = Ok = 0;	


	/*fprintf(fpmpi,"\n AT PrintOutGMRCAFiles_Codon_AncestralMat \n");*/

	fprintf(fpGMRCAancPrint,"%d %d\n", 1, numNuc);
	for (m = 0; m < nextAvailable; m++)
		{
		f = nodes + m;
		/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
		if (f->GMRCA_ancestral == YES)
			{
			rootLabel = f->label;
			/*fprintf(fpmpi,"This node is index = %d, label = %d \n", f->index, f->label);*/
			Ok++;
			}
		}
	if (Ok == 0)
		{
		for (m = 0; m < nextAvailable; m++)
			{
			f = nodes + m;
			/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
			if (f->class == 5)
				{
				rootLabel = f->label;
				/*fprintf(fpmpi,"Sp_This node is index = %d, label = %d \n", f->index, f->label);*/
				Ok++;
				}
			}

		}




	fprintf (fpGMRCAancPrint,"root      ");

	for (j=1; j<=numSites; j++)
		{
		if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
			{
			fprintf (stderr, "\n stop codonX1 \n");
			exit(-1);
			}
		number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
		for (k = 0; k < 3; k++)
			fprintf (fpGMRCAancPrint, "%c", codon[k]);
		}

}












/**** FASTA ****/
/***************************** PrintSequences_FASTA *******************************/
/* Prints sequences to alignment FASTA file */
static void PrintSequences_FASTA (/*int replicate*/) 
{
	int		 i, j, m, dem, outgroupLabel;
	TreeNode	*f;
	
	dem = 0;
	outgroupLabel = 0;


/*	if (thereisOutgroup == YES)
		fprintf(fpAlignment,"%d %d \n", numSequences+1, numNuc);*/
		/*fprintf(fpAlignment,"Dataset_%d %d %d \n", replicate+1, numSequences+1, numNuc);*/
/*	else*/
		/*fprintf(fpAlignment,"Dataset_%d %d %d \n", replicate+1, numSequences, numNuc);*/
/*		fprintf(fpAlignment,"%d %d \n", numSequences, numNuc);*/
	
	if (doMigration == NO)
		{
		if(thereisOutgroup == YES)
			{
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 2)
					{
					outgroupLabel = f->label;
					break;
					}
				}
	
			for (i=0; i<numSequences+1; i++)
				{
				if (i == numSequences)
					fprintf(fpAlignment, ">outgroup\n");
				else 
					fprintf (fpAlignment,">seq%05d\n", i+1);
				
				if (i == numSequences) /* outgroup */
					{
					for (j=1; j<=numNuc; j++)
						{
						fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(outgroupLabel,j,numNuc)]));
						}
					}
				else /* tip nodes */
					{
					for (j=1; j<=numNuc; j++)
						{
						fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
						}
					}
				fprintf (fpAlignment,"\n");
				}	
			fprintf (fpAlignment,"\n");
			/*fprintf(stderr,"\n sequence %d go with nodo of label %d \n", i+1, i+1);*/
			}
		else
			{
			for (i=0; i<numSequences; i++)
				{
				fprintf (fpAlignment,">seq%05d\n", i+1);

				for (j=1; j<=numNuc; j++)
					fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
				fprintf (fpAlignment,"\n");
				}	
			fprintf (fpAlignment,"\n");
			}
		}
	else /* migration */
		{
		if(thereisOutgroup == YES)
			{	

			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 2)
					{
					outgroupLabel = f->label;
					break;
					}
				}
			/*fprintf (stderr, "\n outgroupLabel = %d \n", outgroupLabel);*/


			for (i=0; i<numSequences+1; i++)
				{

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					if ((f->label == i) && (f->index < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}
				
				if (i == numSequences)
					fprintf(fpAlignment, ">outgrp_p0\n");
				else 
					fprintf (fpAlignment,">s%05d_p%d\n", i+1, dem);

				
				if (i == numSequences) /* outgroup */
					{
					for (j=1; j<=numNuc; j++)
						{
						fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(outgroupLabel,j,numNuc)]));
						}
					}
				else /* tip nodes */
					{
					for (j=1; j<=numNuc; j++)
						{
						fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
						}
					}

				fprintf (fpAlignment,"\n");
				}	
			fprintf (fpAlignment,"\n");
			}
		else
			{
			for (i=0; i<numSequences; i++)
				{
				/*for (m = 0; m < numNodex; m++)
					{
					f = nodex + m;
					if ((f->label == i) && (f->NetIndex < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}*/

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					if ((f->label == i) && (f->index < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}

				fprintf (fpAlignment,">s%05d_p%d\n", i+1, dem);
				
				for (j=1; j<=numNuc; j++)
					{
					fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
					}
				
				fprintf (fpAlignment,"\n");
				}	
			fprintf (fpAlignment,"\n");
			}
		}
}
		


/***************************** PrintAncestralSequences_FASTA *******************************/
/* Prints ancestral sequences to alignment FASTA file */
static void PrintAncestralSequences_FASTA (/*int replicate*/)   
{
	int		 i, j, a, m, dem, outgroupLabel, rootLabel;
	TreeNode	*f;
	dem = 0;
	outgroupLabel = rootLabel = 0;

	if (numRE == 0) /* There are NOT recombinations */
		{
		if (thereisOutgroup == YES)
			{
			/*fprintf(fpAlignment,"%d %d\n", 2*numSequences, numNuc);*/
			/*fprintf(fpAlignment,"Dataset_%d %d %d\n", replicate+1, 2*numSequences, numNuc);*/
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 2)
					outgroupLabel = f->label;
				if (f->class == 5)
					rootLabel = f->label;
				}
			}
		else
			{
			/*fprintf(fpAlignment,"%d %d\n", 2*numSequences-1, numNuc);*/
			/*fprintf(fpAlignment,"Dataset_%d %d %d\n",replicate+1, 2*numSequences-1, numNuc);*/
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 5)
					rootLabel = f->label;
				}
			}


		
		if (doMigration == NO)
			{
			if(thereisOutgroup == YES)
				{	
				for (i=0; i<2*numSequences; i++)
					{
					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,">seq%05d\n", i+1);
					else if (i == numSequences)
						fprintf(fpAlignment, ">outgroup\n");
					else if (i == 2*numSequences-1)
						fprintf (fpAlignment,">root\n");
					else
						fprintf (fpAlignment,">anc%05d\n", i+1);
		

					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
							}
						}
					else if (i == numSequences) /* is outgroup */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(outgroupLabel,j,numNuc)]));
							}
						}
					else if (i == 2*numSequences-1) /* is root */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
							}
						}
					else  /* is ancestral */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i-1,j,numNuc)]));
							}
						}

					fprintf (fpAlignment,"\n");
					}
				fprintf (fpAlignment,"\n");
				}
			else
				{
				for (i=0; i<2*numSequences-1; i++)
					{
					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,">seq%05d\n", i+1);
					else if (i == 2*numSequences-2)
						fprintf (fpAlignment,">root\n");
					else
						fprintf (fpAlignment,">anc%05d\n", i+1);


					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
							}
						}
					else if (i == 2*numSequences-2) /* is root */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
							}
						}
					else /* is ancestral */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
							}
						}
					fprintf (fpAlignment,"\n");
					}	
				fprintf (fpAlignment,"\n");
				}
			}
		else /* migration */
			{
			if(thereisOutgroup == YES)
				{	
				for (i=0; i<2*numSequences; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodes + m;
						if ((f->label == i) && (f->index <= numSequences*2-2))
							{
							dem = f->indexOldMigPop;
							break;
							}
						}*/

					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;

						if (i < numSequences) /* tip */
							{
							if ((f->label == i) && (f->index <= numSequences*2-2)) 
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == 2*numSequences-1) /* root */
							{
							if (f->class == 5)
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == numSequences) /* outgroup */
							{
							dem = 0;
							break;
							}
						else /* ancestral */
							{
							if (f->label == i-1)
								{
								dem = f->indexOldMigPop;
								if (f->class != 4)
									{
									fprintf (stderr, "\n Warning in PrintAncestralSequences_C. f->label = %d, f->class = %d \n", f->label, f->class);
									exit(-1);
									}
								break;
								}
							}
						}
					
					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,">s%05d_p%d\n", i+1, dem);
					else if (i == numSequences)
						fprintf(fpAlignment, ">outgrp_p0\n");
					else if (i == 2*numSequences-1)
						fprintf (fpAlignment,">root_p%d\n", dem);
					else
						fprintf (fpAlignment,">a%05d_p%d\n", i+1, dem);
		

					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
							}
						}
					else if (i == numSequences) /* is outgroup */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(outgroupLabel,j,numNuc)]));
							}
						}
					else if (i == 2*numSequences-1) /* is root */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
							}
						}
					else  /* is ancestral */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i-1,j,numNuc)]));
							}
						}
					fprintf (fpAlignment,"\n");
					/*fprintf(stderr," \n\n sequence %d with label(in times file) %d", i+1, i+1);*/
					}	
				fprintf (fpAlignment,"\n");
				}
			else
				{
				for (i=0; i<2*numSequences-1; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodes + m;
						if ((f->label == i) && (f->index <= numSequences*2-2))
							{
							dem = f->indexOldMigPop;
							break;
							}
						}*/

					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;

						if (i < numSequences) /* tip */
							{
							if ((f->label == i) && (f->index <= numSequences*2-2)) 
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == 2*numSequences-2) /* root */
							{
							if (f->class == 5)
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else /* ancestral */
							{
							if (f->label == i)
								{
								dem = f->indexOldMigPop;
								if (f->class != 4)
									{
									fprintf (stderr, "\n Warning in PrintAncestralSequences_C. f->label = %d, f->class = %d \n", f->label, f->class);
									exit(-1);
									}
								break;
								}
							}
						}
						
					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,">s%05d_p%d\n", i+1, dem);
					else if (i == 2*numSequences-2)
						fprintf (fpAlignment,">root_p%d\n", dem);
					else
						fprintf (fpAlignment,">a%05d_p%d\n", i+1, dem);

					
					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
							}
						}
					else if (i == 2*numSequences-2) /* is root */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
							}
						}
					else /* is ancestral */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
							}
						}

					fprintf (fpAlignment,"\n");
					}	
				fprintf (fpAlignment,"\n");
				}
			}
		}
	else /* There are recombinations */
		{
		/*if (thereisOutgroup == YES)
			fprintf(fpAlignment,"%d %d \n", numSequences+2, numNuc);*/
			/*fprintf(fpAlignment,"Dataset_%d %d %d \n", replicate+1, numSequences+1, numNuc);*/
		/*else*/
			/*fprintf(fpAlignment,"Dataset_%d %d %d \n", replicate+1, numSequences, numNuc);*/
			/*fprintf(fpAlignment,"%d %d \n", numSequences+1, numNuc);*/
		
		if (doMigration == NO)
			{
			if(thereisOutgroup == YES)
				{	

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 2)
						outgroupLabel = f->label;
					if (f->class == 5)
						rootLabel = f->label;
					}		


				for (i=0; i<2*numSequences; i++)
					{
					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,">seq%05d\n", i+1);
						a++;
						}
					if (i == numSequences)
						{
						fprintf(fpAlignment, ">outgroup\n");
						a++;
						}
					if (i == 2*numSequences-1)
						{
						fprintf (fpAlignment,">root\n");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == numSequences) /* is outgroup */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(outgroupLabel,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-1) /* is root */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}

						/*fprintf(stderr," \n\n sequence %d with label(in times file) %d", i+1, i+1);*/
						}
					}	
				fprintf (fpAlignment,"\n");
				}
			else
				{
				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 5)
						rootLabel = f->label;
					}


				for (i=0; i<2*numSequences-1; i++)
					{
					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,">seq%05d\n", i+1);
						a++;
						}
					if (i == 2*numSequences-2)
						{
						fprintf (fpAlignment,">root\n");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{

						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-2) /* is root */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}


						}
					}	
				fprintf (fpAlignment,"\n");
				}
			}
		else /* migration */
			{
			if(thereisOutgroup == YES)
				{	
				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 2)
						outgroupLabel = f->label;
					if (f->class == 5)
						rootLabel = f->label;
					}	





				for (i=0; i<2*numSequences; i++)
					{					
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodes + m;
						if (f->index < numSequences*/ /*|| i == 2*numSequences-1*//*)
							if ( f->label == i)
								{
								dem = f->indexOldMigPop;*/
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
							/*	break;
								}
						}*/
					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;
						if (f->index < numSequences /*|| i == 2*numSequences-1*/)
							if (f->label == i)
								{
								dem = f->indexOldMigPop;
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								break;
								}
						}

					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,">s%05d_p%d\n", i+1, dem);
						a++;
						}
					if (i == numSequences)
						{
						fprintf(fpAlignment, ">outgrp_p0\n");
						a++;
						}
					if (i == 2*numSequences-1)
						{
						fprintf (fpAlignment,">root\n");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{

						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == numSequences) /* is outgroup */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(outgroupLabel,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-1) /* is root */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}

						}
					}	
				fprintf (fpAlignment,"\n");
				}
			else
				{

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 5)
						rootLabel = f->label;
					}


				for (i=0; i<2*numSequences-1; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodes + m;
						if (f->index < numSequences*/ /*|| i == 2*numSequences-1*//*)
							if ( f->label == i)
								{
								dem = f->indexOldMigPop;
								break;
								}
						}*/
					
					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;
						if (f->index < numSequences /*|| i == 2*numSequences-1*/)
							if ( f->label == i)
								{
								dem = f->indexOldMigPop;
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								break;
								}
						}


					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,">s%05d_p%d\n", i+1, dem);
						a++;
						}
					if (i == 2*numSequences-2)
						{
						fprintf (fpAlignment,">root\n");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-2) /* is root */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}

						}
					}	
				fprintf (fpAlignment,"\n");
				}

			}
		}
}
		



		

/***************************** PrintSequences for codon Model *******************************/
/* Prints sequences to alignment file in fasta sequential format */

static void PrintSequences_C_FASTA (/*int replicate*/)
{
	int		 i, j, k, m, dem, outgroupLabel;
	char codon[3];
	TreeNode	*f;
	dem = 0;
	outgroupLabel = 0;
	
	/*if (thereisOutgroup == YES)
		fprintf(fpAlignment,"%d %d \n", numSequences+1, numNuc);*/
		/*fprintf(fpAlignment,"Dataset_%d %d %d \n", replicate+1, numSequences+1, numNuc);*/
	/*else
		fprintf(fpAlignment,"%d %d \n", numSequences, numNuc);*/
		/*fprintf(fpAlignment,"Dataset_%d %d %d \n", replicate+1, numSequences, numNuc);*/

	if (doMigration == NO)
		{
		if(thereisOutgroup == YES)
			{
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 2)
					{
					outgroupLabel = f->label;
					break;
					}
				}
	
			for (i=0; i<numSequences+1; i++)
				{
				if (i == numSequences)
					fprintf(fpAlignment, ">outgroup\n");
				else 
					fprintf (fpAlignment,">seq%05d\n", i+1);

				if (i == numSequences) /* outgroup */
					{
					for (j=1; j<=numSites; j++)
						{
						if (matrixC[pos(outgroupLabel,j,numSites)] > 60 || matrixC[pos(outgroupLabel,j,numSites)] > 60)
							{
							fprintf (stderr, "\n stop codon5 \n");
							exit(-1);
							}
						number_to_codon(matrixC[pos(outgroupLabel,j,numSites)], codon);
						for (k = 0; k < 3; k++)
							fprintf (fpAlignment, "%c", codon[k]);
						}
					}
				else /* tip nodes */
					{
					for (j=1; j<=numSites; j++)
						{
						if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
							{
							fprintf (stderr, "\n stop codon5 \n");
							exit(-1);
							}
						number_to_codon(matrixC[pos(i,j,numSites)], codon);
						for (k = 0; k < 3; k++)
							fprintf (fpAlignment, "%c", codon[k]);
						}
					}
				fprintf (fpAlignment,"\n");
				/*fprintf(stderr,"\n sequence %d go with nodo of label %d \n", i+1, i+1);*/
				}	
			fprintf (fpAlignment,"\n");
			}
		else
			{
			for (i=0; i<numSequences; i++)
				{
				fprintf (fpAlignment,">seq%05d\n", i+1);

				for (j=1; j<=numSites; j++)
					{
					if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
						{
						fprintf (stderr, "\n stop codon6 \n");
						exit(-1);
						}
					number_to_codon(matrixC[pos(i,j,numSites)], codon);
					for (k = 0; k < 3; k++)
						fprintf (fpAlignment, "%c", codon[k]);
					}
				fprintf (fpAlignment,"\n");
				}	
			fprintf (fpAlignment,"\n");
			}
		}
	else /* migration */
		{
		if(thereisOutgroup == YES)
			{
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 2)
					{
					outgroupLabel = f->label;
					break;
					}
				}
			/*fprintf (stderr, "\n outgroupLabel = %d \n", outgroupLabel);*/
			

			for (i=0; i<numSequences+1; i++)
				{
				/*for (m = 0; m < numNodex; m++)
					{
					f = nodex + m;
					if ((f->label == i) && (f->NetIndex < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}*/
				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					if ((f->label == i) && (f->index < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}		
		
				if (i == numSequences)
					fprintf(fpAlignment, ">outgrp_p0\n");
				else 
					fprintf (fpAlignment,">s%05d_p%d\n", i+1, dem);
				
				if (i == numSequences) /* outgroup */
					{
					for (j=1; j<=numSites; j++)
						{
						if (matrixC[pos(outgroupLabel,j,numSites)] > 60 || matrixC[pos(outgroupLabel,j,numSites)] > 60)
							{
							fprintf (stderr, "\n stop codon7 \n");
							exit(-1);
							}
						number_to_codon(matrixC[pos(outgroupLabel,j,numSites)], codon);
						for (k = 0; k < 3; k++)
							fprintf (fpAlignment, "%c", codon[k]);
						}
					}
				else /* tip nodes */
					{
					for (j=1; j<=numSites; j++)
						{
						if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
							{
							fprintf (stderr, "\n stop codon7 \n");
							exit(-1);
							}
						number_to_codon(matrixC[pos(i,j,numSites)], codon);
						for (k = 0; k < 3; k++)
							fprintf (fpAlignment, "%c", codon[k]);
						}
					}

				fprintf (fpAlignment,"\n");
				/*fprintf(stderr,"\n sequence %d go with node label(times file label) %d \n", i+1, i+1);*/
				}	
			fprintf (fpAlignment,"\n");
			}
		else
			{
			for (i=0; i<numSequences; i++)
				{
				/*for (m = 0; m < numNodex; m++)
					{
					f = nodex + m;
					if ((f->label == i) && (f->NetIndex < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}*/
				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					if ((f->label == i) && (f->index < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}
					
				fprintf (fpAlignment,">s%05d_p%d\n", i+1, dem);

				for (j=1; j<=numSites; j++)
					{
					if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
						{
						fprintf (stderr, "\n stop codon8 \n");
						exit(-1);
						}
					number_to_codon(matrixC[pos(i,j,numSites)], codon);
					for (k = 0; k < 3; k++)
						fprintf (fpAlignment, "%c", codon[k]);
					}
				fprintf (fpAlignment,"\n");
				}	
			fprintf (fpAlignment,"\n");
			}
		}
}
		


/***************************** PrintAncestralSequences for Codon Model *******************************/
/* Prints ancestral sequences to alignment file for codon Model */
static void PrintAncestralSequences_C_FASTA (/*int replicate*/)
{
	int		 i, j, k, a, m, dem, outgroupLabel, rootLabel;
	char codon[3];
	TreeNode	*f;
	dem = 0;
	outgroupLabel = rootLabel = 0;	

	if (numRE == 0) /* There are NOT recombinations */
		{
		if (thereisOutgroup == YES)
			{
			/*fprintf(fpAlignment,"%d %d\n", 2*numSequences, numNuc);*/
			/*fprintf(fpAlignment,"Dataset_%d %d %d\n", replicate+1, 2*numSequences, numNuc);*/
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 2)
					outgroupLabel = f->label;
				if (f->class == 5)
					rootLabel = f->label;
				}
			}
		else
			{
			/*fprintf(fpAlignment,"%d %d\n", 2*numSequences-1, numNuc);*/
			/*fprintf(fpAlignment,"Dataset_%d %d %d\n",replicate+1, 2*numSequences-1, numNuc);*/
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 5)
					rootLabel = f->label;
				}
			}
		

		if (doMigration == NO)
			{
			if (thereisOutgroup == YES)
				{	
				for (i=0; i<2*numSequences; i++)
					{
					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,">seq%05d\n", i+1);
					else if (i == numSequences)
						fprintf(fpAlignment, ">outgroup\n");
					else if (i == 2*numSequences-1)
						fprintf (fpAlignment,">root\n");
					else
						fprintf (fpAlignment,">anc%05d\n", i+1);


					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else if (i == numSequences) /* is outgroup */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(outgroupLabel,j,numSites)] > 60 || matrixC[pos(outgroupLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(outgroupLabel,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else if (i == 2*numSequences-1) /* is root */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else  /* is ancestral */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i-1,j,numSites)] > 60 || matrixC[pos(i-1,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i-1,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}

					fprintf (fpAlignment,"\n");
					}	
				fprintf (fpAlignment,"\n");
				}
			else
				{
				for (i=0; i<2*numSequences-1; i++)
					{
					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,">seq%05d\n", i+1);
					else if (i == 2*numSequences-2)
						fprintf (fpAlignment,">root\n");
					else
						fprintf (fpAlignment,">anc%05d\n", i+1);

					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon11 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i,j,numSites)], codon);
							for	(k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else if (i == 2*numSequences-2) /* is root */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon11 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
							for	(k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else /* is ancestral */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon11 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i,j,numSites)], codon);
							for	(k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}

					fprintf (fpAlignment,"\n");
					}	
				fprintf (fpAlignment,"\n");
				}
			}
		else /* migration */
			{
			if (thereisOutgroup == YES)
				{	
				for (i=0; i<2*numSequences; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodex + m;
						if ((f->label == i) && (f->NetIndex <= numSequences*2-2))
							{
							dem = f->indexOldMigPop;
							break;
							}
						}*/
					
					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;

						if (i < numSequences) /* tip */
							{
							if ((f->label == i) && (f->index <= numSequences*2-2)) 
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == 2*numSequences-1) /* root */
							{
							if (f->class == 5)
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == numSequences) /* outgroup */
							{
							dem = 0;
							break;
							}
						else /* ancestral */
							{
							if (f->label == i-1)
								{
								dem = f->indexOldMigPop;
								if (f->class != 4)
									{
									fprintf (stderr, "\n Warning in PrintAncestralSequences_C. f->label = %d, f->class = %d \n", f->label, f->class);
									exit(-1);
									}
								break;
								}
							}
						}



					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,">s%05d_p%d\n", i+1, dem);
					else if (i == numSequences)
						fprintf(fpAlignment, ">outgrp_p0\n");
					else if (i == 2*numSequences-1)
						fprintf (fpAlignment,">root_p%d\n",dem);
					else
						fprintf (fpAlignment,">a%05d_p%d\n", i+1, dem);
		


					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else if (i == numSequences) /* is outgroup */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(outgroupLabel,j,numSites)] > 60 || matrixC[pos(outgroupLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(outgroupLabel,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else if (i == 2*numSequences-1) /* is root */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else  /* is ancestral */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i-1,j,numSites)] > 60 || matrixC[pos(i-1,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i-1,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}


					fprintf (fpAlignment,"\n");
					}	
				fprintf (fpAlignment,"\n");
				}
			else
				{
				for (i=0; i<2*numSequences-1; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodex + m;
						if ((f->label == i) && (f->NetIndex <= numSequences*2-2))
							{
							dem = f->indexOldMigPop;
							break;
							}
						}*/
					
					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;

						if (i < numSequences) /* tip */
							{
							if ((f->label == i) && (f->index <= numSequences*2-2)) 
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == 2*numSequences-2) /* root */
							{
							if (f->class == 5)
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else /* ancestral */
							{
							if (f->label == i)
								{
								dem = f->indexOldMigPop;
								if (f->class != 4)
									{
									fprintf (stderr, "\n Warning in PrintAncestralSequences_C. f->label = %d, f->class = %d \n", f->label, f->class);
									exit(-1);
									}
								break;
								}
							}
						}




					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,">s%05d_p%d\n", i+1, dem);
					else if (i == 2*numSequences-2)
						fprintf (fpAlignment,">root_p%d\n", dem);
					else
						fprintf (fpAlignment,">a%05d_p%d\n", i+1, dem);


					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon11 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i,j,numSites)], codon);
							for	(k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else if (i == 2*numSequences-2) /* is root */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon11 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
							for	(k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else /* is ancestral */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon11 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i,j,numSites)], codon);
							for	(k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}

					fprintf (fpAlignment,"\n");

					}	
				}
			fprintf (fpAlignment,"\n");
			}
		}
	else /* There are recombinations */
		{
		/*if (thereisOutgroup == YES)
			fprintf(fpAlignment,"%d %d\n", numSequences+2, numNuc);
		else
			fprintf(fpAlignment,"%d %d\n", numSequences+1, numNuc);*/
		
		if (doMigration == NO)
			{
			if (thereisOutgroup == YES)
				{
				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 2)
						outgroupLabel = f->label;
					if (f->class == 5)
						rootLabel = f->label;
					}				
	
				for (i=0; i<2*numSequences; i++)
					{
					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,">seq%05d\n", i+1);
						a++;
						}
					if (i == numSequences)
						{
						fprintf(fpAlignment, ">outgroup\n");
						a++;
						}
					if (i == 2*numSequences-1)
						{
						fprintf (fpAlignment,">root\n");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon14 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(i,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == numSequences) /* is outgroup */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(outgroupLabel,j,numSites)] > 60 || matrixC[pos(outgroupLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon14 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(outgroupLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-1) /* is root */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon14 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}

						}
					}	
				fprintf (fpAlignment,"\n");
				}
			else
				{

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 5)
						rootLabel = f->label;
					}

				for (i=0; i<2*numSequences-1; i++)
					{
					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,">seq%05d\n", i+1);
						a++;
						}
					if (i == 2*numSequences-2)
						{
						fprintf (fpAlignment,">root\n");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon15 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(i,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-2) /* is root */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon15 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}

						}
					}	
				fprintf (fpAlignment,"\n");
				}
			}
		else /* migration */
			{
			if (thereisOutgroup == YES)
				{	

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 2)
						outgroupLabel = f->label;
					if (f->class == 5)
						rootLabel = f->label;
					}	


				for (i=0; i<2*numSequences; i++)
					{

					/*for (m = 0; m < numNodex; m++)
						{
						f = nodex + m;
						if (f->NetIndex < numSequences*/ /*|| i == 2*numSequences-1*//*
							if ( f->label == i)
								{
								dem = f->indexOldMigPop;*/
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								/*break;
								}
						}*/

					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;
						if (f->index < numSequences /*|| i == 2*numSequences-1*/)
							if (f->label == i)
								{
								dem = f->indexOldMigPop;
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								break;
								}
						}



					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,">s%05d_p%d\n", i+1, dem);
						a++;
						}
					if (i == numSequences)
						{
						fprintf(fpAlignment, ">outgrp_p0\n");
						a++;
						}
					if (i == 2*numSequences-1)
						{
						fprintf (fpAlignment,">root\n");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{

						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon16 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(i,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == numSequences) /* is outgroup */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(outgroupLabel,j,numSites)] > 60 || matrixC[pos(outgroupLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon16 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(outgroupLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-1) /* is root */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon16 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}

						}
					}	
				fprintf (fpAlignment,"\n");
				}
			else
				{

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 5)
						rootLabel = f->label;
					}

				for (i=0; i<2*numSequences-1; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodex + m;
						if (f->NetIndex < numSequences*/ /*|| i == 2*numSequences-1*//*)
							if (f->label == i)
								{
								dem = f->indexOldMigPop;*/
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								/*break;
								}
						}*/
					
					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;
						if (f->index < numSequences /*|| i == 2*numSequences-1*/)
							if ( f->label == i)
								{
								dem = f->indexOldMigPop;
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								break;
								}
						}



					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,">s%05d_p%d\n", i+1, dem);
						a++;
						}
					if (i == 2*numSequences-2)
						{
						fprintf (fpAlignment,">root\n");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon17 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(i,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-2) /* is root */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon17 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}


						}
					}	
				fprintf (fpAlignment,"\n");
				}
			}
		}
}






/**** NEXUS ****/
/***************************** PrintSequences_NEXUS *******************************/
/* Prints sequences to nexus alignment file */

static void PrintSequences_NEXUS (/*int replicate*/) 
{
	int		 i, j, m, dem, outgroupLabel;
	TreeNode	*f;
	
	dem = 0;
	outgroupLabel = 0;

	PrintNEXUS_initial();

	if (thereisOutgroup == YES)
		{
		fprintf(fpAlignment,"Dimensions ntax=%d nchar=%d;\n", numSequences+1, numNuc);
		/*fprintf(fpAlignment,"Dataset_%d %d %d \n", replicate+1, numSequences+1, numNuc);*/
		}
	else
		{
		/*fprintf(fpAlignment,"Dataset_%d %d %d \n", replicate+1, numSequences, numNuc);*/
		fprintf(fpAlignment,"Dimensions ntax=%d nchar=%d;\n", numSequences, numNuc);
		}

	fprintf(fpAlignment,"	Format datatype=nucleotide gap=- missing=? matchchar=.;\n");
	fprintf(fpAlignment,"	Matrix\n");

	
	if (doMigration == NO)
		{
		if(thereisOutgroup == YES)
			{
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 2)
					{
					outgroupLabel = f->label;
					break;
					}
				}
	
			for (i=0; i<numSequences+1; i++)
				{
				if (i == numSequences)
					fprintf(fpAlignment, "outgroup     ");
				else 
					fprintf (fpAlignment,"seq%05d     ", i+1);
				
				if (i == numSequences) /* outgroup */
					{
					for (j=1; j<=numNuc; j++)
						{
						fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(outgroupLabel,j,numNuc)]));
						}
					}
				else /* tip nodes */
					{
					for (j=1; j<=numNuc; j++)
						{
						fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
						}
					}
				fprintf (fpAlignment,"\n");
				}	
			/*fprintf (fpAlignment,"\n");*/
			/*fprintf(stderr,"\n sequence %d go with nodo of label %d \n", i+1, i+1);*/
			}
		else
			{
			for (i=0; i<numSequences; i++)
				{
				fprintf (fpAlignment,"seq%05d     ", i+1);

				for (j=1; j<=numNuc; j++)
					fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
				fprintf (fpAlignment,"\n");
				}	
			/*fprintf (fpAlignment,"\n");*/
			}
		}
	else /* migration */
		{
		if(thereisOutgroup == YES)
			{	

			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 2)
					{
					outgroupLabel = f->label;
					break;
					}
				}
			/*fprintf (stderr, "\n outgroupLabel = %d \n", outgroupLabel);*/


			for (i=0; i<numSequences+1; i++)
				{

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					if ((f->label == i) && (f->index < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}
				
				if (i == numSequences)
					fprintf(fpAlignment, "outgrp_p0    ");
				else 
					fprintf (fpAlignment,"s%05d_p%d    ", i+1, dem);

				
				if (i == numSequences) /* outgroup */
					{
					for (j=1; j<=numNuc; j++)
						{
						fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(outgroupLabel,j,numNuc)]));
						}
					}
				else /* tip nodes */
					{
					for (j=1; j<=numNuc; j++)
						{
						fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
						}
					}

				fprintf (fpAlignment,"\n");
				}	
			/*fprintf (fpAlignment,"\n");*/
			}
		else
			{
			for (i=0; i<numSequences; i++)
				{
				/*for (m = 0; m < numNodex; m++)
					{
					f = nodex + m;
					if ((f->label == i) && (f->NetIndex < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}*/

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					if ((f->label == i) && (f->index < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}

				fprintf (fpAlignment,"s%05d_p%d    ", i+1, dem);
				
				for (j=1; j<=numNuc; j++)
					{
					fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
					}
				
				fprintf (fpAlignment,"\n");
				}	
			/*fprintf (fpAlignment,"\n");*/
			}
		}
	PrintNEXUS_end();
}
		


/***************************** PrintAncestralSequences_NEXUS *******************************/
/* Prints ancestral sequences to nexus alignment file */

static void PrintAncestralSequences_NEXUS (/*int replicate*/)   
{
	int		 i, j, a, m, dem, outgroupLabel, rootLabel;
	TreeNode	*f;
	dem = 0;
	outgroupLabel = rootLabel = 0;

	PrintNEXUS_initial();

	if (numRE == 0) /* There are NOT recombinations */
		{
		if (thereisOutgroup == YES)
			{
			fprintf(fpAlignment,"Dimensions ntax=%d nchar=%d;\n", 2*numSequences, numNuc);
			/*fprintf(fpAlignment,"Dataset_%d %d %d\n", replicate+1, 2*numSequences, numNuc);*/
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 2)
					outgroupLabel = f->label;
				if (f->class == 5)
					rootLabel = f->label;
				}
			}
		else
			{
			fprintf(fpAlignment,"Dimensions ntax=%d nchar=%d;\n", 2*numSequences-1, numNuc);
			/*fprintf(fpAlignment,"Dataset_%d %d %d\n",replicate+1, 2*numSequences-1, numNuc);*/
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 5)
					rootLabel = f->label;
				}
			}
		
		fprintf(fpAlignment,"	Format datatype=nucleotide gap=- missing=? matchchar=.;\n");
		fprintf(fpAlignment,"	Matrix\n");

		
		if (doMigration == NO)
			{
			if(thereisOutgroup == YES)
				{	
				for (i=0; i<2*numSequences; i++)
					{
					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,"seq%05d     ", i+1);
					else if (i == numSequences)
						fprintf(fpAlignment, "outgroup     ");
					else if (i == 2*numSequences-1)
						fprintf (fpAlignment,"root         ");
					else
						fprintf (fpAlignment,"anc%05d     ", i+1);
		

					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
							}
						}
					else if (i == numSequences) /* is outgroup */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(outgroupLabel,j,numNuc)]));
							}
						}
					else if (i == 2*numSequences-1) /* is root */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
							}
						}
					else  /* is ancestral */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i-1,j,numNuc)]));
							}
						}

					fprintf (fpAlignment,"\n");
					}
				/*fprintf (fpAlignment,"\n");*/
				}
			else
				{
				for (i=0; i<2*numSequences-1; i++)
					{
					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,"seq%05d     ", i+1);
					else if (i == 2*numSequences-2)
						fprintf (fpAlignment,"root         ");
					else
						fprintf (fpAlignment,"anc%05d     ", i+1);


					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
							}
						}
					else if (i == 2*numSequences-2) /* is root */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
							}
						}
					else /* is ancestral */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
							}
						}
					fprintf (fpAlignment,"\n");
					}	
				/*fprintf (fpAlignment,"\n");*/
				}
			}
		else /* migration */
			{
			if(thereisOutgroup == YES)
				{	
				for (i=0; i<2*numSequences; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodes + m;
						if ((f->label == i) && (f->index <= numSequences*2-2))
							{
							dem = f->indexOldMigPop;
							break;
							}
						}*/

					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;

						if (i < numSequences) /* tip */
							{
							if ((f->label == i) && (f->index <= numSequences*2-2)) 
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == 2*numSequences-1) /* root */
							{
							if (f->class == 5)
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == numSequences) /* outgroup */
							{
							dem = 0;
							break;
							}
						else /* ancestral */
							{
							if (f->label == i-1)
								{
								dem = f->indexOldMigPop;
								if (f->class != 4)
									{
									fprintf (stderr, "\n Warning in PrintAncestralSequences_C. f->label = %d, f->class = %d \n", f->label, f->class);
									exit(-1);
									}
								break;
								}
							}
						}
					
					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,"s%05d_p%d    ", i+1, dem);
					else if (i == numSequences)
						fprintf(fpAlignment, "outgrp_p0    ");
					else if (i == 2*numSequences-1)
						fprintf (fpAlignment,"root_p%d      ", dem);
					else
						fprintf (fpAlignment,"a%05d_p%d    ", i+1, dem);
		

					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
							}
						}
					else if (i == numSequences) /* is outgroup */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(outgroupLabel,j,numNuc)]));
							}
						}
					else if (i == 2*numSequences-1) /* is root */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
							}
						}
					else  /* is ancestral */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i-1,j,numNuc)]));
							}
						}
					fprintf (fpAlignment,"\n");
					/*fprintf(stderr," \n\n sequence %d with label(in times file) %d", i+1, i+1);*/
					}	
				/*fprintf (fpAlignment,"\n");*/
				}
			else
				{
				for (i=0; i<2*numSequences-1; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodes + m;
						if ((f->label == i) && (f->index <= numSequences*2-2))
							{
							dem = f->indexOldMigPop;
							break;
							}
						}*/

					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;

						if (i < numSequences) /* tip */
							{
							if ((f->label == i) && (f->index <= numSequences*2-2)) 
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == 2*numSequences-2) /* root */
							{
							if (f->class == 5)
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else /* ancestral */
							{
							if (f->label == i)
								{
								dem = f->indexOldMigPop;
								if (f->class != 4)
									{
									fprintf (stderr, "\n Warning in PrintAncestralSequences_C. f->label = %d, f->class = %d \n", f->label, f->class);
									exit(-1);
									}
								break;
								}
							}
						}
						
					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,"s%05d_p%d    ", i+1, dem);
					else if (i == 2*numSequences-2)
						fprintf (fpAlignment,"root_p%d      ", dem);
					else
						fprintf (fpAlignment,"a%05d_p%d    ", i+1, dem);

					
					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
							}
						}
					else if (i == 2*numSequences-2) /* is root */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
							}
						}
					else /* is ancestral */
						{
						for (j=1; j<=numNuc; j++)
							{
							fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
							}
						}

					fprintf (fpAlignment,"\n");
					}	
				/*fprintf (fpAlignment,"\n");*/
				}
			}
		}
	else /* There are recombinations */
		{
		if (thereisOutgroup == YES)
			{
			fprintf(fpAlignment,"Dimensions ntax=%d nchar=%d;\n", numSequences+2, numNuc);
			/*fprintf(fpAlignment,"Dataset_%d %d %d \n", replicate+1, numSequences+1, numNuc);*/
			}
		else
			{
			/*fprintf(fpAlignment,"Dataset_%d %d %d \n", replicate+1, numSequences, numNuc);*/
			fprintf(fpAlignment,"Dimensions ntax=%d nchar=%d;\n", numSequences+1, numNuc);
			}
		
		if (doMigration == NO)
			{
			if(thereisOutgroup == YES)
				{	

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 2)
						outgroupLabel = f->label;
					if (f->class == 5)
						rootLabel = f->label;
					}		


				for (i=0; i<2*numSequences; i++)
					{
					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,"seq%05d     ", i+1);
						a++;
						}
					if (i == numSequences)
						{
						fprintf(fpAlignment, "outgroup     ");
						a++;
						}
					if (i == 2*numSequences-1)
						{
						fprintf (fpAlignment,"root         ");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == numSequences) /* is outgroup */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(outgroupLabel,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-1) /* is root */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}

						/*fprintf(stderr," \n\n sequence %d with label(in times file) %d", i+1, i+1);*/
						}
					}	
				/*fprintf (fpAlignment,"\n");*/
				}
			else
				{
				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 5)
						rootLabel = f->label;
					}


				for (i=0; i<2*numSequences-1; i++)
					{
					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,"seq%05d     ", i+1);
						a++;
						}
					if (i == 2*numSequences-2)
						{
						fprintf (fpAlignment,"root         ");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{

						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-2) /* is root */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}

						}
					}	
				/*fprintf (fpAlignment,"\n");*/
				}
			}
		else /* migration */
			{
			if(thereisOutgroup == YES)
				{	
				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 2)
						outgroupLabel = f->label;
					if (f->class == 5)
						rootLabel = f->label;
					}	



				for (i=0; i<2*numSequences; i++)
					{					
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodes + m;
						if (f->index < numSequences*/ /*|| i == 2*numSequences-1*//*)
							if ( f->label == i)
								{
								dem = f->indexOldMigPop;*/
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
							/*	break;
								}
						}*/
					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;
						if (f->index < numSequences /*|| i == 2*numSequences-1*/)
							if (f->label == i)
								{
								dem = f->indexOldMigPop;
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								break;
								}
						}

					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,"s%05d_p%d    ", i+1, dem);
						a++;
						}
					if (i == numSequences)
						{
						fprintf(fpAlignment, "outgrp_p0    ");
						a++;
						}
					if (i == 2*numSequences-1)
						{
						fprintf (fpAlignment,"root         ");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == numSequences) /* is outgroup */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(outgroupLabel,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-1) /* is root */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}

						}
					}	
				/*fprintf (fpAlignment,"\n");*/
				}
			else
				{

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 5)
						rootLabel = f->label;
					}


				for (i=0; i<2*numSequences-1; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodes + m;
						if (f->index < numSequences*/ /*|| i == 2*numSequences-1*//*)
							if ( f->label == i)
								{
								dem = f->indexOldMigPop;
								break;
								}
						}*/
					
					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;
						if (f->index < numSequences /*|| i == 2*numSequences-1*/)
							if ( f->label == i)
								{
								dem = f->indexOldMigPop;
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								break;
								}
						}


					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,"s%05d_p%d    ", i+1, dem);
						a++;
						}
					if (i == 2*numSequences-2)
						{
						fprintf (fpAlignment,"root         ");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(i,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-2) /* is root */
							{
							for (j=1; j<=numNuc; j++)
								{
								fprintf (fpAlignment, "%c", WhichNuc(matrix[pos(rootLabel,j,numNuc)]));
								}
							fprintf (fpAlignment,"\n");
							}

						}
					}	
				/*fprintf (fpAlignment,"\n");*/
				}

			}
		}
	PrintNEXUS_end();
}
		















/***************************** PrintSequences for codon Model *******************************/
/* Prints sequences to alignment file in nexus sequential format */

static void PrintSequences_C_NEXUS (/*int replicate*/)
{
	int		 i, j, k, m, dem, outgroupLabel;
	char codon[3];
	TreeNode	*f;
	dem = 0;
	outgroupLabel = 0;
	
	PrintNEXUS_initial();

	if (thereisOutgroup == YES)
		fprintf(fpAlignment,"Dimensions ntax=%d nchar=%d;\n", numSequences+1, numNuc);
		/*fprintf(fpAlignment,"Dataset_%d %d %d \n", replicate+1, numSequences+1, numNuc);*/
	else
		fprintf(fpAlignment,"Dimensions ntax=%d nchar=%d;\n", numSequences, numNuc);
		/*fprintf(fpAlignment,"Dataset_%d %d %d \n", replicate+1, numSequences, numNuc);*/
	
	fprintf(fpAlignment,"	Format datatype=nucleotide gap=- missing=? matchchar=.;\n");
	fprintf(fpAlignment,"	Matrix\n");

	if (doMigration == NO)
		{
		if(thereisOutgroup == YES)
			{
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 2)
					{
					outgroupLabel = f->label;
					break;
					}
				}
	
			for (i=0; i<numSequences+1; i++)
				{
				if (i == numSequences)
					fprintf(fpAlignment, "outgroup     ");
				else 
					fprintf (fpAlignment,"seq%05d     ", i+1);

				if (i == numSequences) /* outgroup */
					{
					for (j=1; j<=numSites; j++)
						{
						if (matrixC[pos(outgroupLabel,j,numSites)] > 60 || matrixC[pos(outgroupLabel,j,numSites)] > 60)
							{
							fprintf (stderr, "\n stop codon5 \n");
							exit(-1);
							}
						number_to_codon(matrixC[pos(outgroupLabel,j,numSites)], codon);
						for (k = 0; k < 3; k++)
							fprintf (fpAlignment, "%c", codon[k]);
						}
					}
				else /* tip nodes */
					{
					for (j=1; j<=numSites; j++)
						{
						if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
							{
							fprintf (stderr, "\n stop codon5 \n");
							exit(-1);
							}
						number_to_codon(matrixC[pos(i,j,numSites)], codon);
						for (k = 0; k < 3; k++)
							fprintf (fpAlignment, "%c", codon[k]);
						}
					}
				fprintf (fpAlignment,"\n");
				/*fprintf(stderr,"\n sequence %d go with nodo of label %d \n", i+1, i+1);*/
				}	
			/*fprintf (fpAlignment,"\n");*/
			}
		else
			{
			for (i=0; i<numSequences; i++)
				{
				fprintf (fpAlignment,"seq%05d     ", i+1);

				for (j=1; j<=numSites; j++)
					{
					if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
						{
						fprintf (stderr, "\n stop codon6 \n");
						exit(-1);
						}
					number_to_codon(matrixC[pos(i,j,numSites)], codon);
					for (k = 0; k < 3; k++)
						fprintf (fpAlignment, "%c", codon[k]);
					}
				fprintf (fpAlignment,"\n");
				}	
			/*fprintf (fpAlignment,"\n");*/
			}
		}
	else /* migration */
		{
		if(thereisOutgroup == YES)
			{
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 2)
					{
					outgroupLabel = f->label;
					break;
					}
				}
			/*fprintf (stderr, "\n outgroupLabel = %d \n", outgroupLabel);*/
			

			for (i=0; i<numSequences+1; i++)
				{
				/*for (m = 0; m < numNodex; m++)
					{
					f = nodex + m;
					if ((f->label == i) && (f->NetIndex < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}*/
				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					if ((f->label == i) && (f->index < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}		
		
				if (i == numSequences)
					fprintf(fpAlignment, "outgrp_p0    ");
				else 
					fprintf (fpAlignment,"s%05d_p%d    ", i+1, dem);
				
				if (i == numSequences) /* outgroup */
					{
					for (j=1; j<=numSites; j++)
						{
						if (matrixC[pos(outgroupLabel,j,numSites)] > 60 || matrixC[pos(outgroupLabel,j,numSites)] > 60)
							{
							fprintf (stderr, "\n stop codon7 \n");
							exit(-1);
							}
						number_to_codon(matrixC[pos(outgroupLabel,j,numSites)], codon);
						for (k = 0; k < 3; k++)
							fprintf (fpAlignment, "%c", codon[k]);
						}
					}
				else /* tip nodes */
					{
					for (j=1; j<=numSites; j++)
						{
						if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
							{
							fprintf (stderr, "\n stop codon7 \n");
							exit(-1);
							}
						number_to_codon(matrixC[pos(i,j,numSites)], codon);
						for (k = 0; k < 3; k++)
							fprintf (fpAlignment, "%c", codon[k]);
						}
					}

				fprintf (fpAlignment,"\n");
				/*fprintf(stderr,"\n sequence %d go with node label(times file label) %d \n", i+1, i+1);*/
				}	
			/*fprintf (fpAlignment,"\n");*/
			}
		else
			{
			for (i=0; i<numSequences; i++)
				{
				/*for (m = 0; m < numNodex; m++)
					{
					f = nodex + m;
					if ((f->label == i) && (f->NetIndex < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}*/
				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					if ((f->label == i) && (f->index < numSequences))
						{
						dem = f->indexOldMigPop;
						break;
						}
					}
					
				fprintf (fpAlignment,"s%05d_p%d    ", i+1, dem);

				for (j=1; j<=numSites; j++)
					{
					if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
						{
						fprintf (stderr, "\n stop codon8 \n");
						exit(-1);
						}
					number_to_codon(matrixC[pos(i,j,numSites)], codon);
					for (k = 0; k < 3; k++)
						fprintf (fpAlignment, "%c", codon[k]);
					}
				fprintf (fpAlignment,"\n");
				}	
			/*fprintf (fpAlignment,"\n");*/
			}
		}
	PrintNEXUS_end();
}
		


/***************************** PrintAncestralSequences for Codon Model *******************************/
/* Prints ancestral sequences to nexus alignment file for codon Model */

static void PrintAncestralSequences_C_NEXUS (/*int replicate*/)
{
	int		 i, j, k, a, m, dem, outgroupLabel, rootLabel;
	char codon[3];
	TreeNode	*f;
	dem = 0;
	outgroupLabel = rootLabel = 0;	

	PrintNEXUS_initial();

	if (numRE == 0) /* There are NOT recombinations */
		{
		if (thereisOutgroup == YES)
			{
			fprintf(fpAlignment,"Dimensions ntax=%d nchar=%d;\n", 2*numSequences, numNuc);
			/*fprintf(fpAlignment,"Dataset_%d %d %d\n", replicate+1, 2*numSequences, numNuc);*/
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 2)
					outgroupLabel = f->label;
				if (f->class == 5)
					rootLabel = f->label;
				}
			}
		else
			{
			fprintf(fpAlignment,"Dimensions ntax=%d nchar=%d;\n", 2*numSequences-1, numNuc);
			/*fprintf(fpAlignment,"Dataset_%d %d %d\n",replicate+1, 2*numSequences-1, numNuc);*/
			for (m = 0; m < nextAvailable; m++)
				{
				f = nodes + m;
				/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
				if (f->class == 5)
					rootLabel = f->label;
				}
			}

		fprintf(fpAlignment,"	Format datatype=nucleotide gap=- missing=? matchchar=.;\n");
		fprintf(fpAlignment,"	Matrix\n");

		if (doMigration == NO)
			{
			if (thereisOutgroup == YES)
				{	
				for (i=0; i<2*numSequences; i++)
					{
					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,"seq%05d     ", i+1);
					else if (i == numSequences)
						fprintf(fpAlignment, "outgroup     ");
					else if (i == 2*numSequences-1)
						fprintf (fpAlignment,"root         ");
					else
						fprintf (fpAlignment,"anc%05d     ", i+1);


					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else if (i == numSequences) /* is outgroup */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(outgroupLabel,j,numSites)] > 60 || matrixC[pos(outgroupLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(outgroupLabel,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else if (i == 2*numSequences-1) /* is root */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else  /* is ancestral */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i-1,j,numSites)] > 60 || matrixC[pos(i-1,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i-1,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}

					fprintf (fpAlignment,"\n");
					}	
				/*fprintf (fpAlignment,"\n");*/
				}
			else
				{
				for (i=0; i<2*numSequences-1; i++)
					{
					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,"seq%05d     ", i+1);
					else if (i == 2*numSequences-2)
						fprintf (fpAlignment,"root         ");
					else
						fprintf (fpAlignment,"anc%05d     ", i+1);

					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon11 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i,j,numSites)], codon);
							for	(k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else if (i == 2*numSequences-2) /* is root */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon11 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
							for	(k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else /* is ancestral */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon11 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i,j,numSites)], codon);
							for	(k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}

					fprintf (fpAlignment,"\n");
					}	
				/*fprintf (fpAlignment,"\n");*/
				}
			}
		else /* migration */
			{
			if (thereisOutgroup == YES)
				{	
				for (i=0; i<2*numSequences; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodex + m;
						if ((f->label == i) && (f->NetIndex <= numSequences*2-2))
							{
							dem = f->indexOldMigPop;
							break;
							}
						}*/
					
					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;

						if (i < numSequences) /* tip */
							{
							if ((f->label == i) && (f->index <= numSequences*2-2)) 
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == 2*numSequences-1) /* root */
							{
							if (f->class == 5)
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == numSequences) /* outgroup */
							{
							dem = 0;
							break;
							}
						else /* ancestral */
							{
							if (f->label == i-1)
								{
								dem = f->indexOldMigPop;
								if (f->class != 4)
									{
									fprintf (stderr, "\n Warning in PrintAncestralSequences_C. f->label = %d, f->class = %d \n", f->label, f->class);
									exit(-1);
									}
								break;
								}
							}
						}



					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,"s%05d_p%d    ", i+1, dem);
					else if (i == numSequences)
						fprintf(fpAlignment, "outgrp_p0    ");
					else if (i == 2*numSequences-1)
						fprintf (fpAlignment,"root_p%d      ",dem);
					else
						fprintf (fpAlignment,"a%05d_p%d    ", i+1, dem);
		


					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else if (i == numSequences) /* is outgroup */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(outgroupLabel,j,numSites)] > 60 || matrixC[pos(outgroupLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(outgroupLabel,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else if (i == 2*numSequences-1) /* is root */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else  /* is ancestral */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i-1,j,numSites)] > 60 || matrixC[pos(i-1,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon10 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i-1,j,numSites)], codon);
							for (k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}


					fprintf (fpAlignment,"\n");
					}	
				/*fprintf (fpAlignment,"\n");*/
				}
			else
				{
				for (i=0; i<2*numSequences-1; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodex + m;
						if ((f->label == i) && (f->NetIndex <= numSequences*2-2))
							{
							dem = f->indexOldMigPop;
							break;
							}
						}*/
					
					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;

						if (i < numSequences) /* tip */
							{
							if ((f->label == i) && (f->index <= numSequences*2-2)) 
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else if (i == 2*numSequences-2) /* root */
							{
							if (f->class == 5)
								{
								dem = f->indexOldMigPop;
								break;
								}
							}
						else /* ancestral */
							{
							if (f->label == i)
								{
								dem = f->indexOldMigPop;
								if (f->class != 4)
									{
									fprintf (stderr, "\n Warning in PrintAncestralSequences_C. f->label = %d, f->class = %d \n", f->label, f->class);
									exit(-1);
									}
								break;
								}
							}
						}




					if (i < numSequences) /* is tip */
						fprintf (fpAlignment,"s%05d_p%d    ", i+1, dem);
					else if (i == 2*numSequences-2)
						fprintf (fpAlignment,"root_p%d      ", dem);
					else
						fprintf (fpAlignment,"a%05d_p%d    ", i+1, dem);


					if (i < numSequences) /* is tip */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon11 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i,j,numSites)], codon);
							for	(k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else if (i == 2*numSequences-2) /* is root */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon11 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
							for	(k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}
					else /* is ancestral */
						{
						for (j=1; j<=numSites; j++)
							{
							if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
								{
								fprintf (stderr, "\n stop codon11 \n");
								exit(-1);
								}
							number_to_codon(matrixC[pos(i,j,numSites)], codon);
							for	(k = 0; k < 3; k++)
								fprintf (fpAlignment, "%c", codon[k]);
							}
						}

					fprintf (fpAlignment,"\n");

					}	
				}
			/*fprintf (fpAlignment,"\n");*/
			}
		}
	else /* There are recombinations */
		{
		if (thereisOutgroup == YES)
			fprintf(fpAlignment,"Dimensions ntax=%d nchar=%d;\n", numSequences+2, numNuc);
		else
			fprintf(fpAlignment,"Dimensions ntax=%d nchar=%d;\n", numSequences+1, numNuc);
		
		fprintf(fpAlignment,"	Format datatype=nucleotide gap=- missing=? matchchar=.;\n");
		fprintf(fpAlignment,"	Matrix\n");

		if (doMigration == NO)
			{
			if (thereisOutgroup == YES)
				{
				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 2)
						outgroupLabel = f->label;
					if (f->class == 5)
						rootLabel = f->label;
					}				
	
				for (i=0; i<2*numSequences; i++)
					{
					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,"seq%05d     ", i+1);
						a++;
						}
					if (i == numSequences)
						{
						fprintf(fpAlignment, "outgroup     ");
						a++;
						}
					if (i == 2*numSequences-1)
						{
						fprintf (fpAlignment,"root         ");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon14 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(i,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == numSequences) /* is outgroup */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(outgroupLabel,j,numSites)] > 60 || matrixC[pos(outgroupLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon14 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(outgroupLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-1) /* is root */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon14 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}

						}
					}	
				/*fprintf (fpAlignment,"\n");*/
				}
			else
				{

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 5)
						rootLabel = f->label;
					}

				for (i=0; i<2*numSequences-1; i++)
					{
					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,"seq%05d     ", i+1);
						a++;
						}
					if (i == 2*numSequences-2)
						{
						fprintf (fpAlignment,"root         ");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon15 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(i,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-2) /* is root */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon15 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}

						}
					}	
				/*fprintf (fpAlignment,"\n");*/
				}
			}
		else /* migration */
			{
			if (thereisOutgroup == YES)
				{	

				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 2)
						outgroupLabel = f->label;
					if (f->class == 5)
						rootLabel = f->label;
					}	


				for (i=0; i<2*numSequences; i++)
					{

					/*for (m = 0; m < numNodex; m++)
						{
						f = nodex + m;
						if (f->NetIndex < numSequences*/ /*|| i == 2*numSequences-1*//*
							if ( f->label == i)
								{
								dem = f->indexOldMigPop;*/
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								/*break;
								}
						}*/

					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;
						if (f->index < numSequences /*|| i == 2*numSequences-1*/)
							if (f->label == i)
								{
								dem = f->indexOldMigPop;
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								break;
								}
						}



					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,"s%05d_p%d    ", i+1, dem);
						a++;
						}
					if (i == numSequences)
						{
						fprintf(fpAlignment, "outgrp_p0    ");
						a++;
						}
					if (i == 2*numSequences-1)
						{
						fprintf (fpAlignment,"root         ");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{

						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon16 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(i,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == numSequences) /* is outgroup */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(outgroupLabel,j,numSites)] > 60 || matrixC[pos(outgroupLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon16 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(outgroupLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-1) /* is root */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon16 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}

						}
					}	
				/*fprintf (fpAlignment,"\n");*/
				}
			else
				{


				for (m = 0; m < nextAvailable; m++)
					{
					f = nodes + m;
					/*fprintf (stderr, "\n f->label = %d   f->class = %d \n", f->label, f->class);*/
					if (f->class == 5)
						rootLabel = f->label;
					}

				for (i=0; i<2*numSequences-1; i++)
					{
					/*for (m = 0; m < numNodex; m++)
						{
						f = nodex + m;
						if (f->NetIndex < numSequences*/ /*|| i == 2*numSequences-1*//*)
							if (f->label == i)
								{
								dem = f->indexOldMigPop;*/
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								/*break;
								}
						}*/
					
					for (m = 0; m < nextAvailable; m++)
						{
						f = nodes + m;
						if (f->index < numSequences /*|| i == 2*numSequences-1*/)
							if ( f->label == i)
								{
								dem = f->indexOldMigPop;
								/*if (i == 2*numSequences-1)
									fprintf(stderr,"\n\n ****** f->label = %d,    f->NetIndex = %d,    f->indexOldMigPop = %d     \n", f->label, f->NetIndex, f->indexOldMigPop);*/
								break;
								}
						}



					a = 0;
					if (i < numSequences) /* is tip */
						{
						fprintf (fpAlignment,"s%05d_p%d    ", i+1, dem);
						a++;
						}
					if (i == 2*numSequences-2)
						{
						fprintf (fpAlignment,"root         ");
						a++;
						}
					/*else
						fprintf (fpAlignment,"anc%05d  ", i+1);*/
					if (a != 0)
						{
						if (i < numSequences) /* is tip */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(i,j,numSites)] > 60 || matrixC[pos(i,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon17 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(i,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}
						if (i == 2*numSequences-2) /* is root */
							{
							for (j=1; j<=numSites; j++)
								{
								if (matrixC[pos(rootLabel,j,numSites)] > 60 || matrixC[pos(rootLabel,j,numSites)] > 60)
									{
									fprintf (stderr, "\n stop codon17 \n");
									exit(-1);
									}
								number_to_codon(matrixC[pos(rootLabel,j,numSites)], codon);
								for (k = 0; k < 3; k++)
									fprintf (fpAlignment, "%c", codon[k]);
								}
							fprintf (fpAlignment,"\n");
							}


						}
					}	
				/*fprintf (fpAlignment,"\n");*/
				}
			}
		}
	PrintNEXUS_end();
}










/***************************** PrintNexus1 *******************************/
/* Prints information header nexus output files*/
static void PrintNEXUS_initial ()
{
	char *date;				/* define date */
	time_t now;
	
	now=time(NULL);
	date= ctime(&now);

	fprintf(fpAlignment,"#NEXUS\n\n");
	fprintf(fpAlignment,"[\n\n");
	fprintf(fpAlignment,"Generated by %s on %s\n", PROGRAM_NAME, date);
	fprintf(fpAlignment,"%s - Simulating Codon Sequences with Total Recombination\n", PROGRAM_NAME);
	fprintf(fpAlignment, "Version %s", VERSION_NUMBER);
	fprintf(fpAlignment, "\nWritten by Miguel Arenas (miguelab@uvigo.es) and David Posada (dposada@uvigo.es)");
	fprintf(fpAlignment, "\nArenas, M. and Posada, D. 2010. Coalescent simulation of intracodon recombination. Genetics, 184(2):429-437.\n\n");

	fprintf(fpAlignment,"]\n\n\n\n");
	fprintf(fpAlignment,"Begin data;\n");
}

/***************************** PrintNexus2 *******************************/
/* Prints information at the end of nexus output files*/
static void PrintNEXUS_end ()
{
	fprintf(fpAlignment,"	;\n");
	fprintf(fpAlignment,"End;\n");
}






/***************************** PrintRunSettings *******************************/
/* Prints a summary of run settings */

static void	PrintRunSettings (FILE *filep, long int seed)
	{
		int i, j;
		double stdErrorErep, stdErrGMRCA, stdErrTrep;

		
		stdErrorErep = stdErrGMRCA = stdErrTrep = 0.0;		

		/* A–adido Miguel: Para que escriba por la pantalla del master */
		#ifdef MPI
			if (rank==root)
				{
		#endif
		
		if (doSettingsFile == YES)
			filep = fpSettings;
		/*fprintf(fpTrees, "Dataset %d \n", replicate+1);*/

		fprintf (filep, "\n\nRun settings\n-------------------------");
		fprintf (filep, "\n[Assumptions in brackets]");
		fprintf (filep, "\n\nSimulations");
		fprintf (filep, "\nSeed                                =  %-3ld", seed);
		fprintf (filep, "\nNumber replicate data sets          =  %-3d", numDataSets);
		fprintf (filep, "\nNumber of sequences                 =  %-3d", numSequences);
		fprintf (filep, "\nNumber of sites (bp or codons)      =  %-3d", numSites);
		if (Nscaling == 1)
			{
			fprintf (filep, "\nHaploid data");
			}
		if (Nscaling == 2)
			{
			fprintf (filep, "\nDiploid data");
			}
		/*fprintf (filep, "\nEffective population size           =  %-3ld", N);*/
		if (doDatedTips == YES)
			{
			fprintf (filep, "\nDated tips");
			fprintf (filep, "\n Generation time                    =   %3.2f", generationTime);
			fprintf (filep, "\n Sample    Time    Size   Sequences");
			for (i=0; i<numTipDates; i++)
				fprintf (filep, "\n  %2d      %7.2f   %2d     %d-%d", 
				i+1, datedSample[i].time, datedSample[i].size, datedSample[i].member[0], datedSample[i].member[datedSample[i].size-1]);	
			}


		/* Global events */
		fprintf (filep, "\n\nGlobal Data Sets");
		j = cumNumCA+cumNumRE+cumNumMIG+cumNumCONV;
		for (i = 0; i < numDataSets; i++)
			{
			varEvent[i] = (varEvent[i]-(j/numDataSets))*(varEvent[i]-(j/numDataSets));
			varianceErep = varianceErep + (varEvent[i]*1.0);
			}
		varianceErep = varianceErep/numDataSets;
		stdErrorErep = sqrt (varianceErep/numDataSets);
		if (doMigration == YES)
			{
			if (doConvergDemes == YES)
				fprintf (filep, "\nTotal number of events \n (CA+RE+MI+CONV)                    =  %d \n Mean events per replicate          =  %d (stdErr = %3.2f)", j, j/numDataSets, stdErrorErep);
			else
				fprintf (filep, "\nTotal number of events (CA+RE+MI)   =  %d \n Mean events per replicate          =  %d (stdErr = %3.2f)", j, j/numDataSets, stdErrorErep);
			}
		else
			fprintf (filep, "\nTotal number of events (CA+RE)      =  %d \n Mean events per replicate          =  %d (stdErr = %3.2f)", j, j/numDataSets, stdErrorErep);
		/*fprintf (filep, "\nTotal number of recombinations      =  %d", cumNumRE);*/
		fprintf (filep, "\nMean number of coalescence events   =  %3.2f", meanNumCA);

		/* Recombination */
		fprintf (filep, "\n\nRecombination");
		fprintf (filep, "\nRecombination rate                  =  %2.1e", recombinationRate);
		fprintf (filep, "\nRho (2*Nscaling*NrL)                =  %3.2f", rho);
		fprintf (filep, "\nExpected number of rec events       =  %3.2f", expNumRE);
		fprintf (filep, "\n  [N=con, no mig]");
		if (doCountsForExpNumRec == YES)
			fprintf (filep, "\nCountable recombination events      =  %3.2f", meanNumREtc);
		if (doFixNumRecEvents == NO)
			fprintf (filep, "\n#Reps with 0 rec events             =  %d", zeroRec);
		if (doFixNumRecEvents == YES)
			fprintf (filep, "\nFixed number of recombination events =  %-3d", fixedNumRecEvents);

		fprintf (filep, "\nMean number of recombination events =  %3.2f", meanNumRE);
		if ((meanNumRE > 0.00) && (doCodonModel == YES))
			{
			fprintf (filep, "\n Countable Mean number of recombinations breaking codons              =  %3.2f", meanNumREbreakCod);
			fprintf (filep, "\n Countable Expected mean number of recombinations breaking codons     =  %3.2f", meanNumREtc*2/3);
			fprintf (filep, "\n Countable Mean of stop codons generated/restarted by IntraCodon Rec  =  %3.2f", meanNumStopCodonREC);
			
			/*fprintf (filep, "\n meanNumEqual2                      =  %3.2f", meanNumEqual2);
			fprintf (filep, "\n meanNumEqual1                      =  %3.2f", meanNumEqual1);
			fprintf (filep, "\n meanNumDifCodSameAA                =  %3.2f", meanNumDifCodSameAA);
			fprintf (filep, "\n meanNumDifCodDifAA                 =  %3.2f", meanNumDifCodDifAA);*/
			fprintf (filep, "\n Countable Mean number of recombinations with 0 nonsynonymous changes =  %3.2f (%3.2f%%)", meanNumNonSyn0, 1.00*meanNumNonSyn0*100/meanNumREtc);
			fprintf (filep, "\n Countable Mean number of recombinations with 1 nonsynonymous changes =  %3.2f (%3.2f%%)", meanNumNonSyn1, 1.00*meanNumNonSyn1*100/meanNumREtc);
			fprintf (filep, "\n Countable Mean number of recombinations with 2 nonsynonymous changes =  %3.2f (%3.2f%%)", meanNumNonSyn2, 1.00*meanNumNonSyn2*100/meanNumREtc);
			fprintf (filep, "\n Countable Total number of nonsynonymous changes by intracodon recombination =  %3.2f", 1.0*(meanNumNonSyn1+meanNumNonSyn2+meanNumNonSyn2));
			}
		
		fprintf (filep, "\n\nMean number of mutation events      =  %3.2f", meanNumMU);
		if (doCodonModel == YES)
			{
			fprintf (filep, "\n  synonymous changes =  %3.2f (%3.2f%%)", meanNumMU_S, 1.00*meanNumMU_S*100/meanNumMU);
			fprintf (filep, "\n  nonsynonymous changes =  %3.2f (%3.2f%%)", meanNumMU_NS, 1.00*meanNumMU_NS*100/meanNumMU);
			}
		/* Average tMRCA and tGMRCA inter all replicates */
		for (i = 0; i < numDataSets; i++)
			{			
			#ifdef HUDSON_UNITS
				varTimeGMRCA[i] = (varTimeGMRCA[i]-(counterTime/(2*Nscaling*N*numDataSets)))*(varTimeGMRCA[i]-(counterTime/(2*Nscaling*N*numDataSets))); /* Hudson units(/2*Nscaling*N) */
			#else
				varTimeGMRCA[i] = (varTimeGMRCA[i]-(counterTime/numDataSets))*(varTimeGMRCA[i]-(counterTime/numDataSets)); 
			#endif
			varianceGMRCArep = varianceGMRCArep + (varTimeGMRCA[i]*1.0);
			}
		varianceGMRCArep = varianceGMRCArep/numDataSets;
		stdErrGMRCA = sqrt(varianceGMRCArep/numDataSets);
		for (i = 0; i < numDataSets; i++)
			{
			varTimeT[i] = (varTimeT[i]-(countTMRCAReps/numDataSets))*(varTimeT[i]-(countTMRCAReps/numDataSets));
			varianceTrep = varianceTrep + (varTimeT[i]*1.0);
			}
		varianceTrep = varianceTrep/numDataSets;
		stdErrTrep = sqrt(varianceTrep/numDataSets);
		#ifdef HUDSON_UNITS
			fprintf(filep, "\nAverage time to GMRCA               =  %3.2f (stdErr = %3.2f)  (Hd units, 1/2*Nscaling*N)", (counterTime/(numDataSets*2*Nscaling*N)), stdErrGMRCA);
		#else
			fprintf(filep, "\nAverage time to GMRCA               =  %3.2f (stdErr = %3.2f)", counterTime/numDataSets, stdErrGMRCA);
		#endif
		fprintf(filep, "\nAverage time to MRCA                =  %3.2f (stdErr = %3.2f)", countTMRCAReps/numDataSets, stdErrTrep);
		fprintf(filep, "\nExpected time to the MRCA           =  %3.2f (var = %3.2f)", expTMRCA, expVarTMRCA);
		fprintf(filep, "\n  [N=con, no rec, no mig]");
				
		/* Demographics */
		fprintf (filep, "\n\nDemographics");
		fprintf (filep, "\nEffective population size           =  %-3d", N);
		if (doExponential == YES)
			fprintf (filep, "\nExponential growth rate             =  %+2.1e", growthRate);
		if (doDemographics == YES)
			{
			fprintf (filep, "\nDemographic periods:");
			fprintf (filep, "\n Period  Nbegin    Nend   Duration   Growth rate");
			for (i=1; i<=numPeriods; i++)
				fprintf (filep, "\n  %d   %7d  %7d    %7d      %+2.1e", 
			i, Nbegin[i], Nend[i], cumDuration[i]-cumDuration[i-1], periodGrowth[i]);	
			}
			
		/* Migration */
		if (doMigration == YES)
			{
			fprintf (filep, "\n\nIsland migration model");
			fprintf (filep, "\nMigration rate                      =  %2.1e", migrationRate);
			fprintf (filep, "\nExpected FST = 1/(2*Nscaling*Nm+1)  =  %2.1e", 1/(2.0*Nscaling*N*migrationRate+1));
			fprintf (filep, "\n  [N=con, no rec]");
			fprintf (filep, "\nInitial number of demes             =  %-3d", numPopulations);
			fprintf (filep, "\n Deme         size");
			for (i=1; i<=numPopulations; i++)
				fprintf (filep, "\n  %d            %-3d", i, initPopulation[i]);
			fprintf (filep, "\n\nMean number of migration events     =  %3.2f", meanNumMIG);
			if (doConvergDemes == YES) /* Convergence of demes */
				{
				fprintf (filep, "\n\nConvergence of demes");
				fprintf (filep, "\n Event    Time       Deme 1   Deme 2");
				for (i=1; i<=numConvergDemes; i++)
					{
					if (convDemTimes_old[i] < 100)
						fprintf (filep, "\n   %d      %3.2f        %-4d     %-2d", i, convDemTimes_old[i], deme_a_old[i], deme_b_old[i]);
					if (convDemTimes_old[i] < 1000 && convDemTimes_old[i] >= 100)
						fprintf (filep, "\n   %d      %3.2f       %-4d     %-2d", i, convDemTimes_old[i], deme_a_old[i], deme_b_old[i]);
					if (convDemTimes_old[i] < 10000 && convDemTimes_old[i] >= 1000)
						fprintf (filep, "\n   %d      %3.2f      %-4d     %-2d", i, convDemTimes_old[i], deme_a_old[i], deme_b_old[i]);
					if (convDemTimes_old[i] < 100000 && convDemTimes_old[i] >= 10000)
						fprintf (filep, "\n   %d      %3.2f     %-4d     %-2d", i, convDemTimes_old[i], deme_a_old[i], deme_b_old[i]);
					if (convDemTimes_old[i] >= 100000)
						fprintf (filep, "\n   %d      %3.2f    %-4d     %-2d", i, convDemTimes_old[i], deme_a_old[i], deme_b_old[i]);
					}
				fprintf (filep, "\nMean number of convergence demes\n events                             =  %3.2f", meanNumCONV);
				}
			}

		/* Other Settings */
			fprintf (filep, "\n\nOther Settings");
		if (thereisOutgroup == NO)
			fprintf (filep, "\nSimulate outgroup                   =  NO");
		else
			{
			fprintf (filep, "\nSimulate outgroup                   =  YES");
			fprintf (filep, "\n   branch length                    =  %3.2f", outgroupBranchLength);
			}	
	
		if (doPrintAncestralSequences == YES)
			fprintf (filep, "\nPrint ancestral states              =  YES");
		else
			fprintf (filep, "\nPrint ancestral states              =  NO");

		if (doMRCAFile == YES)
			fprintf (filep, "\nMRCA sequence                       =  specified by the user");
		else
			fprintf (filep, "\nMRCA sequence                       =  simulated");
			
			
		
		/* Substitution model */
		if (doCodonModel == YES)
			fprintf (filep, "\n\n\nCodon replacement model\n------------------------");
		else
			fprintf (filep, "\n\n\nNucleotide substitution model\n-----------------------------");
		fprintf (filep, "\nSubstitution rate                   =  %2.1e", mutationRate);
		fprintf (filep, "\nExpected theta (2*Nscaling*NuL) per site [ISM] =  %5.4f", 2.0 * Nscaling * N * mutationRate * numSites);
		
		if (doCodonModel == YES)
			{
			fprintf (filep, "\nBase frequencies (piA piC piG piT) for codon positions:");
			fprintf (filep, "\n   first  =   %3.2f %3.2f %3.2f %3.2f", p_i_codon[0], p_i_codon[1], p_i_codon[2], p_i_codon[3]);
			fprintf (filep, "\n   second =   %3.2f %3.2f %3.2f %3.2f", p_i_codon[4], p_i_codon[5], p_i_codon[6], p_i_codon[7]);
			fprintf (filep, "\n   third  =   %3.2f %3.2f %3.2f %3.2f", p_i_codon[8], p_i_codon[9], p_i_codon[10], p_i_codon[11]);
			
			if (doOmegaCat == NO && doM8 == NO && doM7 == NO)
				fprintf (filep, "\nOmega (dN/dS ratio)                 =  %3.2f", OmegaInit);
			if (doOmegaCat == YES)
				{
				if (doM1 == YES)
					{
					fprintf (filep, "\n\nM1 codon model:\n  P0 = %3.2f, W0 = %3.2f, P1 = %3.2f\n", M1_P0_omeg0, M1_omega0, M1_P1_omeg1);
					}
				else
					{
					fprintf (filep, "\nOmega (dN/dS ratio) in %d user-categories: \n Category      Omega      Probability", numOmegaCat);
					for (j=1; j<=numOmegaCat; j++)
						fprintf(filep, "\n    %d          %3.2f           %3.2f", j, omegaVal[j], omegaProb[j]);
					}
				}
			if (doOmegaRateHetCont == YES)
				fprintf (filep, "\nOmega variation among sites\n alpha                              =  %3.2f", OmegaRateHet);
			if (doOmegaRateHetDisc == YES)
				{
				fprintf (filep, "\nOmega discrete variation among sites (%d categories)\n alpha (gamma shape)                =  %3.2f", numOmegaCat, OmegaRateHet);
				fprintf (filep, "\n Category      Omega      Probability");
				for (j=1; j<=numOmegaCat; j++)
					fprintf(filep, "\n    %d          %3.2f           %3.2f", j, OmegaInit*gammaRates[j-1], (double)1/numOmegaCat);
				/*fprintf (filep, "\nGamma values obtained               =");*/ /* Active to see the gamma values for omega */
				/*for (j=1; j<=numOmegaCat; j++)
					fprintf(filep, "  %3.3f", gammaRates[j-1]);*/
				}
			if (doM8 == YES)
				{
				fprintf (filep, "\n\nM8 codon model:\n  P0 = %3.2f P1 = %3.2f\n", M8_P0_beta, M8_P1_omega);
				fprintf (filep, "  Beta distribution (P0) p = %3.2f q = %3.2f; Omega (P1) = %3.2f\n", M8_p_beta, M8_q_beta, M8_omegaP1);
				}
			if (doM7 == YES)
				{
				fprintf (filep, "\n\nM7 codon model:\n  Beta distribution, p = %3.2f q = %3.2f\n", M7_p_beta, M7_q_beta);
				}
			}
		else
			{
			fprintf (filep, "\nBase frequencies (piA piC piG piT)  =  %3.2f %3.2f %3.2f %3.2f", p_i[0], p_i[1], p_i[2], p_i[3]);
			}
		if (doHKY == YES || doCodon_HKY == YES)
			fprintf (filep, "\nTransition/transversion ratio       =  %3.2f (kappa = %3.2f)", titv, kappa);
		if (doGTR == YES || doCodon_GTR == YES)
			fprintf (filep, "\nR-matrix                            =  %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f", 
				Rmat[0], Rmat[1], Rmat[2], Rmat[3], Rmat[4], Rmat[5]);
		if (doGTnR == YES || doCodon_NGTR == YES)
			fprintf (filep, "\nnR-matrix                           =  %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f\n                                       %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f ", 
				NRmat[0], NRmat[1], NRmat[2], NRmat[3], NRmat[4], NRmat[5], NRmat[6], NRmat[7], NRmat[8], NRmat[9], NRmat[10], NRmat[11]);

		fprintf (filep, "\nProportion of invariable sites      =  %3.2f ", pinv);
		if (doRateHet == YES)
			fprintf (filep, "\nRate variation among sites\n  alpha (gamma shape)               =  %3.2f ", alpha);
		else
			fprintf (filep, "\nNo rate variation among sites");
		
		if (doCodonModel == YES)
			{
			if (doCodon_HKY == YES)
				{
				if (freqNumber == 12 && doM0 == YES)
					{
					sprintf(model,"GY94xHKY_3x4 M0");
					/*fprintf(filep,"\nGY94xHKY_3x4 M0\n");*/
					}
				if (freqNumber == 12 && doM1 == YES)
					{
					sprintf(model,"GY94xHKY_3x4 M1");
					}
				if (freqNumber == 12 && doM7 == YES)
					{
					sprintf(model,"GY94xHKY_3x4 M7");
					}
				if (freqNumber == 12 && doM8 == YES)
					{
					if (M8_omegaP1 > 1)
						sprintf(model,"GY94xHKY_3x4 M8");
					if (M8_omegaP1 == 1.0)
						sprintf(model,"GY94xHKY_3x4 M8a");
					if (M8_omegaP1 < 1)
						sprintf(model,"GY94xHKY_3x4 M8b");
					}
				if (freqNumber == 12 && doOmegaCat == YES && doM1 == NO)
					{
					sprintf(model,"GY94xHKY_3x4 by categories");
					}
				if (freqNumber == 12 && doOmegaRateHetCont == YES)
					{
					sprintf(model,"GY94xHKY_3x4 with variable omega by continuous alpha (gamma shape)");
					}
				if (freqNumber == 12 && doOmegaRateHetDisc == YES)
					{
					sprintf(model,"GY94xHKY_3x4 with variable omega by discrete alpha (gamma shape)");
					}

				if (freqNumber == 4 && doM0 == YES)
					{
					sprintf(model,"GY94xHKY_1x4 M0");
					}
				if (freqNumber == 4 && doM1 == YES)
					{
					sprintf(model,"GY94xHKY_1x4 M1");
					}
				if (freqNumber == 4 && doM7 == YES)
					{
					sprintf(model,"GY94xHKY_1x4 M7");
					}
				if (freqNumber == 4 && doM8 == YES)
					{
					if (M8_omegaP1 > 1)
						sprintf(model,"GY94xHKY_1x4 M8");
					if (M8_omegaP1 == 1.0)
						sprintf(model,"GY94xHKY_1x4 M8a");
					if (M8_omegaP1 < 1)
						sprintf(model,"GY94xHKY_1x4 M8b");
					}
				if (freqNumber == 4 && doOmegaCat == YES && doM1 == NO)
					{
					sprintf(model,"GY94xHKY_1x4 by categories");
					}
				if (freqNumber == 4 && doOmegaRateHetCont == YES)
					{
					sprintf(model,"GY94xHKY_1x4 with variable omega by continuous alpha (gamma shape)");
					}
				if (freqNumber == 4 && doOmegaRateHetDisc == YES)
					{
					sprintf(model,"GY94xHKY_1x4 with variable omega by discrete alpha (gamma shape)");
					}

				}
			if (doCodon_GTR == YES)
				{
				if (freqNumber == 12 && doM0 == YES)
					{
					sprintf(model,"GY94xGTR_3x4 M0");
					/*fprintf(filep,"\nGY94xGTR_3x4 M0\n");*/
					}
				if (freqNumber == 12 && doM1 == YES)
					{
					sprintf(model,"GY94xGTR_3x4 M1");
					}
				if (freqNumber == 12 && doM7 == YES)
					{
					sprintf(model,"GY94xGTR_3x4 M7");
					}
				if (freqNumber == 12 && doM8 == YES)
					{
					if (M8_omegaP1 > 1)
						sprintf(model,"GY94xGTR_3x4 M8");
					if (M8_omegaP1 == 1.0)
						sprintf(model,"GY94xGTR_3x4 M8a");
					if (M8_omegaP1 < 1)
						sprintf(model,"GY94xGTR_3x4 M8b");
					}
				if (freqNumber == 12 && doOmegaCat == YES && doM1 == NO)
					{
					sprintf(model,"GY94xGTR_3x4 by categories");
					}
				if (freqNumber == 12 && doOmegaRateHetCont == YES)
					{
					sprintf(model,"GY94xGTR_3x4 with variable omega by continuous alpha (gamma shape)");
					}
				if (freqNumber == 12 && doOmegaRateHetDisc == YES)
					{
					sprintf(model,"GY94xGTR_3x4 with variable omega by discrete alpha (gamma shape)");
					}

				if (freqNumber == 4 && doM0 == YES)
					{
					sprintf(model,"GY94xGTR_1x4 M0");
					}
				if (freqNumber == 4 && doM1 == YES)
					{
					sprintf(model,"GY94xGTR_1x4 M1");
					}
				if (freqNumber == 4 && doM7 == YES)
					{
					sprintf(model,"GY94xGTR_1x4 M7");
					}
				if (freqNumber == 4 && doM8 == YES)
					{
					if (M8_omegaP1 > 1)
						sprintf(model,"GY94xGTR_1x4 M8");
					if (M8_omegaP1 == 1.0)
						sprintf(model,"GY94xGTR_1x4 M8a");
					if (M8_omegaP1 < 1)
						sprintf(model,"GY94xGTR_1x4 M8b");
					}
				if (freqNumber == 4 && doOmegaCat == YES && doM1 == NO)
					{
					sprintf(model,"GY94xGTR_1x4 by categories");
					}
				if (freqNumber == 4 && doOmegaRateHetCont == YES)
					{
					sprintf(model,"GY94xGTR_1x4 with variable omega by continuous alpha (gamma shape)");
					}
				if (freqNumber == 4 && doOmegaRateHetDisc == YES)
					{
					sprintf(model,"GY94xGTR_1x4 with variable omega by discrete alpha (gamma shape)");
					}
				}
			if (doCodon_NGTR == YES)
				{
				if (freqNumber == 12 && doM0 == YES)
					{
					sprintf(model,"GY94xGTnR_3x4 M0");
					/*fprintf(filep,"\nGY94xGTnR_3x4 M0\n");*/
					}
				if (freqNumber == 12 && doM1 == YES)
					{
					sprintf(model,"GY94xGTnR_3x4 M1");
					}
				if (freqNumber == 12 && doM7 == YES)
					{
					sprintf(model,"GY94xGTnR_3x4 M7");
					}
				if (freqNumber == 12 && doM8 == YES)
					{
					if (M8_omegaP1 > 1)
						sprintf(model,"GY94xGTnR_3x4 M8");
					if (M8_omegaP1 == 1.0)
						sprintf(model,"GY94xGTnR_3x4 M8a");
					if (M8_omegaP1 < 1)
						sprintf(model,"GY94xGTnR_3x4 M8b");
					}
				if (freqNumber == 12 && doOmegaCat == YES && doM1 == NO)
					{
					sprintf(model,"GY94xGTnR_3x4 by categories");
					}
				if (freqNumber == 12 && doOmegaRateHetCont == YES)
					{
					sprintf(model,"GY94xGTnR_3x4 with variable omega by continuous alpha (gamma shape)");
					}
				if (freqNumber == 12 && doOmegaRateHetDisc == YES)
					{
					sprintf(model,"GY94xGTnR_3x4 with variable omega by discrete alpha (gamma shape)");
					}

				if (freqNumber == 4 && doM0 == YES)
					{
					sprintf(model,"GY94xGTnR_1x4 M0");
					}
				if (freqNumber == 4 && doM1 == YES)
					{
					sprintf(model,"GY94xGTnR_1x4 M1");
					}
				if (freqNumber == 4 && doM7 == YES)
					{
					sprintf(model,"GY94xGTnR_1x4 M7");
					}
				if (freqNumber == 4 && doM8 == YES)
					{
					if (M8_omegaP1 > 1)
						sprintf(model,"GY94xGTnR_1x4 M8");
					if (M8_omegaP1 == 1.0)
						sprintf(model,"GY94xGTnR_1x4 M8a");
					if (M8_omegaP1 < 1)
						sprintf(model,"GY94xGTnR_1x4 M8b");
					}
				if (freqNumber == 4 && doOmegaCat == YES && doM1 == NO)
					{
					sprintf(model,"GY94xGTnR_1x4 by categories");
					}
				if (freqNumber == 4 && doOmegaRateHetCont == YES)
					{
					sprintf(model,"GY94xGTnR_1x4 with variable omega by continuous alpha (gamma shape)");
					}
				if (freqNumber == 4 && doOmegaRateHetDisc == YES)
					{
					sprintf(model,"GY94xGTnR_1x4 with variable omega by discrete alpha (gamma shape)");
					}
				}
			}
		else
			{
			if (doGTR == YES)
				sprintf (model,"GTR");
			else if (doGTnR == YES)
				sprintf (model,"GTnR");
			else
				{
				if (equalBaseFreq == YES)
					{
					if (titv == 0.5)
						sprintf (model,"JC");
					else
						sprintf (model,"K80");
					}
				else
					{
					if (titv == 0.5)
						sprintf (model,"F81");
					else
						sprintf (model,"HKY");
					}
				}
			}
		
		if (doRateHet == NO && pinv == 0)
			fprintf (filep,"\n\nThese settings correspond to the %s model\n", model);
		if (doRateHet == NO && pinv != 0)
			fprintf (filep,"\n\nThese settings correspond to the %s + I model\n", model);
		if (doRateHet == YES && pinv == 0)
			fprintf (filep,"\n\nThese settings correspond to the %s + G model\n", model);
		if (doRateHet == YES && pinv != 0)
			fprintf (filep,"\n\nThese settings correspond to the %s + G + I model\n", model);
	

	#ifdef MPI
			}
	#endif
	}






/******************** ReadParametersFromCommandLine *************************/
static void ReadParametersFromCommandLine (int argc,char **argv)
{
	int		i, j, k, h, l, modelNumber;
	char	flag;
	int		to, from; 
 	float	argument;
	double	sumPi, sumPi_b, sumPi_Cod_first, sumPi_Cod_second, sumPi_Cod_third;
	int		sumPopul;
	
	modelNumber = 0;
	strcpy(alignmentFile, "sequences");
	
	for (i=1; i<argc; i++)
	{
		argv[i]++;
		flag=*argv[i];
		argv[i]++;
		argument = -9999;
		

	/* Used: A B C D E F G H I J K L N M O P R S T U V W X Y Z ? # @ % $ * + = / _      next ideas : ; */

	   switch (toupper(flag))
			{
			case '#':
				argument = atof(argv[i]);
				userSeed = (int) argument;
				if (userSeed < 0) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad user seed for random number (%lu bytes)\n\n", userSeed);
					PrintUsage();
					}
			break;
			case 'S':
				argument = atof(argv[i]);
				numSequences = (int) argument;
				if (numSequences < 1)
				 	{
					fprintf (stderr, "PARAMETER ERROR: Bad sample size (%d)\n\n", numSequences);
					PrintUsage();
					}
			break;
			case 'L':
				argument = atof(argv[i]);
				numSites = (int) argument;
				if (numSites<1) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad number of sites (%d)\n\n", numSites);
					PrintUsage();
					}
			break;
			case 'E':
				argument = atof(argv[i]);
				N = (int) argument;
				if (N < 1) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad effective population size (%d)\n\n", N);
					PrintUsage();
					}
			break;		
			case '=':
				doDatedTips = YES;
	        	argument = atof(argv[i]);
				numTipDates = (int) argument;
				if (numTipDates <= 0) 
					{
					fprintf (stderr, "COMMAND-LINE PARAMETER ERROR: Bad number of sampling dates (%d)\n\n", numTipDates);
					PrintUsage();
					}
				
				datedSample = 	(SampleSt *) calloc(numTipDates, sizeof(SampleSt));
				if (datedSample == NULL)
					{
					fprintf (stderr, "COMMAND-LINE PARAMETER ERROR: Could not allocate sampling dates vectors (%ld)\n", numTipDates * (long) sizeof(SampleSt));
					exit (1);
					}
	
				for (j=0; j<numTipDates; j++)
					{
					argument = atof(argv[++i]);
	        		datedSample[j].time = (float) argument;
					argument = atof(argv[++i]);
        			from = (int) argument;
					argument = atof(argv[++i]);
        			to = (int) argument;
	
					datedSample[j].member = 	(int *) calloc(to-from+1, sizeof(SampleSt));
					if (datedSample[j].member == NULL)
						{
						fprintf (stderr, "COMMAND-LINE PARAMETER ERROR: Could not allocate sampling dates vector[j=%d] (%ld)\n",j, to-from+1 * (long) sizeof(int));
						exit (1);
						}
			
					l = 0;
					for (k=from; k<=to; k++)
						datedSample[j].member[l++] = k;
					
					datedSample[j].size = l;
					}
			break;
			case '/':
				generationTime = atof(argv[i]);
				if (generationTime < 0) 
					{
					fprintf (stderr, "COMMAND-LINE PARAMETER ERROR: Bad generation time (%f)\n\n", generationTime);
					PrintUsage();
					}
			break;
			case 'R':
	        	recombinationRate = atof(argv[i]);
				if (recombinationRate < 0) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad recombination rate (%f)\n\n", recombinationRate );
					PrintUsage();
					}
			break;
			case 'U':
	        	mutationRate = atof(argv[i]);
				if (mutationRate < 0) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad mutation rate (%f)\n\n", mutationRate);
					PrintUsage();
					}
			break;
			case '_':
				argument = atof(argv[i]);
				Nscaling = (int) argument;
				if (Nscaling > 10 || Nscaling < 1) /* max 10 */
					{
					fprintf (stderr, "PARAMETER ERROR: Bad  haplid/diploid chosen (%d)\n\n", Nscaling);
					PrintUsage();
					}
				if (Nscaling < 1 || Nscaling > 2)
					{
					fprintf (stderr, "PARAMETER ERROR: Haploid/diplod option (1-2) (%d)\n\n", Nscaling);
					PrintUsage();
					}
			break;	
			case 'M':
				/*doCodonModel = YES;*/
				sumPi = sumPi_b = 0.0;
				
				modelNumber = atof(argv[i]);
				if (modelNumber > 7 || modelNumber < 1) /* max 10 categories */
					{
					fprintf (stderr, "PARAMETER ERROR: Bad model chosen (1-7) (%d)\n\n", modelNumber);
					PrintUsage();
					}
				
				if (modelNumber == 1)		/** omega constant **/
					{
					doOmegaCat = NO;
					doOmegaRateHetCont = NO;
					doOmegaRateHetDisc = NO;
					doM0 = YES;					

					argument = atof(argv[++i]);
					omega = (double) argument;
					if (omega < 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad omega (%3.2f)\n\n", omega);
						PrintUsage();
						}
					OmegaInit = omega;
					}
				else if (modelNumber == 2) /** omega by categories from the user **/
					{					
					argument = atof(argv[++i]);
					numOmegaCat = (int) argument;
					if (numOmegaCat < 0 || numOmegaCat > 10) /* max 10 categories */
						{
						fprintf (stderr, "PARAMETER ERROR: Bad number of omega categories (%d)\n\n", numOmegaCat);
						PrintUsage();
						}
					omegaVal = (double *) calloc((numOmegaCat+1),(long) sizeof(double));
					if (!omegaVal)
						{
						fprintf (stderr, "PARAMETER ERROR: Could not allocate omega values of categories (%lu bytes)\n", numOmegaCat *(long) sizeof(double));
						exit (1);
						}
					omegaProb = (double *) calloc((numOmegaCat+1),(long) sizeof(double));
					if (!omegaProb)
						{
						fprintf (stderr, "PARAMETER ERROR: Could not allocate omega probabilities of categories1 (%lu bytes)\n", numOmegaCat *(long) sizeof(double));
						exit (1);
						}
						
					for (j=1; j<=numOmegaCat; j++)
						{
						argument = atof(argv[++i]);
						omegaVal[j] = (double) argument;
						
						argument = atof(argv[++i]);
						omegaProb[j] = (double) argument;
						
						if (omegaProb[j] > 1)
							{
							fprintf (stderr, "PARAMETER ERROR: Bad number of probabilities of omega categories (%lf)\n\n", omegaProb[j]);
							PrintUsage();
							}
						
						sumPi = sumPi + omegaProb[j];
						}						
					if (sumPi != 1) /* update probabilities of categories */
						{
						if (numOmegaCat == 1)
							omegaProb[1] = 1.0;
						else
							{
							for (j=1; j<=numOmegaCat; j++)
								{
								omegaProb[j]/=sumPi;
								sumPi_b = sumPi_b + omegaProb[j];
								}
							if ((int)sumPi_b != 1)
								{
								fprintf (stderr, "\n ERROR in the sum of probabilities of omega categories");
								exit (-1);
								}
							}
						}
					doOmegaCat = YES;
					doOmegaProb = YES;
					doOmegaRateHetCont = NO;
					doOmegaRateHetDisc = NO;
					}
				else if (modelNumber == 3) /** omega by discrete heterogeneous rate **/ /* max 10 categories */
					{
					argument = atof(argv[++i]);
					numOmegaCat = (int) argument;
					if (numOmegaCat < 0 || numOmegaCat > 10) 
						{
						fprintf (stderr, "PARAMETER ERROR: Bad number of omega categories (%d)\n\n", numOmegaCat);
						PrintUsage();
						}
					argument = atof(argv[++i]);
					OmegaRateHet = (double) argument;
					if (OmegaRateHet <= 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad omega heterogeneous rate (%3.2f)\n\n", OmegaRateHet);
						PrintUsage();
						}
					argument = atof(argv[++i]);
					omega = (double) argument;
					if (omega < 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad omega (%3.2f)\n\n", omega);
						PrintUsage();
						}
					OmegaInit = omega;
					
					doOmegaCat = NO;
					doOmegaRateHetCont = NO;
					doOmegaRateHetDisc = YES;
					}
				else if (modelNumber == 4) /** omega by continuous heterogeneous rate **/
					{
					doOmegaCat = NO;
					doOmegaRateHetCont = YES;
					doOmegaRateHetDisc = NO;
					
					argument = atof(argv[++i]);
					OmegaRateHet = (double) argument;
					if (OmegaRateHet <= 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad omega heterogeneous rate (%3.2f)\n\n", OmegaRateHet);
						PrintUsage();
						}
					argument = atof(argv[++i]);
					omega = (double) argument;
					if (omega < 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad omega (%3.2f)\n\n", omega);
						PrintUsage();
						}
					OmegaInit = omega;
					}
				else if (modelNumber == 5) /** M1 **/
					{
					doM1 = YES;	

					/*doOmegaCat = NO;*/
					doOmegaRateHetCont = NO;
					doOmegaRateHetDisc = NO;
					
					argument = atof(argv[++i]);
					M1_P0_omeg0 = (double) argument;
					if (M1_P0_omeg0 < 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad P0 in model M1 (%3.2f)\n\n", M1_P0_omeg0);
						PrintUsage();
						}

					argument = atof(argv[++i]);
					M1_omega0 = (double) argument;
					if (M1_omega0 < 0.0 || M1_omega0 > 1)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad omega in M1 (%3.2f)\n\n", M1_omega0);
						PrintUsage();
						}
					

					M1_P1_omeg1 = 1 - M1_P0_omeg0;

					doOmegaCat = YES;
					numOmegaCat = 2;

					omegaVal = (double *) calloc((numOmegaCat+1),(long) sizeof(double));
					if (!omegaVal)
						{
						fprintf (stderr, "PARAMETER ERROR: Could not allocate omega values of categories (%lu bytes)\n", numOmegaCat *(long) sizeof(double));
						exit (1);
						}
					omegaProb = (double *) calloc((numOmegaCat+1),(long) sizeof(double));
					if (!omegaProb)
						{
						fprintf (stderr, "PARAMETER ERROR: Could not allocate omega probabilities of categories1 (%lu bytes)\n", numOmegaCat *(long) sizeof(double));
						exit (1);
						}

					omegaProb[1] = M1_P0_omeg0;
					omegaProb[2] = M1_P1_omeg1;
					omegaVal[1] = M1_omega0;
					omegaVal[2] = 1.0;

					}
				else if (modelNumber == 6) /** M7 **/
					{
					doM7 = YES;

					doOmegaCat = NO;
					doOmegaRateHetCont = NO;
					doOmegaRateHetDisc = NO;
				
					argument = atof(argv[++i]);
					M7_p_beta = (double) argument;
					if (M7_p_beta < 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad p in model M7 (%3.2f)\n\n", M7_p_beta);
						PrintUsage();
						}
					
					argument = atof(argv[++i]);
					M7_q_beta = (double) argument;
					if (M7_q_beta < 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad q in model M7 (%3.2f)\n\n", M7_q_beta);
						PrintUsage();
						}

					}
				else if (modelNumber == 7) /** M8 **/
					{
					doM8 = YES;

					doOmegaCat = NO;
					doOmegaRateHetCont = NO;
					doOmegaRateHetDisc = NO;

					argument = atof(argv[++i]);
					M8_P0_beta = (double) argument;
					if (M8_P0_beta < 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad P0 in model M8 (%3.2f)\n\n", M8_P0_beta);
						PrintUsage();
						}
					M8_P1_omega = 1 - M8_P0_beta;

					argument = atof(argv[++i]);
					M8_p_beta = (double) argument;
					if (M8_p_beta < 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad p in model M8 (%3.2f)\n\n", M8_p_beta);
						PrintUsage();
						}
					
					argument = atof(argv[++i]);
					M8_q_beta = (double) argument;
					if (M8_q_beta < 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad q in model M8 (%3.2f)\n\n", M8_q_beta);
						PrintUsage();
						}
				
					argument = atof(argv[++i]);
					M8_omegaP1 = (double) argument;
					if (M8_omegaP1 < 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad omega for P1 in model M8 (%3.2f)\n\n", M8_omegaP1);
						PrintUsage();
						}
					}
				else
					{
					fprintf (stderr, "PARAMETER ERROR: Bad model chosen (1-7) (%d)\n\n", (int) argument);
					PrintUsage();
					}
				sumPi = sumPi_b = 0.0;
			break;
			case 'F':
				freqNumber = atof(argv[i]);
				if (freqNumber == 4)
					{
					p_i[0] = atof(argv[++i]);
					p_i[1] = atof(argv[++i]);
					p_i[2] = atof(argv[++i]);
					p_i[3] = atof(argv[++i]);

					if ((p_i[0] == p_i[1]) && (p_i[1] == p_i[2]) && (p_i[2] == p_i[3]))
						equalBaseFreq = YES;
					else
						equalBaseFreq = NO;
					
					sumPi = p_i[0] + p_i[1] + p_i[2] + p_i[3];
					if (sumPi !=1.0) 
						{
						p_i[0]/=sumPi;
						p_i[1]/=sumPi;
						p_i[2]/=sumPi;
						p_i[3]/=sumPi;
						}
					
					p_i_codon[0]=p_i[0];
					p_i_codon[1]=p_i[1];
					p_i_codon[2]=p_i[2];
					p_i_codon[3]=p_i[3];
					p_i_codon[4]=p_i[0];
					p_i_codon[5]=p_i[1];
					p_i_codon[6]=p_i[2];
					p_i_codon[7]=p_i[3];
					p_i_codon[8]=p_i[0];
					p_i_codon[9]=p_i[1];
					p_i_codon[10]=p_i[2];
					p_i_codon[11]=p_i[3];					
					}
				else if (freqNumber == 12)
					{
					p_i_codon[0] = atof(argv[++i]);
					p_i_codon[1] = atof(argv[++i]);
					p_i_codon[2] = atof(argv[++i]);
					p_i_codon[3] = atof(argv[++i]);
					
					p_i_codon[4] = atof(argv[++i]);
					p_i_codon[5] = atof(argv[++i]);
					p_i_codon[6] = atof(argv[++i]);
					p_i_codon[7] = atof(argv[++i]);
					
					p_i_codon[8] = atof(argv[++i]);
					p_i_codon[9] = atof(argv[++i]);
					p_i_codon[10] = atof(argv[++i]);
					p_i_codon[11] = atof(argv[++i]);
					
					
						
					for (k = 0; k < 12; k++)
						if (p_i_codon[k] < 0 || p_i_codon[k] > 1)
							{
							fprintf (stderr, "PARAMETER ERROR: Bad number of enter frequencies (it must to be between 0 and 1) (%lf)\n\n", p_i_codon[i]);
							PrintUsage();
							}
					equalBaseFreqCod = YES;
					for (k = 1; k < 12; k++)		
						if (p_i_codon[k] != p_i_codon[k-1])
							{
							equalBaseFreqCod = NO;
							break;
							}		
					
					sumPi_Cod_first = p_i_codon[0] + p_i_codon[1] + p_i_codon[2] + p_i_codon[3];
					sumPi_Cod_second = p_i_codon[4] + p_i_codon[5] + p_i_codon[6] + p_i_codon[7];
					sumPi_Cod_third = p_i_codon[8] + p_i_codon[9] + p_i_codon[10] + p_i_codon[11];
					if (sumPi_Cod_first != 1.0) 
						{
						p_i_codon[0]/=sumPi_Cod_first;
						p_i_codon[1]/=sumPi_Cod_first;
						p_i_codon[2]/=sumPi_Cod_first;
						p_i_codon[3]/=sumPi_Cod_first;
						}
					if (sumPi_Cod_second != 1.0) 
						{
						p_i_codon[4]/=sumPi_Cod_second;
						p_i_codon[5]/=sumPi_Cod_second;
						p_i_codon[6]/=sumPi_Cod_second;
						p_i_codon[7]/=sumPi_Cod_second;
						}
					if (sumPi_Cod_third != 1.0) 
						{
						p_i_codon[8]/=sumPi_Cod_third;
						p_i_codon[9]/=sumPi_Cod_third;
						p_i_codon[10]/=sumPi_Cod_third;
						p_i_codon[11]/=sumPi_Cod_third;
						}
					}
				else
					{
					fprintf (stderr, "PARAMETER ERROR: Bad number of frequencies (4 or 12) (%d)\n\n", freqNumber);
					PrintUsage();
					}
			break;
			case 'T':
	        	titv = atof(argv[i]);
				if (titv < 0) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad ti/tv (%f)\n\n", titv);
					PrintUsage();
					}
				/*else
					{
					if (doCodonModel == NO)
						doHKY = YES;
					else
						{
						doCodon_HKY = YES;
						doCodon_GTR = NO;
						doCodon_NGTR = NO;
						}
					}*/
			break;
			case 'V':
	        	Rmat[0] = atof(argv[i]);
	        	Rmat[1] = atof(argv[++i]);
	        	Rmat[2] = atof(argv[++i]);
	        	Rmat[3] = atof(argv[++i]);
	        	Rmat[4] = atof(argv[++i]);
	        	Rmat[5] = atof(argv[++i]);
				if (Rmat[5]!=1.0) 
					{
					for (j=0; j<5; j++) 
						Rmat[j]/=Rmat[5];
					Rmat[5]=1.0;
					}
				/*doGTR = YES;
				doHKY = NO;
				doGTnR = NO;
				if (doCodonModel == YES)
					{
					doGTR = NO;
					doHKY = NO;
					doCodon_GTR = YES;
					doCodon_HKY = NO;
					doCodon_NGTR = NO;
					}*/
			break;
			case '@':							/*	AC CA AG GA AT TA CG GC CT TC GT=1 TG */
	        	NRmat[0] = atof(argv[i]);
	        	NRmat[1] = atof(argv[++i]);
	        	NRmat[2] = atof(argv[++i]);
	        	NRmat[3] = atof(argv[++i]);
	        	NRmat[4] = atof(argv[++i]);
	        	NRmat[5] = atof(argv[++i]);
	        	NRmat[6] = atof(argv[++i]);
	        	NRmat[7] = atof(argv[++i]);
	        	NRmat[8] = atof(argv[++i]);
	        	NRmat[9] = atof(argv[++i]);
	        	NRmat[10] = atof(argv[++i]); 
	        	NRmat[11] = atof(argv[++i]);
				if (NRmat[10]!=1.0) 
					{
					for (j=0; j<12 && j!=10; j++) 
						NRmat[j]/=NRmat[10];
					NRmat[10]=1.0;
					}
				/*doGTR = NO;
				doGTnR = YES;
				doHKY = NO;
				if (doCodonModel == YES)
					{
					doGTnR = NO;
					doGTR = NO;
					doCodon_GTR = NO;
					doCodon_HKY = NO;
					doCodon_NGTR = YES;
					}*/
			break;	
			case 'I':
	        	pinv = atof(argv[i]);
				if (pinv < 0) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad p-inv (%f)\n\n", pinv);
					PrintUsage();
					}
			break;
			case 'A':
	        	alpha = atof(argv[i]);
				if (alpha < 0) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad gamma shape (alpha) for rate variation among sites(%f)\n\n", alpha);
					PrintUsage();
					}
				doRateHet = YES;
			break;
			case 'N':
				argument = atof(argv[i]);
				numDataSets = (int) argument;
				if (numDataSets <1)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad number of replicates (%d)\n\n", numDataSets);
					PrintUsage();
					}
			break;
			case 'Y':
				argument = atof(argv[i]);
				noisy = (int) argument;
				if (noisy < 0)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad noisy value (%d)\n\n", noisy);
					PrintUsage();
					}
			break;
			case 'O':
	        	outgroupBranchLength = atof(argv[i]);
				if (outgroupBranchLength <= 0) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad outgroup branch lenght value (%f)\n\n", outgroupBranchLength);
					PrintUsage();
					}
				thereisOutgroup = YES;
			break;
			case 'G':
	        	growthRate = atof(argv[i]);
				/*if (growthRate < 0) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad growth rate (%f)\n\n", growthRate);
					PrintUsage();
					}*/
				/*if (growthRate != 0) 
        			{
        			doExponential = YES;
		        	if (doDemographics == YES)
						{
						fprintf (stderr, "PARAMETER ERROR: Cannot have both exponential (-g) and demographics (-p)\n\n");
						exit (1);
						}
					}*/
			break;
			case 'P':
				/*doDemographics = YES;
		       if (doExponential == YES)
					{
					fprintf (stderr, "PARAMETER ERROR: Cannot have both demographics periods (-g) and other demographics (-p)\n\n");
					exit (1);
					}*/
	        	argument = atof(argv[i]);
				numPeriods = (int) argument;
				if (numPeriods < 0) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad number of periods (%d)\n\n", numPeriods);
					PrintUsage();
					}
				/*if (numPeriods == 0)
					doDemographics = NO;*/
				if (numPeriods > 0)
					{
					Nbegin = 	(int *) calloc(numPeriods+1,(long) sizeof(int));
					Nend =		(int *) calloc(numPeriods+1,(long) sizeof(int));
					cumDuration =	(int *) calloc(numPeriods+1,(long) sizeof(int));
					periodGrowth =	(double *) calloc(numPeriods+1,(long) sizeof(double));
					if (Nbegin == NULL || Nend == NULL || cumDuration == NULL)
						{
						fprintf (stderr, "PARAMETER ERROR: Could not allocate demographic vectors (%lu bytes)\n", numPeriods *(long) sizeof(int));
						exit (1);
						}
					for (j=1; j<=numPeriods; j++)
						{
						argument = atof(argv[++i]);
						Nbegin[j] = (int) argument;
						argument = atof(argv[++i]);
						Nend[j] = (int) argument;
						argument = atof(argv[++i]);
						cumDuration[j] = (int) argument + cumDuration[j-1];
						}
					}
			break;
			case 'Q':
				/*doMigration = YES;*/
	        	sumPopul = 0;
				
				argument = atof(argv[i]);
				numPopulations = (int) argument;
				if (numPopulations < 1 || numPopulations > numSequences) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad number of subpopulations (%d)\n\n", numPopulations);
					PrintUsage();
					}
				
				initPopulation = (int *) calloc(numPopulations+1,(long) sizeof(int));
				if (!initPopulation)
					{
					fprintf (stderr, "PARAMETER ERROR: Could not allocate initPopulation of migration model (%lu bytes)\n", (numPopulations+1) *(long) sizeof(int));
					exit (1);
					}
				for (j=1; j<=numPopulations; j++)
					{
					argument = atof(argv[++i]);
					initPopulation[j] = (int) argument;
					if (initPopulation[j] > numSequences || initPopulation[j] <= 0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad initial population (%d)\n\n", initPopulation[j]);
						PrintUsage();
						}					
					sumPopul = sumPopul + initPopulation[j];
					}
				argument = atof(argv[++i]);
				migrationRate = (double) argument;
				if (migrationRate < 0)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad migration rate (%lf)\n\n", migrationRate);
					PrintUsage();
					}
				if (sumPopul != numSequences)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad initial subpopulation. The number of sequences is different that the sum of nodes in demes \n\n");
					PrintUsage();
					}
			break;
			case '%':
				argument = atof(argv[i]);
				numConvergDemes = (int) argument;
				if (numConvergDemes < 0 || numConvergDemes >= numSequences) 
					{
					fprintf(stderr, "PARAMETER ERROR: Bad number of convergencies demes events (%d) \n\n", numConvergDemes);
					PrintUsage();
					}
				/*fprintf(stderr, "\n\n numConvergDemes = %d \n", numConvergDemes);*/
				
				deme_a_old = (int *) calloc(numConvergDemes+1,(long) sizeof(int));
				if (!deme_a_old)
					{
					fprintf (stderr, "PARAMETER ERROR: Could not allocate deme_a of convergencies demes events (%lu bytes)\n", (numConvergDemes+1) *(long) sizeof(int));
					exit (1);
					}
				deme_b_old = (int *) calloc(numConvergDemes+1,(long) sizeof(int));
				if (!deme_b_old)
					{
					fprintf (stderr, "PARAMETER ERROR: Could not allocate deme_b of convergencies demes events (%lu bytes)\n", (numConvergDemes+1) *(long) sizeof(int));
					exit (1);
					}
				convDemTimes_old = (double*) calloc ((numConvergDemes+1), sizeof (double)); 
				if (convDemTimes_old == NULL)
					{
					fprintf (stderr, "PARAMETER ERROR: Could not allocate convDemTimes of convergencies demes events (%lu bytes)\n", (numConvergDemes+1) *(long) sizeof(double));
					exit (1);
					}
				if (numConvergDemes > 0)
					{
					for (j=1; j<=numConvergDemes; j++)
						{
						argument = atof(argv[++i]);
						deme_a_old[j] = (int) argument;
						argument = atof(argv[++i]);
						deme_b_old[j] = (int) argument;
						argument = atof(argv[++i]);
						convDemTimes_old[j] = (double) argument;
						}
					for (j=1; j<=numConvergDemes; j++)
						{
						/*fprintf(stderr, "\n\n deme_a_old[%d] = %d \n", j, deme_a_old[j]);
						fprintf(stderr, "\n\n deme_b_old[%d] = %d \n", j, deme_b_old[j]);
						fprintf(stderr, "\n\n convDemTimes_old[%d] = %lf \n", j, convDemTimes_old[j]);*/
						if (deme_a_old[j] == deme_b_old[j])
							{
							fprintf (stderr, "PARAMETER ERROR: Bad number of deme in convergencies demes events (%d)(%d), they must to be differents \n\n", deme_a_old[j], deme_b_old[j]);
							PrintUsage();
							}
						if (convDemTimes_old[j] <= 0)
							{
							fprintf (stderr, "PARAMETER ERROR: Bad time to convergencies demes events (%lf)\n\n", convDemTimes_old[j]);
							PrintUsage();
							}
						}
					k = h = 0;
					for (j=1; j<=numConvergDemes; j++)
						{
						k = convDemTimes_old[j];
						for (h=1; h<=numConvergDemes; h++)
							{
							if (j != h && convDemTimes_old[j] == convDemTimes_old[h])
								{
								fprintf (stderr, "PARAMETER ERROR: Bad time to convergencies demes events (%lf), it can not have two events at the same time\n\n", convDemTimes_old[h]);
								PrintUsage();
								}
							if (j != h && deme_a_old[j] == deme_a_old[h])
								{
								fprintf (stderr, "PARAMETER ERROR: Bad number of demes to convergencies demes events (%d)(%d), it can not have a same deme converging at two different times\n\n", deme_a_old[j], deme_a_old[h]);
								PrintUsage();
								}
							if (j != h && deme_b_old[j] == deme_b_old[h])
								{
								fprintf (stderr, "PARAMETER ERROR: Bad number of demes to convergencies demes events (%d)(%d), it can not have a same deme converging at two different times\n\n", deme_b_old[j], deme_b_old[h]);
								PrintUsage();
								}
							}
						}
					}
			break;
			case 'J':
				strcpy(treeFile, "trees");
				doPrintTrees = YES;
			break;
			case 'B':
				strcpy(alignmentFile, "sequences");
			break;
			case 'K':
				strcpy(timesFile, "times");
				doPrintTimes = YES;
			break;
			case 'D':
				strcpy(breakpointFile, "breakpoints");
				doPrintBreakpoints = YES;
			break;
			case 'X':
				strcpy(MRCAFile, "seqGMRCA"); 
				doMRCAFile = YES;
			break;
			case '*':
				argument = atof(argv[i]);
				formatNumber = (int) argument;
				if (formatNumber > 3 || formatNumber < 1) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad format of output alignments files chosen (1-3) (%d)\n\n", formatNumber);
					PrintUsage();
					}				
				if (formatNumber == 1) /* phylip */
					{
					doPrintFASTA = NO;
					doPrintNEXUS = NO;
					}
				if (formatNumber == 2) /* fasta */
					{
					doPrintFASTA = YES;
					doPrintNEXUS = NO;
					}
				if (formatNumber == 3) /* nexus */
					{
					doPrintFASTA = NO;
					doPrintNEXUS = YES;
					}
				if (formatNumber > 3 || formatNumber < 1)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad format of output alignments files chosen (1-3) (%d)\n\n", (int) argument);
					PrintUsage();
					}					
			break;
			case 'C':
				doPrintAncestralSequences = YES;
			break;
			case 'Z':
				doSeparatedSequences = YES;
			break;
			case '$':
				doOutMRCAfiles = YES;
			break;
			case '+':
				doPrintOmegasPerSitefiles = YES;
			break;
			case 'W':
				if (recombinationRate <= 0)
					{
					fprintf (stderr, "PARAMETER ERROR: Cannot fix number of recombination events when the recombination rate is 0\n\n");
					PrintUsage();
					}	
				argument = atof(argv[i]);
				fixedNumRecEvents = (int) argument;
				if (fixedNumRecEvents <= 0)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad fixed number of recombination events value (%d)\n\n", fixedNumRecEvents);
					PrintUsage();
					}
			doFixNumRecEvents = YES;
			break;

			case '?':
				PrintUsage();
			break;
			case 'H':
				PrintUsage();
			break;
			default :
				fprintf(stderr, "PARAMETER ERROR: Incorrect parameter: %c\n\n", flag);
				PrintUsage();
			break;
		}
	} 
}



/***************************** ReadParametersFromFile *******************************/
/* Reads parameter values from the parameter file */

void ReadParametersFromFile()
{
	int  	j, i, modelNumber, k, h, l;
	char 	ch;
	int		to, from;
	double	sumPi, sumPi_b, sumPi_Cod_first, sumPi_Cod_second, sumPi_Cod_third;
	int		sumPopul;
	float argument;
	
	/* Used: A B C D E F G H I J K L M N O P Q R S T U V W X Y Z ? # @ $ % * + = / _    next ideas : ;  */

	if(feof(stdin)){
		fprintf(stderr, "PARAMETER ERROR: Unable to read parameters from stdin\n");
		exit(0);
	}

	ch=fgetc(stdin);
	while(isspace(ch))
		ch=fgetc(stdin);
	while(ch=='[')
	{
		ReadUntil(stdin, ']', "closing bracket");
		ch=fgetc(stdin);
		while(isspace(ch))
			ch=fgetc(stdin);
	}

	while(!feof(stdin))
		{
		argument = 0;
		ch=toupper(ch);
		switch (ch) 
			{
			case '#':
				if (fscanf(stdin, "%lu bytes", &userSeed) !=1 || userSeed < 0) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad seed (%d)\n\n", (int) userSeed);
					PrintUsage();
					}
			break;
			case 'S':
				if (fscanf(stdin, "%d", &numSequences)!=1 || numSequences<1)
				 	{
					fprintf(stderr, "PARAMETER ERROR: Bad sample size (%d) \n\n", numSequences);
					PrintUsage();
					}
			break;
			case 'L':
				if (fscanf(stdin, "%d", &numSites)!=1 || numSites<1) 
					{
					fprintf(stderr, "PARAMETER ERROR: Bad sequence number of sites (%d)\n\n", numSites);
					PrintUsage();
					}
			break;
			case 'E':
				if (fscanf(stdin, "%d", &N)!=1 || N<1) 
					{
					fprintf(stderr, "PARAMETER ERROR: Bad effective population size (%d) \n\n", N);
					PrintUsage();
					}
			break;
			case '=':
				doDatedTips = YES;
				if (fscanf(stdin, "%f", &argument) !=1 || argument < 1) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad number of periods (%d)\n\n", (int) argument);
					PrintUsage();
					}
 				numTipDates = (int) argument;
				if (numTipDates <= 0) 
					{
					fprintf (stderr, "COMMAND-LINE PARAMETER ERROR: Bad number of sampling dates (%d)\n\n", numTipDates);
					PrintUsage();
					}
				
				datedSample = 	(SampleSt *) calloc(numTipDates, sizeof(SampleSt));
				if (datedSample == NULL)
					{
					fprintf (stderr, "COMMAND-LINE PARAMETER ERROR: Could not allocate sampling dates vectors (%ld)\n", numTipDates * (long) sizeof(SampleSt));
					exit (1);
					}
	
				for (j=0; j<numTipDates; j++)
					{
					fscanf(stdin, "%f", &argument);
	        		datedSample[j].time = (float) argument;
					fscanf(stdin, "%f", &argument);
					from = (int) argument;
					fscanf(stdin, "%f", &argument);
        			to = (int) argument;
	
					datedSample[j].member = 	(int *) calloc(to-from+1, sizeof(int));
					if (datedSample[j].member == NULL)
						{
						fprintf (stderr, "COMMAND-LINE PARAMETER ERROR: Could not allocate sampling dates vector[j=%d] (%ld)\n",j, to-from+1 * (long) sizeof(int));
						exit (1);
						}

					l = 0;
					for (k=from; k<=to; k++)
						datedSample[j].member[l++] = k;
					
					datedSample[j].size = l;
					}
															
				/* debugging */
				/*for (j=0; j<numTipDates; j++)
					{
					fprintf (stderr, "\nERR: time=%f size=%d sequences=", datedSample[j].time, datedSample[j].size);
					for (k=0; k<datedSample[j].size; k++)
						fprintf (stderr, " %d", datedSample[j].member[k]);	
					}
				*/
			break;
			case '/':
				if (fscanf(stdin, "%lf", &generationTime) !=1) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad generation time (%f)\n\n", generationTime);
					PrintUsage();
					}
			break;
			case 'R':
				if (fscanf(stdin, "%lf", &recombinationRate)!=1) 
				{
					fprintf(stderr, "PARAMETER ERROR: Bad recombination rate (%f) \n\n", recombinationRate);
					PrintUsage();
				}
			break;
			case 'U':
				if (fscanf(stdin, "%lf", &mutationRate)!=1) 
					{
					fprintf(stderr, "PARAMETER ERROR: Bad mutation rate (%f) \n\n", mutationRate);
					PrintUsage();
					}
			break;
			case '_':
				if (fscanf(stdin, "%f", &argument)!=1) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad haplid/diploid chosen (%d)\n\n", (int) argument);
					PrintUsage();
					}
				Nscaling = (int) argument;
				if (Nscaling < 1 || Nscaling > 2)
					{
					fprintf (stderr, "PARAMETER ERROR: Haploid/diplod option (1-2) (%d)\n\n", Nscaling);
					PrintUsage();
					}
			break;
			case 'M':
				/*doCodonModel = YES;*/
				sumPi = sumPi_b = 0.0;
				if (fscanf(stdin, "%f", &argument)!=1) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad model chosen (1-7) (%d)\n\n", (int) argument);
					PrintUsage();
					}
				modelNumber = (int) argument;
				/*fscanf(stdin, "%f", &modelNumber);
				if (modelNumber < 0.0)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad modelNumber (%3.2f)\n\n", modelNumber);
					PrintUsage();
					}*/
				
				if (modelNumber == 1)		/** omega constant (M0) **/
					{
					doOmegaCat = NO;
					doOmegaRateHetCont = NO;
					doOmegaRateHetDisc = NO;
					doM0 = YES;	
					
					fscanf(stdin, "%lf", &omega);
					if (omega < 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad omega (%3.2f)\n\n", omega);
						PrintUsage();
						}
					OmegaInit = omega;
					}
				else if (modelNumber == 2) /** omega by categories from the user **/
					{
					doOmegaCat = YES;
					doOmegaProb = YES;
					doOmegaRateHetCont = NO;
					doOmegaRateHetDisc = NO;
					
					if (fscanf(stdin, "%f", &argument) !=1 || argument > 10) /* max 10 categories */
						{
						fprintf (stderr, "PARAMETER ERROR: Bad number of omega categories (%d)\n\n", (int) argument);
						PrintUsage();
						}
					numOmegaCat = (int) argument;
					omegaVal = (double *) calloc((numOmegaCat+1),(long) sizeof(double));
					if (!omegaVal)
						{
						fprintf (stderr, "PARAMETER ERROR: Could not allocate omega values of categories (%lu bytes)\n", numOmegaCat *(long) sizeof(double));
						exit (1);
						}
					omegaProb = (double *) calloc((numOmegaCat+1),(long) sizeof(double));
					if (!omegaProb)
						{
						fprintf (stderr, "PARAMETER ERROR: Could not allocate omega probabilities of categories1 (%lu bytes)\n", numOmegaCat *(long) sizeof(double));
						exit (1);
						}
					for (j=1; j<=numOmegaCat; j++)
						{
						fscanf(stdin, "%lf", &omegaVal[j]);
						fscanf(stdin, "%lf", &omegaProb[j]);

						if (omegaProb[j] > 1 || omegaProb[j] < 0)
							{
							fprintf (stderr, "PARAMETER ERROR: Bad number of probabilities of omega categories2 (%3.2f)\n\n", omegaProb[j]);
							PrintUsage();
							}
						sumPi = sumPi + omegaProb[j];
						}
					if (sumPi != 1) /* update probabilities of categories */
						{
						if (numOmegaCat == 1)
							omegaProb[1] = 1.0;
						else
							{
							for (j=1; j<=numOmegaCat; j++)
								{
								omegaProb[j]/=sumPi;
								sumPi_b = sumPi_b + omegaProb[j];
								}
							if ((int)sumPi_b != 1)
								{
								fprintf (stderr, "\n ERROR in the sum of probabilities of omega categories");
								exit (-1);
								}
							}
						}
					}
				else if (modelNumber == 3) /** omega by discrete heterogeneous rate **/
					{
					doOmegaCat = NO;
					doOmegaRateHetCont = NO;
					doOmegaRateHetDisc = YES;
					
					fscanf(stdin, "%d", &numOmegaCat);
					if (numOmegaCat < 0 || numOmegaCat > 10)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad number of omega categories (%d)\n\n", numOmegaCat);
						PrintUsage();
						}
					fscanf(stdin, "%lf", &OmegaRateHet);
					if (OmegaRateHet <= 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad omega heterogeneous rate (%3.2f)\n\n", OmegaRateHet);
						PrintUsage();
						}
					fscanf(stdin, "%lf", &omega);
					if (omega < 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad omega (%3.2f)\n\n", omega);
						PrintUsage();
						}
					OmegaInit = omega;
					/*fprintf (stderr, "\n\n %d %lf %lf \n",numOmegaCat, OmegaRateHet, omega);*/
					}
				else if (modelNumber == 4) /** omega by continuous heterogeneous rate **/
					{
					/*fprintf (stderr, "\n lee file model M1 \n\n");*/
					doOmegaCat = NO;
					doOmegaRateHetCont = YES;
					doOmegaRateHetDisc = NO;
					
					fscanf(stdin, "%lf", &OmegaRateHet);
					if (OmegaRateHet <= 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad omega heterogeneous rate (%3.2f)\n\n", OmegaRateHet);
						PrintUsage();
						}
					fscanf(stdin, "%lf", &omega);
					if (omega < 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad omega (%3.2f)\n\n", omega);
						PrintUsage();
						}
					OmegaInit = omega;
					}
				else if (modelNumber == 5) /** M1 **/
					{
					doM1 = YES;	
					/*doOmegaCat = NO;*/
					doOmegaRateHetCont = NO;
					doOmegaRateHetDisc = NO;
				
					fscanf(stdin, "%lf", &M1_P0_omeg0);
					if (M1_P0_omeg0 < 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad P0 in model M1 (%3.2f)\n\n", M1_P0_omeg0);
						PrintUsage();
						}
					fscanf(stdin, "%lf", &M1_omega0);
					if (M1_omega0 < 0.0 || M1_omega0 > 1)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad omega for P0 in model M1 (%3.2f)\n\n", M1_omega0);
						PrintUsage();
						}

					M1_P1_omeg1 = 1 - M1_P0_omeg0;
					doOmegaCat = YES;
					numOmegaCat = 2;
					
					omegaVal = (double *) calloc((numOmegaCat+1),(long) sizeof(double));
					if (!omegaVal)
						{
						fprintf (stderr, "PARAMETER ERROR: Could not allocate omega values of categories (%lu bytes)\n", numOmegaCat *(long) sizeof(double));
						exit (1);
						}
					omegaProb = (double *) calloc((numOmegaCat+1),(long) sizeof(double));
					if (!omegaProb)
						{
						fprintf (stderr, "PARAMETER ERROR: Could not allocate omega probabilities of categories1 (%lu bytes)\n", numOmegaCat *(long) sizeof(double));
						exit (1);
						}

					omegaProb[1] = M1_P0_omeg0;
					omegaProb[2] = M1_P1_omeg1;
					omegaVal[1] = M1_omega0;
					omegaVal[2] = 1.0;

					}
				else if (modelNumber == 6) /** M7 **/
					{
					doM7 = YES;

					doOmegaCat = NO;
					doOmegaRateHetCont = NO;
					doOmegaRateHetDisc = NO;
				
					fscanf(stdin, "%lf", &M7_p_beta);
					if (M7_p_beta < 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad p in model M7 (%3.2f)\n\n", M7_p_beta);
						PrintUsage();
						}
					
					fscanf(stdin, "%lf", &M7_q_beta);
					if (M7_q_beta < 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad q in model M8 (%3.2f)\n\n", M7_q_beta);
						PrintUsage();
						}

					}
				else if (modelNumber == 7) /** M8 **/
					{
					doM8 = YES;

					doOmegaCat = NO;
					doOmegaRateHetCont = NO;
					doOmegaRateHetDisc = NO;

					fscanf(stdin, "%lf", &M8_P0_beta);
					if (M8_P0_beta < 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad P0 in model M8 (%3.2f)\n\n", M8_P0_beta);
						PrintUsage();
						}
					M8_P1_omega = 1 - M8_P0_beta;
					
					fscanf(stdin, "%lf", &M8_p_beta);
					if (M8_p_beta < 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad p in model M8 (%3.2f)\n\n", M8_p_beta);
						PrintUsage();
						}
					
					fscanf(stdin, "%lf", &M8_q_beta);
					if (M8_q_beta < 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad q in model M8 (%3.2f)\n\n", M8_q_beta);
						PrintUsage();
						}
					
					fscanf(stdin, "%lf", &M8_omegaP1);
					if (M8_omegaP1 < 0.0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad omega for P1 in model M8 (%3.2f)\n\n", M8_omegaP1);
						PrintUsage();
						}
					}
				else
					{
					fprintf (stderr, "PARAMETER ERROR: Bad model chosen (1-7) (%d) fromFile\n\n", modelNumber);
					PrintUsage();
					}
				sumPi = sumPi_b = 0.0;
			break;
			case 'F':
				if (fscanf(stdin, "%f", &argument) !=1 || argument > 12) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad number of frequencies (4 or 12) (%d)\n\n", (int) argument);
					PrintUsage();
					}
				freqNumber = (int) argument;
				
				if (freqNumber == 4)
					{
					if (fscanf(stdin, "%lf %lf %lf %lf", &p_i[0], &p_i[1], &p_i[2], &p_i[3])!=4) 
						{
						fprintf(stderr, "PARAMETER ERROR: Bad Base Frequencies (put only 4 base frequencies) \n\n");
						PrintUsage();
						}
					
					equalBaseFreq = YES;
					for (i = 1; i < 4; i++)		
						if (p_i[i] != p_i[i-1])
							{
							equalBaseFreq = NO;
							break;
							}
							
					sumPi = p_i[0] + p_i[1] + p_i[2] + p_i[3];
					if (sumPi !=1.0) 
						{
						p_i[0]/=sumPi;
						p_i[1]/=sumPi;
						p_i[2]/=sumPi;
						p_i[3]/=sumPi;
						}
					p_i_codon[0]=p_i[0];
					p_i_codon[1]=p_i[1];
					p_i_codon[2]=p_i[2];
					p_i_codon[3]=p_i[3];
					p_i_codon[4]=p_i[0];
					p_i_codon[5]=p_i[1];
					p_i_codon[6]=p_i[2];
					p_i_codon[7]=p_i[3];
					p_i_codon[8]=p_i[0];
					p_i_codon[9]=p_i[1];
					p_i_codon[10]=p_i[2];
					p_i_codon[11]=p_i[3];					
					}
				else if (freqNumber == 12)
					{
					fscanf(stdin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &p_i_codon[0], &p_i_codon[1], &p_i_codon[2], &p_i_codon[3], 
					&p_i_codon[4], &p_i_codon[5], &p_i_codon[6], &p_i_codon[7], &p_i_codon[8], &p_i_codon[9], &p_i_codon[10], &p_i_codon[11]);
						
					for (i = 0; i < 12; i++)
						if (p_i_codon[i] < 0 || p_i_codon[i] > 1)
							{
							fprintf (stderr, "PARAMETER ERROR: Bad number of enter frequencies (it must to be between 0 and 1) (%lf)\n\n", p_i_codon[i]);
							PrintUsage();
							}
					equalBaseFreqCod = YES;
					for (i = 1; i < 12; i++)		
						if (p_i_codon[i] != p_i_codon[i-1])
							{
							equalBaseFreqCod = NO;
							break;
							}		
					
					sumPi_Cod_first = p_i_codon[0] + p_i_codon[1] + p_i_codon[2] + p_i_codon[3];
					sumPi_Cod_second = p_i_codon[4] + p_i_codon[5] + p_i_codon[6] + p_i_codon[7];
					sumPi_Cod_third = p_i_codon[8] + p_i_codon[9] + p_i_codon[10] + p_i_codon[11];
					if (sumPi_Cod_first != 1.0) 
						{
						p_i_codon[0]/=sumPi_Cod_first;
						p_i_codon[1]/=sumPi_Cod_first;
						p_i_codon[2]/=sumPi_Cod_first;
						p_i_codon[3]/=sumPi_Cod_first;
						}
					if (sumPi_Cod_second != 1.0) 
						{
						p_i_codon[4]/=sumPi_Cod_second;
						p_i_codon[5]/=sumPi_Cod_second;
						p_i_codon[6]/=sumPi_Cod_second;
						p_i_codon[7]/=sumPi_Cod_second;
						}
					if (sumPi_Cod_third != 1.0) 
						{
						p_i_codon[8]/=sumPi_Cod_third;
						p_i_codon[9]/=sumPi_Cod_third;
						p_i_codon[10]/=sumPi_Cod_third;
						p_i_codon[11]/=sumPi_Cod_third;
						}
					}
				else
					{
					fprintf (stderr, "PARAMETER ERROR: Bad number of frequencies (4 or 12) (%d)\n\n", freqNumber);
					PrintUsage();
					}
			break;
			case 'T':
				if (fscanf(stdin, "%lf", &titv)!=1) 
					{
					fprintf(stderr, "PARAMETER ERROR: Bad ti/tv (%f) \n\n", titv);
					PrintUsage();
					}
				/*else
					{
					if (doCodonModel == NO)
						doHKY = YES;
					else
						{
						doCodon_HKY = YES;
						doCodon_GTR = NO;
						doCodon_NGTR = NO;
						}
					}*/
			break;
			case 'V':
				if (fscanf(stdin, "%lf %lf %lf %lf %lf %lf", &Rmat[0], &Rmat[1], 
												&Rmat[2], &Rmat[3], &Rmat[4], &Rmat[5])!=6) 
					{
					fprintf(stderr, "PARAMETER ERROR: Bad general rate matrix (-rx x x x x x)\n\n");
					PrintUsage();
					}
				if (Rmat[5]!=1.0) 
					{
					for (j=0; j<5; j++) 
						Rmat[j]/=Rmat[5];
					Rmat[5]=1.0;
					}
				/*doGTR = YES;
				doHKY = NO;
				doGTnR = NO;
				if (doCodonModel == YES)
					{
					doGTnR = NO;
					doHKY = NO;
					doGTR = NO;
					doCodon_GTR = YES;
					doCodon_HKY = NO;
					doCodon_NGTR = NO;
					}*/
			break;
			case '@':
				if (fscanf(stdin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &NRmat[0], &NRmat[1], &NRmat[2], &NRmat[3], &NRmat[4], &NRmat[5],
									&NRmat[6], &NRmat[7], &NRmat[8], &NRmat[9], &NRmat[10], &NRmat[11])!=12) 
					{
					fprintf(stderr, "Bad general rate matrix (-rx x x x x x x x x x x) (AC CA AG GA AT TA CG GC CT TC GT=1 TG)\n\n");
					PrintUsage();
					}
				if (NRmat[10]!=1.0) 
					{
					for (j=0; j<12 && j!=10; j++) 
						NRmat[j]/=NRmat[10];
					NRmat[10]=1.0;
					}
				/*for (j=0; j<12; j++)
					printf ("\nR-matrix                            =  %3.2f", NRmat[j]);*/
				
				/*doGTR = NO;
				doGTnR = YES;
				doHKY = NO;
				if (doCodonModel == YES)
					{
					doGTnR = NO;
					doHKY = NO;
					doGTR = NO;
					doCodon_GTR = NO;
					doCodon_HKY = NO;
					doCodon_NGTR = YES;
					}*/
			break;	
			case 'I':
				if (fscanf(stdin, "%lf", &pinv)!=1) 
					{
					fprintf(stderr, "PARAMETER ERROR: Bad p-inv (%f) \n\n", pinv);
					PrintUsage();
					}
			break;
			case 'A':
				if (fscanf(stdin, "%lf", &alpha)!=1 || alpha<=0.0) 
					{
					fprintf(stderr, "PARAMETER ERROR: Bad Gamma Shape (%f) \n", alpha);
					exit(0);
					}
				else
					doRateHet = YES;
			break;
			case 'C':
				doPrintAncestralSequences = YES;
			break;
			case '*':
				if (fscanf(stdin, "%f", &argument) !=1 || argument > 10) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad format of output alignments files chosen (1-3) (%d)\n\n", (int) argument);
					PrintUsage();
					}
				formatNumber = (int) argument;
				
				if (formatNumber == 1) /* phylip */
					{
					doPrintFASTA = NO;
					doPrintNEXUS = NO;
					}
				if (formatNumber == 2) /* fasta */
					{
					doPrintFASTA = YES;
					doPrintNEXUS = NO;
					}
				if (formatNumber == 3) /* nexus */
					{
					doPrintFASTA = NO;
					doPrintNEXUS = YES;
					}
				if (formatNumber > 3 || formatNumber < 1)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad format of output alignments files chosen (1-3) (%d)\n\n", (int) argument);
					PrintUsage();
					}					
			break;

			case 'Z':
				doSeparatedSequences = YES;
			break;
			case '$':
				doOutMRCAfiles = YES;
			break;
			case '+':
				doPrintOmegasPerSitefiles = YES;
			break;
			case 'N':
				if (fscanf(stdin, "%d", &numDataSets)!=1 || numDataSets <1)
					{
					fprintf(stderr, "PARAMETER ERROR: Bad number of replicates (%d)\n\n",numDataSets);
					PrintUsage();
					}
			break;
			case 'Y':
				if (fscanf(stdin, "%d", &noisy)!=1 || noisy <0)
					{
					fprintf(stderr, "PARAMETER ERROR: Bad noisy value (%d)\n\n", noisy);
					PrintUsage();
					}
			break;
			case 'W':
				if (recombinationRate <= 0)
					{
					fprintf (stderr, "PARAMETER ERROR: Cannot fix number of recombination events when the recombination rate is 0\n\n");
					PrintUsage();
					}	
				if (fscanf(stdin, "%d", &fixedNumRecEvents) !=1 || fixedNumRecEvents <= 0)
					{
					fprintf(stderr, "PARAMETER ERROR: Bad fixed number of recombination events (%d)\n\n", fixedNumRecEvents);
					PrintUsage();
					}
				doFixNumRecEvents = YES;
			break;
			case 'O':
				if (fscanf(stdin, "%lf", &outgroupBranchLength) <= 0)
					{
					fprintf(stderr, "PARAMETER ERROR: Bad outgroup branch lenght value (%f)\n\n", outgroupBranchLength);
					PrintUsage();
					}
				else
					thereisOutgroup = YES;
			break;
			case 'G':
				if (fscanf(stdin, "%lf", &growthRate) !=1) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad growth rate (%f)\n\n", growthRate);
					PrintUsage();
					}
				/*if (growthRate != 0) 
	    			{
	    			doExponential = YES;
		        	if (doDemographics == YES)
						{
						fprintf (stderr, "PARAMETER ERROR: Cannot have both exponential (-g) and other demographics(-p)\n\n");
						exit (1);
						}
					}*/
			break;
			case 'P':
				/*doDemographics = YES;
	        	if (doExponential == YES)
					{
					fprintf (stderr, "PARAMETER ERROR: Cannot have both demographics periods (-p) and other demographics (-g)\n\n");
					exit (1);
					}*/
				if (fscanf(stdin, "%f", &argument) !=1 || argument < 0) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad number of periods (%d)\n\n", (int) argument);
					PrintUsage();
					}
				numPeriods = (int) argument;
				/*if (numPeriods == 0)
					doDemographics = NO;*/
				if (numPeriods > 0)
					{
					Nbegin =	(int *) calloc(numPeriods+1,(long) sizeof(int));
					Nend = 		(int *) calloc(numPeriods+1,(long) sizeof(int));
					cumDuration =	(int *) calloc(numPeriods+1,(long) sizeof(int));
					periodGrowth =	(double *) calloc(numPeriods+1,(long) sizeof(double));
					if (Nbegin == NULL || Nend == NULL || cumDuration == NULL)
						{
						fprintf (stderr, "PARAMETER ERROR: Could not allocate demographic vectors (%lu bytes)\n", numPeriods *(long) sizeof(int));
						exit (1);
						}
					for (j=1; j<=numPeriods; j++)
						{
						fscanf(stdin, "%f", &argument);
						Nbegin[j] = (int) argument;
						fscanf(stdin, "%f", &argument);
						Nend[j] = (int) argument;
						fscanf(stdin, "%f", &argument);
						cumDuration[j] = (int) argument + cumDuration[j-1];
						}
					}
			break;
			case 'Q':
				/*doMigration = YES;*/
	        	sumPopul = 0;
				
				if (fscanf(stdin, "%f", &argument) !=1 || argument < 1 || argument > numSequences) 
					{
					fprintf (stderr, "PARAMETER ERROR: Bad number of subpopulations (%d)\n\n", (int) argument);
					PrintUsage();
					}
				numPopulations = (int) argument;
				
				initPopulation = (int *) calloc(numPopulations+1,(long) sizeof(int));
				if (!initPopulation)
					{
					fprintf (stderr, "PARAMETER ERROR: Could not allocate initPopulation of migration model (%lu bytes)\n", (numPopulations+1) *(long) sizeof(int));
					exit (1);
					}
				for (j=1; j<=numPopulations; j++)
					{
					fscanf(stdin, "%d", &initPopulation[j]);
					if (initPopulation[j] > numSequences || initPopulation[j] <= 0)
						{
						fprintf (stderr, "PARAMETER ERROR: Bad initial subpopulation (%d)\n\n", initPopulation[j]);
						PrintUsage();
						}
					sumPopul = sumPopul + initPopulation[j];
					}
				fscanf(stdin, "%lf", &migrationRate);
				if (migrationRate < 0)
					{
					fprintf (stderr, "PARAMETER ERROR: Bad migration rate (%lf)\n\n", migrationRate);
					PrintUsage();
					}
				if (sumPopul != numSequences)
					{
					fprintf (stderr, "\n PARAMETER ERROR: Bad initial population. The sample size (%d) is different to the sum of nodes in demes (%d) \n\n", numSequences, sumPopul);
					PrintUsage();
					}
			break;
			case '%':
				if (fscanf(stdin, "%f", &argument) !=1 || argument < 1 || argument >= numSequences) 
					{
					fprintf(stderr, "PARAMETER ERROR: Bad number of convergencies demes events (%d) \n\n", numConvergDemes);
					PrintUsage();
					}
				numConvergDemes = (int) argument;
				/*fprintf(stderr, "\n\n numConvergDemes = %d \n", numConvergDemes);*/
				
				deme_a_old = (int *) calloc(numConvergDemes+1,(long) sizeof(int));
				if (!deme_a_old)
					{
					fprintf (stderr, "PARAMETER ERROR: Could not allocate deme_a of convergencies demes events (%lu bytes)\n", (numConvergDemes+1) *(long) sizeof(int));
					exit (1);
					}
				deme_b_old = (int *) calloc(numConvergDemes+1,(long) sizeof(int));
				if (!deme_b_old)
					{
					fprintf (stderr, "PARAMETER ERROR: Could not allocate deme_b of convergencies demes events (%lu bytes)\n", (numConvergDemes+1) *(long) sizeof(int));
					exit (1);
					}
				convDemTimes_old = (double*) calloc ((numConvergDemes+1), sizeof (double)); 
				if (convDemTimes_old == NULL)
					{
					fprintf (stderr, "PARAMETER ERROR: Could not allocate convDemTimes of convergencies demes events (%lu bytes)\n", (numConvergDemes+1) *(long) sizeof(double));
					exit (1);
					}
				if (numConvergDemes > 0)
					{
					for (j=1; j<=numConvergDemes; j++)
						{
						fscanf(stdin, "%f", &argument);
						deme_a_old[j] = (int) argument;
						fscanf(stdin, "%f", &argument);
						deme_b_old[j] = (int) argument;
						fscanf(stdin, "%f", &argument);
						convDemTimes_old[j] = (double) argument;
						}
					for (j=1; j<=numConvergDemes; j++)
						{
						/*fprintf(stderr, "\n\n deme_a_old[%d] = %d \n", j, deme_a_old[j]);
						fprintf(stderr, "\n\n deme_b_old[%d] = %d \n", j, deme_b_old[j]);
						fprintf(stderr, "\n\n convDemTimes_old[%d] = %lf \n", j, convDemTimes_old[j]);*/
						if (deme_a_old[j] == deme_b_old[j])
							{
							fprintf (stderr, "PARAMETER ERROR: Bad number of deme in convergencies demes events (%d)(%d), they must to be differents\n\n", deme_a_old[j], deme_b_old[j]);
							PrintUsage();
							}
						if (convDemTimes_old[j] <= 0)
							{
							fprintf (stderr, "PARAMETER ERROR: Bad time to convergencies demes events (%lf)\n\n", convDemTimes_old[j]);
							PrintUsage();
							}
						}
					k = h = 0;
					for (j=1; j<=numConvergDemes; j++)
						{
						k = convDemTimes_old[j];
						for (h=1; h<=numConvergDemes; h++)
							{
							if (j != h && convDemTimes_old[j] == convDemTimes_old[h])
								{
								fprintf (stderr, "PARAMETER ERROR: Bad time to convergencies demes events (%lf), it can not have two events at the same time\n\n", convDemTimes_old[h]);
								PrintUsage();
								}
							if (j != h && deme_a_old[j] == deme_a_old[h])
								{
								fprintf (stderr, "PARAMETER ERROR: Bad number of demes to convergencies demes events (%d)(%d), it can not have a same deme converging at two different times\n\n", deme_a_old[j], deme_a_old[h]);
								PrintUsage();
								}
							if (j != h && deme_b_old[j] == deme_b_old[h])
								{
								fprintf (stderr, "PARAMETER ERROR: Bad number of demes to convergencies demes events (%d)(%d), it can not have a same deme converging at two different times\n\n", deme_b_old[j], deme_b_old[h]);
								PrintUsage();
								}
							}
						}
				}
			break;
			case 'B':
				ch=fgetc(stdin);
				if(isspace(ch))
				{
					strcpy(alignmentFile, "sequences");
				}
				else
					{
					j=0;
					do
						{
						alignmentFile[j]=ch;
						j++;
						ch=fgetc(stdin);
						}
					while(!isspace(ch));
					alignmentFile[j]='\0';
					}
			break;
			case 'J':
				ch=fgetc(stdin);
				if(isspace(ch))
				{
					strcpy(treeFile, "trees");
				}
				else
					{
					j=0;
					do
						{
						treeFile[j]=ch;
						j++;
						ch=fgetc(stdin);
						}
					while(!isspace(ch));
					treeFile[j]='\0';
					}
			doPrintTrees = YES;
			break;
			case 'K':
				ch=fgetc(stdin);
				if(isspace(ch))
				{
					strcpy(timesFile, "times");
				}
				else
					{
					j=0;
					do
						{
						timesFile[j]=ch;
						j++;
						ch=fgetc(stdin);
						}
					while(!isspace(ch));
					timesFile[j]='\0';
					}
			doPrintTimes = YES;
			break;
			case 'D':
				ch=fgetc(stdin);
				if(isspace(ch))
					{
					strcpy(breakpointFile, "breakpoints");
					}
				else
					{
					j=0;
					do
						{
						breakpointFile[j]=ch;
						j++;
						ch=fgetc(stdin);
						}
					while(!isspace(ch));
					breakpointFile[j]='\0';
					}
			doPrintBreakpoints = YES;
			break;
			case 'X':
				ch=fgetc(stdin);
				if(isspace(ch))
					{
					strcpy(MRCAFile, "MRCA");
					}
				else
					{
					j=0;
					do
						{
						MRCAFile[j]=ch;
						j++;
						ch=fgetc(stdin);
						}
					while(!isspace(ch));
					MRCAFile[j]='\0';
					}
			doMRCAFile = YES;
			break;
			case '?':
				PrintUsage();
			break;
			case 'H':
				PrintUsage();
			break;
			default :
				fprintf(stderr, "PARAMETER ERROR: Incorrect parameter: %c\n\n", ch);
				PrintUsage();
			break;
		}
		ch=fgetc(stdin);
		while(isspace(ch) && !feof(stdin))
			ch=fgetc(stdin);
		while(ch=='[')
			{
			ReadUntil(stdin, ']', "closing bracket");
			ch=fgetc(stdin);
			while(isspace(ch))
				ch=fgetc(stdin);
			}
		}
}





#ifdef USER_INPUT
/***************************** UserInput *******************************/
/* Prompts the user for parameter values */

void UserInput (long int *seed)
{
	int i, j, k, n, modelNumber;
	double sumPi, sumPi_b, sumPi_Cod_first, sumPi_Cod_second, sumPi_Cod_third;
	int sumPopul;
	i = j = k = n = 0;
	
	printf ("\n << USER INPUT >> ");
	printf ("\n\n Number of data sets      = ");
	scanf ("%d", &numDataSets);
	printf (" Number of sequences      = ");
	scanf ("%d", &numSequences);
	printf (" Number of sites (bp or codons)     = ");
	scanf ("%d", &numSites);
	printf (" Mutation rate            = ");
	scanf ("%lf", &mutationRate);
	printf (" Effective population size = ");
	scanf ("%d", &N);
	printf (" Recombination rate       = ");
	scanf ("%lf", &recombinationRate);
	
	printf (" Number fixed of recombination events       = ");
	scanf ("%d", &fixedNumRecEvents);
	
	printf (" How many frequencies value?, 4 or 12 (12 only codon model)       = ");
	scanf ("%d", &freqNumber);
	sumPi = sumPi_b = 0.0;
	if (freqNumber == 4)
		{
		printf (" Enter frequencies       = ");
		for (i = 0;i < 4;i++)
			scanf ("%lf", &p_i[i]);	
		equalBaseFreq = YES;
			for (i = 1; i < 4; i++)		
				if (p_i[i] != p_i[i-1])
					{
					equalBaseFreq = NO;
					break;
					}
							
		sumPi = p_i[0] + p_i[1] + p_i[2] + p_i[3];
		if (sumPi !=1.0) 
			{
			p_i[0]/=sumPi;
			p_i[1]/=sumPi;
			p_i[2]/=sumPi;
			p_i[3]/=sumPi;
			}
		p_i_codon[0]=p_i[0];
		p_i_codon[1]=p_i[1];
		p_i_codon[2]=p_i[2];
		p_i_codon[3]=p_i[3];
		p_i_codon[4]=p_i[0];
		p_i_codon[5]=p_i[1];
		p_i_codon[6]=p_i[2];
		p_i_codon[7]=p_i[3];
		p_i_codon[8]=p_i[0];
		p_i_codon[9]=p_i[1];
		p_i_codon[10]=p_i[2];
		p_i_codon[11]=p_i[3];					
		}
	else if (freqNumber == 12)
		{
		printf (" Enter frequencies       = ");
		for (i = 0;i < 12;i++)
			scanf ("%lf", &p_i_codon[i]);
		
		for (i = 0; i < 12; i++)
			if (p_i_codon[i] < 0 || p_i_codon[i] > 1)
				{
				fprintf (stderr, "PARAMETER ERROR: Bad number of enter frequencies (it must to be between 0 and 1) (%lf)\n\n", p_i_codon[i]);
				PrintUsage();
				}
	
		equalBaseFreqCod = YES;
		for (i = 1; i < 12; i++)		
			if (p_i_codon[i] != p_i_codon[i-1])
				{
				equalBaseFreqCod = NO;
				break;
				}		
					
		sumPi_Cod_first = p_i_codon[0] + p_i_codon[1] + p_i_codon[2] + p_i_codon[3];
		sumPi_Cod_second = p_i_codon[4] + p_i_codon[5] + p_i_codon[6] + p_i_codon[7];
		sumPi_Cod_third = p_i_codon[8] + p_i_codon[9] + p_i_codon[10] + p_i_codon[11];
		if (sumPi_Cod_first != 1.0) 
			{
			p_i_codon[0]/=sumPi_Cod_first;
			p_i_codon[1]/=sumPi_Cod_first;
			p_i_codon[2]/=sumPi_Cod_first;
			p_i_codon[3]/=sumPi_Cod_first;
			}
		if (sumPi_Cod_second != 1.0) 
			{
			p_i_codon[4]/=sumPi_Cod_second;
			p_i_codon[5]/=sumPi_Cod_second;
			p_i_codon[6]/=sumPi_Cod_second;
			p_i_codon[7]/=sumPi_Cod_second;
			}
		if (sumPi_Cod_third != 1.0) 
			{
			p_i_codon[8]/=sumPi_Cod_third;
			p_i_codon[9]/=sumPi_Cod_third;
			p_i_codon[10]/=sumPi_Cod_third;
			p_i_codon[11]/=sumPi_Cod_third;
			}
		}
	else
		{
		fprintf (stderr, "PARAMETER ERROR: Bad number of frequencies (4 or 12) (%d)\n\n", freqNumber);
		PrintUsage();
		}

	printf (" titv              = ");
	scanf ("%lf", &titv);
		
   printf (" Do you want rate heterogeneity among sites (0 = no; 1 = yes)? ");
   scanf ("%d", &doRateHet);
   if (doRateHet == YES)
           {
           printf (" Enter a shape parameter for the discrete gamma distribution: ");
           scanf ("%lf", &alpha);
           }
   else
           {
           alpha = 1.0;
           }
	printf (" Random number seed = ");
	scanf ("%lu bytes", &userSeed);
	
	printf (" Substitution Model:\n");
	printf ("  - If Codon Model enter 1\n");
	printf ("  - If Nucleotide Model enter 2\n");
	scanf ("%d", &i);
	while (i != 1 || i != 2)
		{
		printf (" You have to choose an option: Enter 1 to codon Model or 2 to Nucleotide Model \n");
		scanf ("%d", &i);
		}
		
	if (i == 2)
		{
		printf ("\n Nucleotide Model");
		doCodonModel = NO;
		
		printf ("\n Nucleotide Model can to be from types: HKY, GTR or GTR not reversible. If you want GTR Model enter number 1, GTR not reversible enter 2, and for HKY Model enter another number");
		j = 0;
		scanf ("%d", &j);
		if (j == 1)
			{
			doGTR = YES;
			doHKY = NO;
			printf ("\n You have to introduce the relative substitution rates to changes (rate matrix A-C A-G A-T C-G C-T G-T, for GTR models)(e.g. -v1 3 0.2 0.5 5 1)");
			scanf ("%lf", &Rmat[0]);			
			scanf ("%lf", &Rmat[1]);
			scanf ("%lf", &Rmat[2]);
			scanf ("%lf", &Rmat[3]);
			scanf ("%lf", &Rmat[4]);
			scanf ("%lf", &Rmat[5]);
			if (Rmat[5]!=1.0) 
				{
				for (j=0; j<5; j++) 
					Rmat[j]/=Rmat[5];
				Rmat[5]=1.0;
				}
			}
		else if (j == 2)
			{
			doHKY = NO;
			doGTnR = YES;
			printf("\n Introduce the vale rates for changes: AC CA AG GA AT TA CG GC CT TC GT");
			for (n=0;n<12;n++)
				scanf ("%lf", &NRmat[n]);
			if (NRmat[10]!=1.0) 
				{
				for (j=0; j<12 && j!=10; j++) 
					NRmat[j]/=NRmat[10];
				NRmat[10]=1.0;
				}
			}
		else
			{
			doHKY = YES;
			doGTR = NO;
			}
		}
	if (i == 1)
		{
		printf ("\n Codon Model");
		doCodonModel = YES;
		doGTR = NO;
		doHKY = NO;
		
		printf ("\n Omega values: You have to Introduce 1 (omega is constant), or 2 (omega categories), or 3 (discrete hetereogeneous rate), or 4 (continuous hetereogeneous rate)");
		scanf ("%d", &modelNumber);
		if (modelNumber > 4 || modelNumber < 1)
			{
			printf("\n Bad number of type of omega model");
			exit(-1);
			}
		
		if (modelNumber == 1) /* omega constant */
			{
			doOmegaCat = NO;
			doOmegaRateHetCont = NO;
			doOmegaRateHetDisc = NO;
			
			printf ("\n Introduce the omega value");
			scanf("%lf", &omega);
			OmegaInit = omega;
			}
		
		if (modelNumber == 2) /* omega by categories */
			{
			doOmegaProb = YES;
			doOmegaCat = YES;
			doOmegaRateHetCont = NO;
			doOmegaRateHetDisc = NO;
			
			sumPi = sumPi_b = 0.0;
			printf ("\n Introduce the number of categories");
			scanf ("%d", &numOmegaCat);
			if (numOmegaCat < 1 || numOmegaCat > 10)
				{
				printf ("PARAMETER ERROR: Bad number of omega categories (%d)\n\n", numOmegaCat);
				exit (-1);
				}
			omegaVal = (double *) calloc((numOmegaCat+1),(long) sizeof(double));
			if (!omegaVal)
				{
				fprintf (stderr, "PARAMETER ERROR: Could not allocate omega values of categories (%lu bytes)\n", (numOmegaCat+1) *(long) sizeof(double));
				exit (1);
				}
			omegaProb = (double *) calloc((numOmegaCat+1),(long) sizeof(double));
			if (!omegaProb)
				{
				fprintf (stderr, "PARAMETER ERROR: Could not allocate omega probabilities of categories (%lu bytes)\n", (numOmegaCat+1) *(long) sizeof(double));
				exit (1);
				}
			for (j=1; j<=numOmegaCat; j++)
				{
				printf ("\n Introduce the value of omega for the category %d", j);
				scanf("%lf", &omegaVal[j]);
				printf ("\n Introduce the value of omega for the category %d", j);
				scanf("%lf", &omegaProb[j]);
					
				if (omegaProb[j] > 1)
					{
					printf ("PARAMETER ERROR: Bad number of probabilities of omega categories (%lf)\n\n", omegaProb[j]);
					exit (-1);
					}
				sumPi = sumPi + omegaProb[j];
				}
					
			if (sumPi != 1) /* update probabilities of categories */
				{
				if (numOmegaCat == 1)
					omegaProb[1] = 1.0;
				else
					{
					for (j=1; j<=numOmegaCat; j++)
						{
						omegaProb[j]/=sumPi;
						sumPi_b = sumPi_b + omegaProb[j];
						}
					if ((int)sumPi_b != 1)
						{
						fprintf (stderr, "\n ERROR in the sum of probabilities of omega categories");
						exit (-1);
						}
					}
				}
			}
		
		if (modelNumber == 3) /* omega by a discrete hetereogeneous rate */
			{
			doOmegaCat = NO;
			doOmegaRateHetCont = NO;
			doOmegaRateHetDisc = YES;
			
			printf ("\n Introduce the number of categories");
			scanf ("%d", &numOmegaCat);
			if (numOmegaCat < 1 || numOmegaCat > 10)
				{
				printf ("PARAMETER ERROR: Bad number of omega categories (%d)\n\n", numOmegaCat);
				exit (-1);
				}
			printf ("\n Introduce the Omega rate heterogeneous");
			scanf("%lf", &OmegaRateHet);
			if (OmegaRateHet <= 0.0)
				{
				fprintf (stderr, "PARAMETER ERROR: Bad omega heterogeneous rate (%3.2f)\n\n", OmegaRateHet);
				PrintUsage();
				}
			printf ("\n Introduce the omega value");
			scanf("%lf", &omega);
			if (omega < 0.0)
				{
				fprintf (stderr, "PARAMETER ERROR: Bad omega (%3.2f)\n\n", omega);
				PrintUsage();
				}
			OmegaInit = omega;
			}
		
		if (modelNumber == 4)
			{
			doOmegaCat = NO;
			doOmegaRateHetCont = YES;
			doOmegaRateHetDisc = NO;
			
			printf ("\n Introduce the Omega rate heterogeneous");
			scanf("%lf", &OmegaRateHet);
			if (OmegaRateHet <= 0.0)
				{
				fprintf (stderr, "PARAMETER ERROR: Bad omega heterogeneous rate (%3.2f)\n\n", OmegaRateHet);
				PrintUsage();
				}
			printf ("\n Introduce the omega value");
			scanf("%lf", &omega);
			if (omega < 0.0)
				{
				fprintf (stderr, "PARAMETER ERROR: Bad omega (%3.2f)\n\n", omega);
				PrintUsage();
				}
			OmegaInit = omega;
			}
		
		printf ("\n Codon Model can to be from types: HKY, GTR or GTR not reversible. If you want GTR Model enter number 1, GTR not reversible enter 2, and for HKY Model enter another number");
		i = 0;
		scanf ("%d", &i);
		if (i == 1)
			{
			doCodon_NGTR = NO;
			doCodon_GTR = YES;
			doCodon_HKY = NO;
			printf ("\n You have to introduce the relative substitution rates to changes (rate matrix A-C A-G A-T C-G C-T G-T, for GTR models)(e.g. -v1 3 0.2 0.5 5 1)");
			scanf ("%lf", &Rmat[0]);			
			scanf ("%lf", &Rmat[1]);
			scanf ("%lf", &Rmat[2]);
			scanf ("%lf", &Rmat[3]);
			scanf ("%lf", &Rmat[4]);
			scanf ("%lf", &Rmat[5]);
			if (Rmat[5]!=1.0) 
				{
				for (j=0; j<5; j++) 
					Rmat[j]/=Rmat[5];
				Rmat[5]=1.0;
				}
			}
		else if (i == 2)
			{
			doHKY = NO;
			doCodon_NGTR = YES;
			doCodon_GTR = NO;
			doCodon_HKY = NO;
			printf("\n Introduce the vale rates for changes: AC CA AG GA AT TA CG GC CT TC GT");
			for (n=0;n<12;n++)
				scanf ("%lf", &NRmat[n]);
			if (NRmat[10]!=1.0) 
				{
				for (j=0; j<12 && j!=10; j++) 
					NRmat[j]/=NRmat[10];
				NRmat[10]=1.0;
				}
			}
		else
			{
			doCodon_NGTR = NO;
			doCodon_HKY = YES;
			doCodon_GTR = NO;
			}
		}
	
	printf ("\n proportion of invariable sites (e.g. -p0.6)");
	scanf ("%lf", &pinv);

	printf (" Noisy level ");
	printf ("\n   noisy = 0: does not print run information");
	printf ("\n   noisy = 1: + run settings");
	printf ("\n   noisy = 2: + calculation status and event information");
	printf ("\n   noisy = 3: + print ancestral status for each sequence at each event + MRCA");
	printf ("\n   noisy = 4: + gi's and Gi");
	printf ("\n   noisy = 5: + tree for each site");
	printf ("\n Noisy level (0-5) = ");
	scanf ("%d", &noisy);

   printf ("\n Do you want to simulate an outgroup (0 = no; 1 = yes)? ");
	scanf ("%d", &thereisOutgroup);

	if (thereisOutgroup == YES)
    	{       
    	printf ("\n Enter the length of the branch to the outgroup: ");
       scanf ("%lf", &outgroupBranchLength);
		}

	printf ("\n rate of exponential growth (e.g. 1e-6) = ");
	scanf ("%lf", &growthRate);
	fprintf (stderr, "\n If you want to introduce periods demographics enter 1 or enter another number if it doesn«t ");
	i = 0;
	scanf ("%d", &i);
	if (i == 1)
		{
		printf ("\n Enter number of periods");
		scanf ("%d", &numPeriods);
		if (numPeriods == 0)
			doDemographics = NO;
		else
			{
			doDemographics = YES;
			Nbegin =	(int *) calloc(numPeriods+1,(long) sizeof(int));
			Nend = 		(int *) calloc(numPeriods+1,(long) sizeof(int));
			cumDuration =	(int *) calloc(numPeriods+1,(long) sizeof(int));
			periodGrowth =	(double *) calloc(numPeriods+1,(long) sizeof(double));
			if (Nbegin == NULL || Nend == NULL || cumDuration == NULL)
				{
				fprintf (stderr, "PARAMETER ERROR: Could not allocate demographic vectors (%lu bytes)\n", numPeriods *(long) sizeof(int));
				exit (1);
				}
			printf ("\n Now, you have to introduce the Start, End and Duration (1 or 2 * N units) for each period:");
			for (j=1; j<=numPeriods; j++)
				{
				printf ("\n Period %d:", j);
				printf ("\n begin = ");
				scanf("%d", &Nbegin[j]);
				printf ("\n end = ");
				scanf("%d", &Nend[j]);
				printf ("\n duration = ");
				scanf("%d", &cumDuration[j]);
				cumDuration[j] = cumDuration[j] + cumDuration[j-1];
				}
			}
		}
	
   printf ("\n Enter the name of sequences file to be created: ");
	scanf ("%s", alignmentFile);
	i = 0;
	printf ("\n If you want Separated sets of sequences, you can introduce 1. If you do not, you can introduce another number");
	scanf ("%d", &i);
	if (i == 1)
		doSeparatedSequences = YES;
	else
		doSeparatedSequences = NO;
	i = 0;
	printf ("\n Enter number 1 to print the ancestral sequences in the sequences file or enter another number if it doesn«t ");
	scanf ("%d", &i);
	if (i == 1)
		doPrintAncestralSequences = YES;
	else
		doPrintAncestralSequences = NO;
	i = 0;
	
	printf ("\n\n Printed files:");
	printf ("\n Do you want an output Tree file?. Enter 1 to YES or press another number to NO");
	scanf ("%d", &i);
	if (i == 1)
		{
		printf ("\n Enter the name of tree file to be created: ");
		scanf ("%s", treeFile);
		doPrintTrees = YES;
		i = 0;
		}
	printf ("\n Do you want an output Times file?. Enter 1 to YES or press another number to NO");
	scanf ("%d", &i);
	if (i == 1)
		{
		printf ("\n Enter the name of times file to be created: ");
		scanf ("%s", timesFile);
		doPrintTimes = YES;
		i = 0;
		}
	printf ("\n Do you want an output breakpoints file?. Enter 1 to YES or press another number to NO");
	scanf ("%d", &i);
	if (i == 1)
		{	
		printf ("\n Enter the name of breakpoints file to be created: ");
		scanf ("%s", breakpointFile);
		doPrintBreakpoints = YES;
		}
	i = 0;
	
	printf ("\n\n Migration model. Enter 1 to YES or press another number to NO");
	scanf ("%d", &i);
	if (i == 1)
		{
		doMigration = YES;
		sumPopul = 0;
		
		printf ("\n First. You can enter the number of subpopulation/s");		
		scanf ("%d", &numPopulations);		
		if (numPopulations < 1 || numPopulations > numSequences) 
			{
			printf ("PARAMETER ERROR: Bad number of subpopulations (%d)\n\n", numPopulations);
			PrintUsage();
			}
				
		initPopulation = (int *) calloc(numPopulations+1,(long) sizeof(int));
		if (!initPopulation)
			{
			fprintf (stderr, "PARAMETER ERROR: Could not allocate initPopulation of migration model (%lu bytes)\n", (numPopulations+1) *(long) sizeof(int));
			exit (1);
			}
				
		printf ("\n Second. You can enter the proportions of tip sequences to each subpopulation/s. Each proportion will be a number between 0 and 1. The sum of all proportions will be 1\n");
		for (j=1; j<=numPopulations; j++)
			{
			scanf("%d", &initPopulation[j]);
			if (initPopulation[j] > numSequences || initPopulation[j] <= 0)
				{
				fprintf (stderr, "PARAMETER ERROR: Bad subpopulation proportion (%d)\n\n", initPopulation[j]);
				PrintUsage();
				}
			sumPopul = sumPopul + initPopulation[j];
			}
		printf ("\nNow you can enter the migration rate");	
		scanf("%lf", &migrationRate);
		if (migrationRate < 0)
			{
			fprintf (stderr, "PARAMETER ERROR: Bad migration rate (%lf)\n\n", migrationRate);
			PrintUsage();
			}
		if (sumPopul != numSequences)
			{
			fprintf (stderr, "PARAMETER ERROR: Bad initial subpopulation. The number of sequences is different to the sum of nodes in demes \n\n");
			PrintUsage();
			}
		}
	else
		doMigration = NO;
	i = 0;
	
	printf ("\n Enter number 1 to make the MRCA sequence from your input file or enter another number if you want to make the MRCA sequence from nucleotide frequencies");
	scanf ("%d", &i);
	if (i == 1)
		{
		printf ("\n Now, you have to introduce the name of your MRCA file \n (if you have problems with this you have to check the input file name and if the location of your file is correct)");
		scanf ("%s", MRCAFile);
		doMRCAFile = YES;
		}
	else
		doMRCAFile = NO;
	i = 0;
	

	printf ("\n Enter number 1 to print the sumMRCA/GMRCA sequence in separate files");
	scanf ("%d", &i);
	if (i == 1)
		{
		doOutMRCAfiles = YES;
		}
	i = 0;


	printf ("\n Finally, if you want to see more information about input parameters enter 1, to run the program enter another number");
	scanf ("%d", &i);
	if (i == 1)
		PrintUsage();
		
	printf ("\n\n");
}

#endif



void ReadUntil(FILE *fv, char stopChar, char *what)
{
	char ch;
	
	ch=fgetc(fv);
	while (!feof(fv) && ch!=stopChar) 
		ch=fgetc(fv);

	if (feof(fv) || ch!=stopChar) {
		fprintf(stderr, "%s missing", what);
		exit(0);
	}
}




/***************************** PrintTitle *******************************/
/* Prints program header */

static void PrintTitle(FILE *filep)
{
	FILE *file;
	file = filep;
	fprintf(filep,"%s - Simulating Codon Sequences with total Recombination", PROGRAM_NAME);
	fprintf(filep, "\nVersion %s", VERSION_NUMBER);
	fprintf(filep, "\nWritten by Miguel Arenas (miguelab@uvigo.es) and David Posada (dposada@uvigo.es)");
	fprintf(filep, "\nArenas, M. and Posada, D. 2010. Coalescent simulation of intracodon recombination. Genetics, 184(2):429-437.");
	fprintf(filep,"\n______________________________________________________________________________________________________________\n\n");
}







/***************************** PrintUsage *******************************/
/* Prints a short description of program usage */
static void PrintUsage()
{
	
#ifdef MPI
	if (rank==root)
		{
#endif
		fprintf (stderr, "\n\nUsage: %s%s [-n# -s# -=# -/# -l# -e# -_# -g# -p# -q# -%%# -r# -w# -u# -m# -f# -t# -v# -@# -a# -i# -o# -xseqGMRCA -b# -c -*# -z -$ -j# -k# -d# -+ -## -y# -? -h]", PROGRAM_NAME, VERSION_NUMBER);
		fprintf (stderr,"\n-n: number of replicates (e.g. -n1000)");
		fprintf (stderr,"\n-s: sample size (#haplotypes) (e.g. -s8)");
		fprintf (stderr,"\n-=: tip dates. (e.g. for 4 samples: 1995:node 1; 2003:nodes 4 and 6; 1997:nodes 2 and 3; 2001:nodes 7 and 8.   -=4 1995 1 1 2003 4 6 1997 2 3 2001 7 8)");
		fprintf (stderr,"\n-/: generation time (e.g. -/300)");
		fprintf (stderr,"\n-l: number of sites (bp or number of codons) (e.g. -l501)");
		fprintf (stderr,"\n-e: effective population size (e.g. -e1000)");
		fprintf (stderr,"\n-_: haploid/diploid scale. Type 1 is for haploid, and 2 is for diploid (e.g. -_2)");
		fprintf (stderr,"\n-g: rate of exponential growth (e.g. -g1e-6) (default is constant population size)");
		fprintf (stderr,"\n-p: number of demographic periods followed by effective population size at the beginning and at the end of the period, ");
		fprintf (stderr, "\n   and the duration of the period in generations. (e.g. -p1 100 200 30000) (e.g. -p2 100 100 40000 1000 2000 20000)");
		fprintf (stderr,"\n-q: migration model has three sets of numbers: (e.g. -q4 2 2 1 3 0.002) ");
		fprintf (stderr,"\n		the first number is the number of subpopulations that you want on tip nodes, at present time.");
		fprintf (stderr,"\n		the next numbers are the number of nodes for each subpopulation. The sum must to be equal to the number of sequences, sample size.");
		fprintf (stderr,"\n		the last number is the migration rate.");
		fprintf (stderr,"\n-%%: convergence demes events: (e.g. -%%3 1 2 200 3 4 1000 5 6 2000) ");
		fprintf (stderr,"\n		the first number is the number of convergence events. Then each 3 numbers are for one convergence event, where: the first and the second numbers ");
		fprintf (stderr,"\n		are the number of the demes, the third number is the time to the event. When two demes are converging the resulting deme was named as the next number of the demes list. ");
		fprintf (stderr,"\n		So, the number of the convergence demes cannot be repeated in the enter instruction. ");
		fprintf (stderr,"\n-r: homogeneous recombination rate (e.g. -r1e-6)");
		fprintf (stderr,"\n-w: fixed number of recombination events (e.g. -w3)");
		fprintf (stderr,"\n-u: mutation rate (e.g. -u1e-8). It's for nucleotide model or for codon model");
		fprintf (stderr,"\n-m: omega = dN/dS ratio (e.g. -m1 0.2). There are 7 models to intoroduce the omega value:");
		fprintf (stderr,"\n		1. Omega constant. Intoduce 1 and then the value of omega (e.g. -m1 0.9)");
		fprintf (stderr,"\n		2. Categories from the user: the max number of categories is 10, category probabibilities can to be more than 1 or negative. (e.g. -m2 2 0.2 0.4 1.5 0.6)");
		fprintf (stderr,"\n		   the first number is 2, and then the number of categories.");
		fprintf (stderr,"\n		   the next numbers are, for each category, first its omega value and then its probability.");
		fprintf (stderr,"\n		3. Discrete heterogeneous rate. Introduce 3, and then the number of categories, the discrete heterogeneous rate and the omega value (e.g. -m3 4 0.78 1.3)");
		fprintf (stderr,"\n		4. Continuous heterogeneous rate. Introduce 4, and then the continuous heterogeneous rate and the omega value (e.g. -m4 0.5 0.86)");
		fprintf (stderr,"\n		5. M1 model. Introduce 5, and then the proportion of sites for omega = 0 (e.g. -m5 0.5)");
		fprintf (stderr,"\n		6. M7 model. Introduce 6, and then p and q for the D. beta (e.g. -m6 0.4 0.3)");
		fprintf (stderr,"\n		7. M8 models. Introduce 7, and then: The proportion of sites of the D. beta, p for D. beta, q for D. beta, omega for P1 (e.g. -m7 0.4 0.3 0.6 1.4)");
		fprintf (stderr,"\n-f: base frequencies (ACGT) (e.g. -f4 0.4 0.3 0.2 0.1) or only for codons 3x4 (e.g. -f12 0.4 0.3 0.2 0.1 0.2 0.2 0.2 0.4 0.1 0.1 0.3 0.5)");
		fprintf (stderr,"\n-t: transition/transversion ratio. It's for HKY models (e.g. -t2)");
		fprintf (stderr,"\n-v: relative substitution rates (AC AG AT CG CT GT) (e.g. -v1 3 0.2 0.5 5 1)");
		fprintf (stderr,"\n-@: Bad general rate matrix (-rx x x x x x x x x x x) (AC CA AG GA AT TA CG GC CT TC GT=1 TG for GTR non reversible)");
		fprintf (stderr,"\n-a: alpha shape of the gamma distribution for rate heterogeneity (e.g. -a0.05)");
		fprintf (stderr,"\n-i: proportion of invariable sites (e.g. -p0.6)");
		fprintf (stderr,"\n-o: branch length to the outgroup (e.g. -o0.0325) (default is no outgroup)");
		fprintf (stderr,"\n-x: input MRCA file name. It cannot have stop codons in codon models (e.g. -xSeqGMRCA)");
		fprintf (stderr,"\n-b: output file name (e.g. -bsequences)");
		fprintf (stderr,"\n-c: include ancestral sequences in alignment (e.g. -c)");
		fprintf (stderr,"\n-*: print in phylip/fasta/nexus format. 1 is for phylip format, 2 for fasta format and 3 for nexus format (e.g. -*2)");
		fprintf (stderr,"\n-z: multiple aligment files (e.g. -z)");
		fprintf (stderr,"\n-$: print to files sumMRCAs and GMRCAs (e.g. -$)");
		fprintf (stderr,"\n-j: tree file name (e.g. -jtrees)");
		fprintf (stderr,"\n-k: times file name (e.g. -ktimes)");
		fprintf (stderr,"\n-d: breakpoints file name (e.g. -dbreakpoints)");
		fprintf (stderr,"\n-+: print simulate omegas per site (e.g. -+)");
		fprintf (stderr,"\n-y: noisy (e.g. -y2)");
		fprintf (stderr,"\n-#: seed (e.g. -#37864287)");
		fprintf (stderr,"\n-? -h: Print help");
		
	#ifdef MPI
		}
	#endif

	exit(-1);
}




/***************************** PrintDate *******************************/
/* Prints date and time */

static void PrintDate (FILE *filep)
{
	char *date;				/* define date */
	time_t now;
	FILE *file;
	file = filep;		/* define */
	
	now=time(NULL);
	date= ctime(&now);
	fprintf(filep,"%s\n",date);
}



/*********************** RandomizeArray ***********************/
/* Randomizes elements of a given array */

/*static void RandomizeArray (int array[],int first, int last)
{
	int i, r, tmp, length;
	
	length = last - first;
	for (i=first; i<last; i++) */        /* iterate through the array */
		/*{              
		r = Rand(length-(i-first)) + i; *//* get random number */
		
		/*tmp = array[r];*/               /* swap (change) elements.   e.j., array[2](=1) = array[4](goes from =3 to =1) */
		/*array[r] = array[i];
		array[i] = tmp ;
		}
}*/



/********************** Rand ****************************/
/* it returns random uniform in range 0...max-1  */

/*static int Rand(int max)
{
	double	rd;
	rd = Rndu();
	return (floor(rd*max));*/		/* floor returns the entire nearer to rd*max */
/*}*/


	
/********************** Rndu ****************************/
/* it returns random uniform                */

/*static double Rndu (void)
{
  static int x_rndu=11, y_rndu=23, z_rndu=137;
  double r;

  x_rndu = 171*(x_rndu%177) - 2*(x_rndu/177);
  y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
  z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
  if (x_rndu<0) x_rndu+=30269;
  if (y_rndu<0) y_rndu+=30307;
  if (z_rndu<0) z_rndu+=30323;
  r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;
  return (r-(int)r);
}*/





/**************		The next functions are for recombination, they return the new segments of the new ancestral nodes			*************/
/********************************************************************************************************************************************/

/************** recSegmentsGeneratesLeft ****************/
/* It builds the new segments on the left of the breakpoint */
static int		recSegmentsGeneratesLeft(int nodeValue, TreeSegment *s, TreeSegment *n, int numNuc, int whichSite, int *actSegIndex)
{
	int breakp;

	breakp = whichSite;
	if (s->sStart < breakp)
		{
		if (breakp >= s->sEnd)
			{
			if (s->sStart == 1 && s->sEnd == breakp)
				{
				s->before1 = n;
				n->before1 = NULL;
				n->before2 = NULL;
				n->after1 = s;
				n->after2 = NULL;
				n->sIndex = *actSegIndex;
				n->sStart = s->sStart;
				n->sEnd = breakp-1;				/* the breakpoint site will be in the other segment */
				n->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			else if (s->sStart == 1 && s->sEnd != breakp)
				{
				s->before1 = n;
				n->before1 = NULL;
				n->before2 = NULL;
				n->after1 = s;
				n->after2 = NULL;
				n->sIndex = *actSegIndex;
				n->sStart = s->sStart;
				n->sEnd = s->sEnd;	
				n->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			else if (s->sStart != 1 && s->sEnd == breakp)
				{
				s->before1 = n;
				n->before1 = NULL;
				n->before2 = NULL;
				n->after1 = s;
				n->after2 = NULL;
				n->sIndex = *actSegIndex;
				n->sStart = s->sStart;
				n->sEnd = breakp-1;				/* the breakpoint site will be in the other segment */
				n->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			else
				{
				s->before1 = n;
				n->before1 = NULL;
				n->before2 = NULL;
				n->after1 = s;
				n->after2 = NULL;
				n->sIndex = *actSegIndex;
				n->sStart = s->sStart;
				n->sEnd = s->sEnd;	
				n->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			}
		else if (s->sStart == 1)
			{
			if (s->sEnd != numNuc)
				{
				s->before1 = n;
				n->before1 = NULL;
				n->before2 = NULL;
				n->after1 = s;
				n->after2 = NULL;
				n->sIndex = *actSegIndex;
				n->sStart = s->sStart;
				n->sEnd = breakp-1;
				n->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			else
				{
				s->before1 = n;
				n->before1 = NULL;
				n->before2 = NULL;
				n->after1 = s;
				n->after2 = NULL;
				n->sIndex = *actSegIndex;
				n->sStart = s->sStart;
				n->sEnd = breakp-1;
				n->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			}
		else if (s->sEnd == numNuc)
			{
			if (s->sStart != 1)
				{
				s->before1 = n;
				n->before1 = NULL;
				n->before2 = NULL;
				n->after1 = s;
				n->after2 = NULL;
				n->sIndex = *actSegIndex;
				n->sStart = s->sStart;
				n->sEnd = breakp-1;
				n->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			else
				{
				fprintf (fpmpi, "\n\nWarning in recSegmentsGeneratesLeft1, it can't never here");			
				exit(-1);
				}
			}
		else if (s->sStart >= 1 && s->sEnd > breakp)
			{
			s->before1 = n;
			n->before1 = NULL;
			n->before2 = NULL;
			n->after1 = s;
			n->after2 = NULL;
			n->sIndex = *actSegIndex;
			n->sStart = s->sStart;
			n->sEnd = breakp-1;
			n->sIndexNode = nodeValue;
			
			*actSegIndex = *actSegIndex+1;
				
			return (1);
			}
		
		else
			{
			fprintf (fpmpi, "\n\nWarning in recSegmentsGeneratesLeft2, breakp = %d, it can't never here", breakp);			
			exit(-1);
			}
		}
	else		
		{
		fprintf (fpmpi, "\n\nWarning in recSegmentsGeneratesLeft3, it can't never here");			
		exit(-1);
		}

return (0);
}




/************** recSegmentsGeneratesRight ****************/
/* It builds the new segments on the right of the breakpoint */

static int		recSegmentsGeneratesRight(int nodeValue, TreeSegment *s, TreeSegment *m, int numNuc, int whichSite, int *actSegIndex)
{
	int breakp;
	
	breakp = whichSite;
	if (s->sEnd >= breakp)
		{
		if (breakp <= s->sStart)
			{
			if (s->sStart == breakp && s->sEnd == numNuc)
				{
				s->before2 = m;
				m->before1 = NULL;
				m->before2 = NULL;
				m->after1 = s;
				m->after2 = NULL;
				m->sIndex = *actSegIndex;
				m->sStart = s->sStart;
				m->sEnd = s->sEnd;	
				m->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			else if (s->sStart == breakp && s->sEnd != numNuc)
				{
				s->before2 = m;
				m->before1 = NULL;
				m->before2 = NULL;
				m->after1 = s;
				m->after2 = NULL;
				m->sIndex = *actSegIndex;
				m->sStart = s->sStart;
				m->sEnd = s->sEnd;	
				m->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			else if (s->sStart != breakp && s->sEnd == numNuc)
				{
				s->before2 = m;
				m->before1 = NULL;
				m->before2 = NULL;
				m->after1 = s;
				m->after2 = NULL;
				m->sIndex = *actSegIndex;
				m->sStart = s->sStart;
				m->sEnd = s->sEnd;	
				m->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;

				return (1);
				}
			else
				{
				s->before2 = m;
				m->before1 = NULL;
				m->before2 = NULL;
				m->after1 = s;
				m->after2 = NULL;
				m->sIndex = *actSegIndex;
				m->sStart = s->sStart;
				m->sEnd = s->sEnd;	
				m->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;

				return (1);
				}
			}
			
		else if (s->sStart == 1)
			{
			if (s->sEnd != numNuc)
				{
				s->before2 = m;
				m->before1 = NULL;
				m->before2 = NULL;
				m->after1 = s;
				m->after2 = NULL;
				m->sIndex = *actSegIndex;
				m->sStart = breakp;
				m->sEnd = s->sEnd;
				m->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			else
				{
				s->before2 = m;
				m->before1 = NULL;
				m->before2 = NULL;
				m->after1 = s;
				m->after2 = NULL;
				m->sIndex = *actSegIndex;
				m->sStart = breakp;
				m->sEnd = s->sEnd;
				m->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			}
		else if (s->sEnd == numNuc)
			{
			if (s->sStart != 1)
				{
				s->before2 = m;
				m->before1 = NULL;
				m->before2 = NULL;
				m->after1 = s;
				m->after2 = NULL;
				m->sIndex = *actSegIndex;
				m->sStart = breakp;
				m->sEnd = s->sEnd;
				m->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			else
				{
				fprintf (fpmpi, "\n\nWarning in recSegmentsGeneratesRight1, it can't never here");			
				exit(-1);
				}
			}
		else if (s->sStart <= breakp && s->sEnd >= breakp)
			{
			s->before2 = m;
			m->before1 = NULL;
			m->before2 = NULL;
			m->after1 = s;
			m->after2 = NULL;
			m->sIndex = *actSegIndex;
			m->sStart = breakp;
			m->sEnd = s->sEnd;
			m->sIndexNode = nodeValue;
			
			*actSegIndex = *actSegIndex+1;
				
			return (1);
			}
		else
			{
			fprintf (fpmpi, "\n\nWarning in recSegmentsGeneratesRight2, breakp = %d, it can't never here", breakp);			
			exit(-1);
			}
		}
	else
		{
		fprintf (fpmpi, "\n\nWarning in recSegmentsGeneratesRight3, it can't never here");			
		exit(-1);
		}

return (0);
}



/** Broken Codons **/
/************** recSegmentsGeneratesLeftBrokenCodon ****************/
/* It builds the new segments on the left of the breakpoint */
static int		recSegmentsGeneratesLeftBrokenCodon(int nodeValue, TreeSegment *s, TreeSegment *n, int numNuc, int whichSite, int LeftLess, int RightHigh, int *actSegIndex)
{
	int breakp;
	int HereLeftLess, HereRightHigh;

	HereLeftLess = LeftLess;
	HereRightHigh = RightHigh;

	breakp = whichSite;

	if (s->sStart < breakp)
		{
		if (breakp >= s->sEnd)
			{
			if (s->sStart == 1 && s->sEnd == breakp)
				{
				s->before1 = n;
				n->before1 = NULL;
				n->before2 = NULL;
				n->after1 = s;
				n->after2 = NULL;
				n->sIndex = *actSegIndex;
				n->sStart = s->sStart;
				n->sEnd = HereRightHigh;				/* the breakpoint site will be in the other segment */
				n->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			else if (s->sStart == 1 && s->sEnd != breakp)
				{
				s->before1 = n;
				n->before1 = NULL;
				n->before2 = NULL;
				n->after1 = s;
				n->after2 = NULL;
				n->sIndex = *actSegIndex;
				n->sStart = s->sStart;
				n->sEnd = HereRightHigh;	
				n->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			else if (s->sStart != 1 && s->sEnd == breakp)
				{
				s->before1 = n;
				n->before1 = NULL;
				n->before2 = NULL;
				n->after1 = s;
				n->after2 = NULL;
				n->sIndex = *actSegIndex;
				n->sStart = s->sStart;
				n->sEnd = HereRightHigh;				/* the breakpoint site will be in the other segment */
				n->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			else
				{
				s->before1 = n;
				n->before1 = NULL;
				n->before2 = NULL;
				n->after1 = s;
				n->after2 = NULL;
				n->sIndex = *actSegIndex;
				n->sStart = s->sStart;
				n->sEnd = HereRightHigh;	
				n->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			}
		else if (s->sStart == 1)
			{
			if (s->sEnd != numNuc)
				{
				s->before1 = n;
				n->before1 = NULL;
				n->before2 = NULL;
				n->after1 = s;
				n->after2 = NULL;
				n->sIndex = *actSegIndex;
				n->sStart = s->sStart;
				n->sEnd = HereRightHigh;
				n->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			else
				{
				s->before1 = n;
				n->before1 = NULL;
				n->before2 = NULL;
				n->after1 = s;
				n->after2 = NULL;
				n->sIndex = *actSegIndex;
				n->sStart = s->sStart;
				n->sEnd = HereRightHigh;
				n->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			}
		else if (s->sEnd == numNuc)
			{
			if (s->sStart != 1)
				{
				s->before1 = n;
				n->before1 = NULL;
				n->before2 = NULL;
				n->after1 = s;
				n->after2 = NULL;
				n->sIndex = *actSegIndex;
				n->sStart = s->sStart;
				n->sEnd = HereRightHigh;
				n->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			else
				{
				fprintf (fpmpi, "\n\nWarning in recSegmentsGeneratesLeft1, it can't never here");			
				exit(-1);
				}
			}
		else if (s->sStart >= 1 && s->sEnd > breakp)
			{
			s->before1 = n;
			n->before1 = NULL;
			n->before2 = NULL;
			n->after1 = s;
			n->after2 = NULL;
			n->sIndex = *actSegIndex;
			n->sStart = s->sStart;
			n->sEnd = HereRightHigh;
			n->sIndexNode = nodeValue;
			
			*actSegIndex = *actSegIndex+1;
				
			return (1);
			}
		
		else
			{
			fprintf (fpmpi, "\n\nWarning in recSegmentsGeneratesLeft2, breakp = %d, it can't never here", breakp);			
			exit(-1);
			}
		}
	else		
		{
		fprintf (fpmpi, "\n\nWarning in recSegmentsGeneratesLeft3, it can't never here");			
		exit(-1);
		}

return (0);
}




/************** recSegmentsGeneratesRightBrokenCodon ****************/
/* It builds the new segments on the right of the breakpoint */

static int		recSegmentsGeneratesRightBrokenCodon(int nodeValue, TreeSegment *s, TreeSegment *m, int numNuc, int whichSite, int LeftLess, int RightHigh, int *actSegIndex)
{
	int breakp;
	int HereLeftLess, HereRightHigh;
	
	HereLeftLess = LeftLess;
	HereRightHigh = RightHigh;

	breakp = whichSite;

	if (s->sEnd >= breakp)
		{
		if (breakp <= s->sStart)
			{
			if (s->sStart == breakp && s->sEnd == numNuc)
				{
				s->before2 = m;
				m->before1 = NULL;
				m->before2 = NULL;
				m->after1 = s;
				m->after2 = NULL;
				m->sIndex = *actSegIndex;
				m->sStart = HereLeftLess;
				m->sEnd = s->sEnd;	
				m->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			else if (s->sStart == breakp && s->sEnd != numNuc)
				{
				s->before2 = m;
				m->before1 = NULL;
				m->before2 = NULL;
				m->after1 = s;
				m->after2 = NULL;
				m->sIndex = *actSegIndex;
				m->sStart = HereLeftLess;
				m->sEnd = s->sEnd;	
				m->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			else if (s->sStart != breakp && s->sEnd == numNuc)
				{
				s->before2 = m;
				m->before1 = NULL;
				m->before2 = NULL;
				m->after1 = s;
				m->after2 = NULL;
				m->sIndex = *actSegIndex;
				m->sStart = HereLeftLess;
				m->sEnd = s->sEnd;	
				m->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;

				return (1);
				}
			else
				{
				s->before2 = m;
				m->before1 = NULL;
				m->before2 = NULL;
				m->after1 = s;
				m->after2 = NULL;
				m->sIndex = *actSegIndex;
				m->sStart = HereLeftLess;
				m->sEnd = s->sEnd;	
				m->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;

				return (1);
				}
			}
			
		else if (s->sStart == 1)
			{
			if (s->sEnd != numNuc)
				{
				s->before2 = m;
				m->before1 = NULL;
				m->before2 = NULL;
				m->after1 = s;
				m->after2 = NULL;
				m->sIndex = *actSegIndex;
				m->sStart = HereLeftLess;
				m->sEnd = s->sEnd;
				m->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			else
				{
				s->before2 = m;
				m->before1 = NULL;
				m->before2 = NULL;
				m->after1 = s;
				m->after2 = NULL;
				m->sIndex = *actSegIndex;
				m->sStart = HereLeftLess;
				m->sEnd = s->sEnd;
				m->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			}
		else if (s->sEnd == numNuc)
			{
			if (s->sStart != 1)
				{
				s->before2 = m;
				m->before1 = NULL;
				m->before2 = NULL;
				m->after1 = s;
				m->after2 = NULL;
				m->sIndex = *actSegIndex;
				m->sStart = HereLeftLess;
				m->sEnd = s->sEnd;
				m->sIndexNode = nodeValue;
				
				*actSegIndex = *actSegIndex+1;
				
				return (1);
				}
			else
				{
				fprintf (fpmpi, "\n\nWarning in recSegmentsGeneratesRight1, it can't never here");			
				exit(-1);
				}
			}
		else if (s->sStart <= breakp && s->sEnd >= breakp)
			{
			s->before2 = m;
			m->before1 = NULL;
			m->before2 = NULL;
			m->after1 = s;
			m->after2 = NULL;
			m->sIndex = *actSegIndex;
			m->sStart = HereLeftLess;
			m->sEnd = s->sEnd;
			m->sIndexNode = nodeValue;
			
			*actSegIndex = *actSegIndex+1;
				
			return (1);
			}
		else
			{
			fprintf (fpmpi, "\n\nWarning in recSegmentsGeneratesRight2, breakp = %d, it can't never here", breakp);			
			exit(-1);
			}
		}
	else
		{
		fprintf (fpmpi, "\n\nWarning in recSegmentsGeneratesRight3, it can't never here");			
		exit(-1);
		}

return (0);
}






/*********************** coalescence overlap segments ************************/
/* When 2 nodes have a coalescence, this function returns YES 
if the site is contained in some segment of both nodes (reduces MRCA) or NO 
if the site is not in both nodes (does not reduce MRCA)*/
static int overLapSegmentsCoalMRCA (TreeNode *p, TreeNode *q, int sizeNode_p, int sizeNode_q, int site)
{
	int a, b, position;
	TreeSegment *s, *n;
	
	position = a = b = 0;
	
			
	for (position = 0; position < sizeNode_p; position++)
		{
		s = segments + post(position,p->index,distance);
		if (site >= s->sStart && site <= s->sEnd)
			{
			a = 1;
			break;
			}
		}
			
	for (position = 0; position < sizeNode_q; position++)
		{
		n = segments + post(position,q->index,distance);
		if (site >= n->sStart && site <= n->sEnd)
			{
			b = 1;
			break;
			}
		}
			
	if (a == 1 && b == 1)
		return YES;
	else
		return NO;
							
}






/************************************************/
/*************** BUILDING TREES *****************/
/* These functions build the genealogies of the trees from a net. 
These functions belong to MakeCoalescenceTree function */

/************ buildTreeCoal *************/
/* It builds the tree when there aren't recombinations */
static void		buildTreeCoal (TreeNode *p, TreeNodex *f, int numSequences, int *numActNodex)
	{		
	TreeNodex *g, *h;
	
	g = NULL;
	h = NULL;
	
	if (p->left != NULL) 
		{
		*numActNodex = *numActNodex+1; 
		g = nodex + *numActNodex;
						
		g->index = *numActNodex;
		g->NetIndex = p->left->index;
		g->length = p->left->length;
		g->time = p->left->time;
		g->indexOldMigPop = p->left->indexOldMigPop;
		f->left = g;
		g->anc1 = f;
		}
	
	if (p->right != NULL)
		{
		*numActNodex = *numActNodex +1;
		h = nodex + *numActNodex;
				
		h->index = *numActNodex;
		h->NetIndex = p->right->index;
		h->length = p->right->length;
		h->time = p->right->time;
		h->indexOldMigPop = p->right->indexOldMigPop;
		f->right = h;
		h->anc1 = f;
		}
	
	if (p->left != NULL)
		buildTreeCoal (p->left, g, numSequences, numActNodex);
	if (p->right != NULL)
		buildTreeCoal (p->right, h, numSequences, numActNodex);
	}



/* Building with recombination */
/************* buildTreeInit ************/
/* It builds the initial tree of a net */
static void		buildTreeInit (TreeNode *p, TreeNodex *f, int numNuc, int *arrayIndBreakpointsOrd, int whoBreakp, int numSequences, int *numActNodex/*, int *numNodex*/)
	{		/*buildTreeInit(p, f, numNuc, arrayIndBreakpointsOrd, j, numSequences, &numActNodex);*/
	int hind, gind, a, w, j, x, k, step, numSegInit;
	int *arrayIndexTreeSegments;
	int numberNq, numberNr;
	
	TreeSegment *s;
	TreeNode *q, *r;
	TreeNodex *g, *h;
	
	r = NULL;
	q = NULL;
	g = NULL;
	h = NULL;
	w = x = k = a = j = hind = gind = step = numSegInit = numberNq = numberNr = 0;

	arrayIndexTreeSegments = (int *) calloc((p->numSegNode +1),(long) sizeof(int));
	if (!arrayIndexTreeSegments)
		{
		fprintf (fpmpi, "Could not allocate arrayIndexTreeSegments (%lu bytes)\n", (p->numSegNode*2 +1) *(long) sizeof(int));
		exit (-1);
		}
		
	for (w = 0; w < p->numSegNode; w++)
		{
		s = segments + post(w,p->index,distance);
		
		if (s->sStart == 1 && s->sEnd >= (arrayIndBreakpointsOrd[whoBreakp]-1) && s->after1 != NULL) /* first tree */
			{
			arrayIndexTreeSegments[x] = s->sIndex;
			x++;
			}
		}
	numSegInit = x;
	for (w = 0; w < numSegInit; w++) 
		{
		if (arrayIndexTreeSegments[w] < numSequences)
			step++;
		}
	
	if (step < numSegInit)
		{
		for (w = 0; w < p->numSegNode; w++)
			{
			s = segments + post(w,p->index,distance);
	
			if (s->sIndex == arrayIndexTreeSegments[a] && a <= numSegInit) 
				{
				if (k == 0) /* first segment that make a nodex */
					{
					if (s->after1 != NULL) 
						{
						gind++;
						*numActNodex = *numActNodex+1;
						
						g = nodex + *numActNodex;
												
						numberNq = s->after1->sIndexNode;
						q = nodes + numberNq;
												
						g->index = *numActNodex;

						g->MRCAfrom = 1;
						g->MRCAto = arrayIndBreakpointsOrd[whoBreakp]-1;

						g->NetIndex = q->index;
						g->length = q->length;
						g->time = q->time;
						g->indexOldMigPop = q->indexOldMigPop;
						f->left = g;
						g->anc1 = f;
						
						if (s->after2 != NULL)
							{
							if (s->after2->sIndexNode != s->after1->sIndexNode)
								{
								*numActNodex = *numActNodex +1;
								hind++;
								
								h = nodex + *numActNodex;
								numberNr = s->after2->sIndexNode;
								r = nodes + numberNr;

								h->index = *numActNodex;
								h->NetIndex = r->index;

								h->MRCAfrom = 1;
								h->MRCAto = arrayIndBreakpointsOrd[whoBreakp]-1;

								h->length = r->length;
								h->time = r->time;
								h->indexOldMigPop = r->indexOldMigPop;
								f->right = h;
								h->anc1 = f;
								}
							}
						}
					}
	
	
			if ((k > 0) && (gind == 0 || hind == 0)) /* others segments */
				{
				if (s->after1 != NULL) 
					{
					if (s->after2 == NULL && hind == 0 && s->after1->sIndexNode != numberNq)
						{
						hind++;	
						*numActNodex = *numActNodex+1;
								
						h = nodex + *numActNodex;
						numberNr = s->after1->sIndexNode;
						r = nodes + numberNr;
						
						h->index = *numActNodex;
						h->NetIndex = r->index;

						h->MRCAfrom = 1;
						h->MRCAto = arrayIndBreakpointsOrd[whoBreakp]-1;

						h->length = r->length;
						h->time = r->time;
						h->indexOldMigPop = r->indexOldMigPop;
						f->right = h;
						h->anc1 = f;
						}
					if (s->after2 != NULL)
						{
						if (hind == 0)	
							{
							if (s->after1->sIndexNode == numberNq && s->after2->sIndexNode != numberNq)
								{
								hind++;
								*numActNodex=*numActNodex+1;
								
								h = nodex + *numActNodex;
								numberNr = s->after2->sIndexNode;
								r = nodes + numberNr;
								
								h->index = *numActNodex;
								h->NetIndex = r->index;

								h->MRCAfrom = 1;
								h->MRCAto = arrayIndBreakpointsOrd[whoBreakp]-1;

								h->length = r->length;
								h->time = r->time;
								h->indexOldMigPop = r->indexOldMigPop;
								f->right = h;
								h->anc1 = f;
								}
							}
						}
					}
				}
				k++;
				if (a < numSegInit)
					w = 0;
				a++;
			}
			if (a >= numSegInit)
				break;
		}
		
		free (arrayIndexTreeSegments);
				
		if (q != NULL && gind != 0)
			buildTreeInit (q, g, numNuc, arrayIndBreakpointsOrd, whoBreakp, numSequences, numActNodex);
		if (r != NULL && hind != 0)
			buildTreeInit (r, h, numNuc, arrayIndBreakpointsOrd, whoBreakp, numSequences, numActNodex);
		}
	}






/************* buildTreeEnd ************/
/* It builds the final tree of a net */
static void		buildTreeEnd (TreeNode *p, TreeNodex *f, int numNuc, int *arrayIndBreakpointsOrd, int whoBreakp, int numSequences, int *numActNodex)
	{		
	int hind, gind, end_b, end_c, a, w, j, x, k, step, counter, numSegInit;
	int numberNr, numberNq;
	int *arrayIndexTreeSegments;

	TreeSegment *s;
	TreeNode *q, *r;
	TreeNodex *g, *h;
					
	r = NULL;
	q = NULL;
	g = NULL;
	h = NULL;
	w = x = k = a = j = counter = end_b = end_c = hind = gind = step = numSegInit = numberNr = numberNq = 0;
	
	arrayIndexTreeSegments = (int *) calloc((p->numSegNode +1),(long) sizeof(int));
	if (!arrayIndexTreeSegments)
		{
		fprintf (fpmpi, "Could not allocate arrayIndexTreeSegments (%lu bytes)\n", (p->numSegNode*2 +1) *(long) sizeof(int));
		exit (-1);
		}
	
	for (w = 0; w < p->numSegNode; w++) 
		{
		s = segments + post(w,p->index,distance);
	
		if (s->sStart <= arrayIndBreakpointsOrd[whoBreakp] && s->sEnd == numNuc && s->after1 != NULL)		/* last tree */		/* j = indNumRE-1 */
			{
			arrayIndexTreeSegments[x] = s->sIndex;
			x++;
			}
		}
	numSegInit = x;

	for (w = 0; w < numSegInit; w++) 
		{
		if (arrayIndexTreeSegments[w] < numSequences)
			step++;
		}

	
	if (step < numSegInit)
		{
		for (w = 0; w < p->numSegNode; w++)
			{
			s = segments + post(w,p->index,distance);
	
			if (s->sIndex == arrayIndexTreeSegments[a] && a <= numSegInit) 
				{
				if (k == 0) /* first segment that make a nodex */
					{
					if (s->after1 != NULL) 
						{				
						gind++;
						*numActNodex = *numActNodex+1;
						
						g = nodex + *numActNodex;
						numberNq = s->after1->sIndexNode;
						q = nodes + numberNq;

						g->index = *numActNodex;
						g->NetIndex = q->index;
						g->length = q->length;
						g->time = q->time;

						g->MRCAfrom = arrayIndBreakpointsOrd[whoBreakp];
						g->MRCAto = numNuc;

						g->indexOldMigPop = q->indexOldMigPop;
						f->left = g;
						g->anc1 = f;
						
						if (s->after2 != NULL)
							{
							if (s->after2->sIndexNode != s->after1->sIndexNode)
								{						
								*numActNodex = *numActNodex +1;
								hind++;
								
								h = nodex + *numActNodex;
								numberNr = s->after2->sIndexNode;
								r = nodes + numberNr;
								
								h->index = *numActNodex;
								h->NetIndex = r->index;

								h->MRCAfrom = arrayIndBreakpointsOrd[whoBreakp];
								h->MRCAto = numNuc;

								h->length = r->length;
								h->time = r->time;
								h->indexOldMigPop = r->indexOldMigPop;
								f->right = h;
								h->anc1 = f;
								}
							}
						}
					}
	
				if ((k > 0) && (gind == 0 || hind == 0)) /* others segments */
					{
					if (s->after1 != NULL) 
						{
						if (s->after2 == NULL && hind == 0 && s->after1->sIndexNode != numberNq)
							{
							hind++;	
							*numActNodex = *numActNodex+1;
							
							h = nodex + *numActNodex;
							numberNr = s->after1->sIndexNode;
							r = nodes + numberNr;
							
							h->index = *numActNodex;
							h->NetIndex = r->index;

							h->MRCAfrom = arrayIndBreakpointsOrd[whoBreakp];
							h->MRCAto = numNuc;

							h->length = r->length;
							h->time = r->time;
							h->indexOldMigPop = r->indexOldMigPop;
							f->right = h;
							h->anc1 = f;
							}
						if (s->after2 != NULL)
							{
							if (s->after1->sIndexNode == numberNq && s->after2->sIndexNode != numberNq && hind == 0) /* h aun no se ha creado */
								{
								hind++;
								*numActNodex=*numActNodex+1;
								
								h = nodex + *numActNodex;
								numberNr = s->after2->sIndexNode;
								r = nodes + numberNr;
								
								h->index = *numActNodex;
								h->NetIndex = r->index;
								
								h->MRCAfrom = arrayIndBreakpointsOrd[whoBreakp];
								h->MRCAto = numNuc;

								h->length = r->length;
								h->time = r->time;
								h->indexOldMigPop = r->indexOldMigPop;
								f->right = h;
								h->anc1 = f;
								}
								
							}
						}
					}
					k++;
					if (a < numSegInit)
						w = 0;
					a++; 
				}
				if (a >= numSegInit)
					break;
			}
	
		free (arrayIndexTreeSegments);
	
		if (q != NULL && gind == 1)
			buildTreeEnd (q, g, numNuc, arrayIndBreakpointsOrd, whoBreakp, numSequences, numActNodex);
		if (r != NULL && hind == 1)
			buildTreeEnd (r, h, numNuc, arrayIndBreakpointsOrd, whoBreakp, numSequences, numActNodex);
		}
	}








/************ buildTreeIntern *************/
/* It builds the internal trees of a net */
static void		buildTreeIntern (TreeNode *p, TreeNodex *f, int numNuc, int *arrayIndBreakpointsOrd, int whoBreakp, int numSequences, int *numActNodex)
	{		
	int hind, gind, a, w, j, x, k, step, numSegInit;
	int numberNq, numberNr;
	int *arrayIndexTreeSegments;

	TreeSegment *s;
	TreeNode *q, *r;
	TreeNodex *g, *h;
	
	r = NULL;
	q = NULL;
	g = NULL;
	h = NULL;
	w = x = k = a = j = hind = gind = step = numSegInit = numberNq = numberNr = 0;
	
	arrayIndexTreeSegments = (int *) calloc((p->numSegNode +1),(long) sizeof(int));
	if (!arrayIndexTreeSegments)
		{
		fprintf (fpmpi, "Could not allocate arrayIndexTreeSegments (%lu bytes)\n", (p->numSegNode*2 +1) *(long) sizeof(int));
		exit (-1);
		}
	
	for (w = 0; w < p->numSegNode; w++)
		{
		s = segments + post(w,p->index,distance);
	
		if (s->sStart <= (arrayIndBreakpointsOrd[whoBreakp]) && s->sEnd >= (arrayIndBreakpointsOrd[whoBreakp+1]-1) && s->after1 != NULL) /* internal trees */
			{
			arrayIndexTreeSegments[x] = s->sIndex;
			x++;
			}
		}
	numSegInit = x;
	for (w = 0; w < numSegInit; w++) 
		{
		if (arrayIndexTreeSegments[w] < numSequences)
			step++;
		}
	
	if (step < numSegInit)
		{
		for (w = 0; w < p->numSegNode; w++)
			{
			s = segments + post(w,p->index,distance);
	
			if (s->sIndex == arrayIndexTreeSegments[a] && a <= numSegInit)
				{
									
				if (k == 0) /* first segment that make a nodex */
					{
					if (s->after1 != NULL) 
						{
						gind++;
						*numActNodex = *numActNodex+1;
						
						g = nodex + *numActNodex;
						numberNq = s->after1->sIndexNode;
						q = nodes + numberNq;
						
						g->index = *numActNodex;

						g->MRCAfrom = arrayIndBreakpointsOrd[whoBreakp];
						g->MRCAto = arrayIndBreakpointsOrd[whoBreakp+1]-1;

						g->NetIndex = q->index;
						g->length = q->length;
						g->time = q->time;
						g->indexOldMigPop = q->indexOldMigPop;
						f->left = g;
						g->anc1 = f;
							
						if (s->after2 != NULL)
							{
							if (s->after2->sIndexNode != s->after1->sIndexNode)
								{
								*numActNodex = *numActNodex +1;
								hind++;
							
								h = nodex + *numActNodex;
								numberNr = s->after2->sIndexNode;
								r = nodes + numberNr;
								
								h->index = *numActNodex;
								h->NetIndex = r->index;

								h->MRCAfrom = arrayIndBreakpointsOrd[whoBreakp];
								h->MRCAto = arrayIndBreakpointsOrd[whoBreakp+1]-1;								

								h->length = r->length;
								h->time = r->time;
								h->indexOldMigPop = r->indexOldMigPop;
								f->right = h;
								h->anc1 = f;
								}
							}
						}
					}
	

				if ((k > 0) && (gind == 0 || hind == 0)) /* others segments */
					{
					if (s->after1 != NULL) 
						{
						if (s->after2 == NULL && hind == 0 && s->after1->sIndexNode != numberNq)
							{
							hind++;	
							*numActNodex = *numActNodex+1;
							
							h = nodex + *numActNodex;
							numberNr = s->after1->sIndexNode;
							r = nodes + numberNr;
							
							h->index = *numActNodex;
							h->NetIndex = r->index;

							h->MRCAfrom = arrayIndBreakpointsOrd[whoBreakp];
							h->MRCAto = arrayIndBreakpointsOrd[whoBreakp+1]-1;	

							h->length = r->length;
							h->time = r->time;
							h->indexOldMigPop = r->indexOldMigPop;
							f->right = h;
							h->anc1 = f;
							}
						if (s->after2 != NULL)
							{
							if (hind == 0 && s->after1->sIndexNode == numberNq && s->after2->sIndexNode != numberNq)
								{
								hind++;
								*numActNodex=*numActNodex+1;
							
								h = nodex + *numActNodex;
								numberNr = s->after2->sIndexNode;
								r = nodes + numberNr;
								
								h->index = *numActNodex;
								h->NetIndex = r->index;

								h->MRCAfrom = arrayIndBreakpointsOrd[whoBreakp];
								h->MRCAto = arrayIndBreakpointsOrd[whoBreakp+1]-1;	

								h->length = r->length;
								h->time = r->time;
								h->indexOldMigPop = r->indexOldMigPop;
								f->right = h;
								h->anc1 = f;
								}
							}
						}
					}
					k++;
					if (a < numSegInit)
						w = 0;
					a++;
				}
				if (a >= numSegInit)
					break;
			}
			
		free (arrayIndexTreeSegments);
	
		if (q != NULL && gind != 0)
			buildTreeIntern (q, g, numNuc, arrayIndBreakpointsOrd, whoBreakp, numSequences, numActNodex);
			
		if (r != NULL && hind != 0)
			buildTreeIntern (r, h, numNuc, arrayIndBreakpointsOrd, whoBreakp, numSequences, numActNodex);
		}
	}





/**************** Functions for print output files **************/

/****************** Print Segments Trees *****************/
/* Print unrooted trees to treefile in Newick format */
void PrintTreesSeg(int replicate, int indNumRE, int numNuc, int *arrayIndBreakpointsOrd)
	{
	int		step;
	
	if (numRE == 0)			/* there aren't recombinations */
		{
		step = 0;
		fprintf(fpTrees, "Dataset %d \n", replicate+1);
		fprintf (fpTrees, "Tree %05d_%05d [1-%d] = ", replicate+1, step+1, numNuc);
		/*fprintf(fpTrees,"Tree.%05d = ", replicate+1);*/
		WriteTreeSeg (treeRootNodex[0]);
		fprintf(fpTrees,");\n");
		fprintf(fpTrees,"\n\n");
		}
	else					/* there are recombinations */
		{
		fprintf(fpTrees, "Dataset %d \n", replicate+1);
		for (step = 0; step <= indNumRE; step++)
			{
			if (step == 0)	
				fprintf (fpTrees, "Tree %05d_%05d [1-%d] = ", replicate+1, step+1, arrayIndBreakpointsOrd[0]-1);
			if (step == indNumRE)
				fprintf (fpTrees, "Tree %05d_%05d [%d-%d] = ", replicate+1, step+1, arrayIndBreakpointsOrd[indNumRE-1], numNuc);
			if (step > 0 && step < indNumRE)
				fprintf (fpTrees, "Tree %05d_%05d [%d-%d] = ", replicate+1, step+1, arrayIndBreakpointsOrd[step-1], arrayIndBreakpointsOrd[step]-1); 
				
			WriteTreeSeg (treeRootNodex[step]);
			fprintf (fpTrees, ");\n");
			}
		fprintf(fpTrees,"\n\n");
		}
	/*fprintf (fpTrees,"\n");*/
	}



/******************* WriteTree ****************/
/* Writes a given (unrooted) tree from PrintTrees */

void WriteTreeSeg (TreeNodex *f)
	{		
	if (f != NULL)
		{
		if (doMigration == YES)
			{
			if (f->isOutgroup == YES)			/* outgroup */
				fprintf (fpTrees, ",outgroup_pop%d:%8.6f", f->indexOldMigPop, f->length*mutationRate);
			
			else if (f->left == NULL && f->right == NULL)		/* tip of the tree */
				fprintf (fpTrees, "seq%05d_pop%d:%8.6f", LabelSeg(f), f->indexOldMigPop, (f->anc1->time - f->time)*mutationRate);
			else								/* all ancester */
				{
				fprintf (fpTrees, "(");
				WriteTreeSeg (f->left);
				fprintf (fpTrees, ",");
				WriteTreeSeg (f->right);	
				if (f->anc1 != NULL)
					fprintf (fpTrees, "):%8.6f",(f->anc1->time - f->time)*mutationRate);
				WriteTreeSeg (f->outgroup);	
				}
			}
		else
			{
			if (f->isOutgroup == YES)			/* outgroup */
				fprintf (fpTrees, ",outgroup:%8.6f", f->length*mutationRate);
			
			else if (f->left == NULL && f->right == NULL)		/* tip of the tree */
				fprintf (fpTrees, "seq%05d:%8.6f", LabelSeg(f), (f->anc1->time - f->time)*mutationRate);
			else								/* all ancester */
				{
				fprintf (fpTrees, "(");
				WriteTreeSeg (f->left);
				fprintf (fpTrees, ",");
				WriteTreeSeg (f->right);	
				if (f->anc1 != NULL)
					fprintf (fpTrees, "):%8.6f",(f->anc1->time - f->time)*mutationRate);
				WriteTreeSeg (f->outgroup);	
				}
			}
		}
	}




/********************* PrintTimes **********************/
/* Prints to timesfile a detailed description of
the tree: nodes, times, branch lengths */

void PrintTimesSeg(int replicate, int indNumRE, int numNuc, int *arrayIndBreakpointsOrd)
	{

	int		step;
	step = 0;
	
	if (doMigration == NO)
		{
		if (numRE == 0)				/* there aren't recombinations */
			{
			fprintf (fpTimes, "\n\nDataset_%d Fragment_1 [1-%d]", replicate+1, numNuc);
			fprintf (fpTimes, "\n----------------- Nodes -----------------------");
			fprintf (fpTimes, "\n  class  label index (left right anc) NetNode |        time    time length   branch lenght");
			fprintf (fpTimes, "\n-----------------------------------------------------------------------------------------\n");
			ListTimesSeg (treeRootNodex[0]);
			}
		else						/* there are recombinations */
			{
			for (step = 0; step <= indNumRE; step++)
				{
				if (step == 0)
					fprintf (fpTimes, "\n\nDataset_%d Fragment_%d [1-%d]", replicate+1, step+1, arrayIndBreakpointsOrd[0]-1);
				if (step == indNumRE)
					fprintf (fpTimes, "\n\nDataset_%d Fragment_%d [%d-%d]", replicate+1, step+1, arrayIndBreakpointsOrd[indNumRE-1], numNuc);
				if (step > 0 && step < indNumRE)
					fprintf (fpTimes, "\n\nDataset_%d Fragment_%d [%d-%d]", replicate+1, step+1, arrayIndBreakpointsOrd[step-1], arrayIndBreakpointsOrd[step]-1);
		
				fprintf (fpTimes, "\n----------------- Nodes -------------------");
			fprintf (fpTimes, "\n  class  label index (left right anc) NetNode |        time    time length   branch lenght");
				fprintf (fpTimes, "\n-----------------------------------------------------------------------------------------\n");
				
				ListTimesSeg (treeRootNodex[step]);
				}
			}
		}
	else
		{
		if (numRE == 0)				/* there aren't recombinations */
			{
			fprintf (fpTimes, "\n\nDataset_%d Fragment_1 [1-%d]", replicate+1, numNuc);
			fprintf (fpTimes, "\n-------------------- Nodes -------------------------");
			fprintf (fpTimes, "\n  class  label index (left right anc) NetNode deme |        time    time length   branch lenght");
			fprintf (fpTimes, "\n-----------------------------------------------------------------------------------------\n");
			ListTimesSeg (treeRootNodex[0]);
			}
		else						/* there are recombinations */
			{
			for (step = 0; step <= indNumRE; step++)
				{
				if (step == 0)
					fprintf (fpTimes, "\n\nDataset_%d Fragment_%d [1-%d]", replicate+1, step+1, arrayIndBreakpointsOrd[0]-1);
				if (step == indNumRE)
					fprintf (fpTimes, "\n\nDataset_%d Fragment_%d [%d-%d]", replicate+1, step+1, arrayIndBreakpointsOrd[indNumRE-1], numNuc);
				if (step > 0 && step < indNumRE)
					fprintf (fpTimes, "\n\nDataset_%d Fragment_%d [%d-%d]", replicate+1, step+1, arrayIndBreakpointsOrd[step-1], arrayIndBreakpointsOrd[step]-1);
		
				fprintf (fpTimes, "\n-------------------- Nodes -------------------------");
				fprintf (fpTimes, "\n  class  label index (left right anc) NetNode deme |        time    time length   branch lenght");
				fprintf (fpTimes, "\n-----------------------------------------------------------------------------------------\n");
				
				ListTimesSeg (treeRootNodex[step]);
				}
			}
		}
	
	}






/******************* ListTimes ****************/
/* It makes a list of the nodes with information as time, lenght.. */
void ListTimesSeg (TreeNodex *f)
	{
	if (doMigration == NO)
		{
		if (f != NULL)			
			{
			if (f->isOutgroup == YES)			/* Outgroup */
				fprintf (fpTimes, "%8s  %4d  %4d (%4d %4d %4d)  %5d  |  %10.2lf     %10.2lf      %10.4lf\n", 
						"outgroup", LabelSeg(f), IndexSeg(f), IndexSeg(f->left), IndexSeg(f->right), IndexSeg(f->anc1), f->NetIndex, f->time, f->length, f->length*mutationRate);
			else if (f->anc1 != NULL && f->left != NULL && f->right != NULL)				/* No MRCA, no tip (internal ancester) */
				fprintf (fpTimes, "%8s  %4d  %4d (%4d %4d %4d)  %5d  |  %10.2lf     %10.2lf      %10.4lf\n", 
						"internal", LabelSeg(f), IndexSeg(f), IndexSeg(f->left), IndexSeg(f->right), IndexSeg(f->anc1), f->NetIndex, f->time, f->anc1->time - f->time, (f->anc1->time - f->time)*mutationRate);
			else if (f->anc1 != NULL && f->left == NULL && f->right == NULL)				/* tip */
				fprintf (fpTimes, "%8s  %4d  %4d (%4d %4d %4d)  %5d  |  %10.2lf     %10.2lf      %10.4lf\n", 
						"tip", LabelSeg(f), IndexSeg(f), IndexSeg(f->left), IndexSeg(f->right), IndexSeg(f->anc1), f->NetIndex, f->time, f->anc1->time - f->time, (f->anc1->time - f->time)*mutationRate);
			else if (f->anc1 == NULL && f->left != NULL && f->right != NULL)				/* root, MRCA */
				fprintf (fpTimes, "%8s  %4d  %4d (%4d %4d %4d)  %5d  |  %10.2lf     %10.2lf      %10.4lf\n", 
						"root", LabelSeg(f), IndexSeg(f), IndexSeg(f->left), IndexSeg(f->right), IndexSeg(f->anc1), f->NetIndex, f->time, 0.0, 0.0);
			else
				fprintf (fpTimes," ");
		
			ListTimesSeg (f->left);
			ListTimesSeg (f->right);
			ListTimesSeg (f->outgroup);
			}
		}
	else
		{
		if (f != NULL)			
			{
			if (f->isOutgroup == YES)			/* Outgroup */
				fprintf (fpTimes, "%8s %4d %4d   (%4d %4d %4d) %5d   %2d   |  %10.2lf     %10.2lf      %10.4lf\n", 
						"outgroup", LabelSeg(f), IndexSeg(f), IndexSeg(f->left), IndexSeg(f->right), IndexSeg(f->anc1), f->NetIndex,  f->indexOldMigPop, f->time, f->length, f->length*mutationRate);
			else if (f->anc1 != NULL && f->left != NULL && f->right != NULL)				/* No MRCA, no tip (internal ancester) */
				fprintf (fpTimes, "%8s %4d %4d   (%4d %4d %4d) %5d   %2d   |  %10.2lf     %10.2lf      %10.4lf\n", 
						"internal", LabelSeg(f), IndexSeg(f), IndexSeg(f->left), IndexSeg(f->right), IndexSeg(f->anc1), f->NetIndex,  f->indexOldMigPop, f->time, f->anc1->time - f->time, (f->anc1->time - f->time)*mutationRate);
			else if (f->anc1 != NULL && f->left == NULL && f->right == NULL)				/* tip */
				fprintf (fpTimes, "%8s %4d %4d   (%4d %4d %4d) %5d   %2d   |  %10.2lf     %10.2lf      %10.4lf\n", 
						"tip", LabelSeg(f), IndexSeg(f), IndexSeg(f->left), IndexSeg(f->right), IndexSeg(f->anc1), f->NetIndex,  f->indexOldMigPop, f->time, f->anc1->time - f->time, (f->anc1->time - f->time)*mutationRate);
			else if (f->anc1 == NULL && f->left != NULL && f->right != NULL)				/* root, MRCA */
				fprintf (fpTimes, "%8s %4d %4d   (%4d %4d %4d) %5d   %2d   |  %10.2lf     %10.2lf      %10.4lf\n", 
						"root", LabelSeg(f), IndexSeg(f), IndexSeg(f->left), IndexSeg(f->right), IndexSeg(f->anc1), f->NetIndex,  f->indexOldMigPop, f->time, 0.0, 0.0);
			else
				fprintf (fpTimes," ");
		
			ListTimesSeg (f->left);
			ListTimesSeg (f->right);
			ListTimesSeg (f->outgroup);
			}
		}
	}





/***************** Dex ***************/
/* Returns index for a given node */

int IndexSeg (TreeNodex *f)
{
	return (f == NULL) ? -1 : f->index+1; /* If the node haven't got bond => index = -1, else index = index+1 */
}



/***************** Lab ***************/
/* Returns label for a given node */

int LabelSeg (TreeNodex *f)
{
	return (f->anc1 == NULL && f->left == NULL && f->right == NULL) ? -1 : f->label+1; /* If the node haven't got ancester and descendants => label = -1, else label = label+1 */
}





/**************** RelabelNodes **************/
/*	After getting rid of superfluos node, we
	need to relabel those so they are consecutive 
	Use the indexes as labels when there
	is recombination
*/
void RelabelNodesSeg (TreeNodex *f)
	{
	if (f != NULL)
		{
		RelabelNodesSeg (f->left);
		RelabelNodesSeg (f->right);
		
		
		if (f->left == NULL && f->right == NULL) /* is tip */
			f->label = /*tipLabel++;*/f->NetIndex;
		else /* all ancester */
			f->label = intLabel++;/*f->NetIndex*/;

		/*if (f->left != NULL && f->right != NULL && f->anc1 == NULL)
			{
			fprintf (fpmpi, "\n Relabel, MRCA, f->label = %d, f->index = %d \n", f->label, f->index);
			}*/

		/*if (f->NetIndex == 0)
				fprintf(stderr,"\n\n  ********************** f->NetIndex = 0 y f->label = %d \n\n", f->label);*/
			
		/*if (tipLabel > intLabel)
			{
			fprintf (fpmpi, "\n\nWarning in RelabelNodesSeg, intLabel %d and tipLabel %d", intLabel, tipLabel);
			exit (-1);
			}*/
		}
	}
	
	














