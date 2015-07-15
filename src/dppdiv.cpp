/* 
 * DPPDiv version 1.1b source code (https://github.com/trayc7/FDPPDIV)
 * Copyright 2009-2013
 * Tracy Heath(1,2,3) 
 * Mark Holder(1)
 * John Huelsenbeck(2)
 *
 * (1) Department of Ecology and Evolutionary Biology, University of Kansas, Lawrence, KS 66045
 * (2) Integrative Biology, University of California, Berkeley, CA 94720-3140
 * (3) email: tracyh@berkeley.edu
 *
 * Also: T Stadler, D Darriba, AJ Aberer, T Flouri, F Izquierdo-Carrasco, and A Stamatakis
 *
 * DPPDiv is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License (the file gpl.txt included with this
 * distribution or http://www.gnu.org/licenses/gpl.txt for more
 * details.
 *
 * Some of this code is from publicly available source by John Huelsenbeck and Fredrik Ronquist
 *
 */


#include <iostream>
#include <string>
#include <cstring>
#include "Alignment.h"
#include "MbRandom.h"
#include "Mcmc.h"
#include "Model.h"
#include "util.h"

using namespace std;

void printHelp(bool files);

void printHelp(bool files) 
{
	if(files){
		
		cout << "*****\n";
		cout << "\nFile formats: \n\n";
		cout << "Tree file (newick format)\n";
		cout << "-----------------------------------------------\n";
		cout << "(((T1,T2),T3),(T4,T5));\n\n";
		cout << "-----------------------------------------------\n\n";
		cout << "Data file (phylip format)\n";
		cout << "-----------------------------------------------\n";
		cout << "5 10\n";
		cout << "T1\tTTCTTAGATT\n";
		cout << "T2\tTTATTAGATT\n";
		cout << "T3\tTTCCCAGATT\n";
		cout << "T4\tTTGCTAGATT\n";
		cout << "T5\tTTGCTAGATT\n";
		cout << "-----------------------------------------------\n\n";
		cout << "Node Calibration file \n-U = uniform with min and max bounds\n-E = offset exponentail with min age\n";
		cout << "-----------------------------------------------\n";
		cout << "3\n";
		cout << "-E\troot\t25.01\n";
		cout << "-U\tT1\tT3\t14.22\t20.532\n";
		cout << "-E\tT4\tT5\t4.773\n";
		cout << "-----------------------------------------------\n\n";
		cout << "*****\n";
	}
	else{
		cout << "\n\texample:      \n\n\t$ dppdiv -in datafile.in -out file -tre tree.phy -n 10000 -sf 10\n\n";
		cout << "\tHere are the available options that you can change (default values are in []):\n";
		cout << "\t\t-h    : print this menu **\n";
		cout << "\t\t-hf   : display example file formats **\n";
		cout << "\t\t-in   : Input file name **\n";
		cout << "\t\t-out  : output file name prefix **\n";
		cout << "\t\t-tre  : tree file name **\n";
		cout << "\t\t-pm   : prior mean of number of rate categories [= 1.0]\n";
		cout << "\t\t-ra   : shape for gamma of rates [= 2.0]\n";
		cout << "\t\t-rb   : scale for gamma of rates [= 4.0]\n";
		cout << "\t\t-hsh  : shape for gamma hyper prior on alpha concentration parameter [= 2.0]\n";
		cout << "\t\t-hsc  : scale for gamma hyper prior on alpha concentration parameter [calculated from prior mean on categories]\n";
		cout << "\t\t-n    : Number of MCMC cycles [= 1000000]\n";
		cout << "\t\t-pf   : print frequency [= 100] \n";
		cout << "\t\t-sf   : sample frequency [= 100] \n";
		cout << "\t\t-s1   : seed 1 (use this if you only pass in one seed) \n";
		cout << "\t\t-s2   : seed 2 \n";
		cout << "\t\t-ubl  : use input branch lengths \n";
		cout << "\t\t-snm  : single node move is turned on \n";
		cout << "\t\t-rdn  : random-order node moves \n";
		cout << "\t\t-offm : turn off a given move (1,2,3,4,5,6) \n";
		cout << "\t\t-rnp  : return 0.0 for lnl, run under prior \n";
		cout << "\t\t-cal  : file name with internal node calibrations \n";
		cout << "\t\t-npr  : 1=uniform, 2=yule, 3=cbd, 4=cbd fix with given vals\n";
		cout << "\t\t-bdr  : initial diversification rate (lambda - mu)\n";
		cout << "\t\t-bda  : initial relative death rate (mu / lambda)\n"; // what about bds?
        //cout << "\t\t-rho  : extant species sampling (fixed) -fofbd only\n"; //rw: probably should be fixed
		cout << "\t\t-soft : turn on soft bounds on calibrated nodes\n";
		cout << "\t\t-clok : run under strict clock (and estimate substitution rate)\n";
		cout << "\t\t-urg  : run under uncorrelated gamma-distributed rates\n";
		cout << "\t\t-exhp : all calibrated nodes are offset exponential and this also turns on the hyperprior on the exp rates\n";
		cout << "\t\t-dphp : all cal nodes have a DPM hyperprior this also gets a value for the expected # of calibration clusters\n";
		cout << "\t\t-ghp  : hyperprior on calibrations from a gamma\n";
		cout << "\t\t-mup  : modify update probabilities mid run\n";
		cout << "\t\t-fxm  : fix some model params\n";
		cout << "\t\t-ihp  : run under independent hyperprior on exp cals\n";
        cout << "\t\t-po   : print origin to log file\n"; //rw
        //cout << "\t\t-pfat : print fossil attachment times (zf) to log file \n"; //rw: right now this also prints yf times (probably not neccessary)
		cout << "\t\t** required\n\n";
	}
}

int main (int argc, char * const argv[]) {

	// read user settings
	
	seedType s1			= 0;
	seedType s2			= 0;
	string dataFileName	= "";
	string treeFileName = "";
	string calibFN		= "";
	string tipDateFN	= "";
	string outName		= "out";
	double priorMean    = 3.0;		// prior mean number of rate cats
	double rateSh       = 2.0;		// shape param for gamma dist on rates
	double rateSc       = 4.0;		// scale param for gamma dist on rates
	double hyperSh		= -1.0;		// scale hyperparam for gamma dist on concentration param
	double hyperSc		= -1.0;		// scale hyperparam for gamma dist on concentration param
	double netDiv		= -1.0;		// initial diversificaton rate (lambda - mu)
	double relDeath		= -1.0;		// initial relative death rate (mu / lambda)
	double ssbdPrS		= -1.0;
	double fixclokrt	= -1.0;		// fix the clock rate to this
    double rho          = 1.0;      //rw: extant species sampling
	int offmove			= 0;		// used to turn off one particular move
	int printFreq		= 100;
	int sampleFreq		= 100;
	int numCycles		= 1000000;
	int burn			= 10000;	// not used
	int treeNodePrior	= 1;
	bool userBLs		= false;	// initialize tree with user branch lenghts
	bool writeDataFile	= false;	// write moves to info.out file
	bool verbose		= false;	// output to logger this isn't used at the moment
	bool runPrior		= false;	// causes lnl calc to return 0.0 so this is run just under the prior
	bool justTree		= false;	// just print a uniformized tree and exit
	bool rndNdMv		= false;	// turns on randomized node move
	bool moveAllN		= true;		// if false, then node move only moves a single node
	bool rootfix		= true;
	bool printalign		= false;
	bool softbnd		= false;
	bool calibHyP		= false;
	bool dpmExpHyp		= false;
	bool gammaExpHP		= false;
	bool modUpdatePs	= false;
	bool fixModelPs		= false;
	int dpmEHPPrM		= 3;
	int modelType		= 1;		// 1 = DPP, 2 = strict clock, 3 = uncorrelated-gamma
	bool fixClockToR    = false;
	bool indHP			= false;
	bool doAbsRts		= false;
	bool fixTest		= false;
    bool printOrigin    = false;    //rw: print origin to log file
    bool printAttach    = false;     //rw: maybe issue a warning if cal file is large & printAttach = true
	int numIntervals	= 1;
	double originMax		= 0.0;
	
	if(argc > 1){
		for (int i = 1; i < argc; i++){
			char *curArg = argv[i];
			if(strlen(curArg) > 1 && curArg[0] == '-'){
				if(!strcmp(curArg, "-in"))
					dataFileName = argv[i+1];
				else if(!strcmp(curArg, "-out"))
					outName = argv[i+1];
				else if(!strcmp(curArg, "-tre"))
					treeFileName = argv[i+1];
				else if(!strcmp(curArg, "-pm")){
					priorMean = atof(argv[i+1]);
					modelType = 1;
				}
				else if(!strcmp(curArg, "-ra"))
					rateSh = atof(argv[i+1]);
				else if(!strcmp(curArg, "-rb"))
					rateSc = atof(argv[i+1]);
				else if(!strcmp(curArg, "-n"))
					numCycles = atoi(argv[i+1]);
				else if(!strcmp(curArg, "-pf"))
					printFreq = atoi(argv[i+1]);
				else if(!strcmp(curArg, "-sf"))
					sampleFreq = atoi(argv[i+1]);
				else if(!strcmp(curArg, "-s1"))
					s1 = (seedType) atoi(argv[i+1]);
				else if(!strcmp(curArg, "-s2"))
					s2 = (seedType) atoi(argv[i+1]);
				else if(!strcmp(curArg, "-bi"))
					burn = atoi(argv[i+1]);
				else if(!strcmp(curArg, "-v"))
					verbose = true;
				else if(!strcmp(curArg, "-ubl"))
					userBLs = true;
				else if(!strcmp(curArg, "-snm"))
					moveAllN = false;
				else if(!strcmp(curArg, "-rdn"))
					rndNdMv = true;
				else if(!strcmp(curArg, "-vb"))
					writeDataFile = true;
				else if(!strcmp(curArg, "-offm"))
					offmove = atoi(argv[i+1]);
				else if(!strcmp(curArg, "-hsh"))
					hyperSh = atof(argv[i+1]);
				else if(!strcmp(curArg, "-hsc"))
					hyperSc = atof(argv[i+1]);
				else if(!strcmp(curArg, "-rnp"))
					runPrior = true;
				else if(!strcmp(curArg, "-tfu"))
					justTree = true;
				else if(!strcmp(curArg, "-cal"))
					calibFN = argv[i+1];
				else if(!strcmp(curArg, "-tip"))
					tipDateFN = argv[i+1];
				else if(!strcmp(curArg, "-npr"))
					treeNodePrior = atoi(argv[i+1]);
				else if(!strcmp(curArg, "-tgs"))		// set treeNodePrior to do calibrated birth-death
					treeNodePrior = 6;
				else if(!strcmp(curArg, "-tga"))		// set treeNodePrior to do calibrated birth-death with ancestor fossils
					treeNodePrior = 7;
                else if(!strcmp(curArg, "-fbds")) // condition fbd on the age of the origin
                    treeNodePrior = 8;
                else if(!strcmp(curArg, "-fofbd")) // fossil only fossilised birth death process
                    treeNodePrior = 9;
				else if(!strcmp(curArg, "-bdr"))	// (lambda - mu)
					netDiv = atof(argv[i+1]);
				else if(!strcmp(curArg, "-bda"))	// (mu / lambda)
					relDeath = atof(argv[i+1]);
				else if(!strcmp(curArg, "-bds"))	// (mu / lambda) //rw: is this correct, is this not psi? how do these flags work?
					ssbdPrS = atof(argv[i+1]);
                else if(!strcmp(curArg, "-rho"))	//rw: extant species sampling
                    rho = atof(argv[i+1]);
				else if(!strcmp(curArg, "-fix")){	// fix clock
					fixclokrt = atof(argv[i+1]);
					fixClockToR = true;
				}
				else if(!strcmp(curArg, "-res"))
					rootfix = false;
				else if(!strcmp(curArg, "-pal"))
					printalign = false;
				else if(!strcmp(curArg, "-soft"))
					softbnd = true;
				else if(!strcmp(curArg, "-clok")){  // run under strict clock
					if(!fixClockToR) fixclokrt = -1.0;
					modelType = 2;
				}
				else if(!strcmp(curArg, "-urg"))  // run under uncorrelated-gamma rates
					modelType = 3;
				else if(!strcmp(curArg, "-exhp"))
					calibHyP = true;
				else if(!strcmp(curArg, "-dphp")){
					calibHyP = true;
					dpmExpHyp = true;
					dpmEHPPrM = atoi(argv[i+1]);
				}
				else if(!strcmp(curArg, "-ghp"))
					gammaExpHP = true;
				else if(!strcmp(curArg, "-mup"))
					modUpdatePs = true;
				else if(!strcmp(curArg, "-fxm"))
					fixModelPs = true;
				else if(!strcmp(curArg, "-ihp"))
					indHP = true;
				else if(!strcmp(curArg, "-abs")){
					doAbsRts = true;
					rateSh = 1.0;
				}
				else if(!strcmp(curArg, "-fxtr")){
					fixTest = true;
				}
				else if(!strcmp(curArg, "-h")){
					printHelp(false);
					return 0;
				}
				else if(!strcmp(curArg, "-hf")){
					printHelp(true);
					return 0;
				}
                else if(!strcmp(curArg, "-po")){//RW
                    printOrigin = true;
                }
                else if(!strcmp(curArg, "-pfat")){//RW
                    printAttach = true;
                }
                else if(!strcmp(curArg, "-omax")){//RW
                    originMax = atof(argv[i+1]);
                }
				else if(!strcmp(curArg, "-sky"))
					numIntervals = atoi(argv[i+1]);
				else {
					cout << "\n############################ !!! ###########################\n";
					cout << "\n\n\tPerhaps you mis-typed something, here are the \n\tavailable options:\n";
					printHelp(false);
					cout << "\n############################ !!! ###########################\n";
					return 0;
				}
			}
		}
	}
	
	else {
		cout << "\n############################ !!! ###########################\n";
		cout << "\n\n\tPlease specify data and tree files, here are the \n\tavailable options:\n";
		printHelp(false);
		cout << "\n############################ !!! ###########################\n";
		return 0;
	}
	
    //rw: suppress tree requirement
    if(treeNodePrior != 9) {
        if(dataFileName.empty() || treeFileName.empty()){
            cout << "\n############################ !!! ###########################\n";
            cout << "\n\n\tPlease specify data and tree files, here are the \n\tavailable options:\n";
            printHelp(false);
            cout << "\n############################ !!! ###########################\n";
            return 0;
        }
    }
    
    MbRandom myRandom;
    myRandom.setSeed(s1, s2);
	
    if(treeNodePrior == 9){
        Model myModel(&myRandom, calibFN, treeNodePrior, rho, runPrior);
        Mcmc mcmc(&myRandom, &myModel, numCycles, printFreq, sampleFreq, outName, writeDataFile, modUpdatePs, printOrigin, printAttach);
        return 0;
    }
    else {
        // set up the model
        Alignment myAlignment( dataFileName );
        cout << "Reading data from file -- " << dataFileName << endl;
        myAlignment.compress();
        if(printalign)
            myAlignment.print(std::cout);
        cout << "   Number of Site Patterns = " << myAlignment.getNumPatterns() << endl;
        string treeStr = getLineFromFile(treeFileName, 1);

        Model myModel(&myRandom, &myAlignment, treeStr, priorMean, rateSh, rateSc,
                      hyperSh, hyperSc, userBLs, moveAllN, offmove, rndNdMv, calibFN,
                      treeNodePrior, netDiv, relDeath, ssbdPrS, fixclokrt, rootfix, softbnd, calibHyP,
                      dpmExpHyp, dpmEHPPrM, gammaExpHP, modelType, fixModelPs, indHP, tipDateFN, fixTest, numIntervals,
					  originMax, runPrior);
        if(doAbsRts)
            myModel.setEstAbsRates(true);
        if(runPrior)
            myModel.setRunUnderPrior(true);
        if(justTree){
            myModel.writeUnifTreetoFile();
            return 0;
        }
        Mcmc mcmc(&myRandom, &myModel, numCycles, printFreq, sampleFreq, outName, writeDataFile, modUpdatePs, printOrigin, printAttach); //RW:
    }
	
    return 0;
}

