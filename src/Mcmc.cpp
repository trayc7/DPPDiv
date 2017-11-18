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

#include "MbRandom.h"
#include "Mcmc.h"
#include "Model.h"
#include "Parameter.h"
#include "Parameter_basefreq.h"
#include "Parameter_exchangeability.h"
#include "Parameter_expcalib.h"
#include "Parameter_fossilgraph.h"
#include "Parameter_fossilrangegraph.h"
#include "Parameter_fossilrangegraphskyline.h"
#include "Parameter_origin.h"
#include "Parameter_rate.h"
#include "Parameter_tree.h"
#include "Parameter_shape.h"
#include "Parameter_speciaton.h"
#include "Parameter_speciationskyline.h"
#include "Parameter_treescale.h"
#include "util.h"

#include <iomanip>
#include <iostream>
#include <ctime>

#include <time.h>

using namespace std;

Mcmc::Mcmc(MbRandom *rp, Model *mp, int nc, int pf, int sf, string ofp, bool wdf, bool modUpP, bool po, bool pfat) {

	ranPtr          = rp;
	modelPtr        = mp;
	numCycles       = nc;
	printFrequency  = pf;
	sampleFrequency = sf;
	fileNamePref    = ofp;
	writeInfoFile   = wdf;//rw: check where this is switched on
	printratef		= false;
	modUpdateProbs  = modUpP;
    printOrigin     = po;
    printAttach     = pfat;
    treeTimePr      = mp->getTreeTimePriorNum();
    if(treeTimePr < 9)
        runChain();
    else if(treeTimePr == 9)
        runFOFBDChain();
    else if(treeTimePr == 10 || treeTimePr == 11)
        runFRGFBDChain();
}

void Mcmc::runChain(void) {
	
	int expMode = modelPtr->getFBDSExperimentalMode();
	
	string pFile = fileNamePref + ".p"; // parameter file name
	string figTFile = fileNamePref + ".ant.tre"; // write to a file with the nodes colored by their rate classes
	string dFile = fileNamePref + ".info.out";
	string ndFile = fileNamePref + ".nodes.out"; // info about nodes
	string mtxFile = fileNamePref + ".rates.out"; // info about nodes
	ofstream pOut;
	ofstream fTOut;
	if(expMode == 0){
		pOut.open(pFile.c_str(), ios::out);
		fTOut.open(figTFile.c_str(), ios::out);
	}
	ofstream nOut(ndFile.c_str(), ios::out);
	ofstream mxOut;
	ofstream dOut;
	
	if(writeInfoFile)
		dOut.open(dFile.c_str(), ios::out);
	if(printratef)
		mxOut.open(mtxFile.c_str(), ios::out);
		
	double oldLnLikelihood = modelPtr->lnLikelihood();
	
	// verbose logging
	if(writeInfoFile){
		dOut << "Running MCMC with:\n";
		dOut << "   Starting seeds = { " << modelPtr->getStartingSeed1() << " , " << modelPtr->getStartingSeed2() << " } \n";
		dOut << "   # Gens = " << numCycles << "\n";
		dOut << "   Prior mean # groups = " << modelPtr->getPriorMeanV() << "\n";
		dOut << "   lnL = " << oldLnLikelihood << "\n";
		printAllModelParams(dOut);
	}
	
	int timeSt = (int)time(NULL);
	bool testLnL = false;
	int modifyUProbsGen = (int)numCycles * 0.5;
	for (int n=1; n<=numCycles; n++){
		if(modUpdateProbs && n == modifyUProbsGen)
			modelPtr->setUpdateProbabilities(false);

		modelPtr->switchActiveParm();
		Parameter *parm = modelPtr->pickParmToUpdate();
        
		double prevlnl = oldLnLikelihood;
		double lnPriorProposalRatio = parm->update(oldLnLikelihood);
		
		double newLnLikelihood = modelPtr->getMyCurrLnL(); 
		double lnLikelihoodRatio = newLnLikelihood - oldLnLikelihood;
		
		double lnR = lnLikelihoodRatio + lnPriorProposalRatio;
		double r = safeExponentiation(lnR);
		
		bool isAccepted = false;
		if ( ranPtr->uniformRv() < r )
			isAccepted = true;
		
		if ( n % printFrequency == 0 || n == 1){
			cout << setw(6) << n << " -- " << fixed << setprecision(3) << prevlnl << " -> " << newLnLikelihood << endl;
			if(writeInfoFile){
				dOut << setw(6) << n << " -- " << fixed << setprecision(3) << prevlnl << " -> " << newLnLikelihood << endl;
				dOut << n << " -- " << parm->writeParam();
			}
		}
		
		if (isAccepted == true){
			oldLnLikelihood = newLnLikelihood;
			modelPtr->updateAccepted();
		}
		else{
			modelPtr->updateRejected();
			
			// TAH : without this, the lnls were not right for some moves following rejected ones 
			Tree *t = modelPtr->getActiveTree();
			Treescale *ts = modelPtr->getActiveTreeScale();
			t->setTreeScale(ts->getScaleValue());
			t->flipAllCls();
			t->flipAllTis();
			t->upDateAllCls();
			t->upDateAllTis();
			modelPtr->upDateRateMatrix();
			modelPtr->setTiProb();
		}
		
		if(n < 100){
			Tree *t = modelPtr->getActiveTree(); 
			t->setNodeRateValues();
		}
		
		// sample chain
		if ( n % sampleFrequency == 0 || n == 1){
			if(expMode == 0)
				sampleChain(n, pOut, fTOut, nOut, oldLnLikelihood);
			else if(expMode > 0)
				sampleFBDExpChain(n, nOut, oldLnLikelihood);
			//sampleRtsFChain(n, mxOut);
		}
		
		//Logger & logger = Logger::getInstance();
		//std::ostream &o = logger.debugStream();
		//parm->print(o);
		//parm->print(std::cerr);
#		if 0
		if( n % 20 == 0 || n == 1)
			modelPtr->getActiveTree()->verifyTreeDebug(n, parm->getName());
#		endif
	
		if(testLnL){
			Tree *t = modelPtr->getActiveTree(); 
			t->flipAllCls();
			t->flipAllTis();
			t->upDateAllCls();
			t->upDateAllTis();
			modelPtr->upDateRateMatrix();
			modelPtr->setTiProb();
			double chklnl = modelPtr->lnLikelihood();
			cout << "   " << n << "  " << parm->getName() << " -->  OL" << prevlnl << "   ---   NL" << newLnLikelihood << "  (" << chklnl << ", " << oldLnLikelihood << ")\n\n";
		}
	}
	int timeEnd = (int)time(NULL);
	cout << "   Markov chain completed in " << (static_cast<float>(timeEnd - timeSt)) << " seconds" << endl;
	pOut.close();
	fTOut.close();
	dOut.close();
	nOut.close();
	mxOut.close();
}

void Mcmc::runFOFBDChain() {
    
    string dFile = fileNamePref + ".info.out"; // I don't think we use this
    string occFile = fileNamePref + ".occ.out"; // info about occurrences, speciation parameters
    ofstream oOut(occFile.c_str(), ios::out);
    ofstream mxOut; //?
    ofstream dOut; //?
	
	int numMoves = modelPtr->getTotalUpdateWeights();
    
    if(writeInfoFile)
        dOut.open(dFile.c_str(), ios::out);
    
    double oldLnLikelihood = modelPtr->getActiveFossilGraph()->getActiveFossilGraphProb();
    
    // verbose logging
    if(writeInfoFile){
        dOut << "Running MCMC with:\n";
        dOut << "   Starting seeds = { " << modelPtr->getStartingSeed1() << " , " << modelPtr->getStartingSeed2() << " } \n";
        dOut << "   # Gens = " << numCycles << "\n";
        dOut << "   lnL = " << oldLnLikelihood << "\n";
        //printAllModelParams(dOut);
    }
    
    int timeSt = (int)time(NULL);
    int modifyUProbsGen = (int)numCycles * 0.5;
    for (int n=1; n<=numCycles; n++){
        if(modUpdateProbs && n == modifyUProbsGen)
            modelPtr->setUpdateProbabilities(false);
        double prevlnl = oldLnLikelihood;
		double newLnLikelihood = 0.0;
		
		for(int it=0; it<numMoves; it++){
			modelPtr->switchActiveParm();
			Parameter *parm = modelPtr->pickParmToUpdate();
			
			prevlnl = oldLnLikelihood;
			newLnLikelihood = parm->update(oldLnLikelihood); //all of the moves for FOFBD return the new likelihood for reporting
			bool isAccepted = true;
			if (isAccepted == true){
				oldLnLikelihood = newLnLikelihood;
				modelPtr->updateAccepted();
			}
        }
        
        if ( n % printFrequency == 0 || n == 1){
            cout << setw(6) << n << " -- " << fixed << setprecision(3) << prevlnl << " -> " << newLnLikelihood << endl;
//            if(writeInfoFile){
//                dOut << setw(6) << n << " -- " << fixed << setprecision(3) << prevlnl << " -> " << newLnLikelihood << endl;
//                dOut << n << " -- " << parm->writeParam();
//            }
        }
        
        // sample chain
        if ( n % sampleFrequency == 0 || n == 1){
            sampleChain(n, oOut, oldLnLikelihood);
        }
        
    }
    
    int timeEnd = (int)time(NULL);
    cout << "   Markov chain completed in " << (static_cast<float>(timeEnd - timeSt)) << " seconds" << endl;
    dOut.close();
    oOut.close();
    mxOut.close();
    
}

void Mcmc::runFRGFBDChain() {
    
    string dFile = fileNamePref + ".info.out"; // I don't think we use this
    string frFile = fileNamePref + ".FR.out"; // info about occurrences, speciation parameters
    
    ofstream frOut(frFile.c_str(), ios::out);
    //ofstream mxOut;
    ofstream dOut;
    
    int numMoves = modelPtr->getTotalUpdateWeights();
    
    if(writeInfoFile)
        dOut.open(dFile.c_str(), ios::out);
    
    double oldLnLikelihood = 0.0;
    
    if(treeTimePr == 10)
        oldLnLikelihood = modelPtr->getActiveFossilRangeGraph()->getActiveFossilRangeGraphProb();
    else if(treeTimePr == 11)
        oldLnLikelihood  = modelPtr->getActiveFossilRangeGraphSkyline()->getActiveFossilRangeGraphSkylineProb();
        
    // verbose logging
    if(writeInfoFile){
        dOut << "Running MCMC with:\n";
        dOut << "   Starting seeds = { " << modelPtr->getStartingSeed1() << " , " << modelPtr->getStartingSeed2() << " } \n";
        dOut << "   # Gens = " << numCycles << "\n";
        dOut << "   lnL = " << oldLnLikelihood << "\n";
        //printAllModelParams(dOut);
    }
    
    int timeSt = (int)time(NULL);
    int modifyUProbsGen = (int)numCycles * 0.5;
    
    for (int n=1; n<=numCycles; n++){
        
        if(modUpdateProbs && n == modifyUProbsGen) // do we need this here?
            modelPtr->setUpdateProbabilities(false);
        
        double prevlnl = oldLnLikelihood;
        double newLnLikelihood = 0.0;
        
        for(int it=0; it<numMoves; it++){
            modelPtr->switchActiveParm();
            Parameter *parm = modelPtr->pickParmToUpdate();

            prevlnl = oldLnLikelihood;
            newLnLikelihood = parm->update(oldLnLikelihood); // all of the moves for FRGFBD return the new likelihood for reporting
            bool isAccepted = true;
            if (isAccepted == true){
                oldLnLikelihood = newLnLikelihood;
                modelPtr->updateAccepted();
            }
        }

        if ( n % printFrequency == 0 || n == 1){
            cout << setw(6) << n << " -- " << fixed << setprecision(3) << prevlnl << " -> " << newLnLikelihood << endl;
            //if(writeInfoFile){
                //dOut << setw(6) << n << " -- " << fixed << setprecision(3) << prevlnl << " -> " << newLnLikelihood << endl;
                //dOut << n << " -- " << parm->writeParam();
            //}
        }

        // sample chain
        if ( n % sampleFrequency == 0 || n == 1){
            if(treeTimePr == 10)
                sampleChainFR(n, frOut, oldLnLikelihood);
            if(treeTimePr == 11)
                sampleChainFRSkyline(n, frOut, oldLnLikelihood);
        }
    }
    
    int timeEnd = (int)time(NULL);
    cout << "   Markov chain completed in " << (static_cast<float>(timeEnd - timeSt)) << " seconds" << endl;
    dOut.close();
    frOut.close();
    //mxOut.close();
    
}

double Mcmc::safeExponentiation(double lnX) {

	if (lnX < -300.0)
		return 0.0;
	else if (lnX > 0.0)
		return 1.0;
	else
		return exp(lnX);
}

void Mcmc::sampleChain(int gen, ofstream &paraOut, ofstream &figTOut, 
					   ofstream &nodeOut, double lnl) {

	Basefreq *f = modelPtr->getActiveBasefreq();
	Exchangeability *e = modelPtr->getActiveExchangeability();
	NodeRate *nr = modelPtr->getActiveNodeRate();
	Tree *t = modelPtr->getActiveTree();
	Shape *sh = modelPtr->getActiveShape();
	Speciation *sp = modelPtr->getActiveSpeciation();
	Treescale *ts = modelPtr->getActiveTreeScale();
    OriginTime *ot = modelPtr->getActiveOriginTime();
	ExpCalib *hpex;
	sp->setAllBDFossParams();
	bool expHPCal = modelPtr->getExponCalibHyperParm();
	bool dpmHPCal = modelPtr->getExponDPMCalibHyperParm();
	int treePr = modelPtr->getTreeTimePriorNum();
	if(expHPCal)
		hpex = modelPtr->getActiveExpCalib();
	
	if(gen == 1){
		paraOut << "Gen\tlnLikelihood\tf(A)\tf(C)\tf(G)\tf(T)";
//		paraOut << "\tr(AC)\tr(AG)\tr(AT)\tr(CG)\tr(CT)\tr(GT)\tshape\tave rate\tnum rate groups\tconc param\n";
		paraOut << "\tr(AC)\tr(AG)\tr(AT)\tr(CG)\tr(CT)\tr(GT)\tshape\n";
		figTOut << "#NEXUS\nbegin trees;\n";
		nodeOut << "Gen\tlnL";
		nodeOut << "\tNetDiv(b-d)\tRelativeDeath(d/b)";
		if(treePr > 3)
			nodeOut << "\tFBD.psi\tFBD.rho";
		if(treePr == 4)
			nodeOut << "\tbdss.torig";
		if(treePr > 5)
			nodeOut << "\tFBD.lambda\tFBD.mu\tFBD.prsp";
        if(treePr == 8){
            if(printOrigin)
            nodeOut << "\tFBD.OriginTime";
        }
		nodeOut << "\tPr(speciation)\tave.subrate\tnum.DPMgroups\tDPM.conc";
		if(expHPCal){
			if(dpmHPCal)
				nodeOut << "\texpHP.dpmConP\texpHP.dpmNumLs";
			else
				nodeOut << "\texpHP.lambda1\texpHP.lambda2\texpHP.epsilon";
		}
			
		nodeOut << t->getNodeInfoNames();
		if(expHPCal)
			nodeOut << t->getCalNodeInfoNames();
		
		if(treePr > 5){
			nodeOut << t->getCalBDSSNodeInfoParamNames();
		}
		if(treePr >= 7){
//			nodeOut << t->getCalBDSSNodeInfoIndicatorNames();
			nodeOut << "\tnum.tip_fossils";
		}
		
		nodeOut << "\n";
	}
	// then print stuff
	paraOut << gen << "\t" << lnl;
	for(int i=0; i<f->getNumStates(); i++)
		paraOut << "\t" << f->getFreq(i);
	for(int i=0; i<6; i++)
		paraOut << "\t" << e->getRate(i);
	paraOut << "\t" << sh->getAlphaSh();
//	paraOut << "\t" << nr->getAverageRate();
//	paraOut << "\t" << nr->getNumRateGroups();
//	paraOut << "\t" << nr->getConcenParam();
	paraOut << "\n";
		
	figTOut << "  tree t" << gen << " = ";
	figTOut << t->getFigTreeDescription() << "\n";
	
	nodeOut << gen << "\t" << lnl;
	nodeOut << "\t" << sp->getNetDiversification();
	nodeOut << "\t" << sp->getRelativeDeath();
	if(treePr > 3){
		nodeOut << "\t" << sp->getBDSSFossilSampRatePsi();
		nodeOut << "\t" << sp->getBDSSSppSampRateRho();
	}
	if(treePr == 4)
		nodeOut << "\t" << ts->getTreeOriginTime();
	if(treePr > 5){
		nodeOut << "\t" << sp->getBDSSSpeciationRateLambda();
		nodeOut << "\t" << sp->getBDSSExtinctionRateMu();
		nodeOut << "\t" << sp->getBDSSFossilSampProbS();
	}
    if(treePr == 8){
        if(printOrigin)
        nodeOut << "\t" << ot->getOriginTime();
    }
    nodeOut << "\t" << t->getTreeSpeciationProbability();
	nodeOut << "\t" << nr->getAverageRate();
	nodeOut << "\t" << nr->getNumRateGroups();
	nodeOut << "\t" << nr->getConcenParam();
	if(expHPCal){
		if(dpmHPCal){
			nodeOut << "\t" << hpex->getDPMExpHPConcentParam();
			nodeOut << "\t" << hpex->getNumLambdaTables();
		}
		else{
			nodeOut << "\t" << hpex->getCurMajorityLambda();
			nodeOut << "\t" << hpex->getCurOutlieLambda();
			nodeOut << "\t" << hpex->getEpsilonValue();
		}
	}
	nodeOut << t->getNodeInfoList();
	if(expHPCal)
		nodeOut << t->getCalNodeInfoList();

	if(treePr > 5){
		nodeOut << t->getCalBDSSNodeInfoParamList();
	}
	if(treePr >= 7){
//		nodeOut << t->getCalBDSSNodeInfoIndicatorList();
		nodeOut << "\t" << t->getSumIndicatorV();
	}
	nodeOut << "\n";
	
	if(gen == numCycles){
		figTOut << "end;\n";
		figTOut << "\nbegin figtree;\n";
		figTOut << "    set appearance.branchColorAttribute=\"rate_cat\";\n";
		figTOut << "    set appearance.branchLineWidth=2.0;\n";
		figTOut << "    set scaleBar.isShown=false;\n";
		figTOut << "end;\n";
	}
}

void Mcmc::sampleFBDExpChain(int gen, std::ofstream &nodeOut, double lnl) {

	NodeRate *nr = modelPtr->getActiveNodeRate();
	Tree *t = modelPtr->getActiveTree();
	Speciation *sp = modelPtr->getActiveSpeciation();
	Treescale *ts = modelPtr->getActiveTreeScale();
    OriginTime *ot = modelPtr->getActiveOriginTime();
	sp->setAllBDFossParams();
	int treePr = modelPtr->getTreeTimePriorNum();
	
	if(gen == 1){
		nodeOut << "Gen\tlnL";
		nodeOut << "\tNetDiv(b-d)\tRelativeDeath(d/b)";
		if(treePr > 3)
			nodeOut << "\tFBD.psi\tFBD.rho";
		if(treePr == 4)
			nodeOut << "\tbdss.torig";
		if(treePr > 5)
			nodeOut << "\tFBD.lambda\tFBD.mu\tFBD.prsp";
        if(treePr == 8){
            if(printOrigin)
            nodeOut << "\tFBD.OriginTime";
        }
		nodeOut << "\tPr(speciation)\tave.subrate\tnum.DPMgroups\tDPM.conc";
		
		nodeOut << t->getNodeInfoNames();
		
		if(treePr > 5){
			nodeOut << t->getCalBDSSNodeInfoParamNames();
		}
		if(treePr >= 7){
//			nodeOut << t->getCalBDSSNodeInfoIndicatorNames();
			nodeOut << "\tnum.tip_fossils";
		}
		
		nodeOut << "\n";
	}
	// then print stuff
	
	nodeOut << gen << "\t" << lnl;
	nodeOut << "\t" << sp->getNetDiversification();
	nodeOut << "\t" << sp->getRelativeDeath();
	if(treePr > 3){
		nodeOut << "\t" << sp->getBDSSFossilSampRatePsi();
		nodeOut << "\t" << sp->getBDSSSppSampRateRho();
	}
	if(treePr == 4)
		nodeOut << "\t" << ts->getTreeOriginTime();
	if(treePr > 5){
		nodeOut << "\t" << sp->getBDSSSpeciationRateLambda();
		nodeOut << "\t" << sp->getBDSSExtinctionRateMu();
		nodeOut << "\t" << sp->getBDSSFossilSampProbS();
	}
    if(treePr == 8){
        if(printOrigin)
        nodeOut << "\t" << ot->getOriginTime();
    }
    nodeOut << "\t" << t->getTreeSpeciationProbability();
	nodeOut << "\t" << nr->getAverageRate();
	nodeOut << "\t" << nr->getNumRateGroups();
	nodeOut << "\t" << nr->getConcenParam();
	nodeOut << t->getNodeInfoList();


	if(treePr > 5){
		nodeOut << t->getCalBDSSNodeInfoParamList();
	}
	if(treePr >= 7){
//		nodeOut << t->getCalBDSSNodeInfoIndicatorList();
		nodeOut << "\t" << t->getSumIndicatorV();
	}
	nodeOut << "\n";
	
}

void Mcmc::printAllModelParams(ofstream &dOut){
	
	dOut << "\n--------------------------------------------------\n";
	dOut << "Initial: \n";
	dOut << modelPtr->getActiveBasefreq()->writeParam();
	dOut << modelPtr->getActiveExchangeability()->writeParam();
	dOut << modelPtr->getActiveShape()->writeParam();
	dOut << modelPtr->getActiveNodeRate()->writeParam();
	dOut << modelPtr->getActiveTree()->writeParam();
	dOut << "--------------------------------------------------\n\n";
}

void Mcmc::writeCalibrationTree(){
	
	modelPtr->setNodeRateGrpIndxs();
	Tree *t = modelPtr->getActiveTree();
	if(t->getIsCalibratedTree()){
		string ctfn = fileNamePref + ".CALIB.tre";
		ofstream cto(ctfn.c_str(), ios::out);
		cto << "#NEXUS\nbegin trees;\n";
		cto << "    tree calib_init = ";
		cto << t->getCalibInitialTree() << "\nend;\n";
		cto << "\nbegin figtree;\n";
		cto << "    set appearance.branchColorAttribute=\"User selection\";\n";
		cto << "    set appearance.branchLineWidth=3.0;\n";
		cto << "    set nodeBars.isShown=true;\n";
		cto << "    set nodeBars.barWidth=27.0;\n";
		cto << "    set nodeLabels.isShown=true;\n";
		cto << "    set nodeLabels.displayAttribute=\"Node ages\";\n";
		cto << "    set nodeLabels.fontSize=18;\n";
		cto << "    set nodeLabels.fontStyle=2;\n";
		cto << "    set rectilinearLayout.rootLength=0;\n";
		cto << "    set scaleAxis.isShown=true;\n";
		cto << "    set scaleAxis.reverseAxis=true;\n";
		cto << "    set tipLabels.fontSize=20;\n";
		cto << "    set scaleBar.isShown=false;\n";
		cto << "end;\n";
	}
}

void Mcmc::sampleRtsFChain(int gen, std::ofstream &rOut){
	
	NodeRate *nr = modelPtr->getActiveNodeRate();
	Tree *t = modelPtr->getActiveTree();
	if(gen == 1){
		rOut << "Gen\tDPM.conc";
		rOut << t->getDownPNodeInfoNames();
		rOut << "\n";
	}
	
	rOut << gen << "\t" << nr->getConcenParam();
	rOut << t->getDownPNodeInfoList();
	rOut << "\n";
}


void Mcmc::sampleChain(int gen, ofstream &occOut, double lnl) {
    
    FossilGraph *fg = modelPtr->getActiveFossilGraph();
    Speciation *sp = modelPtr->getActiveSpeciation();
    OriginTime *ot = modelPtr->getActiveOriginTime();

    sp->setAllBDFossParams();

    //int treePr = modelPtr->getTreeTimePriorNum(); // probably required for options >9
    
    if(gen == 1){
        
        occOut << "Gen\tlnL";
        occOut << "\tNetDiv(b-d)\tRelativeDeath(d/b)";
        occOut << "\tFBD.psi\tFBD.rho";
        occOut << "\tFBD.lambda\tFBD.mu\tFBD.prsp";
        if(printOrigin)
            occOut << "\tFBD.OriginTime";
        if(printAttach)
            occOut << fg->getOccInfoParamNames(); // cf getNodeInfoNames()
        occOut << "\tnum.tip_fossils";
        occOut << "\n";
    }
    
    // then print stuff
    occOut << gen << "\t" << lnl;
    occOut << "\t" << sp->getNetDiversification();
    occOut << "\t" << sp->getRelativeDeath();
    occOut << "\t" << sp->getBDSSFossilSampRatePsi();
    occOut << "\t" << sp->getBDSSSppSampRateRho();
    occOut << "\t" << sp->getBDSSSpeciationRateLambda();
    occOut << "\t" << sp->getBDSSExtinctionRateMu();
    occOut << "\t" << sp->getBDSSFossilSampProbS();
    if(printOrigin)
        occOut << "\t" << ot->getOriginTime();
    if(printAttach)
        occOut << fg->getOccInfoParamList(); // cf getNodeInfoList()
    occOut << "\t" << fg->getSumIndicatorFG();
    occOut << "\n";
    
}

void Mcmc::sampleChainFR(int gen, ofstream &frOut, double lnl) {
    
    FossilRangeGraph *frg = modelPtr->getActiveFossilRangeGraph();
    Speciation *sp = modelPtr->getActiveSpeciation();
    
    sp->setAllBDFossParams();
    
    if(gen == 1){
        
        frOut << "Gen\tlnL";
        frOut << "\tNetDiv(b-d)\tRelativeDeath(d/b)";
        frOut << "\tFBD.psi\tFBD.rho";
        frOut << "\tFBD.lambda\tFBD.mu\tFBD.prsp";
        if(printOrigin)
            frOut << "\tFBD.OriginTime";
        if(printAttach)
            frOut << frg->getFossilRangeInfoParamNames();
        frOut << "\tnum.fossils(k)";
        frOut << "\n";
    }
    
    // then print stuff
    //frOut << gen << "\t" << fixed << setprecision(3) <<  lnl;
    frOut << gen << "\t" <<  lnl;
    frOut << "\t" << sp->getNetDiversification();
    frOut << "\t" << sp->getRelativeDeath();
    frOut << "\t" << sp->getBDSSFossilSampRatePsi();
    frOut << "\t" << sp->getBDSSSppSampRateRho();
    frOut << "\t" << sp->getBDSSSpeciationRateLambda();
    frOut << "\t" << sp->getBDSSExtinctionRateMu();
    frOut << "\t" << sp->getBDSSFossilSampProbS();
    if(printOrigin)
        frOut << "\t" << frg->getFossilRangeGraphOriginTime();
    if(printAttach)
        frOut << frg->getFossilRangeInfoParamList();
    //occOut << "\t" << fg->getSumIndicatorFG();
    frOut << "\t" << frg->getNumFossils();
    frOut << "\n";
    
}

void Mcmc::sampleChainFRSkyline(int gen, ofstream &frOut, double lnl) {
    
    FossilRangeGraphSkyline *frg = modelPtr->getActiveFossilRangeGraphSkyline();
    SpeciationSkyline *sp = modelPtr->getActiveSpeciationSkyline();
    
    //sp->setAllBDFossParams();
    
    int numIntervals = frg->getNumIntervals();
    
    if(gen == 1){
        
        frOut << "Gen\tlnL";
        
        //for(int i=0; i < numIntervals; i++){
        //    frOut << "\tbirth[" << i << "]\tdeath[" << i << "]\tpsi[" << i << "]";
        //}
        
        int intID = 0;
        for(int i = numIntervals - 1; i >= 0; i--){
            frOut << "\tlambda[" << intID << "]";
            intID += 1;
        }
        intID = 0;
        for(int i = numIntervals - 1; i >= 0; i--){
            frOut << "\tmu[" << intID << "]";
            intID += 1;
        }
        intID = 0;
        for(int i = numIntervals - 1; i >= 0; i--){
            frOut << "\tpsi[" << intID << "]";
            intID += 1;
        }
        
        frOut << "\tFBD.rho[0]\t";
        if(printOrigin)
          frOut << "\tFBD.OriginTime";
        if(printAttach)
            frOut << frg->getFossilRangeSkylineInfoParamNames();
        frOut << "\tnum.fossils(k)";
        frOut << "\n";
    }
    frOut << gen << "\t" <<  lnl;
    
    //for(int i=0; i < numIntervals; i++){
    //frOut << "\t" << sp->getSpeciationRates()[i] << "\t" << sp->getExtinctionRates()[i] << "\t" << sp->getFossilSampRates()[i];
    //}
    
    for(int i = numIntervals - 1; i >= 0; i--){
        frOut << "\t" << sp->getSpeciationRates()[i];
    }
    for(int i = numIntervals - 1; i >= 0; i--){
        frOut << "\t" << sp->getExtinctionRates()[i];
    }
    for(int i = numIntervals - 1; i >= 0; i--){
        frOut << "\t" << sp->getFossilSampRates()[i];
    }
    frOut << "\t" << sp->getSppSampRateRho()[0];
    if(printOrigin)
        frOut << "\t" << frg->getFossilRangeGraphSkylineOriginTime();
    if(printAttach)
        frOut << frg->getFossilRangeSkylineInfoParamList();
    frOut << "\t" << frg->getNumFossils();
    frOut << "\n";
}

// end
