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
#include "Parameter_origin.h"
#include "Parameter_rate.h"
#include "Parameter_tree.h"
#include "Parameter_shape.h"
#include "Parameter_speciaton.h"
#include "Parameter_treescale.h"
#include "util.h"

#include <iomanip>
#include <iostream>
#include <ctime>

#include <time.h>

using namespace std;

Mcmc::Mcmc(MbRandom *rp, Model *mp, int nc, int pf, int sf, string ofp, bool wdf, bool modUpP, bool po) {//RW:

	ranPtr          = rp;
	modelPtr        = mp;
	numCycles       = nc;
	printFrequency  = pf;
	sampleFrequency = sf;
	fileNamePref    = ofp;
	writeInfoFile   = wdf;
	printratef		= false;
	modUpdateProbs  = modUpP;
    printOrigin     = po;
	runChain();
}

void Mcmc::runChain(void) {
	
	string pFile = fileNamePref + ".p"; // parameter file name
	string figTFile = fileNamePref + ".ant.tre"; // write to a file with the nodes colored by their rate classes
	string dFile = fileNamePref + ".info.out";
	string ndFile = fileNamePref + ".nodes.out"; // info about nodes
	string mtxFile = fileNamePref + ".rates.out"; // info about nodes
	ofstream pOut(pFile.c_str(), ios::out);
	ofstream fTOut(figTFile.c_str(), ios::out);
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
			sampleChain(n, pOut, fTOut, nOut, oldLnLikelihood);
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
	Tree *t = modelPtr->getActiveTree();
	int timeEnd = time(NULL);
	cout << "   Markov chain completed in " << (static_cast<float>(timeEnd - timeSt)) << " seconds" << endl;
	pOut.close();
	fTOut.close();
	dOut.close();
	nOut.close();
	mxOut.close();
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
            if(printOrigin)//RW
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
        if(printOrigin)//RW
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






// end
