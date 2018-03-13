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

#include "cpuspec.h"
#include "Alignment.h"
#include "MbRandom.h"
#include "MbTransitionMatrix.h"
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
#include "Parameter_shape.h"
#include "Parameter_speciaton.h"
#include "Parameter_speciationskyline.h"
#include "Parameter_tree.h"
#include "Parameter_cphyperp.h"
#include "Parameter_treescale.h"
#include "Calibration.h"
#include "util.h"
#include <string>
#include <vector>
#include <fstream>
#include <cstring>

using namespace std;

Model::Model(MbRandom *rp, Alignment *ap, string ts, double pm, double ra, double rb, 
			 double hal, double hbe, bool ubl, bool alnm, int offmv, bool rndNo, 
			 string clfn, int nodpr, double bdr, double bda, double bds, double fxclkrt, bool roofix,
			 bool sfb, bool ehpc, bool dphpc, int dphpng, bool gamhp, int rmod, bool fxmod,
			 bool ihp, string tipdfn, bool fxtr, int expmo, bool igfoss, double rh, int bdp, bool fxPsi, double psi,
             int specPr, int psiPr, double bPrRate, double dPrRate, double pPrRate) {
	// remember pointers to important objects...
	ranPtr       = rp;
	alignmentPtr = ap;
	priorMeanN   = pm;
	ranPtr->getSeed(startS1, startS2);
	runUnderPrior = false;
	treeTimePrior = nodpr;
	myCurLnL = 0.0;
	lnLGood = false;
	calibfilen = clfn;
	tipDateFileN = tipdfn;
	fixRootHeight = roofix;
	zeroNodeTimeMove = false;
	exponCalibHyperParm = ehpc;
	exponDPMCalibHyperParm = dphpc;
	turnedOffMove = offmv;
	fixedClockRate = fxclkrt;
	fixSomeModParams = fxmod;
	fixTestRun = fxtr;
	estAbsRts = false;
    originMax = 100000.0;
    rho = rh; //RW: only active in experimental mode
    fbdPar = bdp;
    if(treeTimePrior == 8)
        conditionOnOrigin = true; // some redundancy for now
	double initRootH = 1.0;
	runIndCalHP = ihp;
	rHtY = 0.0;
	rHtO = 0.0;
	int tsPrDist = 1;
	rootNExpRate = -1.0;
	bool rtCalib = false;
    cout << "expmo" << expmo << endl;
    fbdsExperimentalMode = expmo;// experimental mode = 1: fix user specified branch lengths; = 2 fix user specified branch lengths and ignore the calibrations
    ignoreFossils = igfoss;
        
	if(calibfilen.empty() == false){
		initRootH = readCalibFile();
//        initRootH = 100.0; // TAH: DEBUG
//        initOT = 150.0;
		Calibration *rCal = getRootCalibration();
		if(rCal != NULL && fixRootHeight == false){
			rHtY = rCal->getYngTime();
			rHtO = rCal->getOldTime();
			tsPrDist = rCal->getPriorDistributionType();
			rtCalib = true;
			rootNExpRate =  rCal->getCalExponRate();
		}
		else if(rCal == NULL && fixRootHeight == false){
			rHtY = rHtY / 1.1;
			rHtO = -1.0;
			tsPrDist = 2;
		}
	}
    if(conditionOnOrigin)
        rHtO = initOT;
    else rHtO = originMax;
    
	cout << "\nStarting with seeds: { " << startS1 << " , " << startS2 << " } \n\n";
	
	// ...and initialize some important variables
    
	numGammaCats = 4;
	numPatterns  = alignmentPtr->getNumChar();
	
	cpfix = false;
	if(turnedOffMove == 5)
		cpfix = true;
	else if(turnedOffMove == 6)
		cpfix = true;
	if(rmod > 1)
		cpfix = true;
    
	int nn = 2*alignmentPtr->getNumTaxa()-1;
	if(pm > nn - 1){
		cerr << "ERROR: the prior on the mean number of tables cannot exceed the number of nodes in the tree!" << endl;
		exit(1);
	}
	Cphyperp *conp = new Cphyperp(ranPtr, this, hal, hbe, nn, priorMeanN, cpfix);
	ExpCalib *excal = new ExpCalib(ranPtr, this, dphpc, dphpng, initRootH, gamhp, runIndCalHP);
	NodeRate *nr = new NodeRate(ranPtr, this, nn, ra, rb, conp->getCurrentCP(), fxclkrt, rmod);
	for (int i=0; i<2; i++){ 
		parms[i].push_back( new Basefreq(ranPtr, this, 4, fxmod) );					// base frequency parameter
		parms[i].push_back( new Exchangeability(ranPtr, this) );				// rate parameters of the GTR model
		parms[i].push_back( new Shape(ranPtr, this, numGammaCats, 2.0, fxmod) );		// gamma shape parameter for rate variation across sites
		parms[i].push_back( new Tree(ranPtr, this, alignmentPtr, ts, ubl, alnm, rndNo, 
									 calibrs, initRootH, initOT, sfb, ehpc, excal, tipDates, fbdsExperimentalMode) );    // rooted phylogenetic tree
		parms[i].push_back( nr );												// restaurant containing node rates
		parms[i].push_back( conp );												// hyper prior on DPP concentration parameter
		parms[i].push_back( new Treescale(ranPtr, this, initRootH, rHtY, rHtO, tsPrDist, rtCalib, ehpc) ); // the tree scale prior
		parms[i].push_back( new Speciation(ranPtr, this, bdr, bda, bds, initRootH, rho, fbdPar, fxPsi, psi, specPr, psiPr, bPrRate, dPrRate, pPrRate) );												// hyper prior on diversification for cBDP speciation
		parms[i].push_back( excal );											// hyper prior exponential node calibration parameters
        parms[i].push_back( new OriginTime(ranPtr, this, initOT, rHtY, originMax) ); // the origin time parameters
	}
	numParms = (int)parms[0].size();
	activeParm = 0;
	for (int i=0; i<numParms; i++)
		*parms[0][i] = *parms[1][i];

	for (int i=0; i<numParms; i++)
		parms[0][i]->print(std::cout);
	
	updateTreeRates();
	// initialize the probabilities for updating different parameters
	setUpdateProbabilities(true);	
	if(ehpc)
		excal->getAllExpHPCalibratedNodes();
		
	// allocate and initialize conditional likelihoods
	initializeConditionalLikelihoods();
	
	// allocate the transition probability matrices
	initializeTransitionProbabilityMatrices();
	
	// instantiate the transition probability calculator
	tiCalculator = new MbTransitionMatrix( getActiveExchangeability()->getRate(), getActiveBasefreq()->getFreq(), true );
	
	setTiProb();
	myCurLnL = lnLikelihood();
	cout << "lnL = " << myCurLnL << endl;

}

Model::Model(MbRandom *rp, std::string clfn, int nodpr, double rh, bool rnp, int bdp, bool fixFRG, bool estExt, bool fixInd, bool lnSurf, bool fxPsi, double psi, int compS,
             int specPr, int psiPr, double bPrRate, double dPrRate, double pPrRate, int expMode){
    
    ranPtr = rp;
    ranPtr->getSeed(startS1, startS2);
    treeTimePrior = nodpr;
    calibfilen = clfn;
	runUnderPrior = rnp;
    rho = rh;
    fbdPar = bdp;
    numFossils = 0;
    numLineages = 0;
    myCurLnL = 0.0;
    
    if(calibfilen.empty() == false){
        if(treeTimePrior == 9)
            readOccurrenceFile(); // --> this function will read the file, create a Calibration obj for each one, and initialize initOT
        else if (treeTimePrior == 10)
            readFossilRangeFile(); // --> this function will read the file, create a Calibration obj for each fossil range
    }

    cout << "\nStarting with seeds: { " << startS1 << " , " << startS2 << " } \n\n";
    
    if(treeTimePrior == 9){
        FossilGraph *fg = new FossilGraph(ranPtr, this, numFossils, initOT, calibrs, runUnderPrior);
        OriginTime *ot = new OriginTime(ranPtr, this, initOT, rHtY, originMax);
        Speciation *sp = new Speciation(ranPtr, this, -1.0, -1.0, -1.0, 100.0, rho, bdp, fxPsi, psi, specPr, psiPr, bPrRate, dPrRate, pPrRate);
        for (int i=0; i<2; i++){
            parms[i].push_back( ot );
            parms[i].push_back( sp );
            parms[i].push_back( fg );
            
        }
        numParms = (int)parms[0].size();
        activeParm = 0;
        for (int i=0; i<numParms; i++)
            *parms[0][i] = *parms[1][i];
        
        for (int i=0; i<numParms; i++)
            parms[0][i]->print(std::cout);
        
        updateProb.clear();
        updateProb.push_back(0.0); // 1 origin time
        updateProb.push_back(3.0); // 2 speciation
        updateProb.push_back(4.0); // 3 fossil graph
        double sum = 0.0;
        for (unsigned i=0; i<updateProb.size(); i++)
            sum += updateProb[i];
        for (unsigned i=0; i<updateProb.size(); i++)
            updateProb[i] /= sum;
        
        totalUpdateWeights = (int)sum;
        
        myCurLnL = this->getActiveFossilGraph()->getActiveFossilGraphProb();
        cout << "lnL = " << myCurLnL << endl;
    }
    
    if(treeTimePrior == 10){
        
        FossilRangeGraph *frg = new FossilRangeGraph(ranPtr, this, numFossils, numLineages, calibrs, runUnderPrior, fixFRG, estExt, fixInd, compS, expMode);
        Speciation *sp = new Speciation(ranPtr, this, -1.0, -1.0, -1.0, 100.0, rho, fbdPar, fxPsi, psi, specPr, psiPr, bPrRate, dPrRate, pPrRate);
        
        for (int i=0; i<2; i++){
            parms[i].push_back( sp );
            parms[i].push_back( frg );
        }
        numParms = (int)parms[0].size();
        activeParm = 0;
        for (int i=0; i<numParms; i++)
            *parms[0][i] = *parms[1][i];
        
        for (int i=0; i<numParms; i++)
            parms[0][i]->print(std::cout);
        
        updateProb.clear();
        if(expMode == 1)
            updateProb.push_back(0.0); // 1 speciation
        else
            updateProb.push_back(3.0); // 1 speciation
        if(fixFRG)
            updateProb.push_back(0.0); // 2 fossil graph
        else
            updateProb.push_back(4.0); // 2 fossil graph
        
        double sum = 0.0;
        for (unsigned i=0; i<updateProb.size(); i++)
            sum += updateProb[i];
        for (unsigned i=0; i<updateProb.size(); i++)
            updateProb[i] /= sum;
        
        totalUpdateWeights = (int)sum;
        
        if(lnSurf)
            this->getActiveFossilRangeGraph()->lnSurfaceGenerator(clfn);
        
        myCurLnL = this->getActiveFossilRangeGraph()->getActiveFossilRangeGraphProb();
        cout << "lnL = " << myCurLnL << endl;
    }
    
}

Model::Model(MbRandom *rp, std::string clfn, std::string intfn, std::string pafn, int nodpr, double rh, bool rnp, int bdp, bool fixFRG, bool estExt, int expMode, int fbdLk,
             int specPr, int psiPr, double bPrRate, double dPrRate, double pPrRate, bool proxy){
    
    ranPtr = rp;
    ranPtr->getSeed(startS1, startS2);
    calibfilen = clfn;
    intfilen = intfn;
    pafilen = pafn;
    treeTimePrior = nodpr;
    rho = rh;
    runUnderPrior = rnp;
    if(rho < 0.0 || rho > 1.0) {
        cerr << "ERROR: Extant species sampling (-rho) must be > 0 and < 1." << endl;
        exit(1);
    }
    fbdPar = bdp;
    numFossils = 0;
    numLineages = 0;
    userSpecifiedIntervals = 0; //**skyline note user specified intervals
    myCurLnL = 0.0;
    
    if(fbdLk == 1 || fbdLk == 2)
        readFossilRangeFile(); // this function will read the file and create a Calibration obj for each fossil range
    else if(fbdLk == 3)
        readPresenceAbsenceFile(); // this function will read the file and create a Calibration obj for each fossil range, inc presence / absence data
    readIntervalsFile(); // this function will read the file, create a Calibration obj for each interval
    
    cout << "\nStarting with seeds: { " << startS1 << " , " << startS2 << " } \n\n";
    
    FossilRangeGraphSkyline *frg = new FossilRangeGraphSkyline(ranPtr, this, numFossils, numLineages, calibrs, userSpecifiedIntervals, intervals, runUnderPrior, fixFRG, estExt, expMode, fbdLk);
    SpeciationSkyline *sp = new SpeciationSkyline(ranPtr, this, userSpecifiedIntervals, rho, specPr, psiPr, bPrRate, dPrRate, pPrRate, proxy);
    
    for (int i=0; i<2; i++){
        parms[i].push_back( sp );
        parms[i].push_back( frg );
    }
    numParms = (int)parms[0].size();
    activeParm = 0;
    for (int i=0; i<numParms; i++)
        *parms[0][i] = *parms[1][i];
    
    for (int i=0; i<numParms; i++)
        parms[0][i]->print(std::cout);
    
    updateProb.clear();
    
    if(expMode == 1)
        updateProb.push_back(0.0); // 1 speciation
    else
        updateProb.push_back(3.0); // 1 speciation
    if(fixFRG)
        updateProb.push_back(0.0); // 2 fossil graph
    else
        updateProb.push_back(4.0); // 2 fossil graph
    
    double sum = 0.0;
    for (unsigned i=0; i<updateProb.size(); i++)
        sum += updateProb[i];
    for (unsigned i=0; i<updateProb.size(); i++)
        updateProb[i] /= sum;
    
    totalUpdateWeights = (int)sum;

    this->getActiveFossilRangeGraphSkyline()->setAllIntervalConstants();
    myCurLnL = this->getActiveFossilRangeGraphSkyline()->getActiveFossilRangeGraphSkylineProb();
    cout << "lnL = " << myCurLnL << endl;
    
}

Model::~Model(void) {

    if(treeTimePrior < 9){
        delete [] cls;
        delete tiCalculator;
        for (int i=0; i<2; i++){
            delete [] clPtr[i];
            delete [] tis[i][0];
            delete [] tis[i];
        }
    }
}

Basefreq* Model::getActiveBasefreq(void) {

	for (int i=0; i<numParms; i++){
		Parameter *p = parms[activeParm][i];
		Basefreq *derivedPtr = dynamic_cast<Basefreq *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
	}
	return NULL;
}

Tree* Model::getActiveTree(void) {

	for (int i=0; i<numParms; i++){
		Parameter *p = parms[activeParm][i];
		Tree *derivedPtr = dynamic_cast<Tree *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
	}
	return NULL;
}

Treescale* Model::getActiveTreeScale(void) {
	
	for (int i=0; i<numParms; i++){
		Parameter *p = parms[activeParm][i];
		Treescale *derivedPtr = dynamic_cast<Treescale *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
	}
	return NULL;
}


Exchangeability* Model::getActiveExchangeability(void) {

	for (int i=0; i<numParms; i++){
		Parameter *p = parms[activeParm][i];
		Exchangeability *derivedPtr = dynamic_cast<Exchangeability *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
	}
	return NULL;
}

Shape* Model::getActiveShape(void) {

	for (int i=0; i<numParms; i++){
		Parameter *p = parms[activeParm][i];
		Shape *derivedPtr = dynamic_cast<Shape *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
	}
	return NULL;
}

NodeRate* Model::getActiveNodeRate(void) {

	for (int i=0; i<numParms; i++){
		Parameter *p = parms[activeParm][i];
		NodeRate *derivedPtr = dynamic_cast<NodeRate *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
	}
	return NULL;
}

Speciation* Model::getActiveSpeciation(void) {
	
	for (int i=0; i<numParms; i++){
		Parameter *p = parms[activeParm][i];
		Speciation *derivedPtr = dynamic_cast<Speciation *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
	}
	return NULL;
}

Cphyperp* Model::getActiveCphyperp(void) {
	
	for (int i=0; i<numParms; i++){
		Parameter *p = parms[activeParm][i];
		Cphyperp *derivedPtr = dynamic_cast<Cphyperp *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
	}
	return NULL;
}

ExpCalib* Model::getActiveExpCalib(void) {
	
	for (int i=0; i<numParms; i++){
		Parameter *p = parms[activeParm][i];
		ExpCalib *derivedPtr = dynamic_cast<ExpCalib *>(p);
		if ( derivedPtr != 0 )
			return derivedPtr;
	}
	return NULL;
}

OriginTime* Model::getActiveOriginTime(void) {
    
    for (int i=0; i<numParms; i++){
        Parameter *p = parms[activeParm][i];
        OriginTime *derivedPtr = dynamic_cast<OriginTime *>(p);
        if ( derivedPtr != 0 )
            return derivedPtr;
    }
    return NULL;
}

FossilGraph* Model::getActiveFossilGraph(void) {
    
    for (int i=0; i<numParms; i++){
        Parameter *p = parms[activeParm][i];
        FossilGraph *derivedPtr = dynamic_cast<FossilGraph *>(p);
        if ( derivedPtr != 0 )
            return derivedPtr;
    }
    return NULL;
}

FossilRangeGraph* Model::getActiveFossilRangeGraph(void) {
    
    for (int i=0; i<numParms; i++){
        Parameter *p = parms[activeParm][i];
        FossilRangeGraph *derivedPtr = dynamic_cast<FossilRangeGraph *>(p);
        if ( derivedPtr != 0 )
            return derivedPtr;
    }
    return NULL;
}

FossilRangeGraphSkyline* Model::getActiveFossilRangeGraphSkyline(void) {
    
    for (int i=0; i<numParms; i++){
        Parameter *p = parms[activeParm][i];
        FossilRangeGraphSkyline *derivedPtr = dynamic_cast<FossilRangeGraphSkyline *>(p);
        if ( derivedPtr != 0 )
            return derivedPtr;
    }
    return NULL;
}

SpeciationSkyline* Model::getActiveSpeciationSkyline(void) {
    
    for (int i=0; i<numParms; i++){
        Parameter *p = parms[activeParm][i];
        SpeciationSkyline *derivedPtr = dynamic_cast<SpeciationSkyline *>(p);
        if ( derivedPtr != 0 )
            return derivedPtr;
    }
    return NULL;
}

void Model::initializeConditionalLikelihoods(void) {

	// allocate conditional likelihoods
	int nNodes = 2*alignmentPtr->getNumTaxa()-1;
	int nChar  = alignmentPtr->getNumChar();
	int sizeOneNode = nChar * numGammaCats * 4;
	int sizeOneSpace = nNodes * sizeOneNode;
	cls = new double[2 * sizeOneSpace];
	for (int i=0; i<2*sizeOneSpace; i++)
		cls[i] = 0.0;
	for (int i=0; i<2; i++)
		{
		clPtr[i] = new double*[nNodes];
		for (int j=0; j<nNodes; j++)
			clPtr[i][j] = &cls[ i * sizeOneSpace + j * sizeOneNode ];
		}
		
	// initialize the tip conditional likelihoods
	for (int i=0; i<alignmentPtr->getNumTaxa(); i++)
		{
		double *cl0 = clPtr[0][i];
		double *cl1 = clPtr[1][i];
		for (int j=0; j<alignmentPtr->getNumChar(); j++)
			{
			int nucCode = alignmentPtr->getNucleotide(i, j);
			int possibleNucs[4];
			alignmentPtr->getPossibleNucs(nucCode, possibleNucs);
			for (int k=0; k<numGammaCats; k++)
				{
				for (int s=0; s<4; s++)
					{
					cl0[s] = (double)possibleNucs[s];
					cl1[s] = (double)possibleNucs[s];
					}
				cl0 += 4;
				cl1 += 4;
				}
			}
		}
}

void Model::initializeTransitionProbabilityMatrices(void) {

	int nNodes = 2*alignmentPtr->getNumTaxa()-1;
	for (int i=0; i<2; i++)
		{
		tis[i] = new MbMatrix<double>*[nNodes];
		tis[i][0] = new MbMatrix<double>[numGammaCats*nNodes];
		for (int j=1; j<nNodes; j++)
			tis[i][j] = tis[i][j-1] + numGammaCats;
		}
	for (int i=0; i<2; i++)
		for (int j=0; j<nNodes; j++)
			for (int k=0; k<numGammaCats; k++)
				tis[i][j][k] = MbMatrix<double>(4,4);
}


Parameter* Model::pickParmToUpdate(void) {

	double u = ranPtr->uniformRv();
	double sum = 0.0;
	Parameter *parm = NULL;
	for (unsigned i=0; i<updateProb.size(); i++)
		{
		sum += updateProb[i];
		if ( u < sum )
			{
			parm = parms[activeParm][i];
			break;
			}
		}
	return parm;
}

void Model::printTis(std::ostream & o) const {

	int nNodes = 2*alignmentPtr->getNumTaxa()-1;
	for (int j=0; j<nNodes; j++)
		{
		o << "Node " << j << endl;
		for (int a=0; a<4; a++)
			{
			for (int i=0; i<2; i++)
				{
				for (int k=0; k<numGammaCats; k++)
					{
					for (int b=0; b<4; b++)
						{
						o << fixed << setprecision(10) << tis[i][j][k][a][b] << " ";
						}
					}
				}
			o << '\n';
			}
		}
	o.flush();	
}

double Model::safeExponentiation(double lnX) {

	if (lnX < -300.0)
		return 0.0;
	else if (lnX > 0.0)
		return 1.0;
	else
		return exp(lnX);
}

void Model::setTiProb(void) {

	Tree *t     = getActiveTree();
	Shape *s    = getActiveShape();
	NodeRate *r = getActiveNodeRate();
	for (int n=0; n<t->getNumNodes(); n++)
		{
		Node *p = t->getDownPassNode(n);
		if (p->getAnc() != NULL && p->getIsTiDirty() == true) 
			{
			setTiProb(p, s, r);
			p->setIsTiDirty(false);
			}
		}
	// TAH root rate debug. This stuff below is stupid anyway
#	if ASSIGN_ROOT
	Node *roo = t->getRoot();
	t->setRootRateValue(r->getRateForNodeIndexed(roo->getIdx()));
#	endif
}

void Model::setTiProb(Node *p, Shape *s, NodeRate *r) {

	int activeTi = p->getActiveTi();
	int idx      = p->getIdx();
	Treescale *ts     = getActiveTreeScale();
	double sv = 1.0; //ts->getScaleValue();  // FIXPARM
	if(estAbsRts) sv = ts->getScaleValue();
	double branchProportion =  (p->getAnc()->getNodeDepth() - p->getNodeDepth()) * sv;
	double rP = r->getRateForNodeIndexed(p->getIdx());
#	if 0
	double rA = r->getRateForNodeIndexed(p->getAnc()->getIdx());
	double v = branchProportion * (rP + rA) * 0.5;
#	else
	double v = branchProportion * rP;
#	endif
	
	if(rP == 0.0){
		// then this means that this is calculating the lnl for nodes that are not assigned to 
		// RateGroups in currentRateGroups
		cerr << "ERROR: Problem rP = 0" << endl;
		exit(1);
	}
	
#	if 0
	for (int k=0; k<4; k++){
		double rt = s->getRate(k);
		tis[activeTi][idx][k] = tiCalculator->tiProbs( v*rt, tis[activeTi][idx][k] );
	}
#	else
	double rt = s->getRate(0);
	tis[activeTi][idx][0] = tiCalculator->tiProbs( v*rt, tis[activeTi][idx][0] );
	
	rt = s->getRate(1);
	tis[activeTi][idx][1] = tiCalculator->tiProbs( v*rt, tis[activeTi][idx][1] );
	
	rt = s->getRate(2);
	tis[activeTi][idx][2] = tiCalculator->tiProbs( v*rt, tis[activeTi][idx][2] );
	
	rt = s->getRate(3);
	tis[activeTi][idx][3] = tiCalculator->tiProbs( v*rt, tis[activeTi][idx][3] );
#	endif
	// set node info for printing
	//p->setBranchTime(branchProportion);
	//p->setRtGrpVal(rP);
}

void Model::setNodeRateGrpIndxs(void) {
	
	Tree *t     = getActiveTree();
	NodeRate *r = getActiveNodeRate();
#	if ASSIGN_ROOT
	int rtID = t->getNumNodes() + 2;
#	else
	int rtID = t->getRoot()->getIdx();  // TAH root rate debug
#	endif
	//int numRCats = r->getNumRateGroups();
	for(int n=0; n<t->getNumNodes(); n++){
		if(n != rtID){
			Node *p = t->getNodeByIndex(n);
			int tn = r->getTableNumForNodeIndexed(n);
			double srt = r->getRateForNodeIndexed(n);
			p->setRtGrpIdx(tn);
			p->setRtGrpVal(srt);
		}
	}
}

void Model::updateAccepted(void) {

	int from = activeParm, to;
	if (from == 0)
		to = 1;
	else
		to = 0;
	for (int i=0; i<numParms; i++)
		*parms[to][i] = *parms[from][i];
}

void Model::updateRejected(void) {

	int to = activeParm, from;
	if (to == 0)
		from = 1;
	else
		from = 0;
	for (int i=0; i<numParms; i++)
		*parms[to][i] = *parms[from][i];

}


void Model::upDateRateMatrix(void) {

	tiCalculator->updateQ( getActiveExchangeability()->getRate(), getActiveBasefreq()->getFreq() );
}

void Model::writeUnifTreetoFile(void) {
	
	Tree *t = getActiveTree();
	t->checkNodeCalibrationCompatibility();
	t->setAllNodeBranchTimes();
	string ts = t->getTreeDescription();
	ofstream out;
	out.open("uniformized_t.phy");
	out << ts << "\n";
	out.close();
	if(t->getIsCalibratedTree()){
		ofstream xmlo;
		xmlo.open("calib_beast.xml");
		t->writeCalNodeBEASTInfoXML(xmlo);
		xmlo.close();
	}
	else{
		ofstream xmlo;
		xmlo.open("calib_beast.xml");
		t->writeRRTNodeBEASTInfoXML(xmlo);
		xmlo.close();
	}
}

double Model::getMyCurrLnL(void) {
	
	if(lnLGood){
		lnLGood = false;
		return myCurLnL;
	}
	else
		return lnLikelihood();
}


double Model::readCalibFile(void) {
	
	/*
	3
	T1	T3	0.4	0.8
	T6	T7	0.1	0.2
	T10 T9  0.9	0.9
	*/
	
	/*
	TAH
	This also initializes the root height parameter
	I'm not sure how best to do this
	
	*/
	
	cout << "\nCalibrations:" << endl;
	bool rootIs = false;
	Calibration *rooCal;
	string ln = getLineFromFile(calibfilen, 1);
    string tg = "-s"; // total group fossil indicator
    string fx = "-x"; // fix node age
    int nlins = atoi(ln.c_str());
	int nnodes = alignmentPtr->getNumTaxa() - 1;
	string *calList = new string[nlins];
	for(int i=0; i<nlins; i++){
		calList[i] = getLineFromFile(calibfilen, i+2);
        if (calList[i].find(tg) != string::npos) {
            if (treeTimePrior != 8) {
                cerr << "ERROR: Total group fossils (-s) cannot be included with the -tga option!\nTry using -fbds." << endl;
                exit(1);
            }
        }
        
        if (calList[i].find(fx) != string::npos) {
            Calibration *cal = new Calibration(calList[i], 3);
            fixedNodes.push_back(cal);
            //rooCal = cal;
            //rootIs = true; // I dont know if we want to do this here
        }
        else if (fbdsExperimentalMode != 2 or ignoreFossils)  {
            Calibration *cal = new Calibration(calList[i], 0);
            calibrs.push_back(cal);
            if(cal->getIsRootCalib()){
                rooCal = cal;
                rootIs = true;
            }
        }
	}
	delete [] calList;

	double initTScale = 1.0;
	double yb = 0.0;
	double ob = 0.0;
    
    bool rootDone = false;
    
    // initialization of the origin time here
	
	double oldest = 0.0;
	for(vector<Calibration *>::iterator v = calibrs.begin(); v != calibrs.end(); v++){
		double fa = (*v)->getYngTime();
		if(fa > oldest)
			oldest = fa;
	}
	
    if(!fixedNodes.empty()){
        for(vector<Calibration *>::iterator v = fixedNodes.begin(); v != fixedNodes.end(); v++){
            if((*v)->getIsRootCalib()){
                fixRootHeight = true;
                yb = (*v)->getYngTime();
                initTScale = yb;
                ob = yb;
                rootDone = true;
                rHtY = yb;
                rHtO = ob;
                initOT = initTScale * 2.0;
            }
        }
    }
    
    if(!rootDone && fbdsExperimentalMode > 0) {
        cerr << "ERROR: You need to specify the age of the root in experiemental mode 1 or 2. Use -x root <age> in the calibration file." << endl;
        exit(1);
    }
    
    if(!rootDone){
        if(rootIs){
            if(rooCal->getPriorDistributionType() == 1){
                yb = rooCal->getYngTime();
                ob = rooCal->getOldTime();
                if(yb == ob){
                    initTScale = yb;
                    fixRootHeight = true;
                }
                else{
                    initTScale = yb + (ranPtr->uniformRv() * (ob - yb));
                    fixRootHeight = false;
                }
            }
            //else if(rooCal->getPriorDistributionType() == 2){
            else if(rooCal->getPriorDistributionType() > 1){
                fixRootHeight = false;
                yb = rooCal->getYngTime();
                double expMean = yb * 0.2;
                initTScale = yb + ranPtr->exponentialRv(1 / expMean);
            }
            rHtY = yb;
            rHtO = ob;
            initOT = initTScale * 2.0;
        }
        else{
            yb = 0.0;
            fixRootHeight = false;
            for(vector<Calibration *>::iterator v = calibrs.begin(); v != calibrs.end(); v++){
                double tmpv;
                if((*v)->getPriorDistributionType() == 1)
                    tmpv = (*v)->getOldTime();
                //else if((*v)->getPriorDistributionType() == 2)
                else if((*v)->getPriorDistributionType() > 1)
                    tmpv = (*v)->getYngTime() * 1.1;
                if(tmpv > yb)
                    yb = tmpv;
            }
            ob = yb + (yb * 2);
            double tsc = yb + (ranPtr->uniformRv() * (ob - yb));
            initTScale = tsc;
            rHtY = yb;
            rHtO = ob;
            initOT = tsc * 2.0;
        }
    }
	
	if(nlins == nnodes){
		bool fixall = true;
		for(vector<Calibration *>::iterator v = calibrs.begin(); v != calibrs.end(); v++){
			double tmpo = (*v)->getOldTime();
			double tmpy = (*v)->getYngTime();
			if(tmpo != tmpy){
				fixall = false;
				break;
			}
		}
		zeroNodeTimeMove = fixall;
	}
	else if(nlins == nnodes-1 && rootIs == false){
		bool fixall = true;
		for(vector<Calibration *>::iterator v = calibrs.begin(); v != calibrs.end(); v++){
			double tmpo = (*v)->getOldTime();
			double tmpy = (*v)->getYngTime();
			if(tmpo != tmpy){
				fixall = false;
				break;
			}
		}
		zeroNodeTimeMove = fixall;
	}
	
	if(initOT < oldest){
		double oMxTmp = oldest * 4.0;
		initOT = oldest + (ranPtr->uniformRv() * (oMxTmp - oldest));
	}
	
	cout << "\nInitial root height : " << initTScale <<  " [" << yb << ", " << ob << "]" << endl;
	return initTScale;
}

void Model::readTipDateFile(void){
	
	/*
	    3
		X50    20.4
		X59    40.1
		
	*/
	string ln = getLineFromFile(tipDateFileN, 1);
	int nlins = atoi(ln.c_str());
	string *calList = new string[nlins];
	for(int i=0; i<nlins; i++){
		calList[i] = getLineFromFile(tipDateFileN, i+2);
		Calibration *cal = new Calibration(calList[i], 1);
		tipDates.push_back(cal);
	}
    
	delete [] calList;

}

void Model::readOccurrenceFile(void){
    
    /*
	    3
     20.4
     40.1
     90.0
     */
    
    // initialize occurrences
    cout << "\nOccurrences:" << endl;
    string ln = getLineFromFile(calibfilen, 1);
    int nlins = atoi(ln.c_str());
    string *calList = new string[nlins];
    for(int i=0; i<nlins; i++){
        numFossils++;
        calList[i] = getLineFromFile(calibfilen, i+2);
        Calibration *cal = new Calibration(calList[i], 2);
        calibrs.push_back(cal);
    }
    delete [] calList;
    
    // initialize origin time
    double yb = 0.0;
    double ob = 0.0;
    
    for(vector<Calibration *>::iterator v = calibrs.begin(); v != calibrs.end(); v++){
        double tmpv;
        tmpv = (*v)->getYngTime() * 1.1;
        if(tmpv > yb)
            yb = tmpv;
    }
    ob = yb + (yb * 1.2);
    initOT = yb + (ranPtr->uniformRv() * (ob - yb));
    rHtY = yb;
    
    // initialize terminal time
    for(vector<Calibration *>::iterator v = calibrs.begin(); v != calibrs.end(); v++){
        double yt;
        yt = (*v)->getYngTime();
        if(yt < yb)
            yb = yt;
    }
    rHtY = yb;
    
    cout << "\nTotal number of occurrences: " << numFossils << endl;
    cout << "\nInitial origin time: " << initOT <<  " [" << yb << ", " << ob << "]" << endl;
    cout << "\nTerminal time: " << rHtY << endl;
    
}


void Model::readFossilRangeFile(void){
    /*
    n k
    yi oi
    yj oj etc.
    */
    
    /* OR to specify fixed attachment times
     n k
     yi oi bi
     yj oj bj etc.
     */
    
    cout << "\nFossil ranges:" << endl;
    string ln = getLineFromFile(calibfilen, 1);
    int nlins = atoi(ln.c_str());
    
    // fetch the number of lineages and fossils
    stringstream ss;
    string tmp = "";
    ss << ln;
    ss >> tmp;
    numLineages = atoi(tmp.c_str());
    ss >> tmp;
    numFossils = atoi(tmp.c_str());
    
    string *calList = new string[nlins];
    for(int i=0; i<nlins; i++){
        calList[i] = getLineFromFile(calibfilen, i+2);
        Calibration *cal = new Calibration(calList[i], 4);
        calibrs.push_back(cal);
    }
    
    delete [] calList;
    
    //initOT = ot * 2.0; // origin time initiated by class frg

    cout << "\nTotal number of lineages: " << numLineages << endl;
    cout << "\nTotal number of fossils: " << numFossils << endl;

}

void Model::readPresenceAbsenceFile(void){
    
    cout << "\nFossil ranges:" << endl;
    string lnCal = getLineFromFile(calibfilen, 1);
    string lnPa = getLineFromFile(pafilen, 1);
    
    int nlins = atoi(lnCal.c_str());
    int nlinsPa = atoi(lnPa.c_str());
    
    if(nlins != nlinsPa){
        cout << "ERROR: There's a problem with the calibration files \n";
        exit(1);
    }
    
    // fetch the number of lineages and fossils
    stringstream ss;
    string tmp = "";
    ss << lnCal;
    ss >> tmp;
    numLineages = atoi(tmp.c_str());
    ss >> tmp;
    numFossils = atoi(tmp.c_str());
    
    stringstream pp;
    tmp = "";
    pp << lnPa;
    pp >> tmp;
    int lnNum = atoi(tmp.c_str());
    pp >> tmp;
    int intNum = atoi(tmp.c_str());
    
    // do some cross checking here with the pa file
    //if(numLineages != lnNum){
    //    cout << "ERROR: There's a problem with the calibration files \n";
    //    exit(1);
    //}
    
    string *calList = new string[nlins];
    string *paList = new string[nlins];
    
    for(int i=0; i<nlins; i++){
        calList[i] = getLineFromFile(calibfilen, i+2);
        paList[i] = getLineFromFile(pafilen, i+2);
        Calibration *cal = new Calibration(calList[i], paList[i], intNum);
        calibrs.push_back(cal);
    }
    
    delete [] calList;
    
    cout << "\nTotal number of lineages: " << numLineages << endl;
    cout << "\nTotal number of fossils: " << numFossils << endl;
    
}

void Model::readIntervalsFile(void){
    /*
     i
     h1
     h2 etc.
     */
    cout << "\nIntervals:" << endl;
    string ln = getLineFromFile(intfilen, 1);
    int nlins = atoi(ln.c_str());
    
    // fetch the number of intervals
    stringstream ss;
    string tmp = "";
    ss << ln;
    ss >> tmp;
    userSpecifiedIntervals = atoi(tmp.c_str());
    
    string *intList = new string[nlins];
    for(int i=0; i<nlins; i++){
        intList[i] = getLineFromFile(intfilen, i+2);
        Calibration *interval = new Calibration(intList[i], 5);
        intervals.push_back(interval);
    }
    delete [] intList;
    
    cout << "\nTotal number of user defined intervals: " << userSpecifiedIntervals << endl;
    
}


Calibration* Model::getRootCalibration(void) {

	for(vector<Calibration *>::iterator v = calibrs.begin(); v != calibrs.end(); v++){
		if((*v)->getIsRootCalib())
			return (*v);
	}
	return NULL;
}


void Model::setUpdateProbabilities(bool initial) {

	double bfp, srp, shp, ntp, dpp, cpa, tsp, spp, ehp, fcp, otp;
	if(initial){
		bfp = 0.3;
		srp = 0.3;
		shp = 0.3;
		ntp = 0.6;
		dpp = 0.5;
		cpa = 0.3;
		tsp = 0.5;
		spp = 0.4;
		ehp = 0.0;
		fcp = 0.0;
        otp = 0.0;
	}
	else{
		bfp = 0.2;
		srp = 0.2;
		shp = 0.2;
		ntp = 0.4;
		dpp = 0.5;
		cpa = 0.3;
		tsp = 0.4;
		spp = 0.4;
		ehp = 0.0;
        otp = 0.0;
	}
	if(turnedOffMove == 1)
		bfp = 0.0;
	else if(turnedOffMove == 2)
		srp = 0.0;
	else if(turnedOffMove == 3)
		shp = 0.0;
	else if(turnedOffMove == 4){
		ntp = 0.0;
		tsp = 0.0;
		spp = 0.0;
		ehp = 0.0;
	}
	else if(turnedOffMove == 5){
		dpp = 0.0;
		cpa = 0.0;
		setNodeRateGrpIndxs(); // set the rates on the tree
		*parms[0][3] = *parms[1][3]; //make sure both trees are the same
	}
	else if(turnedOffMove == 6 || cpfix == true) // if this move is turned off then the cp is set using the prior mean number of groups
		cpa = 0.0;
	else if(turnedOffMove == 8)
		spp = 0.0;
	if(fixRootHeight)
		tsp = 0.0;
	if(treeTimePrior == 1) // treeTimePrior = 4 is now BDSS || treeTimePrior == 4)
		spp = 0.0;
	if(zeroNodeTimeMove == 1){
		ntp = 0.0;
		cout << "All internal node times are fixed" << endl;
	}
	if(exponCalibHyperParm){
		if(initial)
			ehp = 0.3;
		else
			ehp = 0.4;
	}
	if(fixSomeModParams){
		bfp = 0.0;
		shp = 0.0;
	}
	
	if(treeTimePrior > 3) // might want to change this
		spp = 0.5;
	
	if(treeTimePrior > 5){
		ntp = 0.8;
		fcp = 0.4;
	}
	if(fixedClockRate > 0.0){ // FIXPARM
		cpa = 0.0;
		dpp = 0.0;
		//spp = 0.0;
	}
	
	// TAH: FIXING 6 Feb
	if(fixTestRun){
		spp = 0.0;
		//tsp = 0.0;
		//ntp = 0.0;
	}
    
    if(true){ // TAH: fixing for test
        ntp = 0.5;
        tsp = 0.5;
        spp = 0.5;
        otp = 0.5; // set to 0.0 to fix origin to initOT
    }
    
    if(fbdsExperimentalMode > 0){
        bfp = 0.0;
        srp = 0.0;
        shp = 0.0;
        ntp = 0.4;
        dpp = 0.0;
        cpa = 0.0;
        tsp = 0.0;
        spp = 0.4;
        ehp = 0.0;
        fcp = 0.0;
        otp = 0.0;
    }
	
	updateProb.clear();
	updateProb.push_back(bfp); // 1 basefreq
	updateProb.push_back(srp); // 2 sub rates
	updateProb.push_back(shp); // 3 gamma shape
	updateProb.push_back(ntp); // 4 node times
	updateProb.push_back(dpp); // 5 dpp rates
	updateProb.push_back(cpa); // 6 concentration parameter
	updateProb.push_back(tsp); // 7 tree scale parameter
	updateProb.push_back(spp); // 8 speciation parameters
	updateProb.push_back(ehp); // 9 exponential calibration hyper priors
    updateProb.push_back(otp); // 10 origin time parameter
	double sum = 0.0;
	for (unsigned i=0; i<updateProb.size(); i++)
		sum += updateProb[i];
	for (unsigned i=0; i<updateProb.size(); i++)
		updateProb[i] /= sum;
}

void Model::updateTreeRates(void){
	
	for (int i=0; i<numParms; i++){
		Parameter *p1 = parms[0][i];
		Tree *dp1 = dynamic_cast<Tree *>(p1);
		if ( dp1 != 0 ){
			dp1->setNodeRateValues();
		}
		Parameter *p2 = parms[1][i];
		Tree *dp2 = dynamic_cast<Tree *>(p2);
		if ( dp2 != 0 ){
			dp2->setNodeRateValues();
		}
	}
}




// end model



