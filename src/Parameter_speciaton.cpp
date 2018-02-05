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


#include "Parameter.h"
#include "Parameter_speciaton.h"
#include "Parameter_fossilgraph.h"
#include "Parameter_fossilrangegraph.h"
#include "Parameter_treescale.h"
#include "Parameter_tree.h"
#include "MbRandom.h"
#include "Model.h"
#include <iostream>
#include <iomanip>
#include <sstream>

//#include <unistd.h>

using namespace std;

Speciation::Speciation(MbRandom *rp, Model *mp, double bdr, double bda, double bds, double initRH, double rh, int bdp, bool fxPsi, double psi, int specPr,
                       int psiPr, double bPrRate, double dPrRate, double pPrRate) : Parameter(rp, mp) {
	
	
    maxdivV = 30000.0;
	name = "SP";
    relativeDeath = ranPtr->uniformRv();
    netDiversificaton = ranPtr->uniformRv();
    probSpeciationS = ranPtr->uniformRv();
    fossilRate = psi;
    birthRate = 0.02;
    deathRate = 0.01;
    extantSampleRate = rh;
	treeTimePrior = modelPtr->getTreeTimePriorNum();
	currentFossilGraphLnL = 0.0;
    currentFossilRangeGraphLnL = 0.0;
	parameterization = bdp;
	if(treeTimePrior == 9){
		parameterization = 3;
		extantSampleRate = 0.0;
	}
    fixPsi = fxPsi;
	setAllBDFossParams();
    
    // priors on birth death paras
    deathRatePrior = specPr; // 1 = unifrom prior, 2 = exponential prior
    birthRatePrior = specPr;
    fossilSamplingRatePrior = psiPr;
    netDivRatePrior = 1; // only used by FG and FRG models
    birthRateExpRate = bPrRate;
    deathRateExpRate = dPrRate;
    fossilSamplingRateExpRate = pPrRate;
    netDivRateExpRate = 10.0;
	
	if(mp->getFixTestRun()){
		relativeDeath = bda;
		netDiversificaton = bdr;
		probSpeciationS = bds;
	}
	setAllBDFossParams();	
	
	if(relativeDeath < 0.0 && netDiversificaton < 0.0){
		relativeDeath = ranPtr->uniformRv();
		netDiversificaton = ranPtr->uniformRv() * maxdivV;
	}
	if(relativeDeath >= 1.0){
		cerr << "ERROR: The relative death rate (-bda) is >= 1 " << endl;
		exit(1);
	}
	if(netDiversificaton <= 0.0){
		cerr << "ERROR: The net diversification rate (-bdr) is <= 0 " << endl;
		exit(1);
	}
	if(treeTimePrior == 1){
		relativeDeath = 0.0;
		netDiversificaton = 0.0;
	}	
	else if(treeTimePrior == 2)
		relativeDeath = 0.0;
	if(netDiversificaton >= maxdivV)
		netDiversificaton = maxdivV * maxdivV / netDiversificaton;
		
	if(treeTimePrior > 5){
        cout << "Paramterization: " << parameterization << endl;
		cout << "Speciaton parameters are initialized with: d = " << netDiversificaton << " , r = " << relativeDeath<< " , s = " << probSpeciationS << endl;
		cout << "                                   l = " << birthRate << " , m = " << deathRate << " , psi = " << fossilRate << " , rho = " << extantSampleRate << endl;
	}
    cout << "BD initialized\n";
	
}

Speciation::~Speciation(void) {
	
}

Speciation& Speciation::operator=(const Speciation &c) {
	
	if (this != &c)
		clone(c);
	return *this;
}

void Speciation::clone(const Speciation &c) {
	
	setAllBDFossParams();
	relativeDeath = c.relativeDeath;
	netDiversificaton = c.netDiversificaton;
	probSpeciationS = c.probSpeciationS;
	birthRate = c.birthRate;
	deathRate = c.deathRate;
	fossilRate = c.fossilRate;
	parameterization = c.parameterization;
	extantSampleRate = c.extantSampleRate;
	treeTimePrior = c.treeTimePrior;
	name = "SP";
}

void Speciation::print(std::ostream & o) const {
	if(treeTimePrior > 1){
		o << "Speciaton parameters: d/b = ";
		o << fixed << setprecision(4) << relativeDeath << " , b-d = ";
		o << fixed << setprecision(4) << netDiversificaton << " ";
		o << fixed << setprecision(4) << probSpeciationS << " ";
		o << fixed << setprecision(4) << fossilRate << " ";
		o << fixed << setprecision(4) << birthRate << " ";
		o << fixed << setprecision(4) << deathRate << " ";
		o << endl;
	}
}

double Speciation::update(double &oldLnL) {
	
	double lnR = 0.0;
	if(treeTimePrior == 9){
		return updateFossilGraphBDParams(oldLnL);
	}
    else if(treeTimePrior == 10){
        return updateFossilRangeGraphBDParams(oldLnL);
    }
	else{
		return updateTreeBDParams(oldLnL);
	}
	return lnR;
}

double Speciation::updateTreeBDParams(double &oldLnL) {
	
	double lnR = 0.0;
	//int v;

	Tree *t = modelPtr->getActiveTree();
	double oldtreeprob = getLnTreeProb(t); 
	double lnProposalRatio = 0.0;
	
	if(treeTimePrior == 2){
		lnProposalRatio += updateNetDivRate();
		relativeDeath = 0.0;
	}
	else if(treeTimePrior == 3){
		updateRelDeathRt(t);
		updateNetDivRate(t);
		modelPtr->setLnLGood(true);
		modelPtr->setMyCurrLnl(oldLnL);
		
		return 0.0;
	}
	else if(treeTimePrior > 3){
		
		//v = (int)(ranPtr->uniformRv() * 3);
		
		if(parameterization == 1){
            updateRelDeathRt(t); // r
            updateNetDivRate(t); // d
            updateBDSSFossilProbS(t); // s
		}
		else if(parameterization == 2){
			updateRelDeathRt(t); // r
			updateNetDivRate(t); // d
			updatePsiRate(t); // psi
		}
		// this parameterization is currently problematic
		// because it allows mu >> lambda
		else if(parameterization == 3){
            updateBirthRate(t); // lambda
			updateDeathRate(t); // mu
			updatePsiRate(t); // psi
		}
		else if(parameterization == 4){
			updateNetDivRate(t); // d
            updateDeathRate(t); // mu
			updatePsiRate(t); // psi
		}
		modelPtr->setLnLGood(true);
		modelPtr->setMyCurrLnl(oldLnL);
		return 0.0;
		
	}
	double newtreeprob = getLnTreeProb(t); 
	double lnPriorRatio = (newtreeprob - oldtreeprob);
	lnR = lnPriorRatio + lnProposalRatio;
	
	modelPtr->setLnLGood(true);
	modelPtr->setMyCurrLnl(oldLnL);
	return lnR;
}

double Speciation::updateFossilGraphBDParams(double &oldLnL){

	int v;
	currentFossilGraphLnL = oldLnL;
	FossilGraph *fg = modelPtr->getActiveFossilGraph();
	int numMoves = 3;
	if(parameterization == 3) numMoves = 4;
	v = (int)(ranPtr->uniformRv() * numMoves);
	
	if(parameterization == 1){
		if(v == 0)
			updateRelDeathRt(fg); // r
		else if(v == 1)
			updateNetDivRate(fg); // d
		else
			updateBDSSFossilProbS(fg); // s
	}
	else if(parameterization == 2){
		if(v == 0)
			updateRelDeathRt(fg); // r
		else if(v == 1)
			updateNetDivRate(fg); // d
		else
			updatePsiRate(fg); // psi
	}
	else if(parameterization == 3){  //****//
		if(v == 0)
			updateDeathRate(fg); // mu
		else if(v == 1)
			updateBirthRate(fg); // lambda
		else if(v == 2)
			updatePsiRate(fg); // psi
		else
			updateBirthAndDeath(fg); // mu and lambda
	}
	else if(parameterization == 4){
		if(v == 0)
			updateNetDivRate(fg); // d
		else if(v == 1)
			updateDeathRate(fg); // mu
		else
			updatePsiRate(fg); // psi
	}
	return currentFossilGraphLnL;
}


double Speciation::updateFossilRangeGraphBDParams(double &oldLnL){
    
    int v, numMoves;
    currentFossilRangeGraphLnL = oldLnL;
    FossilRangeGraph *frg = modelPtr->getActiveFossilRangeGraph();
    
    if(fixPsi)
        numMoves = 2;
    else
        numMoves = 3;
    
    //if(parameterization == 3) numMoves = 4;
    v = (int)(ranPtr->uniformRv() * numMoves);
    
    if(parameterization == 1){
        if(v == 0)
            updateRelDeathRt(frg); // r
        else if(v == 1)
            updateNetDivRate(frg); // d
        else
            updateBDSSFossilProbS(frg); // s
    }
    else if(parameterization == 2){
        if(v == 0)
            updateRelDeathRt(frg); // r
        else if(v == 1)
            updateNetDivRate(frg); // d
        else
            updatePsiRate(frg); // psi
    }
    else if(parameterization == 3){  //
        if(v == 0)
            updateDeathRate(frg); // mu
        else if(v == 1)
            updateBirthRate(frg); // lambda
        else if(v == 2)
            updatePsiRate(frg); // psi
        //else
          //  updateBirthAndDeath(fg); // mu and lambda
    }
    
    //cout << "Now = " << netDiversificaton << " , r = " << relativeDeath<< " , s = " << probSpeciationS << endl; // debugging
    
    return currentFossilRangeGraphLnL;
}

double Speciation::updateRelDeathRt(void) {
	
	double rdwindow = 0.2;
	double oldRD = relativeDeath;
	double newRD = getNewValSWindoMv(oldRD, 0.0, 0.99999, rdwindow);
	relativeDeath = newRD;
	return 0.0;
	
}

double Speciation::updateRelDeathRt(Tree *t) {
	
	double oldtreeprob = getLnTreeProb(t);
	double lnPropR = 0.0;
	double rdwindow = 0.2;
	double oldRD = relativeDeath;
	double newRD = getNewValSWindoMv(oldRD, 0.0, 0.99999, rdwindow);
	relativeDeath = newRD;
	double newtreeprob = getLnTreeProb(t); 
	double lnPriorRatio = (newtreeprob - oldtreeprob);
	double lnR = lnPriorRatio + lnPropR;
	double r = modelPtr->safeExponentiation(lnR);
	if(ranPtr->uniformRv() < r){
		setAllBDFossParams();
	}
	else{
		relativeDeath = oldRD;
		setAllBDFossParams();
	}
	return 0.0;
	
}

double Speciation::updateRelDeathRt(FossilGraph *fg) {
    
	double oldfgprob = currentFossilGraphLnL;
	double lnPropR = 0.0;
	double rdwindow = 0.2;
	double oldRD = relativeDeath;
	double newRD = getNewValSWindoMv(oldRD, 0.0, 0.99999, rdwindow);
	relativeDeath = newRD;
	double newfgprob = getLnFossilGraphProb(fg);
	double lnLikeRatio = (newfgprob - oldfgprob);
	double lnR = lnLikeRatio + lnPropR;
	double r = modelPtr->safeExponentiation(lnR);
	if(ranPtr->uniformRv() < r){
		setAllBDFossParams();
		currentFossilGraphLnL = newfgprob;
	}
	else{
		relativeDeath = oldRD;
		setAllBDFossParams();
		currentFossilGraphLnL = oldfgprob;
	}
	return 0.0;
	
}

double Speciation::updateRelDeathRt(FossilRangeGraph *frg) {
    
    double oldfgprob = currentFossilRangeGraphLnL;
    double lnPropR = 0.0;
    double rdwindow = 0.2;
    double oldRD = relativeDeath;
    double newRD = getNewValSWindoMv(oldRD, 0.0, 0.99999, rdwindow);
    relativeDeath = newRD;
    double newfgprob = getLnFossilRangeGraphProb(frg);
    double lnLikeRatio = (newfgprob - oldfgprob);
    double lnR = lnLikeRatio + lnPropR;
    double r = modelPtr->safeExponentiation(lnR);
    if(ranPtr->uniformRv() < r){
        setAllBDFossParams();
        currentFossilRangeGraphLnL = newfgprob;
    }
    else{
        relativeDeath = oldRD;
        setAllBDFossParams();
        currentFossilRangeGraphLnL = oldfgprob;
    }
    return 0.0;
    
}

double Speciation::getNewValScaleMv(double &nv, double ov, double vmin, double vmax, double tv){
	
	double rv = ranPtr->uniformRv();
	double c = tv * (rv - 0.5);
	double newcv = ov * exp(c);
	bool validV = false;
	do{
		if(newcv < vmin)
			newcv = vmin * vmin / newcv;
		else if(newcv > vmax)
			newcv = vmax * vmax / newcv;
		else
			validV = true;
	} while(!validV);
	nv = newcv;
	return c;
}

double Speciation::getNewValSWindoMv(double ov, double vmin, double vmax, double tv){
	double nv = 0.0;
	double u = ranPtr->uniformRv(-0.5,0.5) * (tv);
	nv = ov + u;
	bool validV = false;
	do{
		if(nv < vmin)
			nv = 2.0 * vmin - nv;
		else if(nv > vmax)
			nv = (2.0 * vmax) - nv;
		else
			validV = true;
	}while(!validV);
	return nv;
}

double Speciation::getNewValUpDownScaleMv(double &nv1, double ov1, double &nv2, double ov2, double sf){
	
	double u = ranPtr->uniformRv();
	double c = exp(sf * (u - 0.5));
	nv1 = ov1 * c;
	nv2 = ov2 / c;
	
	return 0.0;
}

double Speciation::updateNetDivRate(void) {
	
	double lpr = 0.0;
	double oldND = netDiversificaton;
	double newND;	
	double tuning = log(2.0);
	double minV = 0.0001;
	double c = getNewValScaleMv(newND, oldND, minV, maxdivV, tuning);
	netDiversificaton = newND;
	lpr = c; 
	return lpr;
}

double Speciation::updateNetDivRate(Tree *t) {
	
	double oldtreeprob = getLnTreeProb(t);
	double lpr = 0.0;
	double oldND = netDiversificaton;
	double newND;	
	double tuning = log(2.0);
	double minV = 0.0001;
	double c = getNewValScaleMv(newND, oldND, minV, maxdivV, tuning);
	netDiversificaton = newND;
	lpr = c;
	double newtreeprob = getLnTreeProb(t); 
	double lnPriorRatio = (newtreeprob - oldtreeprob);
	double lnR = lnPriorRatio + lpr;
	double r = modelPtr->safeExponentiation(lnR);
	if(ranPtr->uniformRv() < r){
		setAllBDFossParams();
	}
	else{
		netDiversificaton = oldND;
		setAllBDFossParams();
	}
	return 0.0;
}

double Speciation::updateNetDivRate(FossilGraph *fg) {
    
    double oldfgprob = currentFossilGraphLnL;
    double lpr = 0.0;
    double oldND = netDiversificaton;
    double newND;
    double tuning = log(2.0);
    double minV = 0.000001;
    double c = getNewValScaleMv(newND, oldND, minV, maxdivV, tuning);
    netDiversificaton = newND;
    lpr = c;
    double newfgprob = getLnFossilGraphProb(fg);
    double lnPriorRat = getExpPriorRatio(oldND, newND, netDivRateExpRate, netDivRatePrior);
    double lnLikeRat = (newfgprob - oldfgprob);
    double lnR = lnLikeRat + lpr + lnPriorRat;
    double r = modelPtr->safeExponentiation(lnR);
    if(ranPtr->uniformRv() < r){
        setAllBDFossParams();
        currentFossilGraphLnL = newfgprob;
    }
    else{
        netDiversificaton = oldND;
        setAllBDFossParams();
        currentFossilGraphLnL = oldfgprob;
    }
    return 0.0;
}

double Speciation::updateRealValNetDiv(FossilGraph *fg) {
    
	// a move that allows d to be a negative value
    double oldfgprob = currentFossilGraphLnL;
    double oldND = netDiversificaton;
    double window = 0.2;
    double minV = -1000.0;
    double newND = getNewValSWindoMv(oldND, minV, maxdivV, window);
    netDiversificaton = newND;
    double newfgprob = getLnFossilGraphProb(fg);
    double lnPriorRat = getExpPriorRatio(oldND, newND, netDivRateExpRate, netDivRatePrior);
    double lnLikeRat = (newfgprob - oldfgprob);
    double lnR = lnLikeRat + lnPriorRat;
    double r = modelPtr->safeExponentiation(lnR);
    if(ranPtr->uniformRv() < r){
        setAllBDFossParams();
        currentFossilGraphLnL = newfgprob;
    }
    else{
        netDiversificaton = oldND;
        setAllBDFossParams();
        currentFossilGraphLnL = oldfgprob;
    }
	
	return 0.0;
}

double Speciation::updateNetDivRate(FossilRangeGraph *frg) {
    
    double oldfgprob = currentFossilRangeGraphLnL;
    double lpr = 0.0;
    double oldND = netDiversificaton;
    double newND;
    double tuning = log(2.0);
    double minV = 0.000001;
    double c = getNewValScaleMv(newND, oldND, minV, maxdivV, tuning);
    netDiversificaton = newND;
    lpr = c;
    double newfgprob = getLnFossilRangeGraphProb(frg);
    double lnPriorRat = getExpPriorRatio(oldND, newND, netDivRateExpRate, netDivRatePrior);
    double lnLikeRat = (newfgprob - oldfgprob);
    double lnR = lnLikeRat + lpr + lnPriorRat;
    double r = modelPtr->safeExponentiation(lnR);
    if(ranPtr->uniformRv() < r){
        setAllBDFossParams();
        currentFossilRangeGraphLnL = newfgprob;
    }
    else{
        netDiversificaton = oldND;
        setAllBDFossParams();
        currentFossilRangeGraphLnL = oldfgprob;
    }
    return 0.0;
}


double Speciation::updateBDSSSampleProbRho(void) {
	
	double rdwindow = 0.2;
	double oldSP = extantSampleRate;
	double u;
	double newSP;
	u = ranPtr->uniformRv(-0.5,0.5) * (rdwindow);
	newSP = oldSP + u;
	bool validV = false;
	do{
		if(newSP < 0.0)
			newSP = 0.0 - newSP;
		else if(newSP > 0.9999)
			newSP = (2 * 0.9999) - newSP;
		else
			validV = true;
	}while(!validV);
	extantSampleRate = newSP;
	return 0.0;
}

double Speciation::updateBDSSSampleProbRho(FossilGraph *fg) {
    
    double oldfgprob = getLnFossilGraphProb(fg);
    double lnPropR = 0.0;
    double rdwindow = 0.2;
    double oldSP = extantSampleRate;
    double newSP = getNewValSWindoMv(oldSP, 0.0, 0.99999, rdwindow);
    extantSampleRate = newSP;
    double newfgprob = getLnFossilGraphProb(fg);
    double lnLikeRatio = (newfgprob - oldfgprob);
    double lnR = lnLikeRatio + lnPropR;
    double r = modelPtr->safeExponentiation(lnR);
    if(ranPtr->uniformRv() < r){
        //setAllBDFossParams();
        // move not yet implemented
        currentFossilGraphLnL = newfgprob;
    }
    else{
        extantSampleRate = oldSP;
        //setAllBDFossParams();
        // move not yet implemented
        currentFossilGraphLnL = oldfgprob;
    }
    return 0.0;
}

double Speciation::updateBDSSFossilProbS(void) {
	
	double swindow = 0.2;
	double oldS = probSpeciationS;
	double u;
	double newS;
	u = ranPtr->uniformRv(-0.5,0.5) * (swindow);
	newS = oldS + u;
	bool validV = false;
	do{
		if(newS < 0.0)
			newS = 0.0 - newS;
		else if(newS > 0.9999)
			newS = (2 * 0.9999) - newS;
		else
			validV = true;
	}while(!validV);
	probSpeciationS = newS;
	return 0.0;
}

double Speciation::updateBDSSFossilProbS(Tree *t) {
	
    double oldtreeprob = getLnTreeProb(t);
	double lnPropR = 0.0;
	double swindow = 0.2;
	double oldS = probSpeciationS;
	double newS = getNewValSWindoMv(oldS, 0.0, 0.99999, swindow);
	probSpeciationS = newS;
	double newtreeprob = getLnTreeProb(t); 
	double lnPriorRatio = (newtreeprob - oldtreeprob);
	double lnR = lnPriorRatio + lnPropR;
	double r = modelPtr->safeExponentiation(lnR);
	if(ranPtr->uniformRv() < r){
		setAllBDFossParams();
	}
	else{
		probSpeciationS = oldS;
		setAllBDFossParams();
	}
	return 0.0;

}

double Speciation::updateBDSSFossilProbS(FossilGraph *fg) {
    
    double oldfgprob = currentFossilGraphLnL;
    double lnPropR = 0.0;
    double swindow = 0.2;
    double oldS = probSpeciationS;
    double newS = getNewValSWindoMv(oldS, 0.0, 0.99999, swindow);
    probSpeciationS = newS;
    double newfgprob = getLnFossilGraphProb(fg);
    double lnPriorRatio = (newfgprob - oldfgprob);
    double lnR = lnPriorRatio + lnPropR;
    double r = modelPtr->safeExponentiation(lnR);
    if(ranPtr->uniformRv() < r){
        setAllBDFossParams();
        currentFossilGraphLnL = newfgprob;
    }
    else{
        probSpeciationS = oldS;
        setAllBDFossParams();
        currentFossilGraphLnL = oldfgprob;
    }
    return 0.0;
    
}

double Speciation::updateBDSSFossilProbS(FossilRangeGraph *frg) {
    
    double oldfgprob = currentFossilRangeGraphLnL;
    double lnPropR = 0.0;
    double swindow = 0.2;
    double oldS = probSpeciationS;
    double newS = getNewValSWindoMv(oldS, 0.0, 0.99999, swindow);
    probSpeciationS = newS;
    double newfgprob = getLnFossilRangeGraphProb(frg);
    double lnPriorRatio = (newfgprob - oldfgprob);
    double lnR = lnPriorRatio + lnPropR;
    double r = modelPtr->safeExponentiation(lnR);
    if(ranPtr->uniformRv() < r){
        setAllBDFossParams();
        currentFossilRangeGraphLnL = newfgprob;
    }
    else{
        probSpeciationS = oldS;
        setAllBDFossParams();
        currentFossilRangeGraphLnL = oldfgprob;
    }
    return 0.0;
    
}

double Speciation::updatePsiRate(Tree *t) {
	
	double lpr = 0.0;
	double oldtreeprob = getLnTreeProb(t);
	double oldPsi = fossilRate;
	double newPsi;	
	double tuning = log(2.0);
	double rv = ranPtr->uniformRv();
	double c = tuning * (rv - 0.5);
	newPsi = oldPsi * exp(c);
	double minV = 0.0001;
	double maxV = 100.00;
	bool validV = false;
	do{
		if(newPsi < minV)
			newPsi = minV * minV / newPsi;
		else if(newPsi > maxV)
			newPsi = newPsi * maxV / newPsi;
		else
			validV = true;
	} while(!validV);
	fossilRate = newPsi;
	probSpeciationS = newPsi / (deathRate + newPsi) ;
	lpr = c;
	double newtreeprob = getLnTreeProb(t); 
	double lnPriorRatio = (newtreeprob - oldtreeprob);
	double lnR = lnPriorRatio + lpr;
	double r = modelPtr->safeExponentiation(lnR);
	if(ranPtr->uniformRv() < r){
		setAllBDFossParams();
	}
	else{
		fossilRate = oldPsi;
		setAllBDFossParams();
	}
	return 0.0;
}

double Speciation::updatePsiRate(FossilGraph *fg) {
    
    double oldfgprob = currentFossilGraphLnL;
    double lpr = 0.0;
    double oldPsi = fossilRate;
    double newPsi;
    double tuning = log(2.0);
    double minV = 0.0001;
    double c = getNewValScaleMv(newPsi, oldPsi, minV, maxdivV, tuning);
    fossilRate = newPsi;
    lpr = c;
    double newfgprob = getLnFossilGraphProb(fg);
    double lnPriorRat = getExpPriorRatio(oldPsi, newPsi, fossilSamplingRateExpRate, fossilSamplingRatePrior);
    double lnLikeRat = (newfgprob - oldfgprob);
    double lnR = lnLikeRat + lpr + lnPriorRat;
    double r = modelPtr->safeExponentiation(lnR);
    
    if(ranPtr->uniformRv() < r){
        currentFossilGraphLnL = newfgprob;
    }
    else{
        fossilRate = oldPsi;
        setAllBDFossParams();
        currentFossilGraphLnL = oldfgprob;
    }
    return 0.0;
}

double Speciation::updatePsiRate(FossilRangeGraph *frg) {
    
    double oldfgprob = currentFossilRangeGraphLnL;
    double lpr = 0.0;
    double oldPsi = fossilRate;
    double newPsi;
    double tuning = log(2.0);
    double minV = 0.0001;
    double maxV = 100000;
    double c = getNewValScaleMv(newPsi, oldPsi, minV, maxV, tuning);
    fossilRate = newPsi;
    lpr = c;
    double newfgprob = getLnFossilRangeGraphProb(frg);
    double lnPriorRat = getExpPriorRatio(oldPsi, newPsi, fossilSamplingRateExpRate, fossilSamplingRatePrior);
    double lnLikeRat = (newfgprob - oldfgprob);
    double lnR = lnLikeRat + lpr + lnPriorRat;
    double r = modelPtr->safeExponentiation(lnR);
    
    if(ranPtr->uniformRv() < r){
        currentFossilRangeGraphLnL = newfgprob;
    }
    else{
        fossilRate = oldPsi;
        setAllBDFossParams();
        currentFossilRangeGraphLnL = oldfgprob;
    }
    return 0.0;
}

double Speciation::updateDeathRate(FossilGraph *fg) {
    
    double oldfgprob = currentFossilGraphLnL;
    double lpr = 0.0;
    double oldMu = deathRate;
    double newMu;
    double tuning = log(2.0);
    double minV = 0.0001;
    double c = getNewValScaleMv(newMu, oldMu, minV, 100.0, tuning);
    deathRate = newMu;
    lpr = c;
    double newfgprob = getLnFossilGraphProb(fg);
    double lnPriorRat = getExpPriorRatio(oldMu, newMu, deathRateExpRate, deathRatePrior);
    double lnLikeRat = (newfgprob - oldfgprob);
    double lnR = lnLikeRat + lpr + lnPriorRat;
    double r = modelPtr->safeExponentiation(lnR);
    
    if(ranPtr->uniformRv() < r){
        currentFossilGraphLnL = newfgprob;
    }
    else{
        deathRate = oldMu;
        setAllBDFossParams();
        currentFossilGraphLnL = oldfgprob;
    }
    return 0.0;
}

double Speciation::updateDeathRate(FossilRangeGraph *frg) {
    
    double oldfgprob = currentFossilRangeGraphLnL;
    double lpr = 0.0;
    double oldMu = deathRate;
    double newMu;
    double tuning = log(2.0);
    double minV = 0.0001;
    double c = getNewValScaleMv(newMu, oldMu, minV, maxdivV, tuning);
    deathRate = newMu;
    lpr = c;
    double newfgprob = getLnFossilRangeGraphProb(frg);
    double lnPriorRat = getExpPriorRatio(oldMu, newMu, deathRateExpRate, deathRatePrior);
    double lnLikeRat = (newfgprob - oldfgprob);
    double lnR = lnLikeRat + lpr + lnPriorRat;
    double r = modelPtr->safeExponentiation(lnR);
    
    if(ranPtr->uniformRv() < r){
        currentFossilRangeGraphLnL = newfgprob;
    }
    else{
        deathRate = oldMu;
        setAllBDFossParams();
        currentFossilRangeGraphLnL = oldfgprob;
    }
    return 0.0;
}

double Speciation::updateDeathRate(Tree *t) {
    
    double oldtreeprob = getLnTreeProb(t);
    double lpr = 0.0;
    double oldMu = deathRate;
    double newMu;
    double tuning = log(2.0);
    double minV = 0.0001;
    double c = getNewValScaleMv(newMu, oldMu, minV, 100.0, tuning);
    deathRate = newMu;
    lpr = c;
    double newtreeprob = getLnTreeProb(t);
    double lnPriorRat = getExpPriorRatio(oldMu, newMu, deathRateExpRate, deathRatePrior);
    double lnLikeRat = (newtreeprob - oldtreeprob);
    double lnR = lnLikeRat + lpr + lnPriorRat;
    double r = modelPtr->safeExponentiation(lnR);
    if(ranPtr->uniformRv() < r){
        setAllBDFossParams();
    }
    else{
        deathRate = oldMu;
        setAllBDFossParams();
    }
    
    return 0.0;
}

double Speciation::updateBirthRate(FossilGraph *fg) {
    
    double oldfgprob = currentFossilGraphLnL;
    double lpr = 0.0;
    double oldLambda = birthRate;
    double newLambda;
    double tuning = log(3.0);
    double minV = 0.0001;
    double c = getNewValScaleMv(newLambda, oldLambda, minV, maxdivV, tuning);
    birthRate = newLambda;
    lpr = c;
    double newfgprob = getLnFossilGraphProb(fg);
    double lnPriorRat = getExpPriorRatio(oldLambda, newLambda, birthRateExpRate, birthRatePrior);
    double lnLikeRat = (newfgprob - oldfgprob);
    double lnR = lnLikeRat + lpr + lnPriorRat;
    double r = modelPtr->safeExponentiation(lnR);
    
    if(ranPtr->uniformRv() < r){
        currentFossilGraphLnL = newfgprob;
    }
    else{
        birthRate = oldLambda;
        setAllBDFossParams();
        currentFossilGraphLnL = oldfgprob;
    }
    return 0.0;
}

double Speciation::updateBirthRate(FossilRangeGraph *frg) {
    
    double oldfgprob = currentFossilRangeGraphLnL;
    double lpr = 0.0;
    double oldLambda = birthRate;
    double newLambda;
    double tuning = log(3.0);
    double minV = 0.0001;
    double c = getNewValScaleMv(newLambda, oldLambda, minV, maxdivV, tuning);
    birthRate = newLambda;
    lpr = c;
    double newfgprob = getLnFossilRangeGraphProb(frg);
    double lnPriorRat = getExpPriorRatio(oldLambda, newLambda, birthRateExpRate, birthRatePrior);
    double lnLikeRat = (newfgprob - oldfgprob);
    double lnR = lnLikeRat + lpr + lnPriorRat;
    double r = modelPtr->safeExponentiation(lnR);
    
    if(ranPtr->uniformRv() < r){
        currentFossilRangeGraphLnL = newfgprob;
    }
    else{
        birthRate = oldLambda;
        setAllBDFossParams();
        currentFossilRangeGraphLnL = oldfgprob;
    }
    return 0.0;
}

double Speciation::updateBirthAndDeath(FossilGraph *fg) {
    
    double oldfgprob = currentFossilGraphLnL;
    double lpr = 0.0;
    double newLambda, newMu;
    double oldLambda = birthRate;
	double oldMu = deathRate;
    double scaleFactor = 0.5;
    getNewValUpDownScaleMv(newLambda, oldLambda, newMu, oldMu, scaleFactor);
    birthRate = newLambda;
	deathRate = newMu;
    double newfgprob = getLnFossilGraphProb(fg);
    double lnPriorRatL = getExpPriorRatio(oldLambda, newLambda, birthRateExpRate, birthRatePrior);
    double lnPriorRatM = getExpPriorRatio(oldMu, newMu, deathRateExpRate, deathRatePrior);
    double lnLikeRat = (newfgprob - oldfgprob);
    double lnR = lnLikeRat + lpr + lnPriorRatL + lnPriorRatM;
    double r = modelPtr->safeExponentiation(lnR);
    
    if(ranPtr->uniformRv() < r){
        currentFossilGraphLnL = newfgprob;
    }
    else{
        birthRate = oldLambda;
		deathRate = oldMu;
        setAllBDFossParams();
        currentFossilGraphLnL = oldfgprob;
    }
    return 0.0;
}

double Speciation::updateBirthRate(Tree *t) {
    
    double oldtreeprob = getLnTreeProb(t);
    double lpr = 0.0;
    double oldLambda = birthRate;
    double newLambda;
    double tuning = log(2.0);
    double minV = 0.0001;
    double c = getNewValScaleMv(newLambda, oldLambda, minV, maxdivV, tuning);
    birthRate = newLambda;
    lpr = c;
    double newtreeprob = getLnTreeProb(t);
    double lnPriorRat = getExpPriorRatio(oldLambda, newLambda, birthRateExpRate, birthRatePrior);
    double lnLikeRat = (newtreeprob - oldtreeprob);
    double lnR = lnLikeRat + lpr + lnPriorRat;
    double r = modelPtr->safeExponentiation(lnR);
    if(ranPtr->uniformRv() < r){
        setAllBDFossParams();
    }
    else{
        birthRate = oldLambda;
        setAllBDFossParams();
    }
    return 0.0;
}

double Speciation::getExpPriorRatio(double oldVal, double newVal, double rate, double prior) {
    if (prior == 2)
        return (ranPtr->lnExponentialPdf(rate, newVal)) - (ranPtr->lnExponentialPdf(rate, oldVal));
    else
        return 0.0;
}

double Speciation::lnPrior(void) {
	
	return 0.0;
}

string Speciation::writeParam(void){
	
	stringstream ss;
	ss << "Speciation parameters: m/l = " << fixed << setprecision(4) << relativeDeath 
	   << " , l-m = " << fixed << setprecision(4) << netDiversificaton
	   << " , s = " << fixed << setprecision(4) << probSpeciationS
	   << " , psi = " << fixed << setprecision(4) << fossilRate
	   << " , l = " << fixed << setprecision(4) << birthRate
	   << " , m = " << fixed << setprecision(4) << deathRate
	   
	   << endl;
	return ss.str();
}

double Speciation::getLnTreeProb(Tree *t) {
	
	if(treeTimePrior == 1)
		return 0.0;
	else if(treeTimePrior == 2){
		int ntax = t->getNumTaxa();
		double c1 = (ntax - 1) * log(netDiversificaton); //  ntax * log(1 - rel) = 0;
		double nps = t->getTreeCBDNodePriorProb(netDiversificaton, relativeDeath);
		return c1 + nps;
	}
	else if(treeTimePrior == 3){
		int ntax = t->getNumTaxa();
		double c1 = (ntax - 1) * log(netDiversificaton) + ntax * log(1 - relativeDeath);
		double nps = t->getTreeCBDNodePriorProb(netDiversificaton, relativeDeath);
		return c1 + nps;
	}
	else if(treeTimePrior == 4 || treeTimePrior == 5){
		Treescale *ts = modelPtr->getActiveTreeScale();
		double nps = t->getTreeBDSSTreeNodePriorProb(netDiversificaton, relativeDeath, fossilRate, extantSampleRate, ts->getTreeOriginTime());
		return nps;
	}
	else if(treeTimePrior == 6){
		setAllBDFossParams();
		double nps = t->getTreeCalBDSSTreeNodePriorProb(birthRate, deathRate, fossilRate, extantSampleRate);
		return nps;
	}
	else if(treeTimePrior == 7){ // FBD conditioned on the root
		setAllBDFossParams();
		double nps = t->getTreeAncCalBDSSTreeNodePriorProb(birthRate, deathRate, fossilRate, extantSampleRate);
		return nps;
	}
    else if(treeTimePrior == 8){ // FBD conditioned on the origin time
        setAllBDFossParams();
        double nps = t->getTreeStemAncCalBDSSTreeNodePriorProb(birthRate, deathRate, fossilRate, extantSampleRate);
        return nps;
    }
	return 0.0;
}

double Speciation::getLnFossilGraphProb(FossilGraph *fg) {
	
	setAllBDFossParams();
	double fgprob = fg->getFossilGraphProb(birthRate, deathRate, fossilRate, extantSampleRate);
	return fgprob;
}

double Speciation::getLnFossilRangeGraphProb(FossilRangeGraph *frg) {
    
    double ot = frg->getFossilRangeGraphOriginTime();
    setAllBDFossParams();
    double fgprob = frg->getFossilRangeGraphProb(birthRate, deathRate, fossilRate, extantSampleRate, ot);
    return fgprob;
}

void Speciation::setAllBDFossParams(){
	
	if(parameterization == 1){
		fossilRate = (probSpeciationS / (1-probSpeciationS)) * ((relativeDeath * netDiversificaton) / (1 - relativeDeath));
		birthRate = netDiversificaton / (1.0 - relativeDeath); 
		deathRate = (relativeDeath * netDiversificaton) / (1 - relativeDeath);
	}
	else if(parameterization == 2){
		birthRate = netDiversificaton / (1.0 - relativeDeath);
		deathRate = (relativeDeath * netDiversificaton) / (1 - relativeDeath);
		probSpeciationS = fossilRate / (deathRate + fossilRate);
	}
	else if(parameterization == 3){
		netDiversificaton = birthRate - deathRate;
		relativeDeath = deathRate / birthRate;
		probSpeciationS = fossilRate / (deathRate + fossilRate);
	}
    else if(parameterization == 4){
        birthRate = netDiversificaton + deathRate;
        relativeDeath = deathRate / birthRate;
        probSpeciationS = fossilRate / (deathRate + fossilRate);
    }
}
