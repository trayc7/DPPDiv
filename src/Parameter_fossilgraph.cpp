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

#include "Calibration.h"
#include "Parameter.h"
#include "Parameter_fossilgraph.h"
#include "Parameter_origin.h"
#include "Parameter_speciaton.h"
#include "MbRandom.h"
#include "Model.h"
#include "util.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <algorithm>

using namespace std;

FossilGraph::FossilGraph(MbRandom *rp, Model *mp, int nf, double initOrigTime, vector<Calibration *> clb, bool rnp) : Parameter(rp, mp){

    name = "FG";
    numFossils = nf;
    originTime = initOrigTime;
    createOccurrenceVector(clb);
    numAncFossilsk=0;
	runUnderPrior = rnp;
	tuningVal = 2.0;
    initializeOccurrenceSpecVariables();
}

FossilGraph::~FossilGraph(void){
    
}

FossilGraph& FossilGraph::operator=(const FossilGraph &fg) {
    
    if (this != &fg)
        clone(fg);
    return *this;
}

void FossilGraph::clone(const FossilGraph &fg){

    numFossils = fg.numFossils;
    
    /* any error checking neccessary here?
    if (numNodes != t.numNodes || numTaxa != t.numTaxa){
        cerr << "ERROR: Attempting to clone trees of unequal size" << endl;
        exit(1);
    */
        
    for(int i=0; i<occurrenceSpecimens.size(); i++){
        Occurrence *fTo = occurrenceSpecimens[i];
        Occurrence *fFrom = fg.occurrenceSpecimens[i];
        fTo->setFossilIndex(fFrom->getFossilIndex());
        fTo->setFossilAge(fFrom->getFossilAge());
        fTo->setFossilSppTime(fFrom->getFossilSppTime());
        fTo->setFossilID(fFrom->getFossilID());
        fTo->setFossilFossBrGamma(fFrom->getFossilFossBrGamma());
        fTo->setFossilIndicatorVar(fFrom->getFossilIndicatorVar());
        fTo->setIsTerminal(fFrom->getIsTerminal());
    }
    
    numAncFossilsk = fg.numAncFossilsk;
    treeTimePrior = fg.treeTimePrior;
    
        
    // TAH: double check this (more bookkeeping is probably needed, for now this is a placeholder)
    OriginTime *ot = modelPtr->getActiveOriginTime();
    originTime = ot->getOriginTime();

}

double FossilGraph::update(double &oldLnL){
    
    OriginTime *ot = modelPtr->getActiveOriginTime();
    originTime = ot->getOriginTime();
    
    //double lppr = 0.0;
	if(ranPtr->uniformRv() < 0.5)
		updateOccurrenceAttachmentTimesPhi();
	else{
		for(int i=0; i < numFossils; i++){
			updateRJMoveAddDelEdge();
		}
	}
    return currentFossilGraphLnL;
}

double FossilGraph::lnPrior(){

    return 0;
}

void FossilGraph::print(ostream & o) const {

}

string FossilGraph::writeParam(void){

    string s="";
    return s;
}


void FossilGraph::createOccurrenceVector(vector<Calibration *> clb){
    
    double et = originTime;
    int oid=1;
    for(int c = 0; c < clb.size(); c++){
        Calibration *p = clb[c];
        double yf = p->getYngTime();
        Occurrence *o = new Occurrence(yf, oid);
        occurrenceSpecimens.push_back(o);
        if( yf < et)
            et = yf;
            oid ++;
    }
    terminalTime = et;
    for(int c = 0; c < numFossils; c++){
        Occurrence *o = occurrenceSpecimens[c];
        if( o->getFossilAge() == terminalTime){
            o->setIsTerminal(true);
            break; // currently no two youngest occurrences can have the same age
        }
    }
}

void FossilGraph::initializeOccurrenceSpecVariables(){
    
    numAncFossilsk = 0;
    for(int f = 0; f < numFossils; f++){
        Occurrence *o = occurrenceSpecimens[f];
        if(o->getIsTerminal()){
            numAncFossilsk += 1;
            o->setFossilSppTime(o->getFossilAge());
            o->setFossilIndicatorVar(0);
        }
        else{
            double u = ranPtr->uniformRv();
            if(u < 0.5){
                numAncFossilsk += 1;
                o->setFossilSppTime(o->getFossilAge());
                o->setFossilIndicatorVar(0);
            }
            else {
                double yf = o->getFossilAge();
                double zf = ranPtr->uniformRv(yf,originTime);
                o->setFossilSppTime(zf);
                o->setFossilIndicatorVar(1);
                o->setFossilFossBrGamma(0);
            }
        }
    }
    recountOccurrenceAttachNums();
    cout << "1. num anc Fossils = " << numAncFossilsk << endl;
    
}

void FossilGraph::recountOccurrenceAttachNums(){

    for(int f = 0; f < numFossils; f++){
        Occurrence *o = occurrenceSpecimens[f];
        double zf = o->getFossilSppTime();
        int g = 1; // start at 1 because the lineage leading to OT is implicit
        for(int j = 0; j < numFossils; j++){
            if(f != j){
                Occurrence *p = occurrenceSpecimens[j];
                double zj = p->getFossilSppTime();
                double yj = p->getFossilAge();
                if(zj > zf && zf > yj)
                    g++;
            }
        }
        //cout << "gamma :" << g << endl; // gamma is correct
        o->setFossilFossBrGamma(g);
    }
}

double FossilGraph::getOldestFossilGraphSpeciationTime(){
    
    double oldestTime = 0.0;
    
    for(int f = 0; f < numFossils; f++){
        Occurrence *o = occurrenceSpecimens[f];
        double zf = o->getFossilSppTime();
        if(zf > oldestTime)
            oldestTime = zf;
    }
    
    return oldestTime;
}

double FossilGraph::getActiveFossilGraphProb(){
    if(runUnderPrior)
		return 0.0;

    double nprb = 0.0;
    Speciation *s = modelPtr->getActiveSpeciation();
    s->setAllBDFossParams();
    double lambda = s->getBDSSSpeciationRateLambda();
    double mu = s->getBDSSExtinctionRateMu();
    double fossRate = s->getBDSSFossilSampRatePsi();
    double sppSampRate = s->getBDSSSppSampRateRho();
//    recountOccurrenceAttachNums(); // debug
    nprb = getFossilGraphProb(lambda, mu, fossRate, sppSampRate);
    return nprb;
    
}


double FossilGraph::getFossilGraphProb(double lambda, double mu, double fossRate, double sppSampRate) {
	if(runUnderPrior)
		return 0.0;
    OriginTime *ot = modelPtr->getActiveOriginTime();
    double nprb = getFossilGraphProb(lambda, mu, fossRate, sppSampRate, ot->getOriginTime());
    return nprb;
}

double FossilGraph::getFossilGraphProb(double lambda, double mu, double fossRate, double sppSampRate, double ot) {
	if(runUnderPrior)
		return 0.0;
    originTime = ot;
    
    // the following has been modified for the fofbd
    double nprb = -(log(2*lambda) + log(1.0 - fbdPFxn(lambda, mu, fossRate, sppSampRate, originTime)));
    
    for(int f=0; f < occurrenceSpecimens.size(); f++){
        Occurrence *o = occurrenceSpecimens[f];
		if(!o->getIsTerminal()){
			nprb += log( fossRate * o->getFossilFossBrGamma() );
			if(o->getFossilIndicatorVar()){
				double yf = o->getFossilAge();
				double zf = o->getFossilSppTime();
				double fossPr = log(2.0 * lambda) + log( fbdPFxn(lambda, mu, fossRate, sppSampRate, yf) );
				fossPr += fbdQFxn(lambda, mu, fossRate, sppSampRate, yf);
				fossPr -= fbdQFxn(lambda, mu, fossRate, sppSampRate, zf);

				nprb += fossPr;
			}
		}
    }
	currentFossilGraphLnL = nprb;
    return nprb;
}

double FossilGraph::fbdPFxn(double b, double d, double psi, double rho, double t){
    
    double c1Val = fbdC1Fxn(b,d,psi);
    double c2Val = fbdC2Fxn(b,d,psi,rho);
	
	double expC1MinusC2 = exp(-c1Val * t) * (1.0 - c2Val);
    
    double eCfrac = (expC1MinusC2 - (1.0 + c2Val)) / (expC1MinusC2 + (1.0 + c2Val));
    double v = 1.0 + ((-(b - d - psi)) + (c1Val * eCfrac)) / (2.0 * b);
    
    return v;
}

double FossilGraph::fbdQHatFxn(double b, double d, double psi, double rho, double t){
    
    double v = log(4.0);
    v -= fbdQFxn(b,d,psi,rho,t);
    return v;
}

double FossilGraph::fbdQFxn(double b, double d, double psi, double rho, double t){
    
    double c1Val = fbdC1Fxn(b,d,psi);
    double c2Val = fbdC2Fxn(b,d,psi,rho);
    
    double vX = c1Val * t + 2.0 * log(exp(-c1Val * t) * (1.0 - c2Val) + (1.0 + c2Val));
    
    // returns log(q)
    return vX;
}

double FossilGraph::fbdC1Fxn(double b, double d, double psi){
    
    double v = fabs( sqrt( ( (b-d-psi) * (b-d-psi) ) + 4*b*psi) );
    return v;
}

double FossilGraph::fbdC2Fxn(double b, double d, double psi, double rho){
    
    double v = -( ( b-d-(2*b*rho)-psi ) / (fbdC1Fxn(b,d,psi)) );
    return v;
}

// the following functions are for sampling/printing from the mcmc
int FossilGraph::getSumIndicatorFG(){
    int niv = 0;
    for(int i=0; i<occurrenceSpecimens.size(); i++){
        Occurrence *o = occurrenceSpecimens[i];
        niv += o->getFossilIndicatorVar();
    }
    return niv;
}

string FossilGraph::getOccInfoParamNames(void){
    
    stringstream ss;
    Occurrence *o = NULL;
    for(int i=0; i<occurrenceSpecimens.size(); i++){
        o = occurrenceSpecimens[i];
        int oID = o->getFossilID();
        // we probably don't want it to print the fossil age forever
        ss << "\ty_f(Oc_" << oID << ")";
        ss << "\tz_f(Oc_" << oID << ")";
    }
    for(int i=0; i<occurrenceSpecimens.size(); i++){
        o = occurrenceSpecimens[i];
        int oID = o->getFossilID();
        ss << "\tgamma_f(Oc_" << oID << ")"; //".nd" << nID << ")";
    }
    string ni = ss.str();
    return ni;
}

string FossilGraph::getOccInfoParamList(void){
    
    stringstream ss;
    Occurrence *o = NULL;

    for(int i=0; i<occurrenceSpecimens.size(); i++){
        o = occurrenceSpecimens[i];
        // we probably don't want it to print the fossil age forever
        ss << "\t" << o->getFossilAge();
        ss << "\t" << o->getFossilSppTime();
    }
    for(int i=0; i<occurrenceSpecimens.size(); i++){
        o = occurrenceSpecimens[i];
        ss << "\t" << o->getFossilFossBrGamma();
    }
    string ni = ss.str();
    return ni;
}


// the following functions are all to do with update moves

double FossilGraph::updateOccurrenceAttachmentTimesPhi() {
    
    OriginTime *ot = modelPtr->getActiveOriginTime();
    originTime = ot->getOriginTime(); // active origin?
    
    Speciation *s = modelPtr->getActiveSpeciation();
    s->setAllBDFossParams();
    double lambda = s->getBDSSSpeciationRateLambda();
    double mu = s->getBDSSExtinctionRateMu();
    double fossRate = s->getBDSSFossilSampRatePsi();
    double sppSampRate = s->getBDSSSppSampRateRho();
	
	// maybe unnecessary here for debugging
	getFossilGraphProb(lambda, mu, fossRate, sppSampRate, originTime);
    
    vector<int> rndOccIDs;
    for(int i=0; i<occurrenceSpecimens.size(); i++)
        rndOccIDs.push_back(i);
    random_shuffle(rndOccIDs.begin(), rndOccIDs.end());
    for(vector<int>::iterator it=rndOccIDs.begin(); it!=rndOccIDs.end(); it++){
        Occurrence *o = occurrenceSpecimens[(*it)];
        if(o->getFossilIndicatorVar()){
            
            double fossDepth = o->getFossilAge();
            
            double oldPhi = o->getFossilSppTime(); // * treeScale;
			
			double oldLike = currentFossilGraphLnL;
            
            double rv = ranPtr->uniformRv();
            double newPhi, c;
            
            /*
            if(nodeProposal == 1){ // used for options < 5
                double delta = 5.0;
                c = doAWindoMove(newPhi, oldPhi, delta, fossDepth, nodeDepth, rv);
            }
            else if(nodeProposal == 2){ // used for options > 5
                double tv = tuningVal;
                c = doAScaleMove(newPhi, oldPhi, tv, fossDepth, nodeDepth, rv);
            }
            else if(nodeProposal == 3){ // never used
                newPhi = fossDepth + rv*(nodeDepth-fossDepth);
                c = 0.0;
            }
             */
            
            // we should spectify the following in the FG initialization
            int nodeProposal = 2; // proposal type 1=window, 2=scale, 3=slide
            
            if(nodeProposal == 2){
                double tv = tuningVal;
                c = doAScaleMove(newPhi, oldPhi, tv, fossDepth, originTime, rv);
            }
            
            o->setFossilSppTime(newPhi);
            
			recountOccurrenceAttachNums();
			double newLike = getFossilGraphProb(lambda, mu, fossRate, sppSampRate, originTime);
			double lnLikeRat2 = newLike - oldLike;
            
            double r = modelPtr->safeExponentiation(lnLikeRat2 + c);
            
            if(ranPtr->uniformRv() > r){
                o->setFossilSppTime(oldPhi);
				recountOccurrenceAttachNums();
				currentFossilGraphLnL = oldLike;
            }
        }
    }
    return 0.0;
}

double FossilGraph::doAScaleMove(double &nv, double cv, double tv, double lb, double hb, double rv){
    
    /* 
    nv = new phi
    cv = old phi
    tv = tuning value
    lb = fossil depth
    hb = origin time (previously, node depth)
    rv = random number
    */
    
    double c = tv * (rv - 0.5);
    double newcv = cv * exp(c);
    bool validV = false;
    do{
        if(newcv < lb)
            newcv = lb * lb / newcv;
        else if(newcv > hb)
            newcv = hb * hb / newcv;
        else
            validV = true;
    } while(!validV);
    nv = newcv;
    return c;
}

double FossilGraph::getSumLogAllAttachNums(){
    
    double sumLog = 0.0;
    recountOccurrenceAttachNums();
    for(int i=0; i<occurrenceSpecimens.size(); i++){
        Occurrence *fi = occurrenceSpecimens[i];
        sumLog += log(fi->getFossilFossBrGamma());
    }
    return sumLog;
}

double FossilGraph::updateRJMoveAddDelEdge() {
	int nAncs = numFossils - getSumIndicatorFG() - 1;
	int nf = numFossils - 1;
    
    double gA = 0.5;
    if(nAncs == nf)
        gA = 1.0;
    else if(nAncs == 0)
        gA = 0.0;
    
    double u = ranPtr->uniformRv();
    if(u < gA)
        doAddEdgeMove(nAncs);
    else
        doDeleteEdgeMove(nAncs);
    
    return 0.0;
}

// double check the maths
void FossilGraph::doAddEdgeMove(int k){
    
	int b = k;
    int w = numFossils-1;
    int bN = b-1;
    Speciation *s = modelPtr->getActiveSpeciation();
    s->setAllBDFossParams();
    double lambda = s->getBDSSSpeciationRateLambda();
    double mu = s->getBDSSExtinctionRateMu();
    double fossRate = s->getBDSSFossilSampRatePsi();
    double sppSampRate = s->getBDSSSppSampRateRho();
    double lnHastings, lnJacobian; //, lnPriorR;
    
    int mvFoss = pickRandAncestorFossil();
    Occurrence *o = occurrenceSpecimens[mvFoss]; // rand fossil
	
    double alA = 1.0;
    if(bN == 0)
        alA = 2.0;
    else if(b == w)
        alA = 0.5;
    
    OriginTime *ot = modelPtr->getActiveOriginTime();
    double cf = ot->getOriginTime();
	double oldLnl = getFossilGraphProb(lambda, mu, fossRate, sppSampRate, cf);
    
    
    double yf = o->getFossilAge();
    double nu = ranPtr->uniformRv() * (cf - yf);
    lnHastings = log(alA) + (log(b) - log(w - b + 1.0));
    lnJacobian = log(cf - yf);
    double newPhi = yf + nu;
    o->setFossilSppTime(newPhi);
    o->setFossilIndicatorVar(1);
    

	recountOccurrenceAttachNums();
	double newLnl = getFossilGraphProb(lambda, mu, fossRate, sppSampRate, cf);
	double lnLikeR = newLnl - oldLnl;
	


    double lpr = lnLikeR + lnHastings + lnJacobian;
    double r = modelPtr->safeExponentiation(lpr);
    
    if(ranPtr->uniformRv() > r){
        o->setFossilIndicatorVar(0);
        o->setFossilSppTime(yf);
		recountOccurrenceAttachNums();
		currentFossilGraphLnL = oldLnl;
    }
}

void FossilGraph::doDeleteEdgeMove(int k){
    
	int b = k;
    int w = numFossils-1;
    int bN = b+1;
    Speciation *s = modelPtr->getActiveSpeciation();
    s->setAllBDFossParams();
    double lambda = s->getBDSSSpeciationRateLambda();
    double mu = s->getBDSSExtinctionRateMu();
    double fossRate = s->getBDSSFossilSampRatePsi();
    double sppSampRate = s->getBDSSSppSampRateRho();
    double lnHastings, lnJacobian; //, lnPriorR;
    
    int mvFoss = pickRandTipFossil();
    Occurrence *o = occurrenceSpecimens[mvFoss];
	double oldLnl = getFossilGraphProb(lambda, mu, fossRate, sppSampRate);

    double alD = 1.0;
    if(b == 0)
        alD = 0.5;
    else if(bN == w)
        alD = 2.0;

    double yf = o->getFossilAge();
	double cf = originTime;
    
    double oldPhi = o->getFossilSppTime();
    o->setFossilIndicatorVar(0);
    o->setFossilSppTime(yf);

    
    lnHastings = log(alD) + (log(w-b) - log(b+1));
    lnJacobian = -(log(cf - yf));
    
    recountOccurrenceAttachNums();
	double newLnl = getFossilGraphProb(lambda, mu, fossRate, sppSampRate);
	double lnLikeR = newLnl - oldLnl;


    double lpr = lnLikeR + lnHastings + lnJacobian;
    double r = modelPtr->safeExponentiation(lpr);
    
    if(ranPtr->uniformRv() > r){
        o->setFossilSppTime(oldPhi);
        o->setFossilIndicatorVar(1);
		recountOccurrenceAttachNums();
		currentFossilGraphLnL = oldLnl;
    }
}

double FossilGraph::doSinglePhiMove() {
    
    OriginTime *ot = modelPtr->getActiveOriginTime();
    originTime = ot->getOriginTime(); // active origin?
    
    Speciation *s = modelPtr->getActiveSpeciation();
    s->setAllBDFossParams();
    double lambda = s->getBDSSSpeciationRateLambda();
    double mu = s->getBDSSExtinctionRateMu();
    double fossRate = s->getBDSSFossilSampRatePsi();
    double sppSampRate = s->getBDSSSppSampRateRho();
	
	int pick = (int)(ranPtr->uniformRv() * occurrenceSpecimens.size());
    Occurrence *o = occurrenceSpecimens[pick];
	while(o->getFossilIndicatorVar() == 0){
		pick = (int)(ranPtr->uniformRv() * occurrenceSpecimens.size());
		o = occurrenceSpecimens[pick];
	}
	
            
	double fossDepth = o->getFossilAge();
	
	double oldPhi = o->getFossilSppTime();
	
	double oldLike = currentFossilGraphLnL;
	
	double rv = ranPtr->uniformRv();
	double newPhi, c;
	
	/*
	if(nodeProposal == 1){ // used for options < 5
		double delta = 5.0;
		c = doAWindoMove(newPhi, oldPhi, delta, fossDepth, nodeDepth, rv);
	}
	else if(nodeProposal == 2){ // used for options > 5
		double tv = tuningVal;
		c = doAScaleMove(newPhi, oldPhi, tv, fossDepth, nodeDepth, rv);
	}
	else if(nodeProposal == 3){ // never used
		newPhi = fossDepth + rv*(nodeDepth-fossDepth);
		c = 0.0;
	}
	 */
	
	// we should spectify the following in the FG initialization
	int nodeProposal = 2; // proposal type 1=window, 2=scale, 3=slide
	
	if(nodeProposal == 2){
		double tv = tuningVal;
		c = doAScaleMove(newPhi, oldPhi, tv, fossDepth, originTime, rv);
	}
	
	o->setFossilSppTime(newPhi);
	
	recountOccurrenceAttachNums();
	double newLike = getFossilGraphProb(lambda, mu, fossRate, sppSampRate, originTime);
	double lnLikeRat = newLike - oldLike;
	
	double r = modelPtr->safeExponentiation(lnLikeRat + c);
	
	if(ranPtr->uniformRv() > r){
		o->setFossilSppTime(oldPhi);
		recountOccurrenceAttachNums();
		currentFossilGraphLnL = oldLike;
	}
     return 0.0;
}



int FossilGraph::pickRandAncestorFossil(){
    
    
    vector<int> af;
    for(int i=0; i<occurrenceSpecimens.size(); i++){
        Occurrence *f = occurrenceSpecimens[i];
		if(!f->getIsTerminal()){ // avoid the terminal fossil
				if(f->getFossilIndicatorVar() == 0)
					af.push_back(i);
		}
    }
    int v = (int)(ranPtr->uniformRv()*af.size());
    
    return af[v];
    
}

int FossilGraph::pickRandTipFossil(){

    vector<int> af;
    for(int i=0; i<occurrenceSpecimens.size(); i++){
        Occurrence *f = occurrenceSpecimens[i];
        if(f->getFossilIndicatorVar() == 1)
            af.push_back(i);
    }
    int v = (int)(ranPtr->uniformRv()*af.size());
    return af[v];
}







// END
