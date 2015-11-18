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
#include "Parameter_fossilrangegraph.h"
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

FossilRangeGraph::FossilRangeGraph(MbRandom *rp, Model *mp, int nf, int nl, vector<Calibration *> clb, bool rnp) : Parameter(rp, mp){
    
    name = "FRG";
    numFossils = nf;
    numLineages = nl;
    runUnderPrior = rnp;
    originTime = 0.0;
    ancientBound = 1000.0;
    printInitialFossilRangeVariables = 1; //rw: this is for debugging
    createFossilRangeVector(clb);
    initializeFossilRangeVariables();
    moves = 1; // 1: update lineage start or stop times; 2: update both
    proposal = 2; // proposal type 1=window, 2=scale, 3=slide
    
    //numAncFossilsk=0;
}

FossilRangeGraph::~FossilRangeGraph(void){
    
}

FossilRangeGraph& FossilRangeGraph::operator=(const FossilRangeGraph &frg) {
    
    if (this != &frg)
        clone(frg);
    return *this;
}

void FossilRangeGraph::clone(const FossilRangeGraph &frg){
    
    numFossils = frg.numFossils;
    numLineages = frg.numLineages;
    
    for(int i=0; i<fossilRanges.size(); i++){
        FossilRange *frTo = fossilRanges[i];
        FossilRange *frFrom = frg.fossilRanges[i];
        frTo->setFossilRangeIndex(frFrom->getFossilRangeIndex());
        frTo->setFirstAppearance(frFrom->getFirstAppearance());
        frTo->setLastAppearance(frFrom->getLastAppearance());
        frTo->setExtant(frFrom->getIsExtant());
        frTo->setExtantOnly(frFrom->getIsExtantOnly());
        frTo->setFossilRangeID(frFrom->getFossilRangeID());
        frTo->setFossilRangeBrGamma(frFrom->getFossilRangeBrGamma());
        frTo->setLineageStart(frFrom->getLineageStart());
        frTo->setLineageStop(frFrom->getLineageStop());
    }
    
    //numAncFossilsk = frg.numAncFossilsk;
    //treeTimePrior = frg.treeTimePrior;
    
}

double FossilRangeGraph::update(double &oldLnL){
    
    if(moves == 1){
        if(ranPtr->uniformRv() < 0.5)
            updateLineageStartTimes();
        //else
          //updateLineageStopTimes();
    }
    else if (moves == 2) {
        updateLineageStartTimes();
        //updateLineageStopTimes();
    }
    
    return currentFossilRangeGraphLnL;
    
}

double FossilRangeGraph::lnPrior(){
    
    return 0;
}

void FossilRangeGraph::print(ostream & o) const {
    
}

string FossilRangeGraph::writeParam(void){

    string s="";
    return s;
}

void FossilRangeGraph::createFossilRangeVector(vector<Calibration *> clb){
    
    //double et = originTime;
    //terminalTime = et;
    
    int frid=1;
    for(int c = 0; c < clb.size(); c++){
        Calibration *p = clb[c];
        double fa = p->getFirstAppearance();
        double la = p->getLastAppearance();
        bool e = p->getIsExtant();
        bool eo = p->getIsExtantOnly();
        
        FossilRange *fr = new FossilRange(fa, la, e, eo, frid);
        fossilRanges.push_back(fr);
        
        frid ++;
    }
}


void FossilRangeGraph::initializeFossilRangeVariables(){
    
    //numAncFossilsk = 0; //rw: do we ever need to know this for the frg?
    double stop, start, la, fa;
    
    for(int f = 0; f < numLineages; f++){
        FossilRange *fr = fossilRanges[f];
        
        // lineage start times = speciation times
        if(fr->getIsExtantOnly()){ //rw: need to check this with T&T
            start = ranPtr->uniformRv(0.0,ancientBound);
            fr->setLineageStart(start);
        }
        else{
            fa = fr->getFirstAppearance(); // yf
            start = ranPtr->uniformRv(fa,ancientBound); // zf
            fr->setLineageStart(start);
        }
        
        // lineage stop times = extinction times
        if(fr->getIsExtant()){
            fr->setLineageStop(0.0);
        }
        else {
            la = fr->getLastAppearance();
            stop = ranPtr->uniformRv(0.0,la);
            fr->setLineageStop(stop);
        }
    }
    
    redefineOriginTime();
    recountFossilRangeAttachNums();
    cout << "\nInitial origin time: " << originTime << endl;
}

void FossilRangeGraph::recountFossilRangeAttachNums(){
    
    for(int f = 0; f < numLineages; f++){
        FossilRange *fr = fossilRanges[f];
        double zf = fr->getLineageStart();
        int g = 0; //rw: should this be 1? cf fossilgraph
        for(int j = 0; j < numLineages; j++){
            if(f != j){
                FossilRange *p = fossilRanges[j];
                double zj = p->getLineageStart();
                double yj = p->getLineageStop(); //rw: check this is correct
                if(zj > zf && zf > yj)
                    g++;
            }
        }
        
        fr->setFossilRangeBrGamma(g);
    }
    
    //rw: potentially place this in a seperate fxn so you can monitor the mcmc
    if(printInitialFossilRangeVariables) {
        cout << "Debugging stuff...\nprintInitialFossilRangeVariables:\n" << endl;
        for(int f = 0; f < numLineages; f++){
            FossilRange *fr = fossilRanges[f];
            cout << "Fossil range ID :" << fr->getFossilRangeID() << endl;
            cout << "First appearance :" << fr->getFirstAppearance() << endl;
            cout << "Last appearance :" << fr->getLastAppearance() << endl;
            cout << "Lineage start :" << fr->getLineageStart() << endl;
            cout << "Lineage end :" << fr->getLineageStop() << endl;
            cout << "Gamma :" << fr->getFossilRangeBrGamma() << endl; //rw: is gamma correct?
        }
    }
}

void FossilRangeGraph::redefineOriginTime(){
    
    double ot = originTime;
    
    for(int f = 0; f < numLineages; f++){
        FossilRange *fr = fossilRanges[f];
        double ls = fr->getLineageStart();
        if(ls > ot)
            ot = ls;
    }
    
    originTime = ot;
}

double FossilRangeGraph::updateLineageStartTimes(){
    
    Speciation *s = modelPtr->getActiveSpeciation();
    s->setAllBDFossParams();
    double lambda = s->getBDSSSpeciationRateLambda();
    double mu = s->getBDSSExtinctionRateMu();
    double fossRate = s->getBDSSFossilSampRatePsi();
    double sppSampRate = s->getBDSSSppSampRateRho();
    
    // maybe unnecessary here, for debugging
    //getFossilGraphProb(lambda, mu, fossRate, sppSampRate, originTime);
    
    vector<int> rndFossilRangeIDs;
    for(int i=0; i<fossilRanges.size(); i++)
        rndFossilRangeIDs.push_back(i);
    random_shuffle(rndFossilRangeIDs.begin(), rndFossilRangeIDs.end());
    
    for(vector<int>::iterator it=rndFossilRangeIDs.begin(); it!=rndFossilRangeIDs.end(); it++){
        
        FossilRange *fr = fossilRanges[(*it)];
        
        // define old values
        double fa = fr->getFirstAppearance(); // yf
        double oldStart = fr->getLineageStart(); // zf
        double oldLike = currentFossilRangeGraphLnL; //rw: where is this initially defined?

        // propose new values
        double newStart = 0.0;
        double c = 0.0;
        
        if(proposal == 1){ //rw: probably not appropriate for this parameter
            double rdwindow = 0.2;
            newStart = getNewValSWindoMv(oldStart, fa, ancientBound, rdwindow);
        }
        
        if(proposal == 2){
            double rv = ranPtr->uniformRv();
            double tv = 2.0; //rw: tuningVal
            c = doAScaleMove(newStart, oldStart, tv, fa, ancientBound, rv);
        }
        
        // redefine values
        fr->setLineageStart(newStart);
        redefineOriginTime();
        recountFossilRangeAttachNums();
        
        // recalculate the FRG probability
        double newLike = getFossilRangeGraphProb(lambda, mu, fossRate, sppSampRate, originTime);
        
        // calculate the likelihood/prior ratio
        double lnLikeRatio = newLike - oldLike;
        double r = modelPtr->safeExponentiation(lnLikeRatio + c);
            
        if(ranPtr->uniformRv() > r){
            fr->setLineageStart(oldStart);
            redefineOriginTime();
            recountFossilRangeAttachNums();
            currentFossilRangeGraphLnL = oldLike;
        }
    }
    return 0.0;
}

double FossilRangeGraph::updateLineageStopTimes(){
 
    return 0.0;
}

double FossilRangeGraph::getFossilRangeGraphProb(double lambda, double mu, double fossRate, double sppSampRate, double ot){
    if(runUnderPrior)
        return 0.0;
        
    currentFossilRangeGraphLnL = 0.0;
    return 0.0;
}

double FossilRangeGraph::getNewValSWindoMv(double ov, double vmin, double vmax, double tv){
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
    } while(!validV);
    return nv;
}

double FossilRangeGraph::doAScaleMove(double &nv, double cv, double tv, double lb, double hb, double rv){
    
    /*
     nv = new start
     cv = old start
     tv = tuning value
     lb = first appearance
     hb = ancient bound  (previously, node depth or origin time)
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

double FossilRangeGraph::getActiveFossilRangeGraphProb(){
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
    
    nprb = getFossilRangeGraphProb(lambda, mu, fossRate, sppSampRate, originTime);
    return nprb;
    
}

//END