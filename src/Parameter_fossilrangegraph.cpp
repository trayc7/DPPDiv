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

FossilRangeGraph::FossilRangeGraph(MbRandom *rp, Model *mp, int nf, int nl, vector<Calibration *> clb, bool rnp, int fxSt, int fxSp) : Parameter(rp, mp){
    
    name = "FRG";
    numFossils = nf;
    numLineages = nl;
    runUnderPrior = rnp;
    originTime = 0.0;
    ancientBound = 1000.0;
    printInitialFossilRangeVariables = 1;
    numExtinctLineages = 0;
    fixFRG = 1; //1: fix start and end range times to FAs and LAs
    fixStart = fxSt;
    fixStop = fxSp;
    createFossilRangeVector(clb);
    initializeFossilRangeVariables();
    currentFossilRangeGraphLnL = 0.0;
    moves = 1; // 1: update lineage start or stop times; 2: update both
    proposal = 2; // proposal type 1=window, 2=scale, 3=slide
    getAltProb=0;
    
    cout << "Number of lineages: " << numLineages << endl;
    cout << "Number of extinct ranges: " << numExtinctLineages << endl;
    cout << "Number of fossils: " << numFossils << endl;
    cout << "\nInitial origin time: " << originTime << endl;
    cout << "Fossil range graph initialized" << endl;
    
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
        else
          updateLineageStopTimes();
    }
    else if (moves == 2) {
        updateLineageStartTimes();
        updateLineageStopTimes();
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
        if(fr->getIsExtant()){
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
            stop = ranPtr->uniformRv(0.0,la); //di
            fr->setLineageStop(stop);
        }
    }

    if(fixFRG){
        for(int f = 0; f < numLineages; f++){
            FossilRange *fr = fossilRanges[f];
            fr->setLineageStart(fr->getFirstAppearance());
            fr->setLineageStop(fr->getLastAppearance());
            fr->setFixStart(1);
            fr->setFixStop(1);
        }
    }
    
    if(fixStart > 0){
        for(int f = 0; f < fixStart; f++){
            FossilRange *fr = fossilRanges[f];
            fr->setLineageStart(fr->getFirstAppearance());
            fr->setFixStart(1);
        }
    }
    
    if(fixStop > 0){
        for(int f = 0; f < fixStop; f++){
            FossilRange *fr = fossilRanges[f];
            fr->setLineageStop(fr->getLastAppearance());
            fr->setFixStop(1);
        }
    }
    
    redefineOriginTime();
    recountFossilRangeAttachNums();
    countExtinctLineages();
    
    if(printInitialFossilRangeVariables) //rw: debugging code
        printFossilRangeVariables();
    
}

//rw: for debugging
void FossilRangeGraph::printFossilRangeVariables(){
    
    for(int f = 0; f < numLineages; f++){
        FossilRange *fr = fossilRanges[f];
        cout << "Fossil range ID: " << fr->getFossilRangeID() << endl;
        cout << "First appearance: " << fr->getFirstAppearance() << endl;
        cout << "Last appearance: " << fr->getLastAppearance() << endl;
        cout << "Lineage start: " << fr->getLineageStart() << endl;
        cout << "Lineage end: " << fr->getLineageStop() << endl;
        cout << "Is extant: " << fr->getIsExtant() << endl;
        cout << "Gamma: " << fr->getFossilRangeBrGamma() << endl; //rw: is gamma correct?
    }
    
    cout << "Origin time: " << originTime << endl;
    
}

void FossilRangeGraph::recountFossilRangeAttachNums(){
    
    for(int f = 0; f < numLineages; f++){
        FossilRange *fr = fossilRanges[f];
        double zf = fr->getLineageStart();
        int g = 0;
        for(int j = 0; j < numLineages; j++){
            if(f != j){
                FossilRange *p = fossilRanges[j];
                double zj = p->getLineageStart();
                double bj = p->getLineageStop();
                if(zj > zf && zf > bj)
                    g++;
            }
        }
        if(zf == originTime)
            g++;
        
        fr->setFossilRangeBrGamma(g);
    }
    
}

void FossilRangeGraph::countExtinctLineages(){
    
    int el = 0;
    
    for(int f = 0; f < numLineages; f++){
        FossilRange *fr = fossilRanges[f];
        if(!fr->getIsExtant())
            el += 1;
    }
    
    numExtinctLineages = el;
    
}

void FossilRangeGraph::redefineOriginTime(){
    
    double ot = 0.0;
    
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
    getFossilRangeGraphProb(lambda, mu, fossRate, sppSampRate, originTime);
    
    vector<int> rndFossilRangeIDs;
    for(int i=0; i<fossilRanges.size(); i++)
        rndFossilRangeIDs.push_back(i);
    random_shuffle(rndFossilRangeIDs.begin(), rndFossilRangeIDs.end());
    
    for(vector<int>::iterator it=rndFossilRangeIDs.begin(); it!=rndFossilRangeIDs.end(); it++){
        
        FossilRange *fr = fossilRanges[(*it)];
        
        if(fr->getIsFixStart())
            continue;
        
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
    
    Speciation *s = modelPtr->getActiveSpeciation();
    s->setAllBDFossParams();
    double lambda = s->getBDSSSpeciationRateLambda();
    double mu = s->getBDSSExtinctionRateMu();
    double fossRate = s->getBDSSFossilSampRatePsi();
    double sppSampRate = s->getBDSSSppSampRateRho();
    
    // maybe unnecessary here, for debugging -> check this does/doesn't make a difference
    //getFossilGraphProb(lambda, mu, fossRate, sppSampRate, originTime);
    
    vector<int> rndFossilRangeIDs;
    for(int i=0; i<fossilRanges.size(); i++)
        rndFossilRangeIDs.push_back(i);
    random_shuffle(rndFossilRangeIDs.begin(), rndFossilRangeIDs.end());
    
    for(vector<int>::iterator it=rndFossilRangeIDs.begin(); it!=rndFossilRangeIDs.end(); it++){
        
        FossilRange *fr = fossilRanges[(*it)];
        
        if(fr->getIsExtant())
            continue;
        
        if(fr->getIsFixStop())
            continue;
        
        // define old values
        double la = fr->getLastAppearance(); // af
        double oldEnd = fr->getLineageStop(); // bf
        double oldLike = currentFossilRangeGraphLnL; //rw: where is this initially defined?
        
        // propose new values
        double newEnd = 0.0;
        double c = 0.0;
        
        if(proposal == 1){ //rw: probably not appropriate for this parameter
            double rdwindow = 0.2;
            newEnd = getNewValSWindoMv(oldEnd, 0.0, la, rdwindow);
        }
        
        if(proposal == 2){ //rw: double check this move w/ TAH
            double rv = ranPtr->uniformRv();
            double tv = log(2.0); //rw: tuningVal
            c = doAScaleMove(newEnd, oldEnd, tv, 0.0, la, rv);
        }
        
        // redefine values
        fr->setLineageStop(newEnd);
        recountFossilRangeAttachNums();
        
        // recalculate the FRG probability
        double newLike = getFossilRangeGraphProb(lambda, mu, fossRate, sppSampRate, originTime);
        
        // calculate the likelihood/prior ratio
        double lnLikeRatio = newLike - oldLike;
        double r = modelPtr->safeExponentiation(lnLikeRatio + c);
        
        if(ranPtr->uniformRv() > r){
            fr->setLineageStop(oldEnd);
            recountFossilRangeAttachNums();
            currentFossilRangeGraphLnL = oldLike;
        }
    }
    
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

string FossilRangeGraph::getFossilRangeInfoParamNames(void){
    
    stringstream ss;
    FossilRange *frg = NULL;
    for(int i=0; i<fossilRanges.size(); i++){
        frg = fossilRanges[i];
        int frgID = frg->getFossilRangeID();
        ss << "\ty_f(FR_" << frgID << ")"; //rw: first appearance
        ss << "\tb_f(FR_" << frgID << ")"; //rw: lineage start
        ss << "\tx_f(FR_" << frgID << ")"; //rw: last appearance
        ss << "\td_f(FR_" << frgID << ")"; //rw: lineage stop
    }
    for(int i=0; i<fossilRanges.size(); i++){
        frg = fossilRanges[i];
        int frgID = frg->getFossilRangeID();
        ss << "\tgamma_f(FR_" << frgID << ")"; //".nd" << nID << ")";
    }
    string ni = ss.str();
    return ni;
}

string FossilRangeGraph::getFossilRangeInfoParamList(void){
    
    stringstream ss;
    FossilRange *frg = NULL;
    
    for(int i=0; i<fossilRanges.size(); i++){
        frg = fossilRanges[i];
        ss << "\t" << frg->getFirstAppearance(); //rw: first appearance -- note this is fixed in this implementation
        ss << "\t" << frg->getLineageStart(); //rw: lineage start
        ss << "\t" << frg->getLastAppearance(); //rw: last appearance
        ss << "\t" << frg->getLineageStop(); //rw: lineage stop
    }
    for(int i=0; i<fossilRanges.size(); i++){
        frg = fossilRanges[i];
        ss << "\t" << frg->getFossilRangeBrGamma();
    }
    string ni = ss.str();
    return ni;
}

//FBD process augmenting the start and end of species

double FossilRangeGraph::getFossilRangeGraphProb(double lambda, double mu, double fossRate, double sppSampRate, double ot){
    if(runUnderPrior)
        return 0.0;
    
    double nprb = 0.0;
    
    if(getAltProb)
        nprb = getFossilRangeGraphAlternativeProb(lambda, mu, fossRate, sppSampRate, ot);
    
    else {
        nprb = numFossils*log(fossRate);
        nprb += numExtinctLineages*log(mu);
        nprb -= log(lambda * (1-fbdPFxn(lambda,mu,fossRate, sppSampRate,ot)) );
        
        for(int f=0; f < fossilRanges.size(); f++){
            
            FossilRange *fr = fossilRanges[f];
            
            double bi = fr->getLineageStart(); //rw: zf or bi
            double di = fr->getLineageStop();  //rw: di
            
            nprb += log( lambda * fr->getFossilRangeBrGamma() );
            
            double rangePr = 0;
            rangePr += fbdQTildaFxnLog(lambda, mu, fossRate, sppSampRate, bi);
            rangePr -= fbdQTildaFxnLog(lambda, mu, fossRate, sppSampRate, di);
            
            nprb += rangePr;
        }
    }
    
    currentFossilRangeGraphLnL = nprb;
    
    return nprb;
}

double FossilRangeGraph::fbdC1Fxn(double b, double d, double psi){
    
    double v = fabs( sqrt( ( (b-d-psi) * (b-d-psi) ) + 4*b*psi) );
    return v;
}

double FossilRangeGraph::fbdC2Fxn(double b, double d, double psi, double rho){
    
    double v = -( ( b-d-(2*b*rho)-psi ) / (fbdC1Fxn(b,d,psi)) );
    return v;
}

double FossilRangeGraph::fbdC3Fxn(double b, double d, double psi, double rho){
    
    double v = b * (-psi + rho * (d + b * (-1 +rho) + psi ) );
    return v;
}

double FossilRangeGraph::fbdC4Fxn(double b, double d, double psi, double rho){
    
    double v = fbdC3Fxn(b,d,psi,rho) / (fbdC1Fxn(b,d,psi) * fbdC1Fxn(b,d,psi));
    return v;
}

double FossilRangeGraph::fbdPFxn(double b, double d, double psi, double rho, double t){
    
    double c1Val = fbdC1Fxn(b,d,psi);
    double c2Val = fbdC2Fxn(b,d,psi,rho);
    
    double expC1MinusC2 = exp(-c1Val * t) * (1.0 - c2Val);
    
    double eCfrac = (expC1MinusC2 - (1.0 + c2Val)) / (expC1MinusC2 + (1.0 + c2Val));
    double v = 1.0 + ((-(b - d - psi)) + (c1Val * eCfrac)) / (2.0 * b);
    
    return v;
}

//rw: this is unsuitable for large values of t
double FossilRangeGraph::fbdQTildaFxn(double b, double d, double psi, double rho, double t){
    
    double c1 = fbdC1Fxn(b,d,psi);
    double c2 = fbdC2Fxn(b,d,psi,rho);
    double c4 = fbdC4Fxn(b,d,psi,rho);
    
    //double v = rho * sqrt( -( exp(t * (b + d + psi + c1) ) / (c4 * ( (1-exp(-t * c1))*(1-exp(-t * c1)) )  -exp(-t * c1) ) * (  ( (1-c2) * (2 * c4 * (exp(-t*c1) -1) ) + exp(-t*c1) * (1-c2*c2) )  /  ( (1+c2) * (2 * c4 * (exp(-t*c1) - 1)) + exp(-t*c1) * (1-c2*c2) )  ) ) );
    double f1a = exp(-t * (b + d + psi + c1) );
    double f1b = (c4 * ( (1-exp(-t * c1))*(1-exp(-t * c1)) )  -exp(-t * c1) );
    double f2a = ( (1-c2) * (2 * c4 * (exp(-t*c1) -1) ) + exp(-t*c1) * (1-c2*c2) );
    double f2b = ( (1+c2) * (2 * c4 * (exp(-t*c1) - 1)) + exp(-t*c1) * (1-c2*c2) );
    double f  = (f1a/f1b) * (f2a/f2b) * -1;
    
    double v = rho * sqrt(f);

    return v;
}

double FossilRangeGraph::fbdQTildaFxnLog(double b, double d, double psi, double rho, double t){
    
    //b=2.5;
    //d=1.7;
    //psi=0.9;
    //rho=0.8;
    //t=1.3;
    
    double c1 = fbdC1Fxn(b,d,psi);
    double c2 = fbdC2Fxn(b,d,psi,rho);
    double c4 = fbdC4Fxn(b,d,psi,rho);
    
    double f1aLog = -t * (b + d + psi + c1) ;
    double f1b = (c4 * ( (1-exp(-t * c1))*(1-exp(-t * c1)) )  -exp(-t * c1) );
    double f2a = ( (1-c2) * (2 * c4 * (exp(-t*c1) -1) ) + exp(-t*c1) * (1-c2*c2) );
    double f2b = ( (1+c2) * (2 * c4 * (exp(-t*c1) - 1)) + exp(-t*c1) * (1-c2*c2) );
    
    double f  = f1aLog + log(-1/f1b * (f2a/f2b));
    
    double v = f*0.5 + log(rho);
    
    //cout << "t " << t << endl;
    //cout << "q_tilda " << setprecision(12) << exp(v) << endl;
    
    return v;
}

//rw: another test; sppSampRate not used
double FossilRangeGraph::getFossilRangeGraphAlternativeProb(double lambda, double mu, double fossRate, double sppSampRate, double ot){
    if(runUnderPrior)
        return 0.0;
    
    double nprb = 0.0;
    
    nprb = numFossils*log(fossRate);
    nprb += numExtinctLineages*log(mu);
    nprb -= log(lambda * (1-fbdPFxn(lambda,mu,fossRate,sppSampRate,ot)) );
    
    for(int f=0; f < fossilRanges.size(); f++){
        
        FossilRange *fr = fossilRanges[f];
        
        double bi = fr->getLineageStart(); //rw: zf
        double di = fr->getLineageStop();  //rw: bf
        
        int gamma = fr->getFossilRangeBrGamma();
        
        nprb += log( lambda * gamma );
        
        double rangePr;
        rangePr = -(lambda+mu+fossRate)*(bi-di);
        
        nprb += rangePr;
        
    }
    
    currentFossilRangeGraphLnL = nprb;//rw: some redundancy for now
    
    return nprb;
}

//END