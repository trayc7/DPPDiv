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

FossilRangeGraph::FossilRangeGraph(MbRandom *rp, Model *mp, int nf, int nl, vector<Calibration *> clb, bool rnp, bool fxFRG, bool estExt, int compS, int expMode) : Parameter(rp, mp){
    
    name = "FRG";
    numFossils = nf;
    numLineages = nl;
    runUnderPrior = rnp;
    conditionOnSurvival = 1;
    numExtinctLineages = 0;
    numExtantSamples = 0;
    originTime = 0.0;
    ancientBound = 1000.0;
    fixOrigin = 0;
    orderStartStopTimes = 0;
    fixFRG = fxFRG; //1: fix start and end range times to FAs and LAs
    estimateExtant = estExt;
    phyloTest = 0;
    if(expMode == 1){
        fixOrigin = 1;
        originTime = 1.233376;
        ancientBound = originTime;
        orderStartStopTimes = 1;
    }
    printInitialFossilRangeVariables = 1;
    createFossilRangeVector(clb);
    initializeFossilRangeVariables();
    currentFossilRangeGraphLnL = 0.0;
    moves = 1; // 1: update lineage start or stop times; 2: update both
    proposal = 2; // proposal type 1=window, 2=scale, 3=slide
    getAltProb = 0;
    completeSampling = compS; // likelihood function - if this is 0 it should produce the same results as 1 given complete sampling
    bool crossValidate = 0;
    if(crossValidate)
        crossValidateFBDfunctions();
    
    cout << "Number of lineages: " << numLineages << endl;
    cout << "Number of extinct ranges: " << numExtinctLineages << endl;
    cout << "Number of extant samples: " << numExtantSamples << endl;
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
    numExtinctLineages = frg.numExtinctLineages;
    numExtantSamples = frg.numExtantSamples;
    originTime = frg.originTime;
    
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
    
}

double FossilRangeGraph::update(double &oldLnL){
    
    currentFossilRangeGraphLnL = oldLnL;
    
    if(moves == 1){
        if(ranPtr->uniformRv() < 0.5)
            updateLineageStartTimes();
        else{
          updateExtinctIndicator();
          updateLineageStopTimes();
        }
    }
    else if (moves == 2) {
        updateLineageStartTimes();
        updateLineageStopTimes();
    }
    
    if(orderStartStopTimes)
        orderRangeAges();
    
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
    
    int frid=1;
    for(int c = 0; c < clb.size(); c++){
        Calibration *p = clb[c];
        double fa = p->getFirstAppearance();
        double la = p->getLastAppearance();
        double at = p->getAttachmentTime();
        double et = p->getEndTime();
        bool e = p->getIsExtant();
        bool eo = p->getIsExtantOnly();
        
        FossilRange *fr = new FossilRange(fa, la, at, et, e, eo, frid);
        fossilRanges.push_back(fr);
        
        frid ++;
    }
}

void FossilRangeGraph::initializeFossilRangeVariables(){
    
    double stop, start, la, fa;
    
    for(int f = 0; f < numLineages; f++){
        FossilRange *fr = fossilRanges[f];
        
        // lineage start times = speciation times
        if(fr->getIsExtant()){
            start = ranPtr->uniformRv(0.0,ancientBound);
            fr->setLineageStart(start);
            fr->setFixStart(0);
        }
        else{
            fa = fr->getFirstAppearance(); // yf
            start = ranPtr->uniformRv(fa,ancientBound); // zf
            fr->setLineageStart(start);
            fr->setFixStart(0);
        }
        
        // lineage stop times = extinction times
        if(fr->getIsExtant()){
            fr->setLineageStop(0.0);
            fr->setFixStop(1);
            fr->setExtinctIndicator(0);
        }
        else {
            la = fr->getLastAppearance();
            stop = ranPtr->uniformRv(0.0,la); //di
            fr->setLineageStop(stop);
            fr->setFixStop(0);
            fr->setExtinctIndicator(1);
        }
    }

    // redefines the variables defined in the above loop
    if(fixFRG){
        for(int f = 0; f < numLineages; f++){
            FossilRange *fr = fossilRanges[f];
            fr->setLineageStart(fr->getAttachmentTime());
            fr->setLineageStop(fr->getEndTime());
            fr->setFixStart(1);
            fr->setFixStop(1);
        }
    }
    
    // this is only used for expmo = 1
    if(fixOrigin){
        int of = 0;
        double ofa = 0;
        for(int f = 0; f < numLineages; f++){
            FossilRange *fr = fossilRanges[f];
            double fa = fr->getFirstAppearance();
            if(fa > ofa){
                ofa = fa;
                of = f;
            }
        }
        FossilRange *fr = fossilRanges[of];
        fr->setLineageStart(originTime);
        fr->setFixStart(1);
    }
    
    redefineOriginTime();
    recountFossilRangeAttachNums();
    countExtinctLineages();
    
    if(printInitialFossilRangeVariables) // debugging code
        printFossilRangeVariables();
    
}

// debugging code
void FossilRangeGraph::orderRangeAges(){
    
    double start, stop;
    std::vector<double> agesStart, agesStop;
    
    for(int f = 0; f < numLineages; f++){
        FossilRange *fr = fossilRanges[f];
        start = fr->getLineageStart();
        agesStart.push_back(start);
        stop = fr->getLineageStop();
        agesStop.push_back(stop);
        
    }
    
    std::sort(agesStart.begin(), agesStart.end(), std::greater<double>());
    std::sort(agesStop.begin(), agesStop.end(), std::greater<double>());
    
    for(int f = 0; f < numLineages; f++){
        FossilRange *fr = fossilRanges[f];
        fr->setLineageStart(agesStart[f]);
        fr->setLineageStop(agesStop[f]);
    }
    
    agesStart.clear();
    agesStop.clear();
}

void FossilRangeGraph::redefineOriginTime(){
    if(fixOrigin)
        return;
    
    double ot = 0.0;
    
    for(int f = 0; f < numLineages; f++){
        FossilRange *fr = fossilRanges[f];
        double ls = fr->getLineageStart();
        if(ls > ot)
            ot = ls;
    }
    
    originTime = ot;
}

void FossilRangeGraph::recountFossilRangeAttachNums(){
    
    for(int f = 0; f < numLineages; f++){
        FossilRange *fr = fossilRanges[f];
        double zf = fr->getLineageStart();
        int g = 0;
        for(int j = 0; j < numLineages; j++){
            if(f != j){
                FossilRange *p = fossilRanges[j];
                double ind = 0.0;
                if(p->getExtinctIndicator())
                    ind = 1.0;
                double zj = p->getLineageStart();
                double bj = p->getLineageStop() * ind; //TODO: double check this!
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
    
    int m = 0;
    int l = 0;
    
    for(int f = 0; f < numLineages; f++){
        FossilRange *fr = fossilRanges[f];
        if(!(fr->getLineageStop() == 0))
            m += 1;
        if(fr->getLastAppearance() == 0)
            l += 1;
    }
    
    numExtinctLineages = m;
    numExtantSamples = l;
    
}

// debugging code
void FossilRangeGraph::printFossilRangeVariables(){
    
    for(int f = 0; f < numLineages; f++){
        FossilRange *fr = fossilRanges[f];
        cout << "Fossil range ID: " << fr->getFossilRangeID() << endl;
        cout << "First appearance: " << fr->getFirstAppearance() << endl;
        cout << "Last appearance: " << fr->getLastAppearance() << endl;
        cout << "Lineage start: " << fr->getLineageStart() << endl;
        cout << "Lineage end: " << fr->getLineageStop() << endl;
        cout << "Is extant: " << fr->getIsExtant() << endl;
        cout << "Gamma: " << fr->getFossilRangeBrGamma() << endl;
    }
    
    cout << "Origin time: " << originTime << endl;
    
}

void FossilRangeGraph::printFossilRangeVariables(int range){
    
    int f = range;
    
    FossilRange *fr = fossilRanges[f];
    cout << "Fossil range ID: " << fr->getFossilRangeID() << endl;
    cout << "First appearance: " << fr->getFirstAppearance() << endl;
    cout << "Last appearance: " << fr->getLastAppearance() << endl;
    cout << "Lineage start: " << fr->getLineageStart() << endl;
    cout << "Lineage end: " << fr->getLineageStop() << endl;
    cout << "Is extant: " << fr->getIsExtant() << endl;
    cout << "Gamma: " << fr->getFossilRangeBrGamma() << endl;
    cout << "Is fix start: " << fr->getIsFixStart() << endl;
    cout << "Is fix stop: " << fr->getIsFixStop() << endl;
    
    cout << "Origin time: " << originTime << endl;
    
}

double FossilRangeGraph::updateLineageStartTimes(){
    
    Speciation *s = modelPtr->getActiveSpeciation();
    s->setAllBDFossParams();
    double lambda = s->getBDSSSpeciationRateLambda();
    double mu = s->getBDSSExtinctionRateMu();
    double fossRate = s->getBDSSFossilSampRatePsi();
    double sppSampRate = s->getBDSSSppSampRateRho();
    
    // debugging
    // getFossilRangeGraphProb(lambda, mu, fossRate, sppSampRate, originTime);
    
    vector<int> rndFossilRangeIDs;
    for(int i=0; i<fossilRanges.size(); i++)
        rndFossilRangeIDs.push_back(i);
    random_shuffle(rndFossilRangeIDs.begin(), rndFossilRangeIDs.end());
    
    for(vector<int>::iterator it=rndFossilRangeIDs.begin(); it!=rndFossilRangeIDs.end(); it++){
        
        FossilRange *fr = fossilRanges[(*it)];
        
        if(fr->getIsFixStart())
            continue;
        
        if(fixOrigin && fr->getLineageStart() == originTime)
            continue;
        
        // define old values
        double fa = fr->getFirstAppearance(); // yf
        double oldStart = fr->getLineageStart(); // zf
        double oldLike = currentFossilRangeGraphLnL;

        // propose new values
        double newStart = 0.0;
        double c = 0.0;
        
        if(proposal == 1){ // probably not appropriate for this parameter
            double rdwindow = 0.2;
            newStart = getNewValSWindoMv(oldStart, fa, ancientBound, rdwindow);
        }
        
        if(proposal == 2){
            double rv = ranPtr->uniformRv();
            double tv = 2.0; // tuningVal
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
    
    // debugging code
    //printFossilRangeVariables(0);
    
    return 0.0;
}

double FossilRangeGraph::updateLineageStopTimes(){
    
    Speciation *s = modelPtr->getActiveSpeciation();
    s->setAllBDFossParams();
    double lambda = s->getBDSSSpeciationRateLambda();
    double mu = s->getBDSSExtinctionRateMu();
    double fossRate = s->getBDSSFossilSampRatePsi();
    double sppSampRate = s->getBDSSSppSampRateRho();
    
    // this shouldn't be necessary here now
    //getFossilGraphProb(lambda, mu, fossRate, sppSampRate, originTime);
    
    vector<int> rndFossilRangeIDs;
    for(int i=0; i<fossilRanges.size(); i++)
        rndFossilRangeIDs.push_back(i);
    random_shuffle(rndFossilRangeIDs.begin(), rndFossilRangeIDs.end());
    
    for(vector<int>::iterator it=rndFossilRangeIDs.begin(); it!=rndFossilRangeIDs.end(); it++){
        
        FossilRange *fr = fossilRanges[(*it)];
        
        if(fr->getIsFixStop())
            continue;
        
        if(!fr->getExtinctIndicator())
            continue;
        
        // define old values
        double la = fr->getLastAppearance(); // af
        double oldEnd = fr->getLineageStop(); // bf
        double oldLike = currentFossilRangeGraphLnL;
        
        // propose new values
        double newEnd = 0.0;
        double c = 0.0;
        
        if(proposal == 1){ // probably not appropriate for this parameter
            double rdwindow = 0.2;
            newEnd = getNewValSWindoMv(oldEnd, 0.0, la, rdwindow);
        }
        
        if(proposal == 2){
            double rv = ranPtr->uniformRv();
            double tv = log(2.0); //tuningVal
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

double FossilRangeGraph::updateExtinctIndicator(){
    
    Speciation *s = modelPtr->getActiveSpeciation();
    s->setAllBDFossParams();
    double lambda = s->getBDSSSpeciationRateLambda();
    double mu = s->getBDSSExtinctionRateMu();
    double fossRate = s->getBDSSFossilSampRatePsi();
    double sppSampRate = s->getBDSSSppSampRateRho();
    
    vector<int> rndFossilRangeIDs;
    for(int i=0; i<fossilRanges.size(); i++)
        rndFossilRangeIDs.push_back(i);
    random_shuffle(rndFossilRangeIDs.begin(), rndFossilRangeIDs.end());
    
    for(vector<int>::iterator it=rndFossilRangeIDs.begin(); it!=rndFossilRangeIDs.end(); it++){
        
        FossilRange *fr = fossilRanges[(*it)];
        
        if(fr->getIsFixStop())
            continue;
        
        // define old values
        bool oldInd = fr->getExtinctIndicator();
        double oldLike = currentFossilRangeGraphLnL;
        
        // propose new indicator value
        bool newInd;
        double rvInd = ranPtr->discreteUniformRv(0,1);
        if(rvInd == 0)
            newInd = 0;
        else
            newInd = 1;
        
        if(newInd == oldInd) continue;
        
        // redefine values
        fr->setExtinctIndicator(newInd);
        recountFossilRangeAttachNums();
            
        // recalculate the FRG probability
        double newLike = getFossilRangeGraphProb(lambda, mu, fossRate, sppSampRate, originTime);
            
        // calculate the likelihood/prior ratio
        double lnLikeRatio = newLike - oldLike;
        double r = modelPtr->safeExponentiation(lnLikeRatio);
            
        if(ranPtr->uniformRv() > r){
           fr->setExtinctIndicator(oldInd);
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
        ss << "\ty_f(FR_" << frgID << ")"; // first appearance
        ss << "\tb_f(FR_" << frgID << ")"; // lineage start
        ss << "\tx_f(FR_" << frgID << ")"; // last appearance
        ss << "\td_f(FR_" << frgID << ")"; // lineage stop
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
        ss << "\t" << frg->getFirstAppearance(); // first appearance -- note this is fixed in this implementation
        ss << "\t" << frg->getLineageStart(); // lineage start
        ss << "\t" << frg->getLastAppearance(); // last appearance
        ss << "\t" << frg->getLineageStop(); // lineage stop
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
    if(runUnderPrior){
        countExtinctLineages();//TODO: something more elegant with this
        return 0.0;
    }
    
    double nprb = 0.0;

    // debugging
    //lambda = 1.0; mu = 0.1; fossRate = 1.368238; sppSampRate = 1.0;
    
    if(getAltProb)
        nprb = getFossilRangeGraphAlternativeProb(lambda, mu, fossRate, sppSampRate, ot);
    
    else if(phyloTest)
        nprb = getPhyloProb(lambda, mu, sppSampRate, ot);
    
    else if(completeSampling == 1) {
        
        nprb = numFossils*log(fossRate);
        nprb += numExtinctLineages*log(mu);
        nprb -= log(lambda * (1-fbdPFxn(lambda,mu,fossRate,sppSampRate,ot)) );
        
        for(int f=0; f < fossilRanges.size(); f++){
            
            FossilRange *fr = fossilRanges[f];
            
            double bi = fr->getLineageStart(); // bi (zf)
            double di = fr->getLineageStop();  // di
            
            nprb += log( lambda * fr->getFossilRangeBrGamma() );
            
            double rangePr = 0;
            rangePr += fbdQTildaFxnLog(lambda, mu, fossRate, sppSampRate, bi);
            rangePr -= fbdQTildaFxnLog(lambda, mu, fossRate, sppSampRate, di);
            
            nprb += rangePr;
        }
    }
    // Keiding 1975 + Poisson fossil samping
    else if(completeSampling == 2) {
        
        int birthEvents = -1; //B. Note the origin is not a birth event
        int deathEvents = 0; //D
        double totalLineageDuration = 0; //S
        
        for(int f=0; f < fossilRanges.size(); f++){
            
            FossilRange *fr = fossilRanges[f];
            
            double bi = fr->getLineageStart(); // bi
            double di = fr->getLineageStop();  // di
            
            totalLineageDuration += bi - di;
            birthEvents += 1;
            if(di > 0)
                deathEvents +=1;
        }
        
        nprb = birthEvents * log(lambda);
        nprb += deathEvents * log(mu);
        nprb += numFossils * log(fossRate);
        nprb += -(lambda + mu + fossRate) * totalLineageDuration;
        
    }
    
    // correctly accounting for incomplete sampling - Stadler et al. 2018, eq. 7
    else {
        
        nprb = numFossils*log(fossRate);
        
        //nprb += numExtinctLineages*log(mu);
        
        if(sppSampRate == 0) conditionOnSurvival = 0;
        
        if(conditionOnSurvival)
            nprb -= log(lambda * (1-fbdPFxn(lambda,mu,fossRate,sppSampRate,ot)) );
        
        numExtinctLineages = 0;
        
        for(int f=0; f < fossilRanges.size(); f++){
            
            FossilRange *fr = fossilRanges[f];
            
            double ind = 0.0;
            if(fr->getExtinctIndicator())
                ind = 1.0;
            
            double bi = fr->getLineageStart(); // bi
            double di = fr->getLineageStop() * ind;  // di
            double oi = fr->getFirstAppearance(); // oi
            
            // extinction
            if(di != 0){
                nprb += log(mu);
                numExtinctLineages += 1;
            }
            
            // speciation
            nprb += log( lambda * fr->getFossilRangeBrGamma() );
            
            double rangePr = 0;
            
            rangePr += fbdQTildaFxnLog(lambda, mu, fossRate, sppSampRate, oi);
            rangePr -= fbdQTildaFxnLog(lambda, mu, fossRate, sppSampRate, di);
            
            rangePr += fbdQFxnLog(lambda, mu, fossRate, sppSampRate, bi);
            rangePr -= fbdQFxnLog(lambda, mu, fossRate, sppSampRate, oi);
            
            nprb += rangePr;
        }
        
        //cout << "numExtinctLineages = " << numExtinctLineages << endl;
        
        // extant species sampling
        if(sppSampRate < 1 & sppSampRate > 0)
            nprb += ( numExtantSamples * log(sppSampRate) ) + ( (numLineages - numExtinctLineages - numExtantSamples) * log(1 - sppSampRate) );
        
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

// this is unsuitable for large values of t
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

double FossilRangeGraph::fbdQTildaFxnLogAlt(double b, double d, double psi, double rho, double t){
    
    double c1 = fbdC1Fxn(b,d,psi);
    double c2 = fbdC2Fxn(b,d,psi,rho);
    double c4 = fbdC4Fxn(b,d,psi,rho);
    
    double f1aLog = -t * (b + d + psi + c1) ;
    double f1b = (c4 * ( (1-exp(-t * c1))*(1-exp(-t * c1)) )  -exp(-t * c1) );
    double f2a = ( (1-c2) * (2 * c4 * (exp(-t*c1) -1) ) + exp(-t*c1) * (1-c2*c2) );
    double f2b = ( (1+c2) * (2 * c4 * (exp(-t*c1) - 1)) + exp(-t*c1) * (1-c2*c2) );
    
    double f  = f1aLog + log(-1/f1b * (f2a/f2b));
    
    double v = f*0.5 + log(rho);
    
    return v;
}

// simplified version
double FossilRangeGraph::fbdQTildaFxnLog(double b, double d, double psi, double rho, double t){
    
    double q = fbdQFxnLog(b, d, psi, rho, t);
    double v = (q + (- (b + d + psi)*t) ) * 0.5;
    
    return v;
}

double FossilRangeGraph::fbdQFxnLog(double b, double d, double psi, double rho, double t){
    
    double c1 = fbdC1Fxn(b,d,psi);
    double c2 = fbdC2Fxn(b,d,psi,rho);
    
    double f1 = log(4) + (-c1 * t);
    double f2 = 2 * (log( (exp(-c1*t) * (1-c2)) + (1+c2) ));
    
    double v = f1 - f2;
    
    return v;
}

// another test; sppSampRate not used
double FossilRangeGraph::getFossilRangeGraphAlternativeProb(double lambda, double mu, double fossRate, double sppSampRate, double ot){
    if(runUnderPrior)
        return 0.0;
    
    double nprb = 0.0;
    
    nprb = numFossils*log(fossRate);
    nprb += numExtinctLineages*log(mu);
    nprb -= log(lambda * (1-fbdPFxn(lambda,mu,fossRate,sppSampRate,ot)) );
    
    for(int f=0; f < fossilRanges.size(); f++){
        
        FossilRange *fr = fossilRanges[f];
        
        double bi = fr->getLineageStart(); // zf
        double di = fr->getLineageStop();  // bf
        
        int gamma = fr->getFossilRangeBrGamma();
        
        nprb += log( lambda * gamma );
        
        double rangePr;
        rangePr = -(lambda+mu+fossRate)*(bi-di);
        
        nprb += rangePr;
        
    }
    
    currentFossilRangeGraphLnL = nprb;
    
    return nprb;
}

// Stadler 2010, equation 2
double FossilRangeGraph::getPhyloProb(double lambda, double mu, double sppSampRate, double ot){
    
    double nprb = 0.0;
    
    // deal with the origin
    nprb = phyloBDP1FxnLog(lambda, mu, sppSampRate, ot);
    nprb -= log(1 - phyloBDP0Fxn(lambda, mu, sppSampRate, ot));
    
    cout << nprb << endl;
    
    // deal with each internal node
    for(int f=0; f < fossilRanges.size(); f++){
        
        FossilRange *fr = fossilRanges[f];
        
        double bi = fr->getLineageStart();
        
        if(bi==ot)
            continue;
        
        nprb += log(lambda) + phyloBDP1FxnLog(lambda, mu, sppSampRate, bi);
    }
    
    return nprb;
    
}

double FossilRangeGraph::phyloBDP0Fxn(double b, double d, double rho, double t){
    
    double p = 1 - ( ( rho * (b-d) ) / ( (rho * b) + (b * (1-rho) - d) * exp(-(b-d)*t) ) );
    
    return p;
}

double FossilRangeGraph::phyloBDP1FxnLog(double b, double d, double rho, double t){
    
    double p = 2*log(rho * (b-d)) + (-(b-d)*t) - 2*log( (rho * b) + (b*(1-rho)-d) * exp(-(b-d)*t) );
    
    return p;
}

void FossilRangeGraph::crossValidateFBDfunctions(){
    
    double lambda = 1.0;
    double mu = 0.1;
    double fossRate = 10.0;
    double sppSampRate = 1.0;
    double ot = originTime;
    
    if(phyloTest){
        
        double p0 = phyloBDP0Fxn(lambda, mu, sppSampRate, ot);
        double p1 = phyloBDP1FxnLog(lambda, mu, sppSampRate, ot);
        
        cout << "Cross validating phylo BD functions..." << endl;
        cout << "p0(ot) = " << p0 << endl;
        cout << "p1(ot) log = " << p1 << endl;
        
        double bdProb = getPhyloProb(lambda, mu, sppSampRate, ot);
        
        cout << "BD probability = " << setprecision(7) << bdProb << endl;
        
    }
    else {
        double c1 = fbdC1Fxn(lambda, mu, fossRate);
        double c2 = fbdC2Fxn(lambda, mu, fossRate, sppSampRate);
        double c3 = fbdC3Fxn(lambda, mu, fossRate, sppSampRate);
        double c4 = fbdC4Fxn(lambda, mu, fossRate, sppSampRate);
        
        cout << "Cross validating FBD functions..." << endl;
        cout << "c1 = " << c1 << endl;
        cout << "c2 = " << c2 << endl;
        cout << "c3 = " << c3 << endl;
        cout << "c4 = " << c4 << endl;
        
        double fbdP = fbdPFxn(lambda, mu, fossRate, sppSampRate, ot);
        double fbdQLog = fbdQFxnLog(lambda, mu, fossRate, sppSampRate, ot);
        double fbdQtildaLog = fbdQTildaFxnLog(lambda, mu, fossRate, sppSampRate, ot);
        
        cout << "P(ot) = " << fbdP << endl;
        cout << "Q(ot) log = " << fbdQLog << endl;
        cout << "Q tilda (ot) log = " << fbdQtildaLog << endl;
        
        double fbdProb = getFossilRangeGraphProb(lambda, mu, fossRate, sppSampRate, ot);

        cout << "FBD probability = " << setprecision(7) << fbdProb << endl;

    }
    
    exit(0);
    
}

// generate likelihood surface plot
void FossilRangeGraph::lnSurfaceGenerator(string outFile){
    
    // output file for likelihood surface data
    string lnSurfFile = outFile + ".lnSurface.out";
    ofstream lnSurfOut(lnSurfFile.c_str(), ios::out);
    
    lnSurfOut << std::setprecision(9) << "l\tm\tlogL\n";
    
    double birthMin = 0.01;
    double birthMax = 10.0;
    double deathMin = 0.01;
    double deathMax = 10.0;
    
    Speciation *s = modelPtr->getActiveSpeciation();
    double fossRate = s->getBDSSFossilSampRatePsi();
    double sppSampRate = s->getBDSSSppSampRateRho();
    
    for(double i = birthMin; i <= birthMax; i += 0.1) {
        
        for(double j = deathMin; j <= deathMax; j += 0.1){
            
            // calculate FBD probability
            double lnL = getFossilRangeGraphProb(i, j, fossRate, sppSampRate, originTime);
            //cout << lnL << endl;
            lnSurfOut << i << "\t" << j << "\t" << lnL << "\n";
        }
    }
    
    lnSurfOut.close();
    exit(0);
    
}

//END
