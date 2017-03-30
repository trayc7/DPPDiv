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
#include "Parameter_fossilrangegraphskyline.h"
#include "Parameter_speciaton.h" //**skyline - you need to change this
#include "MbRandom.h"
#include "Model.h"
#include "util.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <algorithm>

using namespace std;

FossilRangeGraphSkyline::FossilRangeGraphSkyline(MbRandom *rp, Model *mp, int nf, int nl, vector<Calibration *> clb, int ni, vector<Calibration *> ints, bool rnp, bool fxFRG) : Parameter(rp, mp){

    name = "FRGS";
    numFossils = nf;
    numLineages = nl;
    numIntervals = ni + 1; // user defined intervals + the interval incorporating the ancient bound
    runUnderPrior = rnp;
    numExtinctLineages = 0;
    originTime = 0.0;
    originInterval = 0;
    ancientBound = 1000.0;
    fixOrigin = 0;
    fixFRG = fxFRG; //1: fix start and end range times to FAs and LAs
    printInitialFossilRangeSkylineVariables = 1;
    createIntervalsVector(ints);
    createFossilRangeSkylineVector(clb);
    initializeFossilRangeSkylineVariables();
    currentFossilRangeGraphSkylineLnL = 0.0;
    
    //**skynote add cross validate function
    
    cout << "Number of lineages: " << numLineages << endl;
    cout << "Number of extinct ranges: " << numExtinctLineages << endl;
    cout << "Number of fossils: " << numFossils << endl;
    cout << "Number of intervals (inc ancient bound): " << numIntervals << endl;
    cout << "\nInitial origin time: " << originTime << endl;
    cout << "Fossil range graph skyline initialized" << endl;
    
}

FossilRangeGraphSkyline::~FossilRangeGraphSkyline(void){
    
}

FossilRangeGraphSkyline& FossilRangeGraphSkyline::operator=(const FossilRangeGraphSkyline &frgsl) {
    
    if (this != &frgsl)
        clone(frgsl);
    return *this;
    
}

void FossilRangeGraphSkyline::clone(const FossilRangeGraphSkyline &frgsl){

    numFossils = frgsl.numFossils;
    numLineages = frgsl.numLineages;
    numExtinctLineages = frgsl.numExtinctLineages;
    originTime = frgsl.originTime;
    numIntervals = frgsl.numIntervals;
    
    for(int i=0; i<fossilRangesSkyline.size(); i++){
        FossilRangeSkyline *frTo = fossilRangesSkyline[i];
        FossilRangeSkyline *frFrom = frgsl.fossilRangesSkyline[i];
        frTo->setFossilRangeIndex(frFrom->getFossilRangeIndex());
        frTo->setFirstAppearance(frFrom->getFirstAppearance());
        frTo->setLastAppearance(frFrom->getLastAppearance());
        frTo->setExtant(frFrom->getIsExtant());
        frTo->setExtantOnly(frFrom->getIsExtantOnly());
        frTo->setFossilRangeID(frFrom->getFossilRangeID());
        frTo->setFossilRangeBrGamma(frFrom->getFossilRangeBrGamma());
        frTo->setLineageStart(frFrom->getLineageStart());
        frTo->setLineageStop(frFrom->getLineageStop());
        frTo->setFossilRangeFirstAppearanceInterval(frFrom->getFossilRangeFirstAppearanceInterval());
        frTo->setFossilRangeBirthInterval(frFrom->getFossilRangeBirthInterval());
        frTo->setFossilRangeDeathInterval(frFrom->getFossilRangeDeathInterval());
    }
    
}

double FossilRangeGraphSkyline::update(double &oldLnL){
    
    currentFossilRangeGraphSkylineLnL = oldLnL;
    
    if(ranPtr->uniformRv() < 0.5)
        updateLineageStartTimes();
    else
        updateLineageStopTimes();
    
    return currentFossilRangeGraphSkylineLnL;
}

double FossilRangeGraphSkyline::lnPrior(){
    
    return 0;
}

void FossilRangeGraphSkyline::print(ostream & o) const {
    
}

string FossilRangeGraphSkyline::writeParam(void){
    
    string s="";
    return s;
}

void FossilRangeGraphSkyline::createIntervalsVector(vector<Calibration *> ints){
    
    bool printIntVariables = 1;
    double start = 0;
    double end = 0;
    int fossils = 0;
    
    int intid = 1;
    for(int i = 0; i < ints.size(); i++){
        Calibration *h = ints[i];
        start = h->getIntervalStart();
        end = h->getIntervalEnd();
        fossils = h->getIntervalFossils();
        
        Interval *interval = new Interval(start, end, fossils, intid);
        intervals.push_back(interval);
        
        intid ++;
    }
    
    Interval *interval = new Interval(ancientBound, start, 0, intid);
    intervals.push_back(interval);
    
    if(printIntVariables)
        printIntervalVariables();
    
}

void FossilRangeGraphSkyline::printIntervalVariables(){
    
    for(int h = 0; h < numIntervals; h++){
        Interval *interval = intervals[h];
        cout << "Interval ID: " << interval->getIntervalID() << endl;
        cout << "Interval start: " << interval->getIntervalStart() << endl;
        cout << "Interval end: " << interval->getIntervalEnd() << endl;
        cout << "Number of fossils: " << interval->getIntervalFossils() << endl << endl;
    }
}

void FossilRangeGraphSkyline::createFossilRangeSkylineVector(vector<Calibration *> clb){
    
    int frid=1;
    for(int c = 0; c < clb.size(); c++){
        Calibration *p = clb[c];
        double fa = p->getFirstAppearance();
        double la = p->getLastAppearance();
        double at = p->getAttachmentTime();
        bool e = p->getIsExtant();
        bool eo = p->getIsExtantOnly();
        
        FossilRangeSkyline *frsl = new FossilRangeSkyline(fa, la, at, e, eo, frid);
        fossilRangesSkyline.push_back(frsl);
        
        frid ++;
    }
}

void FossilRangeGraphSkyline::initializeFossilRangeSkylineVariables(){
    
    double stop, start, la, fa;
    
    for(int f = 0; f < numLineages; f++){
        FossilRangeSkyline *fr = fossilRangesSkyline[f];
        
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
        }
        else {
            la = fr->getLastAppearance();
            stop = ranPtr->uniformRv(0.0,la); //di
            fr->setLineageStop(stop);
            fr->setFixStop(0);
        }
    }
    
    // redefines the variables defined in the above loop
    if(fixFRG){
        for(int f = 0; f < numLineages; f++){
            FossilRangeSkyline *fr = fossilRangesSkyline[f];
            fr->setLineageStart(fr->getAttachmentTime());
            fr->setLineageStop(fr->getLastAppearance());
            fr->setFixStart(1);
            fr->setFixStop(1);
        }
    }
    
    // this is only used for expmo = 1
    if(fixOrigin){
        int of = 0;
        double ofa = 0;
        for(int f = 0; f < numLineages; f++){
            FossilRangeSkyline *fr = fossilRangesSkyline[f];
            double fa = fr->getFirstAppearance();
            if(fa > ofa){
                ofa = fa;
                of = f;
            }
        }
        FossilRangeSkyline *fr = fossilRangesSkyline[of];
        fr->setLineageStart(originTime);
        fr->setFixStart(1);
    }
    
    redefineOriginTime();
    recountFossilRangeAttachNums();
    countExtinctLineages();
    
    // assign fas, birth and death time to intervals
    for(int f = 0; f < numLineages; f++){
        FossilRangeSkyline *fr = fossilRangesSkyline[f];
        fr->setFossilRangeFirstAppearanceInterval(assignInterval(fr->getFirstAppearance()));
        fr->setFossilRangeBirthInterval(assignInterval(fr->getLineageStart()));
        fr->setFossilRangeDeathInterval(assignInterval(fr->getLineageStop()));
    }
    
    if(printInitialFossilRangeSkylineVariables)
        printFossilRangeSkylineVariables();
    
}

void FossilRangeGraphSkyline::redefineOriginTime(){
    if(fixOrigin)
        return;
    
    double ot = 0.0;
    
    for(int f = 0; f < numLineages; f++){
        FossilRangeSkyline *fr = fossilRangesSkyline[f];
        double ls = fr->getLineageStart();
        if(ls > ot)
            ot = ls;
    }
    
    originTime = ot;
    originInterval = assignInterval(ot);
    
}


void FossilRangeGraphSkyline::recountFossilRangeAttachNums(){
    
    for(int f = 0; f < numLineages; f++){
        FossilRangeSkyline *fr = fossilRangesSkyline[f];
        double zf = fr->getLineageStart();
        int g = 0;
        for(int j = 0; j < numLineages; j++){
            if(f != j){
                FossilRangeSkyline *p = fossilRangesSkyline[j];
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

void FossilRangeGraphSkyline::countExtinctLineages(){
    
    int el = 0;
    
    for(int f = 0; f < numLineages; f++){
        FossilRangeSkyline *fr = fossilRangesSkyline[f];
        if(!fr->getIsExtant())
            el += 1;
    }
    
    numExtinctLineages = el;
    
}

int FossilRangeGraphSkyline::assignInterval(double time){
    
    int assingment = 0;
    for(int i = 0; i < numIntervals; i++){
        Interval *interval = intervals[i];
        double start = interval->getIntervalStart();
        if(time < start){
            assingment = i;
            break;
        }
    }
    
    return assingment;
}

void FossilRangeGraphSkyline::printFossilRangeSkylineVariables(){
    
    for(int f = 0; f < numLineages; f++){
        FossilRangeSkyline *fr = fossilRangesSkyline[f];
        cout << "Fossil range ID: " << fr->getFossilRangeID() << endl;
        cout << "First appearance: " << fr->getFirstAppearance() << " sampled in interval " << fr->getFossilRangeFirstAppearanceInterval() << endl;
        cout << "Last appearance: " << fr->getLastAppearance() << endl;
        cout << "Lineage start: " << fr->getLineageStart() << " sampled in interval " << fr->getFossilRangeBirthInterval() << endl;
        cout << "Lineage end: " << fr->getLineageStop() << " sampled in interval " << fr->getFossilRangeDeathInterval() << endl;
        cout << "Is extant: " << fr->getIsExtant() << endl;
        cout << "Gamma: " << fr->getFossilRangeBrGamma() << endl;
    }
    
    cout << "Origin time: " << originTime << " sampled in interval " << originInterval << endl;
    
}

void FossilRangeGraphSkyline::printFossilRangeVariables(int range){
    
    int f = range;
    
    FossilRangeSkyline *fr = fossilRangesSkyline[f];
    cout << "Fossil range ID: " << fr->getFossilRangeID() << endl;
    cout << "First appearance: " << fr->getFirstAppearance() << " sampled in interval " << fr->getFossilRangeFirstAppearanceInterval() << endl;
    cout << "Last appearance: " << fr->getLastAppearance() << endl;
    cout << "Lineage start: " << fr->getLineageStart() << " sampled in interval " << fr->getFossilRangeBirthInterval() << endl;
    cout << "Lineage end: " << fr->getLineageStop() << " sampled in interval " << fr->getFossilRangeDeathInterval() << endl;
    cout << "Is extant: " << fr->getIsExtant() << endl;
    cout << "Gamma: " << fr->getFossilRangeBrGamma() << endl;
    cout << "Is fix start: " << fr->getIsFixStart() << endl;
    cout << "Is fix stop: " << fr->getIsFixStop() << endl;
    
    cout << "Origin time: " << originTime << endl;
    
}

double FossilRangeGraphSkyline::updateLineageStartTimes(){
    
    //**skyline note - these should come from skyline speciation
    //Speciation *s = modelPtr->getActiveSpeciation();
    //s->setAllBDFossParams();
    //double lambda = s->getBDSSSpeciationRateLambda();
    //double mu = s->getBDSSExtinctionRateMu();
    //double fossRate = s->getBDSSFossilSampRatePsi();
    //double sppSampRate = s->getBDSSSppSampRateRho();
    
    double sppSampRate = 1;
    vector<double> lambda, mu, fossRate;
    for(int i = 0; numIntervals; i++){
        lambda.push_back(1);
        mu.push_back(0.1);
        fossRate.push_back(2);
    }
    
    vector<int> rndFossilRangeIDs;
    for(int i=0; i<fossilRangesSkyline.size(); i++)
        rndFossilRangeIDs.push_back(i);
    random_shuffle(rndFossilRangeIDs.begin(), rndFossilRangeIDs.end());
    
    for(vector<int>::iterator it=rndFossilRangeIDs.begin(); it!=rndFossilRangeIDs.end(); it++){
        
        FossilRangeSkyline *fr = fossilRangesSkyline[(*it)];
        
        if(fr->getIsFixStart())
            continue;
        
        if(fixOrigin && fr->getLineageStart() == originTime)
            continue;
        
        // define old values
        double fa = fr->getFirstAppearance(); // yf
        double oldStart = fr->getLineageStart(); // zf
        double oldLike = currentFossilRangeGraphSkylineLnL;
        
        // propose new values
        double newStart = 0.0;
        double c = 0.0;
        
    
        double rv = ranPtr->uniformRv();
        double tv = 2.0; //rw: tuningVal
        c = doAScaleMove(newStart, oldStart, tv, fa, ancientBound, rv);
        
        // redefine values
        fr->setLineageStart(newStart);
        redefineOriginTime();
        recountFossilRangeAttachNums();
        fr->setFossilRangeBirthInterval(assignInterval(newStart));
        
        // recalculate the FRG probability
        double newLike = getFossilRangeGraphSkylineProb(lambda, mu, fossRate, sppSampRate, originTime);
        
        // calculate the likelihood/prior ratio
        double lnLikeRatio = newLike - oldLike;
        double r = modelPtr->safeExponentiation(lnLikeRatio + c);
        
        if(ranPtr->uniformRv() > r){
            fr->setLineageStart(oldStart);
            redefineOriginTime();
            recountFossilRangeAttachNums();
            fr->setFossilRangeBirthInterval(assignInterval(oldStart));
            currentFossilRangeGraphSkylineLnL = oldLike;
        }
    }
    
    // debugging code
    //printFossilRangeVariables(0);
    
    return 0.0;
}

double FossilRangeGraphSkyline::updateLineageStopTimes(){
    
//    Speciation *s = modelPtr->getActiveSpeciation();
//    s->setAllBDFossParams();
//    double lambda = s->getBDSSSpeciationRateLambda();
//    double mu = s->getBDSSExtinctionRateMu();
//    double fossRate = s->getBDSSFossilSampRatePsi();
//    double sppSampRate = s->getBDSSSppSampRateRho();
    
    double sppSampRate = 1;
    vector<double> lambda, mu, fossRate;
    for(int i = 0; numIntervals; i++){
        lambda.push_back(1);
        mu.push_back(0.1);
        fossRate.push_back(2);
    }
    
    vector<int> rndFossilRangeIDs;
    for(int i=0; i<fossilRangesSkyline.size(); i++)
        rndFossilRangeIDs.push_back(i);
    random_shuffle(rndFossilRangeIDs.begin(), rndFossilRangeIDs.end());
    
    for(vector<int>::iterator it=rndFossilRangeIDs.begin(); it!=rndFossilRangeIDs.end(); it++){
        
        FossilRangeSkyline *fr = fossilRangesSkyline[(*it)];
        
        if(fr->getIsExtant())
            continue;
        
        if(fr->getIsFixStop())
            continue;
        
        // define old values
        double la = fr->getLastAppearance(); // af
        double oldEnd = fr->getLineageStop(); // bf
        double oldLike = currentFossilRangeGraphSkylineLnL; //rw: where is this initially defined?
        
        // propose new values
        double newEnd = 0.0;
        double c = 0.0;
        
        double rv = ranPtr->uniformRv();
        double tv = log(2.0); //rw: tuningVal
        c = doAScaleMove(newEnd, oldEnd, tv, 0.0, la, rv);
        
        // redefine values
        fr->setLineageStop(newEnd);
        recountFossilRangeAttachNums();
        fr->setFossilRangeDeathInterval(assignInterval(newEnd));
        
        // recalculate the FRG probability
        double newLike = getFossilRangeGraphSkylineProb(lambda, mu, fossRate, sppSampRate, originTime);
        
        // calculate the likelihood/prior ratio
        double lnLikeRatio = newLike - oldLike;
        double r = modelPtr->safeExponentiation(lnLikeRatio + c);
        
        if(ranPtr->uniformRv() > r){
            fr->setLineageStop(oldEnd);
            recountFossilRangeAttachNums();
            fr->setFossilRangeDeathInterval(assignInterval(oldEnd));
            currentFossilRangeGraphSkylineLnL = oldLike;
        }
    }
    
    return 0.0;
}

double FossilRangeGraphSkyline::doAScaleMove(double &nv, double cv, double tv, double lb, double hb, double rv){
    
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

double FossilRangeGraphSkyline::getActiveFossilRangeGraphSkylineProb(){
    if(runUnderPrior)
        return 0.0;
    
    double nprb = 0.0;
    
//    Speciation *s = modelPtr->getActiveSpeciation();
//    s->setAllBDFossParams();
//    double lambda = s->getBDSSSpeciationRateLambda();
//    double mu = s->getBDSSExtinctionRateMu();
//    double fossRate = s->getBDSSFossilSampRatePsi();
//    double sppSampRate = s->getBDSSSppSampRateRho();
    
    double sppSampRate = 1;
    vector<double> lambda, mu, fossRate;
    for(int i = 0; numIntervals; i++){
        lambda.push_back(1);
        mu.push_back(0.1);
        fossRate.push_back(2);
    }
    
    nprb = getFossilRangeGraphSkylineProb(lambda, mu, fossRate, sppSampRate, originTime);
    return nprb;
    
}

//FBD process augmenting the start and end of species skyline model

double FossilRangeGraphSkyline::getFossilRangeGraphSkylineProb(std::vector<double> lambda, std::vector<double> mu, std::vector<double> fossRate, double sppSampRate, double ot){
    if(runUnderPrior)
        return 0.0;
    
    double  nprb = 0.0;
    
    return nprb;
}

//END
