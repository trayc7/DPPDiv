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
#include "Parameter_speciationskyline.h"
#include "MbRandom.h"
#include "Model.h"
#include "util.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <algorithm>

using namespace std;

FossilRangeGraphSkyline::FossilRangeGraphSkyline(MbRandom *rp, Model *mp, int nf, int nl, vector<Calibration *> clb, int ni, vector<Calibration *> ints, bool rnp, bool fxFRG, int expMode) : Parameter(rp, mp){

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
    orderStartStopTimes = 0;
    fixFRG = fxFRG; //1: fix start and end range times to FAs and LAs
    //fixFRG = 1;
    fixStart = 0;
    fixStop = 0;
    if(expMode == 1){
        fixOrigin = 1;
        originTime = 1.233376;
        ancientBound = originTime;
        orderStartStopTimes = 1;
    }
    printInitialFossilRangeSkylineVariables = 1;
    createIntervalsVector(ints);
    createFossilRangeSkylineVector(clb);
    initializeFossilRangeSkylineVariables();
    initializeIntervalConstants();
    currentFossilRangeGraphSkylineLnL = 0.0;
    counter = 0; // debugging
    
    bool crossValidate = 0;
    if(crossValidate)
        crossValidateFBDSkylinefunctions();
    
    cout << "Number of lineages: " << numLineages << endl;
    cout << "Number of extinct ranges: " << numExtinctLineages << endl;
    cout << "Number of fossils: " << numFossils << endl;
    cout << "Number of intervals (inc ancient bound): " << numIntervals << endl;
    cout << "\nInitial origin time: " << originTime << endl;
    cout << "Fossil range graph skyline initialized" << endl;
    
}

FossilRangeGraphSkyline::~FossilRangeGraphSkyline(void){
    
}

FossilRangeGraphSkyline& FossilRangeGraphSkyline::operator=(const FossilRangeGraphSkyline &frgs) {
    
    if (this != &frgs)
        clone(frgs);
    return *this;
    
}

void FossilRangeGraphSkyline::clone(const FossilRangeGraphSkyline &frgs){

    numFossils = frgs.numFossils;
    numLineages = frgs.numLineages;
    numExtinctLineages = frgs.numExtinctLineages;
    originTime = frgs.originTime;
    numIntervals = frgs.numIntervals;
    
    for(int f = 0; f < fossilRangesSkyline.size(); f++){
        FossilRangeSkyline *frTo = fossilRangesSkyline[f];
        FossilRangeSkyline *frFrom = frgs.fossilRangesSkyline[f];
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
    
    if(orderStartStopTimes)
        orderRangeAges();
    
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
    
    for(int i = 0; i < numIntervals; i++){
        Interval *interval = intervals[i];
        cout << "Interval index: " << i << endl;
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
        
        FossilRangeSkyline *frs = new FossilRangeSkyline(fa, la, at, e, eo, frid);
        fossilRangesSkyline.push_back(frs);
        
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
//    if(fixStart){
//        for(int f = 0; f < numLineages; f++){
//            FossilRangeSkyline *fr = fossilRangesSkyline[f];
//            fr->setLineageStart(fr->getAttachmentTime());
//            fr->setLineageStop(fr->getLastAppearance());
//            fr->setFixStart(1);
//        }
//    }

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

// debugging code
void FossilRangeGraphSkyline::orderRangeAges(){
    
    double start, stop;
    std::vector<double> agesStart, agesStop;
    
    for(int f = 0; f < numLineages; f++){
        FossilRangeSkyline *fr = fossilRangesSkyline[f];
        start = fr->getLineageStart();
        agesStart.push_back(start);
        stop = fr->getLineageStop();
        agesStop.push_back(stop);
        
    }
    
    std::sort(agesStart.begin(), agesStart.end(), std::greater<double>());
    std::sort(agesStop.begin(), agesStop.end(), std::greater<double>());
    
    for(int f = 0; f < numLineages; f++){
        FossilRangeSkyline *fr = fossilRangesSkyline[f];
        fr->setLineageStart(agesStart[f]);
        fr->setLineageStop(agesStop[f]);
    }
    
    agesStart.clear();
    agesStop.clear();
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
    
    SpeciationSkyline *s = modelPtr->getActiveSpeciationSkyline();
    s->setAllBDFossParams();
    std::vector<double> rho = s->getSppSampRateRho();
    std::vector<double> lambda = s->getSpeciationRates();
    std::vector<double> mu = s->getExtinctionRates();
    std::vector<double> fossRate = s->getFossilSampRates();
    
    vector<int> rndFossilRangeIDs;
    for(int i = 0; i < fossilRangesSkyline.size(); i++)
        rndFossilRangeIDs.push_back(i);
    random_shuffle(rndFossilRangeIDs.begin(), rndFossilRangeIDs.end());
    
    for(vector<int>::iterator it = rndFossilRangeIDs.begin(); it != rndFossilRangeIDs.end(); it++){
        
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
        double tv = 2.0; // tuningVal
        c = doAScaleMove(newStart, oldStart, tv, fa, ancientBound, rv);
        
        // redefine values
        fr->setLineageStart(newStart);
        redefineOriginTime();
        recountFossilRangeAttachNums();
        fr->setFossilRangeBirthInterval(assignInterval(newStart));
        
        // recalculate the FRG probability
        //double newLike = getFossilRangeGraphSkylineProb(lambda, mu, fossRate, rho, originTime);
        double newLike = getFossilRangeGraphSkylineProb();
        
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
    return 0.0;
}

double FossilRangeGraphSkyline::updateLineageStopTimes(){
    
    SpeciationSkyline *s = modelPtr->getActiveSpeciationSkyline();
    s->setAllBDFossParams();
    std::vector<double> rho = s->getSppSampRateRho();
    std::vector<double> lambda = s->getSpeciationRates();
    std::vector<double> mu = s->getExtinctionRates();
    std::vector<double> fossRate = s->getFossilSampRates();
    
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
        //double newLike = getFossilRangeGraphSkylineProb(lambda, mu, fossRate, rho, originTime);
        double newLike = getFossilRangeGraphSkylineProb();
        
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

string FossilRangeGraphSkyline::getFossilRangeSkylineInfoParamNames(void){
    
    stringstream ss;
    FossilRangeSkyline *fr = NULL;
    
    for(int f = 0;  f < fossilRangesSkyline.size(); f++){
        fr = fossilRangesSkyline[f];
        int frgID = fr->getFossilRangeID();
        ss << "\ty_f(FR_" << frgID << ")"; // first appearance
        ss << "\tb_f(FR_" << frgID << ")"; // lineage start
        ss << "\tx_f(FR_" << frgID << ")"; // last appearance
        ss << "\td_f(FR_" << frgID << ")"; // lineage stop
    }
    for(int f = 0;  f < fossilRangesSkyline.size(); f++){
        fr = fossilRangesSkyline[f];
        int frgID = fr->getFossilRangeID();
        ss << "\ty_I(FR_" << frgID << ")"; // first appearance interval
        ss << "\tb_I(FR_" << frgID << ")"; // lineage start interval
        ss << "\td_I(FR_" << frgID << ")"; // lineage stop interval
    }
    for(int f = 0;  f < fossilRangesSkyline.size(); f++){
        fr = fossilRangesSkyline[f];
        int frgID = fr->getFossilRangeID();
        ss << "\tgamma_f(FR_" << frgID << ")"; //".nd" << nID << ")";
    }
    
    string ni = ss.str();
    return ni;
}

string FossilRangeGraphSkyline::getFossilRangeSkylineInfoParamList(void){
    
    stringstream ss;
    FossilRangeSkyline *fr = NULL;
    
    for(int f = 0;  f < fossilRangesSkyline.size(); f++){
        fr = fossilRangesSkyline[f];
        ss << "\t" << fr->getFirstAppearance(); // first appearance -- note this is fixed in this implementation
        ss << "\t" << fr->getLineageStart(); // lineage start
        ss << "\t" << fr->getLastAppearance(); // last appearance
        ss << "\t" << fr->getLineageStop(); // lineage stop
    }
    for(int f = 0;  f < fossilRangesSkyline.size(); f++){
        fr = fossilRangesSkyline[f];
        ss << "\t" << fr->getFossilRangeFirstAppearanceInterval(); // first appearance interval
        ss << "\t" << fr->getFossilRangeBirthInterval(); // lineage start interval
        ss << "\t" << fr->getFossilRangeDeathInterval(); // lineage stop interval
    }
    for(int f = 0;  f < fossilRangesSkyline.size(); f++){
        fr = fossilRangesSkyline[f];
        ss << "\t" << fr->getFossilRangeBrGamma();
    }
    
    string ni = ss.str();
    return ni;
}

double FossilRangeGraphSkyline::getActiveFossilRangeGraphSkylineProb(){
    if(runUnderPrior)
        return 0.0;
    
//    double nprb = 0.0;

//    SpeciationSkyline *s = modelPtr->getActiveSpeciationSkyline();
//    s->setAllBDFossParams();
//    std::vector<double> rho = s->getSppSampRateRho();
//    std::vector<double> lambda = s->getSpeciationRates();
//    std::vector<double> mu = s->getExtinctionRates();
//    std::vector<double> fossRate = s->getFossilSampRates();
    
    double nprb = getFossilRangeGraphSkylineProb();
    return nprb;
    
}

//FBD process augmenting the start and end of species skyline model
//double FossilRangeGraphSkyline::getFossilRangeGraphSkylineProb(std::vector<double> lambda, std::vector<double> mu, std::vector<double> fossRate, std::vector<double> sppSampRate, double ot){
double FossilRangeGraphSkyline::getFossilRangeGraphSkylineProb(){
    if(runUnderPrior)
        return 0.0;
    
    SpeciationSkyline *s = modelPtr->getActiveSpeciationSkyline();
    std::vector<double> sppSampRate = s->getSppSampRateRho();
    std::vector<double> lambda = s->getSpeciationRates();
    std::vector<double> mu = s->getExtinctionRates();
    std::vector<double> fossRate = s->getFossilSampRates();
    
    double  nprb = 0.0;
    
    bool frgProb = 1;
    if(frgProb){
        nprb = getFossilRangeGraphProb(lambda, mu, fossRate, sppSampRate, originTime);
        return nprb;
    }
    
    double rho = sppSampRate[0];
    
    nprb = (numLineages - numExtinctLineages) * log(rho); // double check you can just exclude (1 - rho) ^ (n - m - l);
    nprb -= log(1 - fbdSkylinePfxn(lambda, mu, fossRate, sppSampRate, originInterval, originTime));
    
    // fossils
    for(int i = 0; i < intervals.size(); i++){
        Interval *interval = intervals[i];
        int numFossils = interval->getIntervalFossils();
        if(numFossils > 0)
            nprb += numFossils * log(fossRate[i]);
    }
    
    // for each range
    for(int f = 0; f < fossilRangesSkyline.size(); f++){
        
        FossilRangeSkyline *fr = fossilRangesSkyline[f];
        
        double bi = fr->getLineageStart(); // bi
        double di = fr->getLineageStop();  // di
        double oi = fr->getFirstAppearance(); // oi
        int bint = fr->getFossilRangeBirthInterval();
        int dint = fr->getFossilRangeDeathInterval();
        int oint = fr->getFossilRangeFirstAppearanceInterval();
        
        nprb += log( fr->getFossilRangeBrGamma() );
        
        // speciation events
        if(bi != originTime)
            nprb += log(lambda[bint]);
        
        // extinction events
        if(!fr->getIsExtant())
            nprb += log(mu[dint]);
        
        nprb += fbdSkylineQTildaFxnLog(lambda, mu, fossRate, sppSampRate, oint, oi);
        nprb -= fbdSkylineQTildaFxnLog(lambda, mu, fossRate, sppSampRate, dint, di);
        
        nprb += fbdSkylineQfxnLog(lambda, mu, fossRate, sppSampRate, bint, bi);
        nprb -= fbdSkylineQfxnLog(lambda, mu, fossRate, sppSampRate, oint, oi);
        
//        for(int i = oint; i < bint; i++){
//            Interval *interval = intervals[i];
//            double t = interval->getIntervalStart();
//            nprb += fbdSkylineQfxnLog(lambda, mu, fossRate, sppSampRate, i, t);
//        }
//        
//        for(int i = dint; i < oint ; i++){
//            Interval *interval = intervals[i];
//            double t = interval->getIntervalStart();
//            nprb += fbdSkylineQTildaFxnLog(lambda, mu, fossRate, sppSampRate, i, t);
//        }
        
        for(int i = oint; i < bint; i++){
            nprb += intervalQs[i];
        }
        for(int i = dint; i < oint ; i++){
            nprb += intervalQts[i];
        }
        
    }
    
    currentFossilRangeGraphSkylineLnL = nprb;
    
    return nprb;
}

// Ai is analagous to c1
double FossilRangeGraphSkyline::fbdSkylineAfxn(std::vector<double> b, std::vector<double> d, std::vector<double> psi, int i){
    
    // skyline note: double check this should be abs
    //double v = fabs( sqrt( ( (b-d-psi) * (b-d-psi) ) + 4*b*psi) );
    double Ai = fabs( sqrt( ( (b[i] - d[i] - psi[i]) * (b[i] - d[i] - psi[i]) ) + 4 * b[i] * psi[i]) );
    return Ai;
}

// Bi is analagous to c2
double FossilRangeGraphSkyline::fbdSkylineBfxn(std::vector<double> b, std::vector<double> d, std::vector<double> psi, std::vector<double> rho, int i){
    
    // t =  interval minimum
    //if(i == 0)
      //  return 1;
    
    // the equations simplify when rho = 1, and 0 (but should still produce the same results
    // for fossil data, rho = 0 for all intevals except the present
    // note also that intervals increase towards the present in the maths notation
    // pi+1 (t_i) - this refers to the p for the NEXT interval and t_i is the minimum age of the current interval
    
    Interval *interval = intervals[i];
    double ti = interval->getIntervalEnd();
    
    double Ai = intervalAs[i];
    
    double Bi;
    
    Bi = (1 - 2 * (1 - rho[i]) * fbdSkylinePfxn(b, d, psi, rho, i-1, ti) ) * b[i];
    //Bi = (1 - 2 * (1 - rho[i]) * intervalPs[i-1] ) * b[i];
    Bi += d[i] + psi[i];
    Bi /= Ai;
    
    return Bi;
}

double FossilRangeGraphSkyline::fbdSkylinePfxn(std::vector<double> b, std::vector<double> d, std::vector<double> psi, std::vector<double> rho, int i, double t){
    
    if(i == -1 || t == 0)
        return 1;

    Interval *interval = intervals[i];
    double ti = interval->getIntervalEnd();
    
//    double Ai = fbdSkylineAfxn(b, d, psi, i);
//    double Bi = fbdSkylineBfxn(b, d, psi, rho, i);

    double Ai = intervalAs[i];
    double Bi = intervalBs[i];
    
    double p;
    
    p = b[i] + d[i] + psi[i];
    p -= Ai * ( (1 + Bi - (1 - Bi) * exp(Ai * (ti - t))) / (1 + Bi + (1 - Bi) * exp(Ai * (ti - t))) );
    p /= (2 * b[i]);
    
    return p;
}

double FossilRangeGraphSkyline::fbdSkylineABPfxnInterval(std::vector<double> b, std::vector<double> d, std::vector<double> psi, std::vector<double> rho, int i, double t){

    if(i == -1 || t == 0)
        return 1;
    
    Interval *interval = intervals[i];
    double ti = interval->getIntervalEnd();
    
    double Ai = fabs( sqrt( ( (b[i] - d[i] - psi[i]) * (b[i] - d[i] - psi[i]) ) + 4 * b[i] * psi[i]) );
    
    double Bi = (1 - 2 * (1 - rho[i]) * fbdSkylineABPfxnInterval(b, d, psi, rho, i-1, ti) ) * b[i];
    Bi += d[i] + psi[i];
    Bi /= Ai;
    
    double p;
    
    p = b[i] + d[i] + psi[i];
    p -= Ai * ( (1 + Bi - (1 - Bi) * exp(Ai * (ti - t))) / (1 + Bi + (1 - Bi) * exp(Ai * (ti - t))) );
    p /= (2 * b[i]);
    
    intervalAs[i] = Ai;
    intervalBs[i] = Bi;
    intervalPs[i] = p;
    
    return p;
}

// unsuitable for large values of t
double FossilRangeGraphSkyline::fbdSkylineQfxn(std::vector<double> b, std::vector<double> d, std::vector<double> psi, std::vector<double> rho, int i, double t){
 
    Interval *interval = intervals[i];
    double ti = interval->getIntervalEnd();
    
    double Ai = fbdSkylineAfxn(b, d, psi, i);
    double Bi = fbdSkylineBfxn(b, d, psi, rho, i);
    
    double q;
    
    q = 4 * exp (Ai * (ti - t) );
    q /= (exp(-Ai * (ti - t) ) * (1 - Bi) + (1 + Bi)) * (exp(-Ai * (ti - t) ) * (1 - Bi) + (1 + Bi));
    
    return q;
}

double FossilRangeGraphSkyline::fbdSkylineQfxnLog(std::vector<double> b, std::vector<double> d, std::vector<double> psi, std::vector<double> rho, int i, double t){
    
    Interval *interval = intervals[i];
    double ti = interval->getIntervalEnd();
    
    //double Ai = fbdSkylineAfxn(b, d, psi, i);
    //double Bi = fbdSkylineBfxn(b, d, psi, rho, i);
    double Ai = intervalAs[i];
    double Bi = intervalBs[i];
    
    double q;
    
    q = log(4) + (Ai * (ti - t)) - (2 * (log( (exp(Ai*(ti - t)) * (1-Bi)) + (1+Bi) )));
    
    return q;
}


// unsuitable for large values of t
double FossilRangeGraphSkyline::fbdSkylineQTildaFxn(std::vector<double> b, std::vector<double> d, std::vector<double> psi, std::vector<double> rho, int i, double t){

    double Ai = fbdSkylineAfxn(b, d, psi, i);
    double Bi = fbdSkylineBfxn(b, d, psi, rho, i);
    
    double f1a = 4 * exp(-t * (b[i] + d[i] + psi[i])) * exp(-t * Ai);
    double f1b = 4 * exp(-t * Ai) + (1 - Bi * Bi) * (1 - exp(-t*Ai)) * (1 - exp(-t*Ai));
    double f2a = (1 + Bi) * exp(-t * Ai) + (1 - Bi);
    double f2b = (1 - Bi) * exp(-t * Ai) + (1 + Bi);
    
    double qt = sqrt ( (f1a / f1b) * (f2a / f2b) );
    
    return qt;
}

// this should be okay for large values of t now
double FossilRangeGraphSkyline::fbdSkylineQTildaFxnLog(std::vector<double> b, std::vector<double> d, std::vector<double> psi, std::vector<double> rho, int i, double t){
    
    //double Ai = fbdSkylineAfxn(b, d, psi, i);
    //double Bi = fbdSkylineBfxn(b, d, psi, rho, i);
    
    double Ai = intervalAs[i];
    double Bi = intervalBs[i];
    
    Interval *interval = intervals[i];
    double ti = interval->getIntervalEnd();
    
    t -= ti;
    
    double f1aLog = log(4) + (-t * (b[i] + d[i] + psi[i])) + (-t * Ai);
    double f1b = 4 * exp(-t * Ai) + (1 - Bi * Bi) * (1 - exp(-t*Ai)) * (1 - exp(-t*Ai));
    double f2a = (1 + Bi) * exp(-t * Ai) + (1 - Bi);
    double f2b = (1 - Bi) * exp(-t * Ai) + (1 + Bi);
    
    double qt = 0.5 * (log( (1/f1b) * (f2a/f2b) ) + f1aLog);

    return qt;
}

// I don't actually use this
double FossilRangeGraphSkyline::fbdSkylineQTildaFxnLogSimplified(std::vector<double> b, std::vector<double> d, std::vector<double> psi, std::vector<double> rho, int i, double t, double q){
    
    Interval *interval = intervals[i];
    double ti = interval->getIntervalEnd();
    
    double qt = 0.5 * ( q + (-(b[i] + d[i] + psi[i]) * (t - ti)) );
    
    return qt;
}

double FossilRangeGraphSkyline::exampleRevBayesPfxn(std::vector<double> l, std::vector<double> m, std::vector<double> psi, std::vector<double> rho, int i, double t){
    
    if (t == 0 || i == -1) // rb: doesn't have i == -1
        return 1.0;
    
    // get the parameters
    double b = l[i]; // rb: [i-1];
    double d = m[i];
    double f = psi[i];
    double r = rho[i];
    
    Interval *interval = intervals[i];
    double ti = interval->getIntervalEnd();
    
    double diff = b - d - f;
    double bp   = b * f;
    double dt   = t - ti;
    
    double A = sqrt(diff * diff + 4.0 * bp);
    double B = ( (1.0 - 2.0 * (1.0 - r) * exampleRevBayesPfxn(l, m, psi, rho, i-1, ti) ) * b + d + f ) / A; // rb: i-1 = i+1
    
    double e = exp(A * dt);
    double tmp = b + d + f - A * (e * (1.0 + B) - (1.0 - B))/(e * (1.0 + B) + (1.0 - B));
    
    return tmp / (2.0 * b);
}

//move this
void FossilRangeGraphSkyline::initializeIntervalConstants(){
    
    for(int i = 0; i < numIntervals; i++){
        intervalAs.push_back(0.01);
        intervalBs.push_back(0.01);
        intervalPs.push_back(0.01);
        intervalQs.push_back(0.01);
        intervalQts.push_back(0.01);
    }
}

void FossilRangeGraphSkyline::setAllIntervalConstants(){
    
    SpeciationSkyline *s = modelPtr->getActiveSpeciationSkyline();
    std::vector<double> rho = s->getSppSampRateRho();
    std::vector<double> lambda = s->getSpeciationRates();
    std::vector<double> mu = s->getExtinctionRates();
    std::vector<double> psi = s->getFossilSampRates();
    
    int i = numIntervals - 1;
    Interval *interval = intervals[i];
    double ti = interval->getIntervalEnd();
    fbdSkylineABPfxnInterval(lambda, mu, psi, rho, i, ti);
    
    for(int i = 0; i < numIntervals; i++){
        Interval *interval = intervals[i];
        double s = interval->getIntervalStart();
        intervalQs[i] = ( fbdSkylineQfxnLog(lambda, mu, psi, rho, i, s));
        intervalQts[i] = ( fbdSkylineQTildaFxnLog(lambda, mu, psi, rho, i, s));
    }
}

void FossilRangeGraphSkyline::crossValidateFBDSkylinefunctions(){
    
    std::vector<double> lambda, mu, psi, rho;
    
    for(int i = 0; i < numIntervals; i++){
        lambda.push_back(0.4758981467371014);
        mu.push_back(0.31179500560019729);
        psi.push_back(2.0);
        rho.push_back(0.0);
    }
    rho[0] = 1.0;
    
    for(int i = 0; i < numIntervals; i++){
        Interval *interval = intervals[i];
        double ti = interval->getIntervalEnd();
        intervalAs[i] = ( fbdSkylineAfxn(lambda, mu, psi, i) );
        intervalBs[i] = ( fbdSkylineBfxn(lambda, mu, psi, rho, i) );
        intervalPs[i] = ( fbdSkylinePfxn(lambda, mu, psi, rho, i, ti));
        double s = interval->getIntervalStart();
        intervalQs[i] = ( fbdSkylineQfxnLog(lambda, mu, psi, rho, i, s));
        intervalQts[i] = ( fbdSkylineQTildaFxnLog(lambda, mu, psi, rho, i, s));
    }
    
    for(int i = 0; i < numIntervals; i++){
        cout << "p " << intervalPs[i] << endl;
    }
    
    for(int i = 0; i < numIntervals; i++){
        cout << "i " << exp(intervalQs[i]) << endl;
    }
    
    for(int i = 0; i < numIntervals; i++){
        cout << "i " << intervalQts[i] << endl;
    }
    //
    
    for(int i = 0; i < numIntervals; i++){
        intervalAs[i] = 1;
        intervalBs[i] = 1;
        intervalPs[i] = 1;
        intervalQs[i] = 1;
        intervalQts[i] = 1;
    }
    
//    for(int i = 0; i < numIntervals; i++){
//        cout << "i " << intervalPs[i] << endl;
//    }
//    
//    for(int i = 0; i < numIntervals; i++){
//        cout << "i " << intervalQs[i] << endl;
//    }
//    
//    for(int i = 0; i < numIntervals; i++){
//        cout << "i " << intervalQts[i] << endl;
//    }
    
    double ti = 10;
    int i = 10;
    fbdSkylineABPfxnInterval(lambda, mu, psi, rho, i, ti);
    
    for(int i = 0; i < numIntervals; i++){
        cout << "p " << intervalPs[i] << endl;
    }
    
    // interval constants
//    std::vector<double> intervalConstantsA, intervalConstantsB, intervalPsx, intervalRBPs;
//    for(int i = 0; i < numIntervals; i++){
//        Interval *interval = intervals[i];
//        double ti = interval->getIntervalEnd();
//        intervalConstantsA.push_back( fbdSkylineAfxn(lambda, mu, psi, i) );
//        intervalConstantsB.push_back( fbdSkylineBfxn(lambda, mu, psi, rho, i) );
//        intervalPsx.push_back( fbdSkylinePfxn(lambda, mu, psi, rho, i, ti));
//        intervalRBPs.push_back( exampleRevBayesPfxn(lambda, mu, psi, rho, i, ti));
//    }
    
//    cout << "Cross validating FBD functions..." << endl;
//    cout << "Constant A" << endl;
//    for(int i = 0; i < numIntervals; i++){
//        cout << "A for interval " << i << " = " << intervalConstantsA[i] << endl;
//    }
//    cout << "Constant B" << endl;
//    for(int i = 0; i < numIntervals; i++){
//        cout << "B for interval " << i << " = " << setprecision(15) << intervalConstantsB[i] << endl;
//    }
//    cout << "P(t) function" << endl;
//    for(int i = 0; i < numIntervals; i++){
//        cout << "P(t) for interval " << i << " = " << intervalPs[i] << endl;
//    }
//    
//    cout << "rb P(t) function" << endl;
//    for(int i = 0; i < numIntervals; i++){
//        cout << "rbP(t) for interval " << i << " = " << intervalRBPs[i] << endl;
//    }

    // other functions
//    double q = fbdSkylineQfxnLog(lambda, mu, psi, rho, 0, 2.1);
//    double qt = fbdSkylineQTildaFxnLog(lambda, mu, psi, rho, 0, 2.1);
//    
//    //cout << "Origin time = " << originTime << ", interval " << originInterval << endl;
//    cout << "Log q(t) for the origin = " << q << endl;
//    cout << "Log qt(t) for the origin = " << qt << endl;
//    
//    //double lk = getFossilRangeGraphSkylineProb(lambda, mu, psi, rho, originTime);
//    double lk = getFossilRangeGraphSkylineProb();
//    cout << "FRG skyline likelihood " << setprecision(7) << lk << endl;
//    
//    cout << "0.5 " << assignInterval(0.5) << endl;
//    cout << "1.3 " << assignInterval(1.3) << endl;
//    cout << "2.9 " << assignInterval(2.9) << endl;
//    cout << "4.6 " << assignInterval(4.6) << endl;
    
    exit(0);
}

//FBD process augmenting the start and end of species
double FossilRangeGraphSkyline::getFossilRangeGraphProb(std::vector<double> b, std::vector<double> d, std::vector<double> s, std::vector<double> r, double ot){
    if(runUnderPrior)
        return 0.0;
    
    double nprb = 0.0;
    
    double lambda = b[0];
    double mu = d[0];
    double fossRate = s[0];
    double sppSampRate = r[0];
    
    nprb = numFossils*log(fossRate);
    nprb += numExtinctLineages*log(mu);
    nprb -= log(lambda * (1-fbdPFxn(lambda,mu,fossRate,sppSampRate,ot)) );
    
    for(int f=0; f < fossilRangesSkyline.size(); f++){
        
        FossilRangeSkyline *fr = fossilRangesSkyline[f];
        
        double bi = fr->getLineageStart(); // bi
        double di = fr->getLineageStop();  // di
        double oi = fr->getFirstAppearance(); // oi
        
        nprb += log( lambda * fr->getFossilRangeBrGamma() );
        
        double rangePr = 0;
        
        rangePr += fbdQTildaFxnLog(lambda, mu, fossRate, sppSampRate, oi);
        rangePr -= fbdQTildaFxnLog(lambda, mu, fossRate, sppSampRate, di);
        
        rangePr += fbdQFxnLog(lambda, mu, fossRate, sppSampRate, bi);
        rangePr -= fbdQFxnLog(lambda, mu, fossRate, sppSampRate, oi);
        
        nprb += rangePr;
    }
    
    currentFossilRangeGraphSkylineLnL = nprb;
    
    return nprb;
}

double FossilRangeGraphSkyline::fbdC1Fxn(double b, double d, double psi){
    
    double v = fabs( sqrt( ( (b-d-psi) * (b-d-psi) ) + 4*b*psi) );
    return v;
}

double FossilRangeGraphSkyline::fbdC2Fxn(double b, double d, double psi, double rho){
    
    double v = -( ( b-d-(2*b*rho)-psi ) / (fbdC1Fxn(b,d,psi)) );
    return v;
}

double FossilRangeGraphSkyline::fbdC3Fxn(double b, double d, double psi, double rho){
    
    double v = b * (-psi + rho * (d + b * (-1 +rho) + psi ) );
    return v;
}

double FossilRangeGraphSkyline::fbdC4Fxn(double b, double d, double psi, double rho){
    
    double v = fbdC3Fxn(b,d,psi,rho) / (fbdC1Fxn(b,d,psi) * fbdC1Fxn(b,d,psi));
    return v;
}

double FossilRangeGraphSkyline::fbdPFxn(double b, double d, double psi, double rho, double t){
    
    double c1Val = fbdC1Fxn(b,d,psi);
    double c2Val = fbdC2Fxn(b,d,psi,rho);
    
    double expC1MinusC2 = exp(-c1Val * t) * (1.0 - c2Val);
    
    double eCfrac = (expC1MinusC2 - (1.0 + c2Val)) / (expC1MinusC2 + (1.0 + c2Val));
    double v = 1.0 + ((-(b - d - psi)) + (c1Val * eCfrac)) / (2.0 * b);
    
    return v;
}

double FossilRangeGraphSkyline::fbdQFxnLog(double b, double d, double psi, double rho, double t){
    
    double c1 = fbdC1Fxn(b,d,psi);
    double c2 = fbdC2Fxn(b,d,psi,rho);
    
    double f1 = log(4) + (-c1 * t);
    double f2 = 2 * (log( (exp(-c1*t) * (1-c2)) + (1+c2) ));
    
    double v = f1 - f2;
    
    return v;
}

// rho = 1
double FossilRangeGraphSkyline::fbdQTildaFxnLog(double b, double d, double psi, double rho, double t){
    
    bool useAlt = 0;
    if(useAlt)
        return fbdQTildaFxnLogAlt(b, d, psi, rho, t);
    
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

// rho < 1
double FossilRangeGraphSkyline::fbdQTildaFxnLogAlt(double b, double d, double psi, double rho, double t){
    
    double c1 = fbdC1Fxn(b,d,psi);
    double c2 = fbdC2Fxn(b,d,psi,rho);
    
    double f1a = 4 * exp(-t * c1) * exp(-t * (b + d + psi));
    double f1b = 4 * exp(-t * c1) + (1 - c2 * c2) * (1 - exp(-t * c1)) *(1 - exp(-t * c1));
    double f2a = (1 + c2) * exp(-t * c1) + (1 - c2);
    double f2b = (1 - c2) * exp(-t * c1) + (1 + c2);
    
    double v = 0.5 * log( (f1a/f1b) * (f2a/f2b) );
    
    return v;
}

//END
