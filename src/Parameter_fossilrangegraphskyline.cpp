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

FossilRangeGraphSkyline::FossilRangeGraphSkyline(MbRandom *rp, Model *mp, int nf, int nl, vector<Calibration *> clb, int ni, vector<Calibration *> ints, bool rnp, bool fxFRG, int expMode, int fbdLk) : Parameter(rp, mp){

    name = "FRGS";
    numFossils = nf;
    numLineages = nl;
    numIntervals = ni + 1; // user defined intervals + the interval incorporating the ancient bound
    runUnderPrior = rnp;
    fbdLikelihood = fbdLk;
    conditionOnSurvival = 1;
    numExtinctLineages = 0;
    originTime = 0.0;
    originInterval = 0;
    ancientBound = 1000.0;
    fixOrigin = 0;
    orderStartStopTimes = 0;
    fixFRG = fxFRG; //1: fix start and end range times to FAs and LAs
    fixStart = 0;
    fixStop = 0;
    speedy = 1;
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
    cout << "Using likelihood function " << fbdLikelihood << endl;
    
}

FossilRangeGraphSkyline::~FossilRangeGraphSkyline(void){
    
}

FossilRangeGraphSkyline& FossilRangeGraphSkyline::operator=(const FossilRangeGraphSkyline &frgs) {
    
    if (this != &frgs)
        clone(frgs);
    return *this;
    
}

//TODO: I think you need to add other stuff here
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
    
    if(fbdLikelihood == 3){
        
        int numMoves = 3;
        int v = (int)(ranPtr->uniformRv() * numMoves);
        
        if(v == 0)
           updateLineageBi();
        else if (v == 1)
            updateLineageOi();
        else
            updateLineageDi();
        
    } else{
        if(ranPtr->uniformRv() < 0.5)
            updateLineageStartTimes();
        else
            updateLineageStopTimes();
    }
    
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
    
    double start = 0;
    double end = 0;
    int fossils = 0;
    
    int intid = 1;
    for(int i = 0; i < ints.size(); i++){
        Calibration *h = ints[i];
        start = h->getIntervalStart();
        end = h->getIntervalEnd();
        fossils = h->getIntervalFossils();
        
        Interval *interval = new Interval(start, end, fossils, intid, 0, 0);
        intervals.push_back(interval);
        
        intid ++;
    }
    
    Interval *interval = new Interval(ancientBound, start, 0, intid, 0, 0);
    intervals.push_back(interval);
    
}   

void FossilRangeGraphSkyline::printIntervalVariables(){
    
    for(int i = 0; i < numIntervals; i++){
        Interval *interval = intervals[i];
        cout << "Interval index: " << i << endl;
        cout << "Interval ID: " << interval->getIntervalID() << endl;
        cout << "Interval start: " << interval->getIntervalStart() << endl;
        cout << "Interval end: " << interval->getIntervalEnd() << endl;
        cout << "Number of fossils: " << interval->getIntervalFossils() << endl;
        cout << "Total sum of range durations (Ls): " << interval->getIntervalSumRangeLengths() << endl;
        cout << "Number of FAs & LAs (kappa'): " << interval->getIntervalKappaPrime() << endl << endl;
    }
}

void FossilRangeGraphSkyline::initializeIntervalConstants(){
    
    for(int i = 0; i < numIntervals; i++){
        intervalAs.push_back(0.01);
        intervalBs.push_back(0.01);
        intervalPs.push_back(0.01);
        intervalQs.push_back(0.01);
        intervalQts.push_back(0.01);
    }
}

void FossilRangeGraphSkyline::createFossilRangeSkylineVector(vector<Calibration *> clb){
    
    int frid=1;
    for(int c = 0; c < clb.size(); c++){
        Calibration *p = clb[c];
        double fa = p->getFirstAppearance();
        double la = p->getLastAppearance();
        double at = p->getAttachmentTime();
        double et = p->getEndTime();
        bool e = p->getIsExtant();
        bool eo = p->getIsExtantOnly();
        
        // vector for gamma interactions
        std::vector<bool> gi(numLineages, 0);
        
        // vector for sub branch lengths
        std::vector<double> ls(numIntervals, 0.0);
        
        FossilRangeSkyline *frs = new FossilRangeSkyline(fa, la, at, et, e, eo, frid, gi, ls);
        
        if(fbdLikelihood == 3)
            frs->setKappaFromCal(p->getKappa());
        
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
            //fr->setLineageStop(fr->getLastAppearance());
            fr->setLineageStop(fr->getEndTime());
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
    
    // assign fas, las, birth and death time to intervals
    for(int f = 0; f < numLineages; f++){
        FossilRangeSkyline *fr = fossilRangesSkyline[f];
        fr->setFossilRangeBirthInterval(assignInterval(fr->getLineageStart()));
        fr->setFossilRangeDeathInterval(assignInterval(fr->getLineageStop()));
        fr->setFossilRangeFirstAppearanceInterval(assignInterval(fr->getFirstAppearance()));
        fr->setFossilRangeLastAppearanceInterval(assignInterval(fr->getLastAppearance()));
    }
    
    // calculate L_s & kappa prime
    calculateIntervalSumRanges();
    
    if(fbdLikelihood == 3)
        initializeIntervalSubBranchLengths();
    
    bool printIntVariables = 1;
    
    if(printIntVariables)
        printIntervalVariables();
    
    if(printInitialFossilRangeSkylineVariables)
        printFossilRangeSkylineVariables();
    
}

// calculate L_s & kappa prime
void FossilRangeGraphSkyline::calculateIntervalSumRanges(){
    
    double fa, la, start, end, duration, sumInterval;
    int kappaPrime;
    
    for(int i = 0; i < numIntervals; i++){
        
        Interval *interval = intervals[i];
        
        start = interval->getIntervalStart();
        end = interval->getIntervalEnd();
        duration = start - end;
        
        sumInterval = 0;
        kappaPrime = 0;
        
        for(int f = 0; f < numLineages; f++){
            
            FossilRangeSkyline *fr = fossilRangesSkyline[f];
            fa = fr->getFirstAppearance();
            la = fr->getLastAppearance();
            
            // range through taxa
            if(fa > start && end > la){
                sumInterval += duration;
            } // first and last appearance occur within the interval (singleton)
            else if (fa < start && fa > end && la < start && la > end){
                sumInterval += fa - la;
                if(fa == la)
                    kappaPrime += 1;
                else
                    kappaPrime += 2;
            }
            // first appearance within interval
            else if (fa < start && fa >= end){
                sumInterval += fa - end;
                if(fa != 0)
                    kappaPrime += 1;
            }
            // last appearance within interval
            else if (la <= start && (la > end || (la == 0 && end == 0) )){
                sumInterval += start - la;
                if(la != 0)
                    kappaPrime += 1;
            }
        }
        interval->setIntervalSumRangeLengths(sumInterval);
        interval->setIntervalKappaPrime(kappaPrime);
    }
}

// initialise L_S and changeable o_i
void FossilRangeGraphSkyline::initializeIntervalSubBranchLengths(){
    
    for(int f = 0; f < numLineages; f++){
        
        FossilRangeSkyline *fr = fossilRangesSkyline[f];
        
        double bi = fr->getLineageStart(); // bi
        double di = fr->getLineageStop();  // di
        //double oi = fr->getFirstAppearance(); // oi
        
        int biInt = fr->getFossilRangeBirthInterval();
        int diInt = fr->getFossilRangeDeathInterval();
        //int oiInt = fr->getFossilRangeFirstAppearanceInterval();
        
        std::vector<double> subBrLens;
        
        for(int i = 0; i < numIntervals; i++){
            
            Interval *interval = intervals[i];
            
            double start = interval->getIntervalStart();
            double end = interval->getIntervalEnd();
            
            //TODO: delete this
            //if(i == oiInt){
            //    if(biInt == oiInt)
            //        fr->setChangeableOi(bi);
            //    else
            //        fr->setChangeableOi(start);
            //}
            
            int kappaS = fr->getFossilRangeKappaS(i);
            
            double subBranchLength;
            
            if(kappaS > 0){
                if(biInt == i && diInt == i)
                    subBranchLength = bi - di;
                else if(biInt == i)
                    subBranchLength = bi - end;
                else if (diInt == i)
                    subBranchLength = start - di;
                else
                    subBranchLength = start - end;
                
                fr->setFossilRangeSubBranchLength(i, subBranchLength);
                
            } else { // this is redundant
                
                fr->setFossilRangeSubBranchLength(i, 0.0);
                
            }
        }
    }
    
}

// calculate changeoable o_i (this isn't used)
void FossilRangeGraphSkyline::calculateChangeableOi(int it){
    
    FossilRangeSkyline *fr = fossilRangesSkyline[it];
    
    int i = fr->getFossilRangeFirstAppearanceInterval();
    double start = intervals[i]->getIntervalStart();
    double bi = fr->getLineageStart();
    
    if(bi > start)
        fr->setChangeableOi(start);
    else fr->setChangeableOi(bi);
    
}

// recalculate L_S. This could probably be more efficient.
void FossilRangeGraphSkyline::recalculateIntervalSubBranchLengths(int it){
    
    FossilRangeSkyline *fr = fossilRangesSkyline[it];
    
    double bi = fr->getLineageStart();
    double di = fr->getLineageStop();
    //double oi = fr->getFirstAppearance();
    
    int biInt = fr->getFossilRangeBirthInterval();
    int diInt = fr->getFossilRangeDeathInterval();
    //int oiInt = fr->getFossilRangeFirstAppearanceInterval();
    
    for(int i = 0; i < numIntervals; i++){
        
        Interval *interval = intervals[i];
        
        double start = interval->getIntervalStart();
        double end = interval->getIntervalEnd();
        
        int kappaS = fr->getFossilRangeKappaS(i);
        
        double subBranchLength;
        
        if(kappaS > 0){
            if(biInt == i && diInt == i)
                subBranchLength = bi - di;
            else if(biInt == i)
                subBranchLength = bi - end;
            else if (diInt == i)
                subBranchLength = start - di;
            else
                subBranchLength = start - end;
            
            fr->setFossilRangeSubBranchLength(i, subBranchLength);
            
        } else {
            
            fr->setFossilRangeSubBranchLength(i, 0.0);
            
        }
    }
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
                if(zj > zf && zf > bj){
                    g++;
                    fr->setFossilRangeGammaInteractions(j, 1);
                }
            }
        }
        // note gamma associated with the origin is handled by the likelihood function

        fr->setFossilRangeBrGamma(g);
    }
}

// gamma trick
void FossilRangeGraphSkyline::recountFossilRangeAttachNumsSpeedy(int i){
    
    FossilRangeSkyline *fr = fossilRangesSkyline[i];
    double zf = fr->getLineageStart();
    double bf = fr->getLineageStop();
    
    int gf;
    
    bool newInteraction = 0;
    bool oldInteraction = 0;
    
    for(int j = 0; j < numLineages; j++){
        
        gf = fr->getFossilRangeBrGamma();
        
        if(i != j){
            
            FossilRangeSkyline *o = fossilRangesSkyline[j];
            double zj = o->getLineageStart();
            double bj = o->getLineageStop();
            int gj = o->getFossilRangeBrGamma();
            
            // first deal with range f
            if(zj > zf && zf > bj)
                newInteraction = 1;
            else
                newInteraction = 0;
            
            // what was the previous interaction?
            oldInteraction = fr->getFossilRangeGammaInteractions(j);
            
            if(oldInteraction != newInteraction){
                if(newInteraction == 0)
                    gf = gf - 1;
                else
                    gf = gf + 1;
                fr->setFossilRangeBrGamma(gf);
                fr->setFossilRangeGammaInteractions(j, newInteraction);
            } // else do nothing
            
            // then deal with range j
            if(zf > zj && zj > bf)
                newInteraction = 1;
            else
                newInteraction = 0;
            
            // what was the previous interaction?
            oldInteraction = o->getFossilRangeGammaInteractions(i);
            
            if(oldInteraction != newInteraction){
                if(newInteraction == 0)
                    gj = gj - 1;
                else
                    gj = gj + 1;
                o->setFossilRangeBrGamma(gj);
                o->setFossilRangeGammaInteractions(i, newInteraction);
            } // else do nothing
            
        }
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
    
    //if(time == 0.0)
      //  return -1;
    
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
    
    double gamma = 0;
    
    for(int f = 0; f < numLineages; f++){
        FossilRangeSkyline *fr = fossilRangesSkyline[f];
        cout << "Fossil range ID: " << fr->getFossilRangeID() << endl;
        cout << "First appearance: " << fr->getFirstAppearance() << " sampled in interval " << fr->getFossilRangeFirstAppearanceInterval() << endl;
        cout << "Last appearance: " << fr->getLastAppearance() << " sampled in interval " << fr->getFossilRangeLastAppearanceInterval() << endl;
        cout << "Lineage start: " << fr->getLineageStart() << " sampled in interval " << fr->getFossilRangeBirthInterval() << endl;
        cout << "Lineage end: " << fr->getLineageStop() << " sampled in interval " << fr->getFossilRangeDeathInterval() << endl;
        cout << "Is extant: " << fr->getIsExtant() << endl;
        cout << "Gamma: " << fr->getFossilRangeBrGamma() << endl;
        gamma += fr->getFossilRangeBrGamma();
        
        if(fbdLikelihood == 3){
            cout << "Kappa S: ";
            for(int i = 0; i < numIntervals; i++){
                if(i == (numIntervals - 1))
                    cout << fr->getFossilRangeKappaS(i) << endl;
                else
                    cout << fr->getFossilRangeKappaS(i) << ", ";
            }
            cout << "Sub branch length: ";
            for(int i = 0; i < numIntervals; i++){
                if(i == (numIntervals - 1))
                    cout << fr->getFossilRangeSubBranchLength(i) << endl;
                else
                    cout << fr->getFossilRangeSubBranchLength(i) << ", ";
            }
            cout << endl;
            //cout << "Chageable Oi: " << fr->getChangeableOi() << endl << endl;
        }
    }
    
    cout << "Total  gamma = " << gamma << endl;
    
    bool printGammaInteractions = 0;
    
    if(printGammaInteractions){
        
        double interactions = 0;
        
        for(int f = 0; f < numLineages; f++){
            cout << "Gamma interactions for range : " << f << endl;
            FossilRangeSkyline *fr = fossilRangesSkyline[f];
            for(int j = 0; j < numLineages; j++){
                cout << "Range : " << j << " = " << fr->getFossilRangeGammaInteractions(j) << endl;
                if(fr->getFossilRangeGammaInteractions(j))
                    interactions++;
            }
        }
        cout << "Total  interactions = " << gamma << endl;
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
    
    //TODO: this doesn't need to be here
//    SpeciationSkyline *s = modelPtr->getActiveSpeciationSkyline();
//    s->setAllBDFossParams();
//    std::vector<double> rho = s->getSppSampRateRho();
//    std::vector<double> lambda = s->getSpeciationRates();
//    std::vector<double> mu = s->getExtinctionRates();
//    std::vector<double> fossRate = s->getFossilSampRates();
    
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
        if(speedy)
            recountFossilRangeAttachNumsSpeedy((*it));
        else
            recountFossilRangeAttachNums();
        fr->setFossilRangeBirthInterval(assignInterval(newStart));
        
        // recalculate the FRG probability
        double newLike = getFossilRangeGraphSkylineProb();
        
        // calculate the likelihood/prior ratio
        double lnLikeRatio = newLike - oldLike;
        double r = modelPtr->safeExponentiation(lnLikeRatio + c);
        
        if(ranPtr->uniformRv() > r){
            fr->setLineageStart(oldStart);
            redefineOriginTime();
            if(speedy)
                recountFossilRangeAttachNumsSpeedy((*it));
            else
                recountFossilRangeAttachNums();
            fr->setFossilRangeBirthInterval(assignInterval(oldStart));
            currentFossilRangeGraphSkylineLnL = oldLike;
        }
    }
    
    // debugging gamma code
    //double gamma = 0;
    //for(int f = 0; f < numLineages; f++){
    // FossilRangeSkyline *fr = fossilRangesSkyline[f];
    // gamma += fr->getFossilRangeBrGamma();
    //}
    // cout << "Total  interactions = " << gamma << endl;
    
    return 0.0;
}

double FossilRangeGraphSkyline::updateLineageStopTimes(){
    
    //TODO: this doesn't need to be here
//    SpeciationSkyline *s = modelPtr->getActiveSpeciationSkyline();
//    s->setAllBDFossParams();
//    std::vector<double> rho = s->getSppSampRateRho();
//    std::vector<double> lambda = s->getSpeciationRates();
//    std::vector<double> mu = s->getExtinctionRates();
//    std::vector<double> fossRate = s->getFossilSampRates();
    
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
        double oldLike = currentFossilRangeGraphSkylineLnL;
        
        // propose new values
        double newEnd = 0.0;
        double c = 0.0;
        
        double rv = ranPtr->uniformRv();
        double tv = log(2.0); // tuningVal
        c = doAScaleMove(newEnd, oldEnd, tv, 0.0, la, rv);
        
        // redefine values
        fr->setLineageStop(newEnd);
        if(speedy)
            recountFossilRangeAttachNumsSpeedy((*it));
        else
            recountFossilRangeAttachNums();
        fr->setFossilRangeDeathInterval(assignInterval(newEnd));
        
        // recalculate the FRG probability
        double newLike = getFossilRangeGraphSkylineProb();
        
        // calculate the likelihood/prior ratio
        double lnLikeRatio = newLike - oldLike;
        double r = modelPtr->safeExponentiation(lnLikeRatio + c);
        
        if(ranPtr->uniformRv() > r){
            fr->setLineageStop(oldEnd);
            if(speedy)
                recountFossilRangeAttachNumsSpeedy((*it));
            else
                recountFossilRangeAttachNums();
            fr->setFossilRangeDeathInterval(assignInterval(oldEnd));
            currentFossilRangeGraphSkylineLnL = oldLike;
        }
    }
    
    return 0.0;
}

double FossilRangeGraphSkyline::updateLineageBi(){
    
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
        double fa = fr->getFirstAppearance(); // oi
        double oldStart = fr->getLineageStart(); // bi
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
        
        if(speedy)
            recountFossilRangeAttachNumsSpeedy((*it));
        else
            recountFossilRangeAttachNums();
        
        fr->setFossilRangeBirthInterval(assignInterval(newStart));
        
        recalculateIntervalSubBranchLengths((*it));
        
        // recalculate the FRG probability
        double newLike = getFossilRangeGraphSkylineProb();
        
        // calculate the likelihood/prior ratio
        double lnLikeRatio = newLike - oldLike;
        double r = modelPtr->safeExponentiation(lnLikeRatio + c);
        
        if(ranPtr->uniformRv() > r){
            fr->setLineageStart(oldStart);
            redefineOriginTime();
            if(speedy) recountFossilRangeAttachNumsSpeedy((*it));
            else recountFossilRangeAttachNums();
            fr->setFossilRangeBirthInterval(assignInterval(oldStart));
            recalculateIntervalSubBranchLengths((*it));
            currentFossilRangeGraphSkylineLnL = oldLike;
        }
    }
    
    return 0.0;
}

double FossilRangeGraphSkyline::updateLineageOi(){
    
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
        
        int i = fr->getFossilRangeFirstAppearanceInterval();
        double bi = fr->getLineageStart(); // bi
        double di = fr->getLineageStop(); // di
        
        Interval *interval = intervals[i];
        double x = interval->getIntervalStart();
        double y = interval->getIntervalEnd();
        
        double min, max;
        if(di > y) min = di;
        else min = y;
        if(bi < x) max = bi;
        else max = x;
        
        double oldFa = fr->getFirstAppearance(); //oi
        double oldLike = currentFossilRangeGraphSkylineLnL;
        
        // propose new values
        double newFa = 0.0;
        double c = 0.0;
        
        double rv = ranPtr->uniformRv();
        double tv = 2.0; // tuningVal
        c = doAScaleMove(newFa, oldFa, tv, min, max, rv);
        
        // redefine values
        fr->setFirstAppearance(newFa);
        
        //TODO: delete
        //redefineOriginTime();
        //if(speedy)
        //    recountFossilRangeAttachNumsSpeedy((*it));
        //else
        //    recountFossilRangeAttachNums();
        //fr->setFossilRangeBirthInterval(assignInterval(newStart));
        
        //TODO: I don't think this needs to be here
        //recalculateIntervalSubBranchLengths((*it));
        
        // recalculate the FRG probability
        double newLike = getFossilRangeGraphSkylineProb();
        
        // calculate the likelihood/prior ratio
        double lnLikeRatio = newLike - oldLike;
        double r = modelPtr->safeExponentiation(lnLikeRatio + c);
        
        if(ranPtr->uniformRv() > r){
            
            fr->setFirstAppearance(oldFa);
            //recalculateIntervalSubBranchLengths((*it));
            currentFossilRangeGraphSkylineLnL = oldLike;
            
            //TODO: delete
            //fr->setLineageStart(oldStart);
            //redefineOriginTime();
            //if(speedy)
            //    recountFossilRangeAttachNumsSpeedy((*it));
            //else
            //    recountFossilRangeAttachNums();
            //fr->setFossilRangeBirthInterval(assignInterval(oldStart));
        }
    }
    
    return 0.0;
}

double FossilRangeGraphSkyline::updateLineageDi(){
    
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
        
        double max;
        
        int i = fr->getFossilRangeFirstAppearanceInterval();
        int j = fr->getFossilRangeLastAppearanceInterval();
        
        if(i == j)
            max = fr->getFirstAppearance();
        else max = intervals[j]->getIntervalStart();
        
        // define old values
        double oldEnd = fr->getLineageStop(); // di
        double oldLike = currentFossilRangeGraphSkylineLnL;
        
        // propose new values
        double newEnd = 0.0;
        double c = 0.0;
        
        double rv = ranPtr->uniformRv();
        double tv = log(2.0); // tuningVal
        c = doAScaleMove(newEnd, oldEnd, tv, 0.0, max, rv);
        
        // redefine values
        fr->setLineageStop(newEnd);
        
        if(speedy)
            recountFossilRangeAttachNumsSpeedy((*it));
        else
            recountFossilRangeAttachNums();
        
        fr->setFossilRangeDeathInterval(assignInterval(newEnd));
        
        recalculateIntervalSubBranchLengths((*it));
        
        // recalculate the FRG probability
        double newLike = getFossilRangeGraphSkylineProb();
        
        // calculate the likelihood/prior ratio
        double lnLikeRatio = newLike - oldLike;
        double r = modelPtr->safeExponentiation(lnLikeRatio + c);
        
        if(ranPtr->uniformRv() > r){
            fr->setLineageStop(oldEnd);
            if(speedy) recountFossilRangeAttachNumsSpeedy((*it));
            else recountFossilRangeAttachNums();
            fr->setFossilRangeDeathInterval(assignInterval(oldEnd));
            recalculateIntervalSubBranchLengths((*it));
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
        ss << "\tgamma_f(FR_" << frgID << ")";
    }
    
    string ni = ss.str();
    return ni;
}

string FossilRangeGraphSkyline::getFossilRangeSkylineInfoParamList(void){
    
    stringstream ss;
    FossilRangeSkyline *fr = NULL;
    
    for(int f = 0;  f < fossilRangesSkyline.size(); f++){
        fr = fossilRangesSkyline[f];
        ss << "\t" << fr->getFirstAppearance(); // first appearance
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
    
    double nprb = getFossilRangeGraphSkylineProb();
    return nprb;
    
}

//FBD process augmenting the start and end of species skyline model
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
    
    if(rho > 0)
        nprb += (numLineages - numExtinctLineages) * log(rho); // you can exclude (1 - rho) ^ (n - m - l) -- THIS IS INCORRECT;
    
    if(conditionOnSurvival) //TODO: double check this
        nprb -= log(1 - fbdSkylinePfxn(lambda, mu, fossRate, sppSampRate, originInterval, originTime));
    
    // fossils
    for(int i = 0; i < intervals.size(); i++){
        
        Interval *interval = intervals[i];
        
        if(fbdLikelihood == 1){
            
            int numFossils = interval->getIntervalFossils();
            
            if(numFossils > 0) nprb += numFossils * log(fossRate[i]);
            
        } else if (fbdLikelihood == 2) {
            
            int kappaPrime = interval->getIntervalKappaPrime();
            
            if(kappaPrime > 0) nprb += kappaPrime * log(fossRate[i]);
            
            nprb += fossRate[i] * interval->getIntervalSumRangeLengths();
            
        } else if (fbdLikelihood == 3){
            
            for(int f = 0; f < fossilRangesSkyline.size(); f++){
                
                FossilRangeSkyline *fr = fossilRangesSkyline[f];
                
                int kappaS = fr->getFossilRangeKappaS(i);
                
                if(kappaS > 0){
                    
                    double subBr = fr->getFossilRangeSubBranchLength(i);
                    
                    nprb += log( exp(fossRate[i] * subBr) * (1 - exp(-fossRate[i] * subBr)) );
                    
                }
            }
        }
    }
    
    //TODO: keep track of this elsewhere
    int numUnsampled = 0;
    
    // for each range
    for(int f = 0; f < fossilRangesSkyline.size(); f++){
        
        FossilRangeSkyline *fr = fossilRangesSkyline[f];
        
        double bi, di, oi;
        
        bi = fr->getLineageStart(); // bi
        di = fr->getLineageStop();  // di
        oi = fr->getFirstAppearance(); // oi
        
        //if(fbdLikelihood == 3)
        //    oi = fr->getChangeableOi();
        //else
        //    oi = fr->getFirstAppearance();
        
        int bint = fr->getFossilRangeBirthInterval();
        int dint = fr->getFossilRangeDeathInterval();
        int oint = fr->getFossilRangeFirstAppearanceInterval();
        
       if(bi == originTime)
            nprb += log( 1 );
        else
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
        
        for(int i = oint; i < bint; i++){
            nprb += intervalQs[i];
        }
        for(int i = dint; i < oint ; i++){
            nprb += intervalQts[i];
        }
        
    }
    
    //TODO: take care of this else where
    nprb += numUnsampled * log(1 - rho);
    
    currentFossilRangeGraphSkylineLnL = nprb;
    
    return nprb;
}

// Ai is analagous to c1
double FossilRangeGraphSkyline::fbdSkylineAfxn(std::vector<double> b, std::vector<double> d, std::vector<double> psi, int i){
    
    //TODO: skyline note: double check this should be abs
    //double v = fabs( sqrt( ( (b-d-psi) * (b-d-psi) ) + 4*b*psi) );
    double Ai = fabs( sqrt( ( (b[i] - d[i] - psi[i]) * (b[i] - d[i] - psi[i]) ) + 4 * b[i] * psi[i]) );
    return Ai;
}

// Bi is analagous to c2
double FossilRangeGraphSkyline::fbdSkylineBfxn(std::vector<double> b, std::vector<double> d, std::vector<double> psi, std::vector<double> rho, int i){
    
    //TODO: delete/clarify some of this
    // t =  interval minimum
    //if(i == 0)
      //  return 1;
    
    // the equations simplify when rho = 1, and 0 (but should still produce the same results)
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

// okay for large values of t
double FossilRangeGraphSkyline::fbdSkylineQTildaFxnLog(std::vector<double> b, std::vector<double> d, std::vector<double> psi, std::vector<double> rho, int i, double t){
    
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

// this isn't used
double FossilRangeGraphSkyline::fbdSkylineQTildaFxnLogSimplified(std::vector<double> b, std::vector<double> d, std::vector<double> psi, std::vector<double> rho, int i, double t, double q){
    
    Interval *interval = intervals[i];
    double ti = interval->getIntervalEnd();
    
    double qt = 0.5 * ( q + (-(b[i] + d[i] + psi[i]) * (t - ti)) );
    
    return qt;
}

double FossilRangeGraphSkyline::exampleRevBayesPfxn(std::vector<double> l, std::vector<double> m, std::vector<double> psi, std::vector<double> rho, int i, double t){
    
    if (t == 0 || i == -1) // note revbayes: doesn't have i == -1
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
    
    //TODO: what's all this? interval constants
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

//FBD process augmenting the start and end of species (this is for cross validation functions with numIntervals = 1)
//TODO: change this to void
double FossilRangeGraphSkyline::getFossilRangeGraphProb(std::vector<double> b, std::vector<double> d, std::vector<double> s, std::vector<double> r, double ot){
    if(runUnderPrior)
        return 0.0;
    
    double nprb = 0.0;
    
    double lambda = b[0];
    double mu = d[0];
    double fossRate = s[0];
    double sppSampRate = r[0];
    
    nprb = numFossils * log(fossRate);
    nprb += numExtinctLineages * log(mu);
    nprb -= log(lambda * (1-fbdPFxn(lambda,mu,fossRate,sppSampRate,ot)) );
    
    int numUnsampled = 0;
    
    for(int f=0; f < fossilRangesSkyline.size(); f++){
        
        FossilRangeSkyline *fr = fossilRangesSkyline[f];
        
        double bi = fr->getLineageStart(); // bi
        double di = fr->getLineageStop();  // di
        double oi = fr->getFirstAppearance(); // oi
        
        //TAKE CARE OF THIS ELSEWHERE
        double yi = fr->getLastAppearance(); // yi
        if (di == 0 && yi != 0)
            numUnsampled += 1;
        
        if(bi == originTime)
            nprb += log( lambda * 1.0 );
        else
            nprb += log( lambda * fr->getFossilRangeBrGamma() );
        
        double rangePr = 0;
        
        rangePr += fbdQTildaFxnLog(lambda, mu, fossRate, sppSampRate, oi);
        rangePr -= fbdQTildaFxnLog(lambda, mu, fossRate, sppSampRate, di);
        
        rangePr += fbdQFxnLog(lambda, mu, fossRate, sppSampRate, bi);
        rangePr -= fbdQFxnLog(lambda, mu, fossRate, sppSampRate, oi);
        
        nprb += rangePr;
    }
    
    //TODO: take care of this else where
    nprb += (numLineages - numExtinctLineages - numUnsampled) * log(sppSampRate);
    nprb += numUnsampled * log(1 - sppSampRate);
    
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
    
    bool useAlt = 1;
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
    
    //TODO: add this to a seperate fxn:
    //double q = fbdQFxnLog(b, d, psi, rho, t);
    //double v = (q + (- (b + d + psi)*t) ) * 0.5;
    
    return v;
}

//END
