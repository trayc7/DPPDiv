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
    printInitialFossilRangeVariables = 1;
    numExtinctLineages = 0;
    createFossilRangeVector(clb);
    initializeFossilRangeVariables();
    currentFossilRangeGraphLnL = 0.0;
    moves = 1; // 1: update lineage start or stop times; 2: update both
    proposal = 2; // proposal type 1=window, 2=scale, 3=slide
    
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
            stop = ranPtr->uniformRv(0.0,la);
            fr->setLineageStop(stop);
        }
    }

    //fixed frg test
//    fossilRanges[0]->setLineageStart(5.560851);
//    fossilRanges[0]->setLineageStop(2.469863);
//    fossilRanges[1]->setLineageStart(3.723645);
//    fossilRanges[1]->setLineageStop(0);
//    fossilRanges[2]->setLineageStart(0.7480391);
//    fossilRanges[2]->setLineageStop(0);
//    fossilRanges[3]->setLineageStart(2.891589);
//    fossilRanges[3]->setLineageStop(1.725818);
//    fossilRanges[4]->setLineageStart(1.926918);
//    fossilRanges[4]->setLineageStop(0);
//    fossilRanges[5]->setLineageStart(0.1037843);
//    fossilRanges[5]->setLineageStop(0);
//    fossilRanges[6]->setLineageStart(0.4343739);
//    fossilRanges[6]->setLineageStop(0);
//    fossilRanges[7]->setLineageStart(1.530841);
//    fossilRanges[7]->setLineageStop(0);
//    fossilRanges[8]->setLineageStart(1.061154);
//    fossilRanges[8]->setLineageStop(0.061442);
//    fossilRanges[9]->setLineageStart(3.224778);
//    fossilRanges[9]->setLineageStop(2.3289);
//    fossilRanges[10]->setLineageStart(2.819214);
//    fossilRanges[10]->setLineageStop(2.718085);
//    fossilRanges[11]->setLineageStart(2.938183);
//    fossilRanges[11]->setLineageStop(0.986128);
//    fossilRanges[12]->setLineageStart(4.1891);
//    fossilRanges[12]->setLineageStop(0);
//    fossilRanges[13]->setLineageStart(0.7107883);
//    fossilRanges[13]->setLineageStop(0);
//    fossilRanges[14]->setLineageStart(0.9757994);
//    fossilRanges[14]->setLineageStop(0);
//    fossilRanges[15]->setLineageStart(0.9359346);
//    fossilRanges[15]->setLineageStop(0);
//    fossilRanges[16]->setLineageStart(1.664171);
//    fossilRanges[16]->setLineageStop(1.471591);
//    fossilRanges[17]->setLineageStart(1.532807);
//    fossilRanges[17]->setLineageStop(0.2839254);
//    fossilRanges[18]->setLineageStart(0.3198285);
//    fossilRanges[18]->setLineageStop(0);
//    fossilRanges[19]->setLineageStart(2.431094);
//    fossilRanges[19]->setLineageStop(2.286943);
//    fossilRanges[20]->setLineageStart(4.747997);
//    fossilRanges[20]->setLineageStop(3.766295);
//    fossilRanges[21]->setLineageStart(3.781237);
//    fossilRanges[21]->setLineageStop(2.150407);
//    fossilRanges[22]->setLineageStart(2.938139);
//    fossilRanges[22]->setLineageStop(0);
//    fossilRanges[23]->setLineageStart(0.4388662);
//    fossilRanges[23]->setLineageStop(0);
//    fossilRanges[24]->setLineageStart(0.334802);
//    fossilRanges[24]->setLineageStop(0);
//    fossilRanges[25]->setLineageStart(0.07767761);
//    fossilRanges[25]->setLineageStop(0);
//    fossilRanges[26]->setLineageStart(0.4731897);
//    fossilRanges[26]->setLineageStop(0.4493463);
//    fossilRanges[27]->setLineageStart(0.5304301);
//    fossilRanges[27]->setLineageStop(0);
//    fossilRanges[28]->setLineageStart(0.3094577);
//    fossilRanges[28]->setLineageStop(0);
//    fossilRanges[29]->setLineageStart(1.331281);
//    fossilRanges[29]->setLineageStop(0);
//    fossilRanges[30]->setLineageStart(0.5376007);
//    fossilRanges[30]->setLineageStop(0);
//    fossilRanges[31]->setLineageStart(0.836111);
//    fossilRanges[31]->setLineageStop(0);
//    fossilRanges[32]->setLineageStart(0.5588915);
//    fossilRanges[32]->setLineageStop(0);
//    fossilRanges[33]->setLineageStart(0.2243971);
//    fossilRanges[33]->setLineageStop(0);
//    fossilRanges[34]->setLineageStart(0.2036094);
//    fossilRanges[34]->setLineageStop(0);
//    fossilRanges[35]->setLineageStart(0.2659481);
//    fossilRanges[35]->setLineageStop(0);
//    fossilRanges[36]->setLineageStart(0.4489092);
//    fossilRanges[36]->setLineageStop(0);
//    fossilRanges[37]->setLineageStart(2.990018);
//    fossilRanges[37]->setLineageStop(1.251231);
//    fossilRanges[38]->setLineageStart(4.067044);
//    fossilRanges[38]->setLineageStop(3.692535);
//    fossilRanges[39]->setLineageStart(3.951721);
//    fossilRanges[39]->setLineageStop(3.926339);
//    fossilRanges[40]->setLineageStart(5.560851);
//    fossilRanges[40]->setLineageStop(3.710196);
//    fossilRanges[41]->setLineageStart(5.23413);
//    fossilRanges[41]->setLineageStop(5.016654);
    
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
        ss << "\tz_f(FR_" << frgID << ")"; //rw: lineage start
        ss << "\ta_f(FR_" << frgID << ")"; //rw: last appearance
        ss << "\tb_f(FR_" << frgID << ")"; //rw: lineage stop
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
    
    nprb = numFossils*log(fossRate);
    nprb += numExtinctLineages*log(mu);
    nprb -= log(lambda * (1-fbdPFxn(lambda,mu,fossRate, sppSampRate,ot)) );
    
    for(int f=0; f < fossilRanges.size(); f++){
        
        FossilRange *fr = fossilRanges[f];
        
        double bi = fr->getLineageStart(); //rw: zf
        double di = fr->getLineageStop();  //rw: bf
        
        nprb += log( lambda * fr->getFossilRangeBrGamma() );
        
        double rangePr = 0;
        rangePr += fbdQTildaFxnLog(lambda, mu, fossRate, sppSampRate, bi);
        rangePr -= fbdQTildaFxnLog(lambda, mu, fossRate, sppSampRate, di);
        
        nprb += rangePr;
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
    //cout << "q_tilda " << v << endl;
    
    return v;
}


//END