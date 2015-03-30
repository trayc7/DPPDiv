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

using namespace std;

FossilGraph::FossilGraph(MbRandom *rp, Model *mp, int nf, double initOrigTime, vector<Calibration *> clb) : Parameter(rp, mp){

    name = "FG";
    numFossils = nf;
    originTime = initOrigTime;
    createOccurrenceVector(clb);
    initializeOccurrenceSpecVariables();
}

FossilGraph::~FossilGraph(void){
    
}

FossilGraph& FossilGraph::operator=(const FossilGraph &t) {
    
    if (this != &t)
        clone(t);
    return *this;
}

void FossilGraph::clone(const FossilGraph &t){

    numFossils = t.numFossils;
}

double FossilGraph::update(double &oldLnL){

    return 0;
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
    for(int c = 0; c < clb.size(); c++){
        Calibration *p = clb[c];
        double yf = p->getYngTime();
        Occurrence *o = new Occurrence(yf);
        occurrenceSpecimens.push_back(o);
        if( yf < et)
            et = yf;
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
    
    for(int f = 0; f < numFossils; f++){
        Occurrence *o = occurrenceSpecimens[f];
        if(o->getIsTerminal()){
            o->setFossilSppTime(o->getFossilAge());
            o->setFossilIndicatorVar(0);
        }
        else{
            double u = ranPtr->uniformRv();
            if(u < 0.5){
                o->setFossilSppTime(o->getFossilAge());
                o->setFossilIndicatorVar(0);
            }
            else {
                double yf = o->getFossilAge();
                double zf = ranPtr->uniformRv(yf,originTime);
                o->setFossilSppTime(zf);//rw: this was the problem, setFossilFossBrGamma used instead of setFossilSppTime
                o->setFossilIndicatorVar(1);
                o->setFossilFossBrGamma(0);
            }
        }
    }
    recountOccurrenceAttachNums();
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
        //cout << "gamma :" << g << endl; //rw: gamma is correct
        o->setFossilFossBrGamma(g);
    }
}

double FossilGraph::getActiveFossilGraphProb(){
    
    double nprb = 0.0;
    Speciation *s = modelPtr->getActiveSpeciation();
    //OriginTime *ot = modelPtr->getActiveOriginTime(); I don't think this is needed here
    s->setAllBDFossParams();
    double lambda = s->getBDSSSpeciationRateLambda();
    double mu = s->getBDSSExtinctionRateMu();
    double fossRate = s->getBDSSFossilSampRatePsi();
    double sppSampRate = s->getBDSSSppSampRateRho();
    recountOccurrenceAttachNums();
    nprb = getFossilGraphProb(lambda, mu, fossRate, sppSampRate);
    return nprb;
}


double FossilGraph::getFossilGraphProb(double lambda, double mu, double fossRate, double sppSampRate) {
    OriginTime *ot = modelPtr->getActiveOriginTime();
    originTime = ot->getOriginTime(); // active origin?
    
    // the following has been modified for the fofbd
    double nprb = 1.0 - (log(2*lambda) + log(1.0 - bdssP0Fxn(lambda, mu, fossRate, sppSampRate, originTime)));
    //cout << "1. nprb = " << nprb << endl;
    for(int f=0; f < occurrenceSpecimens.size(); f++){
        Occurrence *o = occurrenceSpecimens[f];
        nprb += log(fossRate * o->getFossilFossBrGamma() );
        //cout << "Fossil rate " << fossRate << endl;
        //cout << "Gamma " << o->getFossilFossBrGamma() << endl;
        //cout << "2. nprb = " << nprb << endl;
        if(o->getFossilIndicatorVar()){
            double fossAge = o->getFossilAge();
            double fossPhi = o->getFossilSppTime(); // n.b. treescale removed from this part of the equation
            double fossPr = log(2.0 * lambda) + log( bdssP0Fxn(lambda, mu, fossRate, sppSampRate, fossAge) );
            fossPr += fbdQHatFxn(lambda, mu, fossRate, sppSampRate, fossPhi);
            fossPr -= fbdQHatFxn(lambda, mu, fossRate, sppSampRate, fossAge);
            nprb += fossPr;
        }
    }
    //cout << "3. FBDS prob = " << nprb << endl;
    return nprb;
}

double FossilGraph::bdssP0Fxn(double b, double d, double psi, double rho, double t){
    
    double c1Val = bdssC1Fxn(b,d,psi);
    double c2Val = bdssC2Fxn(b,d,psi,rho);
    
    double eCfrac = (exp(-c1Val * t) * (1.0 - c2Val) - (1.0 + c2Val)) / (exp(-c1Val * t) * (1.0 - c2Val) + (1.0 + c2Val));
    double v = 1.0 + ((-(b - d - psi)) + (c1Val * eCfrac)) / (2.0 * b);
    
    return v;
}

double FossilGraph::fbdQHatFxn(double b, double d, double psi, double rho, double t){
    
    double v = log(4.0 * rho);
    v -= bdssQFxn(b,d,psi,rho,t);
    return v;
}

double FossilGraph::bdssQFxn(double b, double d, double psi, double rho, double t){
    
    double c1Val = bdssC1Fxn(b,d,psi);
    double c2Val = bdssC2Fxn(b,d,psi,rho);
    
    double vX = c1Val * t + 2.0 * log(exp(-c1Val * t) * (1.0 - c2Val) + (1.0 + c2Val));
    
    // returns log(q)
    return vX;
}

double FossilGraph::bdssC1Fxn(double b, double d, double psi){
    
    double v = abs( sqrt( ( (b-d-psi) * (b-d-psi) ) + 4*b*psi) );
    return v;
}

double FossilGraph::bdssC2Fxn(double b, double d, double psi, double rho){
    
    double v = -( ( b-d-(2*b*rho)-psi ) / (bdssC1Fxn(b,d,psi)) );
    return v;
}


// the following functions are to do with update moves

double FossilGraph::updateOccurrenceAttachmentTimesPhi() {
    
    OriginTime *ot = modelPtr->getActiveOriginTime();
    originTime = ot->getOriginTime(); // active origin?
    
    Speciation *s = modelPtr->getActiveSpeciation();
    s->setAllBDFossParams();
    double lambda = s->getBDSSSpeciationRateLambda();
    double mu = s->getBDSSExtinctionRateMu();
    double fossRate = s->getBDSSFossilSampRatePsi();
    double sppSampRate = s->getBDSSSppSampRateRho();
    
    vector<int> rndOccIDs;
    for(int i=0; i<occurrenceSpecimens.size(); i++)
        rndOccIDs.push_back(i);
    random_shuffle(rndOccIDs.begin(), rndOccIDs.end());
    for(vector<int>::iterator it=rndOccIDs.begin(); it!=rndOccIDs.end(); it++){
        Occurrence *o = occurrenceSpecimens[(*it)];
        if(o->getFossilIndicatorVar()){
            
            /*
            Node *p = &nodes[f->getFossilMRCANodeID()];
            double nodeDepth = p->getNodeDepth() * treeScale;
            if(f->getIsTotalGroupFossil()){
                if(p != root){
                    nodeDepth = p->getAnc()->getNodeDepth()*treeScale;
                }
                else {
                    OriginTime *ot = modelPtr->getActiveOriginTime();
                    nodeDepth = ot->getOriginTime();
                }
            }
            */
            
            double fossDepth = o->getFossilAge();
            
            double oldPhi = o->getFossilSppTime(); // * treeScale;
            double oldSumLogGammas = getSumLogAllAttachNums(); //??
            
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
            
            double newSumLogGammas = getSumLogAllAttachNums();
            
            double lnPriorRat = 0.0;
            double v1 = newSumLogGammas;
            double v2 = bdssQFxn(lambda, mu, fossRate, sppSampRate, newPhi);
            double v3 = bdssQFxn(lambda, mu, fossRate, sppSampRate, oldPhi);
            double v4 = oldSumLogGammas;
            lnPriorRat = (v1 - v2) - (v4 - v3);
            
            
            double r = modelPtr->safeExponentiation(lnPriorRat + c);
            
            if(ranPtr->uniformRv() < r){ 
                o->setFossilSppTime(newPhi);
                //setNodeOldestAttchBranchTime(p);
            }
            else{
                o->setFossilSppTime(oldPhi);
                //setNodeOldestAttchBranchTime(p);
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
