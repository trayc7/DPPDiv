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
#include "Parameter_speciationskyline.h"
#include "Parameter_fossilrangegraphskyline.h"
#include "MbRandom.h"
#include "Model.h"
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

SpeciationSkyline::SpeciationSkyline(MbRandom *rp, Model *mp, int ni, double rh, int specPr, int psiPr, double bPrRate, double dPrRate, double pPrRate, bool proxy, vector<Calibration *> ints, bool fxPsi) : Parameter(rp, mp) {
    
    name = "SSLP";
    numIntervals = ni + 1;
    extantSampleRate = rh;
    useSamplingProxy = proxy;
    proxyInt = 0;
    currentFossilRangeGraphSkylineLnL = 0.0;
    parameterization = 3; // hard coded for the moment
    maxdivV = 30000.0;
    
    fixPsi = fxPsi;
    constantRateModel = 1;
    fixAllPsi = false;
    if(fixAllPsi)
        numParameters = 2;
    else
        numParameters = 3;
    
    numMoves = numIntervals * numParameters;
    
    initializeIntervalVariables(ints);
    setAllBDFossParams();
    
    // priors on birth death paras
    deathRatePrior = specPr; // 1 = unifrom prior, 2 = exponential prior
    birthRatePrior = specPr;
    fossRatePrior = psiPr;
    birthRateExpRate = bPrRate;
    deathRateExpRate = dPrRate;
    fossRateExpRate = pPrRate;
    
    printInitialIntervalVariables();
    
    cout << "BD initialized\n" << endl;
}

SpeciationSkyline::~SpeciationSkyline(void) {
    
}

SpeciationSkyline& SpeciationSkyline::operator=(const SpeciationSkyline &c) {
    
    if (this != &c)
        clone(c);
    return *this;
    
}

void SpeciationSkyline::clone(const SpeciationSkyline &c) {
    
    //setAllBDFossParams();
    //relativeDeath = c.relativeDeath;
    //netDiversificaton = c.netDiversificaton;
    //probSpeciationS = c.probSpeciationS;
    
    numIntervals = c.numIntervals;
    
    extantSampleRate = c.extantSampleRate;
    for(int i = 0; i < numIntervals; i++){
        birthRates[i] = c.birthRates[i];
        deathRates[i] = c.deathRates[i];
        fossilRates[i] = c.fossilRates[i];
    }
    
    name = "SSLP";
}

void SpeciationSkyline::initializeIntervalVariables(vector<Calibration *> ints){
    
    for(int i = 0; i < numIntervals; i++){
        birthRates.push_back(1.0);
        deathRates.push_back(0.5);
        fossilRates.push_back(0.7);
        tipRates.push_back(0.0);
    }
    tipRates[0] = extantSampleRate;
    
    for(int i=0; i < numIntervals; i++){
        netDiversifications.push_back(birthRates[i] - deathRates[i]);
        relativeDeaths.push_back(deathRates[i] / birthRates[i]);
        probObservations.push_back(fossilRates[i] / (deathRates[i] + fossilRates[i]));
    }
    
    // this is experimental
    if(useSamplingProxy){
        
        double proxy = 0.0;
        double proxyMax = 0.0;
        vector<double> proxyData;
        
        for(int i=0; i < ints.size(); i++){
            Calibration *h = ints[i];
            proxy = h->getIntervalProxy();
            proxyData.push_back(proxy);
            cout << "proxy " << proxy << endl;
            if(proxy > proxyMax){
                proxyInt = i;
                proxyMax = proxy;
            }
        }

        // oldest interval
        // proxyData.push_back(1.0); // this won't actually be used
        
        // define scale values and redefine psi
        for(int i=0; i < proxyData.size(); i++){
            proxyScale.push_back(proxyData[i]/proxyMax);
            fossilRates[i] = fossilRates[proxyInt] * proxyScale[i];
        }
    }
    if(fixPsi) fossilRates[0] = 0;
    if(fixPsi) fossilRates[numIntervals-1] = 0;
}

void SpeciationSkyline::print(std::ostream & o) const {
    
}

double SpeciationSkyline::update(double &oldLnL) {
    
    return updateFossilRangeGraphSkylineBDParams(oldLnL);
    
}

double SpeciationSkyline::lnPrior(void) {
    
    return 0.0;
}

string SpeciationSkyline::writeParam(void){
    
    string s="";
    return s;
    
}

void SpeciationSkyline::printInitialIntervalVariables(){
    
    cout << "Speciaton parameters are initialized with: " << endl;
    int j = 1;
    for(int i = 0; i < numIntervals; i++){
        cout << "Interval " << j << ": ";
        cout << "lambda = " << birthRates[i];
        cout << ", mu = " << deathRates[i];
        if(useSamplingProxy)
            cout << ", psi scale = " << proxyScale[i];
        cout << ", psi = " << fossilRates[i];
        cout << ", rho = " << tipRates[i] << endl;
        j++;
    }
}

double SpeciationSkyline::updateFossilRangeGraphSkylineBDParams(double &oldLnL){
    
    currentFossilRangeGraphSkylineLnL = oldLnL;
    FossilRangeGraphSkyline *frgs = modelPtr->getActiveFossilRangeGraphSkyline();
    
    int t, v;
    
    if(useSamplingProxy){
        for(int i=0; i < numMoves; i++){
            // choose random interval
            t = (int)(ranPtr->uniformRv(0.0, numIntervals));
            // choose random parameters
            v = (int)(ranPtr->uniformRv() * numParameters);
            //cout << "t " << t << endl;
            if(v == 0)
                updateDeathRate(frgs, t); // mu
            else if(v == 1)
                updateBirthRate(frgs, t); // lambda
            else if(v == 2){
                if(t == numIntervals - 1)
                    updatePsiRate(frgs, t);
                else
                    updateProxyRate(frgs); // psi
            }
        }
    } else if(constantRateModel){
        for(int i=0; i < numMoves; i++){
            // choose random parameters
            v = (int)(ranPtr->uniformRv() * numMoves);
            if(v == 0)
                updateBirthOneRate(frgs); // lambda
            else if(v == 1)
                updateDeathOneRate(frgs); // mu
            else if(v == 2)
                updatePsiOneRate(frgs); // psi
        }
    } else {
        for(int i=0; i < numMoves; i++){
            // choose random interval
            t = (int)(ranPtr->uniformRv(0.0, numIntervals));
            // choose random parameters (never go to 3 if fixAllPsi = T)
            v = (int)(ranPtr->uniformRv() * numParameters);
            if(v == 0)
                updateDeathRate(frgs, t); // mu
            else if(v == 1)
                updateBirthRate(frgs, t); // lambda
            else if(v == 2)
                updatePsiRate(frgs, t); // psi
        }
    }
    
    return currentFossilRangeGraphSkylineLnL;
}

double SpeciationSkyline::updateBirthRate(FossilRangeGraphSkyline *frgs, int i) {
    
    double oldfgprob = currentFossilRangeGraphSkylineLnL;
    double lpr = 0.0;
    double oldLambda = birthRates[i];
    double newLambda;
    double tuning = log(3.0);
    double minV = 0.0001;
    double c = getNewValScaleMv(newLambda, oldLambda, minV, maxdivV, tuning);
    birthRates[i] = newLambda;
    lpr = c;
    double newfgprob = getLnFossilRangeGraphSkylineProb(frgs);
    double lnPriorRat = getExpPriorRatio(oldLambda, newLambda, birthRateExpRate, birthRatePrior);
    double lnLikeRat = (newfgprob - oldfgprob);
    double lnR = lnLikeRat + lpr + lnPriorRat;
    double r = modelPtr->safeExponentiation(lnR);
    
    if(ranPtr->uniformRv() < r){
        currentFossilRangeGraphSkylineLnL = newfgprob;
    }
    else{
        birthRates[i] = oldLambda;
        setAllBDFossParams();
        frgs->setAllIntervalConstants();
        currentFossilRangeGraphSkylineLnL = oldfgprob;
    }
    return 0.0;
}

double SpeciationSkyline::updateDeathRate(FossilRangeGraphSkyline *frgs, int i) {
    
    double oldfgprob = currentFossilRangeGraphSkylineLnL;
    double lpr = 0.0;
    double oldMu = deathRates[i];
    double newMu;
    double tuning = log(2.0);
    double minV = 0.0001;
    double c = getNewValScaleMv(newMu, oldMu, minV, maxdivV, tuning);
    deathRates[i] = newMu;
    lpr = c;
    double newfgprob = getLnFossilRangeGraphSkylineProb(frgs);
    double lnPriorRat = getExpPriorRatio(oldMu, newMu, deathRateExpRate, deathRatePrior);
    double lnLikeRat = (newfgprob - oldfgprob);
    double lnR = lnLikeRat + lpr + lnPriorRat;
    double r = modelPtr->safeExponentiation(lnR);
    
    if(ranPtr->uniformRv() < r){
        currentFossilRangeGraphSkylineLnL = newfgprob;
    }
    else{
        deathRates[i] = oldMu;
        setAllBDFossParams();
        frgs->setAllIntervalConstants();
        currentFossilRangeGraphSkylineLnL = oldfgprob;
    }
    return 0.0;
}

double SpeciationSkyline::updatePsiRate(FossilRangeGraphSkyline *frgs, int i) {
    
    if( (fixPsi && i == 0) | (fixPsi && i == numIntervals - 1) )
        return 0.0;
    
    double oldfgprob = currentFossilRangeGraphSkylineLnL;
    double lpr = 0.0;
    double oldPsi = fossilRates[i];
    double newPsi;
    double tuning = log(2.0);
    double minV = 0.0001;
    double maxV = 100000;
    double c = getNewValScaleMv(newPsi, oldPsi, minV, maxV, tuning);
    fossilRates[i] = newPsi;
    lpr = c;
    double newfgprob = getLnFossilRangeGraphSkylineProb(frgs);
    double lnPriorRat = getExpPriorRatio(oldPsi, newPsi, fossRateExpRate, fossRatePrior);
    double lnLikeRat = (newfgprob - oldfgprob);
    double lnR = lnLikeRat + lpr + lnPriorRat;
    double r = modelPtr->safeExponentiation(lnR);
    
    if(ranPtr->uniformRv() < r){
        currentFossilRangeGraphSkylineLnL = newfgprob;
    }
    else{
        fossilRates[i] = oldPsi;
        setAllBDFossParams();
        frgs->setAllIntervalConstants();
        currentFossilRangeGraphSkylineLnL = oldfgprob;
    }
    return 0.0;
}

double SpeciationSkyline::updateProxyRate(FossilRangeGraphSkyline *frgs){
    
    // this is experimental
    
    double oldfgprob = currentFossilRangeGraphSkylineLnL;
    vector<double> oldFossilRates = fossilRates;
    
    double lpr = 0.0;
    double oldPsi = fossilRates[proxyInt];
    double newPsi;
    double tuning = log(2.0);
    double minV = 0.0001;
    double maxV = 100000;
    double c = getNewValScaleMv(newPsi, oldPsi, minV, maxV, tuning);
    fossilRates[proxyInt] = newPsi;
    
    // scale psi values using proxy data
    for(int i=0; i < proxyScale.size(); i++){
        fossilRates[i] = fossilRates[proxyInt] * proxyScale[i];
    }
    
    lpr = c;
    double newfgprob = getLnFossilRangeGraphSkylineProb(frgs);
    double lnPriorRat = getExpPriorRatio(oldPsi, newPsi, fossRateExpRate, fossRatePrior);
    double lnLikeRat = (newfgprob - oldfgprob);
    double lnR = lnLikeRat + lpr + lnPriorRat;
    double r = modelPtr->safeExponentiation(lnR);
    
    if(ranPtr->uniformRv() < r){
        currentFossilRangeGraphSkylineLnL = newfgprob;
    } else {
        fossilRates = oldFossilRates;
        setAllBDFossParams();
        frgs->setAllIntervalConstants();
        currentFossilRangeGraphSkylineLnL = oldfgprob;
    }
    
    return 0.0;
}

double SpeciationSkyline::updateBirthOneRate(FossilRangeGraphSkyline *frgs) {
    
    double oldfgprob = currentFossilRangeGraphSkylineLnL;
    double lpr = 0.0;
    double oldLambda = birthRates[0];
    double newLambda;
    double tuning = log(3.0);
    double minV = 0.0001;
    double c = getNewValScaleMv(newLambda, oldLambda, minV, maxdivV, tuning);
    for(int i=0; i < numIntervals; i++){
        birthRates[i] = newLambda;
    }
    lpr = c;
    double newfgprob = getLnFossilRangeGraphSkylineProb(frgs);
    double lnPriorRat = getExpPriorRatio(oldLambda, newLambda, birthRateExpRate, birthRatePrior);
    double lnLikeRat = (newfgprob - oldfgprob);
    double lnR = lnLikeRat + lpr + lnPriorRat;
    double r = modelPtr->safeExponentiation(lnR);
    
    if(ranPtr->uniformRv() < r){
        currentFossilRangeGraphSkylineLnL = newfgprob;
    }
    else{
        for(int i=0; i < numIntervals; i++){
            birthRates[i] = oldLambda;
        }
        setAllBDFossParams();
        frgs->setAllIntervalConstants();
        currentFossilRangeGraphSkylineLnL = oldfgprob;
    }
    return 0.0;
}

double SpeciationSkyline::updateDeathOneRate(FossilRangeGraphSkyline *frgs) {
    
    double oldfgprob = currentFossilRangeGraphSkylineLnL;
    double lpr = 0.0;
    double oldMu = deathRates[0];
    double newMu;
    double tuning = log(2.0);
    double minV = 0.0001;
    double c = getNewValScaleMv(newMu, oldMu, minV, maxdivV, tuning);
    for(int i=0; i < numIntervals; i++){
        deathRates[i] = newMu;
    }
    lpr = c;
    double newfgprob = getLnFossilRangeGraphSkylineProb(frgs);
    double lnPriorRat = getExpPriorRatio(oldMu, newMu, deathRateExpRate, deathRatePrior);
    double lnLikeRat = (newfgprob - oldfgprob);
    double lnR = lnLikeRat + lpr + lnPriorRat;
    double r = modelPtr->safeExponentiation(lnR);
    
    if(ranPtr->uniformRv() < r){
        currentFossilRangeGraphSkylineLnL = newfgprob;
    }
    else{
        for(int i=0; i < numIntervals; i++){
            deathRates[i] = oldMu;
        }
        setAllBDFossParams();
        frgs->setAllIntervalConstants();
        currentFossilRangeGraphSkylineLnL = oldfgprob;
    }
    return 0.0;
}

double SpeciationSkyline::updatePsiOneRate(FossilRangeGraphSkyline *frgs) {
    
    double oldfgprob = currentFossilRangeGraphSkylineLnL;
    double lpr = 0.0;
    double oldPsi = fossilRates[0];
    double newPsi;
    double tuning = log(2.0);
    double minV = 0.0001;
    double maxV = 100000;
    double c = getNewValScaleMv(newPsi, oldPsi, minV, maxV, tuning);
    for(int i=0; i < numIntervals; i++){
        fossilRates[i] = newPsi;
    }
    lpr = c;
    double newfgprob = getLnFossilRangeGraphSkylineProb(frgs);
    double lnPriorRat = getExpPriorRatio(oldPsi, newPsi, fossRateExpRate, fossRatePrior);
    double lnLikeRat = (newfgprob - oldfgprob);
    double lnR = lnLikeRat + lpr + lnPriorRat;
    double r = modelPtr->safeExponentiation(lnR);
    
    if(ranPtr->uniformRv() < r){
        currentFossilRangeGraphSkylineLnL = newfgprob;
    }
    else{
        for(int i=0; i < numIntervals; i++){
            fossilRates[i] = oldPsi;
        }
        setAllBDFossParams();
        frgs->setAllIntervalConstants();
        currentFossilRangeGraphSkylineLnL = oldfgprob;
    }
    return 0.0;
}

double SpeciationSkyline::getNewValScaleMv(double &nv, double ov, double vmin, double vmax, double tv){
    
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

double SpeciationSkyline::getExpPriorRatio(double oldVal, double newVal, double rate, double prior) {
    if (prior == 2)
        return (ranPtr->lnExponentialPdf(rate, newVal)) - (ranPtr->lnExponentialPdf(rate, oldVal));
    else
        return 0.0;
}

double SpeciationSkyline::getLnFossilRangeGraphSkylineProb(FossilRangeGraphSkyline *frgs) {

    setAllBDFossParams();
    frgs->setAllIntervalConstants();    
    double fgprob = frgs->getFossilRangeGraphSkylineProb();
    return fgprob;
}

// check the overall time this function takes
void SpeciationSkyline::setAllBDFossParams(){
    
    if(parameterization == 3){
        for(int i=0; i < numIntervals; i++){
            netDiversifications[i] = birthRates[i] - deathRates[i];
            relativeDeaths[i] = deathRates[i] / birthRates[i];
            probObservations[i] = fossilRates[i] / (deathRates[i] + fossilRates[i]);
        }
    }
    else {
        cerr << "Parameterization 1, 2 & 4 not yet available in skyline mode " << endl;
        exit(1);
    }
}

//END
