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

#include "Parameter.h"
#include "Parameter_speciationskyline.h"
#include "Parameter_fossilrangegraphskyline.h"
#include "MbRandom.h"
#include "Model.h"
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

SpeciationSkyline::SpeciationSkyline(MbRandom *rp, Model *mp, int ni, double rh) : Parameter(rp, mp) {
    
    name = "SSLP";
    numIntervals = ni + 1;
    extantSampleRate = rh;
    initializeIntervalVariables();
    currentFossilRangeGraphSkylineLnL = 0.0;
    parameterization = 3; // hard coded at the moment
    
    //fixPsi = fxPsi;
    setAllBDFossParams();
    
    // priors on birth death paras
    int specPr = 2; // 1 = unifrom prior, 2 = exponential prior
    int psiPr = 1;
    double bPrRate = 1;
    double dPrRate = 1;
    double pPrRate = 1;
    deathRatePrior = specPr;
    birthRatePrior = specPr;
    fossRatePrior = psiPr;
    birthRateExpRate = bPrRate;
    deathRateExpRate = dPrRate; ;
    fossRateExpRate = pPrRate;
    
    printInitialIntervalVariables();
    
    cout << "BD initialized\n";
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

void SpeciationSkyline::initializeIntervalVariables(){
    
    for(int i = 0; i < numIntervals; i++){
        birthRates.push_back(0.02);
        deathRates.push_back(0.01);
        fossilRates.push_back(2); //**skyline note - add user defined flag psi
    }
    for(int i=0; i < numIntervals; i++){
        netDiversifications.push_back(birthRates[i] - deathRates[i]);
        relativeDeaths.push_back(deathRates[i] / birthRates[i]);
        probObservations.push_back(fossilRates[i] / (deathRates[i] + fossilRates[i]));
    }
    
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
        cout << ", psi = " << fossilRates[i] << endl;
        j++;
    }
}

double SpeciationSkyline::updateFossilRangeGraphSkylineBDParams(double &oldLnL){
    
    currentFossilRangeGraphSkylineLnL = oldLnL;
    FossilRangeGraphSkyline *frg = modelPtr->getActiveFossilRangeGraphSkyline();
    
//    if(fixPsi)
//        int numMoves =  numIntervals * 2;
//    else
        int numMoves =  numIntervals * 3;
    
    int I, v;
    
    for(int i=0; i < numMoves; i++){
        // choose random interval
        I = (int)(ranPtr->uniformRv(0.0, numIntervals - 1));
        // choose random parameters (never go to 3 if fixPsi = T)
        v = (int)(ranPtr->uniformRv() * numMoves);
//        if(v == 0)
//            updateDeathRate(frg, I); // mu
//        else if(v == 1)
//            updateBirthRate(frg, I); // lambda
//        else if(v == 2)
//            updatePsiRate(frg, I); // psi // shouldnt reach here if fixpsi = 1
    }
    
    return currentFossilRangeGraphSkylineLnL;
}

// other functions


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
