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
        birthRates.push_back(1);
        deathRates.push_back(0.1);
        fossilRates.push_back(2);
    }
    
}

void SpeciationSkyline::print(std::ostream & o) const {
    
}

double SpeciationSkyline::update(double &oldLnL) {
    
    double lnR = 0.0;
    
    //return updateFossilRangeGraphBDParams(oldLnL);
    
    return lnR;
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

//END
