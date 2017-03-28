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

FossilRangeGraphSkyline::FossilRangeGraphSkyline(MbRandom *rp, Model *mp, int nf, int nl, vector<Calibration *> clb, int ni) : Parameter(rp, mp){

    //**skynote most of this will be the same
    name = "FRGS";
    numFossils = nf;
    numLineages = nl;
    numIntervals = ni;
    
}

FossilRangeGraphSkyline::~FossilRangeGraphSkyline(void){
    
}

FossilRangeGraphSkyline& FossilRangeGraphSkyline::operator=(const FossilRangeGraphSkyline &frgsl) {
    
    if (this != &frgsl)
        clone(frgsl);
    return *this;
}

void FossilRangeGraphSkyline::clone(const FossilRangeGraphSkyline &frgsl){

    //**skyline note most of this will be the same
    numFossils = frgsl.numFossils;
    numLineages = frgsl.numLineages;
    
}

double FossilRangeGraphSkyline::update(double &oldLnL){
    
    currentFossilRangeGraphSkylineLnL = oldLnL;
    
    // this should be the same
    
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


// create fossil range vector
// and intialize fossil range variables
// print range variables
// most of this will be the same, with the exception that the ranges need to know what interval
// their variables correspond to

// recount gamma
// don't think this should be different

// put functions for assigning ranges to intervals here

// count extinct ranges
// don't think this should be different

// redefine origin
// don't think this should be different

// update functions
// most of this will be the same, with the exception that the ranges need to know what interval
// their variables correspond to

// moves
// don't think these should be different

// ignore: getFossilRangeInfoParamList, getFossilRangeInfoParamNames

// put probability functions here

double FossilRangeGraphSkyline::getFossilRangeGraphProbSkyline(std::vector<double> lambda, double mu){
    
    
    return 0.0;
}

//END
