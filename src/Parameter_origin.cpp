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
#include "Parameter_expcalib.h"
#include "Parameter_speciaton.h"
#include "Parameter_treescale.h"
#include "Parameter_origin.h"
#include "Parameter_tree.h"
#include "MbMath.h"
#include "MbRandom.h"
#include "Model.h"
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

OriginTime::OriginTime(MbRandom *rp, Model *mp, double sv, double yb, double ob, int dt, bool calib) : Parameter(rp, mp) {
    oldBound = ob;
    yngBound = yb;
    originTime = 4000.0;
    name = "OT";
}

OriginTime::~OriginTime(void) {
    
}

OriginTime& OriginTime::operator=(const OriginTime &c) {
    
    if (this != &c)
        clone(c);
    return *this;
}

void OriginTime::clone(const OriginTime &c) {
    
    originTime = c.originTime;
    oldBound = c.oldBound;
    yngBound = c.yngBound;
    isBounded = c.isBounded;
    tuning = c.tuning;
    name = "OT";
}

void OriginTime::print(std::ostream & o) const {
    
    o << "Origin time parameter: ";
    o << fixed << setprecision(4) << originTime << " ";
    o << endl;
}

double OriginTime::update(double &oldLnL) {
    
    double lppr = 0.0;
    
    return lppr;
}

double OriginTime::lnPrior(void) {
    
    return 0.0;
}


string OriginTime::writeParam(void){
    
    stringstream ss;
    ss << "Origin time parameter: " << originTime << " [+/- " << tuning << "]" << endl;
    return ss.str();
}

double OriginTime::getLnTreeProb(Tree *t) {
    
    return 0.0;
}






