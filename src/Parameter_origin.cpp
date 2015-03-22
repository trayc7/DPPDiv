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

OriginTime::OriginTime(MbRandom *rp, Model *mp, double sv, double yb, double ob) : Parameter(rp, mp) {
    oldBound = ob;
    yngBound = yb;
    originTime = sv;
    tuning = ((yngBound + (sv * 1.2)) * 0.5) * 0.2;
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
    
    Treescale *ts = modelPtr->getActiveTreeScale();
    
    Speciation *s = modelPtr->getActiveSpeciation();
    
    Tree *t = modelPtr->getActiveTree();
    
    double oldOT = originTime;
//    cout << "OT1 " << originTime << endl;
    
    double oldOTProb = getFBDProbOriginTime(t, s);
    
    double minAge = t->getOldestTreeSpeciationTime();
    double maxAge = oldBound;
    
    double limO = oldOT + tuning;
    double limY = oldOT - tuning;
    
    
    double newOT;
    
    double u = ranPtr->uniformRv(-0.5,0.5) * (limO - limY);
    newOT = oldOT + u;
    while(newOT < minAge || newOT > maxAge){
        if(newOT < minAge)
            newOT = (2 * minAge) - newOT;
        if(newOT > maxAge)
            newOT = (2 * maxAge) - newOT;
    }
    originTime = newOT;
//    cout << "OT2 " << originTime << endl;
    double newOTProb = getFBDProbOriginTime(t, s);
    
//    cout << "old: OT = " << oldOT << ", old Prob = " << oldOTProb << endl;
//    cout << "new: OT = " << newOT << ", new Prob = " << newOTProb << endl;
    
    // Need to double check this

    modelPtr->setLnLGood(true);
    modelPtr->setMyCurrLnl(oldLnL);
    double myPrR = oldOTProb-newOTProb;
    
    return myPrR;
    
}

//double OriginTime::doAScaleMove(double currentV, double lb, double hb, double rv){
//    
//    double tuningV =
//    double c = tv * (rv - 0.5);
//    double newcv = cv * exp(c);
//    bool validV = false;
//    do{
//        if(newcv < lb)
//            newcv = lb * lb / newcv;
//        else if(newcv > hb)
//            newcv = hb * hb / newcv;
//        else
//            validV = true;
//    } while(!validV);
//    nv = newcv;
//    return c;
//}


double OriginTime::getFBDProbOriginTime(Tree *t, Speciation *s){
    
    // this returns log[ 1 / ((lambda*(1-Phat(x_0))) ]
    s->setAllBDFossParams();
    
    double lambda = s->getBDSSSpeciationRateLambda();
    double mu = s->getBDSSExtinctionRateMu();
    double rho = s->getBDSSSppSampRateRho();
    
    double phat = t->bdssP0HatFxn(lambda, mu, rho, originTime);

    double otPrb = log(lambda) + log(1-phat);
    return otPrb;
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






