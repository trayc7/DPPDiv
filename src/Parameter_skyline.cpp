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
#include "Parameter_skyline.h"
#include "Parameter_fossilgraph.h"
#include "Parameter_treescale.h"
#include "Parameter_tree.h"
#include "MbRandom.h"
#include "Model.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

//#include <unistd.h>//rw: remember to remove

using namespace std;

Skyline::Skyline(MbRandom *rp, Model *mp, int nrts) : Parameter(rp, mp) {
	
	
	rho = 1.0;
	numRates = nrts; // this is the number of rates, which is 1 less than the number of time points
	name = "SL";
	
	for(int i=0; i<numRates; i++){
		lambdas.push_back(0.5);
		mus.push_back(0.1);
		psis.push_back(0.05);
	}
	
}

Skyline::~Skyline(void) {
	
}

Skyline& Skyline::operator=(const Skyline &c) {
	
	if (this != &c)
		clone(c);
	return *this;
}

void Skyline::clone(const Skyline &c) {
	
	rho = c.rho;
	name = "SL";
}

void Skyline::print(std::ostream & o) const {
	
}

double Skyline::update(double &oldLnL) {
	
	double lnR = 0.0;
	
	return lnR;
}



double Skyline::lnPrior(void) {
	
	return 0.0;
}

string Skyline::writeParam(void){
	
	stringstream ss;
	ss << "Skyline parameters: m/l = " << fixed << setprecision(4) << rho
	   
	   << endl;
	return ss.str();
}

