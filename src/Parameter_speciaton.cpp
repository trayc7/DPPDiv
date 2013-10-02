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
#include "Parameter_tree.h"
#include "MbRandom.h"
#include "Model.h"
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

Speciation::Speciation(MbRandom *rp, Model *mp, double bdr, double bda, double initRH) : Parameter(rp, mp) {
	
	
	maxdivV = 30000.0;
	name = "SP";
	relativeDeath = bda;
	netDiversificaton = bdr;
	probSpeciationS = ranPtr->uniformRv();
	extantSampleRate = 1.0; 
	setAllBDFossParams();	
	treeTimePrior = modelPtr->getTreeTimePriorNum(); 
	
	
	if(relativeDeath < 0.0 && netDiversificaton < 0.0){
		relativeDeath = ranPtr->uniformRv();
		netDiversificaton = ranPtr->uniformRv() * maxdivV;
	}
	if(relativeDeath >= 1.0){
		cerr << "ERROR: The relative death rate (-bda) is >= 1 " << endl;
		exit(1);
	}
	if(netDiversificaton <= 0.0){
		cerr << "ERROR: The net diversification rate (-bdr) is <= 0 " << endl;
		exit(1);
	}
	if(treeTimePrior == 1){
		relativeDeath = 0.0;
		netDiversificaton = 0.0;
	}	
	else if(treeTimePrior == 2)
		relativeDeath = 0.0;
	if(netDiversificaton >= maxdivV)
		netDiversificaton = maxdivV * maxdivV / netDiversificaton;
		
	if(treeTimePrior > 5){
		cout << "Speciaton parameters are fixed to: d = " << netDiversificaton << " , r = " << relativeDeath<< " , s = " << probSpeciationS << endl;
		cout << "                                   l = " << birthRate << " , m = " << deathRate<< " , psi = " << fossilRate << endl;
	}
	
}

Speciation::~Speciation(void) {
	
}

Speciation& Speciation::operator=(const Speciation &c) {
	
	if (this != &c)
		clone(c);
	return *this;
}

void Speciation::clone(const Speciation &c) {
	
	setAllBDFossParams();
	relativeDeath = c.relativeDeath;
	netDiversificaton = c.netDiversificaton;
	probSpeciationS = c.probSpeciationS;
	birthRate = c.birthRate;
	deathRate = c.deathRate;
	fossilRate = c.fossilRate;
	extantSampleRate = c.extantSampleRate;
	name = "SP";
}

void Speciation::print(std::ostream & o) const {
	if(treeTimePrior > 1){
		o << "Speciaton parameters: d/b = ";
		o << fixed << setprecision(4) << relativeDeath << " , b-d = ";
		o << fixed << setprecision(4) << netDiversificaton << " ";
		o << fixed << setprecision(4) << probSpeciationS << " ";
		o << fixed << setprecision(4) << fossilRate << " ";
		o << fixed << setprecision(4) << birthRate << " ";
		o << fixed << setprecision(4) << deathRate << " ";
		o << endl;
	}
}

double Speciation::update(double &oldLnL) {
	
	Tree *t = modelPtr->getActiveTree();
	double oldtreeprob = getLnTreeProb(t); 
	double lnProposalRatio = 0.0;	
	
	if(treeTimePrior == 2){
		lnProposalRatio += updateNetDivRate();
		relativeDeath = 0.0;
	}
	else if(treeTimePrior == 3){
		if(ranPtr->uniformRv() < 0.5)
			lnProposalRatio += updateRelDeathRt();
		else
			lnProposalRatio += updateNetDivRate();
	}
	else if(treeTimePrior > 3){ 
		double *probMove = new double[4];
		probMove[0] = 0.3;
		probMove[1] = 0.3;
		probMove[2] = 0.4;
		probMove[3] = 0.0;
		int mvType = ranPtr->categoricalRv(probMove,4);
		if(mvType == 0)
			lnProposalRatio += updateRelDeathRt();
		else if(mvType == 1)
			lnProposalRatio += updateNetDivRate();
		else if(mvType == 2)
			lnProposalRatio += updateBDSSFossilProbS();
		else if(mvType == 3)
			lnProposalRatio += updateBDSSSampleProbRho();
		delete [] probMove;
		setAllBDFossParams();
	}
	double newtreeprob = getLnTreeProb(t); 
	double lnPriorRatio = (newtreeprob - oldtreeprob);
	double lnR = lnPriorRatio + lnProposalRatio;
	
	modelPtr->setLnLGood(true);
	modelPtr->setMyCurrLnl(oldLnL);
	return lnR;
}

double Speciation::updateRelDeathRt(void) {
	
	double rdwindow = 0.2;
	double oldRD = relativeDeath;
	double u;
	double newRD;
	u = ranPtr->uniformRv(-0.5,0.5) * (rdwindow);
	newRD = oldRD + u;
	bool validV = false;
	do{
		if(newRD < 0.0)
			newRD = 0.0 - newRD;
		else if(newRD > 0.9999)
			newRD = (2 * 0.9999) - newRD;
		else
			validV = true;
	}while(!validV);
	relativeDeath = newRD;
	return 0.0;
}

double Speciation::updateNetDivRate(void) {
	
	double lpr = 0.0;
	double oldND = netDiversificaton;
	double newND;	
	double tuning = log(2.0);
	double rv = ranPtr->uniformRv();
	double c = tuning * (rv - 0.5);
	newND = oldND * exp(c);
	double minV = 0.0001;
	double maxV = maxdivV;
	bool validV = false;
	do{
		if(newND < minV)
			newND = minV * minV / newND;
		else if(newND > maxV)
			newND = maxV * maxV / newND;
		else
			validV = true;
	} while(!validV);
	netDiversificaton = newND;
	lpr = c; 
	return lpr;
}

double Speciation::updateBDSSSampleProbRho(void) {
	
	double rdwindow = 0.2;
	double oldSP = extantSampleRate;
	double u;
	double newSP;
	u = ranPtr->uniformRv(-0.5,0.5) * (rdwindow);
	newSP = oldSP + u;
	bool validV = false;
	do{
		if(newSP < 0.0)
			newSP = 0.0 - newSP;
		else if(newSP > 0.9999)
			newSP = (2 * 0.9999) - newSP;
		else
			validV = true;
	}while(!validV);
	extantSampleRate = newSP;
	return 0.0;
}

double Speciation::updateBDSSFossilProbS(void) {
	
	double swindow = 0.2;
	double oldS = probSpeciationS;
	double u;
	double newS;
	u = ranPtr->uniformRv(-0.5,0.5) * (swindow);
	newS = oldS + u;
	bool validV = false;
	do{
		if(newS < 0.0)
			newS = 0.0 - newS;
		else if(newS > 0.9999)
			newS = (2 * 0.9999) - newS;
		else
			validV = true;
	}while(!validV);
	probSpeciationS = newS;
	return 0.0;
}



double Speciation::lnPrior(void) {
	
	return 0.0;
}

string Speciation::writeParam(void){
	
	stringstream ss;
	ss << "Speciation parameters: m/l = " << fixed << setprecision(4) << relativeDeath 
	   << " , l-m = " << fixed << setprecision(4) << netDiversificaton
	   << " , s = " << fixed << setprecision(4) << probSpeciationS
	   << " , psi = " << fixed << setprecision(4) << fossilRate
	   << " , l = " << fixed << setprecision(4) << birthRate
	   << " , m = " << fixed << setprecision(4) << deathRate
	   
	   << endl;
	return ss.str();
}

double Speciation::getLnTreeProb(Tree *t) {
	
	if(treeTimePrior == 1)
		return 0.0;
	else if(treeTimePrior == 2){
		int ntax = t->getNumTaxa();
		double c1 = (ntax - 1) * log(netDiversificaton); //  ntax * log(1 - rel) = 0;
		double nps = t->getTreeCBDNodePriorProb(netDiversificaton, relativeDeath);
		return c1 + nps;
	}
	else if(treeTimePrior == 3){
		int ntax = t->getNumTaxa();
		double c1 = (ntax - 1) * log(netDiversificaton) + ntax * log(1 - relativeDeath);
		double nps = t->getTreeCBDNodePriorProb(netDiversificaton, relativeDeath);
		return c1 + nps;
	}
	else if(treeTimePrior == 4 || treeTimePrior == 5){
		Treescale *ts = modelPtr->getActiveTreeScale();
		double nps = t->getTreeBDSSTreeNodePriorProb(netDiversificaton, relativeDeath, fossilRate, extantSampleRate, ts->getTreeOriginTime());
		return nps;
	}
	else if(treeTimePrior == 6){
		setAllBDFossParams();
		double nps = t->getTreeCalBDSSTreeNodePriorProb(birthRate, deathRate, fossilRate, extantSampleRate);
		return nps;
	}
	else if(treeTimePrior == 7){
		setAllBDFossParams();
		double nps = t->getTreeAncCalBDSSTreeNodePriorProb(birthRate, deathRate, fossilRate, extantSampleRate);
		return nps;
	}
	return 0.0;
}


void Speciation::setAllBDFossParams(){
	fossilRate = (probSpeciationS / (1-probSpeciationS)) * ((relativeDeath * netDiversificaton) / (1 - relativeDeath));
	birthRate = netDiversificaton / (1.0 - relativeDeath); 
	deathRate = (relativeDeath * netDiversificaton) / (1 - relativeDeath);

}