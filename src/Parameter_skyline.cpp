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

Skyline::Skyline(MbRandom *rp, Model *mp, int nrts, bool tree) : Parameter(rp, mp) {
	
	isTreeModel = tree;
	rho = 1.0;
	numRates = nrts; // this is the number of rates, which is 1 less than the number of time points
	name = "SL";
	rateMax = 1000.0;
	netDivExpPriorRt = 1.0;
	psiExpPriorRt = 100.0;
	
	for(int i=0; i<numRates; i++){
		lambdas.push_back(ranPtr->exponentialRv(1.0));
		mus.push_back(ranPtr->exponentialRv(10.0));
		psis.push_back(ranPtr->exponentialRv(10.0));
	}
//	lambdas[1] = 0.6;
//	lambdas[2] = 0.8;
//	lambdas[3] = 1.0;
//	lambdas[4] = 0.4;
//
//	mus[1] = 0.2;
//	mus[2] = 0.3;
//	mus[3] = 0.4;
//	mus[4] = 0.1;

	for(int i=0; i<numRates; i++){
		netDivVec.push_back(ranPtr->exponentialRv(5.0));
		turnoverVec.push_back(ranPtr->betaRv(2.0, 5.0));
	}
	
	

	parameterization = 2; // 1 = lambda, mu, psi; 2 = r, d, psi
	
	initializeAsEqual();
	
	setAllSkyFBDParameters();
	
	setUpdateProbs();
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
	for(int i=0; i<numRates; i++){
		lambdas[i] = c.lambdas[i];
		mus[i] = c.mus[i];
		psis[i] = c.psis[i];
		netDivVec[i] = c.netDivVec[i];
		turnoverVec[i] = c.turnoverVec[i];
	}
	name = "SL";
}

void Skyline::print(std::ostream & o) const {
	
}

double Skyline::update(double &oldLnL) {
	
	double lnR = 0.0;
	
	if(isTreeModel)
		lnR = updateTreeSkyParams(oldLnL);
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


string Skyline::writeSkylineParamLabels(void){
	stringstream ss;
	
	for(int i=0; i<numRates; i++)
		ss << "\tFBDSky.lambda[" << i << "]";
	for(int i=0; i<numRates; i++)
		ss << "\tFBDSky.mu[" << i << "]";
	for(int i=0; i<numRates; i++)
		ss << "\tFBDSky.diversification[" << i << "]";
	for(int i=0; i<numRates; i++)
		ss << "\tFBDSky.turnover[" << i << "]";
	for(int i=0; i<numRates; i++)
		ss << "\tFBDSky.psi[" << i << "]";
	ss << "\tFBDSky.rho";
	string pi = ss.str();
	return pi;
}

string Skyline::writeSkylineParamValues(void){
	stringstream ss;
	
	setAllSkyFBDParameters();
	for(int i=0; i<numRates; i++)
		ss << "\t" << lambdas[i];
	for(int i=0; i<numRates; i++)
		ss << "\t" << mus[i];
	for(int i=0; i<numRates; i++)
		ss << "\t" << netDivVec[i];
	for(int i=0; i<numRates; i++)
		ss << "\t" << turnoverVec[i];
	for(int i=0; i<numRates; i++)
		ss << "\t" << psis[i];
	ss << "\t" << rho;
	string pi = ss.str();
	return pi;
}

void Skyline::setAllSkyFBDParameters(void){
	
	if(parameterization == 1){
		for(int i=0; i<numRates; i++){
			netDivVec[i] = lambdas[i] - mus[i];
			turnoverVec[i] = mus[i] / lambdas[i];
		}
	}
	else if(parameterization == 2){
		for(int i=0; i<numRates; i++){
			lambdas[i] = netDivVec[i] / (1.0 - turnoverVec[i]);
			mus[i] = (turnoverVec[i] * netDivVec[i]) / (1.0 - turnoverVec[i]);
		}
	}
}

void Skyline::setUpdateProbs(void){
	
	updateProbs.clear();
	updateProbs.push_back(0.0); // 0 = update all turnovers
	updateProbs.push_back(0.0); // 1 = scale all turnovers
	updateProbs.push_back(0.0); // 2 = update all net divs
	updateProbs.push_back(0.0); // 3 = scale all net divs
	updateProbs.push_back(0.0); // 4 = update all psis
	updateProbs.push_back(0.0); // 5 = scale all psis
	updateProbs.push_back(1.0); // 6 = scale random turnover
	updateProbs.push_back(1.0); // 7 = scale random net div
	updateProbs.push_back(1.0); // 8 = scale random psi
	
	double sum = 0.0;
	for(size_t i=0; i<updateProbs.size(); i++){
		sum += updateProbs[i];
	}
	for(size_t i=0; i<updateProbs.size(); i++){
		updateProbs[i] /= sum;
	}
}

double Skyline::getNewValScaleMv(double &nv, double ov, double vmin, double vmax, double tv){
	
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
	return log(newcv) - log(ov);
}

double Skyline::getNewValSWindoMv(double &nv, double ov, double vmin, double vmax, double tv){
	double newv = 0.0;
	double u = ranPtr->uniformRv(-0.5,0.5) * (tv);
	newv = ov + u;
	bool validV = false;
	do{
		if(newv < vmin)
			newv = 2.0 * vmin - newv;
		else if(newv > vmax)
			newv = (2.0 * vmax) - newv;
		else
			validV = true;
	}while(!validV);
	nv = newv;
	return 0.0;
}

double Skyline::getNewValUpDownScaleMv(double &nv1, double ov1, double &nv2, double ov2, double sf){
	
	double u = ranPtr->uniformRv();
	double c = exp(sf * (u - 0.5));
	nv1 = ov1 * c;
	nv2 = ov2 / c;
	
	return 0.0;
}

double Skyline::updateTreeSkyParams(double &oldLnL) {
	
	double lnR = 0.0;
	
	size_t updateNum = pickUpdate();
	
	Tree *t = modelPtr->getActiveTree();
	if(updateNum == 0)
		udateTurnovers(t);
	else if (updateNum == 1)
		updateVectorScaleTurnover(t);
	else if (updateNum == 2)
		udateNetDivs(t);
	else if (updateNum == 3)
		updateVectorScaleNetDivs(t);
	else if (updateNum == 4)
		udatePsis(t);
	else if (updateNum == 5)
		updateVectorScalePsis(t);
	else if (updateNum == 6)
		updateRandTurnover(t);
	else if (updateNum == 7)
		updateRandNetDiv(t);
	else if (updateNum == 8)
		updateRandPsi(t);
	
	setAllSkyFBDParameters();
	modelPtr->setLnLGood(true);
	modelPtr->setMyCurrLnl(oldLnL);
	return lnR;
}


double Skyline::udateTurnovers(Tree *t){
	
	double delta = 0.3;
	
	vector<int> tmp;
	for(int i=0; i<numRates; i++)
		tmp.push_back(i);
	random_shuffle(tmp.begin(), tmp.end());
	for(int i=0; i<numRates; i++){
		int vecIDX = tmp[i];
		double oldTO = turnoverVec[vecIDX];
		double oldFBDProb = getTreeSpeciationProb(t);
		double newTO;
		double propRto = getNewValScaleMv(newTO, oldTO, 0.000001, 0.99999, delta);
		turnoverVec[vecIDX] = newTO;
		double newFBDProb = getTreeSpeciationProb(t);
		double prR = (newFBDProb - oldFBDProb) + propRto;
		double r = modelPtr->safeExponentiation(prR);
		if(ranPtr->uniformRv() < r){
			lambdas[vecIDX] = netDivVec[vecIDX] / (1.0 - turnoverVec[vecIDX]);
			mus[vecIDX] = (turnoverVec[vecIDX] * netDivVec[vecIDX]) / (1.0 - turnoverVec[vecIDX]);
		}
		else{
			turnoverVec[vecIDX] = oldTO;
			lambdas[vecIDX] = netDivVec[vecIDX] / (1.0 - turnoverVec[vecIDX]);
			mus[vecIDX] = (turnoverVec[vecIDX] * netDivVec[vecIDX]) / (1.0 - turnoverVec[vecIDX]);
		}
	}
	return 0.0;
}


double Skyline::updateVectorScaleTurnover(Tree *t){
	
	double delta = 0.1;
	double oldFBDProb = getTreeSpeciationProb(t);
	vector<double> stored;
	for(int i=0; i<numRates; i++)
		stored.push_back(turnoverVec[i]);
	double u = ranPtr->uniformRv();
	double sf = exp(delta * (u-0.5));
	bool reject = false;
	for(int i=0; i<numRates; i++){
		double ov = turnoverVec[i];
		double nv = ov * sf;
		if(nv > 0.99999 || nv < 0.000001){
			reject = true;
			break;
		}
		turnoverVec[i] = nv;
	}
	if(reject){
		for(int i=0; i<numRates; i++)
			turnoverVec[i] = stored[i];
	}
	else{
		double hr = (double)numRates * log(sf);
		setAllSkyFBDParameters();
		double newFBDProb = getTreeSpeciationProb(t);
		double prR = (newFBDProb - oldFBDProb) + hr;
		double r = modelPtr->safeExponentiation(prR);
		if(ranPtr->uniformRv() > r){ //reject
			for(int i=0; i<numRates; i++)
				turnoverVec[i] = stored[i];
			setAllSkyFBDParameters();
		}
	}
	
	return 0.0;
}

double Skyline::udateNetDivs(Tree *t){
	
	double delta = 0.3;
	
	vector<int> tmp;
	for(int i=0; i<numRates; i++)
		tmp.push_back(i);
	random_shuffle(tmp.begin(), tmp.end());
	for(int i=0; i<numRates; i++){
		int vecIDX = tmp[i];
		double oldND = netDivVec[vecIDX];
		double oldFBDProb = getTreeSpeciationProb(t);
		double newND;
		double propRto = getNewValScaleMv(newND, oldND, 0.000001, rateMax, delta);
		netDivVec[vecIDX] = newND;
		double newFBDProb = getTreeSpeciationProb(t);
		double prR = (newFBDProb - oldFBDProb) + propRto;
		prR += getExponentialPriorRatio(oldND, newND, netDivExpPriorRt);
		double r = modelPtr->safeExponentiation(prR);
		if(ranPtr->uniformRv() < r){
			lambdas[vecIDX] = netDivVec[vecIDX] / (1.0 - turnoverVec[vecIDX]);
			mus[vecIDX] = (turnoverVec[vecIDX] * netDivVec[vecIDX]) / (1.0 - turnoverVec[vecIDX]);
		}
		else{
			netDivVec[vecIDX] = oldND;
			lambdas[vecIDX] = netDivVec[vecIDX] / (1.0 - turnoverVec[vecIDX]);
			mus[vecIDX] = (turnoverVec[vecIDX] * netDivVec[vecIDX]) / (1.0 - turnoverVec[vecIDX]);
		}
	}
	return 0.0;
}


double Skyline::updateVectorScaleNetDivs(Tree *t){
	
	double delta = 0.1;
	double expPrioR = 0.0;
	double oldFBDProb = getTreeSpeciationProb(t);
	vector<double> stored;
	for(int i=0; i<numRates; i++)
		stored.push_back(netDivVec[i]);
	double u = ranPtr->uniformRv();
	double sf = exp(delta * (u-0.5));
	bool reject = false;
	for(int i=0; i<numRates; i++){
		double ov = netDivVec[i];
		double nv = ov * sf;
		if(nv > rateMax || nv < 0.000001){
			reject = true;
			break;
		}
		netDivVec[i] = nv;
		expPrioR += getExponentialPriorRatio(ov, nv, netDivExpPriorRt);
	}
	if(reject){
		for(int i=0; i<numRates; i++)
			netDivVec[i] = stored[i];
	}
	else{
		double hr = (double)numRates * log(sf);
		setAllSkyFBDParameters();
		double newFBDProb = getTreeSpeciationProb(t);
		double prR = (newFBDProb - oldFBDProb) + hr;
		prR += expPrioR;
		double r = modelPtr->safeExponentiation(prR);
		if(ranPtr->uniformRv() > r){ //reject
			for(int i=0; i<numRates; i++)
				netDivVec[i] = stored[i];
			setAllSkyFBDParameters();
		}
	}
	
	return 0.0;
}

double Skyline::udatePsis(Tree *t){
	
	double delta = 0.3;
	
	vector<int> tmp;
	for(int i=0; i<numRates; i++)
		tmp.push_back(i);
	random_shuffle(tmp.begin(), tmp.end());
	for(int i=0; i<numRates; i++){
		int vecIDX = tmp[i];
		double oldPsi = psis[vecIDX];
		double oldFBDProb = getTreeSpeciationProb(t);
		double newPsi;
		double propRto = getNewValScaleMv(newPsi, oldPsi, 0.000001, rateMax, delta);
		psis[vecIDX] = newPsi;
		double newFBDProb = getTreeSpeciationProb(t);
		double prR = (newFBDProb - oldFBDProb) + propRto;
		prR += getExponentialPriorRatio(oldPsi, newPsi, psiExpPriorRt);
		double r = modelPtr->safeExponentiation(prR);
		if(ranPtr->uniformRv() > r){
			psis[vecIDX] = oldPsi;
		}
	}
	return 0.0;
}


double Skyline::updateVectorScalePsis(Tree *t){
	
	double delta = 0.1;
	double expPrioR = 0.0;
	double oldFBDProb = getTreeSpeciationProb(t);
	vector<double> stored;
	for(int i=0; i<numRates; i++)
		stored.push_back(psis[i]);
	double u = ranPtr->uniformRv();
	double sf = exp(delta * (u-0.5));
	bool reject = false;
	for(int i=0; i<numRates; i++){
		double ov = psis[i];
		double nv = ov * sf;
		if(nv > rateMax || nv < 0.000001){
			reject = true;
			break;
		}
		psis[i] = nv;
		expPrioR += getExponentialPriorRatio(ov, nv, psiExpPriorRt);
	}
	if(reject){
		for(int i=0; i<numRates; i++)
			psis[i] = stored[i];
	}
	else{
		double hr = (double)numRates * log(sf);
		setAllSkyFBDParameters();
		double newFBDProb = getTreeSpeciationProb(t);
		double prR = (newFBDProb - oldFBDProb) + hr;
		prR += expPrioR;
		double r = modelPtr->safeExponentiation(prR);
		if(ranPtr->uniformRv() > r){ //reject
			for(int i=0; i<numRates; i++)
				psis[i] = stored[i];
			setAllSkyFBDParameters();
		}
	}
	
	return 0.0;
}

double Skyline::updateRandTurnover(Tree *t){
	
	double delta = 0.3;
	int vecIDX = (int)(ranPtr->uniformRv() * numRates);
	double ov = turnoverVec[vecIDX];
	double oldFBDProb = getTreeSpeciationProb(t);
	double nv;
	double propRto = getNewValScaleMv(nv, ov, 0.000001, 0.99999, delta);
	turnoverVec[vecIDX] = nv;
	double newFBDProb = getTreeSpeciationProb(t);
	double prR = (newFBDProb - oldFBDProb) + propRto;
	double r = modelPtr->safeExponentiation(prR);
	if(ranPtr->uniformRv() < r){
		lambdas[vecIDX] = netDivVec[vecIDX] / (1.0 - turnoverVec[vecIDX]);
		mus[vecIDX] = (turnoverVec[vecIDX] * netDivVec[vecIDX]) / (1.0 - turnoverVec[vecIDX]);
	}
	else{
		turnoverVec[vecIDX] = ov;
		lambdas[vecIDX] = netDivVec[vecIDX] / (1.0 - turnoverVec[vecIDX]);
		mus[vecIDX] = (turnoverVec[vecIDX] * netDivVec[vecIDX]) / (1.0 - turnoverVec[vecIDX]);
	}
	return 0.0;
}

double Skyline::updateRandNetDiv(Tree *t){
	
	double delta = 0.3;
	int vecIDX = (int)(ranPtr->uniformRv() * numRates);
	double ov = netDivVec[vecIDX];
	double oldFBDProb = getTreeSpeciationProb(t);
	double nv;
	double propRto = getNewValScaleMv(nv, ov, 0.000001, rateMax, delta);
	netDivVec[vecIDX] = nv;
	double newFBDProb = getTreeSpeciationProb(t);
	double prR = (newFBDProb - oldFBDProb) + propRto;
	double r = modelPtr->safeExponentiation(prR);
	if(ranPtr->uniformRv() < r){
		lambdas[vecIDX] = netDivVec[vecIDX] / (1.0 - turnoverVec[vecIDX]);
		mus[vecIDX] = (turnoverVec[vecIDX] * netDivVec[vecIDX]) / (1.0 - turnoverVec[vecIDX]);
	}
	else{
		netDivVec[vecIDX] = ov;
		lambdas[vecIDX] = netDivVec[vecIDX] / (1.0 - turnoverVec[vecIDX]);
		mus[vecIDX] = (turnoverVec[vecIDX] * netDivVec[vecIDX]) / (1.0 - turnoverVec[vecIDX]);
	}
	return 0.0;
}

double Skyline::updateRandPsi(Tree *t){
	
	double delta = 0.3;
	int vecIDX = (int)(ranPtr->uniformRv() * numRates);
	double ov = turnoverVec[vecIDX];
	double oldFBDProb = getTreeSpeciationProb(t);
	double nv;
	double propRto = getNewValScaleMv(nv, ov, 0.000001, rateMax, delta);
	psis[vecIDX] = nv;
	double newFBDProb = getTreeSpeciationProb(t);
	double prR = (newFBDProb - oldFBDProb) + propRto;
	double r = modelPtr->safeExponentiation(prR);
	if(ranPtr->uniformRv() > r){
		psis[vecIDX] = ov;
	}
	return 0.0;
}


double Skyline::getExponentialPriorRatio(double ov, double nv, double r) {
    
	return (ranPtr->lnExponentialPdf(r, nv)) - (ranPtr->lnExponentialPdf(r, ov));
}


size_t Skyline::pickUpdate(void){
	
	double u = ranPtr->uniformRv();
	double sum = 0.0;
	size_t n = 10;
	for (size_t i=0; i<updateProbs.size(); i++){
		sum += updateProbs[i];
		if ( u < sum ){
			n = i;
			break;
		}
	}
	return n;
}

double Skyline::getTreeSpeciationProb(Tree *t){

	double pr;
	setAllSkyFBDParameters();
	pr = t->getFBDSkylineProbability(lambdas, mus, psis, rho);
	return pr;
}

void Skyline::initializeAsEqual(void){
	
	// for debugging
	if(parameterization == 1){
		for(int i=0; i<numRates; i++){
			lambdas[i] = 0.5;
			mus[i] = 0.1;
			psis[i] = 0.05;
		}
	}
	else if(parameterization == 2){
		for(int i=0; i<numRates; i++){
			netDivVec[i] = 0.5;
			turnoverVec[i] = 0.1;
			psis[i] = 0.05;
		}
	}
}



