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
#include "Parameter_tree.h"
#include "MbMath.h"
#include "MbRandom.h"
#include "Model.h"
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

Treescale::Treescale(MbRandom *rp, Model *mp, double sv, double yb, double ob, int dt, bool calib, bool exhpc) : Parameter(rp, mp) {
	oldBound = ob;
	yngBound = yb;
	scaleVal = sv;
	name = "RH";
	treeTimePrior = modelPtr->getTreeTimePriorNum();
	retune = false;
	tsPriorD = dt;
	isUserCalib = calib;
	exponHPCalib = exhpc;
	expoRate = 1.0;
	
	tOrigExpHPRate = 1.0 / (sv * 0.5);
	treeOriginTime = ranPtr->exponentialRv(tOrigExpHPRate) + sv;
	
	if (treeTimePrior > 5)
		tsPriorD = 3;
	
	
	tuning = sv * 0.1;
	if(tsPriorD == 1){
		exponHPCalib = false;
		if(yngBound > 0.0 && oldBound > 0.0){
			isBounded = true;
			if(yngBound != oldBound){
				tuning = ((oldBound + yngBound) * 0.5) * 0.1;
				double windowSize = 2 * tuning;
				double calibSize = oldBound - yngBound;
				while(windowSize > calibSize * 0.85){
					tuning = tuning * 0.9;
					windowSize = 2 * tuning;
				}
				cout << "Calib window = " << calibSize << "  Tuning window = " << windowSize << endl;
			}
		}
	}
	else if(tsPriorD == 2){
		isBounded = true;

		if(isUserCalib)
			expoRate = modelPtr->getRootNExpRate();
		else{
			expoRate = 1.0 / (yngBound * 200.0);
			exponHPCalib = false;
		}
		oldBound = 100000.0;
		tuning = ((yngBound + (sv * 1.2)) * 0.5) * 0.1;
		double windowSize = 2 * tuning;
		cout << "Exponential rate = " << expoRate << "  E[TS] = " << (1 / expoRate) + yngBound
			<< "  offset = " << yngBound << "  Tuning window = " << windowSize << endl;
	}
	else if(tsPriorD > 2 && treeTimePrior > 5){ 
		isBounded = true;
		tsPriorD = 1;
		oldBound = 100000.0;
		tuning = ((yngBound + (sv * 1.2)) * 0.5) * 0.1;
		double windowSize = 2.0 * tuning;
		double calibSize = oldBound - yngBound;
		
		cout << "Root min = " << yngBound << " window size = " << calibSize << "  Tuning window = " << windowSize << endl;
	}
	

}

Treescale::~Treescale(void) {
	
}

Treescale& Treescale::operator=(const Treescale &c) {
	
	if (this != &c)
		clone(c);
	return *this;
}

void Treescale::clone(const Treescale &c) {

	scaleVal = c.scaleVal;
	treeOriginTime = c.treeOriginTime;
	tOrigExpHPRate = c.tOrigExpHPRate;
	oldBound = c.oldBound;
	yngBound = c.yngBound;
	isBounded = c.isBounded;
	tuning = c.tuning;
	name = "RH";
}

void Treescale::print(std::ostream & o) const {

	o << "Root height parameter: ";
	o << fixed << setprecision(4) << scaleVal << " ";
	o << endl;
}

double Treescale::update(double &oldLnL) {
	
	double lppr = 0.0;
	if(treeTimePrior == 4 && ranPtr->uniformRv() < 0.3){
		lppr = updateTreeOrigTime(oldLnL);
	}
	else{
		lppr = updateTreeScalePropSE(oldLnL); 
	}
	return lppr;
}

double Treescale::updateTreeScale(double &oldLnL) {
	
	Tree *t = modelPtr->getActiveTree();
	Node *rt = t->getRoot();

	double oldtreeprob = getLnTreeProb(t);
	
	if(retune && tuning > scaleVal)
		tuning = scaleVal * 0.5;
		
	double limO = scaleVal + tuning;
	double limY = scaleVal - tuning;
	double rtLB = t->getNodeLowerBoundTime(rt) * scaleVal;
	double lowBound = rtLB;
	
	double hiBound = limO;
	if(treeTimePrior == 4){
		if(hiBound > treeOriginTime)
			hiBound = treeOriginTime;
	}
	
	if(isBounded){
		if(lowBound < yngBound)
			lowBound = yngBound;
		if(hiBound > oldBound)
			hiBound = oldBound;
	}
	
	double oldRH, newRH;
	oldRH = scaleVal;
		
	double u = ranPtr->uniformRv(-0.5,0.5) * (limO - limY);
	newRH = oldRH + u;
	while(newRH < lowBound || newRH > hiBound){
		if(newRH < lowBound)
			newRH = (2 * lowBound) - newRH;
		if(newRH > hiBound)
			newRH = (2 * hiBound) - newRH;
	}
		
	double scaleRatio = oldRH / newRH;
	int numNodes = t->getNumNodes();
	for(int i=0; i<numNodes; i++){
		Node *p = t->getNodeByIndex(i);
		if(p != rt){
			if(p->getIsLeaf() && p->getIsCalibratedDepth() == false){
				p->setNodeDepth(0.0);
			}
			else{
				double oldP = p->getNodeDepth();
				double newP = oldP * scaleRatio;
				p->setNodeDepth(newP);
				p->setFossAttchTime(0.0);
				if(false){
					double oldPhi = p->getFossAttchTime();
					double fA = p->getNodeYngTime();
					double newPhi = oldPhi * scaleRatio;
					if(newPhi*newRH < fA){
						cout << i << " op = " << oldPhi*oldRH << "  np = " << newPhi*newRH << " oRH = " << oldRH << "  nRH = " 
							<< newRH << "  fA = " << fA << endl;
					}
					int x = t->countDecLinsTimeIntersect(p, newPhi, newP);
					p->setFossAttchTime(newPhi);
					p->setNumFossAttchLins(x); 
				}
			}
		}
	}
	scaleVal = newRH;
	t->setTreeScale(scaleVal);
	t->treeScaleUpdateFossilAttchTimes(scaleRatio, oldRH, newRH);
	t->treeUpdateNodeOldestBoundsAttchTimes();
	t->setAllNodeBranchTimes();
		
	double lnPriorRatio = 0.0; 
	double newtreeprob = getLnTreeProb(t);
	lnPriorRatio += (newtreeprob - oldtreeprob);
	if(tsPriorD == 2){
		if(exponHPCalib)
			expoRate = t->getRootCalibExpRate();
		lnPriorRatio += lnExponentialTSPriorRatio(newRH, oldRH);
	}
	double lnProposalRatio = 0.0;
	
	double jacobian = 0.0;
	if(treeTimePrior < 2)
		jacobian = (log(oldRH) - log(newRH)) * (t->getNumTaxa() - 2);


	t->flipAllCls();
	t->flipAllTis();
	t->upDateAllCls();
	t->upDateAllTis();
	modelPtr->setTiProb();
	
	
	return lnPriorRatio + lnProposalRatio + jacobian;
}

double Treescale::updateTreeScalePropSE(double &oldLnL) {
	
	Tree *t = modelPtr->getActiveTree();
	Node *rt = t->getRoot();
	
	double oldtreeprob = getLnTreeProb(t);
		
	double rtLB = t->getNodeLowerBoundTime(rt) * scaleVal;
	double lowBound = rtLB;
	
	double hiBound = 100000.0;
	if(treeTimePrior == 4){
		if(hiBound > treeOriginTime)
			hiBound = treeOriginTime;
	}
	
	if(isBounded){
		if(lowBound < yngBound)
			lowBound = yngBound;
		if(hiBound > oldBound)
			hiBound = oldBound;
	}
	
	double oldRH, newRH;
	oldRH = scaleVal;
	
	double rv = ranPtr->uniformRv();
	double tv = log(2.0);
	double c = tv * (rv - 0.5);
	newRH = oldRH * exp(c);
	bool validV = false;
	do{
		if(newRH < lowBound)
			newRH = lowBound * lowBound / newRH;
		else if(newRH > hiBound)
			newRH = hiBound * hiBound / newRH;
		else
			validV = true;
	} while(!validV);
	
	double scaleRatio = oldRH / newRH;
	int numNodes = t->getNumNodes();
	for(int i=0; i<numNodes; i++){
		Node *p = t->getNodeByIndex(i);
		if(p != rt){
			if(p->getIsLeaf() && p->getIsCalibratedDepth() == false){
				p->setNodeDepth(0.0);
			}
			else{
				double oldP = p->getNodeDepth();
				double newP = oldP * scaleRatio;
				p->setNodeDepth(newP);
				p->setFossAttchTime(0.0);
			}
		}
	}
	scaleVal = newRH;
	t->setTreeScale(scaleVal);
	t->treeScaleUpdateFossilAttchTimes(scaleRatio, oldRH, newRH);
	t->treeUpdateNodeOldestBoundsAttchTimes();
	t->setAllNodeBranchTimes();
	
	double lnPriorRatio = 0.0;
	double newtreeprob = getLnTreeProb(t);
	lnPriorRatio += (newtreeprob - oldtreeprob);
	if(tsPriorD == 2){
		if(exponHPCalib)
			expoRate = t->getRootCalibExpRate();
		lnPriorRatio += lnExponentialTSPriorRatio(newRH, oldRH);
	}
	double lnProposalRatio = c;
	
	double jacobian = 0.0; 
	if(treeTimePrior < 2)
		jacobian = (log(oldRH) - log(newRH)) * (t->getNumTaxa() - 2);
	
	t->flipAllCls();
	t->flipAllTis();
	t->upDateAllCls();
	t->upDateAllTis();
	modelPtr->setTiProb();
	
	
	return lnPriorRatio + lnProposalRatio + jacobian;
}


double Treescale::updateTreeOrigTime(double &oldLnL) {
	
	Tree *t = modelPtr->getActiveTree();
	
	double oldtreeprob = getLnTreeProb(t);
	
	
	double limO = treeOriginTime + tuning;
	double limY = treeOriginTime - tuning;
	double lowBound = scaleVal;
	double hiBound = limO;

	double oldTO, newTO;
	oldTO = treeOriginTime;
	
	double u = ranPtr->uniformRv(-0.5,0.5) * (limO - limY);
	newTO = oldTO + u;
	while(newTO < lowBound || newTO > hiBound){
		if(newTO < lowBound)
			newTO = (2 * lowBound) - newTO;
		if(newTO > hiBound)
			newTO = (2 * hiBound) - newTO;
	}
		
	treeOriginTime = newTO;
	
	double lnPriorRatio = 0.0; 
	double newtreeprob = getLnTreeProb(t);
	lnPriorRatio += (newtreeprob - oldtreeprob);
	
	double lnProposalRatio = 0.0; 
	
	double lnR = lnPriorRatio + lnProposalRatio;
	
	modelPtr->setLnLGood(true);
	modelPtr->setMyCurrLnl(oldLnL);
	return lnR;

}


double Treescale::lnPrior(void) {
	
	return 0.0;
}

double Treescale::lnExponentialTSPriorRatio(double newTS, double oldTS) {
	
	double offSt = yngBound;
	double numr = 0.0;
	double dnom = 0.0;
	numr = ranPtr->lnExponentialPdf(expoRate, newTS - offSt);
	dnom = ranPtr->lnExponentialPdf(expoRate, oldTS - offSt);
	return numr - dnom;
}

double Treescale::lnExponentialTreeOrigPriorRatio(double newTO, double oldTO) {
	
	double offSt = scaleVal;
	double numr = 0.0;
	double dnom = 0.0;
	numr = ranPtr->lnExponentialPdf(tOrigExpHPRate, newTO - offSt);
	dnom = ranPtr->lnExponentialPdf(tOrigExpHPRate, oldTO - offSt);
	return numr - dnom;
}


string Treescale::writeParam(void){
	
	stringstream ss;
	ss << "Root height parameter: " << scaleVal << " [+/- " << tuning << "]" << endl;
	return ss.str();
}

double Treescale::getLnTreeProb(Tree *t) {
	
	if(treeTimePrior == 1)
		return 0.0;
	else if(treeTimePrior == 2 || treeTimePrior == 3){
		
		return t->getTreeCBDNodePriorProb();
	}
	else if(treeTimePrior == 4 || treeTimePrior == 5){
		Speciation *s = modelPtr->getActiveSpeciation();
		double nps = t->getTreeBDSSTreeNodePriorProb(s->getNetDiversification(), 
													 s->getRelativeDeath(), s->getBDSSFossilSampRatePsi(), 
													 s->getBDSSSppSampRateRho(), treeOriginTime);
		return nps;
	}
	else if(treeTimePrior == 6){
		Speciation *s = modelPtr->getActiveSpeciation();
		s->setAllBDFossParams();
		double nps = t->getTreeCalBDSSTreeNodePriorProb(s->getBDSSSpeciationRateLambda(), 
														s->getBDSSExtinctionRateMu(), 
														s->getBDSSFossilSampRatePsi(), 
														s->getBDSSSppSampRateRho());
		return nps;
	}
	else if(treeTimePrior == 7){
		Speciation *s = modelPtr->getActiveSpeciation();
		s->setAllBDFossParams();
		double nps = t->getTreeAncCalBDSSTreeNodePriorProb(s->getBDSSSpeciationRateLambda(), 
														s->getBDSSExtinctionRateMu(), 
														s->getBDSSFossilSampRatePsi(), 
														s->getBDSSSppSampRateRho());
		return nps;
	}
	return 0.0;
}






