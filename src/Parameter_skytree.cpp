#include "Alignment.h"
#include "Calibration.h"
#include "MbMath.h"
#include "MbRandom.h"
#include "Model.h"
#include "Parameter.h"
#include "Parameter_expcalib.h"
#include "Parameter_origin.h"
#include "Parameter_rate.h"
#include "Parameter_skyline.h"
#include "Parameter_speciaton.h"
#include "Parameter_treescale.h"
#include "Parameter_tree.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;


double Tree::fbdLogQBarFxn(double b, double d, double psi, double rho, double t){
	
	// 4 / q(t)
	double v = log(4.0);
	v -= bdssQFxn(b,d,psi,rho,t);
	return v;
}

void Tree::intitializeRelativeTimePoints(void){
	
	numIntervals = 5;
	relTimePoints.push_back(1.0);
	for(int i=1; i<numIntervals; i++){
		double v = i * (1.0 / numIntervals);
		relTimePoints.push_back(1.0-v);
	}
	relTimePoints.push_back(0.0);
	
}

int Tree::getSkylineIndexForTime(double ot, double tt){
	
	int myID = 0;
	while(relTimePoints[myID]*ot > tt){
		myID++;
	}
	return myID-1;
}


double Tree::getFBDSkylineProbability(void){

	double nprb = 0.0;
	
	OriginTime *ot = modelPtr->getActiveOriginTime();
	originTime = ot->getOriginTime();
	Skyline *sl = modelPtr->getActiveSkyline();
	vector<double> bv = sl->getSkylineBirthVec();
	vector<double> dv = sl->getSkylineDeathVec();
	vector<double> pv = sl->getSkylinePsiVec();
	double rho = sl->getSkylineRhoVal();
	
	int numFossils = (int)fossSpecimens.size();
	
	double qBarOT = fbdLogQBarFxn(bv[0],dv[0],pv[0],rho,originTime);
	nprb = -(MbMath::lnFactorial(numTaxa + numFossils));
	nprb += (log(rho) + qBarOT) - log(1.0-bdssP0HatFxn(bv[0],dv[0],rho,originTime));
	
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(!p->getIsLeaf()){
			double xi = p->getNodeDepth() * treeScale;
			int pIdx = getSkylineIndexForTime(originTime, xi);
			double qBarXi = fbdLogQBarFxn(bv[pIdx],dv[pIdx],pv[pIdx],rho,xi);
			nprb += log(2.0*bv[pIdx]*rho) + qBarXi;
		}
	}
	
	for(int i=0; i<numFossils; i++){
		Fossil *f = fossSpecimens[i];
		double yf = f->getFossilAge();
		int yIdx = getSkylineIndexForTime(originTime, yf);
		nprb += log(pv[yIdx]) + log(f->getFossilFossBrGamma());
		
		if(f->getFossilIndicatorVar()){
			double zf = f->getFossilSppTime() * treeScale;
			int zIdx = getSkylineIndexForTime(originTime, zf);
			double qBarYf = fbdLogQBarFxn(bv[yIdx],dv[yIdx],pv[yIdx],rho,yf);
			double qBarZf = fbdLogQBarFxn(bv[zIdx],dv[zIdx],pv[zIdx],rho,zf);
			double pFxnYf = bdssP0Fxn(bv[yIdx],dv[yIdx],pv[yIdx],rho,yf);
			nprb += log(2.0 * bv[zIdx]) + qBarZf;
			nprb += log(pFxnYf) - qBarYf;
		}
	}
	return nprb;
}

double Tree::updateFBDSkylineTree(double &oldLnL) {
	
	Treescale *ts = modelPtr->getActiveTreeScale();
	setTreeScale(ts->getScaleValue());

	updateFBDSkylineNodeAges(oldLnL);
	
	updateFBDSkylineFossilAttachTimes();
	treeUpdateNodeOldestBoundsAttchTimes();
	
	for(int i=0; i<fossSpecimens.size(); i++){
		updateFBDSkylineRJMoveAddDelEdge();
		treeUpdateNodeOldestBoundsAttchTimes();
	}
	if(sampleFossilAges){
		updateFBDSkylineFossilAges();
	}
	
	modelPtr->setLnLGood(true);
	modelPtr->setMyCurrLnl(oldLnL);
	modelPtr->setTiProb();
	return 0.0;
}

double Tree::updateFBDSkylineNodeAges(double &oldLnL) {
		
	upDateAllCls(); 
	upDateAllTis();
	Skyline *sl = modelPtr->getActiveSkyline();
	vector<double> bv = sl->getSkylineBirthVec();
	vector<double> dv = sl->getSkylineDeathVec();
	vector<double> pv = sl->getSkylinePsiVec();
	double rho = sl->getSkylineRhoVal();

	double oldLike = oldLnL;
	Node *p = NULL;
	vector<int> rndNodeIDs;
	for(int i=0; i<numNodes; i++)
		rndNodeIDs.push_back(i);
	random_shuffle(rndNodeIDs.begin(), rndNodeIDs.end());
	for(vector<int>::iterator it=rndNodeIDs.begin(); it!=rndNodeIDs.end(); it++){
		p = downPassSequence[(*it)];
		if(p != root && !p->getIsLeaf()){
			double currDepth = p->getNodeDepth() * treeScale;
			double oldSumLogGammas = getSumLogAllAttachNums();
			
			double largestTime = getNodeUpperBoundTime(p) * treeScale;
			double smallestTime = getNodeLowerBoundTime(p) * treeScale;
			
			if (largestTime > smallestTime){
				double rv = ranPtr->uniformRv();
				double newNodeDepth, c;
				if(nodeProposal == 1){
					double delta = 5.0;
					c = doAWindoMove(newNodeDepth, currDepth, delta, smallestTime, largestTime, rv);
				}
				else if(nodeProposal == 2){
					double tv = tuningVal;
					c = doAScaleMove(newNodeDepth, currDepth, tv, smallestTime, largestTime, rv);
				}
				else if(nodeProposal == 3){
					newNodeDepth = smallestTime + rv*(largestTime-smallestTime);
					c = 0.0;
				}
				
				p->setNodeDepth(newNodeDepth/treeScale);
				double newSumLogGammas = getSumLogAllAttachNums();
				double lnPrRatio = lnPrRatioNodeAgeMoveFBDSky(newNodeDepth, currDepth, newSumLogGammas, oldSumLogGammas, bv, dv, pv, rho);

				flipToRootClsTis(p);
				updateToRootClsTis(p);
				modelPtr->setTiProb();
				double newLnl = modelPtr->lnLikelihood();
				double lnLRatio = newLnl - oldLike;

				double lnR = lnPrRatio + lnLRatio + c;
				double r = modelPtr->safeExponentiation(lnR);
				
				if(ranPtr->uniformRv() < r){
					oldLike = newLnl;
				}
				else{
					p->setNodeDepth(currDepth/treeScale);
					flipToRootClsTis(p);
					updateToRootClsTis(p);
				}
			}
			recountFossilAttachNums();
		}
	}
	oldLnL = oldLike;
	return 0.0;
}

double Tree::updateFBDSkylineFossilAttachTimes() {
	
	Skyline *sl = modelPtr->getActiveSkyline();
	vector<double> bv = sl->getSkylineBirthVec();
	vector<double> dv = sl->getSkylineDeathVec();
	vector<double> pv = sl->getSkylinePsiVec();
	double rho = sl->getSkylineRhoVal();

	vector<int> rndFossIDs;
	for(int i=0; i<fossSpecimens.size(); i++)
		rndFossIDs.push_back(i);
	random_shuffle(rndFossIDs.begin(), rndFossIDs.end());
	for(vector<int>::iterator it=rndFossIDs.begin(); it!=rndFossIDs.end(); it++){
		Fossil *f = fossSpecimens[(*it)];
		if(f->getFossilIndicatorVar()){
			
			Node *p = &nodes[f->getFossilMRCANodeID()];
			double nodeDepth = p->getNodeDepth() * treeScale;
            if(f->getIsTotalGroupFossil()){
                if(p != root){
                    nodeDepth = p->getAnc()->getNodeDepth()*treeScale;
                }
                else {
                    OriginTime *ot = modelPtr->getActiveOriginTime();
                    nodeDepth = ot->getOriginTime();
                }
            }
			double fossDepth = f->getFossilAge();
			
			double oldPhi = f->getFossilSppTime() * treeScale;
			double oldSumLogGammas = getSumLogAllAttachNums();
			
			double rv = ranPtr->uniformRv();
			double newPhi, c;
			if(nodeProposal == 1){
				double delta = 5.0;
				c = doAWindoMove(newPhi, oldPhi, delta, fossDepth, nodeDepth, rv);
			}
			else if(nodeProposal == 2){
				double tv = tuningVal;
				c = doAScaleMove(newPhi, oldPhi, tv, fossDepth, nodeDepth, rv);
			}
			else if(nodeProposal == 3){
				newPhi = fossDepth + rv*(nodeDepth-fossDepth);
				c = 0.0;
			}


			f->setFossilSppTime(newPhi/treeScale);
						
			double newSumLogGammas = getSumLogAllAttachNums();
			int newIdx = getSkylineIndexForTime(originTime, newPhi);
			int oldIdx = getSkylineIndexForTime(originTime, oldPhi);

			double lnPriorRat = 0.0;
			double v1 = log(bv[newIdx]) + fbdLogQBarFxn(bv[newIdx],dv[newIdx],pv[newIdx],rho,newPhi) + newSumLogGammas;
			double v2 = log(bv[oldIdx]) + fbdLogQBarFxn(bv[oldIdx],dv[oldIdx],pv[oldIdx],rho,oldPhi) + oldSumLogGammas;
			lnPriorRat = (v1 - v2);
			
			
			double r = modelPtr->safeExponentiation(lnPriorRat + c);
			
			if(ranPtr->uniformRv() < r){ 
				f->setFossilSppTime(newPhi/treeScale);
				setNodeOldestAttchBranchTime(p);
			}
			else{
				f->setFossilSppTime(oldPhi/treeScale);
				setNodeOldestAttchBranchTime(p);
			}
		}
	}
	return 0.0;
}

double Tree::updateFBDSkylineFossilAges(void){
	
	Skyline *sl = modelPtr->getActiveSkyline();
	vector<double> bv = sl->getSkylineBirthVec();
	vector<double> dv = sl->getSkylineDeathVec();
	vector<double> pv = sl->getSkylinePsiVec();
	double rho = sl->getSkylineRhoVal();

	vector<int> rndFossIDs;
	for(int i=0; i<fossSpecimens.size(); i++)
		rndFossIDs.push_back(i);
	random_shuffle(rndFossIDs.begin(), rndFossIDs.end());
	for(vector<int>::iterator it=rndFossIDs.begin(); it!=rndFossIDs.end(); it++){
		Fossil *f = fossSpecimens[(*it)];
		if(f->getFossilIndicatorVar()==1 && f->getDoEstFossilAge() == true){
			double oldAge = f->getFossilAge();
			double oldSumLogGammas = getSumLogAllAttachNums();
			double aMin = f->getFossilMinAge();
			double aMax = f->getFossilMaxAge();
			if(f->getFossilSppTime()*treeScale < aMax)
				aMax = f->getFossilSppTime() * treeScale;
			double newAge = ranPtr->uniformRv(aMin, aMax);
			int newIdx = getSkylineIndexForTime(originTime, newAge);
			int oldIdx = getSkylineIndexForTime(originTime, oldAge);
			
			double prOld = log(pv[oldIdx]) + (log(bdssP0Fxn(bv[oldIdx],dv[oldIdx],pv[oldIdx],rho,oldAge)) - fbdLogQBarFxn(bv[oldIdx],dv[oldIdx],pv[oldIdx],rho,oldAge)) + oldSumLogGammas;
			f->setFossilAge(newAge);
			double newSumLogGammas = getSumLogAllAttachNums();

			double prNew = log(pv[newIdx]) + (log(bdssP0Fxn(bv[newIdx],dv[newIdx],pv[newIdx],rho,newAge)) - fbdLogQBarFxn(bv[newIdx],dv[newIdx],pv[newIdx],rho,newAge)) + newSumLogGammas;
			double prRatio = prNew - prOld;
			double r = modelPtr->safeExponentiation(prRatio);
			
			if(ranPtr->uniformRv() < r){ 
				f->setFossilAge(newAge);
			}
			else{
				f->setFossilAge(oldAge);
				recountFossilAttachNums();
			}
		}
	}
	return 0.0;
}

double Tree::updateFBDSkylineRJMoveAddDelEdge(void){
	
	double gA = 0.5;
	if(numAncFossilsk == numCalibNds)
		gA = 1.0;
	else if(numAncFossilsk == 0)
		gA = 0.0;
		
	double u = ranPtr->uniformRv();
	if(u < gA)
		doFBDSkyAddEdgeMove();
	else
		doFBDSkyDeleteEdgeMove();
	return 0.0;
}

void Tree::doFBDSkyAddEdgeMove(void){
	
	int k = numAncFossilsk;
	int m = numCalibNds;
	int kN = k-1;

	Skyline *sl = modelPtr->getActiveSkyline();
	vector<double> bv = sl->getSkylineBirthVec();
	vector<double> dv = sl->getSkylineDeathVec();
	vector<double> pv = sl->getSkylinePsiVec();
	double rho = sl->getSkylineRhoVal();

	double lnHastings, lnJacobian, lnPriorR;

	int mvFoss = pickRandAncestorFossil();
	Fossil *f = fossSpecimens[mvFoss];
	Node *p = &nodes[f->getFossilMRCANodeID()];
	double oldSumLogGammas = getSumLogAllAttachNums(); 
	
	double alA = 1.0;
	if(kN == 0)
		alA = 2.0;
	else if(k == m)
		alA = 0.5;
		
	double cf = p->getNodeDepth() * treeScale;
    if(f->getIsTotalGroupFossil()){
        if(p != root){
            cf = p->getAnc()->getNodeDepth() * treeScale;
        }
        else {
            OriginTime *ot = modelPtr->getActiveOriginTime();
            cf = ot->getOriginTime();
        }
    }
	double yf = f->getFossilAge();
	double nu = ranPtr->uniformRv() * (cf - yf);
	lnHastings = log(alA) + (log(k) - log(m - k + 1.0)); 
	lnJacobian = log(cf - yf);
	double zf = yf + nu;
	double scaleZf = zf / treeScale;
	f->setFossilSppTime(scaleZf);
	f->setFossilIndicatorVar(1);
		
	double newSumLogGammas = getSumLogAllAttachNums(); 
	
	int newYIdx = getSkylineIndexForTime(originTime, yf);
	int newZIdx = getSkylineIndexForTime(originTime, zf);
	double prNew = log(2.0) + log(bv[newZIdx]) + fbdLogQBarFxn(bv[newZIdx],dv[newZIdx],pv[newZIdx],rho,zf);
	prNew += log(bdssP0Fxn(bv[newYIdx],dv[newYIdx],pv[newYIdx],rho,yf)) - fbdLogQBarFxn(bv[newYIdx],dv[newYIdx],pv[newYIdx],rho,yf);
	prNew += newSumLogGammas;
	
	lnPriorR = prNew - oldSumLogGammas;
	
	double lpr = lnPriorR + lnHastings + lnJacobian;
	double r = modelPtr->safeExponentiation(lpr);
	
	if(ranPtr->uniformRv() < r){ 
		f->setFossilIndicatorVar(1);
		f->setFossilSppTime(scaleZf);
		setNodeOldestAttchBranchTime(p);
		numAncFossilsk = kN;
	}
	else{
		f->setFossilIndicatorVar(0);
		f->setFossilSppTime(yf/treeScale);
		numAncFossilsk = k;
		setNodeOldestAttchBranchTime(p);
	}
}

void Tree::doFBDSkyDeleteEdgeMove(void){
	
	int k = numAncFossilsk;
	int m = numCalibNds;
	int kN = k+1;

	Skyline *sl = modelPtr->getActiveSkyline();
	vector<double> bv = sl->getSkylineBirthVec();
	vector<double> dv = sl->getSkylineDeathVec();
	vector<double> pv = sl->getSkylinePsiVec();
	double rho = sl->getSkylineRhoVal();

	double lnHastings, lnJacobian, lnPriorR;

	int mvFoss = pickRandTipFossil();
	Fossil *f = fossSpecimens[mvFoss];
	Node *p = &nodes[f->getFossilMRCANodeID()];
	double oldSumLogGammas = getSumLogAllAttachNums(); 

	double alD = 1.0;
	if(k == 0)
		alD = 0.5;
	else if(kN == m)
		alD = 2.0;

	double cf = p->getNodeDepth() * treeScale;
	double yf = f->getFossilAge();
	
	double oldZf = f->getFossilSppTime() * treeScale;
	double oldZfScale = f->getFossilSppTime();
	double newZfScale = yf / treeScale;
	f->setFossilIndicatorVar(0);
	f->setFossilSppTime(newZfScale);
	
	double newSumLogGammas = getSumLogAllAttachNums(); 

	int oldYIdx = getSkylineIndexForTime(originTime, yf);
	int oldZIdx = getSkylineIndexForTime(originTime, oldZf);

	lnHastings = log(alD) + (log(m-k) - log(k+1));
	lnJacobian = -(log(cf - yf));
	
	double prOld = log(2.0) + log(bv[oldZIdx]) + fbdLogQBarFxn(bv[oldZIdx],dv[oldZIdx],pv[oldZIdx],rho,oldZf);
	prOld += log(bdssP0Fxn(bv[oldYIdx],dv[oldYIdx],pv[oldYIdx],rho,yf)) - fbdLogQBarFxn(bv[oldYIdx],dv[oldYIdx],pv[oldYIdx],rho,yf);
	prOld += oldSumLogGammas;
		
	lnPriorR = newSumLogGammas - prOld;
	
	double lpr = lnPriorR + lnHastings + lnJacobian;
	double r = modelPtr->safeExponentiation(lpr);

	if(ranPtr->uniformRv() < r){ 
		f->setFossilIndicatorVar(0);
		f->setFossilSppTime(newZfScale);
		setNodeOldestAttchBranchTime(p);
		numAncFossilsk = kN;
	}
	else{
		f->setFossilSppTime(oldZfScale);
		f->setFossilIndicatorVar(1);
		setNodeOldestAttchBranchTime(p);
	}
}


double Tree::lnPrRatioNodeAgeMoveFBDSky(double newD, double oldD, double newG, double oldG, vector<double> bv, vector<double> dv, vector<double> pv, double rho){
	
	int newIdx = getSkylineIndexForTime(originTime, newD);
	int oldIdx = getSkylineIndexForTime(originTime, oldD);
	double v1 = log(bv[newIdx]) + fbdLogQBarFxn(bv[newIdx],dv[newIdx],pv[newIdx],rho,newD) + newG;
	double v2 = log(bv[oldIdx]) + fbdLogQBarFxn(bv[oldIdx],dv[oldIdx],pv[oldIdx],rho,oldD) + oldG;
	
	return v1 - v2;
}




