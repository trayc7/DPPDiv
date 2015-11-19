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
	
	relTimePoints.push_back(1.0);
	for(int i=1; i<numIntervals; i++){
		double v = i * (1.0 / numIntervals);
		relTimePoints.push_back(1.0-v);
	}
	relTimePoints.push_back(0.0);
	
}

int Tree::getSkylineIndexForTime(double ot, double tt){
	
	int myID = 0;
	if(tt == ot)
		return 0;
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
//	sl->setAllSkyFBDParameters();
	lambdaV = sl->getSkylineBirthVec();
	muV = sl->getSkylineDeathVec();
	psiV = sl->getSkylinePsiVec();
	rho = sl->getSkylineRhoVal();
	
	nprb = calcFBDSkylineProbability();

	return nprb;
}

double Tree::getFBDSkylineProbability(double ot){

	double nprb = 0.0;
	
	originTime = ot;
	Skyline *sl = modelPtr->getActiveSkyline();
//	sl->setAllSkyFBDParameters();
	lambdaV = sl->getSkylineBirthVec();
	muV = sl->getSkylineDeathVec();
	psiV = sl->getSkylinePsiVec();
	rho = sl->getSkylineRhoVal();
	
	nprb = calcFBDSkylineProbability();

	return nprb;
}



double Tree::getFBDSkylineProbability(vector<double> bv, vector<double> dv, vector<double> pv, double r){

	double nprb = 0.0;
	
	OriginTime *ot = modelPtr->getActiveOriginTime();
	originTime = ot->getOriginTime();
	lambdaV = bv;
	muV = dv;
	psiV = pv;
	rho = r;
	
	nprb = calcFBDSkylineProbability();

	return nprb;
}



double Tree::calcFBDSkylineProbability(void){

	double nprb = 0.0;
		
	int numFossils = (int)fossSpecimens.size();
	
	double qBarOT = fbdLogQBarFxn(lambdaV[0],muV[0],psiV[0],rho,originTime);
	nprb = -(MbMath::lnFactorial(numTaxa + numFossils));
	nprb += (log(rho) + qBarOT) - log(1.0-bdssP0HatFxn(lambdaV[0],muV[0],rho,originTime));
	
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(!p->getIsLeaf()){
			double xi = p->getNodeDepth() * treeScale;
			int pIdx = getSkylineIndexForTime(originTime, xi);
			double qBarXi = fbdLogQBarFxn(lambdaV[pIdx],muV[pIdx],psiV[pIdx],rho,xi);
			nprb += log(2.0*lambdaV[pIdx]*rho) + qBarXi;
		}
	}
	
	for(int i=0; i<numFossils; i++){
		Fossil *f = fossSpecimens[i];
		double yf = f->getFossilAge();
		int yIdx = getSkylineIndexForTime(originTime, yf);
		nprb += log(psiV[yIdx]) + log(f->getFossilFossBrGamma());
		
		if(f->getFossilIndicatorVar()){
			double zf = f->getFossilSppTime() * treeScale;
			int zIdx = getSkylineIndexForTime(originTime, zf);
			double qBarYf = fbdLogQBarFxn(lambdaV[yIdx],muV[yIdx],psiV[yIdx],rho,yf);
			double qBarZf = fbdLogQBarFxn(lambdaV[zIdx],muV[zIdx],psiV[zIdx],rho,zf);
			double pFxnYf = bdssP0Fxn(lambdaV[yIdx],muV[yIdx],psiV[yIdx],rho,yf);
			nprb += log(2.0 * lambdaV[zIdx]) + qBarZf;
			nprb += log(pFxnYf) - qBarYf;
		}
	}
	return nprb;
}


double Tree::updateFBDSkylineTree(double &oldLnL) {
	
	Treescale *ts = modelPtr->getActiveTreeScale();
	setTreeScale(ts->getScaleValue());
	OriginTime *ot = modelPtr->getActiveOriginTime();
	originTime = ot->getOriginTime();
	Skyline *sl = modelPtr->getActiveSkyline();
	sl->setAllSkyFBDParameters();
	lambdaV = sl->getSkylineBirthVec();
	muV = sl->getSkylineDeathVec();
	psiV = sl->getSkylinePsiVec();
	rho = sl->getSkylineRhoVal();
	
	size_t updateNum = pickUpdate();
	if(updateNum == 0)
		updateFBDSkylineNodeAges(oldLnL);
	else if (updateNum == 1){
		updateFBDSkylineFossilAttachTimes();
		treeUpdateNodeOldestBoundsAttchTimes();
	}
	else if (updateNum == 2){
		for(int i=0; i<fossSpecimens.size(); i++){
			updateFBDSkylineRJMoveAddDelEdge();
			treeUpdateNodeOldestBoundsAttchTimes();
		}
	}
	else if (updateNum == 3)
		updateFBDSkylineFossilAges();
	
	modelPtr->setLnLGood(true);
	modelPtr->setMyCurrLnl(oldLnL);
	modelPtr->setTiProb();
	return 0.0;
}

double Tree::updateFBDSkylineNodeAges(double &oldLnL) {
		
	upDateAllCls(); 
	upDateAllTis();

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
				double lnPrRatio = lnPrRatioNodeAgeMoveFBDSky(newNodeDepth, currDepth, newSumLogGammas, oldSumLogGammas);

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
			double v1 = log(lambdaV[newIdx]) + fbdLogQBarFxn(lambdaV[newIdx],muV[newIdx],psiV[newIdx],rho,newPhi) + newSumLogGammas;
			double v2 = log(lambdaV[oldIdx]) + fbdLogQBarFxn(lambdaV[oldIdx],muV[oldIdx],psiV[oldIdx],rho,oldPhi) + oldSumLogGammas;
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
			
			double prOld = log(psiV[oldIdx]) + (log(bdssP0Fxn(lambdaV[oldIdx],muV[oldIdx],psiV[oldIdx],rho,oldAge)) - fbdLogQBarFxn(lambdaV[oldIdx],muV[oldIdx],psiV[oldIdx],rho,oldAge)) + oldSumLogGammas;
			f->setFossilAge(newAge);
			double newSumLogGammas = getSumLogAllAttachNums();

			double prNew = log(psiV[newIdx]) + (log(bdssP0Fxn(lambdaV[newIdx],muV[newIdx],psiV[newIdx],rho,newAge)) - fbdLogQBarFxn(lambdaV[newIdx],muV[newIdx],psiV[newIdx],rho,newAge)) + newSumLogGammas;
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
	double prNew = log(2.0) + log(lambdaV[newZIdx]) + fbdLogQBarFxn(lambdaV[newZIdx],muV[newZIdx],psiV[newZIdx],rho,zf);
	prNew += log(bdssP0Fxn(lambdaV[newYIdx],muV[newYIdx],psiV[newYIdx],rho,yf)) - fbdLogQBarFxn(lambdaV[newYIdx],muV[newYIdx],psiV[newYIdx],rho,yf);
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
	
	double prOld = log(2.0) + log(lambdaV[oldZIdx]) + fbdLogQBarFxn(lambdaV[oldZIdx],muV[oldZIdx],psiV[oldZIdx],rho,oldZf);
	prOld += log(bdssP0Fxn(lambdaV[oldYIdx],muV[oldYIdx],psiV[oldYIdx],rho,yf)) - fbdLogQBarFxn(lambdaV[oldYIdx],muV[oldYIdx],psiV[oldYIdx],rho,yf);
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


double Tree::lnPrRatioNodeAgeMoveFBDSky(double newD, double oldD, double newG, double oldG){
	
	int newIdx = getSkylineIndexForTime(originTime, newD);
	int oldIdx = getSkylineIndexForTime(originTime, oldD);
	double v1 = log(lambdaV[newIdx]) + fbdLogQBarFxn(lambdaV[newIdx],muV[newIdx],psiV[newIdx],rho,newD) + newG;
	double v2 = log(lambdaV[oldIdx]) + fbdLogQBarFxn(lambdaV[oldIdx],muV[oldIdx],psiV[oldIdx],rho,oldD) + oldG;
	
	return v1 - v2;
}


void Tree::setUpdateProbs(void){

	updateProbs.clear();
	if(skylineFBD){
		updateProbs.push_back(1.0); // 0 = updateFBDSkylineNodeAges
		updateProbs.push_back(1.0); // 1 = updateFBDSkylineFossilAttachTimes
		updateProbs.push_back(1.0); // 2 = updateFBDSkylineRJMoveAddDelEdge
		updateProbs.push_back(0.0); // 3 = updateFBDSkylineFossilAges
		
		if(sampleFossilAges)
			updateProbs[3] = 1.0;
	}
	else{
		updateProbs.push_back(1.0); // 0 = updateAllTGSNodes
		updateProbs.push_back(1.0); // 1 = updateFossilBDSSAttachmentTimePhi
		updateProbs.push_back(1.0); // 2 = updateRJMoveAddDelEdge
		updateProbs.push_back(0.0); // 3 = updateFossilAges
		
		if(sampleFossilAges)
			updateProbs[3] = 1.0;
        if(noSampledAncestors == true){
            updateProbs[2] = 0.0;
        }
	}
	
	double sum = 0.0;
	for(size_t i=0; i<updateProbs.size(); i++){
		sum += updateProbs[i];
	}
	for(size_t i=0; i<updateProbs.size(); i++){
		updateProbs[i] /= sum;
	}
	
}

size_t Tree::pickUpdate(void){
	
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



