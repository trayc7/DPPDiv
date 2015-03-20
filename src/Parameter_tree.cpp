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

#include "Alignment.h"
#include "Calibration.h"
#include "MbMath.h"
#include "MbRandom.h"
#include "Model.h"
#include "Parameter.h"
#include "Parameter_expcalib.h"
#include "Parameter_origin.h"
#include "Parameter_rate.h"
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

Node::Node(void) {

	lft = NULL;
	rht = NULL;
	anc = NULL;
	idx = 0;
	name = "";
	activeCl = 0;
	activeTi = 0;
	isClDirty = true;
	isTiDirty = true;
	nodeDepth = 0.0;
	isLeaf = false;
	isCalib = false;
	youngt = 0.0;
	oldt = -1.0;
	branchTime = 0.0;
	rateGVal = 0.0;
	rateGrpIdx = -1;
	userBL = 0.0;
	nodeCalibPrD = 1;
	nodeExpCalRate = -1.0;
	taintFossil = false;
	redFlag = 0;
	nodeAge = 0.0;
	fossAttachTime = -1.0;
	numFossAttachLins = -1;
	numCalFossils = 0;
}

Tree::Tree(MbRandom *rp, Model *mp, Alignment *ap, string ts, bool ubl, bool allnm, 
		   bool rndNods, vector<Calibration *> clb, double rth, double iot, bool sb, bool exhpc,
		   ExpCalib *ec, vector<Calibration *> tdt) : Parameter(rp, mp) {

	alignmentPtr = ap;
	numTaxa = 0;
	numNodes = 0;
	nodes = NULL;
	root = NULL;
	downPassSequence = NULL;
	useInputBLs = ubl;
	moveAllNodes = allnm;
	calibNds = clb;
	numCalibNds = (int)calibNds.size();
	numAncFossilsk = 0;
	datedTips = tdt;
	treeScale = rth;
    originTime = iot;
	treeTimePrior = modelPtr->getTreeTimePriorNum();
    conditionOnOrigin = modelPtr->getOriginCondition(); // TAH: some redundancy for now
	name = "TR";
	randShufNdMv = rndNods;
	softBounds = sb;
	isCalibTree = false;
	isTipCals = false; 
	expHyperPrCal = exhpc;
	buildTreeFromNewickDescription(ts); 
	
	nodeProposal = 2; // proposal type 1=window, 2=scale, 3=slide
	tuningVal = log(8.0);
	//cout << "input root height " << rth << endl;
	if(rth < 5.0)
		nodeProposal = 1;
	
	

	if(numCalibNds > 0 && treeTimePrior < 6){
		isCalibTree = true;
		if(datedTips.size() > 0){
			setTipDateAges();
		}
		setNodeCalibrationPriors(ec);
		initializeCalibratedNodeDepths();
		while(checkTreeForCalibrationCompatibility() > 0){
			zeroNodeRedFlags();
			initializeCalibratedNodeDepths();
		}
	}
	else if(treeTimePrior > 5){ 
		isCalibTree = true;
		setUPTGSCalibrationFossils();
		initializeCalibratedNodeDepths();
		while(checkTreeForCalibrationCompatibility() > 0){
			zeroNodeRedFlags();
			initializeCalibratedNodeDepths();
		}
		initializeFossilSpecVariables();

	}
	else{
		isCalibTree = false;
		if(useInputBLs)
			initializeTipNodeDepthsFromUserBL(); 
		else if(isTipCals)
			initializeTipNodeDepthsFromUserBL();
		else
			initializeNodeDepths();
	}
	setAllNodeBranchTimes();
	
	if(isCalibTree)
		calibNodes = getListOfCalibratedNodes();
	
}

Tree::~Tree(void) {

	delete [] nodes;
	delete [] downPassSequence;
}

void Tree::buildTreeFromNewickDescription(string ts) {

	numTaxa = alignmentPtr->getNumTaxa();
	numNodes = 2 * numTaxa - 1;
	nodes = new Node[numNodes];
	downPassSequence = new Node*[numNodes];
	for (int i=0; i<numNodes; i++)
		nodes[i].setIdx(i);

	bool readingBrlen = false;
	Node *p = NULL, *q = NULL;
	int nextInteriorNode = numTaxa;
	for (unsigned i=0; i<ts.size(); i++){
		char c = ts[i];
		if ( c == '(' ){
			q = &nodes[nextInteriorNode++];
			if (p == NULL){
				p = q;
				root = p;
			}
			else{
				q->setAnc(p);
				if (p->getLft() == NULL)
					p->setLft(q);
				else if (p->getRht() == NULL)
					p->setRht(q);
				else{
					cerr << "ERROR: Problem adding interior node to tree" << endl;
					exit(1);
				}
			}
			p = q;
			readingBrlen = false;
		}
		else if ( c == ')' ){
			if (p->getAnc() == NULL){
				cerr << "ERROR: Problem moving down the tree" << endl;
				exit(1);
			}
			else
				p = p->getAnc();
			readingBrlen = false;
		}
		else if ( c == ',' ){
			if (p->getAnc() == NULL){
				cerr << "ERROR: Problem moving down the tree" << endl;
				exit(1);
			}
			else
				p = p->getAnc();
			readingBrlen = false;
		}
		else if ( c == ':' ){
			readingBrlen = true;
		}
		else if ( c == ';' ){
			// we are finished with the tree...check
		}
		else{
			string s = "";
			while ( isValidChar(ts[i]) ){
//				cout << ts[i] << " ";
				s += ts[i++];
			}
			i--;
			if (readingBrlen == false){
				if ( alignmentPtr->isTaxonPresent(s) == false ){
					cerr << "ERROR: Cannot find taxon in alignment" << endl;
					exit(1);
				}
				int indexForTaxon = alignmentPtr->getIndexForTaxonNamed(s);
				q = &nodes[indexForTaxon];
				if (p == NULL){
					cerr << "ERROR: Problem adding a tip to the tree" << endl;
					exit(1);
				}
				else{
					q->setAnc(p);
					if (p->getLft() == NULL)
						p->setLft(q);
					else if (p->getRht() == NULL)
						p->setRht(q);
					else{
						cerr << "ERROR: Problem adding interior node to tree" << endl;
						exit(1);
					}
				}
				p = q;
				p->setName(s);
				p->setIsLeaf(true);
			}
			else{
				double x;
				istringstream buf(s);
				buf >> x;
				p->setUerBL(x);
			}
		}
	}
		
	getDownPassSequence();
	root->setRtGrpVal(0.5); 
	
}

Tree& Tree::operator=(const Tree &t) {

	if (this != &t)
		clone(t);
	return *this;
}

void Tree::clone(const Tree &t) {

	if (numNodes != t.numNodes || numTaxa != t.numTaxa){
		cerr << "ERROR: Attempting to clone trees of unequal size" << endl;
		exit(1);
	}
		
	for (int i=0; i<numNodes; i++){
		Node *pTo   = &nodes[i];
		Node *pFrom = &t.nodes[i];
				
		pTo->setName( pFrom->getName() );
		pTo->setActiveCl( pFrom->getActiveCl() );
		pTo->setActiveTi( pFrom->getActiveTi() );
		pTo->setIsClDirty( pFrom->getIsClDirty() );
		pTo->setIsTiDirty( pFrom->getIsTiDirty() );
		pTo->setNodeDepth( pFrom->getNodeDepth() );
		pTo->setIsLeaf( pFrom->getIsLeaf() );
		pTo->setRtGrpVal( pFrom->getRateGVal() );
		pTo->setBranchTime( pFrom->getBranchTime() ); 
		pTo->setIsCalibratedDepth( pFrom->getIsCalibratedDepth() );
		pTo->setNodeYngTime( pFrom->getNodeYngTime() );
		pTo->setNodeOldTime( pFrom->getNodeOldTime() );
		pTo->setRtGrpIdx( pFrom->getRateGrpIdx() );
		pTo->setNodeCalibPrDist( pFrom->getNodeCalibPrDist() );
		pTo->setNumDecendantTax( pFrom->getNumDecendantTax() );
		pTo->setNodeExpCalRate( pFrom->getNodeExpCalRate() );
		pTo->setNodeAge( pFrom->getNodeAge() );
		pTo->setIsContaminatedFossil( pFrom->getIsContaminatedFossil() );

		pTo->setFossAttchTime( pFrom->getFossAttchTime() );
		pTo->setNumFossAttchLins( pFrom->getNumFossAttchLins() );
		
        // add originTime at some point
		
		pTo->setNumFCalibratingFossils( pFrom->getNumCalibratingFossils() );
		
		
		if (pFrom->getLft() == NULL)
			pTo->setLft(NULL);
		else
			pTo->setLft( &nodes[pFrom->getLft()->getIdx()] );
		
		if (pFrom->getRht() == NULL)
			pTo->setRht(NULL);
		else
			pTo->setRht( &nodes[pFrom->getRht()->getIdx()] );

		if (pFrom->getAnc() == NULL)
			pTo->setAnc(NULL);
		else
			pTo->setAnc( &nodes[pFrom->getAnc()->getIdx()] );
	}
	
	for(int i=0; i<fossSpecimens.size(); i++){
		Fossil *fTo = fossSpecimens[i];
		Fossil *fFrom = t.fossSpecimens[i];
		fTo->setFossilIndex(fFrom->getFossilIndex());
		fTo->setFossilAge(fFrom->getFossilAge());
		fTo->setFossilSppTime(fFrom->getFossilSppTime());
		fTo->setFossilMRCANodeID(fFrom->getFossilMRCANodeID());
		fTo->setFossilMRCANodeAge(fFrom->getFossilMRCANodeAge());
		fTo->setFossilFossBrGamma(fFrom->getFossilFossBrGamma());
		fTo->setFossilIndicatorVar(fFrom->getFossilIndicatorVar());
	}
	numAncFossilsk = t.numAncFossilsk;
	root = &nodes[t.root->getIdx()];
	treeScale = t.treeScale;
	treeTimePrior = t.treeTimePrior;
    
    // TAH: double check this (more bookkeeping is probably needed, for now this is a placeholder)
    OriginTime *ot = modelPtr->getActiveOriginTime();
    originTime = ot->getOriginTime();
		
	for (int i=0; i<numNodes; i++)
		downPassSequence[i] = &nodes[ t.downPassSequence[i]->getIdx() ];
}

int Tree::dex(const Node *p) {
	return (p == NULL ? -1 : p->getIdx());
}

void Tree::getDownPassSequence(void) {

	int x = 0;
	passDown(root, &x);
}

void Tree::setTipDateAges(void){
	
	numExtinctTips = datedTips.size();
	for(vector<Calibration *>::iterator v = datedTips.begin(); v != datedTips.end(); v++){
		int calbTip = -1;
		calbTip = findTip((*v)->getTxN1());
		Node *p = &nodes[calbTip];
				
		p->setNodeYngTime((*v)->getYngTime());
		p->setNodeOldTime((*v)->getOldTime());
		double test = (*v)->getYngTime();
		cout << test << endl;
		p->setNodeAge((*v)->getYngTime());
		(*v)->setNodeIndex(calbTip);
		p->setIsCalibratedDepth(true);
		p->setNodeDepth(p->getNodeAge() / treeScale);
	}
}

int Tree::findTip(string tn){
	
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsLeaf()){
			string n = p->getName();
			if(n == tn){
				return i;
			}
		}
	}
	return -1;
}

void Tree::initializeNodeDepthsFromUserBL(void) {
	
	double inDepth = 0.0;
	int k = 0;
	Node *p = &nodes[k];
	while(!p->getIsLeaf())
		p = &nodes[k++];
	while(p != root){
		inDepth += p->getUerBL();
		p = p->getAnc();
	}
	double scaler = 1.0 / inDepth;
	Node *q = NULL;
	double nodeDepth = 0.0;
	for(int i=0; i<numNodes; i++){
		p = &nodes[i];
		if(p->getIsLeaf())
			p->setNodeDepth(0.0);
		else if(p == root) 
			p->setNodeDepth(1.0);
		else {
			nodeDepth = p->getUerBL();
			q = p->getAnc();
			while(q != root){
				nodeDepth += q->getUerBL();
				q = q->getAnc();
			}
			p->setNodeDepth(1.0 - (nodeDepth * scaler));
		}
	}
	
}

void Tree::initializeTipNodeDepthsFromUserBL(void) {
	
	double inDepth = getUBLTreeScaleDepths();
	Node *p = NULL;
	double scaler = 1.0 / inDepth;
	Node *q = NULL;
	double nodeDepth = 0.0;
	for(int i=0; i<numNodes; i++){
		p = &nodes[i];
		if(p->getIsLeaf())
			p->setNodeDepth(p->getNodeAge() / inDepth);
		else if(p == root) 
			p->setNodeDepth(1.0);
		else {
			nodeDepth = p->getUerBL();
			q = p->getAnc();
			while(q != root){
				nodeDepth += q->getUerBL();
				q = q->getAnc();
			}
			p->setNodeDepth(1.0 - (nodeDepth * scaler));
		}
	}
	
}

double Tree::getUBLTreeScaleDepths(void){
	
	double rootDepth = 0.0;
	Node *p = NULL;
	double *ds = new double[numNodes];
	for(int i=0; i<numNodes; i++){
		p = &nodes[i];
		double inDepth = 0.0;
		if(p != root){
			inDepth = getNodePathDepth(p);
			if(inDepth > rootDepth)
				rootDepth = inDepth;
			cout << "node: " << i << " -- " << inDepth << endl;
		}
		ds[i] = inDepth;
	}
	for(int i=0; i<numNodes; i++){
		p = &nodes[i];
		p->setNodeAge(rootDepth - ds[i]);
		cout << "node: " << i << " -- " << p->getNodeAge() << endl;
	}
	return rootDepth;

}

double Tree::getNodePathDepth(Node *t){
	
	double inDepth = t->getUerBL();;
	Node *p = t->getAnc();
	while(p != root){
		inDepth += p->getUerBL();
		p = p->getAnc();
	}
	return inDepth;
	
}


void Tree::initializeNodeDepths(void) {

	vector<Node*> potentialSplits;
	vector<double> nodeTimes;
	
	for (int i=0; i<numTaxa-2; i++)
		nodeTimes.push_back( ranPtr->uniformRv() );
	sort( nodeTimes.begin(), nodeTimes.end(), greater<double>() );

	for (unsigned i=0; i<nodeTimes.size(); i++)
		cout << nodeTimes[i] << " ";
	cout << endl;
	root->setNodeDepth(1.0);
	
	if (root->getLft()->getIsLeaf() == false)
		potentialSplits.push_back(root->getLft());
	if (root->getRht()->getIsLeaf() == false)
		potentialSplits.push_back(root->getRht());
	
	int nextTime = 0;
	while (potentialSplits.size() > 0){
		Node *p = potentialSplits[(int)(ranPtr->uniformRv()*potentialSplits.size())];
		p->setNodeDepth(nodeTimes[nextTime++]);

		for (vector<Node *>::iterator n=potentialSplits.begin(); n != potentialSplits.end(); n++){
			if ( (*n) == p ){
				potentialSplits.erase( n );
				break;
			}
		}
			
		if (p->getLft()->getIsLeaf() == false)
			potentialSplits.push_back(p->getLft());
		if (p->getRht()->getIsLeaf() == false)
			potentialSplits.push_back(p->getRht());
	}
}

void Tree::initializeCalibratedNodeDepths(void) {
	root->setNodeDepth(1.0);
	bool goodinit = false;
	while(!goodinit){
		for(vector<Calibration *>::iterator v = calibNds.begin(); v != calibNds.end(); v++){
			Node *p = &nodes[(*v)->getNodeIndex()];
			if(p != root){
				Node *anc, *ldes, *rdes;
				anc = p->getAnc();
				ldes = p->getLft();
				rdes = p->getRht();
				double ltime = ldes->getNodeDepth();
				double rtime = rdes->getNodeDepth();
				double atime = anc->getNodeDepth();
				double oldbound = 0.0;
				double yngbound = p->getNodeYngTime() / treeScale;
				if(p->getNodeCalibPrDist() == 1)
					oldbound = p->getNodeOldTime() / treeScale;
				else
					oldbound = getTemporaryNodeMaxBound(p) / treeScale;
				if(p->getNodeOldTime() == p->getNodeYngTime()){
					double newNodeDepth = yngbound;
					p->setNodeDepth(newNodeDepth);
					goodinit = true;
				}
				else{
					if(ltime > 0.0 || rtime > 0.0){
						if(rtime > ltime){
							if(rtime > yngbound)
								yngbound = rtime;
						}
						else{
							if(ltime > yngbound)
								yngbound = ltime;
						}
					}
					if(atime > 0.0){
						if(atime < oldbound)
							oldbound = atime;
					}
					if(oldbound > root->getNodeDepth())
						oldbound = root->getNodeDepth();
					if(yngbound > oldbound){
						goodinit = false;
						break;
					}
					else{
						double newNodeDepth = yngbound + ranPtr->uniformRv()*(oldbound-yngbound);
						p->setNodeDepth(newNodeDepth);
						goodinit = true;
					}
				}
			}
			else
				goodinit = true;
		}
	}
	setNodesNumberDecendantTaxa(root);
	int nnc = 0;
	recursiveNodeDepthInitialization(root, nnc, 1.0);
}


void Tree::checkNodeInit(){
	
	for(int i=0; i<fossSpecimens.size(); i++){
		Fossil *f = fossSpecimens[i];
		int nID = f->getFossilMRCANodeID();
		Node *p = &nodes[nID];
		double nage = p->getNodeDepth() * treeScale;
		double fage = f->getFossilAge();
		cout << nID << " -- " << nage - fage << " --- " << nage << " --- " << fage << endl;
	}
}

double Tree::getTemporaryNodeMaxBound(Node *p){
	
	double tempmax = treeScale;
	Node *q = p->getAnc();
	while(q != root && tempmax == treeScale){
		if(q->getIsCalibratedDepth())
			tempmax = q->getNodeYngTime();
		else
			q = q->getAnc();
	}
	return tempmax;
}


vector<double> Tree::recursiveNodeDepthInitialization(Node *p, int &nCont, double maxD) {
	
	vector<double> ndTimesC;
	if(p->getIsLeaf()){
		for(int i=0; i< nCont; i++)
			ndTimesC.push_back( ranPtr->uniformRv()*(maxD) );
		if(ndTimesC.size() > 1)
			sort( ndTimesC.begin(), ndTimesC.end() );
		return ndTimesC;
	}
	else{
		bool getDFrmVec = false;
		int nAncN = nCont;
		double myDepth = maxD;
		if(p->getIsCalibratedDepth()){
			myDepth = p->getNodeDepth();
			nCont = 0;
			getDFrmVec = false;
		}
		else{
			nCont += 1;
			getDFrmVec = true;
		}
	
		Node *d1 = p->getLft();
		Node *d2 = p->getRht();
		if(p->getLft()->getNumDecendantTax() < p->getRht()->getNumDecendantTax()){
			if(p->getRht()->getIsCalibratedDepth() == false && p->getLft()->getIsCalibratedDepth() == true){
				d1 = p->getLft();
				d2 = p->getRht();
			}
			else{
				d2 = p->getLft();
				d1 = p->getRht();
			}
		}
		else{
			if(p->getRht()->getIsCalibratedDepth() == true && p->getLft()->getIsCalibratedDepth() == false){
				d2 = p->getLft();
				d1 = p->getRht();
			}
			else{
				d1 = p->getLft();
				d2 = p->getRht();
			}
		}
		ndTimesC = recursiveNodeDepthInitialization(d1, nCont, myDepth);
		if(getDFrmVec){
			if(p->getNodeDepth() > 0.0){
				myDepth = p->getNodeDepth();
				if(ndTimesC.size() > 0)
					ndTimesC.erase(ndTimesC.begin()+0);
			}
			else{
				if(ndTimesC.size() == 0){
					int myID = p->getIdx();
					cerr << "ERROR: node vector is empty at Node " << myID << endl;
					exit(1);
				}
				myDepth = ndTimesC[0];
				p->setNodeDepth(myDepth);
				ndTimesC.erase(ndTimesC.begin()+0);
			}
		}
		nCont = 0;
		recursiveNodeDepthInitialization(d2, nCont, myDepth);
		
		if(getDFrmVec)
			return ndTimesC;
		else{
			if(ndTimesC.size() > 0){
				int myID = p->getIdx();
				cout << "Calibrated " << myID << "  (MAX = " << maxD*treeScale; 
				cout << ")  (MIN = " << myDepth*treeScale << ") ";
				for(vector<double>::iterator v = ndTimesC.begin(); v != ndTimesC.end(); v++){
					cout << "  " << (*v)*treeScale << "  ";
				}
				cout << endl;
				cerr << "ERROR: node vector has stuff in it" << endl;
				exit(1);
			}
			if(myDepth > maxD){

				double badDep = maxD;
				vector<int> idNFix;
				Node *ancs = p->getAnc();
				while(badDep < myDepth){
					badDep = ancs->getNodeDepth();
					if(badDep < myDepth)
						idNFix.push_back(ancs->getIdx());
					ancs = ancs->getAnc();
				}
				for(int i=0; i<idNFix.size(); i++)
					ndTimesC.push_back( myDepth + ranPtr->uniformRv()*(badDep - myDepth) );
				sort( ndTimesC.begin(), ndTimesC.end() );
				for(int i=0; i<idNFix.size(); i++){
					nodes[idNFix[i]].setNodeDepth(ndTimesC[i]);
				}
				ndTimesC.clear();
				return ndTimesC;
			}
			for(int i=0; i< nAncN; i++)
				ndTimesC.push_back( myDepth + ranPtr->uniformRv()*(maxD - myDepth) );
			sort( ndTimesC.begin(), ndTimesC.end() );
			return ndTimesC;
		}
	}
	return ndTimesC;
}

int Tree::setNodesNumberDecendantTaxa(Node *p) {
	
	if(p->getIsLeaf()){
		p->setNumDecendantTax(0);
		return 1;
	}
	else{
		int nds = setNodesNumberDecendantTaxa(p->getLft());
		nds += setNodesNumberDecendantTaxa(p->getRht());
		p->setNumDecendantTax(nds);
		return nds;
	}
	return 0;
}



bool Tree::isValidChar(char c) {

	if ( c == '(' || c == ')' || c == ',' || c == ':' || c == ';' )
		return false;
	return true;
}

void Tree::passDown(Node *p, int *x) {

    
	if ( p != NULL )
		{
		passDown(p->getLft(), x);
		passDown(p->getRht(), x);
		downPassSequence[(*x)++] = p;
		}
}

void Tree::print(std::ostream & o) const {

	o << "Tree:\n";
	showNodes(root, 3, o);
}

double Tree::update(double &oldLnL) {
	
    OriginTime *ot = modelPtr->getActiveOriginTime();
    originTime = ot->getOriginTime();

	double lppr = 0.0;
	double probCalMove = 0.5;
	if(treeTimePrior == 6){
		recountFossilAttachNums();
		if(ranPtr->uniformRv() > probCalMove){
			if(moveAllNodes){ 
				if(randShufNdMv)
					lppr = updateAllNodesRnd(oldLnL);
				else lppr = updateAllNodes(oldLnL);
			}	
			else lppr = updateOneNode();
			setAllNodeBranchTimes();
		}
		else{
			updateFossilBDSSAttachmentTimePhi();
			modelPtr->setLnLGood(true);
			modelPtr->setMyCurrLnl(oldLnL);
			Tree *t = modelPtr->getActiveTree();
			t->upDateAllCls();
			t->upDateAllTis();
			modelPtr->setTiProb();
			return 0.0;
		}
		recountFossilAttachNums();
	}
	else if(treeTimePrior >= 7){
		
		Treescale *ts = modelPtr->getActiveTreeScale();
		setTreeScale(ts->getScaleValue());

		updateAllTGSNodes(oldLnL);
		
		updateFossilBDSSAttachmentTimePhi();
        treeUpdateNodeOldestBoundsAttchTimes();
		modelPtr->setLnLGood(true);
		modelPtr->setMyCurrLnl(oldLnL);
		modelPtr->setTiProb();

		for(int i=0; i<5; i++){
			updateRJMoveAddDelEdge();
            treeUpdateNodeOldestBoundsAttchTimes();
			modelPtr->setLnLGood(true);
			modelPtr->setMyCurrLnl(oldLnL);
			modelPtr->setTiProb();
		}
		return 0.0;
	}
	else{ 
		if(moveAllNodes){ 
			updateAllNodes(oldLnL);
		}	
		else lppr = updateOneNode();
		setAllNodeBranchTimes();
	}
	return lppr;
}


double Tree::updateFossilBDSSAttachmentTimePhi() {
	
	Speciation *s = modelPtr->getActiveSpeciation();
	s->setAllBDFossParams();
	double lambda = s->getBDSSSpeciationRateLambda();	
	double mu = s->getBDSSExtinctionRateMu();
	double fossRate = s->getBDSSFossilSampRatePsi();
	double sppSampRate = s->getBDSSSppSampRateRho();
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


			double lnPriorRat = 0.0;
			double v1 = newSumLogGammas;
			double v2 = bdssQFxn(lambda, mu, fossRate, sppSampRate, newPhi);
			double v3 = bdssQFxn(lambda, mu, fossRate, sppSampRate, oldPhi);
			double v4 = oldSumLogGammas;
			lnPriorRat = (v1 - v2) - (v4 - v3);
			
			
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

double Tree::updateRJMoveAddDelEdge() {
	
	double gA = 0.5; 
	if(numAncFossilsk == numCalibNds)
		gA = 1.0;
	else if(numAncFossilsk == 0)
		gA = 0.0;
	
	double u = ranPtr->uniformRv();
	if(u < gA)
		doAddEdgeMove();
	else
		doDeleteEdgeMove();
	
	return 0.0;
}

void Tree::doAddEdgeMove(){
	
	int k = numAncFossilsk;
	int m = numCalibNds;
	int kN = k-1;
	Speciation *s = modelPtr->getActiveSpeciation();
	s->setAllBDFossParams();
	double lambda = s->getBDSSSpeciationRateLambda();	
	double mu = s->getBDSSExtinctionRateMu();
	double fossRate = s->getBDSSFossilSampRatePsi();
	double sppSampRate = s->getBDSSSppSampRateRho();
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
	double newPhi = yf + nu;
	double scnewPhi = newPhi / treeScale; 
	f->setFossilSppTime(scnewPhi);
	f->setFossilIndicatorVar(1);
		
	double newSumLogGammas = getSumLogAllAttachNums(); 
		
	double v1 = log(bdssP0Fxn(lambda, mu, fossRate, sppSampRate, yf)); 
	double v2 = bdssQFxn(lambda, mu, fossRate, sppSampRate, yf); 
	double v3 = bdssQFxn(lambda, mu, fossRate, sppSampRate, newPhi); 
	
	lnPriorR = (numTaxa - 2 + numCalibNds - kN) * log(lambda);
	lnPriorR += (newSumLogGammas + log(2.0) + v1 + v2 - v3);
	lnPriorR -= (((numTaxa - 2 + numCalibNds - k) * log(lambda)) + oldSumLogGammas);
	
	double lpr = lnPriorR + lnHastings + lnJacobian;
	double r = modelPtr->safeExponentiation(lpr);
	
	if(ranPtr->uniformRv() < r){ 
		f->setFossilIndicatorVar(1);
		f->setFossilSppTime(scnewPhi);
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

void Tree::doDeleteEdgeMove(){

	
	int k = numAncFossilsk;
	int m = numCalibNds;
	int kN = k+1;
	Speciation *s = modelPtr->getActiveSpeciation();
	s->setAllBDFossParams();
	double lambda = s->getBDSSSpeciationRateLambda();	
	double mu = s->getBDSSExtinctionRateMu();
	double fossRate = s->getBDSSFossilSampRatePsi();
	double sppSampRate = s->getBDSSSppSampRateRho();
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
	
	double oldPhi = f->getFossilSppTime() * treeScale;
	double oldScPhi = f->getFossilSppTime();
	double scnewPhi = yf / treeScale;
	f->setFossilIndicatorVar(0);
	f->setFossilSppTime(scnewPhi);
	
	double newSumLogGammas = getSumLogAllAttachNums(); 

	lnHastings = log(alD) + (log(m-k) - log(k+1));
	lnJacobian = -(log(cf - yf));

	double v1 = log(bdssP0Fxn(lambda, mu, fossRate, sppSampRate, yf)); 
	double v2 = bdssQFxn(lambda, mu, fossRate, sppSampRate, yf); 
	double v3 = bdssQFxn(lambda, mu, fossRate, sppSampRate, oldPhi); 
		
	lnPriorR = ((numTaxa - 2 + numCalibNds - kN) * log(lambda)) + newSumLogGammas;
	lnPriorR -= (((numTaxa - 2 + numCalibNds - k) * log(lambda)) + (oldSumLogGammas + log(2.0) + v1 + v2 - v3));
	
	double lpr = lnPriorR + lnHastings + lnJacobian;
	double r = modelPtr->safeExponentiation(lpr);

	if(ranPtr->uniformRv() < r){ 
		f->setFossilIndicatorVar(0);
		f->setFossilSppTime(scnewPhi);
		setNodeOldestAttchBranchTime(p);
		numAncFossilsk = kN;
	}
	else{
		f->setFossilSppTime(oldScPhi);
		f->setFossilIndicatorVar(1);
		setNodeOldestAttchBranchTime(p);
	}
}

double Tree::updateOneNode() {
	
	Node *p = NULL;
	do {
		p = downPassSequence[(int)(ranPtr->uniformRv()*numNodes)];
	} while (p->getLft() == NULL || p->getRht() == NULL || p->getAnc() == NULL);
	
	double largestTime = getNodeUpperBoundTime(p);
	double smallestTime = getNodeLowerBoundTime(p);
	double currDepth = p->getNodeDepth();
	double newNodeDepth = currDepth;
	double lnPrRatio = 0.0;
	if(largestTime > smallestTime){
		newNodeDepth = smallestTime + ranPtr->uniformRv()*(largestTime-smallestTime);
		
		if(treeTimePrior > 5)
			lnPrRatio = lnPriorRatioTGS(newNodeDepth, currDepth, p);
		else lnPrRatio = lnPriorRatio(newNodeDepth, currDepth);

		p->setNodeDepth(newNodeDepth);
		
		flipToRootClsTis(p);
		updateToRootClsTis(p);
		modelPtr->setTiProb();
		if(p->getIsCalibratedDepth()){
			if(softBounds && p->getNodeCalibPrDist() == 1){
				double ycal = p->getNodeYngTime();
				double ocal = p->getNodeOldTime();
				lnPrRatio += lnCalibPriorRatio(newNodeDepth*treeScale, currDepth*treeScale, ycal, ocal);
			}
			if(p->getNodeCalibPrDist() == 2){
				double offst = p->getNodeYngTime();
				double nodeExCR = p->getNodeExpCalRate();
				lnPrRatio += lnExpCalibPriorRatio(newNodeDepth*treeScale, currDepth*treeScale, offst, nodeExCR);
			}
		}
		
	}
	if(treeTimePrior > 5)
		recountFossilAttachNums();
	return lnPrRatio;
}

double Tree::updateAllNodes(double &oldLnL) {
	
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
				
				double lnPrRatio = lnPriorRatio(newNodeDepth/treeScale, currDepth/treeScale);
				p->setNodeDepth(newNodeDepth/treeScale);
				
				flipToRootClsTis(p);
				updateToRootClsTis(p);
				modelPtr->setTiProb();
				
				double newLnl = modelPtr->lnLikelihood();
				double lnLRatio = newLnl - oldLike;
				if(p->getIsCalibratedDepth()){					
					if(softBounds && p->getNodeCalibPrDist() == 1){
						double ycal = p->getNodeYngTime();
						double ocal = p->getNodeOldTime();
						lnPrRatio += lnCalibPriorRatio(newNodeDepth, currDepth, ycal, ocal);
					}
					if(p->getNodeCalibPrDist() == 2){
						double offst = p->getNodeYngTime();
						double nodeExCR = p->getNodeExpCalRate();
						lnPrRatio += lnExpCalibPriorRatio(newNodeDepth, currDepth, offst, nodeExCR);
					}
				}
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
		}
	}
	oldLnL = oldLike;
	return 0.0;
}

double Tree::updateAllTGSNodes(double &oldLnL) {
		
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
				
				double lnPrRatio = lnPriorRatioTGS(newNodeDepth, currDepth, p);
				p->setNodeDepth(newNodeDepth/treeScale);
	
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

double Tree::doAScaleMove(double &nv, double cv, double tv, double lb, double hb, double rv){
	
	double c = tv * (rv - 0.5);
	double newcv = cv * exp(c);
	bool validV = false;
	do{
		if(newcv < lb)
			newcv = lb * lb / newcv;
		else if(newcv > hb)
			newcv = hb * hb / newcv;
		else
			validV = true;
	} while(!validV);
	nv = newcv;
	return c;
}

double Tree::doAWindoMove(double &nv, double cv, double tv, double lb, double hb, double rv){
	
	double newCV = cv + (tv * (rv - 0.5));
	do{
		if(newCV < lb)
			newCV = 2.0 * lb - newCV;
		else if(newCV > hb)
			newCV = 2.0 * hb - newCV;
	}while(newCV < lb || newCV > hb);
	nv = newCV;
	return 0.0;
}


double Tree::updateAllNodesRnd(double &oldLnL) {
	
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
			double currDepth = p->getNodeDepth();
			double largestTime = getNodeUpperBoundTime(p);
			double smallestTime = getNodeLowerBoundTime(p);

			if (largestTime > smallestTime){
				double newNodeDepth = smallestTime + ranPtr->uniformRv()*(largestTime-smallestTime);
				
				double lnPrRatio = 0.0;
				if(treeTimePrior > 5)
					lnPrRatio = lnPriorRatioTGS(newNodeDepth, currDepth, p);
				else lnPrRatio = lnPriorRatio(newNodeDepth, currDepth);

				p->setNodeDepth(newNodeDepth);
				
				flipToRootClsTis(p);
				updateToRootClsTis(p);
				modelPtr->setTiProb();
				
				double newLnl = modelPtr->lnLikelihood();
				double lnLRatio = newLnl - oldLike;
				if(p->getIsCalibratedDepth()){
					if(softBounds && p->getNodeCalibPrDist() == 1){
						double ycal = p->getNodeYngTime();
						double ocal = p->getNodeOldTime();
						lnPrRatio += lnCalibPriorRatio(newNodeDepth*treeScale, currDepth*treeScale, ycal, ocal);
					}
					if(p->getNodeCalibPrDist() == 2){
						double offst = p->getNodeYngTime();
						double nodeExCR = p->getNodeExpCalRate();
						lnPrRatio += lnExpCalibPriorRatio(newNodeDepth*treeScale, currDepth*treeScale, offst, nodeExCR);
					}
				}
				
				double lnR = lnPrRatio + lnLRatio + 0.0;
				double r = modelPtr->safeExponentiation(lnR);
				
				if(ranPtr->uniformRv() < r){
					oldLike = newLnl;
				}
				else{
					p->setNodeDepth(currDepth);
					flipToRootClsTis(p);
					updateToRootClsTis(p);
				}
			}
		}
		else if(p->getIsLeaf() && p->getIsCalibratedDepth() == false)
			p->setNodeDepth(0.0);
		if(treeTimePrior > 5)
			recountFossilAttachNums();
	}
	oldLnL = oldLike;	
	return 0.0;
}


void Tree::setAllNodeBranchTimes(void) {
	
	for (int n=0; n<numNodes; n++){
		Node *p = &nodes[n];
		if(p != root){
			if(p->getIsLeaf() && p->getIsCalibratedDepth() == false)
				p->setNodeDepth(0.0);
			double branchtime = (p->getAnc()->getNodeDepth() - p->getNodeDepth()) * treeScale;
			if(branchtime < 0){
				int myID = p->getIdx();
				cerr << "ERROR: The tree has a negative branch length! " << branchtime << " At Node " << myID 
				<< ", With depth = " << p->getNodeDepth() * treeScale << endl;
				exit(1);
			}
			else
				p->setBranchTime(branchtime);
		}
	}
}

double Tree::lnPrior() {
	
	return 0.0;
}

double Tree::lnPriorRatio(double snh, double soh) {
		
	double nh = snh*treeScale;
	double oh = soh*treeScale;

	
	if(treeTimePrior == 1) 
		return 0.0;
	else if(treeTimePrior == 2){ 
		double diff = modelPtr->getActiveSpeciation()->getNetDiversification();		
		double nator = (-(diff)*nh);
		double dator = (-(diff)*oh);
		
		return nator - dator;
	}
	else if(treeTimePrior == 3){
		Speciation *s = modelPtr->getActiveSpeciation();
		double diff = s->getNetDiversification();	
		double rel = s->getRelativeDeath();			
		
		double zn = log(1 - (rel) * exp(-(diff)*nh));
		double nator = -2 * zn + (-(diff)*nh);
		
		double zd = log(1 - (rel) * exp(-(diff)*oh));
		double dator = -2 * zd + (-(diff)*oh);
		
		return nator - dator;
	}
	else if(treeTimePrior == 4 || treeTimePrior == 5){ 
		Speciation *s = modelPtr->getActiveSpeciation();
		double netDiv = s->getNetDiversification();	
		double rel = s->getRelativeDeath();			
		double lambda = netDiv / (1.0 - rel); 
		double mu = rel * lambda;
		double fossRate = s->getBDSSFossilSampRatePsi();
		double sppSampRate = s->getBDSSSppSampRateRho();
		
		double topV = -(bdssQFxn(lambda,mu,fossRate,sppSampRate,nh));
		double botV = -(bdssQFxn(lambda,mu,fossRate,sppSampRate,oh));
		return topV - botV;
	}
	return 0.0;
}

double Tree::lnPriorRatioTGS(double snh, double soh, Node *p) {

	double nh = snh; 
	double oh = soh; 
	
	p->setNodeDepth(soh/treeScale);
	double oldSumLogGammas = getSumLogAllAttachNums(); 
	
	p->setNodeDepth(snh/treeScale);
	double newSumLogGammas = getSumLogAllAttachNums(); 

	Speciation *s = modelPtr->getActiveSpeciation();
	s->setAllBDFossParams();
	double lambda = s->getBDSSSpeciationRateLambda();	
	double mu = s->getBDSSExtinctionRateMu();			
	double fossRate = s->getBDSSFossilSampRatePsi();
	double sppSampRate = s->getBDSSSppSampRateRho();
	
	double topV = -(bdssQFxn(lambda,mu,fossRate,sppSampRate,nh)) + newSumLogGammas;
	double botV = -(bdssQFxn(lambda,mu,fossRate,sppSampRate,oh)) + oldSumLogGammas;
	double v = topV - botV;
	return v;
		
}


double Tree::bdssC1Fxn(double b, double d, double psi){
	
	double v = abs( sqrt( ( (b-d-psi) * (b-d-psi) ) + 4*b*psi) );
	return v;
}

double Tree::bdssC2Fxn(double b, double d, double psi, double rho){
	
	double v = -( ( b-d-(2*b*rho)-psi ) / (bdssC1Fxn(b,d,psi)) );
	return v;
}

double Tree::bdssQFxn(double b, double d, double psi, double rho, double t){
	
	double c1Val = bdssC1Fxn(b,d,psi);
	double c2Val = bdssC2Fxn(b,d,psi,rho);
	
	double vX = c1Val * t + 2.0 * log(exp(-c1Val * t) * (1.0 - c2Val) + (1.0 + c2Val));
	
	// returns log(q)
	return vX;
}

double Tree::bdssP0Fxn(double b, double d, double psi, double rho, double t){
	
	double c1Val = bdssC1Fxn(b,d,psi);
	double c2Val = bdssC2Fxn(b,d,psi,rho);
	
	double eCfrac = (exp(-c1Val * t) * (1.0 - c2Val) - (1.0 + c2Val)) / (exp(-c1Val * t) * (1.0 - c2Val) + (1.0 + c2Val));
	double v = 1.0 + ((-(b - d - psi)) + (c1Val * eCfrac)) / (2.0 * b);
	
	return v;
}

double Tree::bdssP0HatFxn(double b, double d, double rho, double t){
	
	double v = 1.0 - ((rho * (b - d)) / ((b*rho) + (((b * (1-rho)) - d) * exp(-(b-d)*t))));
	return v;
}

double Tree::fbdQHatFxn(double b, double d, double psi, double rho, double t){

	double v = log(4.0 * rho);
	v -= bdssQFxn(b,d,psi,rho,t);
	return v;
}


double Tree::lnCalibPriorRatio(double nh, double oh, double lb, double ub) {
	
	/*
		Prior ratio for a uniform dist on calibrated nodes with soft bounds
		From Yang and Rannala, MBE (2006)
	*/
	double numr = 0.0;
	double dnom = 0.0;
	double th1 = (0.95 * lb) / (0.025 * (ub - lb));
	double th2 = 0.95 / (0.025 * (ub - lb));
	
	if(nh <= lb)
		numr += log(0.025) + log(th1 / lb) + ((th1 - 1) * log(nh / lb));
	else if(nh > ub)
		numr += log(0.025) + log(th2) + (-th2 * (nh - ub));
	else
		numr += log(0.95) - log(ub - lb);
	
	if(oh <= lb)
		dnom += log(0.025) + log(th1 / lb) + ((th1 - 1) * log(oh / lb));
	else if(oh > ub)
		dnom += log(0.025) + log(th2) + (-th2 * (oh - ub));
	else
		dnom += log(0.95) - log(ub - lb);
	
	return numr - dnom;
}

double Tree::lnExpCalibPriorRatio(double nh, double oh, double offSt, double expRate) {
	
	/*
	 TAH: Prior ratio for offset exponential
	 */
	double numr = 0.0;
	double dnom = 0.0;

	numr = ranPtr->lnExponentialPdf(expRate, nh - offSt);
	dnom = ranPtr->lnExponentialPdf(expRate, oh - offSt);
	return numr - dnom;
}



void Tree::flipAllCls(void) {

	for (int n=0; n<numNodes; n++)
		nodes[n].flipCl();
}

void Tree::flipAllTis(void) {

	for (int n=0; n<numNodes; n++)
		nodes[n].flipTi();
}

void Tree::upDateAllCls(void) {

	for (int n=0; n<numNodes; n++)
		nodes[n].setIsClDirty(true);
}

void Tree::upDateAllTis(void) {

	for (int n=0; n<numNodes; n++)
		nodes[n].setIsTiDirty(true);
}

void Tree::flipToRootClsTis(Node *p){
	
	Node *q = p;
	if(!q->getIsLeaf()){
		q->getLft()->flipCl();
		q->getRht()->flipCl();
		q->getLft()->flipTi();
		q->getRht()->flipTi();
	}
	while(q != NULL){
		q->flipCl();
		q->flipTi();
		q = q->getAnc();
	}
}

void Tree::updateToRootClsTis(Node *p){
	
	Node *q = p;
	if(!q->getIsLeaf()){
		q->getLft()->setIsClDirty(true);
		q->getRht()->setIsClDirty(true);
		q->getLft()->setIsTiDirty(true);
		q->getRht()->setIsTiDirty(true);
	}
	while(q != NULL){
		q->setIsClDirty(true);
		q->setIsTiDirty(true);
		q = q->getAnc();
	}
}

// overloaded function when passing in just the node ID, then this is called by the DPP rate move.
void Tree::flipToRootClsTis(int ndID){

	Node *q = &nodes[ndID];
	while(q != NULL){
		q->flipCl();
		q->flipTi();
		q = q->getAnc();
	}
}

void Tree::updateToRootClsTis(int ndID){

	Node *q = &nodes[ndID];
	while(q != NULL){
		q->setIsClDirty(true);
		q->setIsTiDirty(true);
		q = q->getAnc();
	}
}


string Tree::getTreeDescription(void){ 
	
	stringstream ss;
	writeTree(root, ss);
	ss << ";";
	string treDesc = ss.str();
	return treDesc;
}

void Tree::writeTree(Node *p, stringstream &ss){
	
	if (p != NULL){
		if(p->getLft() == NULL){ 
			ss << p->getName();
		}
		else{
			ss << "(";
			writeTree(p->getLft(), ss);
			ss << ":" << p->getLft()->getBranchTime();
			ss << ",";
			writeTree(p->getRht(), ss);
			ss << ":" << p->getRht()->getBranchTime();
			ss << ")";
		}
	}
}

string Tree::getFigTreeDescription(void){ 
	
	stringstream ss;
	writeFigTree(root, ss);
	ss << ";";
	string treDesc = ss.str();
	return treDesc;
}

void Tree::writeFigTree(Node *p, stringstream &ss){
	
	if (p != NULL){
		if(p->getLft() == NULL){ 
			ss << p->getName();
		}
		else{
			ss << "(";
			writeFigTree(p->getLft(), ss);
			ss << "[&rate=" << p->getLft()->getRateGVal() << ",";
			ss << "rate_cat=" << p->getLft()->getRateGrpIdx() << "]";
			ss << ":" << p->getLft()->getBranchTime();
			ss << ",";
			writeFigTree(p->getRht(), ss);
			ss << "[&rate=" << p->getRht()->getRateGVal() << ",";
			ss << "rate_cat=" << p->getRht()->getRateGrpIdx() << "]";
			ss << ":" << p->getRht()->getBranchTime();
			ss << ")";
		}
	}
}

string Tree::getCalibInitialTree(void){ 
	
	stringstream ss;
	writeCalibrationFigTree(root, ss);
	ss << "[&rate=0.0,";
	ss << "rate_cat=0";
	if(root->getIsCalibratedDepth())
		ss << ",cal_range={" << root->getNodeYngTime() << "," << root->getNodeOldTime() << "}]";
	else
		ss << ",cal_range={" << treeScale << "," << treeScale << "}]";
	
	ss << ";";
	string treDesc = ss.str();
	return treDesc;
}

void Tree::writeCalibrationFigTree(Node *p, stringstream &ss){
	
	if (p != NULL){
		if(p->getLft() == NULL){ 
			ss << p->getName();
		}
		else{
			ss << "(";
			writeCalibrationFigTree(p->getLft(), ss);
			ss << "[&rate=" << p->getLft()->getRateGVal() << ",";
			ss << "rate_cat=" << p->getLft()->getRateGrpIdx();
			if(p->getLft()->getIsCalibratedDepth())
				ss << ",cal_range={" << p->getLft()->getNodeYngTime() << "," << p->getLft()->getNodeOldTime() << "}]";
			else
				ss << ",cal_range={" << p->getLft()->getNodeDepth() * treeScale << "," << p->getLft()->getNodeDepth() * treeScale << "}]";
			ss << ":" << p->getLft()->getBranchTime();
			ss << ",";
			writeCalibrationFigTree(p->getRht(), ss);
			ss << "[&rate=" << p->getRht()->getRateGVal() << ",";
			ss << "rate_cat=" << p->getRht()->getRateGrpIdx();
			if(p->getRht()->getIsCalibratedDepth())
				ss << ",cal_range={" << p->getRht()->getNodeYngTime() << "," << p->getRht()->getNodeOldTime() << "}]";
			else
				ss << ",cal_range={" << p->getRht()->getNodeDepth() * treeScale << "," << p->getRht()->getNodeDepth() * treeScale << "}]";
			ss << ":" << p->getRht()->getBranchTime();
			ss << ")";
		}
	}
}


string Tree::writeParam(void){
	
	stringstream ss;
	ss << "Tree:" << endl;
	showNodes(root, 3, ss);
	string outp = ss.str();
	return outp;
}

void Tree::showNodes(Node *p, int indent, std::ostream &ss) const {
	
	if (p != NULL){
		for (int i=0; i<indent; i++)
			ss << " ";
		ss << dex(p) << " (" << dex(p->getLft()) << ", " << dex(p->getRht()) << ", " << dex(p->getAnc()) 
		   << ") -- " << fixed << setprecision(5) << (p->getNodeDepth() * treeScale) << " -- ";
		if (p->getLft() == NULL && p->getRht() == NULL )
			ss << " (" << p->getName() << ") ";
		if (p == root)
			ss << " <- Root";
		ss << endl;
		showNodes (p->getLft(),  indent + 2, ss);
		showNodes (p->getRht(), indent + 2, ss);
	}
}

string Tree::getNodeInfoNames(void){
	
	stringstream ss;
	Node *p = NULL;
	for(int i=0; i<numNodes; i++){
		p = &nodes[i];
		if(p->getIsLeaf() == false){
			if(p == root)
				ss << "\tRootDepth.Time(N" << p->getIdx() << ")";
			else
				ss << "\tTime(N" << p->getIdx() << ")";
		}
	}
	for(int i=0; i<numNodes; i++){
		p = &nodes[i];
#		if ASSIGN_ROOT
		ss << "\tSRate(N" << p->getIdx() << ")";
#		else
		if(p != root)
			ss << "\tSRate(N" << p->getIdx() << ")";
#		endif
	}
	string ni = ss.str();
	return ni;
}


string Tree::getNodeInfoList(void){
	
	stringstream ss;
	Node *p = NULL;
	for(int i=0; i<numNodes; i++){
		p = &nodes[i];
		if(p->getIsLeaf() == false)
			ss << "\t" << p->getNodeDepth() * treeScale;
	}
	for(int i=0; i<numNodes; i++){
		p = &nodes[i];
#		if ASSIGN_ROOT
		ss << "\t" << p->getRateGVal();
#		else
		if(p != root)
			ss << "\t" << p->getRateGVal();
#		endif
	}
	string ni = ss.str();
	return ni;
}

string Tree::getDownPNodeInfoNames(void){
	
	stringstream ss;
	Node *p = NULL;
	for(int i=0; i<numNodes; i++){
		p = downPassSequence[i];
		if(p->getIsLeaf() == false){
			if(p == root)
				ss << "\tRootDepth.Time(DP" << i << "|N" << p->getIdx() << ")";
			else
				ss << "\tTime(DP" << i << "|N" << p->getIdx() << ")";
		}
	}
	for(int i=0; i<numNodes; i++){
		p = downPassSequence[i];
#		if ASSIGN_ROOT
		ss << "\tSRate(DP" << i << "|N" << p->getIdx() << ")";
#		else
		if(p != root)
			ss << "\tSRate(DP" << i << "|N" << p->getIdx() << ")";
#		endif
	}
	string ni = ss.str();
	return ni;
}


string Tree::getDownPNodeInfoList(void){
	
	stringstream ss;
	Node *p = NULL;
	for(int i=0; i<numNodes; i++){
		p = downPassSequence[i];
		if(p->getIsLeaf() == false)
			ss << "\t" << p->getNodeDepth() * treeScale;
	}
	for(int i=0; i<numNodes; i++){
		p = downPassSequence[i];
#		if ASSIGN_ROOT
		ss << "\t" << p->getRateGVal();
#		else
		if(p != root)
			ss << "\t" << p->getRateGVal();
#		endif
	}
	string ni = ss.str();
	return ni;
}

string Tree::getCalNodeInfoNames(void){
	
	stringstream ss;
	Node *p = NULL;
	for(int i=0; i<numNodes; i++){
		p = &nodes[i];
		if(p->getIsCalibratedDepth()){
			if(p == root)
				ss << "\tRoot.expHP(N" << p->getIdx() << ")";
			else
				ss << "\texpHP(N" << p->getIdx() << ")";
		}
	}
	string ni = ss.str();
	return ni;
}


string Tree::getCalNodeInfoList(void){
	
	stringstream ss;
	Node *p = NULL;
	for(int i=0; i<numNodes; i++){
		p = &nodes[i];
		if(p->getIsCalibratedDepth())
			ss << "\t" << p->getNodeExpCalRate();
	}
	string ni = ss.str();
	return ni;
}

string Tree::getCalBDSSNodeInfoParamNames(void){
	
	stringstream ss;
	Fossil *f = NULL;
//	for(vector<Node *>::iterator v = calibNodes.begin(); v != calibNodes.end(); v++){
//		int idx = (*v)->getIdx();
//		ss << "\tcalib.time(N" << idx << ")";
//	}
//	for(int i=0; i<fossSpecimens.size(); i++){
//		f = fossSpecimens[i];
//		int nID = f->getFossilMRCANodeID();
//		ss << "\tcal.dist(C" << i << ".nd" << nID << ")";
//	}
	for(int i=0; i<fossSpecimens.size(); i++){
		f = fossSpecimens[i];
		int nID = f->getFossilMRCANodeID();
		ss << "\tphi(C" << i << ".nd" << nID << ")";
	}
	for(int i=0; i<fossSpecimens.size(); i++){
		f = fossSpecimens[i];
		int nID = f->getFossilMRCANodeID();
		ss << "\tgamma(C" << i << ".nd" << nID << ")";
	}
	string ni = ss.str();
	return ni;
}


string Tree::getCalBDSSNodeInfoParamList(void){
	
	stringstream ss;
	Fossil *f = NULL;
//	for(vector<Node *>::iterator v = calibNodes.begin(); v != calibNodes.end(); v++){
//		setNodeOldestAttchBranchTime((*v));
//		double t = (*v)->getFossAttchTime() * treeScale;
//		ss << "\t" << t;
//	}
//	for(int i=0; i<fossSpecimens.size(); i++){
//		f = fossSpecimens[i];
//		ss << "\t" << f->getCalibrationDistance();
//	}
	for(int i=0; i<fossSpecimens.size(); i++){
		f = fossSpecimens[i];
		if(1) 
			ss << "\t" << f->getFossilSppTime() * treeScale;
		else ss << "\t" << f->getFossilAge();
	}
	for(int i=0; i<fossSpecimens.size(); i++){
		f = fossSpecimens[i];
		ss << "\t" << f->getFossilFossBrGamma();
	}
	string ni = ss.str();
	return ni;
}

string Tree::getCalBDSSNodeInfoIndicatorNames(void){
	
	stringstream ss;
	Fossil *f = NULL;
	for(int i=0; i<fossSpecimens.size(); i++){
		f = fossSpecimens[i];
		int nID = f->getFossilMRCANodeID();
		ss << "\tI_f(C" << i << ".nd" << nID << ")";
	}
	string ni = ss.str();
	return ni;
}

string Tree::getCalBDSSNodeInfoIndicatorList(void){
	
	stringstream ss;
	Fossil *f = NULL;
	for(int i=0; i<fossSpecimens.size(); i++){
		f = fossSpecimens[i];
		ss << "\t" << f->getFossilIndicatorVar();
	}
	string ni = ss.str();
	return ni;
}



void Tree::setNodeCalibrationPriors(ExpCalib *ec) {
	

	for(vector<Calibration *>::iterator v = calibNds.begin(); v != calibNds.end(); v++){
		int calbNo = -1;
		calbNo = findCalibNode((*v)->getTxN1(), (*v)->getTxN2());
		Node *p = &nodes[calbNo];
		int calDistrib = (*v)->getPriorDistributionType();
		p->setNodeCalibPrDist(calDistrib);
		
		if(calDistrib == 1){
			p->setNodeYngTime((*v)->getYngTime());
			p->setNodeOldTime((*v)->getOldTime());
		}
		else if(calDistrib > 1){
			if(expHyperPrCal){
				double lv = ec->getLambdaForNode();
				p->setNodeYngTime((*v)->getYngTime());
				p->setNodeOldTime(treeScale * 1.1);
				p->setNodeExpCalRate(lv);
				p->setIsContaminatedFossil(ec->getIsLambdaContaminationClass(lv));
			}
			else{
				p->setNodeYngTime((*v)->getYngTime());
				p->setNodeOldTime(treeScale * 1.1);
				p->setNodeExpCalRate((*v)->getCalExponRate());
			}
		}
		(*v)->setNodeIndex(calbNo);
	}
	for(vector<Calibration *>::iterator v = calibNds.begin(); v != calibNds.end(); v++){
		Node *p = &nodes[(*v)->getNodeIndex()];
		int myNdPrD = p->getNodeCalibPrDist();
		if(p != root){
			if(p->getLft()->getIsCalibratedDepth()){
				int yourNdPrD = p->getLft()->getNodeCalibPrDist();
				double pOldT = p->getNodeOldTime(), lfOldT = p->getLft()->getNodeOldTime();
				if(lfOldT > pOldT && myNdPrD == 1 && yourNdPrD == 1){
					cerr << "ERROR: problem with node calibration descendant old bound > anc old bound" << endl;
					cerr << "Node " << p->getIdx() << " conflicts with it's left descendant, node " 
						 << p->getLft()->getIdx() << endl;
					exit(1);
				}
			}
			if(p->getRht()->getIsCalibratedDepth()){
				int yourNdPrD = p->getRht()->getNodeCalibPrDist();
				double pOldT = p->getNodeOldTime(), rtOldT = p->getRht()->getNodeOldTime();
				if(rtOldT > pOldT && myNdPrD == 1 && yourNdPrD == 1){
					cerr << "ERROR: problem with node calibration descendant old bound > anc old bound" << endl;
					cerr << "Node " << p->getIdx() << " conflicts with it's right descendant, node " 
						 << p->getRht()->getIdx() << endl;
					exit(1);
				}
			}
			if(p->getAnc()->getIsCalibratedDepth()){
				double pYngT = p->getNodeYngTime(), anYngT = p->getAnc()->getNodeYngTime();
				if(anYngT < pYngT && myNdPrD == 1){
					cerr << "ERROR: problem with node calibration descendant young bound > anc young bound" << endl;
					cerr << "Node " << p->getIdx() << " conflicts with it's ancestor, node " 
						 << p->getAnc()->getIdx() << endl;
					exit(1);
				}
			}
		}
	}
}

int Tree::findCalibNode(std::string t1, std::string t2) {
	
	if(t1 == "root"){
		root->setIsCalibratedDepth(true);
		return root->getIdx();
	}
	int nidn = -1;
	int chk;
	chk = findTaxLftRht(root, t1, t2, nidn);
	if(chk != 3){
		cerr << "ERROR: problem finding calibration node (findTaxLftRht)" << endl;
		cout << t1 << "\t" << t2 << endl;
		exit(1);
	}
	if(nidn < 0){
		cerr << "ERROR: problem finding calibration node (setting nidn)" << endl;
		exit(1);
	}
	return nidn;
}

int Tree::findTaxLftRht(Node *p, std::string t1, std::string t2, int &setNd) {
	
	if(p->getIsLeaf()){
		string tn = p->getName();
		if(tn == t1)
			return 1;
		else if(tn == t2)
			return 2;
		else
			return 0;
	}
	else{
		int ld = 0;
		int rd = 0;
		ld = findTaxLftRht(p->getLft(), t1, t2, setNd);
		rd = findTaxLftRht(p->getRht(), t1, t2, setNd);
		if(ld > 0 && rd > 0){
			setNd = p->getIdx();
			if(p->getIsCalibratedDepth() == true && treeTimePrior < 6){
				cerr << "ERROR: This node (ID = " << setNd << ") has already been calibrated. WTF?" << endl;
				exit(1);
			}
			else
				p->setIsCalibratedDepth(true);
		}
		return ld + rd;
	}
	return 0;
}


void Tree::verifyTreeDebug(int iter, string pn) {
	
	getTreeDotFormat(iter, pn);
}

void Tree::getTreeDotFormat(int ngen, string pn) {
	
	ofstream out;
	out.open("testdot.dbg.dot");
	out << "digraph T {\n   ";
	out << "g [shape=record, label=\"{" << ngen << "|" << pn << "}\"]\n";
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		out << "   n" << p->getIdx();
		out << " [shape=";
		stringstream rl;
		if(p->getIsLeaf()){
			rl << "{";
			if(p->getAnc() == root)
				rl << "anc = ROOT";
			else
				rl << "anc = N" << p->getAnc()->getIdx();
			rl << "|R = " << setprecision(4) << p->getRateGVal() << "|";
			rl << p->getName() << "}";
			out << "record, label=\"" << rl.str() << "\", fillcolor=lightblue, style=filled]\n";
		}
		else{
			if(p == root){
				rl << "{ROOT|H = " << setprecision(4) << treeScale << "}";
				out << "record, style=filled, label=\"" << rl.str() << "\"]\n";
			}
			else{
				rl << "{";
				if(p->getAnc() == root)
					rl << "anc = ROOT";
				else
					rl << "anc = N" << p->getAnc()->getIdx();
				rl << "|{H = " << setprecision(5) << (p->getNodeDepth() * treeScale);
				rl << "|R = " << setprecision(4) << p->getRateGVal() << "}|";
				rl << "N" << p->getIdx() << "}";
				out << "Mrecord, label=\"" << rl.str() << "\", fillcolor=pink, style=filled]\n";
			}
		}
	}
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsLeaf() == false){
			out << "   n" << p->getIdx() << " -> n" << p->getLft()->getIdx();
			out << " [label=\"" << setprecision(4) << p->getLft()->getBranchTime()*treeScale << "\"]\n";
			out << "   n" << p->getIdx() << " -> n" << p->getRht()->getIdx();
			out << " [label=\"" << setprecision(4) << p->getRht()->getBranchTime()*treeScale << "\"]\n";
		}
	}
	out << "   {rank=same;";
	for(int i=0; i<numTaxa; i++)
		out << " n" << i << ";";
	out << "}\n";
	out << "}\n";
	out.close();
}

double Tree::getTreeCBDNodePriorProb() {
	
	double nprb = 0.0;
	if(treeTimePrior == 1)
		return 0.0;
	else if(treeTimePrior == 2){
		Speciation *s = modelPtr->getActiveSpeciation();
		double diff = s->getNetDiversification();	
		nprb = getTreeCBDNodePriorProb(diff, 0.0);
		return nprb;
	}
	else if(treeTimePrior == 3){ 
		Speciation *s = modelPtr->getActiveSpeciation();
		double diff = s->getNetDiversification();	
		double rel = s->getRelativeDeath();			
		nprb = getTreeCBDNodePriorProb(diff, rel);
		return nprb;
	}
	else if(treeTimePrior == 4 || treeTimePrior == 5){ 
		Speciation *s = modelPtr->getActiveSpeciation();
		Treescale *ts = modelPtr->getActiveTreeScale();
		double div = s->getNetDiversification();	
		double rel = s->getRelativeDeath();			
		double fossRate = s->getBDSSFossilSampRatePsi();
		double sppSampRate = s->getBDSSSppSampRateRho();
		double treeOriginTime = ts->getTreeOriginTime();
		nprb = getTreeBDSSTreeNodePriorProb(div, rel, fossRate, sppSampRate,treeOriginTime);
		return nprb;
	}
	else if(treeTimePrior == 6){
		Speciation *s = modelPtr->getActiveSpeciation();
		s->setAllBDFossParams();
		double lambda = s->getBDSSSpeciationRateLambda();	
		double mu = s->getBDSSExtinctionRateMu();			
		double fossRate = s->getBDSSFossilSampRatePsi();
		double sppSampRate = s->getBDSSSppSampRateRho();
		nprb = getTreeCalBDSSTreeNodePriorProb(lambda, mu, fossRate, sppSampRate);
		return nprb;
	}
	else if(treeTimePrior == 7){ // FBD conditioned on the root
		Speciation *s = modelPtr->getActiveSpeciation();
		s->setAllBDFossParams();
		double lambda = s->getBDSSSpeciationRateLambda();	
		double mu = s->getBDSSExtinctionRateMu();			
		double fossRate = s->getBDSSFossilSampRatePsi();
		double sppSampRate = s->getBDSSSppSampRateRho();
		nprb = getTreeAncCalBDSSTreeNodePriorProb(lambda, mu, fossRate, sppSampRate);
		return nprb;
	}
    else if(treeTimePrior == 8){ // FBD conditioned on the origin time
        Speciation *s = modelPtr->getActiveSpeciation();
        s->setAllBDFossParams();
        double lambda = s->getBDSSSpeciationRateLambda();
        double mu = s->getBDSSExtinctionRateMu();
        double fossRate = s->getBDSSFossilSampRatePsi();
        double sppSampRate = s->getBDSSSppSampRateRho();
        nprb = getTreeStemAncCalBDSSTreeNodePriorProb(lambda, mu, fossRate, sppSampRate);
        return nprb;
    }
	return 0.0;
}

double Tree::getTreeCBDNodePriorProb(double netDiv, double relDeath) {
	
	double nprb = 0.0;
	if(treeTimePrior == 1)
		return 0.0;
	else if(treeTimePrior == 2){
		double diff = netDiv;
		for(int i=0; i<numNodes; i++){
			Node *p = &nodes[i];
			if(p->getIsLeaf() == false){
				double nh = p->getNodeDepth() * treeScale;
				double l = (-(diff)*nh); 
				if(p == root){
					l += (-(diff)*nh);
				}
				nprb += l;
			}
		}
		return nprb;
	}
	else if(treeTimePrior == 3){ 
		double diff = netDiv;	
		double rel = relDeath;
		for(int i=0; i<numNodes; i++){
			Node *p = &nodes[i];
			if(p->getIsLeaf() == false){
				double nh = p->getNodeDepth() * treeScale;
				double zn = log(1 - ((rel) * exp(-(diff)*nh)));
				double l = (-2 * zn) + (-(diff)*nh); 
				if(p == root){
					l += (-(diff)*nh) - zn;
				}
				nprb += l;
			}
		}
		return nprb;
	}
	else if(treeTimePrior == 4 || treeTimePrior == 5){
		Speciation *s = modelPtr->getActiveSpeciation();
		Treescale *ts = modelPtr->getActiveTreeScale();
		double fossRate = s->getBDSSFossilSampRatePsi();
		double sppSampRate = s->getBDSSSppSampRateRho();
		double treeOriginTime = ts->getTreeOriginTime();
		
		nprb = getTreeBDSSTreeNodePriorProb(netDiv, relDeath, fossRate, sppSampRate, treeOriginTime);
		
		return nprb;
	}
	else if(treeTimePrior == 6){
		Speciation *s = modelPtr->getActiveSpeciation();
		s->setAllBDFossParams();
		double lambda = s->getBDSSSpeciationRateLambda();	
		double mu = s->getBDSSExtinctionRateMu();
		double fossRate = s->getBDSSFossilSampRatePsi();
		double sppSampRate = s->getBDSSSppSampRateRho();
		nprb = getTreeCalBDSSTreeNodePriorProb(lambda, mu, fossRate, sppSampRate);
		return nprb;
	}
    else if(treeTimePrior == 7){ // FBD conditioned on the root
        Speciation *s = modelPtr->getActiveSpeciation();
        s->setAllBDFossParams();
        double lambda = s->getBDSSSpeciationRateLambda();
        double mu = s->getBDSSExtinctionRateMu();
        double fossRate = s->getBDSSFossilSampRatePsi();
        double sppSampRate = s->getBDSSSppSampRateRho();
        nprb = getTreeAncCalBDSSTreeNodePriorProb(lambda, mu, fossRate, sppSampRate);
        return nprb;
    }
    else if(treeTimePrior == 8){ // FBD conditioned on the origin time
        Speciation *s = modelPtr->getActiveSpeciation();
        s->setAllBDFossParams();
        double lambda = s->getBDSSSpeciationRateLambda();
        double mu = s->getBDSSExtinctionRateMu();
        double fossRate = s->getBDSSFossilSampRatePsi();
        double sppSampRate = s->getBDSSSppSampRateRho();
        nprb = getTreeStemAncCalBDSSTreeNodePriorProb(lambda, mu, fossRate, sppSampRate);
        return nprb;
    }
	return 0.0;
}

double Tree::getTreeBDSSTreeNodePriorProb(double netDiv, double relDeath, double fossRate, double sppSampRate, double originTime) {
	
	double nprb = 0.0;
	double lambda = netDiv / (1.0 - relDeath); 
	double mu = relDeath * lambda;
	
	int numExtantTips = numTaxa - numExtinctTips;
	if(treeTimePrior == 4){ 
		nprb = -(bdssQFxn(lambda,mu,fossRate,sppSampRate,originTime));
		nprb += numExtantTips * log(4.0 * sppSampRate);
		for(int i=0; i<numNodes; i++){
			Node *p = &nodes[i];
			if(p->getIsLeaf()){
				if(p->getIsCalibratedDepth()){ 
					double nh = p->getNodeDepth() * treeScale;
					nprb += log(fossRate) + bdssQFxn(lambda,mu,fossRate,sppSampRate,nh);
					nprb += log(bdssP0Fxn(lambda,mu,fossRate,sppSampRate,nh));
				}
			}
			else if(!p->getIsLeaf()){
				double nh = p->getNodeDepth() * treeScale;
				nprb += log(lambda) - bdssQFxn(lambda,mu,fossRate,sppSampRate,nh);
			}
		}
	}
	else if(treeTimePrior == 5){
		nprb = ((numTaxa - 2.0) * log(lambda)) - (2.0 * log(1.0 - bdssP0Fxn(lambda, mu, 0.0, sppSampRate, treeScale)));
		nprb += (numExtantTips * log(4.0 * sppSampRate)) - bdssQFxn(lambda,mu,fossRate,sppSampRate,treeScale);
		for(int i=0; i<numNodes; i++){
			Node *p = &nodes[i];
			if(p->getIsLeaf()){
				if(p->getIsCalibratedDepth()){
					double nh = p->getNodeDepth() * treeScale;
					nprb += log(fossRate) + bdssQFxn(lambda,mu,fossRate,sppSampRate,nh);
					nprb += log(bdssP0Fxn(lambda,mu,fossRate,sppSampRate,nh));
				}
			}
			else if(!p->getIsLeaf()){
				double nh = p->getNodeDepth() * treeScale;
				nprb += -bdssQFxn(lambda,mu,fossRate,sppSampRate,nh);
			}
		}
		
	}
	return nprb;
}

double Tree::getTreeCalBDSSTreeNodePriorProb(double lambda, double mu, double fossRate, double sppSampRate) {
	
	double nprb = 0.0;
	recountFossilAttachNums();
	
	nprb = ( ((numTaxa - 2 + numCalibNds) * log(lambda)) + (numCalibNds * log(fossRate)) );
	nprb -= ( 2.0 * log(1.0 - bdssP0Fxn(lambda, mu, 0.0, sppSampRate, treeScale)) );
	nprb += ( numTaxa * log(4.0 * sppSampRate) ) - bdssQFxn(lambda, mu, fossRate, sppSampRate, treeScale);
	
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(!p->getIsLeaf()){
			double myAge = p->getNodeDepth() * treeScale;
			nprb += -bdssQFxn(lambda, mu, fossRate, sppSampRate, myAge);
		}
	}
	for(int i=0; i<fossSpecimens.size(); i++){
		Fossil *f = fossSpecimens[i];
		double fossAge = f->getFossilAge();
		double fossPhi = f->getFossilSppTime() * treeScale;
		nprb += log(2.0 * f->getFossilFossBrGamma()) + bdssQFxn(lambda, mu, fossRate, sppSampRate, fossAge);
		nprb -= bdssQFxn(lambda, mu, fossRate, sppSampRate, fossPhi);
	}
	
	return nprb;
}

double Tree::getTreeAncCalBDSSTreeNodePriorProb(double lambda, double mu, double fossRate, double sppSampRate) {
	
	double nprb = 0.0;
	recountFossilAttachNums();
		
	nprb = -(2.0 * (log(lambda) + log(1.0 - bdssP0HatFxn(lambda, mu, sppSampRate, treeScale))));
	nprb += (log(4.0 * lambda * sppSampRate)) - bdssQFxn(lambda, mu, fossRate, sppSampRate, treeScale);

	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(!p->getIsLeaf()){
			double myAge = p->getNodeDepth() * treeScale;
			nprb += log(2.0 * 4.0 * lambda * sppSampRate) - bdssQFxn(lambda, mu, fossRate, sppSampRate, myAge);
		}
	}

	for(int i=0; i<fossSpecimens.size(); i++){
		Fossil *f = fossSpecimens[i];
		nprb += log(fossRate * f->getFossilFossBrGamma() );
		if(f->getFossilIndicatorVar()){
			double fossAge = f->getFossilAge();
			double fossPhi = f->getFossilSppTime() * treeScale;
			double fossPr = log(2.0 * lambda) + bdssQFxn(lambda, mu, fossRate, sppSampRate, fossAge);
			fossPr += log( bdssP0Fxn(lambda, mu, fossRate, sppSampRate, fossAge) );
			fossPr -= bdssQFxn(lambda, mu, fossRate, sppSampRate, fossPhi);
			nprb += fossPr;
		}
	}
	
	return nprb;
}


double Tree::getTreeStemAncCalBDSSTreeNodePriorProb(double lambda, double mu, double fossRate, double sppSampRate) {
    OriginTime *ot = modelPtr->getActiveOriginTime();
    originTime = ot->getOriginTime();

	double nprb = 1.0 - (log(lambda) + log(1.0 - bdssP0HatFxn(lambda, mu, sppSampRate, originTime)));
	
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(!p->getIsLeaf()){
			double myAge = p->getNodeDepth() * treeScale;
			nprb += log(lambda) + fbdQHatFxn(lambda, mu, fossRate, sppSampRate, myAge);
		}
	}
	
	for(int i=0; i<fossSpecimens.size(); i++){
		Fossil *f = fossSpecimens[i];
		nprb += log(fossRate * f->getFossilFossBrGamma() );
		if(f->getFossilIndicatorVar()){
			double fossAge = f->getFossilAge();
			double fossPhi = f->getFossilSppTime() * treeScale;
			double fossPr = log(2.0 * lambda) + log( bdssP0Fxn(lambda, mu, fossRate, sppSampRate, fossAge) );
			fossPr += fbdQHatFxn(lambda, mu, fossRate, sppSampRate, fossPhi);
			fossPr -= fbdQHatFxn(lambda, mu, fossRate, sppSampRate, fossAge);
			nprb += fossPr;
		}
	}
    //cout << "FBDS prob = " << nprb << endl;
	return nprb;
}


double Tree::getTreeSpeciationProbability() {
	
	if(treeTimePrior == 1)
		return 0.0;
	else if(treeTimePrior == 2){ 
		Speciation *s = modelPtr->getActiveSpeciation();
		double diff = s->getNetDiversification();	
		double c1 = (numTaxa - 1) * log(diff); 
		double nps = getTreeCBDNodePriorProb(diff, 0.0);
		return c1 + nps;
	}
	else if(treeTimePrior == 3){ 
		Speciation *s = modelPtr->getActiveSpeciation();
		double diff = s->getNetDiversification();
		double rel = s->getRelativeDeath();
		double lnC = 0.0;
		double c1 = ((numTaxa - 1) * log(diff)) + (numTaxa * log(1 - rel));
		double nps = getTreeCBDNodePriorProb(diff, rel);
		return lnC + c1 + nps;
	}
	else if(treeTimePrior > 3){
		double lnC = 0.0;
		double nps = getTreeCBDNodePriorProb();
		return lnC + nps;
	}
	return 0.0;
}


double Tree::getSumOfNodeHeights() {
	
	double sumh = 0.0;
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsLeaf() == false){
			double nh = p->getNodeDepth() * treeScale;
			sumh += nh;
		}
	}
	return sumh;
}

void Tree::setNodeRateValues() {
	
	NodeRate *nr = modelPtr->getActiveNodeRate();
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
#	if ASSIGN_ROOT
		int nID = p->getIdx();
		p->setRtGrpVal(nr->getRateForNodeIndexed(nID));
		p->setRtGrpIdx(nr->getTableNumForNodeIndexed(nID));		
#	else
		if(p != root){
			int nID = p->getIdx();
			p->setRtGrpVal(nr->getRateForNodeIndexed(nID));
			p->setRtGrpIdx(nr->getTableNumForNodeIndexed(nID));			
		}
		else{
			p->setRtGrpVal(0.0);
			p->setRtGrpIdx(0.0);
		}
#	endif
	}
}

void Tree::assignMixLambdaHyperPrToNode(Node *p){
	
	ExpCalib *ec = modelPtr->getActiveExpCalib();
	p->setNodeExpCalRate(ec->getLambdaForNode());
	
}

double Tree::getAMixLambdaHyperPrToNode(){
	
	ExpCalib *ec = modelPtr->getActiveExpCalib();
	return ec->getLambdaForNode();
	
}


vector<Node *> Tree::getListOfCalibratedNodes(){
	
	vector<Node *> calibList;
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsCalibratedDepth()) 
			calibList.push_back(p);
	}
	return calibList;
}

void Tree::writeCalNodeBEASTInfoXML(std::ostream &o){
	
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsLeaf() == false && p != root){
			vector<string> myDecs;
			getListOfTaxNamesDecFromNode(p, myDecs);
			int nodeIDX = p->getIdx();
			o << "\t<taxa id=\"N" << nodeIDX << "\">\n";
			for(int t=0; t<myDecs.size(); t++)
				o << "\t\t<taxon idref=\"" << myDecs[t] << "\"/>\n";
			o << "\t</taxa>\n";
			myDecs.clear();
		}
	}
	o << "\nSPLIT\n";
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsCalibratedDepth()){
			int nodeIDX = p->getIdx();
			double lb = p->getNodeYngTime();
			double ub = p->getNodeOldTime();
			o << "\t\t\t\t<uniformPrior lower=\"" << lb << "\" upper=\"" << ub << "\">\n";
			if(p == root)
				o << "\t\t\t\t\t<parameter idref=\"treeModel.rootHeight\"/>\n";
			else
				o << "\t\t\t\t\t<statistic idref=\"tmrca(N" << nodeIDX << ")\"/>\n";
			o << "\t\t\t\t</uniformPrior>\n";
		}
	}
	o << "\nSPLIT\n";

	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsLeaf() == false && p != root){
			int nodeIDX = p->getIdx();
			o << "\t<tmrcaStatistic id=\"tmrca(";
			if(p == root)
				o << "root(N" << nodeIDX << "))\">\n";
			else
				o << "N" << nodeIDX << ")\">\n";
			o << "\t\t<mrca>\n\t\t\t<taxa idref=\"";
			if(p == root)
				o << "root(N" << nodeIDX << ")\"/>\n";
			else
				o << "N" << nodeIDX << "\"/>\n";
			o << "\t\t</mrca>\n\t\t<treeModel idref=\"treeModel\"/>\n";
			o << "\t</tmrcaStatistic>\n";
		}
	}

	o << "\nSPLIT\n";
	
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsLeaf() == false && p != root){
			int nodeIDX = p->getIdx();
			o << "\t\t\t<tmrcaStatistic idref=\"tmrca(";
			o << "N" << nodeIDX << ")\"/>\n";
		}
	}
	
}

void Tree::writeRRTNodeBEASTInfoXML(std::ostream &o){
	
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsLeaf() == false && p != root){
			vector<string> myDecs;
			getListOfTaxNamesDecFromNode(p, myDecs);
			int nodeIDX = p->getIdx();
			o << "\t<taxa id=\"N" << nodeIDX << "\">\n";
			for(int t=0; t<myDecs.size(); t++)
				o << "\t\t<taxon idref=\"" << myDecs[t] << "\"/>\n";
			o << "\t</taxa>\n";
			myDecs.clear();
		}
	}
	o << "\nSPLIT\n";
	
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsLeaf() == false && p != root){
			int nodeIDX = p->getIdx();
			o << "\t<tmrcaStatistic id=\"tmrca(";
			if(p == root)
				o << "root(N" << nodeIDX << "))\">\n";
			else
				o << "N" << nodeIDX << ")\">\n";
			o << "\t\t<mrca>\n\t\t\t<taxa idref=\"";
			if(p == root)
				o << "root(N" << nodeIDX << ")\"/>\n";
			else
				o << "N" << nodeIDX << "\"/>\n";
			o << "\t\t</mrca>\n\t\t<treeModel idref=\"treeModel\"/>\n";
			o << "\t</tmrcaStatistic>\n";
		}
	}
	
	o << "\nSPLIT\n";
	
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsLeaf() == false && p != root){
			int nodeIDX = p->getIdx();
			o << "\t\t\t<tmrcaStatistic idref=\"tmrca(";
			o << "N" << nodeIDX << ")\"/>\n";
		}
	}
	
}

void Tree::getListOfTaxNamesDecFromNode(Node *p, std::vector<std::string> &ndNs){
	
	if(p->getIsLeaf())
		ndNs.push_back(p->getName());
	else{
		getListOfTaxNamesDecFromNode(p->getLft(), ndNs);
		getListOfTaxNamesDecFromNode(p->getRht(), ndNs);
	}
	
}

void Tree::checkNodeCalibrationCompatibility(){

	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsCalibratedDepth()){
			double lb = p->getNodeYngTime();
			double ub = p->getNodeOldTime();
			double myDepth = p->getNodeDepth() * treeScale;
			if(myDepth > ub || myDepth < lb){
				cerr << "Calibrated node: " << i << " has the depth: " << myDepth << " and is NOT in the range {" << lb << ", " << ub << "}" << endl;
				exit(1);
			}
			else
				cout << "Calibrated node: " << i << " has the depth: " << myDepth << " and is in the range {" << lb << ", " << ub << "}" << endl;
		}
	}
}

void Tree::zeroNodeRedFlags(){
	
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		p->setRedFlag(0);
		if(p->getIsLeaf() == false && p != root){
			p->setNodeDepth(0.0);
		}
	}
}


int Tree::checkTreeForCalibrationCompatibility(){
	
	int numIncomp = 0;
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p != root){
			double branchtime = (p->getAnc()->getNodeDepth() - p->getNodeDepth());
			if(branchtime < 0){
				p->setRedFlag(2);
				numIncomp++;
			}
			if(p->getIsCalibratedDepth() && p->getNodeCalibPrDist() == 1){
				double lb = p->getNodeYngTime();
				double ub = p->getNodeOldTime();
				double myDepth = p->getNodeDepth() * treeScale;
				if(lb == ub){
					p->setNodeDepth(ub / treeScale);
					myDepth = ub;
				}
				
				if(myDepth > ub || myDepth < lb){
					p->setRedFlag(1);
					numIncomp++;
					cout << p->getIdx() << " -- " << myDepth << " -- " << lb << " - " << ub << endl;
					if(myDepth > ub)
						cout << p->getIdx() << " -- " << myDepth << " > " << ub << endl;
					if(myDepth < lb)
						cout << p->getIdx() << " -- " << myDepth << " < " << lb << endl;
				}
			}
		}
	}
	return numIncomp;
}


void Tree::initializeTGSCalibVariables(){  // depricated
	
	for(int i=0; i<numNodes; i++){
		Node *p = &nodes[i];
		if(p->getIsCalibratedDepth()){
			double fA = p->getNodeYngTime() / treeScale;
			double currDepth = p->getNodeDepth();
			double phii = fA + ranPtr->uniformRv()*(currDepth-fA);
			int xi = countDecLinsTimeIntersect(p, phii, currDepth);
			p->setFossAttchTime(phii);
			p->setNumFossAttchLins(xi);
		}
	}
}

int Tree::countDecLinsTimeIntersect(Node *p, double t, double ancAge){
	
	int n = 0;
	double myAge = p->getNodeDepth();
	if(myAge < t)
		return 1;
	else{
		n += countDecLinsTimeIntersect(p->getLft(), t, myAge);
		n += countDecLinsTimeIntersect(p->getRht(), t, myAge);
	}
	
	
	return n;
}

double Tree::getNodeLowerBoundTime(Node *p){
	
//	bool q = p->getIsCalibratedDepth();
	double t = p->getLft()->getNodeDepth();
	if (p->getRht()->getNodeDepth() > t)
		t = p->getRht()->getNodeDepth();
	if(softBounds == false){
		double fA = p->getNodeYngTime() / treeScale;
		if(fA > t)
			t = fA;
		if(treeTimePrior > 5){
			double phi = p->getFossAttchTime();
			if(phi > t)
				t = phi;
		}
		
	}
	return t;
}

double Tree::getNodeUpperBoundTime(Node *p){
	
	double t = p->getAnc()->getNodeDepth();
	
	if(p->getIsCalibratedDepth() && softBounds == false){
		double ocal =  t;
		if(p->getNodeCalibPrDist() == 1)
			ocal = p->getNodeOldTime() / treeScale;
		if(ocal < t)
			t = ocal;
	}
	return t;
}

void Tree::setUPTGSCalibrationFossils() { // definitely have to change for FBDS
	
	int calID = 0;
	for(vector<Calibration *>::iterator v = calibNds.begin(); v != calibNds.end(); v++){
		int calbNo = -1;
		calbNo = findCalibNode((*v)->getTxN1(), (*v)->getTxN2());
		Node *p = &nodes[calbNo];
		double fAge = (*v)->getYngTime();
		Fossil *f = new Fossil(fAge, calbNo);
		fossSpecimens.push_back(f);
		int nPFoss = p->getNumCalibratingFossils();
		p->setNumFCalibratingFossils(nPFoss + 1);
        bool sf = (*v)->getIsStemFossil();
        f->setIsTotalGroupFossil(sf);
		
		if(nPFoss > 0 && sf == false){
			if(p->getNodeYngTime() < fAge)
				p->setNodeYngTime(fAge);				
		}
		else p->setNodeYngTime(fAge);
		p->setNodeOldTime(treeScale * 1.1);
		p->setNodeExpCalRate((*v)->getCalExponRate());
		
		(*v)->setNodeIndex(calbNo);
        p->insertFossilID(calID);
        p->insertTGFossilID(calID);
		cout << calbNo << " --- " << p->getNodeYngTime() << endl;
		calID += 1;
	}
}

void Tree::initializeFossilSpecVariables(){
	
	numAncFossilsk = 0;
	for(int i=0; i<fossSpecimens.size(); i++){
		Fossil *f = fossSpecimens[i];
		int nID = f->getFossilMRCANodeID();
		Node *p = &nodes[nID];
		double fA = f->getFossilAge() / treeScale;
        double currDepth = p->getNodeDepth();
        bool isTGFossil = f->getIsTotalGroupFossil();
        if(isTGFossil){
            if(p != root){
                currDepth = p->getAnc()->getNodeDepth();
            }
            else{
                currDepth = originTime/treeScale;
            }
        }
		double prTip = 1.0;
		if(treeTimePrior >= 7)
            prTip = 0.5;
		if(ranPtr->uniformRv() < prTip){
			double phi = fA + ranPtr->uniformRv()*(currDepth-fA);
			f->setFossilSppTime(phi); // initialises the attachment times
			f->setFossilIndicatorVar(1);
//            if(fA == 45/treeScale)
//                f->setFossilSppTime(80.0/treeScale);
//            else if(fA == 17.0/treeScale)
//                f->setFossilSppTime(25.0/treeScale);
//            else if(fA == 70.0/treeScale)
//                f->setFossilSppTime(85.0/treeScale);
//            else if(fA == 60.0/treeScale)
//                f->setFossilSppTime(75.0/treeScale);
//            else if(fA == 27.0/treeScale)
//                f->setFossilSppTime(55.0/treeScale);
//            else if(fA == 30.0/treeScale)
//                f->setFossilSppTime(40.0/treeScale);
//            else if(fA == 26.0/treeScale)
//                f->setFossilSppTime(33.0/treeScale);
//            
           cout << f->getFossilAge() << "  --  " << f->getFossilSppTime()*treeScale << "  cd = " << currDepth << endl;

//			double nodeBound = p->getFossAttchTime();
//			if(phi > nodeBound)
//				p->setFossAttchTime(phi);
		}
		else{
			numAncFossilsk += 1;
			f->setFossilSppTime(fA);
			f->setFossilFossBrGamma(0);
			f->setFossilIndicatorVar(0);
			double nodeBound = p->getFossAttchTime();
			if(fA > nodeBound)
				p->setFossAttchTime(fA);

		}
	}
    treeUpdateNodeOldestBoundsAttchTimes();
//	recountFossilAttachNums();
//    for(int i=0; i<fossSpecimens.size(); i++){
//        Fossil *f = fossSpecimens[i];
//        cout << f->getFossilAge() << " -- g_f = " << f->getFossilFossBrGamma() << " -- z_f = " << f->getFossilSppTime() << endl;
//    }
//    recountFossilAttachNums();

}



int Tree::recountFossilAttachNums(){
	
	int numChanged = 0;
	for(int i=0; i<fossSpecimens.size(); i++){
		Fossil *fi = fossSpecimens[i];
		int fgammCur = fi->getFossilFossBrGamma();
		double fiPhi = fi->getFossilSppTime(); // z_f
		Node *p = &nodes[fi->getFossilMRCANodeID()];
        double xi = p->getNodeDepth() * treeScale;
        int fiGamma = 0;
        if(fi->getIsTotalGroupFossil()){
            if(fiPhi > xi)
                fiGamma = 1;
            else fiGamma = getDecendantFossilAttachBranches(p, fiPhi, i);
            fiGamma += getGammaValueForTGFossil(fi);
            
        }
        else if (fiPhi < xi)
            fiGamma = getDecendantFossilAttachBranches(p, fiPhi, i);
        else{
            cerr << "Error: z_i > x_i for a crown fossil." << endl;
            exit(1);
        }
		if(fiGamma != fgammCur)
			numChanged += 1;
		fi->setFossilFossBrGamma(fiGamma);
		fi->setFossilMRCANodeAge(p->getNodeDepth() * treeScale);
	}
	return numChanged;
}

int Tree::getGammaValueForTGFossil(Fossil *f){
    
    Node *p = &nodes[f->getFossilMRCANodeID()];
    double xi = p->getNodeDepth() * treeScale;
    double zf = f->getFossilSppTime();
    int sumGamma = 0;
    vector<int> myTGFossils = p->getTGFossilIDs();
    // for tg fossils of xi
    for (int i = 0; i < myTGFossils.size(); i++){
        Fossil *tGFossilI = fossSpecimens[i];
        if (tGFossilI != f) {
            double ztgf = tGFossilI->getFossilSppTime();
            double ytgf = tGFossilI->getFossilAge();
            if(ztgf > xi) {
                if(ztgf > zf && ytgf < zf){
                    sumGamma += 1;
                }
            }
        }
    }
    // for fossils of xpi and attached to the stem of xi
    if(p != root){
        Node *xAnc = p->getAnc();
        double ancAge = xAnc->getNodeDepth() * treeScale;
        vector<int> ancFossils = xAnc->getCalibFossilIDs();
        for (int i = 0; i < ancFossils.size(); i++){
            Fossil *ancFossI = fossSpecimens[i];
            double zAncF = ancFossI->getFossilSppTime();
            if (zAncF > xi && zAncF < ancAge){
                double yAncF = ancFossI->getFossilAge();
                if(zAncF > zf && yAncF < zf){
                    sumGamma += 1;
                }
            }
        }
    }
    
    return sumGamma;
}

double Tree::getSumLogAllAttachNums(){
	
	double sumLog = 0.0;
	recountFossilAttachNums();
	for(int i=0; i<fossSpecimens.size(); i++){
		Fossil *fi = fossSpecimens[i];
		sumLog += log(fi->getFossilFossBrGamma());
	}
	return sumLog;
}

int Tree::getSumIndicatorV(){
	if(treeTimePrior < 6)
		return 0;
	int niv = 0;
	for(int i=0; i<fossSpecimens.size(); i++){
		Fossil *fi = fossSpecimens[i];
		niv += fi->getFossilIndicatorVar();
	}
	return niv;
}


//int Tree::getFossilLinAttachNumsForFoss(int fID){ // TAH: this is never called
//	
//	Fossil *fi = fossSpecimens[fID];
//	double fiPhi = fi->getFossilSppTime();
//	Node *p = &nodes[fi->getFossilMRCANodeID()];
//	int g = 0;
//	if(fi->getFossilIndicatorVar())
//		g = getDecendantFossilAttachBranches(p, fiPhi, fID);
//	
//	return g;
//	
//}

int Tree::getDecendantFossilAttachBranches(Node *p, double t, int fID){
	
	int n = 0;
	
	double nodeAge = p->getNodeDepth();
	if(nodeAge < t)
		return 1; 
	else{
		for(int i=0; i<p->getNumCalibratingFossils(); i++){
			int idx = p->getIthFossiID(i);
			if(idx != fID){
				Fossil *f = fossSpecimens[idx];
				if(f->getFossilSppTime() > t && (f->getFossilAge()/treeScale) < t)
					n += 1;
			}
		}
		n += getDecendantFossilAttachBranches(p->getLft(), t, fID);
		n += getDecendantFossilAttachBranches(p->getRht(), t, fID);
	}
	
	return n;
}


double Tree::getProposedBranchAttachFossils(Node *p, double t){
	
	double sumLog = 0.0;
	
	double nodeAge = p->getNodeDepth();
	if(nodeAge < t || p->getIsLeaf())
		return 0.0;
	else{
		for(int i=0; i<p->getNumCalibratingFossils(); i++){
			int idx = p->getIthFossiID(i);
			int a = 0;
			Fossil *f = fossSpecimens[idx];
			a = getDecendantFossilAttachBranches(p,f->getFossilSppTime(),idx);
			sumLog += log(a);
		}
		sumLog += getProposedBranchAttachFossils(p->getLft(), t);
		sumLog += getProposedBranchAttachFossils(p->getRht(), t);
	}
	return sumLog;
}

void Tree::treeScaleUpdateFossilAttchTimes(double sr, double ort, double nrt){
	

	for(int i=0; i<fossSpecimens.size(); i++){
		Fossil *f = fossSpecimens[i];
		int nID = f->getFossilMRCANodeID();
		Node *p = &nodes[nID];
		double oldPhi = f->getFossilSppTime();
		double newPhi = oldPhi * sr;
		f->setFossilSppTime(newPhi);
		double nodeBound = p->getFossAttchTime();
		if(newPhi > nodeBound)
			p->setFossAttchTime(newPhi);
	}
	recountFossilAttachNums();
	
}

void Tree::treeUpdateNodeOldestBoundsAttchTimes(){
	
	// make this for all nodes, until ready to optimize computation
//	for(vector<Node *>::iterator v = calibNodes.begin(); v != calibNodes.end(); v++){
//		setNodeOldestAttchBranchTime((*v));
//	}
    for(int i; i<numNodes; i++){
        Node *p = &nodes[i];
        if(!p->getIsLeaf())
            setNodeOldestAttchBranchTime(p);
    }

}

void Tree::setNodeOldestAttchBranchTime(Node *p){
	
    // add query to daughter nodes for stem fossils
	double t = 0.0;
	for(int i=0; i<p->getNumCalibratingFossils(); i++){
		int idx = p->getIthFossiID(i);
		Fossil *f = fossSpecimens[idx];
		double phi = f->getFossilSppTime();
		if(phi > t)
			t = phi;
	}
    Node *ld = p->getLft();
    Node *rd = p->getRht();
    for(int i=0; i<ld->getNumCalibratingFossils(); i++){
        int idx = ld->getIthFossiID(i);
        Fossil *f = fossSpecimens[idx];
        if(f->getIsTotalGroupFossil()){
            double phi = f->getFossilSppTime();
            if(phi > t)
                t = phi;
        }
    }
    for(int i=0; i<rd->getNumCalibratingFossils(); i++){
        int idx = rd->getIthFossiID(i);
        Fossil *f = fossSpecimens[idx];
        if(f->getIsTotalGroupFossil()){
            double phi = f->getFossilSppTime();
            if(phi > t)
                t = phi;
        }
    }
    

	p->setFossAttchTime(t);
    
}

int Tree::pickRandAncestorFossil(){
	
	vector<int> af;
	for(int i=0; i<fossSpecimens.size(); i++){
		Fossil *f = fossSpecimens[i];
		if(f->getFossilIndicatorVar() == 0)
			af.push_back(i);
	}
	int v = (int)(ranPtr->uniformRv()*af.size());
	return af[v];
}

int Tree::pickRandTipFossil(){
	
	vector<int> af;
	for(int i=0; i<fossSpecimens.size(); i++){
		Fossil *f = fossSpecimens[i];
		if(f->getFossilIndicatorVar() == 1)
			af.push_back(i);
	}
	int v = (int)(ranPtr->uniformRv()*af.size());
	return af[v];
}


// END
