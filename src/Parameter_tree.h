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

#ifndef PARAMETER_TREE_H
#define PARAMETER_TREE_H


#include <string>
#include <sstream>
#include <vector>

class Node {

	public:
						Node(void);
		void			flipCl(void) { (activeCl == 0 ? activeCl = 1 : activeCl = 0); }
		void			flipTi(void) { (activeTi == 0 ? activeTi = 1 : activeTi = 0); }
		Node*			getLft(void) { return lft; }
		Node*			getRht(void) { return rht; }
		Node*			getAnc(void) { return anc; }
		int				getIdx(void) const { return idx; }
		std::string		getName(void) { return name; }
		int				getActiveCl(void) { return activeCl; }
		int				getActiveTi(void) { return activeTi; }
		bool			getIsClDirty(void) { return isClDirty; }
		bool			getIsTiDirty(void) { return isTiDirty; }
		double			getNodeDepth(void) { return nodeDepth; }
		bool			getIsLeaf(void) { return isLeaf; }
		bool			getIsCalibratedDepth(void) { return isCalib; }
		void			setLft(Node *p) { lft = p; }
		void			setRht(Node *p) { rht = p; }
		void			setAnc(Node *p) { anc = p; }
		void			setIdx(int x) { idx = x; }
		void			setName(std::string s) { name = s; }
		void			setActiveCl(int x) { activeCl = x; }
		void			setActiveTi(int x) { activeTi = x; }
		void			setIsClDirty(bool x) { isClDirty = x; }
		void			setIsTiDirty(bool x) { isTiDirty = x; }
		void			setNodeDepth(double x) { nodeDepth = x; }
		void			setIsLeaf(bool x ) { isLeaf = x; }
		void			setIsCalibratedDepth(bool x ) { isCalib = x; }
		double			getBranchTime() { return branchTime; }
		double			getRateGVal() { return rateGVal; }
		double			getUerBL() { return userBL; }
		int				getRateGrpIdx() { return rateGrpIdx; }
		void			setBranchTime(double v) { branchTime = v; }
		void			setRtGrpVal(double v) { rateGVal = v; }
		void			setRtGrpIdx(int i) { rateGrpIdx = i; }
		void			setUerBL(double v) { userBL = v; }
		double			getNodeYngTime() { return youngt; }
		double			getNodeOldTime() { return oldt; }
		void			setNodeYngTime(double v) { youngt = v; }
		void			setNodeOldTime(double v) { oldt = v; }
		int				getNumDecendantTax() { return numDecTax; }
		void			setNumDecendantTax(int i) { numDecTax = i; }
		int				getNodeCalibPrDist() { return nodeCalibPrD; }
		void			setNodeCalibPrDist(int i) { nodeCalibPrD = i; }
		double			getNodeExpCalRate() { return nodeExpCalRate; }
		void			setNodeExpCalRate(double v) { nodeExpCalRate = v; }
		bool			getIsContaminatedFossil() { return taintFossil; }
		void			setIsContaminatedFossil(bool b) { taintFossil = b; }
		int				getRedFlag() { return redFlag; }
		void			setRedFlag(int i) { redFlag = i; }
		double			getNodeAge(void) { return nodeAge; }
		void			setNodeAge(double d) { nodeAge = d; }
		
		double			getFossAttchTime(void) { return fossAttachTime; }
		void			setFossAttchTime(double d) { fossAttachTime = d; }
		int				getNumFossAttchLins(void) { return numFossAttachLins; }
		void			setNumFossAttchLins(int i) { numFossAttachLins = i; }
		
		// New implmentation stuff

		int				getNumCalibratingFossils(void) { return numCalFossils; }
		void			setNumFCalibratingFossils(int i) { numCalFossils = i; }
		void			insertFossilID(int i) { calibFossilIDs.push_back(i); }
		int				getIthFossiID(int i) { return calibFossilIDs[i]; }
		

	private:
		Node			*lft;
		Node			*rht;
		Node			*anc;
		int				idx;
		std::string		name;
		int				activeCl;
		int				activeTi;
		bool			isClDirty;
		bool			isTiDirty;
		double			nodeDepth;
		double			nodeAge;
		bool			isLeaf;
		bool			isCalib;
		double			youngt;
		double			oldt;
		double			branchTime;
		double			rateGVal;
		int				rateGrpIdx;
		double			userBL;
		int				numDecTax;
		int				nodeCalibPrD;
		double			nodeExpCalRate;
		bool			taintFossil;
		int				redFlag;
		
		// If TGS calibration phi and xi
		double			fossAttachTime;		// phi_i
		int				numFossAttachLins;  // xi_i
		int				numCalFossils;
		std::vector<int> calibFossilIDs;
};

class Fossil {
	
	public:
										Fossil(double fa, int nid) :age(fa), nodeID(nid),
																	fossilBrGamma(0), ancFossIndicator(1) {}
		
		int								getFossilIndex(void) { return indx; } 
		double							getFossilAge(void) { return age; }
		double							getFossilSppTime(void) { return phi; }
		int								getFossilMRCANodeID(void) { return nodeID; }
		double							getFossilMRCANodeAge(void) { return nodeAge; } 
		int								getFossilFossBrGamma(void) { return fossilBrGamma; }
		int								getFossilIndicatorVar(void) { return ancFossIndicator; }
		double							getCalibrationDistance(void) { return nodeAge - age; }
		
		void							setFossilIndex(int i) { indx = i; } 
		void							setFossilAge(double d) { age = d; }
		void							setFossilSppTime(double d) { phi = d; }
		void							setFossilMRCANodeID(int i) { nodeID = i; }
		void							setFossilMRCANodeAge(double d) { nodeAge = d; } 
		void							setFossilFossBrGamma(int i) { fossilBrGamma = i; }
		void							setFossilIndicatorVar(int i) { ancFossIndicator = i; }

		
		
	private:
		int								indx;
		double							age;
		double							phi;
		int								nodeID;
		int								gamma;
		double							nodeAge;
		int								fossilBrGamma;
		int								ancFossIndicator; // {\cal I} = 0 if anc fossil, 1 otherwise
		
	
	
};


class Alignment;
class Calibration;
class MbRandom;
class Model;
class ExpCalib;
class Tree : public Parameter {

	public:
										Tree(MbRandom *rp, Model *mp, Alignment *ap, std::string ts, 
											 bool ubl, bool allnm, bool rndNods, std::vector<Calibration *> clb, 
											 double rth, bool sb, bool exhpc, ExpCalib *ec, std::vector<Calibration *> tdt);
										~Tree(void); 
		Tree							&operator=(const Tree &t);
		void							clone(const Tree &t);
		void							getDownPassSequence(void);
		Node*							getDownPassNode(int i) { return downPassSequence[i]; }
		Node*							getNodeByIndex(int i) { return &nodes[i]; }
		int								getNumNodes(void) { return numNodes; }
		int								getNumTaxa(void) { return numTaxa; }
		Node*							getRoot(void) { return root; }
		double							update(double &oldLnL);
		double							updateOneNode();
		double							updateAllNodes(double &oldLnL);
		double							updateAllTGSNodes(double &oldLnL);
		double							updateAllNodesRnd(double &oldLnL);
		double							lnPrior();
		double							lnPriorRatio(double snh, double soh);
		double							lnPriorRatioTGS(double snh, double soh, Node *p);
		double							lnCalibPriorRatio(double nh, double oh, double lb, double ub);
		double							lnExpCalibPriorRatio(double nh, double oh, double offSt, double expRate);
		void							print(std::ostream & o) const;
		void							flipAllCls(void);
		void							flipAllTis(void);
		void							upDateAllCls(void);
		void							upDateAllTis(void);
		void							flipToRootClsTis(Node *p);
		void							updateToRootClsTis(Node *p);
		void							flipToRootClsTis(int ndID);
		void							updateToRootClsTis(int ndID);
		std::string						getTreeDescription(void);
		std::string						getFigTreeDescription(void);
		std::string						getCalibInitialTree(void);
		std::string						writeParam(void);
		std::string						getNodeInfoNames(void);
		std::string						getNodeInfoList(void);
		std::string						getDownPNodeInfoNames(void);
		std::string						getDownPNodeInfoList(void);
		std::string						getCalNodeInfoNames(void);
		std::string						getCalNodeInfoList(void);
		void							setRootRateValue(double v) { root->setRtGrpVal(v); }
		void							setAllNodeBranchTimes(void);
		void							setRndShufNdMv(bool b) { randShufNdMv = b; }
		double							getTreeScale() { return treeScale; }
		void							setTreeScale(double s) { treeScale = s; }
		void							verifyTreeDebug(int iter, std::string pn);
		bool							getIsCalibratedTree() { return isCalibTree; }
		double							getTreeCBDNodePriorProb();
		double							getTreeCBDNodePriorProb(double netDiv, double relDeath);
		double							getTreeBDSSTreeNodePriorProb(double netDiv, double relDeath, double fossRate, double sppSampRate, double originTime);
		double							getTreeSpeciationProbability();
		double							getSumOfNodeHeights();
		void							setNodeRateValues();
		std::vector<Node *>				getListOfCalibratedNodes();
		double							getRootCalibExpRate() { return root->getNodeExpCalRate(); }
		void							writeCalNodeBEASTInfoXML(std::ostream &o);
		void							writeRRTNodeBEASTInfoXML(std::ostream &o);
		void							getListOfTaxNamesDecFromNode(Node *p, std::vector<std::string> &ndNs);
		void							checkNodeCalibrationCompatibility();
		int								checkTreeForCalibrationCompatibility();
		void							zeroNodeRedFlags();
		
		double							getTreeCalBDSSTreeNodePriorProb(double lambda, double mu, double fossRate, double sppSampRate);
		double							getTreeAncCalBDSSTreeNodePriorProb(double lambda, double mu, double fossRate, double sppSampRate);
		double							getTreeStemAncCalBDSSTreeNodePriorProb(double lambda, double mu, double fossRate, double sppSampRate);

		std::string						getCalBDSSNodeInfoParamNames(void);
		std::string						getCalBDSSNodeInfoParamList(void);
		std::string						getCalBDSSNodeInfoIndicatorNames(void);
		std::string						getCalBDSSNodeInfoIndicatorList(void);
		int								countDecLinsTimeIntersect(Node *p, double t, double ancAge);
		double							updateFossilBDSSAttachmentTimePhi(void);
		void							treeScaleUpdateFossilAttchTimes(double sr, double ort, double nrt);
		void							treeUpdateNodeOldestBoundsAttchTimes(void);
		double							getNodeLowerBoundTime(Node *p);
		double							getNodeUpperBoundTime(Node *p);
		
		double							updateRJMoveAddDelEdge(void);
		void							doAddEdgeMove(void);
		void							doDeleteEdgeMove(void);
	
		int								getSumIndicatorV(void);

							
	private:
		void							buildTreeFromNewickDescription(std::string ts);
		void							initializeNodeDepthsFromUserBL(void);
		void							initializeTipNodeDepthsFromUserBL(void);
		void							initializeNodeDepths(void);
		void							initializeCalibratedNodeDepths(void);
		void							initializeTGSCalibratedNodeDepths(void);
		double							getUBLTreeScaleDepths(void);
		double							getNodePathDepth(Node *t);
		std::vector<double>				recursiveNodeDepthInitialization(Node *p, int &nCont, double maxD);
		int								setNodesNumberDecendantTaxa(Node *p);
		static int						dex(const Node *p);
		bool							isValidChar(char c);
		void							passDown(Node *p, int *x);
		void							showNodes(Node *p, int indent, std::ostream &ss) const;
		void							writeTree(Node *p, std::stringstream &ss);
		void							writeFigTree(Node *p, std::stringstream &ss);
		void							writeCalibrationFigTree(Node *p, std::stringstream &ss);
		void							setNodeCalibrationPriors(ExpCalib *ec);
		int								findCalibNode(std::string t1, std::string t2);
		int								findTaxLftRht(Node *p, std::string t1, std::string t2, int &setNd);
		void							getTreeDotFormat(int ngen, std::string pn);
		double							getTemporaryNodeMaxBound(Node *p);
		void							assignMixLambdaHyperPrToNode(Node *p);
		double							getAMixLambdaHyperPrToNode();
		
		int								findTip(std::string tn);
		
		void							setTipDateAges(void);
		
		double							bdssC1Fxn(double b, double d, double psi);
		double							bdssC2Fxn(double b, double d, double psi,double rho);
		
		double							bdssQFxn(double b, double d, double psi, double rho, double t); // on log scale
		double							bdssP0Fxn(double b, double d, double psi, double rho, double t);
		double							bdssP0HatFxn(double b, double d, double rho, double t);
		double							fbdQHatFxn(double b, double d, double psi, double rho, double t);

		void							initializeTGSCalibVariables(void);
		
		void							setUPTGSCalibrationFossils();
		void							initializeFossilSpecVariables();
		void							initializeFossilAncestorSpecVariables();
		int								recountFossilAttachNums();
		int								getFossilLinAttachNumsForFoss(int fID);
		
		int								getDecendantFossilAttachBranches(Node *p, double t, int fID);
		double							getProposedBranchAttachFossils(Node *p, double t);
		void							setNodeOldestAttchBranchTime(Node *p);
		int								pickRandAncestorFossil(void);
		int								pickRandTipFossil(void);
		double							getSumLogAllAttachNums(void);
		double							doAScaleMove(double &nv, double cv, double tv, double lb, double hb, double rv);
		double							doAWindoMove(double &nv, double cv, double tv, double lb, double hb, double rv);
		int								getNumDecFossils(Node *p);
		
		void							checkNodeInit(void);
		
		void							recountFromNodeFossilAttchNums(Node *p);
		int								recursivCreateTempFossVec(std::vector<Fossil *> &v, Node *p);
		
		
		Alignment						*alignmentPtr;
		int								numTaxa;
		int								numNodes;
		Node							*nodes;
		Node							*root;
		Node							**downPassSequence;
		bool							useInputBLs;
		bool							moveAllNodes;
		bool							randShufNdMv;
		std::vector<Calibration*>		calibNds;
		std::vector<Calibration*>		datedTips;
		int								numCalibNds; // if -tgs, then this is actually equal to the number of calibrating fossils, not the number of nodes calibrated
		int								numExtinctTips;
		int								numAncFossilsk; // for -tgs, this is the number of fossils that are ancestors
		double							treeScale;
        double                          originTime;
        bool                            conditionOnOrigin;
		bool							isCalibTree;
		int								treeTimePrior;
		bool							softBounds;
		bool							expHyperPrCal;
		bool							isTipCals;
		int								nodeProposal;
		
		std::vector<Node *>				calibNodes;
		std::vector<Fossil *>			fossSpecimens;
		
		double							tuningVal;
		int								numMoves;
		bool							autotune;
		
		
};


#endif
