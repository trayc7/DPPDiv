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

#ifndef PARAMETER_FOSSILGRAPH_H
#define PARAMETER_FOSSILGRAPH_H

#include <vector>
#include <set>
#include <string>
#include <ostream>

class Occurrence {
    
public:
    Occurrence(double fa, int oid) :age(fa), occurrenceID(oid), fossilBrGamma(0), ancFossIndicator(1), isTerminal(false) {}
    
    int								getFossilIndex(void) { return indx; }
    double							getFossilAge(void) { return age; }
    double							getFossilSppTime(void) { return phi; }
    int								getFossilFossBrGamma(void) { return fossilBrGamma; }
    int								getFossilIndicatorVar(void) { return ancFossIndicator; }
    int								getFossilID(void) { return occurrenceID; } //rw: what does getFossilIndex do?
    bool                            getIsTerminal(void) { return isTerminal; }

    
    void							setFossilIndex(int i) { indx = i; }
    void							setFossilAge(double d) { age = d; }
    void							setFossilSppTime(double d) { phi = d; }
    void							setFossilFossBrGamma(int i) { fossilBrGamma = i; }
    void							setFossilIndicatorVar(int i) { ancFossIndicator = i; }
    void							setFossilID(int i) { occurrenceID = i; }
    void                            setIsTerminal(bool b) { isTerminal = b; }
    
private:
    int								indx;
    double							age;
    double							phi;
    int								gamma;
    int                             occurrenceID;
    //double							nodeAge; //rw: not used
    int								fossilBrGamma;
    int								ancFossIndicator; // {\cal I} = 0 if anc fossil, 1 otherwise
    bool                            isTerminal;
    
};

class MbRandom;
class Model;
class Calibration;

class FossilGraph : public Parameter {
    
public:
                                    FossilGraph(MbRandom *rp, Model *mp, int nf, double initOrigTime, std::vector<Calibration *> clb, bool rnp);
    
									~FossilGraph(void);
    
    FossilGraph						&operator=(const FossilGraph &fg);
    void							clone(const FossilGraph &fg);
    double							update(double &oldLnL);
    double							lnPrior();
    void							print(std::ostream & o) const;
    std::string						writeParam(void);
    int                             getNumFossils(void) { return numFossils; }
    
    double                          getFGFBDLogProbability(void); // don't think this will be used for anything
    
    // maybe neccessary
    
    double							getFossilGraphProb(double lambda, double mu, double fossRate, double sppSampRate); // cf getTreeAncCalBDSSTreeNodePriorProb
    double							getFossilGraphProb(double lambda, double mu, double fossRate, double sppSampRate, double ot); // cf getTreeAncCalBDSSTreeNodePriorProb
    double							getActiveFossilGraphProb(); // cf getTreeCBDNodePriorProb, called by getTreeSpeciationProbability
    
    double							bdssC1Fxn(double b, double d, double psi);
    double							bdssC2Fxn(double b, double d, double psi,double rho);
    
    double							bdssQFxn(double b, double d, double psi, double rho, double t); // on log scale
    double							bdssP0Fxn(double b, double d, double psi, double rho, double t);
//    double							bdssP0HatFxn(double b, double d, double rho, double t); // not required for the fofbd
    double							fbdQHatFxn(double b, double d, double psi, double rho, double t);
    
    
    double                          getOldestFossilGraphSpeciationTime(void);
    int								pickRandAncestorFossil(void);
    int                             pickRandTipFossil(void);
    
    int								getSumIndicatorFG(void);
    
    std::string						getOccInfoParamNames(void);
    std::string						getOccInfoParamList(void);
	double							getCurrentFossilGraphLnL(void) { return currentFossilGraphLnL; }

    
private:
 
    void                            createOccurrenceVector(std::vector<Calibration *> clb);
    void							initializeOccurrenceSpecVariables(); // cf initializeFossilSpecVariables
    void							recountOccurrenceAttachNums(); // cf recountFossilAttachNums();

    double							updateOccurrenceAttachmentTimesPhi(void); //
    double							doAScaleMove(double &nv, double cv, double tv, double lb, double hb, double rv); //c.f Tree::doAScaleMove, in this case hb is not node depth, but the origin time
    double							getSumLogAllAttachNums(void);
    
    double							updateRJMoveAddDelEdge(void);
    void							doAddEdgeMove(int k);
    void							doDeleteEdgeMove(int k);
	
	double							doSinglePhiMove(void);
    
    
    int                             numFossils;
    int                             numAncFossilsk;
    double                          originTime;
    double                          terminalTime;
    int								treeTimePrior; // this will always be > 9, but this class might also be used for other options 10+
    std::vector<Occurrence *>		occurrenceSpecimens;
    double							tuningVal;
	double							currentFossilGraphLnL;
	bool							runUnderPrior;

    
    
};

#endif

