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

#ifndef PARAMETER_FOSSILRANGEGRAPH_H
#define PARAMETER_FOSSILRANGEGRAPH_H

#include <vector>
#include <set>
#include <string>
#include <ostream>

class FossilRange {
    
public:
    FossilRange(double fa, double la, bool e, bool eo, int frid) :firstAppearance(fa), lastAppearance(la), extant(e), extantOnly(eo), fossilRangeID(frid), fossilBrGamma(0) {}
    
        bool                              getIsExtant(void) { return extant; }
        bool                              getIsExtantOnly(void) { return extantOnly; }
    
        double                            getFirstAppearance(void) { return firstAppearance; }
        double                            getLastAppearance(void) { return lastAppearance; }
        double                            getLineageStart(void) { return lineageStart; }
        double                            getLineageStop(void) { return lineageStop; }
        int                               getFossilRangeBrGamma(void) { return fossilBrGamma; }
    
        void                              setLineageStop(double d) { lineageStop = d; }
        void                              setLineageStart(double d) { lineageStart = d; }
        void                              setFossilRangeBrGamma(int i) { fossilBrGamma = i; }
    
        //rw: fxns required for cloning only (double check)
    
        int								getFossilRangeIndex(void) { return indx; }
        int								getFossilRangeID(void) { return fossilRangeID; } //rw: what does getFossilIndex do?

        void							setFossilRangeIndex(int i) { indx = i; }
        void							setFirstAppearance(double d) { firstAppearance = d; }
        void							setLastAppearance(double d) { lastAppearance = d; }
        void                            setExtant(bool b) { extant = b; }
        void                            setExtantOnly(bool b) { extantOnly = b; }
        void							setFossilRangeID(int i) { fossilRangeID = i; }
    
    //    int								getFossilIndicatorVar(void) { return ancFossIndicator; }

private:
    int								indx;
    double							firstAppearance; // yf
    double							lastAppearance;
    bool                            extant;
    bool                            extantOnly;
    int                             fossilRangeID;
    int								fossilBrGamma;
    double                          lineageStart; // note this is the speciation time = zf or phi
    double                          lineageStop;
    
//    int								ancFossIndicator; // {\cal I} = 0 if anc fossil, 1 otherwise

};

class MbRandom;
class Model;
class Calibration;

class FossilRangeGraph : public Parameter {
    
public:
    FossilRangeGraph(MbRandom *rp, Model *mp, int nf, int nl, std::vector<Calibration *> clb, bool rnp);
    
    ~FossilRangeGraph(void);
    
    FossilRangeGraph				&operator=(const FossilRangeGraph &frg);
    void							clone(const FossilRangeGraph &frg);
    double							update(double &oldLnL);
    double							lnPrior();
    void							print(std::ostream & o) const;
    std::string						writeParam(void);
    
//    int                             getNumFossils(void) { return numFossils; }
    
    double							getFossilRangeGraphProb(double lambda, double mu, double fossRate, double sppSampRate, double ot); // cf getTreeAncCalBDSSTreeNodePriorProb or getFossilGraphProb
    double							getActiveFossilRangeGraphProb();
    
    double                          getFossilRangeGraphOriginTime(void) { return originTime; }
    
    std::string						getFossilRangeInfoParamNames(void);
    std::string						getFossilRangeInfoParamList(void);
    
    double							fbdC1Fxn(double b, double d, double psi);
    double							fbdC2Fxn(double b, double d, double psi,double rho);
    double							fbdC3Fxn(double b, double d, double psi,double rho);
    double							fbdC4Fxn(double b, double d, double psi,double rho);
    double							fbdPFxn(double b, double d, double psi, double rho, double t);
    double							fbdQTildaFxn(double b, double d, double psi, double rho, double t); // on log scale?
    double							fbdQTildaFxnLog(double b, double d, double psi, double rho, double t);

//    double							fbdQHatFxn(double b, double d, double psi, double rho, double t);
    
//    int								getSumIndicatorFG(void);


//    double							getCurrentFossilGraphLnL(void) { return currentFossilGraphLnL; }

private:
    
    void                            createFossilRangeVector(std::vector<Calibration *> clb);
    void							initializeFossilRangeVariables(); // cf initializeFossilSpecVariables or initializeOccurrenceSpecVariables
    void							recountFossilRangeAttachNums(); // cf recountFossilAttachNums() or recountOccurrenceAttachNums()
    void                            printFossilRangeVariables(); //rw: for debugging
    void                            redefineOriginTime();
    double                          updateLineageStartTimes();
    double                          updateLineageStopTimes();
    //double                          getFossilRangeGraphProb();
    
    int                             numFossils;
    int                             numLineages;
    double                          originTime;
    double                          ancientBound;
    bool							runUnderPrior;
    bool							printInitialFossilRangeVariables;
    int                             moves; //rw:
    int                             proposal; //rw:
    std::vector<FossilRange *>		fossilRanges;
    
    double							currentFossilRangeGraphLnL; //rw:
    double                          getNewValSWindoMv(double ov, double vmin, double vmax, double tv);
    double							doAScaleMove(double &nv, double cv, double tv, double lb, double hb, double rv); //c.f FossilGraph::doAScaleMove
    
//    double							doSinglePhiMove(void);
//    int                               numAncFossilsk;
//    int								treeTimePrior; // this should always be 10
    
};



#endif
