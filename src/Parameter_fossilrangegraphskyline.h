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

#ifndef PARAMETER_FOSSILRANGEGRAPHSKYLINE_H
#define PARAMETER_FOSSILRANGEGRAPHSKYLINE_H

#include <vector>
#include <set>
#include <string>
#include <ostream>

class Interval {
    
public:
    Interval(double start, double end, int fossils, int intid, double durations, int kP) : intervalStart(start), intervalEnd(end), intervalFossils(fossils), intervalID(intid),
    intervalSumRangeLengths(durations), intervalKappaPrime(kP) {}

    double          getIntervalStart(void){ return intervalStart; }
    double          getIntervalEnd(void){ return intervalEnd; }
    int			    getIntervalFossils(void){ return intervalFossils; }
    int			    getIntervalID(void){ return intervalID; }
    double          getIntervalSumRangeLengths(void) { return intervalSumRangeLengths; }
    int             getIntervalKappaPrime(void) {return intervalKappaPrime; }
    
    void            setIntervalStart(double d) { intervalStart = d; };
    void            setIntervalEnd(double d) { intervalEnd = d; };
    void            setIntervalSumRangeLengths(double d) { intervalSumRangeLengths = d; }
    void            setIntervalKappaPrime(int i) { intervalKappaPrime = i; }
    
private:
    double          intervalStart;
    double          intervalEnd;
    int             intervalFossils;
    int             intervalID;
    double          intervalSumRangeLengths;
    int             intervalKappaPrime;
};

class FossilRangeSkyline {
    
public:
    FossilRangeSkyline(double fa, double la, double at, double et, bool e, bool eo, int frid, std::vector<bool> gi) :firstAppearance(fa), lastAppearance(la), attachmentTime(at), endTime(et), extant(e), extantOnly(eo), fossilRangeID(frid), fossilBrGamma(0), gammaInteractions(gi) {}
    
    bool                              getIsExtant(void) { return extant; }
    bool                              getIsExtantOnly(void) { return extantOnly; }
    
    double                            getFirstAppearance(void) { return firstAppearance; }
    double                            getLastAppearance(void) { return lastAppearance; }
    double                            getAttachmentTime(void) { return attachmentTime; } // only used when the FRG is fixed
    double                            getEndTime(void) {return endTime; } // only used when the FRG is fixed
    double                            getLineageStart(void) { return lineageStart; }
    double                            getLineageStop(void) { return lineageStop; }
    int                               getFossilRangeBrGamma(void) { return fossilBrGamma; }
    
    bool                              getIsFixStart(void) { return fixStart; }
    bool                              getIsFixStop(void) { return fixStop; }
    
    void                              setLineageStop(double d) { lineageStop = d; }
    void                              setLineageStart(double d) { lineageStart = d; }
    void                              setFossilRangeBrGamma(int i) { fossilBrGamma = i; }
    
    // skyline parameters
    int                               getFossilRangeBirthInterval(void) { return birthInterval; }
    int                               getFossilRangeDeathInterval(void) { return deathInterval; }
    int                               getFossilRangeFirstAppearanceInterval(void) { return firstAppearanceInterval; }
    
    void                              setFossilRangeBirthInterval(int i) { birthInterval = i; }
    void                              setFossilRangeDeathInterval(int i) { deathInterval = i; }
    void                              setFossilRangeFirstAppearanceInterval(int i) { firstAppearanceInterval = i; }
    
    // fxns required for cloning only
    int								getFossilRangeIndex(void) { return indx; }   // this isn't used
    int								getFossilRangeID(void) { return fossilRangeID; }
    
    void							setFossilRangeIndex(int i) { indx = i; }
    void							setFirstAppearance(double d) { firstAppearance = d; }
    void							setLastAppearance(double d) { lastAppearance = d; }
    void                            setExtant(bool b) { extant = b; }
    void                            setExtantOnly(bool b) { extantOnly = b; }
    void							setFossilRangeID(int i) { fossilRangeID = i; }
    void                            setFixStart(bool b) { fixStart = b; }
    void                            setFixStop(bool b) { fixStop = b; }
    
    //fxns required for gamma trick
    bool                            getFossilRangeGammaInteractions(int i) { return gammaInteractions[i]; }
    void                            setFossilRangeGammaInteractions(int i, bool b) { gammaInteractions[i] = b; }
    
private:
    int								indx;
    double							firstAppearance; // oi
    double							lastAppearance; // yi
    double							attachmentTime; // note this is for fixed ranges
    double							endTime; // note this is for fixed ranges
    bool                            extant;
    bool                            extantOnly;
    int                             fossilRangeID;
    int								fossilBrGamma;
    double                          lineageStart; // note this is the speciation time = zf or phi
    double                          lineageStop;
    bool                            fixStart;
    bool                            fixStop;
    
    // skyline parameters
    int                             birthInterval;
    int                             deathInterval;
    int                             firstAppearanceInterval;
    
    std::vector<bool>               gammaInteractions;
    
};

class MbRandom;
class Model;
class Calibration;

class FossilRangeGraphSkyline : public Parameter {
    
public:
    FossilRangeGraphSkyline(MbRandom *rp, Model *mp, int nf, int nl, std::vector<Calibration *> clb, int ni, std::vector<Calibration *> ints, bool rnp, bool fxFRG, int expMode, int fbdLk);

    ~FossilRangeGraphSkyline(void);
    
    FossilRangeGraphSkyline			&operator=(const FossilRangeGraphSkyline &frgs);
    void							clone(const FossilRangeGraphSkyline &frgs);
    double                          update(double &oldLnL);
    double                          lnPrior();
    void							print(std::ostream & o) const;
    std::string                     writeParam(void);
    
    //double							getFossilRangeGraphSkylineProb(std::vector<double> lambda, std::vector<double> mu, std::vector<double> fossRate, std::vector<double> sppSampRate, double ot);
    double							getFossilRangeGraphSkylineProb();
    double							getFossilRangeGraphProb(std::vector<double> b, std::vector<double> d, std::vector<double> s, std::vector<double> r, double ot); // for cross validation
    double							getActiveFossilRangeGraphSkylineProb();
    
    double                          getFossilRangeGraphSkylineOriginTime(void) { return originTime; }
    int                             getNumFossils(void) { return numFossils; }
    int                             getNumIntervals(void) { return numIntervals; }
    
    std::string						getFossilRangeSkylineInfoParamNames(void);
    std::string						getFossilRangeSkylineInfoParamList(void);
    
    // frg skyline functions
    void                            setAllIntervalConstants(void);
    double                          fbdSkylineABPfxnInterval(std::vector<double> b, std::vector<double> d, std::vector<double> psi, std::vector<double> rho, int i, double t);
    double							fbdSkylineAfxn(std::vector<double> b, std::vector<double> d, std::vector<double> psi, int i);
    double							fbdSkylineBfxn(std::vector<double> b, std::vector<double> d, std::vector<double> psi, std::vector<double> rho, int i);
    double                          fbdSkylinePfxn(std::vector<double> b, std::vector<double> d, std::vector<double> psi, std::vector<double> rho, int i, double t);
    double                          fbdSkylineQfxn(std::vector<double> b, std::vector<double> d, std::vector<double> psi, std::vector<double> rho, int i, double t);
    double                          fbdSkylineQfxnLog(std::vector<double> b, std::vector<double> d, std::vector<double> psi, std::vector<double> rho, int i, double t);
    double                          fbdSkylineQTildaFxn(std::vector<double> b, std::vector<double> d, std::vector<double> psi, std::vector<double> rho, int i, double t);
    double                          fbdSkylineQTildaFxnLog(std::vector<double> b, std::vector<double> d, std::vector<double> psi, std::vector<double> rho, int i, double t);
    double                          fbdSkylineQTildaFxnLogSimplified(std::vector<double> b, std::vector<double> d, std::vector<double> psi, std::vector<double> rho, int i, double t, double q);

    double                          exampleRevBayesPfxn(std::vector<double> l, std::vector<double> m, std::vector<double> psi, std::vector<double> rho, int i, double t);
    
    // frg functions
    double							fbdC1Fxn(double b, double d, double psi);
    double							fbdC2Fxn(double b, double d, double psi,double rho);
    double							fbdC3Fxn(double b, double d, double psi,double rho);
    double							fbdC4Fxn(double b, double d, double psi,double rho);
    double							fbdPFxn(double b, double d, double psi, double rho, double t);
    double							fbdQFxnLog(double b, double d, double psi, double rho, double t);
    double							fbdQTildaFxnLog(double b, double d, double psi, double rho, double t);
    double							fbdQTildaFxnLogAlt(double b, double d, double psi, double rho, double t);
    
    void                            crossValidateFBDSkylinefunctions();
    
private:
    
    void                            createIntervalsVector(std::vector<Calibration *> ints);
    void                            createFossilRangeSkylineVector(std::vector<Calibration *> clb);
    void							initializeFossilRangeSkylineVariables();
    void                            initializeIntervalConstants();
    void							recountFossilRangeAttachNums();
    void							recountFossilRangeAttachNumsSpeedy(int i);
    void                            redefineOriginTime();
    void                            countExtinctLineages();
    int                             assignInterval(double d);
    double                          updateLineageStartTimes();
    double                          updateLineageStopTimes();
    void                            orderRangeAges();
    void                            calculateIntervalSumRanges(); // L_S
    
    int                             numFossils;
    int                             numLineages;
    int                             numExtinctLineages;
    int                             numIntervals;
    int                             originInterval;
    double                          originTime;
    double                          ancientBound;
    bool							runUnderPrior;
    bool                            printInitialFossilRangeSkylineVariables;
    int                             fbdLikelihood;
    
    std::vector<Interval *>         intervals;
    std::vector<FossilRangeSkyline *>		fossilRangesSkyline;
    std::vector<double>           intervalAs;
    std::vector<double>           intervalBs;
    std::vector<double>           intervalPs;
    std::vector<double>           intervalQs;
    std::vector<double>           intervalQts;
    
    void                            printIntervalVariables();
    void                            printFossilRangeSkylineVariables(); //debugging code
    void                            printFossilRangeVariables(int range); //debugging code
    
    double							currentFossilRangeGraphSkylineLnL;
    double							doAScaleMove(double &nv, double cv, double tv, double lb, double hb, double rv);
    
//  int								treeTimePrior; // this should always be 11
    
    bool                            fixFRG;
    bool                            fixOrigin;
    bool                            fixStart;
    bool                            fixStop;
    int                             counter; // debugging
    bool                            orderStartStopTimes;
    
    bool                            speedy;
    
};


#endif
