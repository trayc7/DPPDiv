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
    Interval(double start, double end, int fossils, int intid) : intervalStart(start), intervalEnd(end), intervalFossils(fossils), intervalID(intid) {}

    double          getIntervalStart(void){ return intervalStart; }
    double          getIntervalEnd(void){ return intervalEnd; }
    int			    getIntervalFossils(void){ return intervalFossils; }
    int			    getIntervalID(void){ return intervalID; }
    
    void            setIntervalStart(double d) { intervalStart = d; };
    void            setIntervalEnd(double d) { intervalEnd = d; };
    
private:
    double          intervalStart;
    double          intervalEnd;
    int             intervalFossils;
    int             intervalID;
};

class FossilRangeSkyline {
    
public:
    FossilRangeSkyline(double fa, double la, double at, bool e, bool eo, int frid) :firstAppearance(fa), lastAppearance(la), attachmentTime(at), extant(e), extantOnly(eo), fossilRangeID(frid), fossilBrGamma(0) {}
    
    bool                              getIsExtant(void) { return extant; }
    bool                              getIsExtantOnly(void) { return extantOnly; }
    
    double                            getFirstAppearance(void) { return firstAppearance; }
    double                            getLastAppearance(void) { return lastAppearance; }
    double                            getAttachmentTime(void) { return attachmentTime; } // only used when the FRG is fixed
    double                            getLineageStart(void) { return lineageStart; }
    double                            getLineageStop(void) { return lineageStop; }
    int                               getFossilRangeBrGamma(void) { return fossilBrGamma; }
    
    void                              setLineageStop(double d) { lineageStop = d; }
    void                              setLineageStart(double d) { lineageStart = d; }
    void                              setFossilRangeBrGamma(int i) { fossilBrGamma = i; }
    
    // skyline parameters
    int                               getFossilRangeBirthInterval(void) { return birthInterval; }
    int                               getFossilRangeDeathInterval(void) { return deathInterval; }
    int                               getFossilRangeFirstAppearanceInterval(void) { return firstAppearanceInterval; }
    
    void                               setFossilRangeBirthInterval(int i) { birthInterval = i; }
    void                              setFossilRangeDeathInterval(int i) { deathInterval = i; }
    void                              setFossilRangeFirstAppearanceInterval(int i) { firstAppearanceInterval = i; }
    
    //rw: fxns required for cloning only (double check)
    
    int								getFossilRangeIndex(void) { return indx; }
    int								getFossilRangeID(void) { return fossilRangeID; } //rw: what does getFossilIndex do?
    
    void							setFossilRangeIndex(int i) { indx = i; }
    void							setFirstAppearance(double d) { firstAppearance = d; }
    void							setLastAppearance(double d) { lastAppearance = d; }
    void                            setExtant(bool b) { extant = b; }
    void                            setExtantOnly(bool b) { extantOnly = b; }
    void							setFossilRangeID(int i) { fossilRangeID = i; }
    void                            setFixStart(bool b) { fixStart = b; }
    void                            setFixStop(bool b) { fixStop = b; }
    
    bool                             getIsFixStart(void) { return fixStart; }
    bool                             getIsFixStop(void) { return fixStop; }
    
    //    int								getFossilIndicatorVar(void) { return ancFossIndicator; }
    
private:
    int								indx;
    double							firstAppearance; // oi
    double							lastAppearance; // yi
    double							attachmentTime; // note this is for fixed ranges
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
    
    //    int								ancFossIndicator; // {\cal I} = 0 if anc fossil, 1 otherwise
    
};

class MbRandom;
class Model;
class Calibration;

class FossilRangeGraphSkyline : public Parameter {
    
public:
    FossilRangeGraphSkyline(MbRandom *rp, Model *mp, int nf, int nl, std::vector<Calibration *> clb, int ni, std::vector<Calibration *> ints, bool rnp, bool fxFRG);

    ~FossilRangeGraphSkyline(void);
    
    FossilRangeGraphSkyline			&operator=(const FossilRangeGraphSkyline &frgsl);
    void							clone(const FossilRangeGraphSkyline &frgsl);
    double                          update(double &oldLnL);
    double                          lnPrior();
    void							print(std::ostream & o) const;
    std::string                     writeParam(void);
    
    double							getFossilRangeGraphProbSkyline(std::vector<double> lambda, double mu); //, double fossRate, double sppSampRate, double ot);
    
    // get functions
    
    // probability functions
    // constants
    //
    // fbdSkylineAfxn
    // fbdSkylineBfxn
    // p, q and q tilda

private:
    
    void                            createIntervalsVector(std::vector<Calibration *> ints);
    
    // you need both of these
    //void                            createFossilRangeVector(std::vector<Calibration *> clb);
    //void							initializeFossilRangeVariables(); // cf initializeFossilSpecVariables or initializeOccurrenceSpecVariables
    
    int                             numFossils;
    int                             numLineages;
    int                             numIntervals;
    
    bool							runUnderPrior;
    bool                            fixFRG;
    
    std::vector<Interval *>		intervals;
    void                            printIntervalVariables();
    
    double							currentFossilRangeGraphSkylineLnL;
    
};


#endif
