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

#ifndef PARAMETER_SPECIATION_SKYLINE_H
#define PARAMETER_SPECIATION_SKYLINE_H

#include <vector>

class MbRandom;
class Model;
class FossilRangeGraphSkyline;
class SpeciationSkyline : public Parameter {
    
public:
                            SpeciationSkyline(MbRandom *rp, Model *mp, int ni, double rh);
                            ~SpeciationSkyline(void);
    SpeciationSkyline       &operator=(const SpeciationSkyline &c);
    void                    clone(const SpeciationSkyline &c);
    double                  update(double &oldLnL);
    void                    print(std::ostream & o) const;
    double                  lnPrior(void);
    std::string             writeParam(void);
    double                  getLnFossilRangeGraphSkylineProb(FossilRangeGraphSkyline *frgsl);
    
    double                  getSppSampRateRho(void) { return extantSampleRate; }
    
    std::vector<double>     getSpeciationRates() { return birthRates; }
    std::vector<double>     getExtinctionRates() { return deathRates; }
    std::vector<double>     getFossilSampRates() { return fossilRates; }
    
    double                  getSpeciationRate(int i) { return birthRates[i]; }
    double                  getExtinctionRate(int i) { return deathRates[i]; }
    double                  getFossilSampRate(int i) { return fossilRates[i]; }
    
    void                    setAllBDFossParams(void);
    
private:

    double                  extantSampleRate;
    std::vector<double>		birthRates;
    std::vector<double>		deathRates;
    std::vector<double>		fossilRates;
    std::vector<double>		netDiversifications;
    std::vector<double>		relativeDeaths; // turnover
    std::vector<double>		probObservations;
    double                  currentFossilRangeGraphSkylineLnL;
    
    int                     parameterization; // 1=d,r,s,rho; 2=d,r,psi,rho; 3=lambda,mu,psi,rho;
    
    int                     numIntervals;
    void                    initializeIntervalVariables();
    void                    printInitialIntervalVariables();
    
    // add moves
    
    double                  updateFossilRangeGraphSkylineBDParams(double &oldLnL);
    double                  updateBirthRate(FossilRangeGraphSkyline *frg, int i);
    double                  updateDeathRate(FossilRangeGraphSkyline *frg, int i);
    double                  updatePsiRate(FossilRangeGraphSkyline *frg, int i);
    
    int                     birthRatePrior;
    int                     deathRatePrior;
    int                     fossRatePrior;
    double                  birthRateExpRate;
    double                  deathRateExpRate;
    double                  fossRateExpRate;

};

#endif

