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

class MbRandom;
class Model;
class Calibration;

class FossilRangeGraphSkyline : public Parameter {
    
public:
    FossilRangeGraphSkyline(MbRandom *rp, Model *mp, int nf, int nl, std::vector<Calibration *> clb, int ni);

    ~FossilRangeGraphSkyline(void);
    
    FossilRangeGraphSkyline			&operator=(const FossilRangeGraphSkyline &frgsl);
    void							clone(const FossilRangeGraphSkyline &frgsl);
    double                          update(double &oldLnL);
    double                          lnPrior();
    void							print(std::ostream & o) const;
    std::string                     writeParam(void);
    
    double							getFossilRangeGraphProbSkyline(std::vector<double> lambda, double mu); //, double fossRate, double sppSampRate, double ot);
    
    // probability functions
    // constants
    //
    // fbdSkylineAfxn
    // fbdSkylineBfxn
    // p, q and q tilda

private:
    
    int                             numFossils;
    int                             numLineages;
    int                             numIntervals;
    
    double							currentFossilRangeGraphSkylineLnL;
    
};


#endif
