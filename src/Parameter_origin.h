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

#ifndef PARAMETER_ORIGIN_H
#define PARAMETER_ORIGIN_H

class MbRandom;
class Model;
class ExpCalib;
class FossilGraph; 
class Tree;
class Speciation;
class OriginTime : public Parameter {
    
public:
    OriginTime(MbRandom *rp, Model *mp, double sv, double yb, double ob, bool sky=false);
    ~OriginTime(void);
    OriginTime			&operator=(const OriginTime &c);
    void				clone(const OriginTime &c);
    double				update(double &oldLnL);
    double              getFBDProbOriginTime(Tree *t, Speciation *s);
    double              getFOFBDProbOriginTime(FossilGraph *fg, Speciation *s);
	double				getFossilGraphLnLikelihood(FossilGraph *fg, Speciation *s);
    void				print(std::ostream & o) const;
    double				lnPrior(void);
    std::string			writeParam(void);
    double				getOriginTime() { return originTime; }
    void				setOriginTime(double s) { originTime = s; }
    double				getLnTreeProb(Tree *t);
    double              lnExpOriginTimePriorRatio(double nOT, double oOT, double offSt, double expRate);
    
private:
    double				originTime;
    double				oldBound;
    double				yngBound;
    double				tuning;
    bool				isBounded;
    bool				retune;
    double				numAccepted;
    double				numTried;
	int					treeTimePrior;
	double				currentFossilGraphLnL;
	


	
	
	double				updateFOFBD(void);
	double				updateDPPDiv(void);
	double				updateSkylineDPPDiv(void);
	
    int                 otProposal;
    int                 otPrior;
    bool                moveOnDiff;
    double              expRate;
	
	bool				isSkyline;
    
};

#endif