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

/*

The class for the parameters of the conditioned birth-death process model (Gernhard, 2008)
The code for priors on node times using this model is all based on the implementation in BEAST (https://code.google.com/p/beast-mcmc/)

This also defines parameters for the FBD model (Heath, Huelsenbeck, Stadler. unpub)

*/

#ifndef PARAMETER_SKYLINE_H
#define PARAMETER_SKYLINE_H

#include <vector>

class MbRandom;
class Model;
class FossilGraph;
class Tree;
class Skyline : public Parameter {
	
	public:
								Skyline(MbRandom *rp, Model *mp, int nrts, bool tree);
								~Skyline(void);
		Skyline					&operator=(const Skyline &c);
		void					clone(const Skyline &c);
		double					update(double &oldLnL);
		void					print(std::ostream & o) const;
		double					lnPrior(void);
		std::string				writeParam(void);
		
		double					getSkylineRhoVal(void) { return rho; }
		std::vector<double>		getSkylineBirthVec(void) { return lambdas; }
		std::vector<double>		getSkylineDeathVec(void) { return mus; }
		std::vector<double>		getSkylinePsiVec(void) { return psis; }
		
		std::string				writeSkylineParamLabels(void);
		std::string				writeSkylineParamValues(void);
		
		void					setAllSkyFBDParameters(void);
		
		
	private:
		double					rho;
		int						numRates;
		std::vector<double>		lambdas;
		std::vector<double>		mus;
		std::vector<double>		psis;
		std::vector<double>		netDivVec;
		std::vector<double>		turnoverVec;
		int						parameterization;
		
		double					rateMax;
		double					netDivExpPriorRt;
		double					psiExpPriorRt;
		std::vector<double>		updateProbs;
		
		bool					isTreeModel;
		
		void					setUpdateProbs(void);
		
		double					getNewValScaleMv(double &nv, double ov, double vmin, double vmax, double tv);
		double					getNewValSWindoMv(double &nv, double ov, double vmin, double vmax, double tv);
		double					getNewValUpDownScaleMv(double &nv1, double ov1, double &nv2, double ov2, double sf);
		
		double					updateTreeSkyParams(double &oldLnL);
		double					udateTurnovers(Tree *t);
		double					updateVectorScaleTurnover(Tree *t);
		double					udateNetDivs(Tree *t);
		double					updateVectorScaleNetDivs(Tree *t);
		double					udatePsis(Tree *t);
		double					updateVectorScalePsis(Tree *t);
		double					updateRandTurnover(Tree *t);
		double					updateRandNetDiv(Tree *t);
		double					updateRandPsi(Tree *t);
		
		double					getExponentialPriorRatio(double ov, double nv, double r);
		
		size_t					pickUpdate(void);
		double					getTreeSpeciationProb(Tree *t);
		void					initializeAsEqual(void);
};


#endif