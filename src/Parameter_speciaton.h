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

#ifndef PARAMETER_SPECIATION_H
#define PARAMETER_SPECIATION_H

class MbRandom;
class Model;
class Tree;
class Speciation : public Parameter {
	
	public:
							Speciation(MbRandom *rp, Model *mp, double bdr, double bda, double bds, double initRH);
							~Speciation(void);
		Speciation			&operator=(const Speciation &c);
		void				clone(const Speciation &c);
		double				update(double &oldLnL);
		void				print(std::ostream & o) const;
		double				lnPrior(void);
		std::string			writeParam(void);
		double				getRelativeDeath() { return relativeDeath; }
		void				setRelativeDeath(double v) { relativeDeath = v; }
		double				getNetDiversification() { return netDiversificaton; }
		void				setNetDiversification(double v) { netDiversificaton = v; }
		double				getLnTreeProb(Tree *t);
		
		double				getBDSSSpeciationRateLambda(void) { return birthRate; }
		double				getBDSSExtinctionRateMu(void) { return deathRate; }
		double				getBDSSFossilSampRatePsi(void) { return fossilRate; }
		double				getBDSSSppSampRateRho(void) { return extantSampleRate; }
		double				getBDSSFossSampProbOmega(void) { return fossilStratSampleProb; }
		double				getBDSSFossilSampProbS(void) { return probSpeciationS; }
		void				setAllBDFossParams(void);
		
	private:
		double				relativeDeath;		// mu / lambda
		double				netDiversificaton; // lambda - mu
		int					treeTimePrior; // 1=unif, 2=yule, 3=cbd, 4=ssbd, 5=ssbd2, 6=calib ssbd
		double				maxdivV;
		double				birthRate;  // lambda
		double				deathRate;  // mu
		double				fossilRate; // psi
		double				extantSampleRate; // rho = 1
		double				fossilStratSampleProb; // omega = ? This is for later
		double				probSpeciationS;  // \psi / (\mu + \psi) // Need hyperprior
		
		double				getNewValScaleMv(double &nv, double ov, double vmin, double vmax, double tv);
		double				getNewValSWindoMv(double ov, double vmin, double vmax, double tv);


		
		double				updateRelDeathRt();
		double				updateNetDivRate();
		double				updateBDSSFossilProbS();
		double				updateBDSSSampleProbRho();

		double				updateRelDeathRt(Tree *t);
		double				updateNetDivRate(Tree *t);
		double				updateBDSSFossilProbS(Tree *t);
		double				updatePsiRate(Tree *t);

};


#endif