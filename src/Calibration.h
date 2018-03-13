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


#ifndef CALIBRATION_H
#define CALIBRATION_H

#include <string>
#include <vector>


class Calibration {
	private:
		std::string		txn1, txn2;
		double			youngtime, oldtime;
        double			firstAppearance, lastAppearance, attachmentTime, endTime;
        double			intervalStart, intervalEnd, intervalProxy;
        int             intervalFossils;
		int				nodeIDX;
		bool			isRootCal;
		int				prDistType; // 1 = Uniform, 2 = offset Exp, 3 = TGS
		double			exponRate, exponMean;
        bool            isStem;
        bool            isOrigin;
		bool			isSampledAge;
        bool            isFixedAge; // for -x in the calibration file
        bool            isExtant; // fossil range taxa that are also extant
        bool            isExtantOnly; // lineages only sampled at the present
		
	public:
						Calibration(std::string calstr, int tip);
                        Calibration(std::string calstr, std::string pastr, int intervals); // overload fxn for reading presence absence file
		std::string		getTxN1() { return txn1; }
		std::string		getTxN2() { return txn2; }
		double			getYngTime() { return youngtime; }
		double			getOldTime() { return oldtime; }
		void			setNodeIndex(int i) { nodeIDX = i; }
		int				getNodeIndex() { return nodeIDX; }
		void			setIsRootCalib(bool b) { isRootCal = b; }
		bool			getIsRootCalib() { return isRootCal; }
		int				getPriorDistributionType() { return prDistType; }
		double			getCalExponRate() { return exponRate; }
		double			getCalExponMean() { return exponMean; }
        void			setIsStemFossil(bool b) { isStem = b; }
        bool			getIsStemFossil() { return isStem; }
        void			setIsOriginFossil(bool b) { isOrigin = b; }
        bool			setIsOriginFossil() { return isOrigin; }
		bool			getIsSampledAge(void) { return isSampledAge; }
        bool            getIsFixedAge(void) { return isFixedAge; }
		
		void			initializeNodeCalibration(std::string calstr);
		void			initialzeTipCalibration(std::string calstr);
        void            initializeOccurrence(std::string calstr);
        void            initializeFixedNodeAge(std::string calstr);
        void            initializeFossilRange(std::string calstr);
        void            initializeInterval(std::string calstr);
        void            initializePresenceAbsenceRange(std::string calstr, std::string pastr, int intervalNum);
    
        // for the fossil range graph model (constant rate and skyline)
        double			getFirstAppearance() { return firstAppearance; }
        double			getLastAppearance() { return lastAppearance; }
        double			getAttachmentTime() { return attachmentTime; }
        double          getEndTime() { return endTime; }
        bool			getIsExtant(void) { return isExtant; }
        bool            getIsExtantOnly(void) { return isExtantOnly; }
    
       // for the intervals
       //void             setIntervalEnd(double d) { intervalEnd = d; }
       double			getIntervalStart() { return intervalStart; }
       double			getIntervalEnd() { return intervalEnd; }
       double           getIntervalProxy() { return intervalProxy; }
       int              getIntervalFossils() { return intervalFossils; }
    
      // for presence absence data
      std::vector<int>  kappa;
      std::vector<int>  getKappa() { return kappa; }
    
};


#endif
