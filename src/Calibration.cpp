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

#include "Calibration.h"
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <cstdlib>

using namespace std;

Calibration::Calibration(string calstr, int tip){
	
	if(tip == 0)
		initializeNodeCalibration(calstr);
	else if(tip == 1)
		initialzeTipCalibration(calstr);
    else if (tip==2)
        initializeOccurrence(calstr);
    else if (tip == 3)
        initializeFixedNodeAge(calstr);
    else if (tip == 4)
        initializeFossilRange(calstr);
}


void Calibration::initializeNodeCalibration(string calstr){

	isRootCal = false;
	nodeIDX = -1;
	prDistType = 1;
	exponRate = -1.0;
	exponMean = -1.0;
    isStem = false;
	isSampledAge = false;
	stringstream ss;
	string tmp = "";
	ss << calstr;
	ss >> tmp;
	if(tmp[0] == '-'){
		if(tmp[1] == 'U' || tmp[1] == 'u')
			prDistType = 1;
		else if(tmp[1] == 'E' || tmp[1] == 'e')
			prDistType = 2;
		else if(tmp[1] == 'T' || tmp[1] == 't') // fossil belongs to crown group
			prDistType = 3;
        else if(tmp[1] == 'S' || tmp[1] == 's') { // fossil belongs to total group (stem or crown)
            prDistType = 4;
            isStem=true;
        }
		else{
			cerr << "ERROR: There's a problem with the calibration file " << endl;
			exit(1);
		}
		cout << tmp[2] << endl;
		if(tmp[2] == 'E' || tmp[2] == 'e'){
			isSampledAge = true;
		}
		
		ss >> txn1;
		if(txn1 != "root")
			ss >> txn2;
		else { 
			isRootCal = true;
			txn2 = "root";
            if(prDistType == 4)
                isOrigin = true;
		}
		ss >> tmp;
		if(prDistType == 1){
			youngtime = atof(tmp.c_str());
			ss >> tmp;
			oldtime = atof(tmp.c_str());
			cout << "   Uniform calibration on MRCA[" << txn1 << ", " << txn2 << "] --> ("
			<< youngtime << ", " << oldtime << ")" << endl;
		}
		else if(prDistType == 2){
			youngtime = atof(tmp.c_str());
			ss >> tmp;
			if(tmp[0] == '-'){
				if(tmp == "-r"){
					ss >> tmp;
					exponRate = atof(tmp.c_str());
					exponMean = 1.0 / exponRate;
				}
				else if(tmp == "-m"){
					ss >> tmp;
					exponMean = atof(tmp.c_str()) - youngtime;
					exponRate = 1.0 / exponMean;
				}
			}
			else{
				if(isRootCal){
					exponRate = 1.0 / (youngtime * 0.4);
					exponMean = 1.0 / exponRate;
				}
				else{
					exponRate = 1.0 / (youngtime * 0.25);
					exponMean = 1.0 / exponRate;
				}
			}
			cout << "   Offset-exponential calibration on MRCA[" << txn1 << ", " << txn2 << "] --> (offset="
			<< youngtime << ", lambda=" << exponRate << ", mean=" << exponMean + youngtime << ")" << endl;
		}
        else if(prDistType == 3){
            youngtime = atof(tmp.c_str());
            ss >> tmp;
			if(isSampledAge){
				oldtime = atof(tmp.c_str());
			}
            cout << "   Calibrate Birth-Death Prior on MRCA[" << txn1 << ", " << txn2 << "] --> (offset="
            << youngtime << ")" << endl;
        }
        else if(prDistType == 4){
            youngtime = atof(tmp.c_str());
            ss >> tmp;
			if(isSampledAge){
				oldtime = atof(tmp.c_str());
			}
            cout << "   Calibrate Birth-Death Prior on MRCA[" << txn1 << ", " << txn2 << "] --> (offset="
            << youngtime << ")" << endl;
        }
	}
	else{
		txn1 = tmp;
		isRootCal = false;
		if(txn1 != "root")
			ss >> txn2;
		else { 
			isRootCal = true;
			txn2 = "root";
		}
		ss >> tmp;
		youngtime = atof(tmp.c_str());
		ss >> tmp;
		oldtime = atof(tmp.c_str());
		cout << "   Uniform calibration on MRCA[" << txn1 << ", " << txn2 << "] --> ("
		<< youngtime << ", " << oldtime << ")" << endl;
	}

}

void Calibration::initializeFixedNodeAge(string calstr){
    
    isRootCal = false;
    isFixedAge = false;
    
    stringstream ss;
    string tmp = "";
    ss << calstr;
    ss >> tmp;
    if(tmp[0] == '-'){
        if(tmp[1] == 'X' || tmp[1] == 'x') { // fix the age of the node
            isFixedAge=true;
        }
        else{
            cerr << "ERROR: There's a problem with the calibration file " << endl;
            exit(1);
        }
        cout << tmp[2] << endl;
        
        ss >> txn1;
        if(txn1 != "root")
            ss >> txn2;
        else {
            isRootCal = true;
            txn2 = "root";
        }
        ss >> tmp;
        
        youngtime = atof(tmp.c_str());
        ss >> tmp;
        oldtime = youngtime;
        if(isRootCal)
            cout << "   Fix node age on MRCA[" << txn1 << "] --> (age="<< youngtime << ")" << endl;
        else
            cout << "   Fix node age on MRCA[" << txn1 << ", " << txn2 << "] --> (age="<< youngtime << ")" << endl;
    }
    
}

void Calibration::initialzeTipCalibration(string calstr){

	stringstream ss;
	string tmp = "";
	ss << calstr;
	ss >> txn1;
	txn2 = txn1;
	ss >> tmp;
	youngtime = atof(tmp.c_str());
	oldtime = youngtime;
	cout << "   Tip calibration on " << txn1 << " --> ("
	<< youngtime << ", " << oldtime << ")" << endl;
}

void Calibration::initializeOccurrence(string calstr){
    
    stringstream ss;
    string tmp = "";
    ss << calstr;
    ss >> tmp;
    youngtime = atof(tmp.c_str());
    cout << "Occurrence age " << youngtime << endl;
    
}

void Calibration::initializeFossilRange(string calstr){
    
    isExtant = false;
    isExtantOnly = false;
    
    stringstream ss;
    string tmp = "";
    ss << calstr;
    
    ss >> tmp;
    lastAppearance = atof(tmp.c_str());
    ss >> tmp;
    firstAppearance = atof(tmp.c_str());
    ss >> tmp;
    attachmentTime = atof(tmp.c_str()); // consider giving this a more informative name
    
    if(firstAppearance == 0)
        isExtantOnly = true;
    if(lastAppearance == 0)
         isExtant = true;
    
    cout << firstAppearance << " - " << lastAppearance << endl;
    
}


/*
3
root	70	80
T1	T3	8	12
T8	T9	30	40

or

2
T1	T3	0.04	0.14
T8	T9	0.45	0.65

or

3
-E	root	min	mean
-U	T1	T6	min	max


*/