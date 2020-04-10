//
//  Trace.cpp
//  dppdiv
//
//  Created by Rachel Warnock on 23/03/2020.
//  Copyright © 2020 Rachel Warnock. All rights reserved.
//

#include "MbRandom.h"
#include "Model.h"
#include "Calibration.h"
#include "Parameter.h"
#include "Parameter_speciationskyline.h"
#include "Parameter_fossilrangegraphskyline.h"
#include "Trace.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <istream>
#include <vector>
#include <string>

//#include <cstdlib>

// todo: which ones are needed?

using namespace std;

Trace::Trace(MbRandom *rp, Model *mp, string fn){
    
    ranPtr          = rp; /* not technically needed */
    modelPtr        = mp;

    /* code taken from stackoverflow.com/questions/2221065/ */
    std::fstream ifs;

    /*  open file  */
    ifs.open(fn);
    if(!ifs){
      cerr << "\nCan not open trace file!" << endl;
      exit(1);
    } else
      cout << "\nReading trace file..." << endl;
    
    int lineNum = 0;
    
    while (true) {
        
        std::string line;
        getline(ifs, line);

        std::stringstream ss(line, std::ios_base::out|std::ios_base::in|std::ios_base::binary);

        if (!ifs)
            // mainly catch EOF
            break;

        if (line[0] == '#' || line.empty())
            // catch empty lines or comment lines
            continue;

        if(lineNum == 0){
            
            string buf;
            
            while (ss >> buf)
                parameters.push_back(buf);
        } else {
            
            double buf;
            
            std::vector<double> row;

            while (ss >> buf)
                row.push_back(buf);

            table.push_back(row);
            
        }
        lineNum ++;
    }

    ifs.close();
    
    cout << "Finished reading trace file" << endl;
    
    // calculate how many intervals & ranges you have
    numIntervals = 0;
    numLineages = 0;
    
    for(int i = 0; i < parameters.size(); i++){
        
        string p = parameters[i];
        
        if (p.find("lambda") != string::npos)
            numIntervals ++;
        
        if (p.find("b_f") != string::npos)
            numLineages ++;
    }
    
    // throw an error if out file was generated without using the -pfat option
    if(numLineages == 0){
        cerr << "Not enough information provided in the trace file!" << endl;
        exit(1);
    }
     
    cout << "Number of intervals " << numIntervals << endl;
    cout << "Number of ranges " << numLineages << endl;
    
    SpeciationSkyline *sp = modelPtr->getActiveSpeciationSkyline();
    FossilRangeGraphSkyline *frgs = modelPtr->getActiveFossilRangeGraphSkyline();
    
    // check cal and int files are compatible with the tracefile
    if(frgs->fossilRangesSkyline.size() != numLineages || frgs->getNumIntervals() != numIntervals){
        cerr << "Incompatible input and output files!" << endl;
        exit(1);
    }
    
    for(int step = 0; step <lineNum-1; step++){
        
        // define FBD parameters based on the trace file
        for(int i = 1; i < numIntervals+1; i++){
            
            //create a string for lambda, mu & psi
            stringstream ss;
            ss << "lambda[" << i << "]";
            string lm = ss.str();
            ss.str(string());
            
            ss << "mu[" << i << "]";
            string mu = ss.str();
            ss.str(string());
            
            ss << "psi[" << i << "]";
            string ps = ss.str();
            ss.str(string());
            
            for(int x = 0; x < parameters.size(); x++){
                
                string p = parameters[x];
                
                if(p.find(lm) != string::npos){
                    sp->setSpeciationRate(table[step][x], i-1);
                }
                
                if(p.find(mu) != string::npos){
                    sp->setExtinctionRate(table[step][x], i-1);
                }
                
                if(p.find(ps) != string::npos){
                    sp->setFossilSampRate(table[step][x], i-1);
                }
                
                if(p.find("FBD.OriginTime") != string::npos && i==1){
                    frgs->setOriginTime(table[step][x]);
                    frgs->setOriginInterval(frgs->assignInterval(table[step][x]));
                }
            }
        }
        
        for(int i = 1; i < numLineages+1; i++){
            
            FossilRangeSkyline *fr = frgs->fossilRangesSkyline[i-1];
            
            //create a string for FA, bi, di & gamma
            stringstream ss;
            ss << "y_f(FR_" << i << ")";
            string fa = ss.str();
            ss.str(string());
            
            ss << "b_f(FR_" << i << ")";
            string bi = ss.str();
            ss.str(string());
            
            ss << "d_f(FR_" << i << ")";
            string di = ss.str();
            ss.str(string());
            
            ss << "gamma_f(FR_" << i << ")";
            string gm = ss.str();
            ss.str(string());
            
            for(int x = 0; x < parameters.size(); x++){
                
                string p = parameters[x];
                
                if(p.find(fa) != string::npos) {
                    fr->setFirstAppearance(table[step][x]);
                    fr->setFossilRangeFirstAppearanceInterval(frgs->assignInterval(table[step][x]));
                }
                
                if(p.find(bi) != string::npos) {
                    fr->setLineageStart(table[step][x]);
                    fr->setFossilRangeBirthInterval(frgs->assignInterval(table[step][x]));
                }
                
                if(p.find(di) != string::npos) {
                    fr->setLineageStop(table[step][x]);
                    fr->setFossilRangeDeathInterval(frgs->assignInterval(table[step][x]));
                }
                
                if(p.find(gm) != string::npos)
                    fr->setFossilRangeBrGamma(table[step][x]);
                
            }
        }
        
        sp->setAllBDFossParams();
        frgs->setAllIntervalConstants();
        
        double lk = frgs->getFossilRangeGraphSkylineProb();
        cout << "lk " << lk << endl;
        
    }
    
}
