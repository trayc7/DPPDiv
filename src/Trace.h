//
//  Trace.hpp
//  dppdiv
//
//  Created by Rachel Warnock on 23/03/2020.
//  Copyright © 2020 Rachel Warnock. All rights reserved.
//

#ifndef TRACE_H
#define TRACE_H

#include <string>
#include <vector>
#include <iostream>


class MbRandom;
class Model;
class FossilRangeGraphSkyline;

class Trace {
    
public:
                        Trace(MbRandom *rp, Model *mp, std::string fn);

private:
    
    std::vector< std::vector<double> > table;
    std::vector< std::string > parameters;
    
    // //todo: don't think you need any of this
    
    // FBDR skyline model parameters
    std::vector< double > birthRates;
    std::vector< double > deathRates;
    std::vector< double > fossilRates;
    std::vector< double > tipRates;
    double rho;

    int numFossils; //k
    int numLineages;
    int numIntervals;
    
    MbRandom        *ranPtr;
    Model            *modelPtr;
    
    double                            getFossilRangeGraphSkylineProb();
    
};

#endif
