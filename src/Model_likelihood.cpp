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

#include "cpuspec.h"
#include "Alignment.h"
#include "MbRandom.h"
#include "MbTransitionMatrix.h"
#include "Model.h"
#include "Parameter.h"
#include "Parameter_basefreq.h"
#include "Parameter_exchangeability.h"
#include "Parameter_expcalib.h"
#include "Parameter_rate.h"
#include "Parameter_shape.h"
#include "Parameter_speciaton.h"
#include "Parameter_tree.h"
#include "Parameter_cphyperp.h"
#include "Parameter_treescale.h"
#include "Calibration.h"
#include "util.h"
#include <string>
#include <vector>
#include <fstream>



#include <omp.h>


using namespace std;
double Model::lnLikelihood(void) {
        
        double * clP;
        double * clL;
        double * clR;

	
	if(runUnderPrior){
		myCurLnL = 0.0;
		return 0.0;
	}

	Tree *t = getActiveTree();
	MbMatrix<double> *tL = new MbMatrix<double>[numGammaCats];
	MbMatrix<double> *tR = new MbMatrix<double>[numGammaCats];

	for (int n=0; n<t->getNumNodes(); n++) {
		Node *p = t->getDownPassNode(n);
		if (p->getLft() != NULL && p->getRht() != NULL && p->getIsClDirty() == true) {
			clL = clPtr[p->getLft()->getActiveCl()][p->getLft()->getIdx()];
			clR = clPtr[p->getRht()->getActiveCl()][p->getRht()->getIdx()];
			clP = clPtr[p->getActiveCl()          ][p->getIdx()          ];
			for (int k=0; k<numGammaCats; k++) {
				tL[k] = tis[p->getLft()->getActiveTi()][p->getLft()->getIdx()][k];
				tR[k] = tis[p->getRht()->getActiveTi()][p->getRht()->getIdx()][k];
			}

// parallelisation
			#pragma omp parallel for
			for (int c=0; c<numPatterns; c++) {
// parallelisation                      
                                int p = c * numGammaCats * 4;
				for (int k=0; k<numGammaCats; k++) {

					double sumL = 0.0, sumR = 0.0;
#ifdef _TOM_SSE3
                                        #ifdef _ASM_DEBUG
                                        __asm__ ( "int $0x3" );
                                        #endif
                                        __m128d
                                                cll0, clr0,
                                                cll2, clr2,
                                                p1, p2,
                                                s1, s2,
                                                sr, sl;

                                        cll0 = _mm_load_pd ( clL + p );
                                        clr0 = _mm_load_pd ( clR + p );
                                        cll2 = _mm_load_pd ( clL + p + 2 );
                                        clr2 = _mm_load_pd ( clR + p + 2 );

                                        /* Compute clP[p + 0] and clP[p + 1] */
                                        p1 = _mm_mul_pd ( _mm_load_pd ( tL[k][0] ), cll0 );       // tL[k][0][0] * clL[p + 0], tL[k][0][1] * clL[p + 1]
                                        p2 = _mm_mul_pd ( _mm_load_pd ( tL[k][1] ), cll0 );       // tL[k][1][0] * clL[p + 0], tL[k][1][1] * clL[p + 1]
                                        s1 = _mm_hadd_pd ( p1, p2 );                                                     // tL[k][0][0] * clL[p + 0] + tL[k][0][1] * clL[p + 1], tL[k][1][0] * clL[p + 0] + tL[k][1][1] * clL[p + 1]
                                        
                                        
                                        p1 = _mm_mul_pd ( _mm_load_pd ( tL[k][0] + 2 ), cll2 );   // tL[k][0][2] * clL[p + 2], tL[k][0][3] * clL[p + 3]
                                        p2 = _mm_mul_pd ( _mm_load_pd ( tL[k][1] + 2 ), cll2 );   // tL[k][1][2] * clL[p + 2], tL[k][1][3] * clL[p + 3]
                                        s2 = _mm_hadd_pd ( p1, p2 );                                                     // tL[k][0][2] * clL[p + 2] + tL[k][0][3] * clL[p + 3], tL[k][1][2] * clL[p + 2] + tL[k][1][3] * clL[p + 3]
                                        
                                        /*  tL[k][0][0] * clL[p + 0] + tL[k][0][1] * clL[p + 1] + tL[k][0][2] * clL[p + 2] + tL[k][0][3] * clL[p + 3],
                                         *  tL[k][1][0] * clL[p + 0] + tL[k][1][1] * clL[p + 1] + tL[k][1][2] * clL[p + 2] + tL[k][1][3] * clL[p + 3]
                                         */
                                        sl = _mm_add_pd ( s1, s2 );                                                  
                                        
                                        p1 = _mm_mul_pd ( _mm_load_pd ( tR[k][0] ), clr0 );
                                        p2 = _mm_mul_pd ( _mm_load_pd ( tR[k][1] ), clr0 );
                                        s1 = _mm_hadd_pd ( p1, p2 );
                                        
                                        p1 = _mm_mul_pd ( _mm_load_pd ( tR[k][0] + 2 ), clr2 );
                                        p2 = _mm_mul_pd ( _mm_load_pd ( tR[k][1] + 2 ), clr2 );
                                        s2 = _mm_hadd_pd ( p1, p2 );
                                        
                                        sr = _mm_add_pd ( s1, s2 );
                                        
                                        _mm_store_pd ( clP + p, _mm_mul_pd ( sl, sr ) );
                                        
                                        /* Compute clP[p + 2] and clP[p + 3] */
                                        p1 = _mm_mul_pd ( _mm_load_pd ( tL[k][2] ), cll0 );       // tL[k][0][0] * clL[p + 0], tL[k][0][1] * clL[p + 1]
                                        p2 = _mm_mul_pd ( _mm_load_pd ( tL[k][3] ), cll0 );       // tL[k][1][0] * clL[p + 0], tL[k][1][1] * clL[p + 1]
                                        s1 = _mm_hadd_pd ( p1, p2 );
                                        
                                        p1 = _mm_mul_pd ( _mm_load_pd ( tL[k][2] + 2 ), cll2 );   // tL[k][0][2] * clL[p + 2], tL[k][0][3] * clL[p + 3]
                                        p2 = _mm_mul_pd ( _mm_load_pd ( tL[k][3] + 2 ), cll2 );   // tL[k][1][2] * clL[p + 2], tL[k][1][3] * clL[p + 3]
                                        s2 = _mm_hadd_pd ( p1, p2 );
                                        
                                        sl = _mm_add_pd ( s1, s2 );
                                        
                                        p1 = _mm_mul_pd ( _mm_load_pd ( tR[k][2] ), clr0 );
                                        p2 = _mm_mul_pd ( _mm_load_pd ( tR[k][3] ), clr0 );
                                        s1 = _mm_hadd_pd ( p1, p2 );
                                        
                                        p1 = _mm_mul_pd ( _mm_load_pd ( tR[k][2] + 2 ), clr2 );
                                        p2 = _mm_mul_pd ( _mm_load_pd ( tR[k][3] + 2 ), clr2 );
                                        s2 = _mm_hadd_pd ( p1, p2 );
                                        
                                        sr = _mm_add_pd ( s1, s2 );
                                        
                                        _mm_store_pd ( clP + p + 2, _mm_mul_pd ( sl, sr ) );
                                        #ifdef _ASM_DEBUG
                                        __asm__ ( "int $0x3" );
                                        #endif

#elif _TOM_AVX
                                        #ifdef _ASM_DEBUG
                                        __asm__ ( "int $0x3" );
                                        #endif
                                        __m256d
                                                cll,
                                                clr,
                                                r1, r2, r3, r4, r12, r34, r1234,
                                                l1, l2, l3, l4, l12, l34, l1234,
                                                r, perm, blnd;


                                       // __asm__ ("int $0x3");
//                                        printf ( "OUT %d!!\n", BYTE_ALIGNMENT );
                                        cll = _mm256_load_pd ( clL + p );
                                        clr = _mm256_load_pd ( clR + p );

                                        // Compute sumL rows

                                        l1 = _mm256_mul_pd ( _mm256_load_pd ( tL[k][0] ), cll );
                                        l2 = _mm256_mul_pd ( _mm256_load_pd ( tL[k][1] ), cll );
                                        l3 = _mm256_mul_pd ( _mm256_load_pd ( tL[k][2] ), cll );
                                        l4 = _mm256_mul_pd ( _mm256_load_pd ( tL[k][3] ), cll );

                                        l12 = _mm256_hadd_pd ( l1, l2 );
                                        l34 = _mm256_hadd_pd ( l3, l4 );

                                        blnd = _mm256_blend_pd ( l12, l34, 0b1100 );
                                        perm = _mm256_permute2f128_pd ( l12, l34, 0x21 );
                                        l1234 = _mm256_add_pd ( perm, blnd ); 

                                        // Compute sumR rows

                                        r1 = _mm256_mul_pd ( _mm256_load_pd ( tR[k][0] ), clr );
                                        r2 = _mm256_mul_pd ( _mm256_load_pd ( tR[k][1] ), clr );
                                        r3 = _mm256_mul_pd ( _mm256_load_pd ( tR[k][2] ), clr );
                                        r4 = _mm256_mul_pd ( _mm256_load_pd ( tR[k][3] ), clr );

                                        r12 = _mm256_hadd_pd ( r1, r2 );
                                        r34 = _mm256_hadd_pd ( r3, r4 );

                                        blnd = _mm256_blend_pd ( r12, r34, 0b1100 );
                                        perm = _mm256_permute2f128_pd ( r12, r34, 0x21 );
                                        r1234 = _mm256_add_pd ( perm, blnd ); 

                                        r = _mm256_mul_pd ( l1234, r1234 );


                                        _mm256_store_pd ( clP + p, r );
                                        #ifdef _ASM_DEBUG
                                        __asm__ ( "int $0x3" );
                                        #endif
#else
                                        #ifdef _ASM_DEBUG
                                        __asm__ ( "int $0x3" );
                                        #endif
					sumL += tL[k][0][0] * clL[p + 0] + tL[k][0][1] * clL[p + 1] + tL[k][0][2] * clL[p + 2] +tL[k][0][3] * clL[p + 3];
					sumR += tR[k][0][0] * clR[p + 0] + tR[k][0][1] * clR[p + 1] + tR[k][0][2] * clR[p + 2] +tR[k][0][3] * clR[p + 3];
					clP[p + 0] = sumL * sumR;
					
					sumL = 0.0, sumR = 0.0;
					sumL += tL[k][1][0] * clL[p + 0] + tL[k][1][1] * clL[p + 1] + tL[k][1][2] * clL[p + 2] + tL[k][1][3] * clL[p + 3];
					sumR += tR[k][1][0] * clR[p + 0] + tR[k][1][1] * clR[p + 1] + tR[k][1][2] * clR[p + 2] + tR[k][1][3] * clR[p + 3];
					clP[p + 1] = sumL * sumR;
					
					sumL = 0.0, sumR = 0.0;
					sumL += tL[k][2][0] * clL[p + 0] + tL[k][2][1] * clL[p + 1] + tL[k][2][2] * clL[p + 2] +tL[k][2][3] * clL[p + 3];
					sumR += tR[k][2][0] * clR[p + 0] + tR[k][2][1] * clR[p + 1] + tR[k][2][2] * clR[p + 2] +tR[k][2][3] * clR[p + 3];
					clP[p + 2] = sumL * sumR;
					
					sumL = 0.0, sumR = 0.0;
					sumL += tL[k][3][0] * clL[p + 0] + tL[k][3][1] * clL[p + 1] + tL[k][3][2] * clL[p + 2] +tL[k][3][3] * clL[p + 3];
					sumR += tR[k][3][0] * clR[p + 0] + tR[k][3][1] * clR[p + 1] + tR[k][3][2] * clR[p + 2] +tR[k][3][3] * clR[p + 3];
					clP[p + 3] = sumL * sumR;
                                        #ifdef _ASM_DEBUG
                                        __asm__ ( "int $0x3" );
                                        #endif
#endif
/*
#					endif
*/
// parallelisation
                                        p += 4;
					
				}
			}
			p->setIsClDirty(false);
		}
	}
		
	Node *r = t->getRoot();
	MbVector<double> f = getActiveBasefreq()->getFreq();
	//double *clP = clPtr[r->getActiveCl()][r->getIdx()];
	clP = clPtr[r->getActiveCl()][r->getIdx()];
	double catProb = 1.0 / numGammaCats;
	double lnL = 0.0;
        #pragma omp parallel for reduction ( + : lnL )
	for (int c=0; c<numPatterns; c++){
                int p = c * 16;
		double siteProb = 0.0;
                #ifdef _TOM_AVX
                ALIGNED ( double siteProbTmp[4] );
                #endif

#ifdef _TOM_SSE3
                                        #ifdef _ASM_DEBUG
                                        __asm__ ( "int $0x3" );
                                        #endif
                __m128d
                        p1,
                        p2,
                        v1,
                        v2,
                        v3,
                        v4,
                        m1,
                        m2;

                m1 = _mm_set_pd ( f[1], f[0] );
                m2 = _mm_set_pd ( f[3], f[2] );


                p1 = _mm_mul_pd ( _mm_load_pd ( clP + p + 0 ), m1 );           // clP[p+0] * f[0], clP[p+1] * f[1]
                p2 = _mm_mul_pd ( _mm_load_pd ( clP + p + 2 ), m2 );           // clP[p+2] * f[2], clP[p+3] * f[3]
                v1 = _mm_hadd_pd ( p1, p2 );                                                      // clP[p+0] * f[0] + clP[p+2] * f[2], clP[p+1] * f[1] + clP[p + 3] * f[3]

                p += 4;

                p1 = _mm_mul_pd ( _mm_load_pd ( clP + p + 0 ), m1 );
                p2 = _mm_mul_pd ( _mm_load_pd ( clP + p + 2 ), m2 );           // clP[p+2] * f[2], clP[p+3] * f[3]
                v2 = _mm_hadd_pd ( p1, p2 );                                                      // clP[p+0] * f[0] + clP[p+2] * f[2], clP[p+1] * f[1] + clP[p + 3] * f[3]

                p += 4;

                p1 = _mm_mul_pd ( _mm_load_pd ( clP + p + 0 ), m1 );
                p2 = _mm_mul_pd ( _mm_load_pd ( clP + p + 2 ), m2 );           // clP[p+2] * f[2], clP[p+3] * f[3]
                v3 = _mm_hadd_pd ( p1, p2 );                                                      // clP[p+0] * f[0] + clP[p+2] * f[2], clP[p+1] * f[1] + clP[p + 3] * f[3]

                p += 4;

                p1 = _mm_mul_pd ( _mm_load_pd ( clP + p + 0 ), m1 );
                p2 = _mm_mul_pd ( _mm_load_pd ( clP + p + 2 ), m2 );           // clP[p+2] * f[2], clP[p+3] * f[3]
                v4 = _mm_hadd_pd ( p1, p2 );                                                      // clP[p+0] * f[0] + clP[p+2] * f[2], clP[p+1] * f[1] + clP[p + 3] * f[3]

                p += 4;

                p1 = _mm_hadd_pd ( v1, v2 );
                p2 = _mm_hadd_pd ( v3, v4 );

                v1 = _mm_hadd_pd ( p1, p2 );
                _mm_storel_pd ( &siteProb, _mm_hadd_pd ( v1, v1 ) );
                                        #ifdef _ASM_DEBUG
                                        __asm__ ( "int $0x3" );
                                        #endif

#elif _TOM_AVX
                                        #ifdef _ASM_DEBUG
                                        __asm__ ( "int $0x3" );
                                        #endif
                __m256d
                        m,
                        p1, p2, p3, p4;

                m = _mm256_set_pd ( f[3], f[2], f[1], f[0] );
                
                p1 = _mm256_mul_pd ( _mm256_load_pd ( clP + p ), m );
                p2 = _mm256_mul_pd ( _mm256_load_pd ( clP + p + 4 ), m );
                p3 = _mm256_mul_pd ( _mm256_load_pd ( clP + p + 8 ), m );
                p4 = _mm256_mul_pd ( _mm256_load_pd ( clP + p + 12 ), m );
                p += 16;

                p1 = _mm256_hadd_pd ( p1, p2 );
                p2 = _mm256_hadd_pd ( p3, p4 );
                
                p1 = _mm256_hadd_pd ( p1, p2 );

                p1 = _mm256_add_pd ( p1, _mm256_permute2f128_pd ( p1 , p1 , 1)  );

                p1 = _mm256_hadd_pd ( p1, p1 );

                //_mm_storel_pd ( &siteProb, _mm256_extractf128_pd ( p1, 0 ) );
                _mm256_store_pd ( siteProbTmp, p1 );

                siteProb = siteProbTmp[0];
                                        #ifdef _ASM_DEBUG
                                        __asm__ ( "int $0x3" );
                                        #endif

#else
                                        #ifdef _ASM_DEBUG
                                        __asm__ ( "int $0x3" );
                                        #endif

		siteProb += clP[p + 0] * f[0] + clP[p + 1] * f[1] + clP[p + 2] * f[2] + clP[p + 3] * f[3] ;
                p += 4;
		//clP += 4;
		siteProb += clP[p + 0] * f[0] + clP[p + 1] * f[1] + clP[p + 2] * f[2] + clP[p + 3] * f[3] ;
                p += 4;
		//clP += 4;
		siteProb += clP[p + 0] * f[0] + clP[p + 1] * f[1] + clP[p + 2] * f[2] + clP[p + 3] * f[3] ;
		//clP += 4;
		p += 4;
                siteProb += clP[p + 0] * f[0] + clP[p + 1] * f[1] + clP[p + 2] * f[2] + clP[p + 3] * f[3] ;
		//clP += 4;
                p += 4;
                                        #ifdef _ASM_DEBUG
                                        __asm__ ( "int $0x3" );
                                        #endif
#endif
//#		endif
		siteProb *= catProb;
		lnL += alignmentPtr->getNumSitesOfPattern(c) * log(siteProb);
	}

	delete [] tL;
	delete [] tR;
	myCurLnL = lnL;
	return lnL;
}
