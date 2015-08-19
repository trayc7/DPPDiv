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

#ifndef AlIGNMENT_H
#define AlIGNMENT_H

#include <string>
#include <vector>
#include <iostream>


class Alignment {

	public:
									Alignment(std::string fn);
									~Alignment(void);
		void						compress(void);
		int							getNumTaxa(void) { return numTaxa; }
		int							getNumChar(void);
		std::string					getTaxonName(int i) { return taxonNames[i]; }
		int							getIndexForTaxonNamed(std::string nm);
		int							getNucleotide(int i, int j);
		int							getNumSitesOfPattern(int i);
		void						getPossibleNucs (int nucCode, int *nuc);
		bool						isTaxonPresent(std::string nm);
		void						print(std::ostream &) const;
		int							getNumPatterns(void) { return numPatterns; }

	private:
		int							nucID(char nuc);
		int							numTaxa;
		int							numChar;
		int							**matrix;
		int							**compressedMatrix;
		int							*patternCount;
		int							numPatterns;
		bool						isCompressed;
		std::vector<std::string>	taxonNames;
};

// Class for holding partitioned data
// For general MK model for multi-state morph data need partition for each stateset size
class PartitionedMatrix {
	public:
									PartitionedMatrix(unsigned numTaxa,);
									~PartitionedMatrix();
	  unsigned			getNParts() (return nPartitions);
		unsigned 			getNumCharsVec() {return numCharsPerPart};
		void  				fillPartitions(unsigned index, const unsigned **mat);

	private:
		unsigned 									nPartitions;
		std::vector<unsigned**>		matrices;
		std::vector<unsigned**>		compressedMatrices;
		unsigned									*numCharsPerPart;
}

// Class to hold partitioned data and related variables for morphological data
// Analogous to Alignment class for nucleotide data
class MorphData: public Alignment {

	public:
									MorphData();
									~MorphData(void);
		unsigned 			getNumPartitions() {return numPartitions);

	private:
		unsigned 									numPartitions;
		PartitionedMatrix					partMat;
		std::vector<std::string>  rawData;

#endif
