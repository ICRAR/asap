//#---------------------------------------------------------------------------
//# SDMath.h: A collection of single dish mathematical operations
//#---------------------------------------------------------------------------
//# Copyright (C) 2004
//# ATNF
//#
//# This program is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by the Free
//# Software Foundation; either version 2 of the License, or (at your option)
//# any later version.
//#
//# This program is distributed in the hope that it will be useful, but
//# WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
//# Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with this program; if not, write to the Free Software Foundation, Inc.,
//# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning this software should be addressed as follows:
//#        Internet email: Malte.Marquarding@csiro.au
//#        Postal address: Malte Marquarding,
//#                        Australia Telescope National Facility,
//#                        P.O. Box 76,
//#                        Epping, NSW, 2121,
//#                        AUSTRALIA
//#
//# $Id:
//#---------------------------------------------------------------------------
#ifndef SDMATH_H
#define SDMATH_H

#include <string>
#include <vector>
#include <casa/aips.h>
#include <casa/Utilities/CountedPtr.h>

namespace asap {

class SDMemTable;

class SDMath {

 public:

// Default constructor 
   SDMath();

// Copy Constructor (copy semantics)
   SDMath (const SDMath& other);

// Assignment  (copy semantics)
   SDMath &operator=(const SDMath& other);

// Quotient
   casa::CountedPtr<SDMemTable> quotient(const casa::CountedPtr<SDMemTable>& on, 
					 const casa::CountedPtr<SDMemTable>& off);

// Average in time
   casa::CountedPtr<SDMemTable>  average(const casa::Block<casa::CountedPtr<SDMemTable> >& in,
                                         const casa::Vector<casa::Bool>& mask,
                                         casa::Bool scanAverage, const std::string& weightStr);

// Statistics
   std::vector<float> statistic(const casa::CountedPtr<SDMemTable>& in, 
    		                const std::vector<bool>& mask, const casa::String& which);

// Bin up spectra
   SDMemTable* bin(const SDMemTable& in, casa::Int width);

// Smooth
   SDMemTable* smooth (const SDMemTable& in, const casa::String& kernel,
                       casa::Float width, casa::Bool doAll);

// Simple mathematical operations.  what=0 (mul) or 1 (add)
   SDMemTable* simpleOperate(const SDMemTable& in, casa::Float offset, 
                             casa::Bool doAll, casa::uInt what);

// Average polarizations
   SDMemTable* averagePol(const SDMemTable& in, const casa::Vector<casa::Bool>& mask);

 private:

// Weighting type for time averaging

  enum WeightType {NONE,VAR,TSYS};

// Function to use accumulate data during time averaging

  void accumulate (casa::Double& timeSum, casa::Double& intSum, casa::Int& nAccum,
                   casa::MaskedArray<casa::Float>& sum, casa::Array<casa::Float>& sumSq,
                   casa::Array<casa::Float>& nPts, casa::Array<casa::Float>& tSysSum,
                   const casa::Array<casa::Float>& tSys,  const casa::Array<casa::Float>& nInc,
                   const casa::Vector<casa::Bool>& mask, casa::Double time, casa::Double interval,
                   const casa::Block<casa::CountedPtr<SDMemTable> >& in,
                   casa::uInt iTab, casa::uInt iRow, casa::uInt axis, casa::uInt nAxesSub,
                   casa::Bool useMask, WeightType wtType);

// Function to fill Scan Container when averaging in time

  void fillSDC (SDContainer& sc, const casa::Array<casa::Bool>& mask,
                const casa::Array<casa::Float>& data,
                const casa::Array<casa::Float>& tSys,
                casa::Int scanID, casa::Double timeStamp,
                casa::Double interval, const casa::String& sourceName,
                const casa::Vector<casa::uInt>& freqID);

// Put the data and mask into the SDContainer
   void putDataInSDC (SDContainer& sc, const casa::Array<casa::Float>& data,
                      const casa::Array<casa::Bool>& mask);

// Function to normalize data when averaging in time

  void normalize (casa::MaskedArray<casa::Float>& data,
                  const casa::Array<casa::Float>& sumSq,
                  const casa::Array<casa::Float>& nPts,
                  WeightType wtType, casa::Int axis, casa::Int nAxes);

// Function to get the current cursor location
   void getCursorLocation (casa::IPosition& start, casa::IPosition& end,
                           const SDMemTable& in);

// Convert weight string to enum value

   void convertWeightString (WeightType& wt, const std::string& weightStr);
};

} // namespace

#endif
