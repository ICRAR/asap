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

namespace SDMath {

// Quotient

  casa::CountedPtr<SDMemTable> quotient(const casa::CountedPtr<SDMemTable>& on, 
					 const casa::CountedPtr<SDMemTable>& off);

// Multiply

  void multiplyInSitu(SDMemTable* in, casa::Float factor, casa::Bool all);
  casa::CountedPtr<SDMemTable> multiply(const casa::CountedPtr<SDMemTable>& in, 
                                        casa::Float factor, casa::Bool all);

// Addition

  void addInSitu (SDMemTable* in, casa::Float offset, casa::Bool all);
  casa::CountedPtr<SDMemTable> add(const casa::CountedPtr<SDMemTable>& in, 
                                   casa::Float offset, casa::Bool all);

//  Hanning

  casa::CountedPtr<SDMemTable> hanning(const casa::CountedPtr<SDMemTable>& in);

// Bin up

  casa::CountedPtr<SDMemTable> bin(const casa::CountedPtr<SDMemTable>& in, 
			     casa::Int width);

// Average in time

  casa::CountedPtr<SDMemTable>
  average (const casa::Block<casa::CountedPtr<SDMemTable> >& in,
           const casa::Vector<casa::Bool>& mask,
           bool scanAverage, const std::string& weightStr);

// Average polarizations

  casa::CountedPtr<SDMemTable> 
  averagePol(const casa::CountedPtr<SDMemTable>& in, const casa::Vector<casa::Bool>& mask);

// Statistics

  std::vector<float> statistic(const casa::CountedPtr<SDMemTable>& in, 
   		                const std::vector<bool>& mask, const std::string& which);

// private (not actually...)

// Weighting type for time averaging

  enum weightType {NONE,VAR,TSYS};

// Function to use accumulate data during time averaging

  void accumulate (casa::Double& timeSum, casa::Double& intSum, casa::Int& nAccum,
                   casa::MaskedArray<casa::Float>& sum, casa::Array<casa::Float>& sumSq,
                   casa::Array<casa::Float>& nPts, casa::Array<casa::Float>& tSysSum,
                   const casa::Array<casa::Float>& tSys,  const casa::Array<casa::Float>& nInc,
                   const casa::Vector<casa::Bool>& mask, casa::Double time, casa::Double interval,
                   const casa::Block<casa::CountedPtr<SDMemTable> >& in,
                   casa::uInt iTab, casa::uInt iRow, casa::uInt axis, casa::uInt nAxesSub,
                   casa::Bool useMask, weightType wtType);

// Function to fill Scan Container when averaging in time

  void fillSDC (SDContainer& sc, const casa::Array<casa::Bool>& mask,
                const casa::Array<casa::Float>& data,
                const casa::Array<casa::Float>& tSys,
                casa::Int scanID, casa::Double timeStamp,
                casa::Double interval, const casa::String& sourceName,
                const casa::Vector<casa::uInt>& freqID);

// Function to normalize data when averaging in time

  void normalize (casa::MaskedArray<casa::Float>& data,
                  const casa::Array<casa::Float>& sumSq,
                  const casa::Array<casa::Float>& nPts,
                  weightType wtType, casa::Int axis, casa::Int nAxes);

// Functions for simple mathematical operations.  what=0 (mul) or 1 (add)

  SDMemTable* localOperate (const SDMemTable& in, casa::Float offset, 
                            casa::Bool doAll, casa::uInt what);

// Function to get the current cursor location
   void getCursorLocation (casa::IPosition& start, casa::IPosition& end,
                           const SDMemTable& in);
};

} // namespace

#endif





