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
  //public:
  casa::CountedPtr<SDMemTable> average(const casa::CountedPtr<SDMemTable>& in);
  casa::CountedPtr<SDMemTable> quotient(const casa::CountedPtr<SDMemTable>& on, 
					 const casa::CountedPtr<SDMemTable>& off);
  void multiplyInSitu(SDMemTable* in, casa::Float factor);

  casa::CountedPtr<SDMemTable> multiply(const casa::CountedPtr<SDMemTable>& in, 
				  casa::Float factor);
  casa::CountedPtr<SDMemTable> add(const casa::CountedPtr<SDMemTable>& in, 
			     casa::Float offset);
  
  casa::CountedPtr<SDMemTable> hanning(const casa::CountedPtr<SDMemTable>& in);

  casa::CountedPtr<SDMemTable> 
  averages(const casa::Block<casa::CountedPtr<SDMemTable> >& in,
	   const casa::Vector<casa::Bool>& mask);

  casa::CountedPtr<SDMemTable> 
  averagePol(const casa::CountedPtr<SDMemTable>& in, const casa::Vector<casa::Bool>& mask);

  std::vector<float> statistic(const casa::CountedPtr<SDMemTable>& in, 
   		                const std::vector<bool>& mask, const std::string& which);
  
  casa::CountedPtr<SDMemTable> bin(const casa::CountedPtr<SDMemTable>& in, 
			     casa::Int width);

// private (not actually...)

  float theStatistic(const std::string& which,  const casa::MaskedArray<casa::Float>& data);
  
};

} // namespace

#endif
