//#---------------------------------------------------------------------------
//# MathUtilities.cc: General math operations
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

#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/MaskedArray.h>
#include <casa/Arrays/MaskArrMath.h>
#include <casa/BasicSL/String.h>

#include "MathUtils.h"

using namespace casa;

float mathutil::statistics(const String& which,  
			   const MaskedArray<Float>& data)
{
   String str(which);
   str.upcase();
   if (str.contains(String("MIN"))) {
      return min(data); 
   } else if (str.contains(String("MAX"))) {
      return max(data);
   } else if (str.contains(String("SUMSQ"))) {
      return sumsquares(data);
   } else if (str.contains(String("SUM"))) {
      return sum(data);
   } else if (str.contains(String("MEAN"))) {
      return mean(data);
   } else if (str.contains(String("VAR"))) {
      return variance(data); 
   } else if (str.contains(String("STDDEV"))) {
      return stddev(data);
   } else if (str.contains(String("AVDEV"))) {
      return avdev(data);
   } else if (str.contains(String("RMS"))) {
      uInt n = data.nelementsValid();
      return sqrt(sumsquares(data)/n);
   } else if (str.contains(String("MED"))) {
      return median(data);
   }
}
  

void mathutil::replaceMaskByZero(Vector<Float>& data, const Vector<Bool>& mask)
{
   for (uInt i=0; i<data.nelements(); i++) {
      if (!mask[i]) data[i] = 0.0;
   }
}

