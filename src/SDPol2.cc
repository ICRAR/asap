//#---------------------------------------------------------------------------
//# SDPol2.cc: Templated polarimetric functionality
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

#include "SDPol.h"
#include "SDDefs.h"

#include <casa/Arrays/Array.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/BasicSL/String.h>
#include <casa/Utilities/Assert.h>
#include <casa/Utilities/DataType.h>
#include <casa/Utilities/Assert.h>


using namespace casa;
using namespace asap;


template <class T>
Array<T> SDPolUtil::stokesData (Array<T>& rawData, Bool doLinear)
//
// Generate data for each Stokes parameter from the
// raw flags.  This is a lot of computational work and may
// not be worth the effort.
//
// Designed for use 
//   Bool for mask
//   Float for TSys
//
{
   T* t;
   DataType type = whatType(t);
   AlwaysAssert(type==TpFloat || type==TpBool, AipsError);
//
   IPosition shapeIn = rawData.shape();
   uInt nPol = shapeIn(asap::PolAxis);
   const uInt nDim = shapeIn.nelements();
   Array<T> stokesData;
//
   IPosition start(nDim,0);
   IPosition end(shapeIn-1);
   IPosition shapeOut = shapeIn;
//
   if (doLinear) {
      if (nPol==1) {
         stokesData.resize(shapeOut);
         stokesData = rawData;
      } else if (nPol==2 || nPol==4) {

// Set shape of output array

         if (nPol==2) {
            shapeOut(asap::PolAxis) = 1;
         } else {
            shapeOut(asap::PolAxis) = 4;
         }
         stokesData.resize(shapeOut);

// Get reference slices and assign/compute

         start(asap::PolAxis) = 0;
         end(asap::PolAxis) = 0;
         Array<T> M1In = rawData(start,end);
//
         start(asap::PolAxis) = 1;
         end(asap::PolAxis) = 1;
         Array<T> M2In = rawData(start,end);
//
         start(asap::PolAxis) = 0;
         end(asap::PolAxis) = 0;
         Array<T> M1Out = stokesData(start,end);         
         M1Out = SDPolUtil::andArrays (M1In, M2In);         // I
//
         if (nPol==4) {   
            start(asap::PolAxis) = 2;
            end(asap::PolAxis) = 2;
            Array<T> M3In = rawData(start,end);
//
            start(asap::PolAxis) = 3;
            end(asap::PolAxis) = 3;
            Array<T> M4In = rawData(start,end);
//
            start(asap::PolAxis) = 1;
            end(asap::PolAxis) = 1;
            Array<T> M2Out = stokesData(start,end);
            M2Out = M1Out;                                  // Q
//
            start(asap::PolAxis) = 2;
            end(asap::PolAxis) = 2;
            Array<T> M3Out = stokesData(start,end);
            M3Out = M3In;                                   // U
//
            start(asap::PolAxis) = 3;
            end(asap::PolAxis) = 3;
            Array<T> M4Out = stokesData(start,end);
            M4Out = M4In;                                   // V
         }
      } else {
         throw(AipsError("Can only handle 1,2 or 4 polarizations"));
      }
   } else {
      throw (AipsError("Only implemented for Linear polarizations"));
   }
//
   return stokesData;
}

