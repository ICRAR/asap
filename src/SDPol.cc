//#---------------------------------------------------------------------------
//# SDPol.cc: Polarimetric functionality
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

#include <casa/Arrays/Array.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Containers/Record.h>
#include <casa/BasicSL/Constants.h>
#include <casa/BasicSL/String.h>
#include <casa/Utilities/Assert.h>

#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ColumnDesc.h>
#include <tables/Tables/TableRecord.h>
#include <tables/Tables/DataManError.h>


using namespace casa;
using namespace asap;



SDStokesEngine::SDStokesEngine (const String& sourceColumnName,
			   const String& targetColumnName)

: BaseMappedArrayEngine<Float,Float> (sourceColumnName, targetColumnName)
{}


SDStokesEngine::SDStokesEngine (const Record& spec)
: BaseMappedArrayEngine<Float,Float> ()
{
    if (spec.isDefined("SOURCENAME")  &&  spec.isDefined("TARGETNAME")) {
        setNames (spec.asString("SOURCENAME"), spec.asString("TARGETNAME"));
    }
}

SDStokesEngine::SDStokesEngine (const SDStokesEngine& that)
: BaseMappedArrayEngine<Float,Float> (that)
{}

SDStokesEngine::~SDStokesEngine()
{}


DataManager* SDStokesEngine::clone() const
{
    DataManager* dmPtr = new SDStokesEngine (*this);
    return dmPtr;
}


String SDStokesEngine::dataManagerType() const
{
    return className();
}

String SDStokesEngine::className()
{
    return "SDStokesEngine";
}

String SDStokesEngine::dataManagerName() const
{
    return sourceName();
}

Record SDStokesEngine::dataManagerSpec() const
{
    Record spec;
    spec.define ("SOURCENAME", sourceName());
    spec.define ("TARGETNAME", targetName());
    return spec;
}

DataManager* SDStokesEngine::makeObject (const String&,
						 const Record& spec)
{
    DataManager* dmPtr = new SDStokesEngine(spec);
    return dmPtr;
}


void SDStokesEngine::registerClass()
{
    DataManager::registerCtor (className(), makeObject);
}


void SDStokesEngine::create (uInt initialNrrow)
{
    BaseMappedArrayEngine<Float,Float>::create (initialNrrow);
}

void SDStokesEngine::prepare()
{
    BaseMappedArrayEngine<Float,Float>::prepare();
}

Bool SDStokesEngine::canAccessArrayColumnCells (Bool& reask) const
{
    reask = False;
    return True;
}


void SDStokesEngine::getArray (uInt rownr, Array<Float>& array)
{
    Array<Float> target(array.shape());
    roColumn().get(rownr, target);
//
    computeOnGet (array, target);
}

void SDStokesEngine::putArray (uInt rownr, const Array<Float>& array)
{
    throw(AipsError("This Virtual Column is not writable"));
}



void SDStokesEngine::computeOnGet(Array<Float>& array,
                   		 const Array<Float>& target)
//
// array of shape (nBeam,nIF,nPol,nChan)
//
{
   DebugAssert(target.ndim()==4,AipsError);
   DebugAssert(array.ndim()==4,AipsError);

// The silly Array slice operator does not give me back
// a const reference so have to caste it away

   Array<Float>& target2 = const_cast<Array<Float>&>(target);

// uInt polAxis = asap::PolAxis;
   uInt polAxis = 2;
//
   const uInt nDim = target.ndim();
   IPosition start(nDim,0);
   IPosition end(target.shape()-1);

// Input slices

   start(polAxis) = 0;
   end(polAxis) = 0;
   Array<Float> C1 = target2(start,end);
//
   start(polAxis) = 1;
   end(polAxis) = 1;
   Array<Float> C2 = target2(start,end);
//
   start(polAxis) = 2;
   end(polAxis) = 2;
   Array<Float> C3 = target2(start,end);
//
   start(polAxis) = 3;
   end(polAxis) = 3;
   Array<Float> C4 = target2(start,end);

// Compute Output slices

   start(polAxis) = 0;
   end(polAxis) = 0;
   Array<Float> I = array(start,end); 
   I = Float(0.5)*(C1 + C2);
//
   start(polAxis) = 1;
   end(polAxis) = 1;
   Array<Float> Q = array(start,end); 
   Q = Float(0.5)*(C1 - C2);
//
   start(polAxis) = 2;
   end(polAxis) = 2;
   Array<Float> U = array(start,end); 
   U = C3;
//
   start(polAxis) = 3;
   end(polAxis) = 3;
   Array<Float> V = array(start,end); 
   V = C4;
}




// SDPolUtil

Array<Float> SDPolUtil::polarizedIntensity (const Array<Float>& Q,
                                            const Array<Float>& U)
{
   Array<Float> t1 = pow(Q,Double(2.0));
   Array<Float> t2 = pow(U,Double(2.0));
   return sqrt(t1+t2);
}


Array<Float> SDPolUtil::positionAngle (const Array<Float>& Q,
                                       const Array<Float>& U)
{
   return Float(180.0/C::pi/2.0)*atan2(Q,U);       // Degrees
}


void SDPolUtil::rotateXYPhase (Array<Float>& C3,
                               Array<Float>& C4,
                               Float phase)
{
   Float cosVal = cos(C::pi/180.0*phase);
   Float sinVal = sin(C::pi/180.0*phase);
//
   C3 = C3*cosVal - C4*sinVal;
   C4 = C3*sinVal + C4*cosVal;
}


 
