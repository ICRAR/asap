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
#include "SDDefs.h"

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



SDStokesEngine::SDStokesEngine (const String& outputColumnName,
			   const String& inputColumnName)
: BaseMappedArrayEngine<Float,Float> (outputColumnName, inputColumnName)
{}


SDStokesEngine::SDStokesEngine (const Record& spec)
: BaseMappedArrayEngine<Float,Float> ()
{
    if (spec.isDefined("OUTPUTNAME")  &&  spec.isDefined("INPUTNAME")) {
        setNames (spec.asString("OUTPUTNAME"), spec.asString("INPUTNAME"));
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
    spec.define ("OUTPUTNAME", sourceName());    // Ger uses opposite meaning for source/target
    spec.define ("INPUTNAME", targetName());
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


void SDStokesEngine::getArray (uInt rownr, Array<Float>& output)
{
    IPosition inputShape = findInputShape (output.shape());
    Array<Float> input(inputShape);  
    roColumn().get(rownr, input);
//
    computeOnGet (output, input);
}

void SDStokesEngine::putArray (uInt rownr, const Array<Float>& input)
{
    throw(AipsError("This Virtual Column is not writable"));
}



void SDStokesEngine::setShape (uInt rownr, const IPosition& outputShape)   
{   
    BaseMappedArrayEngine<Float,Float>::setShape (rownr, findInputShape(outputShape));
}
  
void SDStokesEngine::setShapeColumn (const IPosition& outputShape)
{
    BaseMappedArrayEngine<Float,Float>::setShapeColumn (findInputShape(outputShape));
}
      
    
IPosition SDStokesEngine::shape (uInt rownr)
{
   IPosition inputShape = roColumn().shape (rownr);
   return findOutputShape(inputShape);
}



void SDStokesEngine::computeOnGet(Array<Float>& output,
                   		 const Array<Float>& input)
//
// array of shape (nBeam,nIF,nPol,nChan)
//
{

// Checks

   const uInt nDim = input.ndim();
   DebugAssert(nDim==4,AipsError);
   DebugAssert(array.ndim()==4,AipsError);
   const IPosition inputShape = input.shape();
   const uInt polAxis = asap::PolAxis;
   const uInt nPol = inputShape(polAxis);
   DebugAssert(nPol==1 || nPol==2 || nPol==3, AipsError);

// The silly Array slice operator does not give me back
// a const reference so have to caste it away

   Array<Float>& input2 = const_cast<Array<Float>&>(input);

// Slice coordnates

   IPosition start(nDim,0);
   IPosition end(input.shape()-1);

// Generate Slices

   start(polAxis) = 0;
   end(polAxis) = 0;
   Array<Float> C1 = input2(start,end);          // Input : C1
//
   start(polAxis) = 0;
   end(polAxis) = 0;
   Array<Float> I = output(start,end);           // Output : I
//
   if (nPol==1) {
      I = C1;
      return;
   }
//
   start(polAxis) = 1;
   end(polAxis) = 1;
   Array<Float> C2 = input2(start,end);          // Input : C1
//
   I = Float(0.5)*(C1 + C2);
   if (nPol <= 2) return;
//
   start(polAxis) = 2;
   end(polAxis) = 2;
   Array<Float> C3 = input2(start,end);          // Input : C3
//
   start(polAxis) = 3;
   end(polAxis) = 3;
   Array<Float> C4 = input2(start,end);          // Input : C4
//
   start(polAxis) = 1;
   end(polAxis) = 1;
   Array<Float> Q = output(start,end);           // Output : Q
   Q = Float(0.5)*(C1 - C2);
//
   start(polAxis) = 2;
   end(polAxis) = 2;
   Array<Float> U = output(start,end);           // Output : U
   U = C3;
//
   start(polAxis) = 3;
   end(polAxis) = 3;
   Array<Float> V = output(start,end);           // Output : V
   V = C4;
}



IPosition SDStokesEngine::findInputShape (const IPosition& outputShape) const
//
// Don't know how to handle the degeneracy that both
// XX    -> I
// XX,YY -> I
// 
{
   uInt axis = asap::PolAxis;
   uInt nPol = outputShape(axis);
   IPosition inputShape = outputShape;
   if (nPol==1) {
      inputShape(axis) = 2;            // XX YY -> I
   } else if (nPol==4) {
      inputShape(axis) = 4;            // XX YY R(XY) I(XY) -> I Q U V
   }
   return inputShape;
}

IPosition SDStokesEngine::findOutputShape (const IPosition& inputShape) const
{
   uInt axis = 2;
   uInt nPol = inputShape(axis);
   IPosition outputShape = inputShape;
   if (nPol==1) {
      outputShape(axis) = 1;            // XX -> I
   } else if (nPol==2) {
      outputShape(axis) = 1;            // XX YY -> I
   } else if (nPol==4) {
      outputShape(axis) = 4;            // XX YY R(XY) I(XY) -> I Q U V
   }
   return outputShape;
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


 
