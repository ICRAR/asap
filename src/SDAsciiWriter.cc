//#---------------------------------------------------------------------------
//# SDAsciiWriter.cc: ASAP class to write out single dish spectra as FITS images
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
//# $Id$
//#---------------------------------------------------------------------------

#include <casa/aips.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/VectorIter.h>
#include <casa/Utilities/CountedPtr.h>
#include <casa/Quanta/Quantum.h>
#include <casa/Quanta/MVAngle.h>


#include <casa/iostream.h>
#include <casa/fstream.h>

#include <coordinates/Coordinates/CoordinateUtil.h>
#include <coordinates/Coordinates/SpectralCoordinate.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>
#include <coordinates/Coordinates/StokesCoordinate.h>

#include <measures/Measures/MEpoch.h>

#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>

#include "SDDefs.h"
#include "SDContainer.h"
#include "SDMemTable.h"
#include "SDAsciiWriter.h"

using namespace casa;
using namespace asap;


SDAsciiWriter::SDAsciiWriter()
{;}

SDAsciiWriter::~SDAsciiWriter()
{;}


Bool SDAsciiWriter::write(const SDMemTable& sdTable, const String& fileName, Bool toStokes)
{

// Get global Header from Table

   SDHeader header = sdTable.getSDHeader();
   MEpoch::Ref timeRef(MEpoch::UTC);              // Should be in header   
   MDirection::Types dirRef(MDirection::J2000);   // Should be in header   

// Column keywords

   Table tab = sdTable.table();
   ROArrayColumn<Double> dir(tab, String("DIRECTION"));
   ROScalarColumn<Double> time(tab, "TIME");
   ROArrayColumn<uInt> freqid(tab, "FREQID");
   ROScalarColumn<String> src(tab, "SRCNAME");

// Temps

   Vector<Int> whichStokes(1,1);
   Array<Double> whichDir;
   Vector<Double> lonLat(2);
   IPosition posDir(2,0);
   const Unit RAD(String("rad"));

// Open file

   String fName(fileName);
   if (fileName.length()==0) fName = String("ascii.txt");
   ofstream of(fName.chars(), ios::trunc);

// Write header

   of << "row beam IF pol source longitude latitude time nchan spectrum mask" 
      << endl;
   
// Loop over rows

   const uInt nRows = sdTable.nRow();
   for (uInt iRow=0; iRow<nRows; iRow++) {

// Get data

      const MaskedArray<Float>& dataIn(sdTable.rowAsMaskedArray(iRow,False,toStokes));
      const Array<Float>& values = dataIn.getArray();
      const Array<Bool>& mask = dataIn.getMask();

// Epoch

      Double dTmp;
      time.get(iRow, dTmp);
      MVEpoch tmp2(Quantum<Double>(dTmp, Unit(String("d"))));
      MEpoch epoch(tmp2, timeRef);

// Iterate through data in this row by spectra

      ReadOnlyVectorIterator<Float> itData(values, asap::ChanAxis);
      ReadOnlyVectorIterator<Bool> itMask(mask, asap::ChanAxis);
      while (!itData.pastEnd()) {
         const IPosition& pos = itData.pos();

// FreqID

         Vector<uInt> iTmp;
         freqid.get(iRow, iTmp);

// Direction
 
         dir.get(iRow, whichDir);
         posDir(0) = pos(asap::BeamAxis);
         posDir(1) = 0;
         lonLat[0] = whichDir(posDir);
//
         posDir(0) = pos(asap::BeamAxis);
         posDir(1) = 1;
         lonLat[1] = whichDir(posDir);

// Write preamble

         of << iRow << "  " << pos(asap::BeamAxis) << " " <<  pos(asap::IFAxis) << " " << 
               pos(asap::PolAxis) << " " <<
               src(iRow) <<  " " << formatDirection(lonLat) << " " << 
               dTmp << " " << itData.vector().nelements() << " ";
// Write data

         const Vector<Float>& data = itData.vector();
         const Vector<Bool>& mask = itMask.vector();
         const uInt n = data.nelements();
         for (uInt i=0; i<n; i++) {
            of << data[i] << " ";
         }
// Write mask

         for (uInt i=0; i<n; i++) {
            of << mask[i] << " ";
         }
         of << endl;

// Next spectrum

         itData.next();
         itMask.next();
      }
   }
//
   of.close();
   cout << "Wrote " << nRows << " rows into file " << fileName << endl;
//   
   return True;
}



String SDAsciiWriter::formatDirection(const Vector<Double>& lonLat)
{ 
   MVAngle x1(lonLat(0));
   String s1 = x1.string(MVAngle::TIME, 12);
// 
   MVAngle x2(lonLat(1));
   String s2 = x2.string(MVAngle::ANGLE, 12);
//
   String ss = s1 + String(" ") + s2;
   return ss;
}

