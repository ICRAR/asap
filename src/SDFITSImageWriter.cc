//#---------------------------------------------------------------------------
//# SDFITSImageWriter.cc: ASAP class to write out single dish spectra as FITS images
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

#include "SDFITSImageWriter.h"

#include <casa/aips.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/VectorIter.h>
#include <casa/Utilities/CountedPtr.h>

#include <coordinates/Coordinates/CoordinateUtil.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <coordinates/Coordinates/SpectralCoordinate.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>
#include <coordinates/Coordinates/StokesCoordinate.h>
#include <coordinates/Coordinates/ObsInfo.h>

#include <images/Images/ImageFITSConverter.h>
#include <images/Images/TempImage.h>

#include <lattices/Lattices/ArrayLattice.h>

#include <measures/Measures/MEpoch.h>

#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>

#include "Definitions.h"
#include "SDContainer.h"
#include "SDMemTable.h"

using namespace casa;
using namespace asap;


SDFITSImageWriter::SDFITSImageWriter()
{;}

SDFITSImageWriter::~SDFITSImageWriter()
{;}


Bool SDFITSImageWriter::write(const SDMemTable& sdTable, 
                              const String& rootName, Bool verbose)
{

// Get global Header from Table

   SDHeader header = sdTable.getSDHeader();
   MEpoch::Ref timeRef(MEpoch::UTC);              // Should be in header   
   MDirection::Types dirRef(MDirection::J2000);   // Should be in header   

// Prepare initial ObsInfo

   ObsInfo oi;
   oi.setTelescope(String(header.antennaname));
   oi.setObserver(String(header.observer));

// Column keywords

   Table tab = sdTable.table();
   ROArrayColumn<Double> dir(tab, String("DIRECTION"));
   ROScalarColumn<Double> time(tab, "TIME");
   ROArrayColumn<uInt> freqid(tab, "FREQID");
   ROScalarColumn<String> src(tab, "SRCNAME");

// Output Image Shape; spectral axis to be updated

   IPosition shapeOut(4,1);                 // spectral, direction (2), stokes

// Axes (should be in header)

   const uInt beamAxis = asap::BeamAxis;
   const uInt ifAxis = asap::IFAxis;
   const uInt stokesAxis = asap::PolAxis;
   const uInt chanAxis = asap::ChanAxis;
   const Unit RAD(String("rad"));

// Temps

   Vector<Int> whichStokes(1,1);
   Array<Double> whichDir;
   Vector<Double> lonLat(2);
   IPosition posDir(2,0);
   Projection proj(Projection::SIN);                   // What should we use ?
   Matrix<Double> xForm(2,2);
   xForm = 0.0; xForm.diagonal() = 1.0;
   Vector<Double> incLonLat(2,0.0);
   Vector<Double> refPixLonLat(2,0.0);

// Loop over rows

   uInt maxMem = 128;
   Bool preferVelocity = False;
   Bool opticalVelocity = False;
   Int bitPix = -32;                     // FLoating point
   Float minPix = 1.0;
   Float maxPix = -1.0;
   Bool overWrite = True;
   Bool degLast = False;
   Bool reallyVerbose = False;
//
   String errMsg;
   const uInt nRows = sdTable.nRow();
   for (uInt iRow=0; iRow<nRows; iRow++) {

// Get data

      const MaskedArray<Float>& dataIn(sdTable.rowAsMaskedArray(iRow));
      const Array<Float>& values = dataIn.getArray();
      const Array<Bool>& mask = dataIn.getMask();
      const IPosition& shapeIn = dataIn.shape();
      shapeOut(0) = shapeIn(chanAxis);
      TiledShape tShapeOut(shapeOut);

// Update ObsInfo (time changes per integration)

      Double dTmp;
      time.get(iRow, dTmp);
      MVEpoch tmp2(Quantum<Double>(dTmp, Unit(String("d"))));
      MEpoch epoch(tmp2, timeRef);
      oi.setObsDate(epoch);

// Iterate through data in this row by spectra

      ReadOnlyVectorIterator<Float> itData(values, chanAxis);
      ReadOnlyVectorIterator<Bool> itMask(mask, chanAxis);
      while (!itData.pastEnd()) {
         const IPosition& pos = itData.pos();

// Form SpectralCoordinate 

         Vector<uInt> iTmp;
         freqid.get(iRow, iTmp);
         SpectralCoordinate sC = sdTable.getCoordinate(iTmp(pos(ifAxis)));

// Form DirectionCoordinate
 
         dir.get(iRow, whichDir);
         posDir(0) = pos(beamAxis);
         posDir(1) = 0;
         lonLat[0] = whichDir(posDir);
//
         posDir(0) = pos(beamAxis);
         posDir(1) = 1;
         lonLat[1] = whichDir(posDir);
         DirectionCoordinate dC(dirRef, proj, lonLat[0], lonLat[1],
                               incLonLat[0], incLonLat[1], xForm, 
                               refPixLonLat[0], refPixLonLat[1]);

// Form Stokes Coordinate (no true Stokes info yet);

         whichStokes(0) = convertStokes(pos(stokesAxis));
         StokesCoordinate stC(whichStokes);

// Create CoordinateSystem

         CoordinateSystem cSys;
         cSys.addCoordinate(sC);
         cSys.addCoordinate(dC);
         cSys.addCoordinate(stC);
         cSys.setObsInfo(oi);

// Reform data into correct shape for output

         Array<Float> t1(itData.array().reform(shapeOut));
         Array<Bool> m1(itMask.array().reform(shapeOut));

// Create aips++ Image

     
         TempImage<Float> tIm(tShapeOut, cSys);
         tIm.put(t1);
//
         ArrayLattice<Bool> latMask(m1);
         tIm.attachMask(latMask);

// Write out as FITS Image file

         ostringstream oss;
         oss << "row" << iRow << "_beam" << pos(0) << "_if" 
	     << pos(1) << "_pol" << pos(2) << "_" << src(iRow) << ".fits";
         String fileName;
         if (rootName.length()>0) {
            fileName = rootName + "_" + String(oss);
         } else {
            fileName = String(oss);
         }
         if (verbose) cerr << "Writing row " << iRow 
			   << " into file " << fileName << endl;
//
         Bool ok = ImageFITSConverter::ImageToFITS(errMsg, tIm, fileName, 
						   maxMem, preferVelocity,
						   opticalVelocity, bitPix, 
						   minPix, maxPix, overWrite,
						   degLast, reallyVerbose);
         if (!ok) {
            cerr << "Error writing fits - " << errMsg << endl;
         }

// Next spectrum

         itData.next();
         itMask.next();
      }
   }
//
   if (!verbose) {
      cerr << "Wrote " << nRows << " into individual FITS files" << endl;
   }
//   
   return True;
}


Int SDFITSImageWriter::convertStokes(Int val)
{
   Stokes::StokesTypes stokes = Stokes::RR;
   if (val==0) {
      stokes = Stokes::RR;
   } else if (val==1) {
      stokes = Stokes::LL;
   } else if (val==2) {
      stokes = Stokes::RL;
   } else if (val==3) {
      stokes = Stokes::LR;
   } else {
      stokes = Stokes::Undefined;
   }
//
   return Int(stokes);
}
