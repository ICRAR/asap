//#---------------------------------------------------------------------------
//# SDMath.cc: A collection of single dish mathematical operations
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
#include <vector>

#include <casa/aips.h>
#include <casa/BasicSL/String.h>
#include <casa/Arrays/IPosition.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/ArrayIter.h>
#include <casa/Arrays/VectorIter.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Arrays/MaskedArray.h>
#include <casa/Arrays/MaskArrMath.h>
#include <casa/Arrays/MaskArrLogi.h>
#include <casa/BasicMath/Math.h>
#include <casa/Containers/Block.h>
#include <casa/Quanta/QC.h>
#include <casa/Utilities/Assert.h>
#include <casa/Exceptions.h>

#include <scimath/Mathematics/VectorKernel.h>
#include <scimath/Mathematics/Convolver.h>
#include <scimath/Mathematics/InterpolateArray1D.h>
#include <scimath/Functionals/Polynomial.h>

#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ReadAsciiTable.h>

#include <lattices/Lattices/LatticeUtilities.h>
#include <lattices/Lattices/RebinLattice.h>
#include <coordinates/Coordinates/SpectralCoordinate.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <coordinates/Coordinates/CoordinateUtil.h>
#include <coordinates/Coordinates/VelocityAligner.h>

#include "MathUtils.h"
#include "SDDefs.h"
#include "SDContainer.h"
#include "SDMemTable.h"

#include "SDMath.h"

using namespace casa;
using namespace asap;


SDMath::SDMath()
{;}

SDMath::SDMath(const SDMath& other)
{

// No state

}

SDMath& SDMath::operator=(const SDMath& other) 
{
  if (this != &other) {
// No state
  }
  return *this;
}

SDMath::~SDMath()
{;}


CountedPtr<SDMemTable> SDMath::average(const Block<CountedPtr<SDMemTable> >& in,
				       const Vector<Bool>& mask, Bool scanAv,
				       const String& weightStr) const
//Bool alignVelocity)
//
// Weighted averaging of spectra from one or more Tables.
//
{
   Bool alignVelocity = False;

// Convert weight type
  
  WeightType wtType = NONE;
  convertWeightString(wtType, weightStr);

// Create output Table by cloning from the first table

  SDMemTable* pTabOut = new SDMemTable(*in[0],True);

// Setup

  IPosition shp = in[0]->rowAsMaskedArray(0).shape();      // Must not change
  Array<Float> arr(shp);
  Array<Bool> barr(shp);
  const Bool useMask = (mask.nelements() == shp(asap::ChanAxis));

// Columns from Tables

  ROArrayColumn<Float> tSysCol;
  ROScalarColumn<Double> mjdCol;
  ROScalarColumn<String> srcNameCol;
  ROScalarColumn<Double> intCol;
  ROArrayColumn<uInt> fqIDCol;

// Create accumulation MaskedArray. We accumulate for each channel,if,pol,beam
// Note that the mask of the accumulation array will ALWAYS remain ALL True.
// The MA is only used so that when data which is masked Bad is added to it,
// that data does not contribute.

  Array<Float> zero(shp);
  zero=0.0;
  Array<Bool> good(shp);
  good = True;
  MaskedArray<Float> sum(zero,good);

// Counter arrays

  Array<Float> nPts(shp);             // Number of points
  nPts = 0.0;
  Array<Float> nInc(shp);             // Increment
  nInc = 1.0;

// Create accumulation Array for variance. We accumulate for
// each if,pol,beam, but average over channel.  So we need
// a shape with one less axis dropping channels.

  const uInt nAxesSub = shp.nelements() - 1;
  IPosition shp2(nAxesSub);
  for (uInt i=0,j=0; i<(nAxesSub+1); i++) {
     if (i!=asap::ChanAxis) {
       shp2(j) = shp(i);
       j++;
     }
  }
  Array<Float> sumSq(shp2);
  sumSq = 0.0;
  IPosition pos2(nAxesSub,0);                        // For indexing

// Time-related accumulators

  Double time;
  Double timeSum = 0.0;
  Double intSum = 0.0;
  Double interval = 0.0;

// To get the right shape for the Tsys accumulator we need to
// access a column from the first table.  The shape of this
// array must not change

  Array<Float> tSysSum;
  {
    const Table& tabIn = in[0]->table();
    tSysCol.attach(tabIn,"TSYS");
    tSysSum.resize(tSysCol.shape(0));
  }
  tSysSum =0.0;
  Array<Float> tSys;

// Scan and row tracking

  Int oldScanID = 0;
  Int outScanID = 0;
  Int scanID = 0;
  Int rowStart = 0;
  Int nAccum = 0;
  Int tableStart = 0;

// Source and FreqID

  String sourceName, oldSourceName, sourceNameStart;
  Vector<uInt> freqID, freqIDStart, oldFreqID;

// Velocity Aligner. We need an aligner for each Direction and FreqID
// combination.  I don't think there is anyway to know how many
// directions there are.
// For now, assume all Tables have the same Frequency Table

/*
  {
     MEpoch::Ref timeRef(MEpoch::UTC);              // Should be in header
     MDirection::Types dirRef(MDirection::J2000);   // Should be in header
//
     SDHeader sh = in[0].getSDHeader();
     const uInt nChan = sh.nchan;
//
     const SDFrequencyTable freqTab = in[0]->getSDFreqTable();
     const uInt nFreqID = freqTab.length();
     PtrBlock<const VelocityAligner<Float>* > vA(nFreqID);

// Get first time from first table

     const Table& tabIn0 = in[0]->table();
     mjdCol.attach(tabIn0, "TIME");
     Double dTmp;
     mjdCol.get(0, dTmp);
     MVEpoch tmp2(Quantum<Double>(dTmp, Unit(String("d"))));
     MEpoch epoch(tmp2, timeRef);
//
     for (uInt freqID=0; freqID<nFreqID; freqID++) {
        SpectralCoordinate sC = in[0]->getCoordinate(freqID);
        vA[freqID] = new VelocityAligner<Float>(sC, nChan, epoch, const MDirection& dir, 
                                                const MPosition& pos, const String& velUnit, 
                                                MDoppler::Types velType, MFrequency::Types velFreqSystem)
     }
  }
*/

// Loop over tables

  Float fac = 1.0;
  const uInt nTables = in.nelements();
  for (uInt iTab=0; iTab<nTables; iTab++) {

// Should check that the frequency tables don't change if doing VelocityAlignment

// Attach columns to Table

     const Table& tabIn = in[iTab]->table();
     tSysCol.attach(tabIn, "TSYS");
     mjdCol.attach(tabIn, "TIME");
     srcNameCol.attach(tabIn, "SRCNAME");
     intCol.attach(tabIn, "INTERVAL");
     fqIDCol.attach(tabIn, "FREQID");

// Loop over rows in Table

     const uInt nRows = in[iTab]->nRow();
     for (uInt iRow=0; iRow<nRows; iRow++) {

// Check conformance

        IPosition shp2 = in[iTab]->rowAsMaskedArray(iRow).shape();
        if (!shp.isEqual(shp2)) {
           throw (AipsError("Shapes for all rows must be the same"));
        }

// If we are not doing scan averages, make checks for source and
// frequency setup and warn if averaging across them

// Get copy of Scan Container for this row

        SDContainer sc = in[iTab]->getSDContainer(iRow);
        scanID = sc.scanid;

// Get quantities from columns

        srcNameCol.getScalar(iRow, sourceName);
        mjdCol.get(iRow, time);
        tSysCol.get(iRow, tSys);
        intCol.get(iRow, interval);
        fqIDCol.get(iRow, freqID);

// Initialize first source and freqID

        if (iRow==0 && iTab==0) {
          sourceNameStart = sourceName;
          freqIDStart = freqID;
        }

// If we are doing scan averages, see if we are at the end of an
// accumulation period (scan).  We must check soutce names too,
// since we might have two tables with one scan each but different
// source names; we shouldn't average different sources together

        if (scanAv && ( (scanID != oldScanID)  ||
                        (iRow==0 && iTab>0 && sourceName!=oldSourceName))) {

// Normalize data in 'sum' accumulation array according to weighting scheme

           normalize(sum, sumSq, nPts, wtType, asap::ChanAxis, nAxesSub);

// Fill scan container. The source and freqID come from the
// first row of the first table that went into this average (
// should be the same for all rows in the scan average)

           Float nR(nAccum);
           fillSDC(sc, sum.getMask(), sum.getArray(), tSysSum/nR, outScanID,
                    timeSum/nR, intSum, sourceNameStart, freqIDStart);

// Write container out to Table

           pTabOut->putSDContainer(sc);

// Reset accumulators

           sum = 0.0;
           sumSq = 0.0;
           nAccum = 0;
//
           tSysSum =0.0;
           timeSum = 0.0;
           intSum = 0.0;
	   nPts = 0.0;

// Increment

           rowStart = iRow;              // First row for next accumulation
           tableStart = iTab;            // First table for next accumulation
           sourceNameStart = sourceName; // First source name for next accumulation
           freqIDStart = freqID;         // First FreqID for next accumulation
//
           oldScanID = scanID;
           outScanID += 1;               // Scan ID for next accumulation period
        }

// Accumulate

        accumulate(timeSum, intSum, nAccum, sum, sumSq, nPts, tSysSum, 
                    tSys, nInc, mask, time, interval, in, iTab, iRow, asap::ChanAxis, 
                    nAxesSub, useMask, wtType);
//
       oldSourceName = sourceName;
       oldFreqID = freqID;
     }
  }

// OK at this point we have accumulation data which is either
//   - accumulated from all tables into one row
// or
//   - accumulated from the last scan average
//
// Normalize data in 'sum' accumulation array according to weighting scheme
  normalize(sum, sumSq, nPts, wtType, asap::ChanAxis, nAxesSub);

// Create and fill container.  The container we clone will be from
// the last Table and the first row that went into the current
// accumulation.  It probably doesn't matter that much really...

  Float nR(nAccum);
  SDContainer sc = in[tableStart]->getSDContainer(rowStart);
  fillSDC(sc, sum.getMask(), sum.getArray(), tSysSum/nR, outScanID,
           timeSum/nR, intSum, sourceNameStart, freqIDStart);
  pTabOut->putSDContainer(sc);
//
  return CountedPtr<SDMemTable>(pTabOut);
}



CountedPtr<SDMemTable> SDMath::quotient(const CountedPtr<SDMemTable>& on,
                                        const CountedPtr<SDMemTable>& off,
                                        Bool preserveContinuum)  const
{
  const uInt nRowOn = on->nRow();
  const uInt nRowOff = off->nRow();
  Bool ok = (nRowOff==1&&nRowOn>0) ||
            (nRowOff>=1&&nRowOn==nRowOff);
  if (!ok) {
     throw (AipsError("The reference Scan Table can have one row or the same number of rows as the source Scan Table"));
  }

// Input Tables and columns

  Table tabOn = on->table();
  Table tabOff = off->table();
  ROArrayColumn<Float> tSysOn(tabOn, "TSYS");
  ROArrayColumn<Float> tSysOff(tabOff, "TSYS");

// Output Table cloned from input

  SDMemTable* pTabOut = new SDMemTable(*on, True);

// Loop over rows

  MaskedArray<Float>* pMOff = new MaskedArray<Float>(off->rowAsMaskedArray(0));
  IPosition shpOff = pMOff->shape();
//
  Array<Float> tSysOnArr, tSysOffArr;
  tSysOn.get(0, tSysOnArr);
  tSysOff.get(0, tSysOffArr);
//
  for (uInt i=0; i<nRowOn; i++) {
     MaskedArray<Float> mOn(on->rowAsMaskedArray(i));
     IPosition shpOn = mOn.shape(); 
     tSysOn.get(i, tSysOnArr);
//
     if (nRowOff>1) {
        delete pMOff;
        pMOff = new MaskedArray<Float>(off->rowAsMaskedArray(i));
        shpOff = pMOff->shape();
        tSysOff.get(i, tSysOffArr);
     }

// Conformance

     if (!shpOn.isEqual(shpOff)) {
        throw(AipsError("on/off data are not conformant"));
     }
     if (!tSysOnArr.shape().isEqual(tSysOffArr.shape())) {
        throw(AipsError("on/off Tsys data are not conformant"));
     }
     if (!shpOn.isEqual(tSysOnArr.shape())) {
        throw(AipsError("Correlation and Tsys data are not conformant"));
     }

// Get container

     SDContainer sc = on->getSDContainer(i);

// Compute and put quotient into container

     if (preserveContinuum) {     
        MaskedArray<Float> tmp = (tSysOffArr * mOn / *pMOff) - tSysOffArr;
        putDataInSDC(sc, tmp.getArray(), tmp.getMask());
     } else {
        MaskedArray<Float> tmp = (tSysOffArr * mOn / *pMOff) - tSysOnArr;
        putDataInSDC(sc, tmp.getArray(), tmp.getMask());
     }
     sc.putTsys(tSysOffArr);
     sc.scanid = i;

// Put new row in output Table
 
     pTabOut->putSDContainer(sc);
  }
  if (pMOff) delete pMOff;
//
  return CountedPtr<SDMemTable>(pTabOut);
}


CountedPtr<SDMemTable> SDMath::simpleBinaryOperate (const CountedPtr<SDMemTable>& left,
                                                    const CountedPtr<SDMemTable>& right,
                                                    const String& op)  const
//
// Simple binary Table operators. add, subtract, multiply, divide (what=0,1,2,3)
//
{

// CHeck operator

  String op2(op);
  op2.upcase();
  uInt what = 0;
  if (op2=="ADD") {
     what = 0;
  } else if (op2=="SUB") {
     what = 1;
  } else if (op2=="MUL") {
     what = 2;
  } else if (op2=="DIV") {
     what = 3;
  } else {
    throw AipsError("Unrecognized operation");
  }

// Check rows

  const uInt nRows = left->nRow();
  if (right->nRow() != nRows) {
     throw (AipsError("Input Scan Tables must have the same number of rows"));
  }

// Input Tables and columns

  const Table& tLeft = left->table();
  const Table& tRight = right->table();
//
  ROArrayColumn<Float> tSysLeft(tLeft, "TSYS");
  ROArrayColumn<Float> tSysRight(tRight, "TSYS");

// Output Table cloned from input

  SDMemTable* pTabOut = new SDMemTable(*left, True);

// Loop over rows

  for (uInt i=0; i<nRows; i++) {

// Get data 
     MaskedArray<Float> mLeft(left->rowAsMaskedArray(i));
     MaskedArray<Float> mRight(right->rowAsMaskedArray(i));
//
     IPosition shpLeft = mLeft.shape();
     IPosition shpRight = mRight.shape();
     if (!shpLeft.isEqual(shpRight)) {
       throw(AipsError("left/right Scan Tables are not conformant"));
     }

// Get TSys

     Array<Float> tSysLeftArr, tSysRightArr;
     tSysLeft.get(i, tSysLeftArr);
     tSysRight.get(i, tSysRightArr);

// Make container

     SDContainer sc = left->getSDContainer(i);

// Operate on data and TSys

     if (what==0) {                               
        MaskedArray<Float> tmp = mLeft + mRight;
        putDataInSDC(sc, tmp.getArray(), tmp.getMask());
        sc.putTsys(tSysLeftArr+tSysRightArr);
     } else if (what==1) {
        MaskedArray<Float> tmp = mLeft - mRight;
        putDataInSDC(sc, tmp.getArray(), tmp.getMask());
        sc.putTsys(tSysLeftArr-tSysRightArr);
     } else if (what==2) {
        MaskedArray<Float> tmp = mLeft * mRight;
        putDataInSDC(sc, tmp.getArray(), tmp.getMask());
        sc.putTsys(tSysLeftArr*tSysRightArr);
     } else if (what==3) {
        MaskedArray<Float> tmp = mLeft / mRight;
        putDataInSDC(sc, tmp.getArray(), tmp.getMask());
        sc.putTsys(tSysLeftArr/tSysRightArr);
     }

// Put new row in output Table

     pTabOut->putSDContainer(sc);
  }
//
  return CountedPtr<SDMemTable>(pTabOut);
}



std::vector<float> SDMath::statistic(const CountedPtr<SDMemTable>& in,
				     const Vector<Bool>& mask,
				     const String& which, Int row) const
//
// Perhaps iteration over pol/beam/if should be in here
// and inside the nrow iteration ?
//
{
  const uInt nRow = in->nRow();

// Specify cursor location

  IPosition start, end;
  getCursorLocation(start, end, *in);

// Loop over rows

  const uInt nEl = mask.nelements();
  uInt iStart = 0;
  uInt iEnd = in->nRow()-1;
//  
  if (row>=0) {
     iStart = row;
     iEnd = row;
  }
//
  std::vector<float> result(iEnd-iStart+1);
  for (uInt ii=iStart; ii <= iEnd; ++ii) {

// Get row and deconstruct

     MaskedArray<Float> marr(in->rowAsMaskedArray(ii));
     Array<Float> arr = marr.getArray();
     Array<Bool> barr = marr.getMask();

// Access desired piece of data

     Array<Float> v((arr(start,end)).nonDegenerate());
     Array<Bool> m((barr(start,end)).nonDegenerate());

// Apply OTF mask

     MaskedArray<Float> tmp;
     if (m.nelements()==nEl) {
       tmp.setData(v,m&&mask);
     } else {
       tmp.setData(v,m);
     }

// Get statistic

     result[ii-iStart] = mathutil::statistics(which, tmp);
  }
//
  return result;
}


SDMemTable* SDMath::bin(const SDMemTable& in, Int width) const
{
  SDHeader sh = in.getSDHeader();
  SDMemTable* pTabOut = new SDMemTable(in, True);

// Bin up SpectralCoordinates

  IPosition factors(1);
  factors(0) = width;
  for (uInt j=0; j<in.nCoordinates(); ++j) {
    CoordinateSystem cSys;
    cSys.addCoordinate(in.getCoordinate(j));
    CoordinateSystem cSysBin = 
      CoordinateUtil::makeBinnedCoordinateSystem(factors, cSys, False);
//
    SpectralCoordinate sCBin = cSysBin.spectralCoordinate(0);
    pTabOut->setCoordinate(sCBin, j);
  }

// Use RebinLattice to find shape

  IPosition shapeIn(1,sh.nchan);
  IPosition shapeOut = RebinLattice<Float>::rebinShape(shapeIn, factors);
  sh.nchan = shapeOut(0);
  pTabOut->putSDHeader(sh);


// Loop over rows and bin along channel axis
  
  for (uInt i=0; i < in.nRow(); ++i) {
    SDContainer sc = in.getSDContainer(i);
//
    Array<Float> tSys(sc.getTsys());                           // Get it out before sc changes shape

// Bin up spectrum

    MaskedArray<Float> marr(in.rowAsMaskedArray(i));
    MaskedArray<Float> marrout;
    LatticeUtilities::bin(marrout, marr, asap::ChanAxis, width);

// Put back the binned data and flags

    IPosition ip2 = marrout.shape();
    sc.resize(ip2);
//
    putDataInSDC(sc, marrout.getArray(), marrout.getMask());

// Bin up Tsys.  

    Array<Bool> allGood(tSys.shape(),True);
    MaskedArray<Float> tSysIn(tSys, allGood, True);
//
    MaskedArray<Float> tSysOut;    
    LatticeUtilities::bin(tSysOut, tSysIn, asap::ChanAxis, width);
    sc.putTsys(tSysOut.getArray());
//
    pTabOut->putSDContainer(sc);
  }
  return pTabOut;
}

SDMemTable* SDMath::simpleOperate(const SDMemTable& in, Float val, Bool doAll,
				  uInt what) const
//
// what = 0   Multiply
//        1   Add
{
   SDMemTable* pOut = new SDMemTable(in,False);
   const Table& tOut = pOut->table();
   ArrayColumn<Float> spec(tOut,"SPECTRA");  
//
   if (doAll) {
      for (uInt i=0; i < tOut.nrow(); i++) {

// Get

         MaskedArray<Float> marr(pOut->rowAsMaskedArray(i));

// Operate

         if (what==0) {
            marr *= val;
         } else if (what==1) {
            marr += val;
         }

// Put

         spec.put(i, marr.getArray());
      }
   } else {

// Get cursor location

      IPosition start, end;
      getCursorLocation(start, end, in);
//
      for (uInt i=0; i < tOut.nrow(); i++) {

// Get

         MaskedArray<Float> dataIn(pOut->rowAsMaskedArray(i));

// Modify. More work than we would like to deal with the mask

         Array<Float>& values = dataIn.getRWArray();
         Array<Bool> mask(dataIn.getMask());
//
         Array<Float> values2 = values(start,end);
         Array<Bool> mask2 = mask(start,end);
         MaskedArray<Float> t(values2,mask2);
         if (what==0) {
            t *= val;
         } else if (what==1) {
            t += val;
         }
         values(start, end) = t.getArray();     // Write back into 'dataIn'

// Put
         spec.put(i, dataIn.getArray());
      }
   }
//
   return pOut;
}



SDMemTable* SDMath::averagePol(const SDMemTable& in, const Vector<Bool>& mask) const
//
// Average all polarizations together, weighted by variance
//
{
//   WeightType wtType = NONE;
//   convertWeightString(wtType, weight);

   const uInt nRows = in.nRow();
   const uInt polAxis = asap::PolAxis;                     // Polarization axis
   const uInt chanAxis = asap::ChanAxis;                    // Spectrum axis

// Create output Table and reshape number of polarizations

  Bool clear=True;
  SDMemTable* pTabOut = new SDMemTable(in, clear);
  SDHeader header = pTabOut->getSDHeader();
  header.npol = 1;
  pTabOut->putSDHeader(header);

// Shape of input and output data 

  const IPosition& shapeIn = in.rowAsMaskedArray(0u, False).shape();
  IPosition shapeOut(shapeIn);
  shapeOut(polAxis) = 1;                          // Average all polarizations
//
  const uInt nChan = shapeIn(chanAxis);
  const IPosition vecShapeOut(4,1,1,1,nChan);     // A multi-dim form of a Vector shape
  IPosition start(4), end(4);

// Output arrays

  Array<Float> outData(shapeOut, 0.0);
  Array<Bool> outMask(shapeOut, True);
  const IPosition axes(2, 2, 3);              // pol-channel plane
// 
  const Bool useMask = (mask.nelements() == shapeIn(chanAxis));

// Loop over rows

   for (uInt iRow=0; iRow<nRows; iRow++) {

// Get data for this row

      MaskedArray<Float> marr(in.rowAsMaskedArray(iRow));
      Array<Float>& arr = marr.getRWArray();
      const Array<Bool>& barr = marr.getMask();

// Make iterators to iterate by pol-channel planes

      ReadOnlyArrayIterator<Float> itDataPlane(arr, axes);
      ReadOnlyArrayIterator<Bool> itMaskPlane(barr, axes);

// Accumulations

      Float fac = 1.0;
      Vector<Float> vecSum(nChan,0.0);

// Iterate through data by pol-channel planes

      while (!itDataPlane.pastEnd()) {

// Iterate through plane by polarization  and accumulate Vectors

        Vector<Float> t1(nChan); t1 = 0.0;
        Vector<Bool> t2(nChan); t2 = True;
        MaskedArray<Float> vecSum(t1,t2);
        Float varSum = 0.0;
        {
           ReadOnlyVectorIterator<Float> itDataVec(itDataPlane.array(), 1);
           ReadOnlyVectorIterator<Bool> itMaskVec(itMaskPlane.array(), 1);
           while (!itDataVec.pastEnd()) {     

// Create MA of data & mask (optionally including OTF mask) and  get variance 

              if (useMask) {
                 const MaskedArray<Float> spec(itDataVec.vector(),mask&&itMaskVec.vector());
                 fac = 1.0 / variance(spec);
              } else {
                 const MaskedArray<Float> spec(itDataVec.vector(),itMaskVec.vector());
                 fac = 1.0 / variance(spec);
              }

// Normalize spectrum (without OTF mask) and accumulate

              const MaskedArray<Float> spec(fac*itDataVec.vector(), itMaskVec.vector());
              vecSum += spec;
              varSum += fac;

// Next

              itDataVec.next();
              itMaskVec.next();
           }
        }

// Normalize summed spectrum

        vecSum /= varSum;

// FInd position in input data array.  We are iterating by pol-channel
// plane so all that will change is beam and IF and that's what we want.

        IPosition pos = itDataPlane.pos();

// Write out data. This is a bit messy. We have to reform the Vector 
// accumulator into an Array of shape (1,1,1,nChan)

        start = pos;
        end = pos; 
        end(chanAxis) = nChan-1;
        outData(start,end) = vecSum.getArray().reform(vecShapeOut);
        outMask(start,end) = vecSum.getMask().reform(vecShapeOut);

// Step to next beam/IF combination

        itDataPlane.next();
        itMaskPlane.next();
      }

// Generate output container and write it to output table

      SDContainer sc = in.getSDContainer();
      sc.resize(shapeOut);
//
      putDataInSDC(sc, outData, outMask);
      pTabOut->putSDContainer(sc);
   }
//
  return pTabOut;
}


SDMemTable* SDMath::smooth(const SDMemTable& in, 
			   const casa::String& kernelType,
			   casa::Float width, Bool doAll) const
{

// Number of channels

   const uInt chanAxis = asap::ChanAxis;  // Spectral axis
   SDHeader sh = in.getSDHeader();
   const uInt nChan = sh.nchan;

// Generate Kernel

   VectorKernel::KernelTypes type = VectorKernel::toKernelType(kernelType);
   Vector<Float> kernel = VectorKernel::make(type, width, nChan, True, False);

// Generate Convolver

   IPosition shape(1,nChan);
   Convolver<Float> conv(kernel, shape);

// New Table

   SDMemTable* pTabOut = new SDMemTable(in,True);

// Get cursor location
         
  IPosition start, end;
  getCursorLocation(start, end, in);
//
  IPosition shapeOut(4,1);

// Output Vectors

  Vector<Float> valuesOut(nChan);
  Vector<Bool> maskOut(nChan);

// Loop over rows in Table

  for (uInt ri=0; ri < in.nRow(); ++ri) {

// Get copy of data
    
    const MaskedArray<Float>& dataIn(in.rowAsMaskedArray(ri));
    AlwaysAssert(dataIn.shape()(chanAxis)==nChan, AipsError);
//
    Array<Float> valuesIn = dataIn.getArray();
    Array<Bool> maskIn = dataIn.getMask();

// Branch depending on whether we smooth all locations or just
// those pointed at by the current selection cursor

    if (doAll) {
       uInt axis = asap::ChanAxis;
       VectorIterator<Float> itValues(valuesIn, axis);
       VectorIterator<Bool> itMask(maskIn, axis);
       while (!itValues.pastEnd()) {

// Smooth
          if (kernelType==VectorKernel::HANNING) {
             mathutil::hanning(valuesOut, maskOut, itValues.vector(), itMask.vector());
             itMask.vector() = maskOut;
          } else {
             mathutil::replaceMaskByZero(itValues.vector(), itMask.vector());
             conv.linearConv(valuesOut, itValues.vector());
          }
//
          itValues.vector() = valuesOut;
//
          itValues.next();
          itMask.next();
       }
    } else {

// Set multi-dim Vector shape

       shapeOut(chanAxis) = valuesIn.shape()(chanAxis);

// Stuff about with shapes so that we don't have conformance run-time errors

       Vector<Float> valuesIn2 = valuesIn(start,end).nonDegenerate();
       Vector<Bool> maskIn2 = maskIn(start,end).nonDegenerate();

// Smooth

       if (kernelType==VectorKernel::HANNING) {
          mathutil::hanning(valuesOut, maskOut, valuesIn2, maskIn2);
          maskIn(start,end) = maskOut.reform(shapeOut);
       } else {
          mathutil::replaceMaskByZero(valuesIn2, maskIn2);
          conv.linearConv(valuesOut, valuesIn2);
       }
//
       valuesIn(start,end) = valuesOut.reform(shapeOut);
    }

// Create and put back

    SDContainer sc = in.getSDContainer(ri);
    putDataInSDC(sc, valuesIn, maskIn);
//
    pTabOut->putSDContainer(sc);
  }
//
  return pTabOut;
}


SDMemTable* SDMath::convertFlux (const SDMemTable& in, Float a, Float eta, Bool doAll) const
// 
// As it is, this function could be implemented with 'simpleOperate'
// However, I anticipate that eventually we will look the conversion
// values up in a Table and apply them in a frequency dependent way,
// so I have implemented it fully here
//
{
  SDHeader sh = in.getSDHeader();
  SDMemTable* pTabOut = new SDMemTable(in, True);

// FInd out how to convert values into Jy and K (e.g. units might be mJy or mK)
// Also automatically find out what we are converting to according to the
// flux unit

  Unit fluxUnit(sh.fluxunit); 
  Unit K(String("K"));
  Unit JY(String("Jy"));
//
  Bool toKelvin = True;
  Double inFac = 1.0;
  if (fluxUnit==JY) {
     cerr << "Converting to K" << endl;
//
     Quantum<Double> t(1.0,fluxUnit);
     Quantum<Double> t2 = t.get(JY);
     inFac = (t2 / t).getValue();
//
     toKelvin = True;
     sh.fluxunit = "K";
  } else if (fluxUnit==K) {
     cerr << "Converting to Jy" << endl;
//
     Quantum<Double> t(1.0,fluxUnit);
     Quantum<Double> t2 = t.get(K);
     inFac = (t2 / t).getValue();
//
     toKelvin = False;
     sh.fluxunit = "Jy";
  } else {
     throw AipsError("Unrecognized brightness units in Table - must be consistent with Jy or K");
  }
  pTabOut->putSDHeader(sh);

// Compute conversion factor. 'a' and 'eta' are really frequency, time and  
// telescope dependent and should be looked// up in a table

  Float factor = 2.0 * inFac * 1.0e-7 * 1.0e26 * 
                 QC::k.getValue(Unit(String("erg/K"))) / a / eta;
  if (toKelvin) {
    factor = 1.0 / factor;
  }
  cerr << "Applying conversion factor = " << factor << endl;

// For operations only on specified cursor location

  IPosition start, end;
  getCursorLocation(start, end, in);

// Loop over rows and apply factor to spectra
  
  const uInt axis = asap::ChanAxis;
  for (uInt i=0; i < in.nRow(); ++i) {

// Get data

    MaskedArray<Float> dataIn(in.rowAsMaskedArray(i));
    Array<Float>& valuesIn = dataIn.getRWArray();              // writable reference
    const Array<Bool>& maskIn = dataIn.getMask();  

// Need to apply correct conversion factor (frequency and time dependent)
// which should be sourced from a Table. For now we just apply the given 
// factor to everything

    if (doAll) {
       VectorIterator<Float> itValues(valuesIn, asap::ChanAxis);
       while (!itValues.pastEnd()) {
          itValues.vector() *= factor;                            // Writes back into dataIn
//
          itValues.next();
       }
    } else {
       Array<Float> valuesIn2 = valuesIn(start,end);
       valuesIn2 *= factor;
       valuesIn(start,end) = valuesIn2;
    }

// Write out

    SDContainer sc = in.getSDContainer(i);
    putDataInSDC(sc, valuesIn, maskIn);
//
    pTabOut->putSDContainer(sc);
  }
  return pTabOut;
}



SDMemTable* SDMath::gainElevation (const SDMemTable& in, const Vector<Float>& coeffs,
                                   const String& fileName,
                                   const String& methodStr, Bool doAll) const
{

// Get header and clone output table

  SDHeader sh = in.getSDHeader();
  SDMemTable* pTabOut = new SDMemTable(in, True);

// Get elevation data from SDMemTable and convert to degrees

  const Table& tab = in.table();
  ROScalarColumn<Float> elev(tab, "ELEVATION");
  Vector<Float> x = elev.getColumn();
  x *= Float(180 / C::pi);
//
  const uInt nC = coeffs.nelements();
  if (fileName.length()>0 && nC>0) {
     throw AipsError("You must choose either polynomial coefficients or an ascii file, not both");
  }

// Correct

  if (nC>0 || fileName.length()==0) {

// Find instrument

     Bool throwIt = True;
     Instrument inst = SDMemTable::convertInstrument (sh.antennaname, throwIt);
     
// Set polynomial

     Polynomial<Float>* pPoly = 0;
     Vector<Float> coeff;
     String msg;
     if (nC>0) {
        pPoly = new Polynomial<Float>(nC);
        coeff = coeffs;
        msg = String("user");
     } else {
        if (inst==PKSMULTIBEAM) {
        } else if (inst==PKSSINGLEBEAM) {
        } else if (inst==TIDBINBILLA) {
           pPoly = new Polynomial<Float>(3);
           coeff.resize(3);
           coeff(0) = 3.58788e-1;
           coeff(1) = 2.87243e-2;
           coeff(2) = -3.219093e-4;
        } else if (inst==MOPRA) {
        }
        msg = String("built in");
     }
//
     if (coeff.nelements()>0) {
        pPoly->setCoefficients(coeff);
     } else {
        throw AipsError("There is no known gain-el polynomial known for this instrument");
     }
//
     cerr << "Making polynomial correction with " << msg << " coefficients" << endl;
     const uInt nRow = in.nRow();
     Vector<Float> factor(nRow);
     for (uInt i=0; i<nRow; i++) {
        factor[i] = (*pPoly)(x[i]);
     }
     delete pPoly;
//
     correctFromVector (pTabOut, in, doAll, factor);
  } else {

// Indicate which columns to read from ascii file

     String col0("ELEVATION");
     String col1("FACTOR");

// Read and correct

     cerr << "Making correction from ascii Table" << endl;
     correctFromAsciiTable (pTabOut, in, fileName, col0, col1, 
                            methodStr, doAll, x);
   }
//
   return pTabOut;
}

 

SDMemTable* SDMath::opacity (const SDMemTable& in, Float tau, Bool doAll) const
{

// Get header and clone output table

  SDHeader sh = in.getSDHeader();
  SDMemTable* pTabOut = new SDMemTable(in, True);

// Get elevation data from SDMemTable and convert to degrees

  const Table& tab = in.table();
  ROScalarColumn<Float> elev(tab, "ELEVATION");
  Vector<Float> zDist = elev.getColumn();
  zDist = Float(C::pi_2) - zDist;

// Generate correction factor

  const uInt nRow = in.nRow();
  Vector<Float> factor(nRow);
  Vector<Float> factor2(nRow);
  for (uInt i=0; i<nRow; i++) {
     factor[i] = exp(tau)/cos(zDist[i]);
  }

// Correct

  correctFromVector (pTabOut, in, doAll, factor);
//
  return pTabOut;
}




// 'private' functions

void SDMath::fillSDC(SDContainer& sc,
		     const Array<Bool>& mask,
		     const Array<Float>& data,
		     const Array<Float>& tSys,
		     Int scanID, Double timeStamp,
		     Double interval, const String& sourceName,
		     const Vector<uInt>& freqID) const
{
// Data and mask

  putDataInSDC(sc, data, mask);

// TSys

  sc.putTsys(tSys);

// Time things

  sc.timestamp = timeStamp;
  sc.interval = interval;
  sc.scanid = scanID;
//
  sc.sourcename = sourceName;
  sc.putFreqMap(freqID);
}

void SDMath::normalize(MaskedArray<Float>& sum,
                        const Array<Float>& sumSq,
                        const Array<Float>& nPts,
                        WeightType wtType, Int axis,
                        Int nAxesSub) const
{
   IPosition pos2(nAxesSub,0);
//
   if (wtType==NONE) {

// We just average by the number of points accumulated.
// We need to make a MA out of nPts so that no divide by
// zeros occur

      MaskedArray<Float> t(nPts, (nPts>Float(0.0)));
      sum /= t;
   } else if (wtType==VAR) {

// Normalize each spectrum by sum(1/var) where the variance
// is worked out for each spectrum

      Array<Float>& data = sum.getRWArray();
      VectorIterator<Float> itData(data, axis);
      while (!itData.pastEnd()) {
         pos2 = itData.pos().getFirst(nAxesSub);
         itData.vector() /= sumSq(pos2);
         itData.next();
      }
   } else if (wtType==TSYS) {
   }
}


void SDMath::accumulate(Double& timeSum, Double& intSum, Int& nAccum,
			MaskedArray<Float>& sum, Array<Float>& sumSq, 
			Array<Float>& nPts, Array<Float>& tSysSum, 
			const Array<Float>& tSys, const Array<Float>& nInc, 
			const Vector<Bool>& mask, Double time, Double interval,
			const Block<CountedPtr<SDMemTable> >& in,
			uInt iTab, uInt iRow, uInt axis, 
			uInt nAxesSub, Bool useMask,
			WeightType wtType) const
{

// Get data

   MaskedArray<Float> dataIn(in[iTab]->rowAsMaskedArray(iRow));
   Array<Float>& valuesIn = dataIn.getRWArray();           // writable reference
   const Array<Bool>& maskIn = dataIn.getMask();          // RO reference
//
   if (wtType==NONE) {
      const MaskedArray<Float> n(nInc,dataIn.getMask());
      nPts += n;                               // Only accumulates where mask==T
   } else if (wtType==VAR) {

// We are going to average the data, weighted by the noise for each pol, beam and IF.
// So therefore we need to iterate through by spectrum (axis 3)

      VectorIterator<Float> itData(valuesIn, axis);
      ReadOnlyVectorIterator<Bool> itMask(maskIn, axis);
      Float fac = 1.0;
      IPosition pos(nAxesSub,0);  
//
      while (!itData.pastEnd()) {

// Make MaskedArray of Vector, optionally apply OTF mask, and find scaling factor

        if (useMask) {
           MaskedArray<Float> tmp(itData.vector(),mask&&itMask.vector());
           fac = 1.0/variance(tmp);
        } else {
           MaskedArray<Float> tmp(itData.vector(),itMask.vector());
           fac = 1.0/variance(tmp);
        }

// Scale data

        itData.vector() *= fac;     // Writes back into 'dataIn'
//
// Accumulate variance per if/pol/beam averaged over spectrum
// This method to get pos2 from itData.pos() is only valid
// because the spectral axis is the last one (so we can just
// copy the first nAXesSub positions out)

        pos = itData.pos().getFirst(nAxesSub);
        sumSq(pos) += fac;
//
        itData.next();
        itMask.next();
      }
   } else if (wtType==TSYS) {
   }

// Accumulate sum of (possibly scaled) data

   sum += dataIn;

// Accumulate Tsys, time, and interval

   tSysSum += tSys;
   timeSum += time;
   intSum += interval;
   nAccum += 1;
}




void SDMath::getCursorLocation(IPosition& start, IPosition& end,
			       const SDMemTable& in) const
{
  const uInt nDim = 4;
  const uInt i = in.getBeam();
  const uInt j = in.getIF();
  const uInt k = in.getPol();
  const uInt n = in.nChan();
//
  start.resize(nDim);
  start(0) = i;
  start(1) = j;
  start(2) = k;
  start(3) = 0;
//
  end.resize(nDim);
  end(0) = i;
  end(1) = j;
  end(2) = k;
  end(3) = n-1;
}


void SDMath::convertWeightString(WeightType& wtType, const String& weightStr) const
{
  String tStr(weightStr);
  tStr.upcase();
  if (tStr.contains(String("NONE"))) {
     wtType = NONE;
  } else if (tStr.contains(String("VAR"))) {
     wtType = VAR;
  } else if (tStr.contains(String("TSYS"))) {
     wtType = TSYS;
     throw(AipsError("T_sys weighting not yet implemented"));
  } else {
    throw(AipsError("Unrecognized weighting type"));
  }
}

void SDMath::convertInterpString(Int& type, const String& interp) const
{
  String tStr(interp);
  tStr.upcase();
  if (tStr.contains(String("NEAR"))) {
     type = InterpolateArray1D<Float,Float>::nearestNeighbour;
  } else if (tStr.contains(String("LIN"))) {
     type = InterpolateArray1D<Float,Float>::linear;
  } else if (tStr.contains(String("CUB"))) {
     type = InterpolateArray1D<Float,Float>::cubic;
  } else if (tStr.contains(String("SPL"))) {
     type = InterpolateArray1D<Float,Float>::spline;
  } else {
    throw(AipsError("Unrecognized interpolation type"));
  }
}

void SDMath::putDataInSDC(SDContainer& sc, const Array<Float>& data,
			  const Array<Bool>& mask) const
{
    sc.putSpectrum(data);
//
    Array<uChar> outflags(data.shape());
    convertArray(outflags,!mask);
    sc.putFlags(outflags);
}

Table SDMath::readAsciiFile (const String& fileName) const
{
   String formatString;
   Table tbl = readAsciiTable (formatString, Table::Memory, fileName, "", "", False);
   return tbl;
}



void SDMath::correctFromAsciiTable(SDMemTable* pTabOut,
                                   const SDMemTable& in, const String& fileName,
                                   const String& col0, const String& col1,
                                   const String& methodStr, Bool doAll,
                                   const Vector<Float>& xOut) const
{

// Read gain-elevation ascii file data into a Table.

  Table geTable = readAsciiFile (fileName);
//
  correctFromTable (pTabOut, in, geTable, col0, col1, methodStr, doAll, xOut);
}

void SDMath::correctFromTable(SDMemTable* pTabOut, const SDMemTable& in, 
                              const Table& tTable, const String& col0, 
                              const String& col1,
                              const String& methodStr, Bool doAll,
                              const Vector<Float>& xOut) const
{

// Get data from Table

  ROScalarColumn<Float> geElCol(tTable, col0);
  ROScalarColumn<Float> geFacCol(tTable, col1);
  Vector<Float> xIn = geElCol.getColumn();
  Vector<Float> yIn = geFacCol.getColumn();
  Vector<Bool> maskIn(xIn.nelements(),True);

// Interpolate (and extrapolate) with desired method

   Int method = 0;
   convertInterpString(method, methodStr);
//
   Vector<Float> yOut;
   Vector<Bool> maskOut;
   InterpolateArray1D<Float,Float>::interpolate(yOut, maskOut, xOut, 
                                                xIn, yIn, maskIn, method,
                                                True, True);
// Apply 

   correctFromVector (pTabOut, in, doAll, yOut);
}


void SDMath::correctFromVector (SDMemTable* pTabOut, const SDMemTable& in, 
                                Bool doAll, const Vector<Float>& factor) const
{
// For operations only on specified cursor location

  IPosition start, end;
  getCursorLocation(start, end, in);

// Loop over rows and interpolate correction factor
  
  const uInt axis = asap::ChanAxis;
  for (uInt i=0; i < in.nRow(); ++i) {

// Get data

    MaskedArray<Float> dataIn(in.rowAsMaskedArray(i));
    Array<Float>& valuesIn = dataIn.getRWArray();  
    const Array<Bool>& maskIn = dataIn.getMask();  

// Apply factor

    if (doAll) {
       VectorIterator<Float> itValues(valuesIn, asap::ChanAxis);
       while (!itValues.pastEnd()) {
          itValues.vector() *= factor(i);
          itValues.next();
       }
    } else {
       Array<Float> valuesIn2 = valuesIn(start,end);
       valuesIn2 *= factor(i);
       valuesIn(start,end) = valuesIn2;
    }

// Write out

    SDContainer sc = in.getSDContainer(i);
    putDataInSDC(sc, valuesIn, maskIn);
//
    pTabOut->putSDContainer(sc);
  }
}


