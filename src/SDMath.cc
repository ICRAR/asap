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
#include <casa/Utilities/Assert.h>
#include <casa/Exceptions.h>

#include <scimath/Mathematics/VectorKernel.h>
#include <scimath/Mathematics/Convolver.h>

#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>

#include <lattices/Lattices/LatticeUtilities.h>
#include <lattices/Lattices/RebinLattice.h>
#include <coordinates/Coordinates/SpectralCoordinate.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <coordinates/Coordinates/CoordinateUtil.h>

#include "MathUtils.h"
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
				       const std::string& weightStr)
//
// Weighted averaging of spectra from one or more Tables.
//
{

// Convert weight type
  
  WeightType wtType = NONE;
  convertWeightString(wtType, weightStr);

// Create output Table by cloning from the first table

  SDMemTable* pTabOut = new SDMemTable(*in[0],True);

// Setup

  const uInt axis = 3;                                     // Spectral axis
  IPosition shp = in[0]->rowAsMaskedArray(0).shape();      // Must not change
  Array<Float> arr(shp);
  Array<Bool> barr(shp);
  const Bool useMask = (mask.nelements() == shp(axis));

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
     if (i!=axis) {
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

// Loop over tables

  Float fac = 1.0;
  const uInt nTables = in.nelements();
  for (uInt iTab=0; iTab<nTables; iTab++) {

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

           normalize(sum, sumSq, nPts, wtType, axis, nAxesSub);

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
	   nPts = 0.0; // reset this too!!!

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
                    tSys, nInc, mask, time, interval, in, iTab, iRow, axis, 
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
  normalize(sum, sumSq, nPts, wtType, axis, nAxesSub);

// Create and fill container.  The container we clone will be from
// the last Table and the first row that went into the current
// accumulation.  It probably doesn't matter that much really...

  Float nR(nAccum);
  SDContainer sc = in[tableStart]->getSDContainer(rowStart);
  fillSDC(sc, sum.getMask(), sum.getArray(), tSysSum/nR, outScanID,
           timeSum/nR, intSum, sourceNameStart, freqIDStart);
//
  pTabOut->putSDContainer(sc);
/*
   cout << endl;
   cout << "Last accumulation for output scan ID " << outScanID << endl;
   cout << "   The first row in this accumulation is " << rowStart << endl;
   cout << "   The number of rows accumulated is " << nAccum << endl;
   cout << "   The first table in this accumulation is " << tableStart << endl;
   cout << "   The first source in this accumulation is " << sourceNameStart << endl;
   cout << "   The first freqID in this accumulation is " << freqIDStart << endl;
   cout  << "   Average time stamp = " << timeSum/nR << endl;
   cout << "   Integrated time = " << intSum << endl;
*/
  return CountedPtr<SDMemTable>(pTabOut);
}



CountedPtr<SDMemTable>
SDMath::quotient(const CountedPtr<SDMemTable>& on,
		 const CountedPtr<SDMemTable>& off) 
{
//
// Compute quotient spectrum
//
  const uInt nRows = on->nRow();
  if (off->nRow() != nRows) {
     throw (AipsError("Input Scan Tables must have the same number of rows"));
  }

// Input Tables and columns

  Table ton = on->table();
  Table toff = off->table();
  ROArrayColumn<Float> tsys(toff, "TSYS");
  ROScalarColumn<Double> mjd(ton, "TIME");
  ROScalarColumn<Double> integr(ton, "INTERVAL");
  ROScalarColumn<String> srcn(ton, "SRCNAME");
  ROArrayColumn<uInt> freqidc(ton, "FREQID");

// Output Table cloned from input

  SDMemTable* pTabOut = new SDMemTable(*on, True);

// Loop over rows

  for (uInt i=0; i<nRows; i++) {
     MaskedArray<Float> mon(on->rowAsMaskedArray(i));
     MaskedArray<Float> moff(off->rowAsMaskedArray(i));
     IPosition ipon = mon.shape();
     IPosition ipoff = moff.shape();
//
     Array<Float> tsarr;  
     tsys.get(i, tsarr);
     if (ipon != ipoff && ipon != tsarr.shape()) {
       throw(AipsError("on/off not conformant"));
     }

// Compute quotient

     MaskedArray<Float> tmp = (mon-moff);
     Array<Float> out(tmp.getArray());
     out /= moff;
     out *= tsarr;
     Array<Bool> outflagsb = mon.getMask() && moff.getMask();

// Fill container for this row

     SDContainer sc = on->getSDContainer(i);
//
     putDataInSDC(sc, out, outflagsb);
     sc.putTsys(tsarr);
     sc.scanid = i;

// Put new row in output Table

     pTabOut->putSDContainer(sc);
  }
//
  return CountedPtr<SDMemTable>(pTabOut);
}



std::vector<float> SDMath::statistic(const CountedPtr<SDMemTable>& in,
				     const std::vector<bool>& mask,
				     const String& which)
//
// Perhaps iteration over pol/beam/if should be in here
// and inside the nrow iteration ?
//
{
  const uInt nRow = in->nRow();
  std::vector<float> result(nRow);
  Vector<Bool> msk(mask);

// Specify cursor location

  IPosition start, end;
  getCursorLocation(start, end, *in);

// Loop over rows

  const uInt nEl = msk.nelements();
  for (uInt ii=0; ii < in->nRow(); ++ii) {

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
       tmp.setData(v,m&&msk);
     } else {
       tmp.setData(v,m);
     }

// Get statistic

     result[ii] = mathutil::statistics(which, tmp);
  }
//
  return result;
}


SDMemTable* SDMath::bin(const SDMemTable& in, Int width)
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
  
  const uInt axis = 3;
  for (uInt i=0; i < in.nRow(); ++i) {
    SDContainer sc = in.getSDContainer(i);
//
    Array<Float> tSys(sc.getTsys());                           // Get it out before sc changes shape

// Bin up spectrum

    MaskedArray<Float> marr(in.rowAsMaskedArray(i));
    MaskedArray<Float> marrout;
    LatticeUtilities::bin(marrout, marr, axis, width);

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
    LatticeUtilities::bin(tSysOut, tSysIn, axis, width);
    sc.putTsys(tSysOut.getArray());
//
    pTabOut->putSDContainer(sc);
  }
  return pTabOut;
}

SDMemTable* SDMath::simpleOperate(const SDMemTable& in, Float val, Bool doAll,
				  uInt what)
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



SDMemTable* SDMath::averagePol(const SDMemTable& in, const Vector<Bool>& mask)
//
// Average all polarizations together, weighted by variance
//
{
//   WeightType wtType = NONE;
//   convertWeightString(wtType, weight);

   const uInt nRows = in.nRow();
   const uInt polAxis = 2;                     // Polarization axis
   const uInt chanAxis = 3;                    // Spectrum axis

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
			   casa::Float width, Bool doAll)
{

// Number of channels

   const uInt chanAxis = 3;                                     // Spectral axis
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
       uInt axis = 3;
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





// 'private' functions

void SDMath::fillSDC(SDContainer& sc,
		     const Array<Bool>& mask,
		     const Array<Float>& data,
		     const Array<Float>& tSys,
		     Int scanID, Double timeStamp,
		     Double interval, const String& sourceName,
		     const Vector<uInt>& freqID)
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
                        Int nAxesSub)
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
			WeightType wtType)
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
			       const SDMemTable& in)
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


void SDMath::convertWeightString(WeightType& wtType, const std::string& weightStr)
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

void SDMath::putDataInSDC(SDContainer& sc, const Array<Float>& data,
			  const Array<Bool>& mask)
{
    sc.putSpectrum(data);
//
    Array<uChar> outflags(data.shape());
    convertArray(outflags,!mask);
    sc.putFlags(outflags);
}
