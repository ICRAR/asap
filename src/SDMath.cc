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
#include <casa/Exceptions.h>

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
//using namespace asap::SDMath;

CountedPtr<SDMemTable> SDMath::average(const CountedPtr<SDMemTable>& in) 
//
// Average all rows in Table in time
//
{
  Table t = in->table();
  ROArrayColumn<Float> tsys(t, "TSYS");
  ROScalarColumn<Double> mjd(t, "TIME");
  ROScalarColumn<String> srcn(t, "SRCNAME");
  ROScalarColumn<Double> integr(t, "INTERVAL");
  ROArrayColumn<uInt> freqidc(t, "FREQID");
  IPosition ip = in->rowAsMaskedArray(0).shape();
  Array<Float> outarr(ip); outarr =0.0;
  Array<Float> narr(ip);narr = 0.0;
  Array<Float> narrinc(ip);narrinc = 1.0;

  Array<Float> tsarr(tsys.shape(0));
  Array<Float> outtsarr(tsys.shape(0));
  outtsarr =0.0;
  tsys.get(0, tsarr);// this is probably unneccessary as tsys should
  Double tme = 0.0;
  Double inttime = 0.0;

// Loop over rows

  for (uInt i=0; i < t.nrow(); i++) {

// Get data and accumulate sums

    MaskedArray<Float> marr(in->rowAsMaskedArray(i));
    outarr += marr;
    MaskedArray<Float> n(narrinc,marr.getMask());
    narr += n;

// Accumulkate Tsys

    tsys.get(i, tsarr);// this is probably unneccessary as tsys should
    outtsarr += tsarr; // be constant
    Double tmp;
    mjd.get(i,tmp);
    tme += tmp;// average time
    integr.get(i,tmp);
    inttime += tmp;
  }

// Average

  MaskedArray<Float> nma(narr,(narr > Float(0)));
  outarr /= nma;

// Create container and put 

  Array<Bool> outflagsb = !(nma.getMask());
  Array<uChar> outflags(outflagsb.shape());
  convertArray(outflags,outflagsb);
  SDContainer sc = in->getSDContainer();
  Int n = t.nrow();
  outtsarr /= Float(n);
  sc.timestamp = tme/Double(n);
  sc.interval = inttime;
  String tstr; srcn.getScalar(0,tstr);// get sourcename of "mid" point
  sc.sourcename = tstr;
  Vector<uInt> tvec;
  freqidc.get(0,tvec);
  sc.putFreqMap(tvec);
  sc.putTsys(outtsarr);
  sc.scanid = 0;
  sc.putSpectrum(outarr);
  sc.putFlags(outflags);
  SDMemTable* sdmt = new SDMemTable(*in,True);
  sdmt->putSDContainer(sc);
  return CountedPtr<SDMemTable>(sdmt);
}

CountedPtr<SDMemTable>
SDMath::quotient(const CountedPtr<SDMemTable>& on,
                 const CountedPtr<SDMemTable>& off) 
//
// Compute quotient spectrum
//
{
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

  SDMemTable* sdmt = new SDMemTable(*on, True);

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
     Array<Bool> outflagsb = !(mon.getMask() && moff.getMask());
     Array<uChar> outflags(outflagsb.shape());
     convertArray(outflags,outflagsb);

// Fill container for this row

     SDContainer sc = on->getSDContainer();
     sc.putTsys(tsarr);
     sc.scanid = 0;
     sc.putSpectrum(out);
     sc.putFlags(outflags);

// Put new row in output Table

     sdmt->putSDContainer(sc);
  }
//
  return CountedPtr<SDMemTable>(sdmt);
}

CountedPtr<SDMemTable>
SDMath::multiply(const CountedPtr<SDMemTable>& in, Float factor) 
//
// Multiply values by factor
//
{
  SDMemTable* sdmt = new SDMemTable(*in);
  Table t = sdmt->table();
  ArrayColumn<Float> spec(t,"SPECTRA");

  for (uInt i=0; i < t.nrow(); i++) {
    MaskedArray<Float> marr(sdmt->rowAsMaskedArray(i));
    marr *= factor;
    spec.put(i, marr.getArray());
  }
  return CountedPtr<SDMemTable>(sdmt);
}

CountedPtr<SDMemTable>
SDMath::add(const CountedPtr<SDMemTable>& in, Float offset) 
//
// Add offset to values
//
{
  SDMemTable* sdmt = new SDMemTable(*in);

  Table t = sdmt->table();
  ArrayColumn<Float> spec(t,"SPECTRA");

  for (uInt i=0; i < t.nrow(); i++) {
    MaskedArray<Float> marr(sdmt->rowAsMaskedArray(i));
    marr += offset;
    spec.put(i, marr.getArray());
  }
  return CountedPtr<SDMemTable>(sdmt);
}


CountedPtr<SDMemTable>
SDMath::hanning(const CountedPtr<SDMemTable>& in) 
//
// Hanning smooth each row
// Should Tsys be smoothed ?
//
{
  SDMemTable* sdmt = new SDMemTable(*in,True);

// Loop over rows in Table

  for (uInt ri=0; ri < in->nRow(); ++ri) {

// Get data
    
    const MaskedArray<Float>& marr(in->rowAsMaskedArray(ri));
    Array<Float> arr = marr.getArray();
    Array<Bool> barr = marr.getMask();

// Smooth along the channels axis

    uInt axis = 3;
    VectorIterator<Float> itData(arr, axis);
    VectorIterator<Bool> itMask(barr, axis);
    Vector<Float> outv;
    Vector<Bool> outm;
    while (!itData.pastEnd()) {
       mathutil::hanning(outv, outm, itData.vector(), itMask.vector());
       itData.vector() = outv;                                                 
       itMask.vector() = outm;
//
       itData.next();
       itMask.next();
    }

// Create and put back

    Array<uChar> outflags(barr.shape());
    convertArray(outflags,!barr);
    SDContainer sc = in->getSDContainer(ri);
    sc.putSpectrum(arr);
    sc.putFlags(outflags);
    sdmt->putSDContainer(sc);
  }
  return CountedPtr<SDMemTable>(sdmt);
}



CountedPtr<SDMemTable> SDMath::averages(const Block<CountedPtr<SDMemTable> >& in,
                                        const Vector<Bool>& mask) 
//
// Noise weighted averaging of spectra from many Tables.  Tables can have different
// number of rows.  
//
{

// Setup

  const uInt axis = 3;                                 // Spectral axis
  IPosition shp = in[0]->rowAsMaskedArray(0).shape();
  Array<Float> arr(shp);
  Array<Bool> barr(shp);
  Double sumInterval = 0.0;
  const Bool useMask = (mask.nelements() == shp(axis));

// Create data accumulation MaskedArray. We accumulate for each
// channel,if,pol,beam

  Array<Float> zero(shp); zero=0.0;
  Array<Bool> good(shp); good = True;
  MaskedArray<Float> sum(zero,good);

// Create accumulation Array for variance. We accumulate for 
// each if,pol,beam, but average over channel

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
  IPosition pos2(nAxesSub,0);                        // FOr indexing
//  
  Float fac = 1.0;
  const uInt nTables = in.nelements();
  for (uInt iTab=0; iTab<nTables; iTab++) {
     const uInt nRows = in[iTab]->nRow();
     sumInterval += nRows * in[iTab]->getInterval();   // Sum of time intervals
//
     for (uInt iRow=0; iRow<nRows; iRow++) {

// Check conforms

        IPosition shp2 = in[iTab]->rowAsMaskedArray(iRow).shape();
        if (!shp.isEqual(shp2)) {
           throw (AipsError("Shapes for all rows must be the same"));
        }

// Get data and deconstruct
 
        MaskedArray<Float> marr(in[iTab]->rowAsMaskedArray(iRow));
        Array<Float>& arr = marr.getRWArray();                     // writable reference
        const Array<Bool>& barr = marr.getMask();                  // RO reference

// We are going to average the data, weighted by the noise for each
// pol, beam and IF. So therefore we need to iterate through by 
// spectra (axis 3)

        VectorIterator<Float> itData(arr, axis);
        ReadOnlyVectorIterator<Bool> itMask(barr, axis);
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

          itData.vector() *= fac;

// Accumulate variance per if/pol/beam averaged over spectrum
// This method to get pos2 from itData.pos() is only valid
// because the spectral axis is the last one (so we can just
// copy the first nAXesSub positions out)

          pos2 = itData.pos().getFirst(nAxesSub);
          sumSq(pos2) += fac;
//
          itData.next();
          itMask.next();
        }    

// Accumulate sums

       sum += marr;
    }
  }

// Normalize by the sum of the 1/var.  

  Array<Float>& data = sum.getRWArray();    
  VectorIterator<Float> itData(data, axis);
  while (!itData.pastEnd()) {
     pos2 = itData.pos().getFirst(nAxesSub);           // See comments above
     itData.vector() /= sumSq(pos2);
     itData.next();
  }    

// Create and fill output

  Array<uChar> outflags(shp);
  convertArray(outflags,!(sum.getMask()));
//
  SDContainer sc = in[0]->getSDContainer();     // CLone from first container of first Table
  sc.putSpectrum(data);
  sc.putFlags(outflags);
  sc.interval = sumInterval;
//
  SDMemTable* sdmt = new SDMemTable(*in[0],True);  // CLone from first Table
  sdmt->putSDContainer(sc);
  return CountedPtr<SDMemTable>(sdmt);
}


CountedPtr<SDMemTable>
SDMath::averagePol(const CountedPtr<SDMemTable>& in,
                   const Vector<Bool>& mask) 
{
   const uInt nRows = in->nRow();
   const uInt axis = 3;                        // Spectrum
   const IPosition axes(2, 2, 3);              // pol-channel plane

// Create output Table

  SDMemTable* sdmt = new SDMemTable(*in, True);

// Loop over rows

   for (uInt iRow=0; iRow<nRows; iRow++) {

// Get data for this row

      MaskedArray<Float> marr(in->rowAsMaskedArray(iRow));
      Array<Float>& arr = marr.getRWArray();
      const Array<Bool>& barr = marr.getMask();
//
      IPosition shp = marr.shape();
      const Bool useMask = (mask.nelements() == shp(axis));
      const uInt nChan = shp(axis);

// Make iterators to iterate by pol-channel planes

     ArrayIterator<Float> itDataPlane(arr, axes);
     ReadOnlyArrayIterator<Bool> itMaskPlane(barr, axes);

// Accumulations

     Float fac = 0.0;
     Vector<Float> vecSum(nChan,0.0);

// Iterate by plane

     while (!itDataPlane.pastEnd()) {

// Iterate through pol-channel plane by spectrum

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

// We have formed the weighted averaged spectrum from all polarizations
// for this beam and IF.  Now replicate the spectrum to all polarizations 

        {
           VectorIterator<Float> itDataVec(itDataPlane.array(), 1);  // Writes back into 'arr'
           const Vector<Float>& vecSumData = vecSum.getArray();      // It *is* a Vector
//
           while (!itDataVec.pastEnd()) {     
              itDataVec.vector() = vecSumData;
              itDataVec.next();
           }
        }

// Step to next beam/IF combination

        itDataPlane.next();
        itMaskPlane.next();
     }

// Generate output container and write it to output table

     SDContainer sc = in->getSDContainer();
     Array<uChar> outflags(barr.shape());
     convertArray(outflags,!barr);
     sc.putSpectrum(arr);
     sc.putFlags(outflags);
     sdmt->putSDContainer(sc);
   }
//
  return CountedPtr<SDMemTable>(sdmt);
}


CountedPtr<SDMemTable> SDMath::bin(const CountedPtr<SDMemTable>& in,
                                   Int width) 
{
  SDHeader sh = in->getSDHeader();
  SDMemTable* sdmt = new SDMemTable(*in,True);

// Bin up SpectralCoordinates

  IPosition factors(1);
  factors(0) = width;
  for (uInt j=0; j<in->nCoordinates(); ++j) {
    CoordinateSystem cSys;
    cSys.addCoordinate(in->getCoordinate(j));
    CoordinateSystem cSysBin = 
      CoordinateUtil::makeBinnedCoordinateSystem (factors, cSys, False);
//
    SpectralCoordinate sCBin = cSysBin.spectralCoordinate(0);
    sdmt->setCoordinate(sCBin, j);
  }

// Use RebinLattice to find shape

  IPosition shapeIn(1,sh.nchan);
  IPosition shapeOut = RebinLattice<Float>::rebinShape (shapeIn, factors);
  sh.nchan = shapeOut(0);
  sdmt->putSDHeader(sh);


// Loop over rows and bin along channel axis
  
  const uInt axis = 3;
  for (uInt i=0; i < in->nRow(); ++i) {
    SDContainer sc = in->getSDContainer(i);
//
    Array<Float> tSys(sc.getTsys());                           // Get it out before sc changes shape

// Bin up spectrum

    MaskedArray<Float> marr(in->rowAsMaskedArray(i));
    MaskedArray<Float> marrout;
    LatticeUtilities::bin(marrout, marr, axis, width);

// Put back the binned data and flags

    IPosition ip2 = marrout.shape();
    sc.resize(ip2);
    sc.putSpectrum(marrout.getArray());
//
    Array<uChar> outflags(ip2);
    convertArray(outflags,!(marrout.getMask()));
    sc.putFlags(outflags);

// Bin up Tsys.  

    Array<Bool> allGood(tSys.shape(),True);
    MaskedArray<Float> tSysIn(tSys, allGood, True);
//
    MaskedArray<Float> tSysOut;    
    LatticeUtilities::bin(tSysOut, tSysIn, axis, width);
    sc.putTsys(tSysOut.getArray());
    sdmt->putSDContainer(sc);
  }
  return CountedPtr<SDMemTable>(sdmt);
}



std::vector<float> SDMath::statistic (const CountedPtr<SDMemTable>& in,
                                       const std::vector<bool>& mask,
                                       const std::string& which)
//
// Perhaps iteration over pol/beam/if should be in here
// and inside the nrow iteration ?
//
{
  const uInt nRow = in->nRow();
  std::vector<float> result(nRow);
  Vector<Bool> msk(mask);

// Specify cursor location

  uInt i = in->getBeam();
  uInt j = in->getIF();
  uInt k = in->getPol();
  IPosition start(4,i,j,k,0);
  IPosition end(4,i,j,k,in->nChan()-1);

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

     result[ii] = SDMath::theStatistic(which, tmp);
  }
//
  return result;
}


float SDMath::theStatistic(const std::string& which,  const casa::MaskedArray<Float>& data)
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
