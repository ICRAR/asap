//#---------------------------------------------------------------------------
//# SDContainer.cc: A container class for single dish integrations
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
#include <casa/iostream.h>
#include <casa/iomanip.h>
#include <casa/Exceptions.h>
#include <casa/Utilities/Assert.h>
#include <tables/Tables/Table.h>
#include <casa/Arrays/IPosition.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/VectorIter.h>
#include <casa/Arrays/ArrayAccessor.h>
#include <casa/BasicMath/Math.h>
#include <casa/Quanta/MVTime.h>



#include "SDDefs.h"
#include "SDContainer.h"

using namespace casa;
using namespace asap;

void SDHeader::print() const {
  MVTime mvt(this->utc);
  mvt.setFormat(MVTime::YMD);
  cout << "Observer: " << this->observer << endl 
       << "Project: " << this->project << endl
       << "Obstype: " << this->obstype << endl
       << "Antenna: " << this->antennaname << endl
       << "Ant. Position: " << this->antennaposition << endl
       << "Equinox: " << this->equinox << endl
       << "Freq. ref.: " << this->freqref << endl
       << "Ref. frequency: " << this->reffreq << endl
       << "Bandwidth: "  << this->bandwidth << endl
       << "Time (utc): " 
       << mvt
       << endl;
  //setprecision(10) << this->utc << endl;
}


SDContainer::SDContainer(uInt nBeam, uInt nIF, uInt nPol, uInt nChan) 
  : nBeam_(nBeam),
    nIF_(nIF),
    nPol_(nPol),
    nChan_(nChan),
    spectrum_(IPosition(4,nBeam,nIF,nPol,nChan)),
    flags_(IPosition(4,nBeam,nIF,nPol,nChan)),
    tsys_(IPosition(4,nBeam,nIF,nPol,nChan)),
    freqidx_(nIF),
    restfreqidx_(nIF),
    direction_(IPosition(2,nBeam,2)) {
  uChar x = 0;
  flags_ = ~x;
  tcal.resize(2);
}

SDContainer::SDContainer(IPosition shp) 
  : nBeam_(shp(0)),
    nIF_(shp(1)),
    nPol_(shp(2)),
    nChan_(shp(3)),
    spectrum_(shp),
    flags_(shp),
    tsys_(shp),
    freqidx_(shp(1)),
    restfreqidx_(shp(1)) {
  IPosition ip(2,shp(0),2);
  direction_.resize(ip);
  uChar x = 0;
  flags_ = ~x;
  tcal.resize(2);
}

SDContainer::~SDContainer() {
}

Bool SDContainer::resize(IPosition shp) {
  nBeam_ = shp(0);
  nIF_ = shp(1);
  nPol_ = shp(2);
  nChan_ = shp(3);
  spectrum_.resize(shp);
  flags_.resize(shp);
  tsys_.resize(shp);
  freqidx_.resize(shp(1));
  restfreqidx_.resize(shp(1));
  IPosition ip(2,shp(0),2);
  direction_.resize(ip);
}

Bool SDContainer::putSpectrum(const Array<Float>& spec) {
  spectrum_ = spec;
}
Bool SDContainer::putFlags(const Array<uChar>& flag) {
  flags_ = flag;
}
Bool SDContainer::putTsys(const Array<Float>& tsys) {
  tsys_ = tsys;
}

Bool SDContainer::setSpectrum(const Matrix<Float>& spec,
			      uInt whichBeam, uInt whichIF) 
//
// spec is [nChan,nPol] 
// spectrum_ is [,,,nChan]
// How annoying.
//
{

// Get slice and check dim

  IPosition start, end;
  setSlice (start, end, spec.shape(), spectrum_.shape(),
            whichBeam, whichIF, False);

// Get a reference to the Pol/Chan slice we are interested in

  Array<Float> subArr = spectrum_(start,end);

// Iterate through it and fill

  ReadOnlyVectorIterator<Float> inIt(spec,0);
  VectorIterator<Float> outIt(subArr,asap::ChanAxis);
  while (!inIt.pastEnd()) {
     outIt.vector() = inIt.vector();
//
     inIt.next();
     outIt.next();
  }

  // unset flags for this spectrum, they might be set again by the
  // setFlags method

  Array<uChar> arr(flags_(start,end));
  arr = uChar(0);
//
  return True;
}


Bool SDContainer::setFlags (const Matrix<uChar>& flags,
                            uInt whichBeam, uInt whichIF) 
//
// flags is [nChan,nPol] 
// flags_ is [,,,nChan]
// How annoying.
//
{
// Get slice and check dim

  IPosition start, end;
  setSlice (start, end, flags.shape(), flags_.shape(),
            whichBeam, whichIF, False);

// Get a reference to the Pol/Chan slice we are interested in

  Array<uChar> subArr = flags_(start,end);

// Iterate through it and fill

  ReadOnlyVectorIterator<uChar> inIt(flags,0);
  VectorIterator<uChar> outIt(subArr,asap::ChanAxis);
  while (!inIt.pastEnd()) {
     outIt.vector() = inIt.vector();
//
     inIt.next();
     outIt.next();
  }
//
  return True;
}


Bool SDContainer::setTsys(const Vector<Float>& tsys,
			  uInt whichBeam, uInt whichIF) 
//
// Tsys does not depend upon channel but is replicated
// for simplicity of use
//
{

// Get slice and check dim

  IPosition start, end;
  setSlice (start, end, tsys.shape(), tsys_.shape(),
            whichBeam, whichIF, True);

// Get a reference to the Pol/Chan slice we are interested in

  Array<Float> subArr = tsys_(start,end);

// Iterate through it and fill

  VectorIterator<Float> outIt(subArr,asap::ChanAxis);
  uInt i=0;
  while (!outIt.pastEnd()) {
     outIt.vector() = tsys(i++);
     outIt.next();
  }
}

Array<Float> SDContainer::getSpectrum(uInt whichBeam, uInt whichIF) const 
{
  Matrix<Float> spectra(nChan_, nPol_);

  // Beam.
  ArrayAccessor<Float, Axis<asap::BeamAxis> > i0(spectrum_);
  i0.reset(i0.begin(whichBeam));

  // IF.
  ArrayAccessor<Float, Axis<asap::IFAxis> > i1(i0);
  i1.reset(i1.begin(whichIF));

  // Polarization.
  ArrayAccessor<Float, Axis<asap::PolAxis> > i2(i1);
  ArrayAccessor<Float, Axis<asap::IFAxis> > o1(spectra);

  while (i2 != i2.end()) {
    // Channel.
    ArrayAccessor<Float, Axis<asap::ChanAxis> > i3(i2);
    ArrayAccessor<Float, Axis<asap::BeamAxis> > o0(o1);

    while (i3 != i3.end()) {
      *o0 = *i3;

      i3++;
      o0++;
    }

    i2++;
    o1++;
  }

  return spectra.copy();
}


Array<uChar> SDContainer::getFlags(uInt whichBeam, uInt whichIF) const
{
  Matrix<uChar> flagtra(nChan_, nPol_);

  // Beam.
  ArrayAccessor<uChar, Axis<asap::BeamAxis> > i0(flags_);
  i0.reset(i0.begin(whichBeam));

  // IF.
  ArrayAccessor<uChar, Axis<asap::IFAxis> > i1(i0);
  i1.reset(i1.begin(whichIF));

  // Polarization.
  ArrayAccessor<uChar, Axis<asap::PolAxis> > i2(i1);
  ArrayAccessor<uChar, Axis<asap::IFAxis> > o1(flagtra);

  while (i2 != i2.end()) {
    // Channel.
    ArrayAccessor<uChar, Axis<asap::ChanAxis> > i3(i2);
    ArrayAccessor<uChar, Axis<asap::BeamAxis> > o0(o1);

    while (i3 != i3.end()) {
      *o0 = *i3;

      i3++;
      o0++;
    }

    i2++;
    o1++;
  }

  return flagtra.copy();
}

Array<Float> SDContainer::getTsys(uInt whichBeam, uInt whichIF) const
{
  Vector<Float> tsys(nPol_);

  // Beam.
  ArrayAccessor<Float, Axis<asap::BeamAxis> > i0(tsys_);
  i0.reset(i0.begin(whichBeam));

  // IF.
  ArrayAccessor<Float, Axis<asap::IFAxis> > i1(i0);
  i1.reset(i1.begin(whichIF));

  // Channel.
  ArrayAccessor<Float, Axis<asap::ChanAxis> > i3(i1);

  // Polarization.
  ArrayAccessor<Float, Axis<asap::PolAxis> > i2(i3);
  ArrayAccessor<Float, Axis<asap::BeamAxis> > o0(tsys);

  while (i2 != i2.end()) {
    *o0 = *i2;

    i2++;
    o0++;
  }
  return tsys.copy();
}

Array<Double> SDContainer::getDirection(uInt whichBeam) const {
  Vector<Double> direct(2);
  ArrayAccessor<Double, Axis<asap::BeamAxis> > i0(direction_);
  i0.reset(i0.begin(whichBeam));
  ArrayAccessor<Double, Axis<asap::BeamAxis> > o0(direct);
  ArrayAccessor<Double, Axis<asap::IFAxis> > i1(i0);
  while (i1 != i1.end()) {
    *o0 = *i1;
    i1++;
    o0++;
  }  
  return direct.copy();
}


Bool SDContainer::setFrequencyMap(uInt freqID, uInt whichIF) {
  freqidx_[whichIF] = freqID;
  return True;
}

Bool SDContainer::putFreqMap(const Vector<uInt>& freqs) {
  freqidx_.resize();
  freqidx_ = freqs;
  return True;
}

Bool SDContainer::setRestFrequencyMap(uInt freqID, uInt whichIF) {
  restfreqidx_[whichIF] = freqID;
  return True;
}

Bool SDContainer::putRestFreqMap(const Vector<uInt>& freqs) {
  restfreqidx_.resize();
  restfreqidx_ = freqs;
  return True;
}

Bool SDContainer::setDirection(const Vector<Double>& point, uInt whichBeam) {
  if (point.nelements() != 2) return False;
  ArrayAccessor<Double, Axis<asap::BeamAxis> > aa0(direction_);
  aa0.reset(aa0.begin(whichBeam));
  ArrayAccessor<Double, Axis<asap::BeamAxis> > jj(point);
  for (ArrayAccessor<Double, Axis<asap::IFAxis> > i(aa0);i != i.end(); ++i) {
    
    (*i) = (*jj);
    jj++;
  }
  return True;
}

Bool SDContainer::putDirection(const Array<Double>& dir) {
  direction_.resize();
  direction_ = dir;
  return True;
}

Bool SDContainer::appendHistory(const String& hist)
{
  history_.resize(history_.nelements()+1,True);
  history_[history_.nelements()-1] = hist;
  return True;
}
Bool SDContainer::putHistory(const Vector<String>& hist)
{
  history_.resize();
  history_ = hist;
  return True;
}

void SDContainer::setSlice (IPosition& start, IPosition& end,
                            const IPosition& shpIn, const IPosition& shpOut,
                            uInt whichBeam, uInt whichIF, Bool tSys) const
//
// tSYs
//   shpIn [nPol]
// else
//   shpIn [nCHan,nPol]
//
{
  AlwaysAssert(asap::nAxes==4,AipsError);
  if (tSys) {
     AlwaysAssert(shpOut(asap::PolAxis)==shpIn(0),AipsError);     // pol
  } else {
     AlwaysAssert(shpOut(asap::ChanAxis)==shpIn(0),AipsError);    // chan
     AlwaysAssert(shpOut(asap::PolAxis)==shpIn(1),AipsError);     // pol
  }
//
  start.resize(asap::nAxes);
  start = 0;
  start(asap::BeamAxis) = whichBeam;
  start(asap::IFAxis) = whichIF;
//
  end.resize(asap::nAxes);
  end = shpOut-1;
  end(asap::BeamAxis) = whichBeam;
  end(asap::IFAxis) = whichIF;
}


uInt SDFrequencyTable::addFrequency(Double refPix, Double refVal, Double inc) 
{
  if (length() > 0) {
    for (uInt i=0; i< length();i++) {
      if (near(refVal,refVal_[i]) && 
          near(refPix,refPix_[i]) && 
          near(inc,increment_[i])) {
         return i;
      }
    }
  }

// Not found - add it

  nFreq_ += 1;
  refPix_.resize(nFreq_,True);
  refVal_.resize(nFreq_,True);
  increment_.resize(nFreq_,True);
  refPix_[nFreq_-1] = refPix;
  refVal_[nFreq_-1] = refVal;
  increment_[nFreq_-1] = inc;
  return nFreq_-1;
}

uInt SDFrequencyTable::addRestFrequency(Double val)
{
  uInt nFreq = restFreqs_.nelements();
  if (nFreq>0) {
    for (uInt i=0; i<nFreq;i++) {
      if (near(restFreqs_[i],val)) {
         return i;
      }
    }
  }

// Not found - add it

  nFreq += 1;
  restFreqs_.resize(nFreq,True);
  restFreqs_[nFreq-1] = val;
  return nFreq-1;
}


void SDFrequencyTable::restFrequencies(Vector<Double>& rfs, 
				       String& rfunit ) const
{
  rfs.resize(restFreqs_.nelements());
  rfs = restFreqs_;
  rfunit = restFreqUnit_;
}



// SDDataDesc

uInt SDDataDesc::addEntry (const String& source, uInt freqID, const MDirection& dir)
{

// See if already exists

  if (n_ > 0) {
    for (uInt i=0; i<n_; i++) {
      if (source==source_[i] && freqID==freqID_[i]) {
         return i;
      }
    }
  }

// Not found - add it

  n_ += 1;
  source_.resize(n_,True);
  freqID_.resize(n_,True);
  dir_.resize(n_,True,True);
//
  source_[n_-1] = source;
  freqID_[n_-1] = freqID;
  dir_[n_-1] = dir;
//
  return n_-1;
}


void SDDataDesc::summary() const
{
   cerr << "Source    FreqID" << endl;
   for (uInt i=0; i<n_; i++) {
      cerr << setw(11) << source_(i) << freqID_(i) << endl;
   }
}


