//#---------------------------------------------------------------------------
//# SDContainer.cc: A container class for single dish integrations
//#---------------------------------------------------------------------------
//# Copyright (C) 2004
//# Malte Marquarding, ATNF
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
#include <casa/Exceptions.h>
#include <tables/Tables/Table.h>
#include <casa/Arrays/IPosition.h>
#include <casa/Arrays/ArrayAccessor.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Quanta/MVTime.h>

#include "SDContainer.h"

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
    freqidx_(shp(1)) {
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
			      uInt whichBeam, uInt whichIF) {

  ArrayAccessor<Float, Axis<0> > aa0(spectrum_);
  aa0.reset(aa0.begin(whichBeam));
  ArrayAccessor<Float, Axis<1> > aa1(aa0);
  aa1.reset(aa1.begin(whichIF));
  
  //Vector<Float> pols(nPol);
  ArrayAccessor<Float, Axis<1> > j(spec);
  IPosition shp0 = spectrum_.shape();
  IPosition shp1 = spec.shape();
  if ( (shp0(2) != shp1(1)) || (shp0(3) != shp1(0)) ) {
    throw(AipsError("Arrays not conformant"));
    return False;
  }
  // assert dimensions are the same....
  for (ArrayAccessor<Float, Axis<2> > i(aa1);i != i.end(); ++i) {
    ArrayAccessor<Float, Axis<0> > jj(j);
    for (ArrayAccessor<Float, Axis<3> > ii(i);ii != ii.end(); ++ii) {
      (*ii) = (*jj);
      jj++;
    }
    j++;
  }
  // unset flags for this spectrum, they might be set again by the
  // setFlags method

  IPosition shp = flags_.shape();
  IPosition start(4,whichBeam,whichIF,0,0);
  IPosition end(4,whichBeam,whichIF,shp(2)-1,shp(3)-1);
  Array<uChar> arr(flags_(start,end));
  arr = uChar(0);
}

Bool SDContainer::setFlags(const Matrix<uChar>& flag,
			   uInt whichBeam, uInt whichIF) {

  ArrayAccessor<uChar, Axis<0> > aa0(flags_);
  aa0.reset(aa0.begin(whichBeam));
  ArrayAccessor<uChar, Axis<1> > aa1(aa0);
  aa1.reset(aa1.begin(whichIF));
  
  ArrayAccessor<uChar, Axis<1> > j(flag);
  IPosition shp0 = flags_.shape();
  IPosition shp1 = flag.shape();
  if ( (shp0(2) != shp1(1)) || (shp0(3) != shp1(0)) ) {
    cerr << "Arrays not conformant" << endl;      
    return False;
  }

  // assert dimensions are the same....
  for (ArrayAccessor<uChar, Axis<2> > i(aa1);i != i.end(); ++i) {
    ArrayAccessor<uChar, Axis<0> > jj(j);
    for (ArrayAccessor<uChar, Axis<3> > ii(i);ii != ii.end(); ++ii) {
      (*ii) = (*jj);
      jj++;
    }
    j++;
  }
  return True;
}

Bool SDContainer::setTsys(const Vector<Float>& tsys,
			  uInt whichBeam, uInt whichIF) {
  ArrayAccessor<Float, Axis<0> > aa0(tsys_);
  aa0.reset(aa0.begin(whichBeam));
  ArrayAccessor<Float, Axis<1> > aa1(aa0);
  aa1.reset(aa1.begin(whichIF));
  // assert dimensions are the same....
  for (ArrayAccessor<Float, Axis<3> > i(aa1);i != i.end(); ++i) {    
    ArrayAccessor<Float, Axis<0> > j(tsys);
    for (ArrayAccessor<Float, Axis<2> > ii(i);ii != ii.end(); ++ii) {
      (*ii) = (*j);
      j++;
    }
  }
}

Array<Float> SDContainer::getSpectrum(uInt whichBeam, uInt whichIF) const {
  Matrix<Float> spectra(nChan_, nPol_);

  // Beam.
  ArrayAccessor<Float, Axis<0> > i0(spectrum_);
  i0.reset(i0.begin(whichBeam));

  // IF.
  ArrayAccessor<Float, Axis<1> > i1(i0);
  i1.reset(i1.begin(whichIF));

  // Polarization.
  ArrayAccessor<Float, Axis<2> > i2(i1);
  ArrayAccessor<Float, Axis<1> > o1(spectra);

  while (i2 != i2.end()) {
    // Channel.
    ArrayAccessor<Float, Axis<3> > i3(i2);
    ArrayAccessor<Float, Axis<0> > o0(o1);

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
  ArrayAccessor<uChar, Axis<0> > i0(flags_);
  i0.reset(i0.begin(whichBeam));

  // IF.
  ArrayAccessor<uChar, Axis<1> > i1(i0);
  i1.reset(i1.begin(whichIF));

  // Polarization.
  ArrayAccessor<uChar, Axis<2> > i2(i1);
  ArrayAccessor<uChar, Axis<1> > o1(flagtra);

  while (i2 != i2.end()) {
    // Channel.
    ArrayAccessor<uChar, Axis<3> > i3(i2);
    ArrayAccessor<uChar, Axis<0> > o0(o1);

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
  ArrayAccessor<Float, Axis<0> > i0(tsys_);
  i0.reset(i0.begin(whichBeam));

  // IF.
  ArrayAccessor<Float, Axis<1> > i1(i0);
  i1.reset(i1.begin(whichIF));

  // Channel.
  ArrayAccessor<Float, Axis<3> > i3(i1);

  // Polarization.
  ArrayAccessor<Float, Axis<2> > i2(i3);
  ArrayAccessor<Float, Axis<0> > o0(tsys);

  while (i2 != i2.end()) {
    *o0 = *i2;

    i2++;
    o0++;
  }
  return tsys.copy();
}

Array<Double> SDContainer::getDirection(uInt whichBeam) const {
  Vector<Double> direct(2);
  ArrayAccessor<Double, Axis<0> > i0(direction_);
  i0.reset(i0.begin(whichBeam));
  ArrayAccessor<Double, Axis<0> > o0(direct);
  ArrayAccessor<Double, Axis<1> > i1(i0);
  while (i1 != i1.end()) {
    *o0 = *i1;
    i1++;
    o0++;
  }  
  return direct.copy();
}


Bool SDContainer::setFrequencyMap(uInt freqslot, uInt whichIF) {
  freqidx_[whichIF] = freqslot;
  return True;
}

Bool SDContainer::putFreqMap(const Vector<uInt>& freqs) {
  freqidx_.resize();
  freqidx_ = freqs;
  return True;
}

Bool SDContainer::setDirection(const Vector<Double>& point, uInt whichBeam) {
  if (point.nelements() != 2) return False;
  ArrayAccessor<Double, Axis<0> > aa0(direction_);
  aa0.reset(aa0.begin(whichBeam));
  ArrayAccessor<Double, Axis<0> > jj(point);
  for (ArrayAccessor<Double, Axis<1> > i(aa0);i != i.end(); ++i) {
    
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

Int SDFrequencyTable::addFrequency(Int refPix, Double refVal, Double inc) {
  Int idx = -1;
  Bool addit = False;
  if (length() > 0) {
    for (uInt i=0; i< length();++i) {
      if ( refVal == refVal_[i] ) { // probably check with tolerance
	if ( refPix == refPix_[i] )
	  if ( inc == increment_[i] )
	    idx = Int(i);
      }
    }
    if (idx >= 0) {
      return idx;
    }
  }
  nFreq_ += 1;
  refPix_.resize(nFreq_,True);
  refVal_.resize(nFreq_,True);
  increment_.resize(nFreq_,True);
  refPix_[nFreq_-1] = refPix;
  refVal_[nFreq_-1] = refVal;
  increment_[nFreq_-1] = inc;
  idx = nFreq_-1;
  return idx;
}

