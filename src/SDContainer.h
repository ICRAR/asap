//#---------------------------------------------------------------------------
//# SDContainer.h: A container class for single dish integrations
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
#ifndef _SDCONTAINER_H
#define _SDCONTAINER_H

#include <vector>

#include <casa/aips.h>
#include <casa/BasicSL/String.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>

template<class T> class Matrix;

namespace atnf_sd {


struct SDHeader {
  Int nchan;
  Int npol;
  Int nif;
  Int nbeam;
  String observer;
  String project;
  String obstype;
  String antennaname;
  Vector<Double> antennaposition;
  Float equinox;
  String freqref;
  Double reffreq;
  Double bandwidth;
  Double utc;
  void print() const ;
};

class SDFrequencyTable {

public:

  SDFrequencyTable() : nFreq_(0) {;}
  virtual ~SDFrequencyTable() {;}
  
  Int length() const { return nFreq_;};// # of stored Frequencies

  Double referencePixel(uInt which) const { return refPix_[which];}
  Double referenceValue(uInt which) const { return refVal_[which];}
  Double increment(uInt which) const { return increment_[which];}
  Float equinox() const { return equinox_; }
  String refFrame() const { return refFrame_; }

  // returns the index into the table
  // this creates a new one or returns an existing one
  Int addFrequency(Int refPix, Double refVal, Double inc);
  void setEquinox(Float eq) { equinox_ = eq; }
  void setRefFrame(const String& reff) { refFrame_ = reff; }
  
private:
  Int nFreq_;
  Vector<Double> refPix_;
  Vector<Double> refVal_;
  Vector<Double> increment_;
  Float equinox_;
  String refFrame_;
};


class SDContainer {

public:
  SDContainer(uInt nBeam, uInt nIF, uInt nPol, uInt nChan);
  SDContainer(IPosition shp);

  virtual ~SDContainer();

  Bool resize(IPosition shp);

  Bool setSpectrum(const Matrix<Float>& spec,
		   uInt whichBeam, uInt whichIF);
  Bool putSpectrum(const Array<Float>& spec);

  Bool setFlags(const Matrix<uChar>& flgs,
		uInt whichBeam, uInt whichIF);
  Bool putFlags(const Array<uChar>& spec);

  Bool setTsys(const Vector<Float>& ts, 
	       uInt whichBeam, uInt whichIF);
  Bool putTsys(const Array<Float>& spec);

  Bool setDirection(const Vector<Double>& point, uInt whichBeam);
  Bool putDirection(const Array<Double>& dir);

  Bool setFrequencyMap(uInt freqslot, uInt whichIF);
  Bool putFreqMap(const Vector<uInt>& freqs);
  
  Array<Float> getSpectrum(uInt whichBeam, uInt whichIF) const;
  Array<uChar> getFlags(uInt whichBeam, uInt whichIF) const;
  Array<Float> getTsys(uInt whichBeam, uInt whichIF) const;
  Array<Double> getDirection(uInt whichBeam) const;

  const Array<Float>& getSpectrum() const { return spectrum_; }
  const Array<uChar>& getFlags() const { return flags_; }
  const Array<Float>& getTsys() const { return tsys_; }
  const Array<Double>& getDirection() const { return direction_; }

  const Vector<uInt>& getFreqMap() const { return freqidx_; }
  
  Double timestamp;
  String sourcename;
  String fieldname;
  Double interval;
  Int scanid;
  String tcaltime;
  
private:
  uInt nBeam_,nIF_,nPol_,nChan_;

  // (nBeam,nIF,nPol,nChannel)
  Array<Float>    spectrum_;  
  Array<uChar>    flags_;
  // (nBeam,nIF,nPol,[nChannel]) Tsys is not really a function of
  // channel, but this makes it easier to work with at the expense of
  // a little memory
  Array<Float>    tsys_;
  Array<Float>    tcal_;

  //(nIF) indx into "global" frequency table
  Vector<uInt>            freqidx_;
  //(nBeam,2) maybe use Measures here...
  Array<Double>  direction_;

};

} // namespace
#endif
