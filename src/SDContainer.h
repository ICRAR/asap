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
#ifndef _SDCONTAINER_H_
#define _SDCONTAINER_H_

#include <vector>

#include <aips/aips.h>
#include <aips/Utilities/String.h>
#include <aips/Arrays/Array.h>
#include <aips/Arrays/Vector.h>

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

  SDFrequencyTable() {;}
  // returns the index into the table
  // this creates a new one or returns an existing one
  Int addFrequency(Double refPix, Double refVal, Double inc) {;}
  
  Int length() const { return nFreq_;};// # of stored Frequencies
  // returns a Table with nRows == nFreq, and three cols
  
private:
  Int nFreq_;
  std::vector<double> refPix_;
  std::vector<double> revVal_;
  std::vector<double> increment_;
};


class SDContainer {

public:
  SDContainer(uInt nBeam, uInt nIF, uInt nPol, uInt nChan);
  SDContainer(IPosition shp);

  virtual ~SDContainer();

  Bool setSpectrum(const Matrix<Float>& spec,
		   uInt whichBeam, uInt whichIF);
  Bool putSpectrum(const Array<Float>& spec);

  Bool setFlags(const Matrix<uChar>& flgs,
		uInt whichBeam, uInt whichIF);
  Bool putFlags(const Array<uChar>& spec);

  Bool setTsys(const Vector<Float>& ts, 
	       uInt whichBeam, uInt whichIF);
  Bool putTsys(const Array<Float>& spec);

  Bool setPointing(const Vector<Double>& point, uInt whichBeam) {;}

  Bool setFrequencyMap(uInt freqslot, uInt whichIF) {;}
  
  const Array<Float>& getSpectrum() const { return spectrum_; }
  const Array<uChar>& getFlags() const { return flags_; }
  const Array<Float>& getTsys() const { return tsys_; }

  Double timestamp;
  String sourcename;
  Double interval;
  uInt scanid;
  
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

  //(nBeam) maybe use Measures here...
  //*** Vector<Vector<Double>>  pointing_;
  //(nIF) indx into "global" frequency table
  Vector<uInt>            freqidx_;

};

} // namespace
#endif
