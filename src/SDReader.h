//#---------------------------------------------------------------------------
//# SDReader.h: A class to read single dish spectra from SDFITS, RPFITS
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
#ifndef _SDREADER_H_
#define _SDREADER_H_

#include <vector>
#include <string>

#include <aips/aips.h>
#include <aips/iostream.h>
#include <aips/iomanip.h>
#include <aips/Utilities/CountedPtr.h>
#include <aips/Utilities/String.h>
#include <aips/Arrays/Vector.h>

#include "SDMemTable.h"

class PKSreader;

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

class SDReader {
public:
  SDReader();
  SDReader(CountedPtr<SDMemTable> tbl);
  virtual ~SDReader();

  void open(const std::string& filename);
  int read(const std::vector<int>& seq);

  CountedPtr<SDMemTable> getTable() const { return table_;}

  void reset();

  std::vector<int> pseudoHeader() const {
    std::vector<int> v;
    v.push_back(nBeam_);v.push_back(nIF_);
    v.push_back(nChan_);v.push_back(nPol_);
    return v;
  }

protected:

private:
  Int nBeam_,nIF_,nPol_,nChan_;
  Bool getHeader();
  PKSreader* reader_;  
  SDHeader header_;
  CountedPtr<SDMemTable> table_;
  String filename_;
  uInt cursor_;
  Double timestamp_;
};

}// namespace
#endif

