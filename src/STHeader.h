//#---------------------------------------------------------------------------
//# STHeader.h: A container class for single dish integrations
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
#ifndef STHEADER_H
#define STHEADER_H

#include <vector>

#include <casa/aips.h>
#include <casa/BasicSL/String.h>
#include <casa/Arrays/Vector.h>
#include <casa/Containers/Block.h>
#include <measures/Measures/MDirection.h>

namespace casacore {
  template<class T> class Matrix;
}

namespace asap {


struct STHeader {

  bool conformant(const STHeader& other);
  casacore::String diff( const STHeader& other );


  casacore::Int nchan;
  casacore::Int npol;
  casacore::Int nif;
  casacore::Int nbeam;
  casacore::String observer;
  casacore::String project;
  casacore::String obstype;
  casacore::String antennaname;
  casacore::Vector<casacore::Double> antennaposition;
  casacore::Float equinox;
  casacore::String freqref;
  casacore::Double reffreq;
  casacore::Double bandwidth;
  casacore::Double utc;
  casacore::String fluxunit;
  casacore::String epoch;
  casacore::String poltype;
  void print() const ;
};

class SDDataDesc {

public:

  // Constructor
  SDDataDesc() : n_(0) {;}
  ~SDDataDesc() {;}

  // Add an entry if source name and Integer ID (can be anything you
  // like, such as FreqID) are unique.  You can add secondary entries
  // direction and another integer index which are just stored along
  // with the the primary entries
  casacore::uInt addEntry(const casacore::String& source, casacore::uInt ID,
                      const casacore::MDirection& secDir, casacore::uInt secID);

  // Number of entries
  casacore::Int length() const { return n_;}

  // Get attributes
  casacore::String source(casacore::uInt which) const {return source_[which];}
  casacore::uInt ID(casacore::uInt which) const {return ID_[which];}
  casacore::uInt secID(casacore::uInt which) const {return secID_[which];}
  casacore::MDirection secDir(casacore::uInt which) const {return secDir_[which];}

  // Summary
  void summary() const;

private:
  casacore::uInt n_;
  casacore::Vector<casacore::String> source_;
  casacore::Vector<casacore::uInt> ID_, secID_;
  casacore::Block<casacore::MDirection> secDir_;

  SDDataDesc(const SDDataDesc& other);

};


} // namespace
#endif
