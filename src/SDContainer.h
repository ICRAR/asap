//#---------------------------------------------------------------------------
//# SDContainer.h: A container class for single dish integrations
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
#ifndef SDCONTAINER_H
#define SDCONTAINER_H

#include <vector>

#include <casa/aips.h>
#include <casa/BasicSL/String.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>

template<class T> class casa::Matrix;

namespace asap {


struct SDHeader {
  casa::Int nchan;
  casa::Int npol;
  casa::Int nif;
  casa::Int nbeam;
  casa::String observer;
  casa::String project;
  casa::String obstype;
  casa::String antennaname;
  casa::Vector<casa::Double> antennaposition;
  casa::Float equinox;
  casa::String freqref;
  casa::Double reffreq;
  casa::Double bandwidth;
  casa::Double utc;
  casa::String fluxunit;
  casa::String epoch;
  void print() const ;
};

class SDFrequencyTable {

public:

  SDFrequencyTable() : nFreq_(0) {;}
  virtual ~SDFrequencyTable() {;}
  
  casa::Int length() const { return nFreq_;};// # of stored Frequencies

  casa::Double referencePixel(casa::uInt which) const { return refPix_[which];}
  casa::Double referenceValue(casa::uInt which) const { return refVal_[which];}
  casa::Double increment(casa::uInt which) const { return increment_[which];}
  casa::Float equinox() const { return equinox_; }
  casa::String refFrame() const { return refFrame_; }
  void restFrequencies(casa::Vector<casa::Double>& rfs, 
		       casa::String& rfunit ) const ;

  // returns the index into the table
  // this creates a new one or returns an existing one
  casa::Int addFrequency(casa::Int refPix, casa::Double refVal, 
			 casa::Double inc);
  void setEquinox(casa::Float eq) { equinox_ = eq; }
  void setRefFrame(const casa::String& reff) { refFrame_ = reff; }
  void addRestFrequency(casa::Double);
  void setRestFrequencyUnit(const casa::String& theunit) {
    restFreqUnit_ = theunit;
  };

private:
  casa::Int nFreq_;
  casa::Vector<casa::Double> refPix_;
  casa::Vector<casa::Double> refVal_;
  casa::Vector<casa::Double> increment_;
  casa::Float equinox_;
  casa::String refFrame_;
  casa::Vector<casa::Double> restFreqs_;
  casa::String restFreqUnit_;
};


class SDContainer {

public:
  SDContainer(casa::uInt nBeam, casa::uInt nIF, casa::uInt nPol, 
	      casa::uInt nChan);
  SDContainer(casa::IPosition shp);

  virtual ~SDContainer();

  casa::Bool resize(casa::IPosition shp);

  casa::Bool setSpectrum(const casa::Matrix<casa::Float>& spec,
		   casa::uInt whichBeam, casa::uInt whichIF);
  casa::Bool putSpectrum(const casa::Array<casa::Float>& spec);

  casa::Bool setFlags(const casa::Matrix<casa::uChar>& flgs,
		      casa::uInt whichBeam, casa::uInt whichIF);
  casa::Bool putFlags(const casa::Array<casa::uChar>& spec);

  casa::Bool setTsys(const casa::Vector<casa::Float>& ts, 
	       casa::uInt whichBeam, casa::uInt whichIF);
  casa::Bool putTsys(const casa::Array<casa::Float>& spec);

  casa::Bool setDirection(const casa::Vector<casa::Double>& point, 
			  casa::uInt whichBeam);
  casa::Bool putDirection(const casa::Array<casa::Double>& dir);

  casa::Bool setFrequencyMap(casa::uInt freqslot, casa::uInt whichIF);
  casa::Bool putFreqMap(const casa::Vector<casa::uInt>& freqs);
  
  casa::Array<casa::Float> getSpectrum(casa::uInt whichBeam, 
				       casa::uInt whichIF) const;
  casa::Array<casa::uChar> getFlags(casa::uInt whichBeam, 
				    casa::uInt whichIF) const;
  casa::Array<casa::Float> getTsys(casa::uInt whichBeam, 
				   casa::uInt whichIF) const;
  casa::Array<casa::Double> getDirection(casa::uInt whichBeam) const;

  const casa::Array<casa::Float>& getSpectrum() const { return spectrum_; }
  const casa::Array<casa::uChar>& getFlags() const { return flags_; }
  const casa::Array<casa::Float>& getTsys() const { return tsys_; }
  const casa::Array<casa::Double>& getDirection() const { return direction_; }

  const casa::Vector<casa::uInt>& getFreqMap() const { return freqidx_; }
  
  const casa::Vector<casa::String>& getHistory() const { return history_; }
  casa::Bool putHistory(const casa::Vector<casa::String>& hist);
  casa::Bool appendHistory(const casa::String& hist);

  casa::Double timestamp;
  //Double bandwidth;
  casa::String sourcename;
  casa::String fieldname;
  casa::Double interval;
  casa::Int scanid;
  casa::Vector<casa::Float> tcal;
  casa::String tcaltime;
  casa::Float azimuth;
  casa::Float elevation;
  casa::Float parangle;
  casa::Int refbeam;

private:
  casa::uInt nBeam_,nIF_,nPol_,nChan_;

  // (nBeam,nIF,nPol,nChannel)
  casa::Array<casa::Float>    spectrum_;  
  casa::Array<casa::uChar>    flags_;
  // (nBeam,nIF,nPol,[nChannel]) Tsys is not really a function of
  // channel, but this makes it easier to work with at the expense of
  // a little memory
  casa::Array<casa::Float>    tsys_;
  casa::Array<casa::Float>    tcal_;

  //(nIF) indx into "global" frequency table
  casa::Vector<casa::uInt>    freqidx_;
  //(nBeam,2) maybe use Measures here...
  casa::Array<casa::Double>   direction_;
  casa::Vector<casa::String> history_;

};

} // namespace
#endif
