//#---------------------------------------------------------------------------
//# SDMemTableWrapper.h: Wrapper classes to use CountedPtr
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
#ifndef _SDMEMTABLEWRAPPER_H_
#define _SDMEMTABLEWRAPPER_H_

#include <vector>
#include <string>

#include "SDMemTable.h"
#include "SDReader.h"
#include "SDMath.h"

namespace atnf_sd {

class SDMemTableWrapper {

public:
  SDMemTableWrapper(const std::string& name = "SDinput.tbl") : 
    table_(new SDMemTable(name)) {;}
  SDMemTableWrapper(CountedPtr<SDMemTable> cp) : table_(cp) {;}
  SDMemTableWrapper(SDMemTable* sdmt) : table_(sdmt) {;}
  
  SDMemTableWrapper(const SDMemTableWrapper& mt, int scan) :
    table_(new SDMemTable(mt.getCP()->table(), scan)) {;}
  
  SDMemTableWrapper getScan(int scan) {
    return SDMemTableWrapper(*this, scan);
  }
  std::vector<float> getSpectrum(int whichRow) const {
    return table_->getSpectrum(whichRow);
  }
  float getTsys(int whichRow) {return table_->getTsys(whichRow);}
  double getTime(int whichRow) {return table_->getTime(whichRow);}

  std::vector<bool> getMask(int whichRow) const { 
    return table_->getMask(whichRow); 
  }
  bool setMask(const std::vector<int> mvals) const { 
    return table_->setMask(mvals); 
  }

  std::string getSourceName(int whichRow) {
    return table_->getSourceName(whichRow);
  }

  bool setIF(int whichIF=0) {return table_->setIF(whichIF);}
  bool setBeam(int whichBeam=0) {return table_->setBeam(whichBeam);}
  bool setPol(int whichPol=0) {return table_->setPol(whichPol);} 

  int getIF() {return table_->getIF();}
  int getBeam() {return table_->getBeam();}
  int getPol() {return table_->getPol();} 


  //sets the mask
  bool setChannels(const std::vector<int>& whichChans) {
    return setChannels(whichChans); 
  }
  void makePersistent(const std::string& fname) {
    table_->makePersistent(fname);
  }
  CountedPtr<SDMemTable> getCP() const {return table_;}
  void summary() { table_->summary(); }
  std::string name() { return table_->name(); }
  
private:
  CountedPtr<SDMemTable> table_;
};

class SDReaderWrapper : public SDReader {
public:
  SDReaderWrapper() {;}
  SDReaderWrapper(SDMemTableWrapper tbl) :
    SDReader(tbl.getCP()){;}
  SDMemTableWrapper getSDMemTable() const {
    return SDMemTableWrapper(getTable());
  }
};

class SDMathWrapper {
public:
  SDMemTableWrapper average(const SDMemTableWrapper& sdt) {
    return SDMemTableWrapper(SDMath::average(sdt.getCP()));
  }
  SDMemTableWrapper quotient(const SDMemTableWrapper& on,
			     const SDMemTableWrapper& off) {
    return SDMemTableWrapper(SDMath::quotient(on.getCP(),
					     off.getCP()));
  }
  SDMemTableWrapper multiply(const SDMemTableWrapper& in,
			     Float factor) {
    return SDMemTableWrapper(SDMath::multiply(in.getCP(),factor));
  }
};

} // namespace

#endif

