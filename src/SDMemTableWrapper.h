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
#ifndef _SDMEMTABLEWRAPPER_H
#define _SDMEMTABLEWRAPPER_H

#include <vector>
#include <string>

#include "SDMemTable.h"

namespace asap {

class SDMemTableWrapper {

public:
  SDMemTableWrapper(const std::string& name) : 
    table_(new SDMemTable(name)) {;}
  SDMemTableWrapper() : 
    table_(new SDMemTable()) {;}
  
  SDMemTableWrapper(CountedPtr<SDMemTable> cp) : table_(cp) {;}
  SDMemTableWrapper(SDMemTable* sdmt) : table_(sdmt) {;}
  
  SDMemTableWrapper(const SDMemTableWrapper& mt, const std::string& expr) :
    table_(new SDMemTable(mt.getCP()->table(), expr)) {;}
  
  SDMemTableWrapper getScan(int scan) {
    String cond("SELECT * from $1 WHERE SCANID == ");
    cond += String::toString(scan);
    return SDMemTableWrapper(*this, cond);
  }

  SDMemTableWrapper getSource(const std::string& source) {
    String cond("SELECT * from $1 WHERE SRCNAME == '");
    cond += source;cond += "'";
    return SDMemTableWrapper(*this, cond);
  }

  std::vector<float> getSpectrum(int whichRow=0) const {
    return table_->getSpectrum(whichRow);
  }

  std::vector<double> getAbscissa(int whichRow, 
				  const std::string& unit = "GHz",
				  const std::string& frame = "TOPO",
				  double restfreq=0.0) const {
    return table_->getAbscissa(whichRow,unit,frame,restfreq);
  }

  float getTsys(int whichRow=0) {return table_->getTsys(whichRow);}
  std::string getTime(int whichRow=0) {return table_->getTime(whichRow);}

  std::vector<bool> getMask(int whichRow=0) const { 
    return table_->getMask(whichRow); 
  }
  bool setMask(const std::vector<int> mvals) const { 
    return table_->setMask(mvals); 
 }

  std::string getSourceName(int whichRow=0) {
    return table_->getSourceName(whichRow);
  }

  bool setIF(int whichIF=0) {return table_->setIF(whichIF);}
  bool setBeam(int whichBeam=0) {return table_->setBeam(whichBeam);}
  bool setPol(int whichPol=0) {return table_->setPol(whichPol);} 

  int getIF() {return table_->getIF();}
  int getBeam() {return table_->getBeam();}
  int getPol() {return table_->getPol();} 

  int nIF() {return table_->nIF();}
  int nBeam() {return table_->nBeam();}
  int nPol() {return table_->nPol();} 
  int nChan() {return table_->nChan();} 
  int nScans() {return table_->nScans();} 

  //sets the mask
  bool setChannels(const std::vector<int>& whichChans) {
    return setChannels(whichChans); 
  }
  void makePersistent(const std::string& fname) {
    table_->makePersistent(fname);
  }
  CountedPtr<SDMemTable> getCP() const {return table_;}
  void summary() { table_->summary(); }
  
private:
  CountedPtr<SDMemTable> table_;
};

} // namespace
#endif
