//#---------------------------------------------------------------------------
//# SDMemTableWrapper.h: Wrapper classes to use CountedPtr
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
#ifndef SDMEMTABLEWRAPPER_H
#define SDMEMTABLEWRAPPER_H

#include <vector>
#include <string>
#include <casa/Arrays/Vector.h>

#include "SDMemTable.h"

namespace asap {

class SDMemTableWrapper {

public:
  SDMemTableWrapper(const std::string& name) :
    table_(new SDMemTable(name)) {;}
  SDMemTableWrapper() :
    table_(new SDMemTable()) {;}

  SDMemTableWrapper(casa::CountedPtr<SDMemTable> cp) : table_(cp) {;}
  //SDMemTableWrapper(SDMemTable* sdmt) : table_(sdmt) {;}
  SDMemTableWrapper(const SDMemTableWrapper& mt) :
    table_(mt.getCP()) {;}

  SDMemTableWrapper(const SDMemTableWrapper& mt, const std::string& expr) :
    table_(new SDMemTable(mt.getCP()->table(), expr)) {;}

  SDMemTableWrapper copy() {
    //CountedPtr<SDMemTable> cp = new SDMemTable(*this, False);
    return SDMemTableWrapper(new SDMemTable(*(this->getCP()), casa::False));
  }

  //SDMemTableWrapper getScan(int scan) {
  SDMemTableWrapper getScan(std::vector<int> scan) {
    casa::String cond("SELECT FROM $1 WHERE SCANID IN ");
    casa::Vector<casa::Int> v(scan);
    casa::ostringstream oss;
    oss << v;
    cond += casa::String(oss);
    return SDMemTableWrapper(*this, cond);
  }

  SDMemTableWrapper getSource(const std::string& source) {
    casa::String cond("SELECT * from $1 WHERE SRCNAME == pattern('");
    cond += source;cond += "')";
    return SDMemTableWrapper(*this, cond);
  }

  std::vector<float> getSpectrum(int whichRow=0) const {
    return table_->getSpectrum(whichRow);
  }

  std::vector<float> getStokesSpectrum(int whichRow=0, bool doPol=false, float paOff=0.0) const {
    return table_->getStokesSpectrum(whichRow, doPol, paOff);
  }

  std::vector<float> getCircularSpectrum(int whichRow=0, bool doRR=true) const {
    return table_->getCircularSpectrum(whichRow, doRR);
  }

  std::vector<double> getAbcissa(int whichRow=0) const {
    return table_->getAbcissa(whichRow);
  }
  std::string getAbcissaString(int whichRow=0) const {
    return table_->getAbcissaString(whichRow);
  }

  std::vector<float> getTsys() {
     int nRow = table_->nRow();
     std::vector<float> result(nRow);
     for (uint i=0; i<nRow; i++) {
        result[i] = table_->getTsys(i);
     }
     return result;
  }

  std::string getTime(int whichRow=0) {return table_->getTime(whichRow);}

  std::string getFluxUnit() const {return table_->getFluxUnit();}
  void setFluxUnit(const std::string& unit) {table_->setFluxUnit(unit);}

  void setInstrument(const std::string& name) {table_->setInstrument(name);}

  std::vector<bool> getMask(int whichRow=0) const {
    return table_->getMask(whichRow);
  }
  bool setMask(const std::vector<int> mvals) const {
    return table_->setMask(mvals);
  }

  void flag(int whichRow=-1) {
    table_->flag(whichRow);
  }
  std::string getSourceName(int whichRow=0) {
    return table_->getSourceName(whichRow);
  }

  void setSpectrum(std::vector<float> spectrum, int whichRow=0) {
      table_->setSpectrum(spectrum, whichRow);
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
  int nScan() {return table_->nScan();}
  int nRow() {return table_->nRow();}

  //sets the mask
  bool setChannels(const std::vector<int>& whichChans) {
    return setChannels(whichChans);
  }
  void makePersistent(const std::string& fname) {
    table_->makePersistent(fname);
  }

  bool setRestFreqs(const std::vector<double>& restfreqs, 
                    const std::string& unit, 
                    const std::vector<std::string>& lines,
                    const std::string& source, int whichIF) {
    casa::Vector<casa::Double> restFreqs2(restfreqs);
    return table_->setRestFreqs(restFreqs2, casa::String(unit),
                                lines, casa::String(source),
                                casa::Int(whichIF));
  }

  void spectralLines() const {table_->spectralLines();}

  std::vector<double> getRestFreqs() {
    return table_->getRestFreqs();
  }

  void setCoordInfo(std::vector<string> theinfo) {
    table_->setCoordInfo(theinfo);
  }
  std::vector<string> getCoordInfo() const {
    return table_->getCoordInfo();
  }

  casa::CountedPtr<SDMemTable> getCP() const {return table_;}
  SDMemTable* getPtr() {return &(*table_);}

  std::string summary(bool verbose=false) const { 
    return table_->summary(verbose); 
  }
  
  std::vector<std::string> history(int whichRow=0) { 
    return table_->history(whichRow); 
  }

  void rotateXYPhase (float value) {
      table_->rotateXYPhase(value);
  }

private:
  casa::CountedPtr<SDMemTable> table_;
};

} // namespace
#endif

