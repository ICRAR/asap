//#---------------------------------------------------------------------------
//# SDMemTable.cc: A MemoryTable container for single dish integrations
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

#include <aips/iostream.h>
#include <aips/Arrays/Array.h>
#include <aips/Arrays/ArrayMath.h>
#include <aips/Arrays/MaskArrMath.h>
#include <aips/Arrays/ArrayLogical.h>
#include <aips/Arrays/ArrayAccessor.h>

#include <aips/Tables/TableParse.h>
#include <aips/Tables/TableDesc.h>
#include <aips/Tables/SetupNewTab.h>
#include <aips/Tables/ScaColDesc.h>
#include <aips/Tables/ArrColDesc.h>

#include <aips/Tables/ExprNode.h>
#include <aips/Tables/ScalarColumn.h>
#include <aips/Tables/ArrayColumn.h>

#include "SDMemTable.h"
#include "SDContainer.h"

using namespace atnf_sd;

SDMemTable::SDMemTable(const std::string& name) :
  IFSel_(0),
  beamSel_(0),
  polSel_(0) {
  name_ = String(name);
  setup();
}

SDMemTable::SDMemTable(const SDMemTable& other) {
  this->IFSel_= other.IFSel_;
  this->beamSel_= other.beamSel_;
  this->polSel_= other.polSel_;
  this->chanMask_ = other.chanMask_;
  this->name_ = String("dummy");
  this->table_ = other.table_.copyToMemoryTable(String("dummy"));
  // clear all rows()
  this->table_.removeRow(this->table_.rowNumbers());
}

SDMemTable::SDMemTable(const Table& tab, Int scanID) :
  IFSel_(0),
  beamSel_(0),
  polSel_(0) {
  name_ = String("SDMemTable");
  String exprs = String("select * from $1 where SCANID == ")
    +String::toString(scanID);
  cerr << exprs << endl;
  Table t = tableCommand(exprs,tab);
  table_ = t.copyToMemoryTable(name_);
}

SDMemTable::~SDMemTable(){
  cerr << "goodbye from SDMemTable @ " << this << endl;
}

SDMemTable SDMemTable::getScan(Int scanID) {
  return SDMemTable(table_, scanID);
}

void SDMemTable::setup() {
  TableDesc td("", "1", TableDesc::Scratch);
  td.comment() = "A SDMemTable";
  //td.rwKeywordSet().define("VERSION",Float(0.1));
  td.addColumn(ScalarColumnDesc<Double>("TIME"));
  td.addColumn(ScalarColumnDesc<String>("SRCNAME"));
  td.addColumn(ArrayColumnDesc<Float>("SPECTRA"));
  td.addColumn(ArrayColumnDesc<uChar>("FLAGTRA"));
  td.addColumn(ArrayColumnDesc<Float>("TSYS"));  
  td.addColumn(ScalarColumnDesc<Int>("SCANID"));  
  // Now create a new table from the description.
  SetupNewTable aNewTab(name_, td, Table::New);
  table_ = Table(aNewTab, Table::Memory, 0);  
}

std::string SDMemTable::name() const {
  return name_;
}

std::string SDMemTable::getSourceName(Int whichRow) const {
  ROScalarColumn<String> src(table_, "SRCNAME");
  String name;
  src.get(whichRow, name);
  return name;
}

Double SDMemTable::getTime(Int whichRow) const {
  ROScalarColumn<Double> src(table_, "TIME");
  Double tm;
  src.get(whichRow, tm);
  return tm;
}

bool SDMemTable::setIF(Int whichIF) {
  //if ( whichIF >= 0 && whichIF < nIF_) {
    IFSel_ = whichIF;
    return true;
    //}
    //return false;
}
bool SDMemTable::setBeam(Int whichBeam) {
  //if ( whichBeam >= 0 && whichBeam < nBeam_) {
    beamSel_ = whichBeam;
    return true;
    //}
    //return false;

}
bool SDMemTable::setPol(Int whichPol) {
  //if ( whichPol >= 0 && whichPol < nPol_) {
    polSel_ = whichPol;
    return true;
    //}
    //return false;
}

bool SDMemTable::setChannels(const std::vector<int>& whichChans) {
  //std::vector<bool>::iterator it;
  //std::fill(chanMask_.begin(), chanMask_.end(), false);
  //for (it = whichChans.begin(); it != whichChans.end(); it++) {
  //chanMask_[*it] = true;
  //}
  cout << "SDMemTable::setChannels() disabled" << endl;
  return true;
}

std::vector<bool> SDMemTable::getMask() const {
  return chanMask_;
}
std::vector<float> SDMemTable::getSpectrum(Int whichRow) const {

  std::vector<float> spectrum;
  ROArrayColumn<Float> spec(table_, "SPECTRA");
  Array<Float> arr;
  spec.get(whichRow, arr);
  ArrayAccessor<Float, Axis<0> > aa0(arr);
  aa0.reset(aa0.begin(uInt(beamSel_)));//go to beam
  ArrayAccessor<Float, Axis<1> > aa1(aa0);
  aa1.reset(aa1.begin(uInt(IFSel_)));// go to IF
  for (ArrayAccessor<Float, Axis<2> > i(aa1); i != i.end(); ++i) {
    ArrayAccessor<Float, Axis<3> > aa2(i);
    aa2.reset(aa2.begin(uInt(polSel_)));// go to pol
    spectrum.push_back(*aa2);
  }
  return spectrum;
}

MaskedArray<Float> SDMemTable::rowAsMaskedArray(uInt whichRow) {
  ROArrayColumn<Float> spec(table_, "SPECTRA");
  Array<Float> arr;
  ROArrayColumn<uChar> flag(table_, "FLAGTRA");
  Array<uChar> farr;
  spec.get(whichRow, arr);
  flag.get(whichRow, farr);
  IPosition test(4,0,0,100,0);
  Array<Bool> barr(farr.shape());convertArray(barr, farr);
  MaskedArray<Float> marr;
  marr.setData(arr,!barr);
  return marr;
}


Float SDMemTable::getTsys(Int whichRow) const {
  ROArrayColumn<Float> ts(table_, "TSYS");
  Array<Float> arr;
  ts.get(whichRow, arr);
  Float out;
  IPosition ip(arr.shape());
  ip(0) = beamSel_;ip(1) = IFSel_;ip(2)=0;ip(3) = polSel_;
  out = arr(ip);
  return out;
}

bool SDMemTable::putSDContainer(const SDContainer& sdc) {
  ScalarColumn<Double> mjd(table_, "TIME");
  ScalarColumn<String> srcn(table_, "SRCNAME");
  ArrayColumn<Float> spec(table_, "SPECTRA");
  ArrayColumn<uChar> flags(table_, "FLAGTRA");
  ArrayColumn<Float> ts(table_, "TSYS");
  ScalarColumn<Int> scan(table_, "SCANID");

  uInt rno = table_.nrow();
  table_.addRow();
  
  mjd.put(rno, sdc.timestamp);
  srcn.put(rno, sdc.sourcename);
  spec.put(rno, sdc.getSpectrum());
  flags.put(rno, sdc.getFlags());
  ts.put(rno, sdc.getTsys());
  scan.put(rno, sdc.scanid);
  return true;
}
void SDMemTable::makePersistent(const std::string& filename) {
  table_.deepCopy(filename,Table::New);
}

void SDMemTable::summary() const {
  cerr << "SDMemTable::summary()" << endl;
  ROScalarColumn<Int> scans(table_, "SCANID");
  ROScalarColumn<String> srcs(table_, "SRCNAME");
  uInt count = 0;
  String name;
  Int previous = -1;Int current=0;
  cout << "Scan\tSource" << endl;
  for (uInt i=0; i< scans.nrow();i++) {
    scans.getScalar(i,current);
    if (previous != current) {
      srcs.getScalar(i,name);
      previous = current;     
      count++;
      cout << count << "\t" << name << endl;
    }
  }
  cout << "Table contains " << table_.nrow() << "integrations." << endl;
  cout << "in " << count << "scans." << endl;
}

