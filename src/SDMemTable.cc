//#---------------------------------------------------------------------------
//# SDMemTable.cc: A MemoryTable container for single dish integrations
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

#include <casa/aips.h>
#include <casa/iostream.h>
#include <casa/iomanip.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/MaskArrMath.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Arrays/ArrayAccessor.h>

#include <tables/Tables/TableParse.h>
#include <tables/Tables/TableDesc.h>
#include <tables/Tables/SetupNewTab.h>
#include <tables/Tables/ScaColDesc.h>
#include <tables/Tables/ArrColDesc.h>

#include <tables/Tables/ExprNode.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/TableRecord.h>
#include <measures/Measures/MFrequency.h>
#include <measures/Measures/MeasTable.h>
#include <coordinates/Coordinates/CoordinateUtil.h>
#include <casa/Quanta/MVTime.h>

#include "SDMemTable.h"
#include "SDContainer.h"

using namespace casa;
using namespace asap;

SDMemTable::SDMemTable() :
  IFSel_(0),
  beamSel_(0),
  polSel_(0) {
  setup();
}
SDMemTable::SDMemTable(const std::string& name) :
  IFSel_(0),
  beamSel_(0),
  polSel_(0) {
  Table tab(name);
  table_ = tab.copyToMemoryTable("dummy");
}

SDMemTable::SDMemTable(const SDMemTable& other, Bool clear) {
  this->IFSel_= other.IFSel_;
  this->beamSel_= other.beamSel_;
  this->polSel_= other.polSel_;
  this->chanMask_ = other.chanMask_;
  this->table_ = other.table_.copyToMemoryTable(String("dummy"));
  // clear all rows()
  if (clear) {
    this->table_.removeRow(this->table_.rowNumbers());
  } else {
    this->IFSel_ = other.IFSel_;
    this->beamSel_ = other.beamSel_;
    this->polSel_ = other.polSel_;
  }
}

SDMemTable::SDMemTable(const Table& tab, const std::string& exprs) :
  IFSel_(0),
  beamSel_(0),
  polSel_(0) {
  Table t = tableCommand(exprs,tab);
  if (t.nrow() == 0)
      throw(AipsError("Query unsuccessful."));
  table_ = t.copyToMemoryTable("dummy");
}

SDMemTable::~SDMemTable(){
  //cerr << "goodbye from SDMemTable @ " << this << endl;
}

SDMemTable SDMemTable::getScan(Int scanID) {
  String cond("SELECT * from $1 WHERE SCANID == ");
  cond += String::toString(scanID);
  return SDMemTable(table_, cond);
}

SDMemTable SDMemTable::getSource(const std::string& source) {
  String cond("SELECT * from $1 WHERE SRCNAME == ");
  cond += source;
  return SDMemTable(table_, cond);
}

void SDMemTable::setup() {
  TableDesc td("", "1", TableDesc::Scratch);
  td.comment() = "A SDMemTable";
  td.addColumn(ScalarColumnDesc<Double>("TIME"));
  td.addColumn(ScalarColumnDesc<String>("SRCNAME"));
  td.addColumn(ArrayColumnDesc<Float>("SPECTRA"));
  td.addColumn(ArrayColumnDesc<uChar>("FLAGTRA"));
  td.addColumn(ArrayColumnDesc<Float>("TSYS"));
  td.addColumn(ScalarColumnDesc<Int>("SCANID"));
  td.addColumn(ScalarColumnDesc<Double>("INTERVAL"));
  td.addColumn(ArrayColumnDesc<uInt>("FREQID"));
  td.addColumn(ArrayColumnDesc<Double>("DIRECTION"));
  td.addColumn(ScalarColumnDesc<String>("FIELDNAME"));
  td.addColumn(ScalarColumnDesc<String>("TCALTIME"));
  td.addColumn(ArrayColumnDesc<Float>("TCAL"));
  td.addColumn(ScalarColumnDesc<Float>("AZIMUTH"));
  td.addColumn(ScalarColumnDesc<Float>("ELEVATION"));
  td.addColumn(ScalarColumnDesc<Float>("PARANGLE"));
  td.addColumn(ScalarColumnDesc<Int>("REFBEAM"));

  // Now create a new table from the description.

  SetupNewTable aNewTab("dummy", td, Table::New);
  table_ = Table(aNewTab, Table::Memory, 0);
}

std::string SDMemTable::getSourceName(Int whichRow) const {
  ROScalarColumn<String> src(table_, "SRCNAME");
  String name;
  src.get(whichRow, name);
  return name;
}

std::string SDMemTable::getTime(Int whichRow) const {
  ROScalarColumn<Double> src(table_, "TIME");
  Double tm;
  src.get(whichRow, tm);
  MVTime mvt(tm);
  mvt.setFormat(MVTime::TIME);
  ostringstream oss;
  oss << mvt;
  String str(oss);
  return str;
}
double SDMemTable::getInterval(Int whichRow) const {
  ROScalarColumn<Double> src(table_, "INTERVAL");
  Double intval;
  src.get(whichRow, intval);
  return intval;
}

bool SDMemTable::setIF(Int whichIF) {
  if ( whichIF >= 0 && whichIF < nIF()) {
    IFSel_ = whichIF;
    return true;
  }
  return false;
}

bool SDMemTable::setBeam(Int whichBeam) {
  if ( whichBeam >= 0 && whichBeam < nBeam()) {
    beamSel_ = whichBeam;
    return true;
  }
  return false;
}

bool SDMemTable::setPol(Int whichPol) {
  if ( whichPol >= 0 && whichPol < nPol()) {
    polSel_ = whichPol;
    return true;
  }
  return false;
}

bool SDMemTable::setMask(std::vector<int> whichChans) {
  ROArrayColumn<uChar> spec(table_, "FLAGTRA");
  std::vector<int>::iterator it;
  uInt n = spec.shape(0)(3);
  if (whichChans.empty()) {
    chanMask_ = std::vector<bool>(n,true);
    return true;      
  }
  chanMask_.resize(n,true);
  for (it = whichChans.begin(); it != whichChans.end(); ++it) {
    if (*it < n) {
      chanMask_[*it] = false;
    }
  }
  return true;
}

std::vector<bool> SDMemTable::getMask(Int whichRow) const {
  std::vector<bool> mask;
  ROArrayColumn<uChar> spec(table_, "FLAGTRA");
  Array<uChar> arr;
  spec.get(whichRow, arr);
  ArrayAccessor<uChar, Axis<0> > aa0(arr);
  aa0.reset(aa0.begin(uInt(beamSel_)));//go to beam
  ArrayAccessor<uChar, Axis<1> > aa1(aa0);
  aa1.reset(aa1.begin(uInt(IFSel_)));// go to IF
  ArrayAccessor<uChar, Axis<2> > aa2(aa1);
  aa2.reset(aa2.begin(uInt(polSel_)));// go to pol

  Bool useUserMask = ( chanMask_.size() == arr.shape()(3) );

  std::vector<bool> tmp;
  tmp = chanMask_; // WHY the fxxx do I have to make a copy here
  std::vector<bool>::iterator miter;
  miter = tmp.begin();

  for (ArrayAccessor<uChar, Axis<3> > i(aa2); i != i.end(); ++i) {
    bool out =!static_cast<bool>(*i);
    if (useUserMask) {
      out = out && (*miter);
      miter++;
    }
    mask.push_back(out);
  }
  return mask;
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
  ArrayAccessor<Float, Axis<2> > aa2(aa1);
  aa2.reset(aa2.begin(uInt(polSel_)));// go to pol
  for (ArrayAccessor<Float, Axis<3> > i(aa2); i != i.end(); ++i) {
    spectrum.push_back(*i);
  }
  return spectrum;
}
std::vector<string> SDMemTable::getCoordInfo() const {
  String un;
  Table t = table_.keywordSet().asTable("FREQUENCIES");
  String sunit;
  t.keywordSet().get("UNIT",sunit);
  String dpl;
  t.keywordSet().get("DOPPLER",dpl);
  if (dpl == "") dpl = "RADIO";
  String rfrm;
  t.keywordSet().get("REFFRAME",rfrm);
  std::vector<string> inf;
  inf.push_back(sunit);
  inf.push_back(rfrm);
  inf.push_back(dpl);
  return inf;
}

void SDMemTable::setCoordInfo(std::vector<string> theinfo) {

  std::vector<string>::iterator it;
  String un,rfrm,dpl;
  un = theinfo[0];
  rfrm = theinfo[1];
  dpl = theinfo[2];

  //String un(theunit);
  Table t = table_.rwKeywordSet().asTable("FREQUENCIES");
  Vector<Double> rstf;
  t.keywordSet().get("RESTFREQS",rstf);
  Bool canDo = True;
  Unit u1("km/s");Unit u2("Hz");
  if (Unit(un) == u1) {
    Vector<Double> rstf;
    t.keywordSet().get("RESTFREQS",rstf);
    if (rstf.nelements() == 0) {
        throw(AipsError("Can't set unit to km/s if no restfrequencies are specified"));
    }
  } else if (Unit(un) != u2 && un != "") {
        throw(AipsError("Unit not conformant with Spectral Coordinates"));
  }
  t.rwKeywordSet().define("UNIT", un);

  MFrequency::Types mdr;
  if (!MFrequency::getType(mdr, rfrm)) {
    
    Int a,b;const uInt* c;
    const String* valid = MFrequency::allMyTypes(a, b, c);
    String pfix = "Please specify a legal frame type. Types are\n";
    throw(AipsError(pfix+(*valid)));
  } else {
    t.rwKeywordSet().define("REFFRAME",rfrm);
  }

}

std::vector<double> SDMemTable::getAbscissa(Int whichRow) {
  std::vector<double> absc(nChan());
  Vector<Double> absc1(nChan());
  indgen(absc1);
  ROArrayColumn<uInt> fid(table_, "FREQID");
  Vector<uInt> v;
  fid.get(whichRow, v);
  uInt specidx = v(IFSel_);
  SpectralCoordinate spc = getCoordinate(specidx);
  Table t = table_.keywordSet().asTable("FREQUENCIES");
  String rf;
  //t.keywordSet().get("EQUINOX",rf);
  MDirection::Types mdr;
  //if (!MDirection::getType(mdr, rf)) {
  mdr = MDirection::J2000;
  //cout << "Unknown equinox using J2000" << endl;
  //}
  ROArrayColumn<Double> dir(table_, "DIRECTION");
  Array<Double> posit;
  dir.get(whichRow,posit);
  Vector<Double> wpos(2);
  wpos[0] = posit(IPosition(2,beamSel_,0));
  wpos[1] = posit(IPosition(2,beamSel_,1));
  Quantum<Double> lon(wpos[0],Unit(String("rad")));
  Quantum<Double> lat(wpos[1],Unit(String("rad")));
  MDirection direct(lon, lat, mdr);
  ROScalarColumn<Double> tme(table_, "TIME");
  Double obstime;
  tme.get(whichRow,obstime);
  MVEpoch tm2(Quantum<Double>(obstime, Unit(String("d"))));
  MEpoch epoch(tm2);

  Vector<Double> antpos;
  table_.keywordSet().get("AntennaPosition", antpos);
  MVPosition mvpos(antpos(0),antpos(1),antpos(2));
  MPosition pos(mvpos);
  String sunit;
  t.keywordSet().get("UNIT",sunit);
  if (sunit == "") sunit = "pixel";
  Unit u(sunit);
  String frm;
  t.keywordSet().get("REFFRAME",frm);
  if (frm == "") frm = "TOPO";
  String dpl;
  t.keywordSet().get("DOPPLER",dpl);
  if (dpl == "") dpl = "RADIO";
  MFrequency::Types mtype;
  if (!MFrequency::getType(mtype, frm)) {
    cout << "Frequency type unknown assuming TOPO" << endl;
    mtype = MFrequency::TOPO;
  }
  
  if (!spc.setReferenceConversion(mtype,epoch,pos,direct)) {
    throw(AipsError("Couldn't convert frequency frame."));
  }

  if ( u == Unit("km/s") ) {
    Vector<Double> rstf;
    t.keywordSet().get("RESTFREQS",rstf);
    if (rstf.nelements() > 0) {
      if (rstf.nelements() >= nIF())
        spc.selectRestFrequency(uInt(IFSel_));
      spc.setVelocity(u.getName());
      Vector<Double> wrld;
      spc.pixelToVelocity(wrld,absc1);
      std::vector<double>::iterator it;
      uInt i = 0;
      for (it = absc.begin(); it != absc.end(); ++it) {
        (*it) = wrld[i];
        i++;
      }
    }
  } else if (u == Unit("Hz")) {
    Vector<String> wau(1); wau = u.getName();
    spc.setWorldAxisUnits(wau);
    std::vector<double>::iterator it;
    Double tmp;
    uInt i = 0;
    for (it = absc.begin(); it != absc.end(); ++it) {
      spc.toWorld(tmp,absc1[i]);
      (*it) = tmp;
      i++;
    }

  } else {
    // assume channels/pixels
    std::vector<double>::iterator it;
    uInt i=0;
    for (it = absc.begin(); it != absc.end(); ++it) {
      (*it) = Double(i++);
    }
  }
  return absc;
}

std::string SDMemTable::getAbscissaString(Int whichRow)
{
  ROArrayColumn<uInt> fid(table_, "FREQID");
  Table t = table_.keywordSet().asTable("FREQUENCIES");
  String sunit;
  t.keywordSet().get("UNIT",sunit);
  if (sunit == "") sunit = "pixel";
  Unit u(sunit);
  Vector<uInt> v;
  fid.get(whichRow, v);
  uInt specidx = v(IFSel_);
  SpectralCoordinate spc = getCoordinate(specidx);
  String frm;
  t.keywordSet().get("REFFRAME",frm);
  MFrequency::Types mtype;
  if (!MFrequency::getType(mtype, frm)) {
    cout << "Frequency type unknown assuming TOPO" << endl;
    mtype = MFrequency::TOPO;
  }
  spc.setFrequencySystem(mtype);
  String s = "Channel";
  if (u == Unit("km/s")) { 
    spc.setVelocity(u.getName());
    s = CoordinateUtil::axisLabel(spc,0,True,True,True);
  } else if (u == Unit("Hz")) {
    Vector<String> wau(1);wau = u.getName();
    spc.setWorldAxisUnits(wau);
    s = CoordinateUtil::axisLabel(spc);
  }
  return s;
}

void SDMemTable::setSpectrum(std::vector<float> spectrum, int whichRow) {
  ArrayColumn<Float> spec(table_, "SPECTRA");
  Array<Float> arr;
  spec.get(whichRow, arr);
  if (spectrum.size() != arr.shape()(3)) {
    throw(AipsError("Attempting to set spectrum with incorrect length."));
  }

  ArrayAccessor<Float, Axis<0> > aa0(arr);
  aa0.reset(aa0.begin(uInt(beamSel_)));//go to beam
  ArrayAccessor<Float, Axis<1> > aa1(aa0);
  aa1.reset(aa1.begin(uInt(IFSel_)));// go to IF
  ArrayAccessor<Float, Axis<2> > aa2(aa1);
  aa2.reset(aa2.begin(uInt(polSel_)));// go to pol

  std::vector<float>::iterator it = spectrum.begin();
  for (ArrayAccessor<Float, Axis<3> > i(aa2); i != i.end(); ++i) {
    (*i) = Float(*it);
    it++;
  }
  spec.put(whichRow, arr);
}

void SDMemTable::getSpectrum(Vector<Float>& spectrum, Int whichRow) {
  ROArrayColumn<Float> spec(table_, "SPECTRA");
  Array<Float> arr;
  spec.get(whichRow, arr);
  spectrum.resize(arr.shape()(3));
  ArrayAccessor<Float, Axis<0> > aa0(arr);
  aa0.reset(aa0.begin(uInt(beamSel_)));//go to beam
  ArrayAccessor<Float, Axis<1> > aa1(aa0);
  aa1.reset(aa1.begin(uInt(IFSel_)));// go to IF
  ArrayAccessor<Float, Axis<2> > aa2(aa1);
  aa2.reset(aa2.begin(uInt(polSel_)));// go to pol

  ArrayAccessor<Float, Axis<0> > va(spectrum);
  for (ArrayAccessor<Float, Axis<3> > i(aa2); i != i.end(); ++i) {
    (*va) = (*i);
    va++;
  }
}
/*
void SDMemTable::getMask(Vector<Bool>& mask, Int whichRow) const {
  ROArrayColumn<uChar> spec(table_, "FLAGTRA");
  Array<uChar> arr;
  spec.get(whichRow, arr);
  mask.resize(arr.shape()(3));

  ArrayAccessor<uChar, Axis<0> > aa0(arr);
  aa0.reset(aa0.begin(uInt(beamSel_)));//go to beam
  ArrayAccessor<uChar, Axis<1> > aa1(aa0);
  aa1.reset(aa1.begin(uInt(IFSel_)));// go to IF
  ArrayAccessor<uChar, Axis<2> > aa2(aa1);
  aa2.reset(aa2.begin(uInt(polSel_)));// go to pol

  Bool useUserMask = ( chanMask_.size() == arr.shape()(3) );

  ArrayAccessor<Bool, Axis<0> > va(mask);
  std::vector<bool> tmp;
  tmp = chanMask_; // WHY the fxxx do I have to make a copy here. The
                   // iterator should work on chanMask_??
  std::vector<bool>::iterator miter;
  miter = tmp.begin();

  for (ArrayAccessor<uChar, Axis<3> > i(aa2); i != i.end(); ++i) {
    bool out =!static_cast<bool>(*i);
    if (useUserMask) {
      out = out && (*miter);
      miter++;
    }
    (*va) = out;
    va++;
  }
}
*/
MaskedArray<Float> SDMemTable::rowAsMaskedArray(uInt whichRow,
                                                Bool useSelection) {
  ROArrayColumn<Float> spec(table_, "SPECTRA");
  Array<Float> arr;
  ROArrayColumn<uChar> flag(table_, "FLAGTRA");
  Array<uChar> farr;
  spec.get(whichRow, arr);
  flag.get(whichRow, farr);
  Array<Bool> barr(farr.shape());convertArray(barr, farr);
  MaskedArray<Float> marr;
  if (useSelection) {
    ArrayAccessor<Float, Axis<0> > aa0(arr);
    aa0.reset(aa0.begin(uInt(beamSel_)));//go to beam
    ArrayAccessor<Float, Axis<1> > aa1(aa0);
    aa1.reset(aa1.begin(uInt(IFSel_)));// go to IF
    ArrayAccessor<Float, Axis<2> > aa2(aa1);
    aa2.reset(aa2.begin(uInt(polSel_)));// go to pol

    ArrayAccessor<Bool, Axis<0> > baa0(barr);
    baa0.reset(baa0.begin(uInt(beamSel_)));//go to beam
    ArrayAccessor<Bool, Axis<1> > baa1(baa0);
    baa1.reset(baa1.begin(uInt(IFSel_)));// go to IF
    ArrayAccessor<Bool, Axis<2> > baa2(baa1);
    baa2.reset(baa2.begin(uInt(polSel_)));// go to pol

    Vector<Float> a(arr.shape()(3));
    Vector<Bool> b(barr.shape()(3));
    ArrayAccessor<Float, Axis<0> > a0(a);
    ArrayAccessor<Bool, Axis<0> > b0(b);

    ArrayAccessor<Bool, Axis<3> > j(baa2);
    for (ArrayAccessor<Float, Axis<3> > i(aa2); i != i.end(); ++i) {
      (*a0) = (*i);
      (*b0) = !(*j);
      j++;
      a0++;
      b0++;
    }
    marr.setData(a,b);
  } else {
    marr.setData(arr,!barr);
  }
  return marr;
}

Float SDMemTable::getTsys(Int whichRow) const {
  ROArrayColumn<Float> ts(table_, "TSYS");
  Array<Float> arr;
  ts.get(whichRow, arr);
  Float out;
  IPosition ip(arr.shape());
  ip(0) = beamSel_;ip(1) = IFSel_;ip(2) = polSel_;ip(3)=0;
  out = arr(ip);
  return out;
}

SpectralCoordinate SDMemTable::getCoordinate(uInt whichIdx)  const {

  Table t = table_.keywordSet().asTable("FREQUENCIES");
  if (whichIdx > t.nrow() ) {
    cerr << "SDMemTable::getCoordinate - whichIdx out of range" << endl;
    return SpectralCoordinate();
  }

  Double rp,rv,inc;
  String rf;
  Vector<Double> vec;
  ROScalarColumn<Double> rpc(t, "REFPIX");
  ROScalarColumn<Double> rvc(t, "REFVAL");
  ROScalarColumn<Double> incc(t, "INCREMENT");
  t.keywordSet().get("RESTFREQS",vec);
  t.keywordSet().get("BASEREFFRAME",rf);

  MFrequency::Types mft;
  if (!MFrequency::getType(mft, rf)) {
    cerr << "Frequency type unknown assuming TOPO" << endl;
    mft = MFrequency::TOPO;
  }
  rpc.get(whichIdx, rp);
  rvc.get(whichIdx, rv);
  incc.get(whichIdx, inc);
  SpectralCoordinate spec(mft,rv,inc,rp);
  if (vec.nelements() > 0)
    spec.setRestFrequencies(vec);
  return spec;
}

Bool SDMemTable::setCoordinate(const SpectralCoordinate& speccord,
                               uInt whichIdx) {
  Table t = table_.rwKeywordSet().asTable("FREQUENCIES");
  if (whichIdx > t.nrow() ) {
    throw(AipsError("SDMemTable::setCoordinate - coord no out of range"));
  }
  ScalarColumn<Double> rpc(t, "REFPIX");
  ScalarColumn<Double> rvc(t, "REFVAL");
  ScalarColumn<Double> incc(t, "INCREMENT");

  rpc.put(whichIdx, speccord.referencePixel()[0]);
  rvc.put(whichIdx, speccord.referenceValue()[0]);
  incc.put(whichIdx, speccord.increment()[0]);

  return True;
}

Int SDMemTable::nCoordinates() const
{
  return table_.keywordSet().asTable("FREQUENCIES").nrow();
}

void SDMemTable::setRestFreqs(std::vector<double> freqs, const std::string& theunit)
{
  Vector<Double> tvec(freqs);
  Quantum<Vector<Double> > q(tvec, String(theunit));
  tvec.resize();
  tvec = q.getValue("Hz");
  Table t = table_.keywordSet().asTable("FREQUENCIES");
  t.rwKeywordSet().define("RESTFREQS",tvec);
}

bool SDMemTable::putSDFreqTable(const SDFrequencyTable& sdft) {
  TableDesc td("", "1", TableDesc::Scratch);
  td.addColumn(ScalarColumnDesc<Double>("REFPIX"));
  td.addColumn(ScalarColumnDesc<Double>("REFVAL"));
  td.addColumn(ScalarColumnDesc<Double>("INCREMENT"));
  SetupNewTable aNewTab("freqs", td, Table::New);
  Table aTable (aNewTab, Table::Memory, sdft.length());
  ScalarColumn<Double> sc0(aTable, "REFPIX");
  ScalarColumn<Double> sc1(aTable, "REFVAL");
  ScalarColumn<Double> sc2(aTable, "INCREMENT");
  for (uInt i=0; i < sdft.length(); ++i) {
    sc0.put(i,sdft.referencePixel(i));
    sc1.put(i,sdft.referenceValue(i));
    sc2.put(i,sdft.increment(i));
  }
  String rf = sdft.refFrame();
  if (rf.contains("TOPO")) rf = "TOPO";

  aTable.rwKeywordSet().define("BASEREFFRAME", rf);
  aTable.rwKeywordSet().define("REFFRAME", rf);
  aTable.rwKeywordSet().define("EQUINOX", sdft.equinox());
  aTable.rwKeywordSet().define("UNIT", String(""));
  aTable.rwKeywordSet().define("DOPPLER", String("RADIO"));
  Vector<Double> rfvec;
  aTable.rwKeywordSet().define("RESTFREQS", rfvec);
  table_.rwKeywordSet().defineTable ("FREQUENCIES", aTable);
  return True;
}

SDFrequencyTable SDMemTable::getSDFreqTable() const  {
  SDFrequencyTable sdft;

  return sdft;
}

bool SDMemTable::putSDContainer(const SDContainer& sdc) {
  ScalarColumn<Double> mjd(table_, "TIME");
  ScalarColumn<String> srcn(table_, "SRCNAME");
  ScalarColumn<String> fldn(table_, "FIELDNAME");
  ArrayColumn<Float> spec(table_, "SPECTRA");
  ArrayColumn<uChar> flags(table_, "FLAGTRA");
  ArrayColumn<Float> ts(table_, "TSYS");
  ScalarColumn<Int> scan(table_, "SCANID");
  ScalarColumn<Double> integr(table_, "INTERVAL");
  ArrayColumn<uInt> freqid(table_, "FREQID");
  ArrayColumn<Double> dir(table_, "DIRECTION");
  ScalarColumn<Int> rbeam(table_, "REFBEAM");
  ArrayColumn<Float> tcal(table_, "TCAL");
  ScalarColumn<String> tcalt(table_, "TCALTIME");
  ScalarColumn<Float> az(table_, "AZIMUTH");
  ScalarColumn<Float> el(table_, "ELEVATION");
  ScalarColumn<Float> para(table_, "PARANGLE");

  uInt rno = table_.nrow();
  table_.addRow();

  mjd.put(rno, sdc.timestamp);
  srcn.put(rno, sdc.sourcename);
  fldn.put(rno, sdc.fieldname);
  spec.put(rno, sdc.getSpectrum());
  flags.put(rno, sdc.getFlags());
  ts.put(rno, sdc.getTsys());
  scan.put(rno, sdc.scanid);
  integr.put(rno, sdc.interval);
  freqid.put(rno, sdc.getFreqMap());
  dir.put(rno, sdc.getDirection());
  rbeam.put(rno, sdc.refbeam);
  tcal.put(rno, sdc.tcal);
  tcalt.put(rno, sdc.tcaltime);
  az.put(rno, sdc.azimuth);
  el.put(rno, sdc.elevation);
  para.put(rno, sdc.parangle);

  return true;
}

SDContainer SDMemTable::getSDContainer(uInt whichRow) const {
  ROScalarColumn<Double> mjd(table_, "TIME");
  ROScalarColumn<String> srcn(table_, "SRCNAME");
  ROScalarColumn<String> fldn(table_, "FIELDNAME");
  ROArrayColumn<Float> spec(table_, "SPECTRA");
  ROArrayColumn<uChar> flags(table_, "FLAGTRA");
  ROArrayColumn<Float> ts(table_, "TSYS");
  ROScalarColumn<Int> scan(table_, "SCANID");
  ROScalarColumn<Double> integr(table_, "INTERVAL");
  ROArrayColumn<uInt> freqid(table_, "FREQID");
  ROArrayColumn<Double> dir(table_, "DIRECTION");
  ROScalarColumn<Int> rbeam(table_, "REFBEAM");
  ROArrayColumn<Float> tcal(table_, "TCAL");
  ROScalarColumn<String> tcalt(table_, "TCALTIME");
  ROScalarColumn<Float> az(table_, "AZIMUTH");
  ROScalarColumn<Float> el(table_, "ELEVATION");
  ROScalarColumn<Float> para(table_, "PARANGLE");

  SDContainer sdc(nBeam(),nIF(),nPol(),nChan());
  mjd.get(whichRow, sdc.timestamp);
  srcn.get(whichRow, sdc.sourcename);
  integr.get(whichRow, sdc.interval);
  scan.get(whichRow, sdc.scanid);
  fldn.get(whichRow, sdc.fieldname);
  rbeam.get(whichRow, sdc.refbeam);
  az.get(whichRow, sdc.azimuth);
  el.get(whichRow, sdc.elevation);
  para.get(whichRow, sdc.parangle);
  Vector<Float> tc;
  tcal.get(whichRow, tc);
  sdc.tcal[0] = tc[0];sdc.tcal[1] = tc[1];
  tcalt.get(whichRow, sdc.tcaltime);
  Array<Float> spectrum;
  Array<Float> tsys;
  Array<uChar> flagtrum;
  Vector<uInt> fmap;
  Array<Double> direction;
  spec.get(whichRow, spectrum);
  sdc.putSpectrum(spectrum);
  flags.get(whichRow, flagtrum);
  sdc.putFlags(flagtrum);
  ts.get(whichRow, tsys);
  sdc.putTsys(tsys);
  freqid.get(whichRow, fmap);
  sdc.putFreqMap(fmap);
  dir.get(whichRow, direction);
  sdc.putDirection(direction);
  return sdc;
}

bool SDMemTable::putSDHeader(const SDHeader& sdh) {
  table_.lock();
  table_.rwKeywordSet().define("nIF", sdh.nif);
  table_.rwKeywordSet().define("nBeam", sdh.nbeam);
  table_.rwKeywordSet().define("nPol", sdh.npol);
  table_.rwKeywordSet().define("nChan", sdh.nchan);
  table_.rwKeywordSet().define("Observer", sdh.observer);
  table_.rwKeywordSet().define("Project", sdh.project);
  table_.rwKeywordSet().define("Obstype", sdh.obstype);
  table_.rwKeywordSet().define("AntennaName", sdh.antennaname);
  table_.rwKeywordSet().define("AntennaPosition", sdh.antennaposition);
  table_.rwKeywordSet().define("Equinox", sdh.equinox);
  table_.rwKeywordSet().define("FreqRefFrame", sdh.freqref);
  table_.rwKeywordSet().define("FreqRefVal", sdh.reffreq);
  table_.rwKeywordSet().define("Bandwidth", sdh.bandwidth);
  table_.rwKeywordSet().define("UTC", sdh.utc);
  table_.unlock();
  return true;
}

SDHeader SDMemTable::getSDHeader() const {
  SDHeader sdh;
  table_.keywordSet().get("nBeam",sdh.nbeam);
  table_.keywordSet().get("nIF",sdh.nif);
  table_.keywordSet().get("nPol",sdh.npol);
  table_.keywordSet().get("nChan",sdh.nchan);
  table_.keywordSet().get("Observer", sdh.observer);
  table_.keywordSet().get("Project", sdh.project);
  table_.keywordSet().get("Obstype", sdh.obstype);
  table_.keywordSet().get("AntennaName", sdh.antennaname);
  table_.keywordSet().get("AntennaPosition", sdh.antennaposition);
  table_.keywordSet().get("Equinox", sdh.equinox);
  table_.keywordSet().get("FreqRefFrame", sdh.freqref);
  table_.keywordSet().get("FreqRefVal", sdh.reffreq);
  table_.keywordSet().get("Bandwidth", sdh.bandwidth);
  table_.keywordSet().get("UTC", sdh.utc);
  return sdh;
}
void SDMemTable::makePersistent(const std::string& filename) {
  table_.deepCopy(filename,Table::New);
}

Int SDMemTable::nScan() const {
  Int n = 0;
  ROScalarColumn<Int> scans(table_, "SCANID");
  Int previous = -1;Int current=0;
  for (uInt i=0; i< scans.nrow();i++) {
    scans.getScalar(i,current);
    if (previous != current) {
      previous = current;
      n++;
    }
  }
  return n;
}

String SDMemTable::formatSec(Double x) {
  Double xcop = x;
  MVTime mvt(xcop/24./3600.);  // make days
  if (x < 59.95)
    return  String("   ") + mvt.string(MVTime::TIME_CLEAN_NO_HM, 7)+"s";
  return mvt.string(MVTime::TIME_CLEAN_NO_H, 7)+" ";
};

std::string SDMemTable::summary()  {
  ROScalarColumn<Int> scans(table_, "SCANID");
  ROScalarColumn<String> srcs(table_, "SRCNAME");
  ostringstream oss;
  oss << endl;
  oss << "--------------------------------------------------" << endl;
  oss << " Scan Table Summary" << endl;
  oss << "--------------------------------------------------" << endl;
  oss.flags(std::ios_base::left);
  oss << setw(15) << "Beams:" << setw(4) << nBeam() << endl
      << setw(15) << "IFs:" << setw(4) << nIF() << endl
      << setw(15) << "Polarisations:" << setw(4) << nPol() << endl
      << setw(15) << "Channels:"  << setw(4) << nChan() << endl;
  oss << endl;
  String tmp;
  table_.keywordSet().get("Observer", tmp);
  oss << setw(15) << "Observer:" << tmp << endl;
  table_.keywordSet().get("Project", tmp);
  oss << setw(15) << "Project:" << tmp << endl;
  table_.keywordSet().get("Obstype", tmp);
  oss << setw(15) << "Obs. Type:" << tmp << endl;
  table_.keywordSet().get("AntennaName", tmp);
  oss << setw(15) << "Antenna Name:" << tmp << endl;
  Table t = table_.rwKeywordSet().asTable("FREQUENCIES");
  Vector<Double> vec;
  t.keywordSet().get("RESTFREQS",vec);
  oss << setw(15) << "Rest Freqs:";
  if (vec.nelements() > 0) {
      oss << setprecision(0) << vec << " [Hz]" << endl;
  } else {
      oss << "None set" << endl;
  }
  oss << setw(15) << "Abscissa:" << getAbscissaString() << endl;
  oss << setw(15) << "Cursor:" << "Beam[" << getBeam() << "] "
      << "IF[" << getIF() << "] " << "Pol[" << getPol() << "]" << endl;
  oss << endl;
  uInt count = 0;
  String name;
  Int previous = -1;Int current=0;
  Int integ = 0;
  oss << setw(6) << "Scan"
      << setw(12) << "Source"
      << setw(21) << "Time"
      << setw(11) << "Integration" << endl;
  oss << "--------------------------------------------------" << endl;
  for (uInt i=0; i< scans.nrow();i++) {
    scans.getScalar(i,current);
    if (previous != current) {
      srcs.getScalar(i,name);
      previous = current;
      String t = formatSec(Double(getInterval(i)));
      oss << setw(6) << count << setw(12) << name << setw(21) << getTime(i)
          << setw(2) << setprecision(1)
          << t << endl;
      count++;
    } else {
      integ++;
    }
  }
  oss << endl;
  oss << "Table contains " << table_.nrow() << " integration(s)." << endl;
  oss << "in " << count << " scan(s)." << endl;
  oss << "--------------------------------------------------";
  return String(oss);
}

Int SDMemTable::nBeam() const {
  Int n;
  table_.keywordSet().get("nBeam",n);
  return n;
}
Int SDMemTable::nIF() const {
  Int n;
  table_.keywordSet().get("nIF",n);
  return n;
}
Int SDMemTable::nPol() const {
  Int n;
  table_.keywordSet().get("nPol",n);
  return n;
}
Int SDMemTable::nChan() const {
  Int n;
  table_.keywordSet().get("nChan",n);
  return n;
}
/*
void SDMemTable::maskChannels(const std::vector<Int>& whichChans ) {

  std::vector<int>::iterator it;
  ArrayAccessor<uChar, Axis<2> > j(flags_);
  for (it = whichChans.begin(); it != whichChans.end(); it++) {
    j.reset(j.begin(uInt(*it)));
    for (ArrayAccessor<uChar, Axis<0> > i(j); i != i.end(); ++i) {
      for (ArrayAccessor<uChar, Axis<1> > ii(i); ii != ii.end(); ++ii) {
        for (ArrayAccessor<uChar, Axis<3> > iii(ii);
             iii != iii.end(); ++iii) {
          (*iii) =
        }
      }
    }
  }

}
*/
void SDMemTable::flag(int whichRow) {
  ArrayColumn<uChar> spec(table_, "FLAGTRA");
  Array<uChar> arr;
  spec.get(whichRow, arr);

  ArrayAccessor<uChar, Axis<0> > aa0(arr);
  aa0.reset(aa0.begin(uInt(beamSel_)));//go to beam
  ArrayAccessor<uChar, Axis<1> > aa1(aa0);
  aa1.reset(aa1.begin(uInt(IFSel_)));// go to IF
  ArrayAccessor<uChar, Axis<2> > aa2(aa1);
  aa2.reset(aa2.begin(uInt(polSel_)));// go to pol

  for (ArrayAccessor<uChar, Axis<3> > i(aa2); i != i.end(); ++i) {
    (*i) = uChar(True);
  }

  spec.put(whichRow, arr);
}
