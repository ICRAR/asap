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

#include <map>

#include <casa/aips.h>
#include <casa/iostream.h>
#include <casa/iomanip.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/MaskArrMath.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Arrays/ArrayAccessor.h>
#include <casa/Arrays/Vector.h>
#include <casa/Quanta/MVAngle.h>

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
#include <casa/Quanta/MVAngle.h>

#include "SDDefs.h"
#include "SDMemTable.h"
#include "SDContainer.h"


using namespace casa;
using namespace asap;

SDMemTable::SDMemTable() :
  IFSel_(0),
  beamSel_(0),
  polSel_(0) 
{
  setup();
}

SDMemTable::SDMemTable(const std::string& name) :
  IFSel_(0),
  beamSel_(0),
  polSel_(0)
{
  Table tab(name);
  table_ = tab.copyToMemoryTable("dummy");
  //cerr << "hello from C SDMemTable @ " << this << endl;
}

SDMemTable::SDMemTable(const SDMemTable& other, Bool clear)
{
  IFSel_= other.IFSel_;
  beamSel_= other.beamSel_;
  polSel_= other.polSel_;
  chanMask_ = other.chanMask_;
  table_ = other.table_.copyToMemoryTable(String("dummy"));
  // clear all rows()
  if (clear) {
    table_.removeRow(this->table_.rowNumbers());
  } else {
    IFSel_ = other.IFSel_;
    beamSel_ = other.beamSel_;
    polSel_ = other.polSel_;
  }
  //cerr << "hello from CC SDMemTable @ " << this << endl;
}

SDMemTable::SDMemTable(const Table& tab, const std::string& exprs) :
  IFSel_(0),
  beamSel_(0),
  polSel_(0)
{
  Table t = tableCommand(exprs,tab);
  if (t.nrow() == 0)
      throw(AipsError("Query unsuccessful."));
  table_ = t.copyToMemoryTable("dummy");
}

SDMemTable::~SDMemTable()
{
  //cerr << "goodbye from SDMemTable @ " << this << endl;
}

SDMemTable SDMemTable::getScan(Int scanID) const
{
  String cond("SELECT * from $1 WHERE SCANID == ");
  cond += String::toString(scanID);
  return SDMemTable(table_, cond);
}

SDMemTable &SDMemTable::operator=(const SDMemTable& other)
{
  if (this != &other) {
     IFSel_= other.IFSel_;
     beamSel_= other.beamSel_;
     polSel_= other.polSel_;
     chanMask_.resize(0);
     chanMask_ = other.chanMask_;
     table_ = other.table_.copyToMemoryTable(String("dummy"));
  }
  //cerr << "hello from ASS SDMemTable @ " << this << endl;
  return *this;
}

SDMemTable SDMemTable::getSource(const std::string& source) const
{
  String cond("SELECT * from $1 WHERE SRCNAME == ");
  cond += source;
  return SDMemTable(table_, cond);
}

void SDMemTable::setup()
{
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
  td.addColumn(ArrayColumnDesc<String>("HISTORY"));

  // Now create a new table from the description.

  SetupNewTable aNewTab("dummy", td, Table::New);
  table_ = Table(aNewTab, Table::Memory, 0);
}

std::string SDMemTable::getSourceName(Int whichRow) const
{
  ROScalarColumn<String> src(table_, "SRCNAME");
  String name;
  src.get(whichRow, name);
  return name;
}

std::string SDMemTable::getTime(Int whichRow, Bool showDate) const
{
  Double tm;
  if (whichRow > -1) {
    ROScalarColumn<Double> src(table_, "TIME");
    src.get(whichRow, tm);
  } else {
    table_.keywordSet().get("UTC",tm);
  }
  MVTime mvt(tm);
  if (showDate)
    mvt.setFormat(MVTime::YMD);
  else 
    mvt.setFormat(MVTime::TIME);
  ostringstream oss;
  oss << mvt;
  return String(oss);
}

double SDMemTable::getInterval(Int whichRow) const
{
  ROScalarColumn<Double> src(table_, "INTERVAL");
  Double intval;
  src.get(whichRow, intval);
  return intval;
}

bool SDMemTable::setIF(Int whichIF)
{
  if ( whichIF >= 0 && whichIF < nIF()) {
    IFSel_ = whichIF;
    return true;
  }
  return false;
}

bool SDMemTable::setBeam(Int whichBeam)
{
  if ( whichBeam >= 0 && whichBeam < nBeam()) {
    beamSel_ = whichBeam;
    return true;
  }
  return false;
}

bool SDMemTable::setPol(Int whichPol)
{
  if ( whichPol >= 0 && whichPol < nPol()) {
    polSel_ = whichPol;
    return true;
  }
  return false;
}

bool SDMemTable::setMask(std::vector<int> whichChans)
{
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
  ArrayAccessor<uChar, Axis<asap::BeamAxis> > aa0(arr);
  aa0.reset(aa0.begin(uInt(beamSel_)));//go to beam
  ArrayAccessor<uChar, Axis<asap::IFAxis> > aa1(aa0);
  aa1.reset(aa1.begin(uInt(IFSel_)));// go to IF
  ArrayAccessor<uChar, Axis<asap::PolAxis> > aa2(aa1);
  aa2.reset(aa2.begin(uInt(polSel_)));// go to pol

  Bool useUserMask = ( chanMask_.size() == arr.shape()(3) );

  std::vector<bool> tmp;
  tmp = chanMask_; // WHY the fxxx do I have to make a copy here
  std::vector<bool>::iterator miter;
  miter = tmp.begin();

  for (ArrayAccessor<uChar, Axis<asap::ChanAxis> > i(aa2); i != i.end(); ++i) {
    bool out =!static_cast<bool>(*i);
    if (useUserMask) {
      out = out && (*miter);
      miter++;
    }
    mask.push_back(out);
  }
  return mask;
}
std::vector<float> SDMemTable::getSpectrum(Int whichRow) const 
{
  std::vector<float> spectrum;
  ROArrayColumn<Float> spec(table_, "SPECTRA");
  Array<Float> arr;
  spec.get(whichRow, arr);
  ArrayAccessor<Float, Axis<asap::BeamAxis> > aa0(arr);
  aa0.reset(aa0.begin(uInt(beamSel_)));//go to beam
  ArrayAccessor<Float, Axis<asap::IFAxis> > aa1(aa0);
  aa1.reset(aa1.begin(uInt(IFSel_)));// go to IF
  ArrayAccessor<Float, Axis<asap::PolAxis> > aa2(aa1);
  aa2.reset(aa2.begin(uInt(polSel_)));// go to pol
  for (ArrayAccessor<Float, Axis<asap::ChanAxis> > i(aa2); i != i.end(); ++i) {
    spectrum.push_back(*i);
  }
  return spectrum;
}
std::vector<string> SDMemTable::getCoordInfo() const
{
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
  t.keywordSet().get("BASEREFFRAME",rfrm);
  inf.push_back(rfrm);
  return inf;
}

void SDMemTable::setCoordInfo(std::vector<string> theinfo)
{
  std::vector<string>::iterator it;
  String un,rfrm,dpl;
  un = theinfo[0];
  rfrm = theinfo[1];
  dpl = theinfo[2];

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
//
  MFrequency::Types mdr;
  if (!MFrequency::getType(mdr, rfrm)) {
    
    Int a,b;const uInt* c;
    const String* valid = MFrequency::allMyTypes(a, b, c);
    String pfix = "Please specify a legal frame type. Types are\n";
    throw(AipsError(pfix+(*valid)));
  } else {
    t.rwKeywordSet().define("REFFRAME",rfrm);
  }
//
  MDoppler::Types dtype;
  dpl.upcase();
  if (!MDoppler::getType(dtype, dpl)) {
    throw(AipsError("Doppler type unknown"));
  } else {
    t.rwKeywordSet().define("DOPPLER",dpl);
  }
}


std::vector<double> SDMemTable::getAbcissa(Int whichRow) const
{
  std::vector<double> abc(nChan());

// Get header units keyword

  Table t = table_.keywordSet().asTable("FREQUENCIES");
  String sunit;
  t.keywordSet().get("UNIT",sunit);
  if (sunit == "") sunit = "pixel";
  Unit u(sunit);

// Easy if just wanting pixels

  if (sunit==String("pixel")) {
    // assume channels/pixels
    std::vector<double>::iterator it;
    uInt i=0;
    for (it = abc.begin(); it != abc.end(); ++it) {
      (*it) = Double(i++);
    }
//
    return abc;
  }

// Continue with km/s or Hz.  Get FreqID

  ROArrayColumn<uInt> fid(table_, "FREQID");
  Vector<uInt> v;
  fid.get(whichRow, v);
  uInt specidx = v(IFSel_);

// Get SpectralCoordinate, set reference frame conversion,
// velocity conversion, and rest freq state

  SpectralCoordinate spc = getSpectralCoordinate(specidx, whichRow);
//
  Vector<Double> pixel(nChan());
  indgen(pixel);
//
  if (u == Unit("km/s")) {
     Vector<Double> world;
     spc.pixelToVelocity(world,pixel);
     std::vector<double>::iterator it;
     uInt i = 0;
     for (it = abc.begin(); it != abc.end(); ++it) {
       (*it) = world[i];
       i++;
     }
  } else if (u == Unit("Hz")) {

// Set world axis units

    Vector<String> wau(1); wau = u.getName();
    spc.setWorldAxisUnits(wau);
//
    std::vector<double>::iterator it;
    Double tmp;
    uInt i = 0;
    for (it = abc.begin(); it != abc.end(); ++it) {
      spc.toWorld(tmp,pixel[i]);
      (*it) = tmp;
      i++;
    }
  }
  return abc;
}

std::string SDMemTable::getAbcissaString(Int whichRow) const
{
  ROArrayColumn<uInt> fid(table_, "FREQID");
  Table t = table_.keywordSet().asTable("FREQUENCIES");
//
  String sunit;
  t.keywordSet().get("UNIT",sunit);
  if (sunit == "") sunit = "pixel";
  Unit u(sunit);
//
  Vector<uInt> v;
  fid.get(whichRow, v);
  uInt specidx = v(IFSel_);

// Get SpectralCoordinate, with frame, velocity, rest freq state set

  SpectralCoordinate spc = getSpectralCoordinate(specidx, whichRow);
//
  String s = "Channel";
  if (u == Unit("km/s")) { 
    s = CoordinateUtil::axisLabel(spc,0,True,True,True);
  } else if (u == Unit("Hz")) {
    Vector<String> wau(1);wau = u.getName();
    spc.setWorldAxisUnits(wau);
//
    s = CoordinateUtil::axisLabel(spc,0,True,True,False);
  }
  return s;
}

void SDMemTable::setSpectrum(std::vector<float> spectrum, int whichRow)
{
  ArrayColumn<Float> spec(table_, "SPECTRA");
  Array<Float> arr;
  spec.get(whichRow, arr);
  if (spectrum.size() != arr.shape()(3)) {
    throw(AipsError("Attempting to set spectrum with incorrect length."));
  }

  ArrayAccessor<Float, Axis<asap::BeamAxis> > aa0(arr);
  aa0.reset(aa0.begin(uInt(beamSel_)));//go to beam
  ArrayAccessor<Float, Axis<asap::IFAxis> > aa1(aa0);
  aa1.reset(aa1.begin(uInt(IFSel_)));// go to IF
  ArrayAccessor<Float, Axis<asap::PolAxis> > aa2(aa1);
  aa2.reset(aa2.begin(uInt(polSel_)));// go to pol

  std::vector<float>::iterator it = spectrum.begin();
  for (ArrayAccessor<Float, Axis<asap::ChanAxis> > i(aa2); i != i.end(); ++i) {
    (*i) = Float(*it);
    it++;
  }
  spec.put(whichRow, arr);
}

void SDMemTable::getSpectrum(Vector<Float>& spectrum, Int whichRow) const
{
  ROArrayColumn<Float> spec(table_, "SPECTRA");
  Array<Float> arr;
  spec.get(whichRow, arr);
  spectrum.resize(arr.shape()(3));
  ArrayAccessor<Float, Axis<asap::BeamAxis> > aa0(arr);
  aa0.reset(aa0.begin(uInt(beamSel_)));//go to beam
  ArrayAccessor<Float, Axis<asap::IFAxis> > aa1(aa0);
  aa1.reset(aa1.begin(uInt(IFSel_)));// go to IF
  ArrayAccessor<Float, Axis<asap::PolAxis> > aa2(aa1);
  aa2.reset(aa2.begin(uInt(polSel_)));// go to pol

  ArrayAccessor<Float, Axis<asap::BeamAxis> > va(spectrum);
  for (ArrayAccessor<Float, Axis<asap::ChanAxis> > i(aa2); i != i.end(); ++i) {
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

  ArrayAccessor<uChar, Axis<asap::BeamAxis> > aa0(arr);
  aa0.reset(aa0.begin(uInt(beamSel_)));//go to beam
  ArrayAccessor<uChar, Axis<asap::IFAxis> > aa1(aa0);
  aa1.reset(aa1.begin(uInt(IFSel_)));// go to IF
  ArrayAccessor<uChar, Axis<asap::PolAxis> > aa2(aa1);
  aa2.reset(aa2.begin(uInt(polSel_)));// go to pol

  Bool useUserMask = ( chanMask_.size() == arr.shape()(3) );

  ArrayAccessor<Bool, Axis<asap::BeamAxis> > va(mask);
  std::vector<bool> tmp;
  tmp = chanMask_; // WHY the fxxx do I have to make a copy here. The
                   // iterator should work on chanMask_??
  std::vector<bool>::iterator miter;
  miter = tmp.begin();

  for (ArrayAccessor<uChar, Axis<asap::ChanAxis> > i(aa2); i != i.end(); ++i) {
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
                                                Bool useSelection) const 
{
  ROArrayColumn<Float> spec(table_, "SPECTRA");
  Array<Float> arr;
  ROArrayColumn<uChar> flag(table_, "FLAGTRA");
  Array<uChar> farr;
  spec.get(whichRow, arr);
  flag.get(whichRow, farr);
  Array<Bool> barr(farr.shape());convertArray(barr, farr);
  MaskedArray<Float> marr;
  if (useSelection) {
    ArrayAccessor<Float, Axis<asap::BeamAxis> > aa0(arr);
    aa0.reset(aa0.begin(uInt(beamSel_)));//go to beam
    ArrayAccessor<Float, Axis<asap::IFAxis> > aa1(aa0);
    aa1.reset(aa1.begin(uInt(IFSel_)));// go to IF
    ArrayAccessor<Float, Axis<asap::PolAxis> > aa2(aa1);
    aa2.reset(aa2.begin(uInt(polSel_)));// go to pol

    ArrayAccessor<Bool, Axis<asap::BeamAxis> > baa0(barr);
    baa0.reset(baa0.begin(uInt(beamSel_)));//go to beam
    ArrayAccessor<Bool, Axis<asap::IFAxis> > baa1(baa0);
    baa1.reset(baa1.begin(uInt(IFSel_)));// go to IF
    ArrayAccessor<Bool, Axis<asap::PolAxis> > baa2(baa1);
    baa2.reset(baa2.begin(uInt(polSel_)));// go to pol

    Vector<Float> a(arr.shape()(3));
    Vector<Bool> b(barr.shape()(3));
    ArrayAccessor<Float, Axis<asap::BeamAxis> > a0(a);
    ArrayAccessor<Bool, Axis<asap::BeamAxis> > b0(b);

    ArrayAccessor<Bool, Axis<asap::ChanAxis> > j(baa2);
    for (ArrayAccessor<Float, Axis<asap::ChanAxis> > i(aa2); 
	 i != i.end(); ++i) {
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

Float SDMemTable::getTsys(Int whichRow) const
{
  ROArrayColumn<Float> ts(table_, "TSYS");
  Array<Float> arr;
  ts.get(whichRow, arr);
  Float out;
  IPosition ip(arr.shape());
  ip(0) = beamSel_;ip(1) = IFSel_;ip(2) = polSel_;ip(3)=0;
  out = arr(ip);
  return out;
}

MDirection SDMemTable::getDirection(Int whichRow, Bool refBeam) const
{
  MDirection::Types mdr = getDirectionReference();
  ROArrayColumn<Double> dir(table_, "DIRECTION");
  Array<Double> posit;
  dir.get(whichRow,posit);
  Vector<Double> wpos(2);
  wpos[0] = posit(IPosition(2,beamSel_,0));
  wpos[1] = posit(IPosition(2,beamSel_,1));
  Quantum<Double> lon(wpos[0],Unit(String("rad")));
  Quantum<Double> lat(wpos[1],Unit(String("rad")));
  return MDirection(lon, lat, mdr);
}

MEpoch SDMemTable::getEpoch (Int whichRow) const
{
  MEpoch::Types met = getTimeReference();
//
  ROScalarColumn<Double> tme(table_, "TIME");
  Double obstime;
  tme.get(whichRow,obstime);
  MVEpoch tm2(Quantum<Double>(obstime, Unit(String("d"))));
  return MEpoch(tm2, met);
}

MPosition SDMemTable::getAntennaPosition () const
{
  Vector<Double> antpos;
  table_.keywordSet().get("AntennaPosition", antpos);
  MVPosition mvpos(antpos(0),antpos(1),antpos(2));
  return MPosition(mvpos);
}


SpectralCoordinate SDMemTable::getSpectralCoordinate(uInt whichIdx) const
{
  
  Table t = table_.keywordSet().asTable("FREQUENCIES");
  if (whichIdx > t.nrow() ) {
    throw(AipsError("SDMemTable::getSpectralCoordinate - whichIdx out of range"));
  }

  Double rp,rv,inc;
  String rf;
  ROScalarColumn<Double> rpc(t, "REFPIX");
  ROScalarColumn<Double> rvc(t, "REFVAL");
  ROScalarColumn<Double> incc(t, "INCREMENT");
  t.keywordSet().get("BASEREFFRAME",rf);

// Create SpectralCoordinate (units Hz)

  MFrequency::Types mft;
  if (!MFrequency::getType(mft, rf)) {
    cerr << "Frequency type unknown assuming TOPO" << endl;
    mft = MFrequency::TOPO;
  }
  rpc.get(whichIdx, rp);
  rvc.get(whichIdx, rv);
  incc.get(whichIdx, inc);
//
  SpectralCoordinate spec(mft,rv,inc,rp);
//
  return spec;
}


SpectralCoordinate SDMemTable::getSpectralCoordinate(uInt whichIdx, uInt whichRow) const
{
// Create basic SC

  SpectralCoordinate spec = getSpectralCoordinate (whichIdx);
//
  Table t = table_.keywordSet().asTable("FREQUENCIES");

// Set rest frequencies

  Vector<Double> vec;
  t.keywordSet().get("RESTFREQS",vec);
  if (vec.nelements() > 0) {
    spec.setRestFrequencies(vec);

// Select rest freq

    if (vec.nelements() >= nIF()) {
       spec.selectRestFrequency(uInt(IFSel_));
    }
  }

// Set up frame conversion layer

  String frm;
  t.keywordSet().get("REFFRAME",frm);
  if (frm == "") frm = "TOPO";
  MFrequency::Types mtype;
  if (!MFrequency::getType(mtype, frm)) {
    cout << "Frequency type unknown assuming TOPO" << endl;       // SHould never happen
    mtype = MFrequency::TOPO;
  }

// Set reference frame conversion  (requires row)

  MDirection direct = getDirection(whichRow);
  MEpoch epoch = getEpoch(whichRow);
  MPosition pos = getAntennaPosition();
  if (!spec.setReferenceConversion(mtype,epoch,pos,direct)) {
    throw(AipsError("Couldn't convert frequency frame."));
  }

// Now velocity conversion if appropriate

  String unitStr;
  t.keywordSet().get("UNIT",unitStr);
//
  String dpl;
  t.keywordSet().get("DOPPLER",dpl);
  MDoppler::Types dtype;
  MDoppler::getType(dtype, dpl);

// Only set velocity unit if non-blank and non-Hz

  if (!unitStr.empty()) {
     Unit unitU(unitStr);
     if (unitU==Unit("Hz")) {
     } else {
        spec.setVelocity(unitStr, dtype);
     }
  }
//
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

void SDMemTable::setRestFreqs(std::vector<double> freqs, 
			      const std::string& theunit)
{
  Vector<Double> tvec(freqs);
  Quantum<Vector<Double> > q(tvec, String(theunit));
  tvec.resize();
  tvec = q.getValue("Hz");
  Table t = table_.keywordSet().asTable("FREQUENCIES");
  t.rwKeywordSet().define("RESTFREQS",tvec);
}

std::vector<double> SDMemTable::getRestFreqs() const
{
  Table t = table_.keywordSet().asTable("FREQUENCIES");
  Vector<Double> tvec;
  t.keywordSet().get("RESTFREQS",tvec);
  std::vector<double> stlout;
  tvec.tovector(stlout);
  return stlout;  
}

bool SDMemTable::putSDFreqTable(const SDFrequencyTable& sdft)
{
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
  String rfunit;
  sdft.restFrequencies(rfvec,rfunit);
  Quantum<Vector<Double> > q(rfvec, rfunit);
  rfvec.resize();
  rfvec = q.getValue("Hz");
  aTable.rwKeywordSet().define("RESTFREQS", rfvec);
  table_.rwKeywordSet().defineTable ("FREQUENCIES", aTable);
  return True;
}

SDFrequencyTable SDMemTable::getSDFreqTable() const
{
  // TODO !!!!! implement this properly USE with care
  const Table& t = table_.keywordSet().asTable("FREQUENCIES");
  SDFrequencyTable sdft;
  sdft.setLength(t.nrow());
  return sdft;
}

bool SDMemTable::putSDContainer(const SDContainer& sdc)
{
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
  ArrayColumn<String> hist(table_, "HISTORY");

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
  hist.put(rno, sdc.getHistory());

  return true;
}

SDContainer SDMemTable::getSDContainer(uInt whichRow) const
{
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
  ROArrayColumn<String> hist(table_, "HISTORY");

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
  Vector<String> histo;
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
  hist.get(whichRow, histo);
  sdc.putHistory(histo);
  return sdc;
}

bool SDMemTable::putSDHeader(const SDHeader& sdh)
{
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
  table_.rwKeywordSet().define("FluxUnit", sdh.fluxunit);
  table_.rwKeywordSet().define("Epoch", sdh.epoch);
  return true;
}

SDHeader SDMemTable::getSDHeader() const
{
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
  table_.keywordSet().get("FluxUnit", sdh.fluxunit);
  table_.keywordSet().get("Epoch", sdh.epoch);
  return sdh;
}
void SDMemTable::makePersistent(const std::string& filename)
{
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

String SDMemTable::formatSec(Double x) const
{
  Double xcop = x;
  MVTime mvt(xcop/24./3600.);  // make days
  
  if (x < 59.95)
    return  String("      ") + mvt.string(MVTime::TIME_CLEAN_NO_HM, 7)+"s";
  else if (x < 3599.95) 
    return String("   ") + mvt.string(MVTime::TIME_CLEAN_NO_H,7)+" ";
  else {
    ostringstream oss;
    oss << setw(2) << std::right << setprecision(1) << mvt.hour();
    oss << ":" << mvt.string(MVTime::TIME_CLEAN_NO_H,7) << " ";
    return String(oss);
  }   
};

String SDMemTable::formatDirection(const MDirection& md) const
{
  Vector<Double> t = md.getAngle(Unit(String("rad"))).getValue();
  Int prec = 7;

  MVAngle mvLon(t[0]);
  String sLon = mvLon.string(MVAngle::TIME,prec);
  MVAngle mvLat(t[1]);
  String sLat = mvLat.string(MVAngle::ANGLE+MVAngle::DIG2,prec);

   return sLon + String(" ") + sLat;
}


std::string SDMemTable::getFluxUnit() const
{
  String tmp;
  table_.keywordSet().get("FluxUnit", tmp);
  return tmp;
}

void SDMemTable::setFluxUnit(const std::string& unit)
{
  String tmp(unit);
  Unit tU(tmp);
  if (tU==Unit("K") || tU==Unit("Jy")) {
     table_.rwKeywordSet().define(String("FluxUnit"), tmp);
  } else {
     throw AipsError("Illegal unit - must be compatible with Jy or K");
  }
}


void SDMemTable::setInstrument(const std::string& name)
{
  Bool throwIt = True;
  Instrument ins = convertInstrument (name, throwIt);
  String nameU(name);
  nameU.upcase();
  table_.rwKeywordSet().define(String("AntennaName"), nameU);
}

std::string SDMemTable::summary() const  {
  ROScalarColumn<Int> scans(table_, "SCANID");
  ROScalarColumn<String> srcs(table_, "SRCNAME");

  // get number of integrations per scan
  int cIdx = 0;
  int idx = 0;
  int scount = 0;
  std::vector<int> cycles;
  for (uInt i=0; i<scans.nrow();++i) {
    while (idx == cIdx && i<scans.nrow()) {
      scans.getScalar(++i,cIdx);
      ++scount;
    }
    idx = cIdx;
    cycles.push_back(scount);
    scount=0;
    --i;
  }  


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
  oss << setw(15) << "Obs Date:" << getTime(-1,True) << endl;
  table_.keywordSet().get("Project", tmp);
  oss << setw(15) << "Project:" << tmp << endl;
  table_.keywordSet().get("Obstype", tmp);
  oss << setw(15) << "Obs. Type:" << tmp << endl;
  table_.keywordSet().get("AntennaName", tmp);
  oss << setw(15) << "Antenna Name:" << tmp << endl;
  table_.keywordSet().get("FluxUnit", tmp);
  oss << setw(15) << "Flux Unit:" << tmp << endl;
  Table t = table_.keywordSet().asTable("FREQUENCIES");
  Vector<Double> vec;
  t.keywordSet().get("RESTFREQS",vec);
  oss << setw(15) << "Rest Freqs:";
  if (vec.nelements() > 0) {
      oss << setprecision(0) << vec << " [Hz]" << endl;
  } else {
      oss << "None set" << endl;
  }
  oss << setw(15) << "Abcissa:" << getAbcissaString() << endl;
  oss << setw(15) << "Cursor:" << "Beam[" << getBeam() << "] "
      << "IF[" << getIF() << "] " << "Pol[" << getPol() << "]" << endl;
  oss << endl;
  uInt count = 0;
  String name;
  Int previous = -1;Int current=0;
  Int integ = 0;
  String dirtype ="Position ("+
    MDirection::showType(getDirectionReference())+
    ")";
  oss << setw(6) << "Scan"
      << setw(15) << "Source"
      << setw(26) << dirtype
      << setw(10) << "Time"
      << setw(13) << "Integration" << endl;
  oss << "-------------------------------------------------------------------------------" << endl;
  
  std::vector<int>::iterator it = cycles.begin();
  for (uInt i=0; i< scans.nrow();i++) {
    scans.getScalar(i,current);
    if (previous != current) {
      srcs.getScalar(i,name);
      previous = current;
      String t = formatSec(Double(getInterval(i)));
      String posit = formatDirection(getDirection(i,True));
      oss << setw(6) << count 
	  << setw(15) << name
	  << setw(26) << posit
	  << setw(10) << getTime(i,False)
	  << setw(3) << std::right << *it << " x "
          << setw(10) 
          << t << std::left << endl;
      count++;
      it++;
    } else {
      integ++;
    }
  }
  oss << endl;
  oss << "Table contains " << table_.nrow() << " integration(s)." << endl;
  oss << "in " << count << " scan(s)." << endl;
  oss << "-------------------------------------------------------------------------------";
  return String(oss);
}

Int SDMemTable::nBeam() const
{
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
bool SDMemTable::appendHistory(const std::string& hist, int whichRow)
{
  ArrayColumn<String> histo(table_, "HISTORY");
  Vector<String> history;
  histo.get(whichRow, history);
  history.resize(history.nelements()+1,True);
  history[history.nelements()-1] = hist;
  histo.put(whichRow, history);
}

std::vector<std::string> SDMemTable::history(int whichRow) const
{
  ROArrayColumn<String> hist(table_, "HISTORY");
  Vector<String> history;
  hist.get(whichRow, history);
  std::vector<std::string> stlout;
  // there is no Array<String>.tovector(std::vector<std::string>), so
  // do it by hand
  for (uInt i=0; i<history.nelements(); ++i) {
    stlout.push_back(history[i]);
  }
  return stlout;
}
/*
void SDMemTable::maskChannels(const std::vector<Int>& whichChans ) {

  std::vector<int>::iterator it;
  ArrayAccessor<uChar, Axis<asap::PolAxis> > j(flags_);
  for (it = whichChans.begin(); it != whichChans.end(); it++) {
    j.reset(j.begin(uInt(*it)));
    for (ArrayAccessor<uChar, Axis<asap::BeamAxis> > i(j); i != i.end(); ++i) {
      for (ArrayAccessor<uChar, Axis<asap::IFAxis> > ii(i); ii != ii.end(); ++ii) {
        for (ArrayAccessor<uChar, Axis<asap::ChanAxis> > iii(ii);
             iii != iii.end(); ++iii) {
          (*iii) =
        }
      }
    }
  }

}
*/
void SDMemTable::flag(int whichRow)
{
  ArrayColumn<uChar> spec(table_, "FLAGTRA");
  Array<uChar> arr;
  spec.get(whichRow, arr);

  ArrayAccessor<uChar, Axis<asap::BeamAxis> > aa0(arr);
  aa0.reset(aa0.begin(uInt(beamSel_)));//go to beam
  ArrayAccessor<uChar, Axis<asap::IFAxis> > aa1(aa0);
  aa1.reset(aa1.begin(uInt(IFSel_)));// go to IF
  ArrayAccessor<uChar, Axis<asap::PolAxis> > aa2(aa1);
  aa2.reset(aa2.begin(uInt(polSel_)));// go to pol

  for (ArrayAccessor<uChar, Axis<asap::ChanAxis> > i(aa2); i != i.end(); ++i) {
    (*i) = uChar(True);
  }

  spec.put(whichRow, arr);
}

MDirection::Types SDMemTable::getDirectionReference() const
{  
  Float eq;
  table_.keywordSet().get("Equinox",eq);
  std::map<float,string> mp;
  mp[2000.0] = "J2000";
  mp[1950.0] = "B1950";
  MDirection::Types mdr;
  if (!MDirection::getType(mdr, mp[eq])) {   
    mdr = MDirection::J2000;
    cerr  << "Unknown equinox using J2000" << endl;
  }
//
  return mdr;
}

MEpoch::Types SDMemTable::getTimeReference () const
{
  MEpoch::Types met;
  String ep;
  table_.keywordSet().get("Epoch",ep);
  if (!MEpoch::getType(met, ep)) {
    cerr << "Epoch type uknown - using UTC" << endl;
    met = MEpoch::UTC;
  }
//
  return met;
}


Instrument SDMemTable::convertInstrument (const String& instrument,
                                          Bool throwIt)
{
   String t(instrument);
   t.upcase();

// The strings are what SDReader returns, after cunning interrogation
// of station names... :-(

   Instrument inst = asap::UNKNOWN;
   if (t==String("DSS-43")) {                
      inst = TIDBINBILLA;
   } else if (t==String("ATPKSMB")) {
      inst = ATPKSMB;
   } else if (t==String("ATPKSHOH")) {
      inst = ATPKSHOH;
   } else if (t==String("ATMOPRA")) {
      inst = ATMOPRA;
   } else if (t==String("CEDUNA")) {
      inst = CEDUNA;
   } else if (t==String("HOBART")) {
      inst = HOBART;
   } else {
      if (throwIt) {
         throw AipsError("Unrecognized instrument - use function scan.set_instrument to set");
      }
   }
   return inst;
}

