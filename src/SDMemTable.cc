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
#include <casa/BasicMath/Math.h>
#include <casa/Quanta/MVAngle.h>

#include <tables/Tables/TableParse.h>
#include <tables/Tables/TableDesc.h>
#include <tables/Tables/SetupNewTab.h>
#include <tables/Tables/ScaColDesc.h>
#include <tables/Tables/ArrColDesc.h>

#include <tables/Tables/ExprNode.h>
#include <tables/Tables/TableRecord.h>
#include <measures/Measures/MFrequency.h>
#include <measures/Measures/MeasTable.h>
#include <coordinates/Coordinates/CoordinateUtil.h>
#include <casa/Quanta/MVTime.h>
#include <casa/Quanta/MVAngle.h>

#include "SDDefs.h"
#include "SDMemTable.h"
#include "SDContainer.h"
#include "MathUtils.h"
#include "SDPol.h"


using namespace casa;
using namespace asap;

SDMemTable::SDMemTable() :
  IFSel_(0),
  beamSel_(0),
  polSel_(0) 
{
  setup();
  attach();
}

SDMemTable::SDMemTable(const std::string& name) :
  IFSel_(0),
  beamSel_(0),
  polSel_(0)
{
  Table tab(name);
  table_ = tab.copyToMemoryTable("dummy");
  //cerr << "hello from C SDMemTable @ " << this << endl;
  attach();
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
//
  attach();
  //cerr << "hello from CC SDMemTable @ " << this << endl;
}

SDMemTable::SDMemTable(const Table& tab, const std::string& exprs) :
  IFSel_(0),
  beamSel_(0),
  polSel_(0)
{
  cout << exprs << endl;
  Table t = tableCommand(exprs,tab);
  if (t.nrow() == 0)
      throw(AipsError("Query unsuccessful."));
  table_ = t.copyToMemoryTable("dummy");
  attach();
  renumber();
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
     attach();
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
//
  td.addColumn(ScalarColumnDesc<Double>("TIME"));
  td.addColumn(ScalarColumnDesc<String>("SRCNAME"));
  td.addColumn(ArrayColumnDesc<Float>("SPECTRA"));
  td.addColumn(ArrayColumnDesc<uChar>("FLAGTRA"));
  td.addColumn(ArrayColumnDesc<Float>("TSYS"));
  td.addColumn(ArrayColumnDesc<Float>("STOKES"));
  td.addColumn(ScalarColumnDesc<Int>("SCANID"));
  td.addColumn(ScalarColumnDesc<Double>("INTERVAL"));
  td.addColumn(ArrayColumnDesc<uInt>("FREQID"));
  td.addColumn(ArrayColumnDesc<uInt>("RESTFREQID"));
  td.addColumn(ArrayColumnDesc<Double>("DIRECTION"));
  td.addColumn(ScalarColumnDesc<String>("FIELDNAME"));
  td.addColumn(ScalarColumnDesc<String>("TCALTIME"));
  td.addColumn(ArrayColumnDesc<Float>("TCAL"));
  td.addColumn(ScalarColumnDesc<Float>("AZIMUTH"));
  td.addColumn(ScalarColumnDesc<Float>("ELEVATION"));
  td.addColumn(ScalarColumnDesc<Float>("PARANGLE"));
  td.addColumn(ScalarColumnDesc<Int>("REFBEAM"));
  td.addColumn(ArrayColumnDesc<String>("HISTORY"));

  // Now create Table SetUp from the description.

  SetupNewTable aNewTab("dummy", td, Table::New);

// Bind the Stokes Virtual machine to the STOKES column
// Because we don't know how many polarizations will
// be in the data at this point, we must bind the
// Virtual Engine regardless.  The STOKES column won't
// be accessed if not appropriate (nPol=4)

   SDStokesEngine stokesEngine(String("STOKES"), String("SPECTRA"));
   aNewTab.bindColumn ("STOKES", stokesEngine);

// Create Table

  table_ = Table(aNewTab, Table::Memory, 0);
}


void SDMemTable::attach()
{
  timeCol_.attach(table_, "TIME");
  srcnCol_.attach(table_, "SRCNAME");
  specCol_.attach(table_, "SPECTRA");
  flagsCol_.attach(table_, "FLAGTRA");
  tsCol_.attach(table_, "TSYS");
  stokesCol_.attach(table_, "STOKES");
  scanCol_.attach(table_, "SCANID");
  integrCol_.attach(table_, "INTERVAL");
  freqidCol_.attach(table_, "FREQID");
  restfreqidCol_.attach(table_, "RESTFREQID");
  dirCol_.attach(table_, "DIRECTION");
  fldnCol_.attach(table_, "FIELDNAME");
  tcaltCol_.attach(table_, "TCALTIME");
  tcalCol_.attach(table_, "TCAL");
  azCol_.attach(table_, "AZIMUTH");
  elCol_.attach(table_, "ELEVATION");
  paraCol_.attach(table_, "PARANGLE");
  rbeamCol_.attach(table_, "REFBEAM");
  histCol_.attach(table_, "HISTORY");
}


std::string SDMemTable::getSourceName(Int whichRow) const
{
  String name;
  srcnCol_.get(whichRow, name);
  return name;
}

std::string SDMemTable::getTime(Int whichRow, Bool showDate) const
{
  Double tm;
  if (whichRow > -1) {
    timeCol_.get(whichRow, tm);
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
  Double intval;
  integrCol_.get(whichRow, intval);
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

void SDMemTable::resetCursor() 
{
   polSel_ = 0;
   IFSel_ = 0;
   beamSel_ = 0;
}

bool SDMemTable::setMask(std::vector<int> whichChans)
{
  std::vector<int>::iterator it;
  uInt n = flagsCol_.shape(0)(3);
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
  Array<uChar> arr;
  flagsCol_.get(whichRow, arr);
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
  Array<Float> arr;
  specCol_.get(whichRow, arr);
//
  return getFloatSpectrum (arr);
}

std::vector<float> SDMemTable::getStokesSpectrum(Int whichRow, Bool doPol, Float paOffset) const 
//
// Gets
//  doPol=False  : I,Q,U,V
//  doPol=True   : I,P,PA,V   ; P = sqrt(Q**2+U**2), PA = 0.5*atan2(Q,U)
//
{
  if (nPol()!=1 && nPol()!=2 && nPol()!=4) {
     throw (AipsError("You must have 1,2 or 4 polarizations to get the Stokes parameters"));
  }
  Array<Float> arr;
  stokesCol_.get(whichRow, arr);
//
  if (doPol && (polSel_==1 || polSel_==2)) {
     const IPosition& shape = arr.shape();
     IPosition start(asap::nAxes,0);
     IPosition end(shape-1);
//
     start(asap::PolAxis) = 1;                       // Q
     end (asap::PolAxis) = 1;
     Array<Float> Q = arr(start,end);
//
     start(asap::PolAxis) = 2;                       // U
     end (asap::PolAxis) = 2;
     Array<Float> U = arr(start,end);
//
     Array<Float> out;
     if (polSel_==1) {                                        // P
        out = SDPolUtil::polarizedIntensity(Q,U);
     } else if (polSel_==2) {                                 // P.A.
        out = SDPolUtil::positionAngle(Q,U) + paOffset;
     }
//
     IPosition vecShape(1,shape(asap::ChanAxis));
     Vector<Float> outV = out.reform(vecShape);
     std::vector<float> spectrum(out.nelements());
     for (uInt i=0; i<out.nelements(); i++) {
        spectrum[i] = outV[i];
     }
     return spectrum;
  } else {
     return getFloatSpectrum (arr);
  }
}

std::vector<float> SDMemTable::getFloatSpectrum (const Array<Float>& arr) const
{

// Iterate and extract

  ArrayAccessor<Float, Axis<asap::BeamAxis> > aa0(arr);
  aa0.reset(aa0.begin(uInt(beamSel_)));//go to beam
  ArrayAccessor<Float, Axis<asap::IFAxis> > aa1(aa0);
  aa1.reset(aa1.begin(uInt(IFSel_)));// go to IF
  ArrayAccessor<Float, Axis<asap::PolAxis> > aa2(aa1);
  aa2.reset(aa2.begin(uInt(polSel_)));// go to pol
//
  std::vector<float> spectrum;
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
  String un,rfrm, brfrm,dpl;
  un = theinfo[0];              // Abcissa unit
  rfrm = theinfo[1];            // Active (or conversion) frame
  dpl = theinfo[2];             // Doppler
  brfrm = theinfo[3];           // Base frame
//
  Table t = table_.rwKeywordSet().asTable("FREQUENCIES");
//
  Vector<Double> rstf;
  t.keywordSet().get("RESTFREQS",rstf);
//
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
//
  if (!MFrequency::getType(mdr, brfrm)) {
     Int a,b;const uInt* c;
     const String* valid = MFrequency::allMyTypes(a, b, c);
     String pfix = "Please specify a legal frame type. Types are\n";
     throw(AipsError(pfix+(*valid)));
   } else {
    t.rwKeywordSet().define("BASEREFFRAME",brfrm);
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

// Continue with km/s or Hz.  Get FreqIDs

  Vector<uInt> freqIDs;
  freqidCol_.get(whichRow, freqIDs);
  uInt freqID = freqIDs(IFSel_);
  restfreqidCol_.get(whichRow, freqIDs);
  uInt restFreqID = freqIDs(IFSel_);

// Get SpectralCoordinate, set reference frame conversion,
// velocity conversion, and rest freq state

  SpectralCoordinate spc = getSpectralCoordinate(freqID, restFreqID, whichRow);
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
  Table t = table_.keywordSet().asTable("FREQUENCIES");
//
  String sunit;
  t.keywordSet().get("UNIT",sunit);
  if (sunit == "") sunit = "pixel";
  Unit u(sunit);
//
  Vector<uInt> freqIDs;
  freqidCol_.get(whichRow, freqIDs);
  uInt freqID = freqIDs(IFSel_);  
  restfreqidCol_.get(whichRow, freqIDs);
  uInt restFreqID = freqIDs(IFSel_);

// Get SpectralCoordinate, with frame, velocity, rest freq state set

  SpectralCoordinate spc = getSpectralCoordinate(freqID, restFreqID, whichRow);
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
  Array<Float> arr;
  specCol_.get(whichRow, arr);
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
  specCol_.put(whichRow, arr);
}

void SDMemTable::getSpectrum(Vector<Float>& spectrum, Int whichRow) const
{
  Array<Float> arr;
  specCol_.get(whichRow, arr);

// Iterate and extract

  spectrum.resize(arr.shape()(3));
  ArrayAccessor<Float, Axis<asap::BeamAxis> > aa0(arr);
  aa0.reset(aa0.begin(uInt(beamSel_)));//go to beam
  ArrayAccessor<Float, Axis<asap::IFAxis> > aa1(aa0);
  aa1.reset(aa1.begin(uInt(IFSel_)));// go to IF
  ArrayAccessor<Float, Axis<asap::PolAxis> > aa2(aa1);
  aa2.reset(aa2.begin(uInt(polSel_)));// go to pol
//
  ArrayAccessor<Float, Axis<asap::BeamAxis> > va(spectrum);
  for (ArrayAccessor<Float, Axis<asap::ChanAxis> > i(aa2); i != i.end(); ++i) {
    (*va) = (*i);
    va++;
  }
}


/*
void SDMemTable::getMask(Vector<Bool>& mask, Int whichRow) const {
  Array<uChar> arr;
  flagsCol_.get(whichRow, arr);
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
  Array<Float> arr;
  Array<uChar> farr;
  specCol_.get(whichRow, arr);
  flagsCol_.get(whichRow, farr);
  Array<Bool> barr(farr.shape());
  convertArray(barr, farr);
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
  Array<Float> arr;
  tsCol_.get(whichRow, arr);
  Float out;
//
  IPosition ip(arr.shape());
  ip(asap::BeamAxis) = beamSel_;
  ip(asap::IFAxis) = IFSel_;
  ip(asap::PolAxis) = polSel_;
  ip(asap::ChanAxis)=0;               // First channel only
//
  out = arr(ip);
  return out;
}

MDirection SDMemTable::getDirection(Int whichRow, Bool refBeam) const
{
  MDirection::Types mdr = getDirectionReference();
  Array<Double> posit;
  dirCol_.get(whichRow,posit);
  Vector<Double> wpos(2);
  Int rb;
  rbeamCol_.get(whichRow,rb);
  wpos[0] = posit(IPosition(2,beamSel_,0));
  wpos[1] = posit(IPosition(2,beamSel_,1));
  if (refBeam && rb > -1) {  // use refBeam instead if it exists
    wpos[0] = posit(IPosition(2,rb,0));
    wpos[1] = posit(IPosition(2,rb,1));
  }

  Quantum<Double> lon(wpos[0],Unit(String("rad")));
  Quantum<Double> lat(wpos[1],Unit(String("rad")));
  return MDirection(lon, lat, mdr);
}

MEpoch SDMemTable::getEpoch(Int whichRow) const
{
  MEpoch::Types met = getTimeReference();
//
  Double obstime;
  timeCol_.get(whichRow,obstime);
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


SpectralCoordinate SDMemTable::getSpectralCoordinate(uInt freqID) const
{
  
  Table t = table_.keywordSet().asTable("FREQUENCIES");
  if (freqID> t.nrow() ) {
    throw(AipsError("SDMemTable::getSpectralCoordinate - freqID out of range"));
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
  rpc.get(freqID, rp);
  rvc.get(freqID, rv);
  incc.get(freqID, inc);
//
  SpectralCoordinate spec(mft,rv,inc,rp);
//
  return spec;
}


SpectralCoordinate SDMemTable::getSpectralCoordinate(uInt freqID, uInt restFreqID,
                                                     uInt whichRow) const
{

// Create basic SC

  SpectralCoordinate spec = getSpectralCoordinate (freqID);
//
  Table t = table_.keywordSet().asTable("FREQUENCIES");

// Set rest frequency

  Vector<Double> restFreqIDs;
  t.keywordSet().get("RESTFREQS",restFreqIDs);
  if (restFreqID < restFreqIDs.nelements()) {                  // SHould always be true
    spec.setRestFrequency(restFreqIDs(restFreqID),True);
  }

// Set up frame conversion layer

  String frm;
  t.keywordSet().get("REFFRAME",frm);
  if (frm == "") frm = "TOPO";
  MFrequency::Types mtype;
  if (!MFrequency::getType(mtype, frm)) {
    cerr << "Frequency type unknown assuming TOPO" << endl;       // SHould never happen
    mtype = MFrequency::TOPO;
  }

// Set reference frame conversion  (requires row)

  MDirection dir = getDirection(whichRow);
  MEpoch epoch = getEpoch(whichRow);
  MPosition pos = getAntennaPosition();
//
  if (!spec.setReferenceConversion(mtype,epoch,pos,dir)) {
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
                               uInt freqID) {
  Table t = table_.rwKeywordSet().asTable("FREQUENCIES");
  if (freqID > t.nrow() ) {
    throw(AipsError("SDMemTable::setCoordinate - coord no out of range"));
  }
  ScalarColumn<Double> rpc(t, "REFPIX");
  ScalarColumn<Double> rvc(t, "REFVAL");
  ScalarColumn<Double> incc(t, "INCREMENT");

  rpc.put(freqID, speccord.referencePixel()[0]);
  rvc.put(freqID, speccord.referenceValue()[0]);
  incc.put(freqID, speccord.increment()[0]);

  return True;
}

Int SDMemTable::nCoordinates() const
{
  return table_.keywordSet().asTable("FREQUENCIES").nrow();
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
  const Table& t = table_.keywordSet().asTable("FREQUENCIES");
  SDFrequencyTable sdft;

// Add refpix/refval/incr.  What are the units ? Hz I suppose
// but it's nowhere described...

  Vector<Double> refPix, refVal, incr;
  ScalarColumn<Double> refPixCol(t, "REFPIX");
  ScalarColumn<Double> refValCol(t, "REFVAL");
  ScalarColumn<Double> incrCol(t, "INCREMENT");
  refPix = refPixCol.getColumn();
  refVal = refValCol.getColumn();
  incr = incrCol.getColumn();
//
  uInt n = refPix.nelements();
  for (uInt i=0; i<n; i++) {
     sdft.addFrequency(refPix[i], refVal[i], incr[i]);
  }

// Frequency reference frame.  I don't know if this
// is the correct frame.  It might be 'REFFRAME'
// rather than 'BASEREFFRAME' ?

  String baseFrame;
  t.keywordSet().get("BASEREFFRAME",baseFrame);
  sdft.setRefFrame(baseFrame);

// Equinox

  Float equinox;
  t.keywordSet().get("EQUINOX", equinox);
  sdft.setEquinox(equinox);

// Rest Frequency

  Vector<Double> restFreqs;
  t.keywordSet().get("RESTFREQS", restFreqs);
  for (uInt i=0; i<restFreqs.nelements(); i++) {
     sdft.addRestFrequency(restFreqs[i]);
  }
  sdft.setRestFrequencyUnit(String("Hz"));
//
  return sdft;
}

bool SDMemTable::putSDContainer(const SDContainer& sdc)
{
  uInt rno = table_.nrow();
  table_.addRow();
//
  timeCol_.put(rno, sdc.timestamp);
  srcnCol_.put(rno, sdc.sourcename);
  fldnCol_.put(rno, sdc.fieldname);
  specCol_.put(rno, sdc.getSpectrum());
  flagsCol_.put(rno, sdc.getFlags());
  tsCol_.put(rno, sdc.getTsys());
  scanCol_.put(rno, sdc.scanid);
  integrCol_.put(rno, sdc.interval);
  freqidCol_.put(rno, sdc.getFreqMap());
  restfreqidCol_.put(rno, sdc.getRestFreqMap());
  dirCol_.put(rno, sdc.getDirection());
  rbeamCol_.put(rno, sdc.refbeam);
  tcalCol_.put(rno, sdc.tcal);
  tcaltCol_.put(rno, sdc.tcaltime);
  azCol_.put(rno, sdc.azimuth);
  elCol_.put(rno, sdc.elevation);
  paraCol_.put(rno, sdc.parangle);
  histCol_.put(rno, sdc.getHistory());
//
  return true;
}

SDContainer SDMemTable::getSDContainer(uInt whichRow) const
{
  SDContainer sdc(nBeam(),nIF(),nPol(),nChan());
  timeCol_.get(whichRow, sdc.timestamp);
  srcnCol_.get(whichRow, sdc.sourcename);
  integrCol_.get(whichRow, sdc.interval);
  scanCol_.get(whichRow, sdc.scanid);
  fldnCol_.get(whichRow, sdc.fieldname);
  rbeamCol_.get(whichRow, sdc.refbeam);
  azCol_.get(whichRow, sdc.azimuth);
  elCol_.get(whichRow, sdc.elevation);
  paraCol_.get(whichRow, sdc.parangle);
  Vector<Float> tc;
  tcalCol_.get(whichRow, tc);
  sdc.tcal[0] = tc[0];sdc.tcal[1] = tc[1];
  tcaltCol_.get(whichRow, sdc.tcaltime);
//
  Array<Float> spectrum;
  Array<Float> tsys;
  Array<uChar> flagtrum;
  Vector<uInt> fmap;
  Array<Double> direction;
  Vector<String> histo;
//
  specCol_.get(whichRow, spectrum);
  sdc.putSpectrum(spectrum);
  flagsCol_.get(whichRow, flagtrum);
  sdc.putFlags(flagtrum);
  tsCol_.get(whichRow, tsys);
  sdc.putTsys(tsys);
  freqidCol_.get(whichRow, fmap);
  sdc.putFreqMap(fmap);
  restfreqidCol_.get(whichRow, fmap);
  sdc.putRestFreqMap(fmap);
  dirCol_.get(whichRow, direction);
  sdc.putDirection(direction);
  histCol_.get(whichRow, histo);
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
  Int previous = -1;Int current=0;
  for (uInt i=0; i< scanCol_.nrow();i++) {
    scanCol_.getScalar(i,current);
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

std::string SDMemTable::summary(bool verbose) const  {

// Format header info

  ostringstream oss;
  oss << endl;
  oss << "--------------------------------------------------------------------------------" << endl;
  oss << " Scan Table Summary" << endl;
  oss << "--------------------------------------------------------------------------------" << endl;
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
//
  String dirtype ="Position ("+ MDirection::showType(getDirectionReference()) + ")";
  oss << setw(5) << "Scan"
      << setw(15) << "Source"
      << setw(24) << dirtype
      << setw(10) << "Time"
      << setw(18) << "Integration" 
      << setw(7) << "FreqIDs" << endl;
  oss << "--------------------------------------------------------------------------------" << endl;
  
// Generate list of scan start and end integrations

  Vector<Int> scanIDs = scanCol_.getColumn();
  Vector<uInt> startInt, endInt;
  mathutil::scanBoundaries(startInt, endInt, scanIDs);
//
  const uInt nScans = startInt.nelements();
  String name;
  Vector<uInt> freqIDs, listFQ;
  uInt nInt;
//
  for (uInt i=0; i<nScans; i++) {

// Get things from first integration of scan

      String time = getTime(startInt(i),False);
      String tInt = formatSec(Double(getInterval(startInt(i))));
      String posit = formatDirection(getDirection(startInt(i),True));
      srcnCol_.getScalar(startInt(i),name);

// Find all the FreqIDs in this scan

      listFQ.resize(0);      
      for (uInt j=startInt(i); j<endInt(i)+1; j++) {
         freqidCol_.get(j, freqIDs);
         for (uInt k=0; k<freqIDs.nelements(); k++) {
            mathutil::addEntry(listFQ, freqIDs(k));
         }
      }
//
      nInt = endInt(i) - startInt(i) + 1;
      oss << setw(3) << std::right << i << std::left << setw(2) << "  "
          << setw(15) << name
	  << setw(24) << posit
	  << setw(10) << time 
	  << setw(3) << std::right << nInt  << setw(3) << " x " << std::left
	  << setw(6) <<  tInt
          << " " << listFQ << endl;
  }
  oss << endl;
  oss << "Table contains " << table_.nrow() << " integration(s) in " << nScans << " scan(s)." << endl;

// Frequency Table

  if (verbose) {
    std::vector<string> info = getCoordInfo();
    SDFrequencyTable sdft = getSDFreqTable();
    oss << endl << endl;
    oss << "FreqID  Frame   RefFreq(Hz)     RefPix   Increment(Hz)" << endl;
    oss << "--------------------------------------------------------------------------------" << endl;
    for (uInt i=0; i<sdft.length(); i++) {
      oss << setw(8) << i << setw(8)
	  << info[3] << setw(16) << setprecision(8)
	  << sdft.referenceValue(i) << setw(10)
	  << sdft.referencePixel(i) << setw(12)
	  << sdft.increment(i) << endl;
    }
    oss << "--------------------------------------------------------------------------------" << endl;
  }
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
  Vector<String> history;
  histCol_.get(whichRow, history);
  history.resize(history.nelements()+1,True);
  history[history.nelements()-1] = hist;
  histCol_.put(whichRow, history);
}

std::vector<std::string> SDMemTable::history(int whichRow) const
{
  Vector<String> history;
  histCol_.get(whichRow, history);
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
  Array<uChar> arr;
  flagsCol_.get(whichRow, arr);

  ArrayAccessor<uChar, Axis<asap::BeamAxis> > aa0(arr);
  aa0.reset(aa0.begin(uInt(beamSel_)));//go to beam
  ArrayAccessor<uChar, Axis<asap::IFAxis> > aa1(aa0);
  aa1.reset(aa1.begin(uInt(IFSel_)));// go to IF
  ArrayAccessor<uChar, Axis<asap::PolAxis> > aa2(aa1);
  aa2.reset(aa2.begin(uInt(polSel_)));// go to pol

  for (ArrayAccessor<uChar, Axis<asap::ChanAxis> > i(aa2); i != i.end(); ++i) {
    (*i) = uChar(True);
  }

  flagsCol_.put(whichRow, arr);
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

MEpoch::Types SDMemTable::getTimeReference() const
{
  MEpoch::Types met;
  String ep;
  table_.keywordSet().get("Epoch",ep);
  if (!MEpoch::getType(met, ep)) {
    cerr << "Epoch type unknown - using UTC" << endl;
    met = MEpoch::UTC;
  }
//
  return met;
}


Instrument SDMemTable::convertInstrument(const String& instrument,
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

Bool SDMemTable::setRestFreqs (const Vector<Double>& restFreqsIn, 
                                     const String& sUnit,
                                     const vector<string>& lines,
                                     const String& source,
                                     Int whichIF)
{
   const Int nIFs = nIF();
   if (whichIF>=nIFs) {
      throw(AipsError("Illegal IF index"));
   }

// FInd vector of restfrequencies
// Double takes precedence over String

   Unit unit;
   Vector<Double> restFreqs;
   if (restFreqsIn.nelements()>0) {
      restFreqs.resize(restFreqsIn.nelements());
      restFreqs = restFreqsIn;
      unit = Unit(sUnit);
   } else if (lines.size()>0) {
      const uInt nLines = lines.size();
      unit = Unit("Hz");
      restFreqs.resize(nLines);
      MFrequency lineFreq;
      for (uInt i=0; i<nLines; i++) {
         String tS(lines[i]);
         tS.upcase();
         if (MeasTable::Line(lineFreq, tS)) {
            restFreqs[i] = lineFreq.getValue().getValue();          // Hz
         } else {
            String s = String(lines[i]) + String(" is an unrecognized spectral line");
            throw(AipsError(s));
         }
      }
   } else {
      throw(AipsError("You have not specified any rest frequencies or lines"));
   }

// If multiple restfreqs, must be length nIF. In this
// case we will just replace the rest frequencies
// 

   const uInt nRestFreqs = restFreqs.nelements();
   Int idx = -1;
   SDFrequencyTable sdft = getSDFreqTable();

   if (nRestFreqs>1) {

// Replace restFreqs, one per IF

      if (nRestFreqs != nIFs) {
         throw (AipsError("Number of rest frequencies must be equal to the number of IFs"));
      }
      cout << "Replacing rest frequencies with given list, one per IF" << endl;
//
      sdft.deleteRestFrequencies();
      for (uInt i=0; i<nRestFreqs; i++) {
         Quantum<Double> rf(restFreqs[i], unit);
         sdft.addRestFrequency(rf.getValue("Hz"));
      }
   } else {

// Add new rest freq

      Quantum<Double> rf(restFreqs[0], unit);
      idx = sdft.addRestFrequency(rf.getValue("Hz"));
      cout << "Selecting given rest frequency" << endl;
   }

// Replace

   Bool empty = source.empty();
   Bool ok = False;
   if (putSDFreqTable(sdft)) {
      const uInt nRow = table_.nrow();
      String srcName;
      Vector<uInt> restFreqIDs;
      for (uInt i=0; i<nRow; i++) {
         srcnCol_.get(i, srcName);
         restfreqidCol_.get(i,restFreqIDs);       
// 
         if (idx==-1) {

// Replace vector of restFreqs; one per IF. 
// No selection possible

            for (uInt i=0; i<nIFs; i++) restFreqIDs[i] = i;
         } else {

// Set RestFreqID for selected data

            if (empty || source==srcName) {
               if (whichIF<0) {
                  restFreqIDs = idx;
               } else {              
                  restFreqIDs[whichIF] = idx;
               }
            }
         }
//
         restfreqidCol_.put(i,restFreqIDs);       
      }
      ok = True;
   } else {
      ok = False;
   }
//
   return ok;
}

void SDMemTable::spectralLines() const
{
   Vector<String> lines = MeasTable::Lines();
   MFrequency lineFreq;
   Double freq;
//
   cout.flags(std::ios_base::left);
   cout << "Line      Frequency (Hz)" << endl;
   cout << "-----------------------" << endl;
   for (uInt i=0; i<lines.nelements(); i++) {
     MeasTable::Line(lineFreq, lines[i]);
     freq = lineFreq.getValue().getValue();          // Hz
//
     cout << setw(11) << lines[i] << setprecision(10) << freq << endl;
   }
}

void SDMemTable::renumber()
{
  uInt nRow = scanCol_.nrow();
  Int newscanid = 0;
  Int cIdx;// the current scanid
  // get the first scanid
  scanCol_.getScalar(0,cIdx);
  Int pIdx = cIdx;// the scanid of the previous row
  for (uInt i=0; i<nRow;++i) {
    scanCol_.getScalar(i,cIdx);
    if (pIdx == cIdx) {
      // renumber
      scanCol_.put(i,newscanid);
    } else { 
      ++newscanid;
      pIdx = cIdx;   // store scanid
      --i;           // don't increment next loop
    }
  }
}


void SDMemTable::rotateXYPhase (Float value) 
//
// phase in degrees
//
{
   if (nPol() != 4) {
      throw(AipsError("You must have 4 polarizations to run this function"));
   }
//
   IPosition start(asap::nAxes,0);
   IPosition end(asap::nAxes);
//
   uInt nRow = specCol_.nrow();
   Array<Float> data;
   for (uInt i=0; i<nRow;++i) {
      specCol_.get(i,data);
      end = data.shape()-1;

// Get polarization slice references

      start(asap::PolAxis) = 2;
      end(asap::PolAxis) = 2;
      Array<Float> C3 = data(start,end);
//
      start(asap::PolAxis) = 3;
      end(asap::PolAxis) = 3;
      Array<Float> C4 = data(start,end);

// Rotate

      SDPolUtil::rotateXYPhase(C3, C4, value);

// Put

      specCol_.put(i,data);
   }
}
