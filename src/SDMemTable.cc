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
#include <casa/Quanta/MVTime.h>

#include "SDMemTable.h"
#include "SDContainer.h"

using namespace atnf_sd;

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
  //cerr << exprs << endl;
  Table t = tableCommand(exprs,tab);
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
  mvt.setFormat(MVTime::YMD);
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
  chanMask_.resize(n,true);
  for (it = whichChans.begin(); it != whichChans.end(); ++it) {
    if (*it < n)
      chanMask_[*it] = false;
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

std::vector<double> SDMemTable::getAbscissa(Int whichRow, 
					    const std::string& whichUnit,
					    const std::string& whichFrame,
					    double restfreq) {
  std::vector<double> absc(nChan());
  Vector<Double> absc1(nChan()); 
  indgen(absc1);
  ROArrayColumn<uInt> fid(table_, "FREQID");
  Vector<uInt> v;
  fid.get(whichRow, v);
  uInt specidx = v(IFSel_);
  Unit u;
  if (whichUnit == "") {
    // get unit from table
  } else {
    u = String(whichUnit);
  }
  SpectralCoordinate spc = getCoordinate(specidx);
  Table t = table_.keywordSet().asTable("FREQUENCIES");
  String rf;
  t.keywordSet().get("REFFRAME",rf);
  MDirection::Types mdr;
  MDirection::getType(mdr, rf);    
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
  Quantum<Double> tm(obstime, Unit(String("d")));
  MVEpoch tm2(tm);
  MEpoch epoch(tm2);

  Vector<Double> antpos;
  table_.keywordSet().get("AntennaPosition", antpos);
  //MPosition pos;
  MVPosition mvpos(antpos(0),antpos(1),antpos(2));
  MPosition pos(mvpos);
  //MeasTable::Observatory(pos, String("ATCA"));

  MFrequency::Types mtype;
  if (!MFrequency::getType(mtype, whichFrame)) {
    cerr << "Frequency type unknown assuming TOPO" << endl;
    mtype = MFrequency::TOPO;
  }
  spc.setReferenceConversion(mtype,epoch,pos,direct);  

  if ( u == Unit("km/s") ) {
    if (Double(restfreq) >  Double(0.000001)) {
      cerr << "converting to velocities"<< endl;
      spc.setRestFrequency(Double(restfreq));
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
    cerr << " converting in frequency" << endl;
    std::vector<double>::iterator it;
    Double tmp;
    uInt i = 0;
    for (it = absc.begin(); it != absc.end(); ++it) {      
      spc.toWorld(tmp,absc1[i]);
      (*it) = tmp;
      i++;
    }
  }
  return absc;
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
  ROScalarColumn<Double> rpc(t, "REFPIX");
  ROScalarColumn<Double> rvc(t, "REFVAL");
  ROScalarColumn<Double> incc(t, "INCREMENT");
  t.keywordSet().get("REFFRAME",rf);
  
  MFrequency::Types mft;
  if (!MFrequency::getType(mft, rf)) {
    cerr << "Frequency type unknown assuming TOPO" << endl;
    mft = MFrequency::TOPO;    
  }    
  rpc.get(whichIdx, rp);
  rvc.get(whichIdx, rv);
  incc.get(whichIdx, inc);
  //cerr << "creating speccord from " << whichIdx << ": "
  //     << rp <<", " << rv << ", " << inc << ", " << mft <<endl;
  SpectralCoordinate spec(mft,rv,inc,rp);
  return spec;
}

Bool SDMemTable::setCoordinate(const SpectralCoordinate& speccord, 
			       uInt whichIdx) {
  Table t = table_.rwKeywordSet().asTable("FREQUENCIES");
  if (whichIdx > t.nrow() ) {
    cerr << "SDMemTable::setCoordinate - whichIdx out of range" << endl;
    return False;
  }
  ScalarColumn<Double> rpc(t, "REFPIX");
  ScalarColumn<Double> rvc(t, "REFVAL");
  ScalarColumn<Double> incc(t, "INCREMENT");
  
  rpc.put(whichIdx, speccord.referencePixel()[0]);
  rvc.put(whichIdx, speccord.referenceValue()[0]);
  incc.put(whichIdx, speccord.increment()[0]);

  return True;
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
  aTable.rwKeywordSet().define("REFFRAME", sdft.refFrame());
  aTable.rwKeywordSet().define("EQUINOX", sdft.equinox());
  aTable.rwKeywordSet().define("Unit", String("kms-1"));
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
  ArrayColumn<Float> spec(table_, "SPECTRA");
  ArrayColumn<uChar> flags(table_, "FLAGTRA");
  ArrayColumn<Float> ts(table_, "TSYS");
  ScalarColumn<Int> scan(table_, "SCANID");
  ScalarColumn<Double> integr(table_, "INTERVAL");
  ArrayColumn<uInt> freqid(table_, "FREQID");
  ArrayColumn<Double> dir(table_, "DIRECTION");

  uInt rno = table_.nrow();
  table_.addRow();
  
  mjd.put(rno, sdc.timestamp);
  srcn.put(rno, sdc.sourcename);
  spec.put(rno, sdc.getSpectrum());
  flags.put(rno, sdc.getFlags());
  ts.put(rno, sdc.getTsys());
  scan.put(rno, sdc.scanid);
  integr.put(rno, sdc.interval);
  freqid.put(rno, sdc.getFreqMap());
  dir.put(rno, sdc.getDirection());
  
  return true;
}

SDContainer SDMemTable::getSDContainer(uInt whichRow) const {
  ROScalarColumn<Double> mjd(table_, "TIME");
  ROScalarColumn<String> srcn(table_, "SRCNAME");
  ROArrayColumn<Float> spec(table_, "SPECTRA");
  ROArrayColumn<uChar> flags(table_, "FLAGTRA");
  ROArrayColumn<Float> ts(table_, "TSYS");
  ROScalarColumn<Int> scan(table_, "SCANID");
  ROScalarColumn<Double> integr(table_, "INTERVAL");
  ROArrayColumn<uInt> freqid(table_, "FREQID");
  ROArrayColumn<Double> dir(table_, "DIRECTION");

  SDContainer sdc(nBeam(),nIF(),nPol(),nChan());
  mjd.get(whichRow, sdc.timestamp);
  srcn.get(whichRow, sdc.sourcename);
  integr.get(whichRow, sdc.interval);
  scan.get(whichRow, sdc.scanid);
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

Int SDMemTable::nScans() const {
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

void SDMemTable::summary() const {
  ROScalarColumn<Int> scans(table_, "SCANID");
  ROScalarColumn<String> srcs(table_, "SRCNAME");
  cout << "*************** Header ***************" << endl;
  cout << "nBeam = " << nBeam() << "\t"
       << "nIF   = " << nIF() << endl
       << "nPol  = " << nPol() << "\t"
       << "nChan = " << nChan() << "\t" << endl;
  cout << "*************** Header ***************" << endl;
  uInt count = 0;
  String name;
  Int previous = -1;Int current=0;
  cout << "Scan\tSource\t\tTime\t\tIntegration" << endl;
  for (uInt i=0; i< scans.nrow();i++) {
    scans.getScalar(i,current);
    if (previous != current) {
      srcs.getScalar(i,name);
      previous = current;
      Double t = getInterval();
      String unit("sec");
      if (t/60.0 > 1.0) {
	t/=60.0;unit = "min";
      }
      cout << count << "\t" 
	   << name << "\t"
	   << getTime() << "\t"
	   << setprecision(2) <<  setiosflags(std::ios_base::fixed)
	   << t << " " << unit << endl
	   << endl;
      count++;
    }
  }
  cout << "Table contains " << table_.nrow() << " integration(s)." << endl;
  cout << "in " << count << " scan(s)." << endl;
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
