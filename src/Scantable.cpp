//
// C++ Implementation: Scantable
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <map>

#include <casa/aips.h>
#include <casa/iostream.h>
#include <casa/iomanip.h>
#include <casa/OS/Path.h>
#include <casa/OS/File.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/MaskArrMath.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Arrays/ArrayAccessor.h>
#include <casa/Arrays/VectorSTLIterator.h>
#include <casa/Arrays/Vector.h>
#include <casa/BasicMath/Math.h>
#include <casa/BasicSL/Constants.h>
#include <casa/Quanta/MVAngle.h>
#include <casa/Containers/RecordField.h>

#include <tables/Tables/TableParse.h>
#include <tables/Tables/TableDesc.h>
#include <tables/Tables/TableCopy.h>
#include <tables/Tables/SetupNewTab.h>
#include <tables/Tables/ScaColDesc.h>
#include <tables/Tables/ArrColDesc.h>
#include <tables/Tables/TableRow.h>
#include <tables/Tables/TableVector.h>
#include <tables/Tables/TableIter.h>

#include <tables/Tables/ExprNode.h>
#include <tables/Tables/TableRecord.h>
#include <measures/Measures/MFrequency.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MeasTable.h>
#include <measures/Measures/MeasRef.h>
#include <measures/TableMeasures/TableMeasRefDesc.h>
#include <measures/TableMeasures/TableMeasValueDesc.h>
#include <measures/TableMeasures/TableMeasDesc.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>
#include <coordinates/Coordinates/CoordinateUtil.h>
#include <casa/Quanta/MVTime.h>
#include <casa/Quanta/MVAngle.h>

#include "Scantable.h"
#include "STPolLinear.h"
#include "STAttr.h"

using namespace casa;

namespace asap {

std::map<std::string, STPol::STPolFactory *> Scantable::factories_;

void Scantable::initFactories() {
  if ( factories_.empty() ) {
    Scantable::factories_["linear"] = &STPolLinear::myFactory;
  }
}

Scantable::Scantable(Table::TableType ttype) :
  type_(ttype)
{
  initFactories();
  setupMainTable();
  freqTable_ = STFrequencies(*this);
  table_.rwKeywordSet().defineTable("FREQUENCIES", freqTable_.table());
  weatherTable_ = STWeather(*this);
  table_.rwKeywordSet().defineTable("WEATHER", weatherTable_.table());
  focusTable_ = STFocus(*this);
  table_.rwKeywordSet().defineTable("FOCUS", focusTable_.table());
  tcalTable_ = STTcal(*this);
  table_.rwKeywordSet().defineTable("TCAL", tcalTable_.table());
  moleculeTable_ = STMolecules(*this);
  table_.rwKeywordSet().defineTable("MOLECULES", moleculeTable_.table());
  historyTable_ = STHistory(*this);
  table_.rwKeywordSet().defineTable("HISTORY", historyTable_.table());
  setupFitTable();
  fitTable_ = table_.keywordSet().asTable("FITS");
  originalTable_ = table_;
  attach();
}

Scantable::Scantable(const std::string& name, Table::TableType ttype) :
  type_(ttype)
{
  initFactories();
  Table tab(name, Table::Update);
  Int version;
  tab.keywordSet().get("VERSION", version);
  if (version != version_) {
    throw(AipsError("Unsupported version of ASAP file."));
  }
  if ( type_ == Table::Memory )
    table_ = tab.copyToMemoryTable(generateName());
  else
    table_ = tab;
  attachSubtables();
  originalTable_ = table_;
  attach();
}

Scantable::Scantable( const Scantable& other, bool clear )
{
  // with or without data
  String newname = String(generateName());
  type_ = other.table_.tableType();
  if ( other.table_.tableType() == Table::Memory ) {
      if ( clear ) {
        table_ = TableCopy::makeEmptyMemoryTable(newname,
                                                 other.table_, True);
      } else
        table_ = other.table_.copyToMemoryTable(newname);
  } else {
      other.table_.deepCopy(newname, Table::New, False, Table::AipsrcEndian,
                            Bool(clear));
      table_ = Table(newname, Table::Update);
      copySubtables(other);
      table_.markForDelete();
  }

  attachSubtables();
  originalTable_ = table_;
  attach();
}

void Scantable::copySubtables(const Scantable& other) {
  Table t = table_.rwKeywordSet().asTable("FREQUENCIES");
  TableCopy::copyRows(t, other.freqTable_.table());
  t = table_.rwKeywordSet().asTable("FOCUS");
  TableCopy::copyRows(t, other.focusTable_.table());
  t = table_.rwKeywordSet().asTable("WEATHER");
  TableCopy::copyRows(t, other.weatherTable_.table());
  t = table_.rwKeywordSet().asTable("TCAL");
  TableCopy::copyRows(t, other.tcalTable_.table());
  t = table_.rwKeywordSet().asTable("MOLECULES");
  TableCopy::copyRows(t, other.moleculeTable_.table());
  t = table_.rwKeywordSet().asTable("HISTORY");
  TableCopy::copyRows(t, other.historyTable_.table());
}

void Scantable::attachSubtables()
{
  freqTable_ = STFrequencies(table_);
  focusTable_ = STFocus(table_);
  weatherTable_ = STWeather(table_);
  tcalTable_ = STTcal(table_);
  moleculeTable_ = STMolecules(table_);
  historyTable_ = STHistory(table_);
}

Scantable::~Scantable()
{
  cout << "~Scantable() " << this << endl;
}

void Scantable::setupMainTable()
{
  TableDesc td("", "1", TableDesc::Scratch);
  td.comment() = "An ASAP Scantable";
  td.rwKeywordSet().define("VERSION", Int(version_));

  // n Cycles
  td.addColumn(ScalarColumnDesc<uInt>("SCANNO"));
  // new index every nBeam x nIF x nPol
  td.addColumn(ScalarColumnDesc<uInt>("CYCLENO"));

  td.addColumn(ScalarColumnDesc<uInt>("BEAMNO"));
  td.addColumn(ScalarColumnDesc<uInt>("IFNO"));
  td.rwKeywordSet().define("POLTYPE", String("linear"));
  td.addColumn(ScalarColumnDesc<uInt>("POLNO"));

  td.addColumn(ScalarColumnDesc<uInt>("FREQ_ID"));
  td.addColumn(ScalarColumnDesc<uInt>("MOLECULE_ID"));
  // linear, circular, stokes [I Q U V], stokes1 [I Plinear Pangle V]
  td.addColumn(ScalarColumnDesc<Int>("REFBEAMNO"));

  td.addColumn(ScalarColumnDesc<Double>("TIME"));
  TableMeasRefDesc measRef(MEpoch::UTC); // UTC as default
  TableMeasValueDesc measVal(td, "TIME");
  TableMeasDesc<MEpoch> mepochCol(measVal, measRef);
  mepochCol.write(td);

  td.addColumn(ScalarColumnDesc<Double>("INTERVAL"));

  td.addColumn(ScalarColumnDesc<String>("SRCNAME"));
  // Type of source (on=0, off=1, other=-1)
  td.addColumn(ScalarColumnDesc<Int>("SRCTYPE", Int(-1)));
  td.addColumn(ScalarColumnDesc<String>("FIELDNAME"));

  //The actual Data Vectors
  td.addColumn(ArrayColumnDesc<Float>("SPECTRA"));
  td.addColumn(ArrayColumnDesc<uChar>("FLAGTRA"));
  td.addColumn(ArrayColumnDesc<Float>("TSYS"));

  td.addColumn(ArrayColumnDesc<Double>("DIRECTION",
                                       IPosition(1,2),
                                       ColumnDesc::Direct));
  TableMeasRefDesc mdirRef(MDirection::J2000); // default
  TableMeasValueDesc tmvdMDir(td, "DIRECTION");
  // the TableMeasDesc gives the column a type
  TableMeasDesc<MDirection> mdirCol(tmvdMDir, mdirRef);
  // writing create the measure column
  mdirCol.write(td);
  td.addColumn(ScalarColumnDesc<Double>("AZIMUTH"));
  td.addColumn(ScalarColumnDesc<Double>("ELEVATION"));
  td.addColumn(ScalarColumnDesc<Float>("PARANGLE"));

  td.addColumn(ScalarColumnDesc<uInt>("TCAL_ID"));
  td.addColumn(ScalarColumnDesc<uInt>("FIT_ID"));

  td.addColumn(ScalarColumnDesc<uInt>("FOCUS_ID"));
  td.addColumn(ScalarColumnDesc<uInt>("WEATHER_ID"));

  td.rwKeywordSet().define("OBSMODE", String(""));

  // Now create Table SetUp from the description.
  SetupNewTable aNewTab(generateName(), td, Table::Scratch);
  table_ = Table(aNewTab, type_, 0);
  originalTable_ = table_;

}

void Scantable::setupFitTable()
{
  TableDesc td("", "1", TableDesc::Scratch);
  td.addColumn(ScalarColumnDesc<uInt>("FIT_ID"));
  td.addColumn(ArrayColumnDesc<String>("FUNCTIONS"));
  td.addColumn(ArrayColumnDesc<Int>("COMPONENTS"));
  td.addColumn(ArrayColumnDesc<Double>("PARAMETERS"));
  td.addColumn(ArrayColumnDesc<Bool>("PARMASK"));
  td.addColumn(ArrayColumnDesc<String>("FRAMEINFO"));
  SetupNewTable aNewTab("fits", td, Table::Scratch);
  Table aTable(aNewTab, Table::Memory);
  table_.rwKeywordSet().defineTable("FITS", aTable);
}

void Scantable::attach()
{
  timeCol_.attach(table_, "TIME");
  srcnCol_.attach(table_, "SRCNAME");
  specCol_.attach(table_, "SPECTRA");
  flagsCol_.attach(table_, "FLAGTRA");
  tsysCol_.attach(table_, "TSYS");
  cycleCol_.attach(table_,"CYCLENO");
  scanCol_.attach(table_, "SCANNO");
  beamCol_.attach(table_, "BEAMNO");
  ifCol_.attach(table_, "IFNO");
  polCol_.attach(table_, "POLNO");
  integrCol_.attach(table_, "INTERVAL");
  azCol_.attach(table_, "AZIMUTH");
  elCol_.attach(table_, "ELEVATION");
  dirCol_.attach(table_, "DIRECTION");
  paraCol_.attach(table_, "PARANGLE");
  fldnCol_.attach(table_, "FIELDNAME");
  rbeamCol_.attach(table_, "REFBEAMNO");

  mfitidCol_.attach(table_,"FIT_ID");
  //fitidCol_.attach(fitTable_,"FIT_ID");

  mfreqidCol_.attach(table_, "FREQ_ID");

  mtcalidCol_.attach(table_, "TCAL_ID");

  mfocusidCol_.attach(table_, "FOCUS_ID");

  mmolidCol_.attach(table_, "MOLECULE_ID");
}

void Scantable::setHeader(const SDHeader& sdh)
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
}

SDHeader Scantable::getHeader() const
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

bool Scantable::conformant( const Scantable& other )
{
  return this->getHeader().conformant(other.getHeader());
}


int Scantable::nscan() const {
  int n = 0;
  int previous = -1; int current = 0;
  for (uInt i=0; i< scanCol_.nrow();i++) {
    scanCol_.getScalar(i,current);
    if (previous != current) {
      previous = current;
      n++;
    }
  }
  return n;
}

std::string Scantable::formatSec(Double x) const
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

std::string Scantable::formatDirection(const MDirection& md) const
{
  Vector<Double> t = md.getAngle(Unit(String("rad"))).getValue();
  Int prec = 7;

  MVAngle mvLon(t[0]);
  String sLon = mvLon.string(MVAngle::TIME,prec);
  MVAngle mvLat(t[1]);
  String sLat = mvLat.string(MVAngle::ANGLE+MVAngle::DIG2,prec);
  return sLon + String(" ") + sLat;
}


std::string Scantable::getFluxUnit() const
{
  return table_.keywordSet().asString("FluxUnit");
}

void Scantable::setFluxUnit(const std::string& unit)
{
  String tmp(unit);
  Unit tU(tmp);
  if (tU==Unit("K") || tU==Unit("Jy")) {
     table_.rwKeywordSet().define(String("FluxUnit"), tmp);
  } else {
     throw AipsError("Illegal unit - must be compatible with Jy or K");
  }
}

void Scantable::setInstrument(const std::string& name)
{
  bool throwIt = true;
  Instrument ins = STAttr::convertInstrument(name, throwIt);
  String nameU(name);
  nameU.upcase();
  table_.rwKeywordSet().define(String("AntennaName"), nameU);
}

MPosition Scantable::getAntennaPosition () const
{
  Vector<Double> antpos;
  table_.keywordSet().get("AntennaPosition", antpos);
  MVPosition mvpos(antpos(0),antpos(1),antpos(2));
  return MPosition(mvpos);
}

void Scantable::makePersistent(const std::string& filename)
{
  String inname(filename);
  Path path(inname);
  inname = path.expandedName();
  table_.deepCopy(inname, Table::New);
}

int Scantable::nbeam( int scanno ) const
{
  if ( scanno < 0 ) {
    Int n;
    table_.keywordSet().get("nBeam",n);
    return int(n);
  } else {
    // take the first POLNO,IFNO,CYCLENO as nbeam shouldn't vary with these
    Table t = table_(table_.col("SCANNO") == scanno);
    ROTableRow row(t);
    const TableRecord& rec = row.get(0);
    Table subt = t( t.col("IFNO") == Int(rec.asuInt("IFNO"))
                    && t.col("POLNO") == Int(rec.asuInt("POLNO"))
                    && t.col("CYCLENO") == Int(rec.asuInt("CYCLENO")) );
    ROTableVector<uInt> v(subt, "BEAMNO");
    return int(v.nelements());
  }
  return 0;
}

int Scantable::nif( int scanno ) const
{
  if ( scanno < 0 ) {
    Int n;
    table_.keywordSet().get("nIF",n);
    return int(n);
  } else {
    // take the first POLNO,BEAMNO,CYCLENO as nbeam shouldn't vary with these
    Table t = table_(table_.col("SCANNO") == scanno);
    ROTableRow row(t);
    const TableRecord& rec = row.get(0);
    Table subt = t( t.col("BEAMNO") == Int(rec.asuInt("BEAMNO"))
                    && t.col("POLNO") == Int(rec.asuInt("POLNO"))
                    && t.col("CYCLENO") == Int(rec.asuInt("CYCLENO")) );
    if ( subt.nrow() == 0 ) return 0;
    ROTableVector<uInt> v(subt, "IFNO");
    return int(v.nelements());
  }
  return 0;
}

int Scantable::npol( int scanno ) const
{
  if ( scanno < 0 ) {
    Int n;
    table_.keywordSet().get("nPol",n);
    return n;
  } else {
    // take the first POLNO,IFNO,CYCLENO as nbeam shouldn't vary with these
    Table t = table_(table_.col("SCANNO") == scanno);
    ROTableRow row(t);
    const TableRecord& rec = row.get(0);
    Table subt = t( t.col("BEAMNO") == Int(rec.asuInt("BEAMNO"))
                    && t.col("IFNO") == Int(rec.asuInt("IFNO"))
                    && t.col("CYCLENO") == Int(rec.asuInt("CYCLENO")) );
    if ( subt.nrow() == 0 ) return 0;
    ROTableVector<uInt> v(subt, "POLNO");
    return int(v.nelements());
  }
  return 0;
}

int Scantable::ncycle( int scanno ) const
{
  if ( scanno < 0 ) {
    Block<String> cols(2);
    cols[0] = "SCANNO";
    cols[1] = "CYCLENO";
    TableIterator it(table_, cols);
    int n = 0;
    while ( !it.pastEnd() ) {
      ++n;
    }
    return n;
  } else {
    Table t = table_(table_.col("SCANNO") == scanno);
    ROTableRow row(t);
    const TableRecord& rec = row.get(0);
    Table subt = t( t.col("BEAMNO") == Int(rec.asuInt("BEAMNO"))
                    && t.col("POLNO") == Int(rec.asuInt("POLNO"))
                    && t.col("IFNO") == Int(rec.asuInt("IFNO")) );
    if ( subt.nrow() == 0 ) return 0;
    return int(subt.nrow());
  }
  return 0;
}


int Scantable::nrow( int scanno ) const
{
  return int(table_.nrow());
}

int Scantable::nchan( int ifno ) const
{
  if ( ifno < 0 ) {
    Int n;
    table_.keywordSet().get("nChan",n);
    return int(n);
  } else {
    // take the first SCANNO,POLNO,BEAMNO,CYCLENO as nbeam shouldn't vary with these
    Table t = table_(table_.col("IFNO") == ifno);
    if ( t.nrow() == 0 ) return 0;
    ROArrayColumn<Float> v(t, "SPECTRA");
    return v(0).nelements();
  }
  return 0;
}


int Scantable::getBeam(int whichrow) const
{
  return beamCol_(whichrow);
}

int Scantable::getIF(int whichrow) const
{
  return ifCol_(whichrow);
}

int Scantable::getPol(int whichrow) const
{
  return polCol_(whichrow);
}

std::string Scantable::formatTime(const MEpoch& me, bool showdate) const
{
  MVTime mvt(me.getValue());
  if (showdate)
    mvt.setFormat(MVTime::YMD);
  else
    mvt.setFormat(MVTime::TIME);
  ostringstream oss;
  oss << mvt;
  return String(oss);
}

void Scantable::calculateAZEL()
{
  MPosition mp = getAntennaPosition();
  MEpoch::ROScalarColumn timeCol(table_, "TIME");
  ostringstream oss;
  oss << "Computed azimuth/elevation using " << endl
      << mp << endl;
  for (uInt i=0; i<nrow(); ++i) {
    MEpoch me = timeCol(i);
    MDirection md = dirCol_(i);
    dirCol_.get(i,md);
    oss  << " Time: " << formatTime(me,False) << " Direction: " << formatDirection(md)
         << endl << "     => ";
    MeasFrame frame(mp, me);
    Vector<Double> azel =
        MDirection::Convert(md, MDirection::Ref(MDirection::AZEL,
                                                frame)
                            )().getAngle("rad").getValue();
    azCol_.put(i,azel[0]);
    elCol_.put(i,azel[1]);
    oss << "azel: " << azel[0]/C::pi*180.0 << " "
        << azel[1]/C::pi*180.0 << " (deg)" << endl;
  }
  pushLog(String(oss));
}

void Scantable::flag()
{
  if ( selector_.empty() )
    throw(AipsError("Trying to flag whole scantable. Aborted."));
  TableVector<uChar> tvec(table_, "FLAGTRA");
  uChar userflag = 1 << 7;
  tvec = userflag;
}

std::vector<bool> Scantable::getMask(int whichrow) const
{
  Vector<uChar> flags;
  flagsCol_.get(uInt(whichrow), flags);
  Vector<Bool> bflag(flags.shape());
  convertArray(bflag, flags);
  bflag = !bflag;
  std::vector<bool> mask;
  bflag.tovector(mask);
  return mask;
}

std::vector<float> Scantable::getSpectrum( int whichrow,
                                           const std::string& poltype) const
{
  std::vector<float> out;
  Vector<Float> arr;
  uInt requestedpol = polCol_(whichrow);
  String basetype = getPolType();
  if ( String(poltype) == basetype) {
    specCol_.get(whichrow, arr);
  } else {
    STPol* stpol = 0;
    stpol =STPol::getPolClass(Scantable::factories_, basetype);
    try {
      uInt row = uInt(whichrow);
      stpol->setSpectra(getPolMatrix(row));
      Float frot,fang,ftan;
      focusTable_.getEntry(frot, fang, ftan, row);
      stpol->setPhaseCorrections(frot, fang, ftan);
      arr = stpol->getSpectrum(requestedpol, poltype);
      delete stpol;
    } catch (AipsError& e) {
      delete stpol;
      throw(e);
    }
  }
  arr.tovector(out);
  return out;
}

void asap::Scantable::setSpectrum( const std::vector<float>& spec,
                                   int whichrow )
{
  Vector<Float> spectrum(spec);
  Vector<Float> arr;
  specCol_.get(whichrow, arr);
  if ( spectrum.nelements() != arr.nelements() )
    throw AipsError("The spectrum has incorrect number of channels.");
  specCol_.put(whichrow, spectrum);
}


String Scantable::generateName()
{
  return (File::newUniqueName("./","temp")).baseName();
}

const casa::Table& Scantable::table( ) const
{
  return table_;
}

casa::Table& Scantable::table( )
{
  return table_;
}

std::string Scantable::getPolType() const
{
  return table_.keywordSet().asString("POLTYPE");
}


std::string Scantable::getPolarizationLabel(bool linear, bool stokes,
                                            bool linpol, int polidx) const
{
  uInt idx = 0;
  if (polidx >=0) idx = polidx;
  return "";
  //return SDPolUtil::polarizationLabel(idx, linear, stokes, linpol);
}

void Scantable::unsetSelection()
{
  table_ = originalTable_;
  attach();
  selector_.reset();
}

void Scantable::setSelection( const STSelector& selection )
{
  Table tab = const_cast<STSelector&>(selection).apply(originalTable_);
  if ( tab.nrow() == 0 ) {
    throw(AipsError("Selection contains no data. Not applying it."));
  }
  table_ = tab;
  attach();
  selector_ = selection;
}

std::string Scantable::summary( bool verbose )
{
  // Format header info
  ostringstream oss;
  oss << endl;
  oss << asap::SEPERATOR << endl;
  oss << " Scan Table Summary" << endl;
  oss << asap::SEPERATOR << endl;
  oss.flags(std::ios_base::left);
  oss << setw(15) << "Beams:" << setw(4) << nbeam() << endl
      << setw(15) << "IFs:" << setw(4) << nif() << endl
      << setw(15) << "Polarisations:" << setw(4) << npol()
      << "(" << getPolType() << ")" << endl
      << setw(15) << "Channels:"  << setw(4) << nchan() << endl;
  oss << endl;
  String tmp;
  oss << setw(15) << "Observer:"
      << table_.keywordSet().asString("Observer") << endl;
  oss << setw(15) << "Obs Date:" << getTime(-1,true) << endl;
  table_.keywordSet().get("Project", tmp);
  oss << setw(15) << "Project:" << tmp << endl;
  table_.keywordSet().get("Obstype", tmp);
  oss << setw(15) << "Obs. Type:" << tmp << endl;
  table_.keywordSet().get("AntennaName", tmp);
  oss << setw(15) << "Antenna Name:" << tmp << endl;
  table_.keywordSet().get("FluxUnit", tmp);
  oss << setw(15) << "Flux Unit:" << tmp << endl;
  Vector<Float> vec;
  oss << setw(15) << "Rest Freqs:";
  if (vec.nelements() > 0) {
      oss << setprecision(10) << vec << " [Hz]" << endl;
  } else {
      oss << "none" << endl;
  }
  oss << setw(15) << "Abcissa:" << "channel" << endl;
  oss << selector_.print() << endl;
  oss << endl;
  // main table
  String dirtype = "Position ("
                  + MDirection::showType(dirCol_.getMeasRef().getType())
                  + ")";
  oss << setw(5) << "Scan"
      << setw(15) << "Source"
//      << setw(24) << dirtype
      << setw(10) << "Time"
      << setw(18) << "Integration" << endl
      << setw(5) << "" << setw(10) << "Beam" << dirtype << endl
      << setw(15) << "" << setw(5) << "IF"
      << setw(8) << "Frame" << setw(16)
      << "RefVal" << setw(10) << "RefPix" << setw(12) << "Increment" <<endl;
  oss << asap::SEPERATOR << endl;
  TableIterator iter(table_, "SCANNO");
  while (!iter.pastEnd()) {
    Table subt = iter.table();
    ROTableRow row(subt);
    MEpoch::ROScalarColumn timeCol(subt,"TIME");
    const TableRecord& rec = row.get(0);
    oss << setw(4) << std::right << rec.asuInt("SCANNO")
        << std::left << setw(1) << ""
        << setw(15) << rec.asString("SRCNAME")
        << setw(10) << formatTime(timeCol(0), false);
    // count the cycles in the scan
    TableIterator cyciter(subt, "CYCLENO");
    int nint = 0;
    while (!cyciter.pastEnd()) {
      ++nint;
      ++cyciter;
    }
    oss << setw(3) << std::right << nint  << setw(3) << " x " << std::left
        << setw(6) <<  formatSec(rec.asFloat("INTERVAL")) << endl;

    TableIterator biter(subt, "BEAMNO");
    while (!biter.pastEnd()) {
      Table bsubt = biter.table();
      ROTableRow brow(bsubt);
      MDirection::ROScalarColumn bdirCol(bsubt,"DIRECTION");
      const TableRecord& brec = brow.get(0);
      oss << setw(6) << "" <<  setw(10) << brec.asuInt("BEAMNO");
      oss  << setw(24) << formatDirection(bdirCol(0)) << endl;
      TableIterator iiter(bsubt, "IFNO");
      while (!iiter.pastEnd()) {
        Table isubt = iiter.table();
        ROTableRow irow(isubt);
        const TableRecord& irec = irow.get(0);
        oss << std::right <<setw(8) << "" << std::left << irec.asuInt("IFNO");
        oss << frequencies().print(irec.asuInt("FREQ_ID"));

        ++iiter;
      }
      ++biter;
    }
    ++iter;
  }
  /// @todo implement verbose mode
  return String(oss);
}

std::string Scantable::getTime(int whichrow, bool showdate) const
{
  MEpoch::ROScalarColumn timeCol(table_, "TIME");
  MEpoch me;
  if (whichrow > -1) {
    me = timeCol(uInt(whichrow));
  } else {
    Double tm;
    table_.keywordSet().get("UTC",tm);
    me = MEpoch(MVEpoch(tm));
  }
  return formatTime(me, showdate);
}

std::vector< double > asap::Scantable::getAbcissa( int whichrow ) const
{
  if ( whichrow > table_.nrow() ) throw(AipsError("Illegal ro number"));
  std::vector<double> stlout;
  int nchan = specCol_(whichrow).nelements();
  String us = freqTable_.getUnitString();
  if ( us == "" || us == "pixel" || us == "channel" ) {
    for (int i=0; i<nchan; ++i) {
      stlout.push_back(double(i));
    }
    return stlout;
  }

  const MPosition& mp = getAntennaPosition();
  const MDirection& md = dirCol_(whichrow);
  const MEpoch& me = timeCol_(whichrow);
  Double rf = moleculeTable_.getRestFrequency(mmolidCol_(whichrow));
  SpectralCoordinate spc =
    freqTable_.getSpectralCoordinate(md, mp, me, rf, mfreqidCol_(whichrow));
  Vector<Double> pixel(nchan);
  Vector<Double> world;
  indgen(pixel);
  if ( Unit(us) == Unit("Hz") ) {
    for ( int i=0; i < nchan; ++i) {
      Double world;
      spc.toWorld(world, pixel[i]);
      stlout.push_back(double(world));
    }
  } else if ( Unit(us) == Unit("km/s") ) {
    Vector<Double> world;
    spc.pixelToVelocity(world, pixel);
    world.tovector(stlout);
  }
  return stlout;
}

std::string Scantable::getAbcissaLabel( int whichrow ) const
{
  if ( whichrow > table_.nrow() ) throw(AipsError("Illegal ro number"));
  const MPosition& mp = getAntennaPosition();
  const MDirection& md = dirCol_(whichrow);
  const MEpoch& me = timeCol_(whichrow);
  const Double& rf = mmolidCol_(whichrow);
  SpectralCoordinate spc =
    freqTable_.getSpectralCoordinate(md, mp, me, rf, mfreqidCol_(whichrow));

  String s = "Channel";
  Unit u = Unit(freqTable_.getUnitString());
  if (u == Unit("km/s")) {
    s = CoordinateUtil::axisLabel(spc,0,True,True,True);
  } else if (u == Unit("Hz")) {
    Vector<String> wau(1);wau = u.getName();
    spc.setWorldAxisUnits(wau);

    s = CoordinateUtil::axisLabel(spc,0,True,True,False);
  }
  return s;

}

void asap::Scantable::setRestFrequencies( double rf, const std::string& unit )
{
  ///@todo lookup in line table
  Unit u(unit);
  Quantum<Double> urf(rf, u);
  uInt id = moleculeTable_.addEntry(urf.getValue("Hz"), "", "");
  TableVector<uInt> tabvec(table_, "MOLECULE_ID");
  tabvec = id;
}

void asap::Scantable::setRestFrequencies( const std::string& name )
{
  throw(AipsError("setRestFrequencies( const std::string& name ) NYI"));
  ///@todo implement
}

std::vector< unsigned int > asap::Scantable::rownumbers( ) const
{
  std::vector<unsigned int> stlout;
  Vector<uInt> vec = table_.rowNumbers();
  vec.tovector(stlout);
  return stlout;
}


Matrix<Float> asap::Scantable::getPolMatrix( uInt whichrow ) const
{
  ROTableRow row(table_);
  const TableRecord& rec = row.get(whichrow);
  Table t =
    originalTable_( originalTable_.col("SCANNO") == Int(rec.asuInt("SCANNO"))
                    && originalTable_.col("BEAMNO") == Int(rec.asuInt("BEAMNO"))
                    && originalTable_.col("IFNO") == Int(rec.asuInt("IFNO"))
                    && originalTable_.col("CYCLENO") == Int(rec.asuInt("CYCLENO")) );
  ROArrayColumn<Float> speccol(t, "SPECTRA");
  return speccol.getColumn();
}


}//namespace asap
