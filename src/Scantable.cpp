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
#include "SDAttr.h"

using namespace casa;

namespace asap {

Scantable::Scantable(Table::TableType ttype) :
  type_(ttype),
  freqTable_(ttype),
  focusTable_(ttype),
  weatherTable_(ttype),
  tcalTable_(ttype),
  moleculeTable_(ttype)
{
  setupMainTable();
  table_.rwKeywordSet().defineTable("FREQUENCIES", freqTable_.table());
  table_.rwKeywordSet().defineTable("WEATHER", weatherTable_.table());
  table_.rwKeywordSet().defineTable("FOCUS", focusTable_.table());
  table_.rwKeywordSet().defineTable("TCAL", tcalTable_.table());
  table_.rwKeywordSet().defineTable("MOLECULES", moleculeTable_.table());
  setupHistoryTable();
  setupFitTable();

  originalTable_ = table_;
  attach();
}

Scantable::Scantable(const std::string& name, Table::TableType ttype) :
  type_(ttype),
  freqTable_(ttype)
{
  Table tab(name);
  Int version;
  tab.keywordSet().get("VERSION", version);
  if (version != version_) {
    throw(AipsError("Unsupported version of ASAP file."));
  }
  if ( type_ == Table::Memory )
    table_ = tab.copyToMemoryTable("dummy");
  else
    table_ = tab;

  originalTable_ = table_;
  attach();
}


Scantable::Scantable( const Scantable& other, bool clear )
{
  // with or without data
  if (clear) {
    table_ = TableCopy::makeEmptyMemoryTable(String("dummy"), other.table_,
                                             True);
  } else {
    table_ = other.table_.copyToMemoryTable(String("dummy"));
  }

  originalTable_ = table_;
  attach();
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
  SetupNewTable aNewTab("dummy", td, Table::New);
  table_ = Table(aNewTab, Table::Memory, 0);

  originalTable_ = table_;

}

void asap::Scantable::setupHistoryTable( )
{
  TableDesc tdh("", "1", TableDesc::Scratch);
  tdh.addColumn(ScalarColumnDesc<String>("ITEM"));
  SetupNewTable histtab("history", tdh, Table::New);
  Table histTable(histtab, Table::Memory);
  table_.rwKeywordSet().defineTable("HISTORY", histTable);
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
  SetupNewTable aNewTab("fits", td, Table::New);
  Table aTable(aNewTab, Table::Memory);
  table_.rwKeywordSet().defineTable("FITS", aTable);
}

void Scantable::attach()
{
  historyTable_ = table_.keywordSet().asTable("HISTORY");
  fitTable_ = table_.keywordSet().asTable("FITS");

  timeCol_.attach(table_, "TIME");
  srcnCol_.attach(table_, "SRCNAME");
  specCol_.attach(table_, "SPECTRA");
  flagsCol_.attach(table_, "FLAGTRA");
  tsCol_.attach(table_, "TSYS");
  cycleCol_.attach(table_,"CYCLENO");
  scanCol_.attach(table_, "SCANNO");
  beamCol_.attach(table_, "BEAMNO");
  integrCol_.attach(table_, "INTERVAL");
  azCol_.attach(table_, "AZIMUTH");
  elCol_.attach(table_, "ELEVATION");
  dirCol_.attach(table_, "DIRECTION");
  paraCol_.attach(table_, "PARANGLE");
  fldnCol_.attach(table_, "FIELDNAME");
  rbeamCol_.attach(table_, "REFBEAMNO");

  mfitidCol_.attach(table_,"FIT_ID");
  fitidCol_.attach(fitTable_,"FIT_ID");

  mfreqidCol_.attach(table_, "FREQ_ID");

  mtcalidCol_.attach(table_, "TCAL_ID");

  mfocusidCol_.attach(table_, "FOCUS_ID");

  mmolidCol_.attach(table_, "MOLECULE_ID");
}

void Scantable::putSDHeader(const SDHeader& sdh)
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

SDHeader Scantable::getSDHeader() const
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

int Scantable::rowToScanIndex( int therow )
{
  int therealrow = -1;

  return therealrow;
}

int Scantable::nScan() const {
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
  String tmp;
  table_.keywordSet().get("FluxUnit", tmp);
  return tmp;
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
  Instrument ins = SDAttr::convertInstrument(name, throwIt);
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

int asap::Scantable::nbeam( int scanno ) const
{
  if ( scanno < 0 ) {
    Int n;
    table_.keywordSet().get("nBeam",n);
    return int(n);
  } else {
    // take the first POLNO,IFNO,CYCLENO as nbeam shouldn't vary with these
    Table tab = table_(table_.col("SCANNO") == scanno
                       && table_.col("POLNO") == 0
                       && table_.col("IFNO") == 0
                       && table_.col("CYCLENO") == 0 );
    ROTableVector<uInt> v(tab, "BEAMNO");
    return int(v.nelements());
  }
  return 0;
}

int asap::Scantable::nif( int scanno ) const
{
  if ( scanno < 0 ) {
    Int n;
    table_.keywordSet().get("nIF",n);
    return int(n);
  } else {
    // take the first POLNO,BEAMNO,CYCLENO as nbeam shouldn't vary with these
    Table tab = table_(table_.col("SCANNO") == scanno
                       && table_.col("POLNO") == 0
                       && table_.col("BEAMNO") == 0
                       && table_.col("CYCLENO") == 0 );
    ROTableVector<uInt> v(tab, "IFNO");
    return int(v.nelements());
  }
  return 0;
}

int asap::Scantable::npol( int scanno ) const
{
  if ( scanno < 0 ) {
    Int n;
    table_.keywordSet().get("nPol",n);
    return n;
  } else {
    // take the first POLNO,IFNO,CYCLENO as nbeam shouldn't vary with these
    Table tab = table_(table_.col("SCANNO") == scanno
                       && table_.col("IFNO") == 0
                       && table_.col("BEAMNO") == 0
                       && table_.col("CYCLENO") == 0 );
    ROTableVector<uInt> v(tab, "POLNO");
    return int(v.nelements());
  }
  return 0;
}

int asap::Scantable::nrow( int scanno ) const
{
  if ( scanno < 0 ) {
    return int(table_.nrow());
  } else {
    // take the first POLNO,IFNO,CYCLENO as nbeam shouldn't vary with these
    Table tab = table_(table_.col("SCANNO") == scanno
                       && table_.col("BEAMNO") == 0
                       && table_.col("IFNO") == 0
                       && table_.col("POLNO") == 0 );
    return int(tab.nrow());
  }
  return 0;
}


int asap::Scantable::nchan( int scanno, int ifno ) const
{
  if ( scanno < 0 && ifno < 0 ) {
    Int n;
    table_.keywordSet().get("nChan",n);
    return int(n);
  } else {
    // take the first POLNO,IFNO,CYCLENO as nbeam shouldn't vary with these
    Table tab = table_(table_.col("SCANNO") == scanno
                       && table_.col("IFNO") == ifno
                       && table_.col("BEAMNO") == 0
                       && table_.col("POLNO") == 0
                       && table_.col("CYCLENO") == 0 );
    ROArrayColumn<Float> v(tab, "SPECTRA");
    cout << v.shape(0) << endl;
    return 0;
  }
  return 0;
}

Table Scantable::getHistoryTable() const
{
  return table_.keywordSet().asTable("HISTORY");
}

void Scantable::appendToHistoryTable(const Table& otherHist)
{
  Table t = table_.rwKeywordSet().asTable("HISTORY");

  addHistory(asap::SEPERATOR);
  TableCopy::copyRows(t, otherHist, t.nrow(), 0, otherHist.nrow());
  addHistory(asap::SEPERATOR);
}

void Scantable::addHistory(const std::string& hist)
{
  Table t = table_.rwKeywordSet().asTable("HISTORY");
  uInt nrow = t.nrow();
  t.addRow();
  ScalarColumn<String> itemCol(t, "ITEM");
  itemCol.put(nrow, hist);
}

std::vector<std::string> Scantable::getHistory() const
{
  Vector<String> history;
  const Table& t = table_.keywordSet().asTable("HISTORY");
  uInt nrow = t.nrow();
  ROScalarColumn<String> itemCol(t, "ITEM");
  std::vector<std::string> stlout;
  String hist;
  for (uInt i=0; i<nrow; ++i) {
    itemCol.get(i, hist);
    stlout.push_back(hist);
  }
  return stlout;
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

double Scantable::getInterval(int whichrow) const
{
  if (whichrow < 0) return 0.0;
  Double intval;
  integrCol_.get(Int(whichrow), intval);
  return intval;
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

std::vector<float> Scantable::getSpectrum(int whichrow) const
{
  Vector<Float> arr;
  specCol_.get(whichrow, arr);
  std::vector<float> out;
  arr.tovector(out);
  return out;
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
  selector_.reset();
}

void Scantable::setSelection( const STSelector& selection )
{
  Table tab = const_cast<STSelector&>(selection).apply(originalTable_);
  if ( tab.nrow() == 0 ) {
    throw(AipsError("Selection contains no data. Not applying it."));
  }
  table_ = tab;
  selector_ = selection;
}

std::string asap::Scantable::summary( bool verbose )
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
      << setw(15) << "Polarisations:" << setw(4) << npol() << endl
      << setw(15) << "Channels:"  << setw(4) << nchan() << endl;
  oss << endl;
  String tmp;
  //table_.keywordSet().get("Observer", tmp);
  oss << setw(15) << "Observer:" << table_.keywordSet().asString("Observer") << endl;
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

}//namespace asap
