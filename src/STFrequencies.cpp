//
// C++ Implementation: STFrequencies
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <casa/iostream.h>
#include <casa/iomanip.h>
#include <casa/Exceptions/Error.h>
#include <casa/Containers/RecordField.h>
#include <casa/Arrays/IPosition.h>

#include <tables/Tables/TableDesc.h>
#include <tables/Tables/SetupNewTab.h>
#include <tables/Tables/ScaColDesc.h>
#include <tables/Tables/TableRecord.h>
#include <tables/Tables/TableParse.h>
#include <tables/Tables/TableRow.h>

#include <coordinates/Coordinates/CoordinateSystem.h>
#include <coordinates/Coordinates/CoordinateUtil.h>

#include "STFrequencies.h"


using namespace casa;

namespace asap {

const String STFrequencies::name_ = "FREQUENCIES";

STFrequencies::STFrequencies(const Scantable& parent) :
  STSubTable(parent, name_)
{
  setup();
}

asap::STFrequencies::STFrequencies( casa::Table tab ) :
  STSubTable(tab, name_)
{
  refpixCol_.attach(table_,"REFPIX");
  refvalCol_.attach(table_,"REFVAL");
  incrCol_.attach(table_,"INCREMENT");

}

STFrequencies::~STFrequencies()
{
}

STFrequencies & asap::STFrequencies::operator =( const STFrequencies & other )
{
  if ( this != &other ) {
    static_cast<STSubTable&>(*this) = other;
    refpixCol_.attach(table_,"REFPIX");
    refvalCol_.attach(table_,"REFVAL");
    incrCol_.attach(table_,"INCREMENT");
  }
  return *this;
}

void STFrequencies::setup( )
{
  // add to base class table
  table_.addColumn(ScalarColumnDesc<Double>("REFPIX"));
  table_.addColumn(ScalarColumnDesc<Double>("REFVAL"));
  table_.addColumn(ScalarColumnDesc<Double>("INCREMENT"));

  table_.rwKeywordSet().define("FRAME", String("TOPO"));
  table_.rwKeywordSet().define("BASEFRAME", String("TOPO"));
  table_.rwKeywordSet().define("EQUINOX",String( "J2000"));
  table_.rwKeywordSet().define("UNIT", String("Hz"));
  table_.rwKeywordSet().define("DOPPLER", String("RADIO"));

  // new cached columns
  refpixCol_.attach(table_,"REFPIX");
  refvalCol_.attach(table_,"REFVAL");
  incrCol_.attach(table_,"INCREMENT");
}

uInt STFrequencies::addEntry( Double refpix, Double refval, Double inc )
{
  // test if this already exists
  Table result = table_( near(table_.col("REFVAL"), refval)
                    && near(table_.col("REFPIX"), refpix)
                    && near(table_.col("INCREMENT"), inc) );
  uInt resultid = 0;
  if ( result.nrow() > 0) {
    ROScalarColumn<uInt> c(result, "ID");
    c.get(0, resultid);

  } else {
    uInt rno = table_.nrow();
    table_.addRow();
    // get last assigned freq_id and increment
    if ( rno > 0 ) {
      idCol_.get(rno-1, resultid);
      resultid++;
    }
    refpixCol_.put(rno, refpix);
    refvalCol_.put(rno, refval);
    incrCol_.put(rno, inc);
    idCol_.put(rno, resultid);
  }
  return resultid;
}



void STFrequencies::getEntry( Double& refpix, Double& refval, Double& inc,
                              uInt id )
{
  Table t = table_(table_.col("ID") == Int(id) );
  if (t.nrow() == 0 ) {
    throw(AipsError("STFrequencies::getEntry - freqID out of range"));
  }
  ROTableRow row(t);
  // get first row - there should only be one matching id
  const TableRecord& rec = row.get(0);
  refpix = rec.asDouble("REFPIX");
  refval = rec.asDouble("REFVAL");
  inc = rec.asDouble("INCREMENT");
}

SpectralCoordinate STFrequencies::getSpectralCoordinate( uInt id ) const
{
  Table t = table_(table_.col("ID") == Int(id) );

  if (t.nrow() == 0 ) {
    throw(AipsError("STFrequencies::getSpectralCoordinate - ID out of range"));
  }

  // get the data
  ROTableRow row(t);
  // get first row - there should only be one matching id
  const TableRecord& rec = row.get(0);

  return SpectralCoordinate( getFrame(true), rec.asDouble("REFVAL"),
                             rec.asDouble("INCREMENT"),
                             rec.asDouble("REFPIX"));
}

SpectralCoordinate
  asap::STFrequencies::getSpectralCoordinate( const MDirection& md,
                                              const MPosition& mp,
                                              const MEpoch& me,
                                              Double restfreq, uInt id ) const
{
  SpectralCoordinate spc = getSpectralCoordinate(id);
  spc.setRestFrequency(restfreq, True);
  if ( !spc.setReferenceConversion(getFrame(), me, mp, md) ) {
    throw(AipsError("Couldn't convert frequency frame."));
  }
  String unitstr = getUnitString();
  if ( !unitstr.empty() ) {
    Unit unitu(unitstr);
    if ( unitu == Unit("Hz") ) {
      Vector<String> wau(1); wau = unitu.getName();
      spc.setWorldAxisUnits(wau);
    } else {
      spc.setVelocity(unitstr, getDoppler());
    }
  }
  return spc;
}


void STFrequencies::rescale( Float factor, const std::string& mode )
{
  TableRow row(table_);
  TableRecord& outrec = row.record();
  RecordFieldPtr<Double> rv(outrec, "REFVAL");
  RecordFieldPtr<Double> rp(outrec, "REFPIX");
  RecordFieldPtr<Double> inc(outrec, "INCREMENT");
  for (uInt i=0; i<table_.nrow(); ++i) {

    const TableRecord& rec = row.get(i);

    SpectralCoordinate sc ( getFrame(true), rec.asDouble("REFVAL"),
                            rec.asDouble("INCREMENT"), rec.asDouble("REFPIX") );

    SpectralCoordinate scout;
    if (mode == "BIN") {
      scout = binCsys(sc, Int(factor));
    } else if (mode == "RESAMPLE") {
      scout = resampleCsys(sc, factor);
    }
    *rv = scout.referenceValue()[0];
    *rp = scout.referencePixel()[0];
    *inc = scout.increment()[0];
    row.put(i);
  }
}

SpectralCoordinate STFrequencies::binCsys(const SpectralCoordinate& sc,
                                          Int factor)
{
  CoordinateSystem csys;
  csys.addCoordinate(sc);
  IPosition factors(1, factor);
  CoordinateSystem binnedcs =
    CoordinateUtil::makeBinnedCoordinateSystem(factors, csys, False);
  return binnedcs.spectralCoordinate(0);
}

SpectralCoordinate STFrequencies::resampleCsys(const SpectralCoordinate& sc,
                                               Float width)
{
  Vector<Float> offset(1,0.0);
  Vector<Float> factors(1,1.0/width);
  Vector<Int> newshape;
  CoordinateSystem csys;
  csys.addCoordinate(sc);
  CoordinateSystem csys2 = csys.subImage(offset, factors, newshape);
  return csys2.spectralCoordinate(0);
}


MFrequency::Types STFrequencies::getFrame(bool base) const
{
  // get the ref frame
  String rf = table_.keywordSet().asString("BASEFRAME");

  // Create SpectralCoordinate (units Hz)
  MFrequency::Types mft;
  if (!MFrequency::getType(mft, rf)) {
    ostringstream oss;
    pushLog("WARNING: Frequency type unknown assuming TOPO");
    mft = MFrequency::TOPO;
  }

  return mft;
}

std::string asap::STFrequencies::getFrameString( bool base ) const
{
  if ( base ) return table_.keywordSet().asString("BASEFRAME");
  else return table_.keywordSet().asString("FRAME");
}

std::string asap::STFrequencies::getUnitString( ) const
{
  return table_.keywordSet().asString("UNIT");
}

Unit asap::STFrequencies::getUnit( ) const
{
  return Unit(table_.keywordSet().asString("UNIT"));
}

std::string asap::STFrequencies::getDopplerString( ) const
{
  return table_.keywordSet().asString("DOPPLER");
}

MDoppler::Types asap::STFrequencies::getDoppler( ) const
{
  String dpl = table_.keywordSet().asString("DOPPLER");

  // Create SpectralCoordinate (units Hz)
  MDoppler::Types mdt;
  if (!MDoppler::getType(mdt, dpl)) {
    throw(AipsError("Doppler type unknown"));
  }
  return mdt;
}

std::string asap::STFrequencies::print( int id )
{
  Table t;
  ostringstream oss;
  if ( id < 0 ) t = table_;
  else  t = table_(table_.col("ID") == Int(id) );
  ROTableRow row(t);
  for (uInt i=0; i<t.nrow(); ++i) {
    const TableRecord& rec = row.get(i);
    oss <<  setw(8)
    << t.keywordSet().asString("FRAME") << setw(16) << setprecision(8)
    << rec.asDouble("REFVAL") << setw(10)
    << rec.asDouble("REFPIX") << setw(12)
    << rec.asDouble("INCREMENT") << endl;
  }
  return String(oss);
}

float STFrequencies::getRefFreq( uInt id, uInt channel )
{
  Table t = table_(table_.col("ID") == Int(id) );
  if ( t.nrow() == 0 ) throw(AipsError("Selected Illegal frequency id"));
  ROTableRow row(t);
  const TableRecord& rec = row.get(0);
  return (Double(channel/2) - rec.asDouble("REFPIX"))
          * rec.asDouble("INCREMENT") + rec.asDouble("REFVAL");
}

bool asap::STFrequencies::conformant( const STFrequencies& other ) const
{
  const Record& r = table_.keywordSet();
  const Record& ro = other.table_.keywordSet();
  return ( r.asString("FRAME") == ro.asString("FRAME") &&
           r.asString("EQUINOX") == ro.asString("EQUINOX") &&
           r.asString("UNIT") == ro.asString("UNIT") &&
           r.asString("DOPPLER") == ro.asString("DOPPLER")
          );
}

std::vector< std::string > asap::STFrequencies::getInfo( ) const
{
  const Record& r = table_.keywordSet();
  std::vector<std::string> out;
  out.push_back(r.asString("UNIT"));
  out.push_back(r.asString("FRAME"));
  out.push_back(r.asString("DOPPLER"));
  return out;
}

void asap::STFrequencies::setInfo( const std::vector< std::string >& theinfo )
{
  if ( theinfo.size() != 3 ) throw(AipsError("setInfo needs three parameters"));
  try {
    setUnit(theinfo[0]);
    setFrame(theinfo[1]);
    setDoppler(theinfo[2]);
  } catch (AipsError& e) {
    throw(e);
  }
}

void asap::STFrequencies::setUnit( const std::string & unit )
{
  if (unit == "" || unit == "pixel" || unit == "channel" ) {
    table_.rwKeywordSet().define("UNIT", "");
  } else {
    Unit u(unit);
    if ( u == Unit("km/s") || u == Unit("Hz") )
      table_.rwKeywordSet().define("UNIT", unit);
    else {
      throw(AipsError("Illegal spectral unit."));
    }
  }
}

void asap::STFrequencies::setFrame( const std::string & frame )
{
  MFrequency::Types mdr;
  if (!MFrequency::getType(mdr, frame)) {
    Int a,b;const uInt* c;
    const String* valid = MFrequency::allMyTypes(a, b, c);
    Vector<String> ftypes(IPosition(1,a), valid);
    ostringstream oss;
    oss <<  String("Please specify a legal frequency type. Types are\n");
    oss << ftypes;
    String msg(oss);
    throw(AipsError(msg));
  } else {
    table_.rwKeywordSet().define("FRAME", frame);
  }
}

void asap::STFrequencies::setDoppler( const std::string & doppler )
{
  MDoppler::Types mdt;
  if (!MDoppler::getType(mdt, doppler)) {
    Int a,b;const uInt* c;
    const String* valid = MDoppler::allMyTypes(a, b, c);
    Vector<String> ftypes(IPosition(1,a), valid);
    ostringstream oss;
    oss <<  String("Please specify a legal doppler type. Types are\n");
    oss << ftypes;
    String msg(oss);
    throw(AipsError(msg));
  } else {
    table_.rwKeywordSet().define("DOPPLER", doppler);
  }
}


} // namespace
