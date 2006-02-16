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

const casa::String STFrequencies::name_ = "FREQUENCIES";

STFrequencies::STFrequencies(casa::Table::TableType tt) :
  STSubTable( name_, tt )
{
  setup();
}


STFrequencies::~STFrequencies()
{
}

void STFrequencies::setup( )
{
  // add to base class table
  table_.addColumn(ScalarColumnDesc<Double>("REFPIX"));
  table_.addColumn(ScalarColumnDesc<Double>("REFVAL"));
  table_.addColumn(ScalarColumnDesc<Double>("INCREMENT"));

  table_.rwKeywordSet().define("REFFRAME", String("TOPO"));
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


SpectralCoordinate STFrequencies::getSpectralCoordinate( uInt freqID )
{
  Table t = table_(table_.col("ID") == Int(freqID) );

  if (t.nrow() == 0 ) {
    throw(AipsError("STFrequencies::getSpectralCoordinate - freqID out of range"));
  }

  // get the data
  ROTableRow row(t);
  // get first row - there should only be one matching id
  const TableRecord& rec = row.get(0);

  return SpectralCoordinate( getFrame(), rec.asDouble("REFVAL"),
                             rec.asDouble("INCREMENT"),
                             rec.asDouble("REFPIX"));
}

void STFrequencies::rescale( casa::Float factor, const std::string& mode )
{
  TableRow row(table_);
  TableRecord& outrec = row.record();
  RecordFieldPtr<Double> rv(outrec, "REFVAL");
  RecordFieldPtr<Double> rp(outrec, "REFPIX");
  RecordFieldPtr<Double> inc(outrec, "INCREMENT");
  for (uInt i=0; i<table_.nrow(); ++i) {

    const TableRecord& rec = row.get(i);

    SpectralCoordinate sc ( getFrame(), rec.asDouble("REFVAL"),
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


casa::MFrequency::Types STFrequencies::getFrame( ) const
{
  // get the ref frame
  String rf;
  table_.keywordSet().get("REFFRAME", rf);

  // Create SpectralCoordinate (units Hz)
  MFrequency::Types mft;
  if (!MFrequency::getType(mft, rf)) {
    ostringstream oss;
    pushLog("WARNING: Frequency type unknown assuming TOPO");
    mft = MFrequency::TOPO;
  }

  return mft;
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
    << "frame" << setw(16) << setprecision(8)
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

} // namespace
