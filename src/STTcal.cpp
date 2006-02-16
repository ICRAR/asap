//
// C++ Implementation: STTcal
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <casa/Exceptions/Error.h>
#include <tables/Tables/TableDesc.h>
#include <tables/Tables/SetupNewTab.h>
#include <tables/Tables/ScaColDesc.h>
#include <tables/Tables/ArrColDesc.h>
#include <tables/Tables/TableRecord.h>
#include <tables/Tables/TableParse.h>
#include <tables/Tables/TableRow.h>
#include <casa/Containers/RecordField.h>

#include "STTcal.h"


using namespace casa;

namespace asap {

const casa::String STTcal::name_ = "TCAL";

STTcal::STTcal(casa::Table::TableType tt) :
  STSubTable( name_, tt )
{
  setup();
}


STTcal::~STTcal()
{
}

void asap::STTcal::setup( )
{
  // add to base class table
  table_.addColumn(ScalarColumnDesc<String>("TIME"));
  table_.addColumn(ArrayColumnDesc<Float>("TCAL"));

  // new cached columns
  timeCol_.attach(table_,"TIME");
  tcalCol_.attach(table_,"TCAL");
}

uInt STTcal::addEntry( const String& time, const Vector<Float>& cal)
{
  // test if this already exists
  Table result = table_( table_.col("TIME") == time );
  uInt resultid = 0;
  if ( result.nrow() > 0) {
    ROScalarColumn<uInt> c(result, "ID");
    c.get(0, resultid);
  } else {
    uInt rno = table_.nrow();
    table_.addRow();
    // get last assigned tcal_id and increment
    if ( rno > 0 ) {
      idCol_.get(rno-1, resultid);
      resultid++;
    }
    tcalCol_.put(rno, cal);
    timeCol_.put(rno, time);
    idCol_.put(rno, resultid);
  }
  return resultid;
}

}
