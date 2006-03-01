//
// C++ Implementation: STFocus
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
#include <tables/Tables/TableRecord.h>
#include <tables/Tables/TableParse.h>
#include <tables/Tables/TableRow.h>
#include <casa/Containers/RecordField.h>

#include "STFocus.h"


using namespace casa;

namespace asap {

const casa::String STFocus::name_ = "FOCUS";

STFocus::STFocus(const Scantable& parent ) :
  STSubTable( parent, name_ )
{
  setup();
}

asap::STFocus::STFocus( casa::Table tab ) : STSubTable(tab, name_)
{
  rotationCol_.attach(table_,"ROTATION");
  angleCol_.attach(table_,"ANGLE");
  tanCol_.attach(table_,"TAN");
}

STFocus::~STFocus()
{
}

STFocus & asap::STFocus::operator =( const STFocus & other )
{
  if (this != &other) {
    static_cast<STSubTable&>(*this) = other;
    rotationCol_.attach(table_,"ROTATION");
    angleCol_.attach(table_,"ANGLE");
    tanCol_.attach(table_,"TAN");
  }
  return *this;
}
void asap::STFocus::setup( )
{
  // add to base class table
  table_.addColumn(ScalarColumnDesc<Float>("ROTATION"));
  table_.addColumn(ScalarColumnDesc<Float>("ANGLE"));
  table_.addColumn(ScalarColumnDesc<Float>("TAN"));

  // new cached columns
  rotationCol_.attach(table_,"ROTATION");
  angleCol_.attach(table_,"ANGLE");
  tanCol_.attach(table_,"TAN");
}

uInt STFocus::addEntry( Float rotation, Float angle, Float ftan)
{
  Table result = table_( near(table_.col("ROTATION"), rotation)
                    && near(table_.col("ANGLE"), angle)
                    && near(table_.col("TAN"), ftan) );
  uInt resultid = 0;
  if ( result.nrow() > 0) {
    ROScalarColumn<uInt> c(result, "ID");
    c.get(0, resultid);
  } else {
    uInt rno = table_.nrow();
    table_.addRow();
    // get last assigned _id and increment
    if ( rno > 0 ) {
      idCol_.get(rno-1, resultid);
      resultid++;
    }
    rotationCol_.put(rno, rotation);
    angleCol_.put(rno, angle);
    tanCol_.put(rno, ftan);
    idCol_.put(rno, resultid);
  }
  return resultid;
}

void asap::STFocus::getEntry( Float& rotation, Float& angle, Float& ftan,
                              uInt id)
{
  Table t = table_(table_.col("ID") == Int(id) );
  if (t.nrow() == 0 ) {
    throw(AipsError("STFocus::getEntry - id out of range"));
  }
  ROTableRow row(t);
  // get first row - there should only be one matching id
  const TableRecord& rec = row.get(0);
  rotation = rec.asFloat("ROTATION");
  angle = rec.asFloat("ANGLE");
  ftan = rec.asFloat("TAN");
}


}
