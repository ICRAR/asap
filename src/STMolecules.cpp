//
// C++ Implementation: STMolecules
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

#include "STMolecules.h"


using namespace casa;

namespace asap {

const casa::String STMolecules::name_ = "MOLECULES";

STMolecules::STMolecules(casa::Table::TableType tt) :
  STSubTable( name_, tt )
{
  setup();
}


STMolecules::~STMolecules()
{
}

void asap::STMolecules::setup( )
{
  // add to base class table
  table_.addColumn(ScalarColumnDesc<Double>("RESTFREQUENCY"));
  table_.addColumn(ScalarColumnDesc<String>("NAME"));
  table_.addColumn(ScalarColumnDesc<String>("FORMATTEDNAME"));
  table_.rwKeywordSet().define("UNIT", String("Hz"));
  // new cached columns
  restfreqCol_.attach(table_,"RESTFREQUENCY");
  nameCol_.attach(table_,"NAME");
  formattednameCol_.attach(table_,"FORMATTEDNAME");
}

uInt STMolecules::addEntry( Double restfreq, const String& name,
                            const String& formattedname )
{

  Table result =
    table_( near(table_.col("RESTFREQUENCY"), restfreq) );
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
    restfreqCol_.put(rno, restfreq);
    nameCol_.put(rno, name);
    formattednameCol_.put(rno, formattedname);
    idCol_.put(rno, resultid);
  }
  return resultid;
}

}