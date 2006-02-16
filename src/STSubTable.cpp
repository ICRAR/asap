//
// C++ Implementation: STSubTable
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

#include "STSubTable.h"

using namespace casa;

namespace asap {

STSubTable::STSubTable(const casa::String& name, casa::Table::TableType tt) :
  type_(tt)

{
  TableDesc td("", "1", TableDesc::Scratch);
  td.addColumn(ScalarColumnDesc<uInt>("ID"));
  SetupNewTable aNewTab(name, td, Table::New);
  table_ = Table(aNewTab, tableType());
  idCol_.attach(table_,"ID");

}

STSubTable::~STSubTable()
{
}

}
