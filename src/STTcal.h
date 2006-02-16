//
// C++ Interface: STTcal
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSTTCAL_H
#define ASAPSTTCAL_H

#include <casa/aips.h>
#include <casa/BasicSL/String.h>
#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>

#include "STSubTable.h"

namespace asap {

/**
The Tcal subtable of the Scantable

@author Malte Marquarding
*/
class STTcal : public STSubTable {
public:
  STTcal( casa::Table::TableType tt = casa::Table::Memory);

  virtual ~STTcal();

  casa::uInt addEntry( const casa::String& time,
                       const casa::Vector<casa::Float>& tcal);

private:
  void setup();
  static const casa::String name_;
  //casa::Table table_;
  casa::ArrayColumn<casa::Float> tcalCol_;
  casa::ScalarColumn<casa::String> timeCol_;
};

}

#endif