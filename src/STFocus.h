//
// C++ Interface: STFocus
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSTFOCUS_H
#define ASAPSTFOCUS_H

#include <casa/aips.h>
#include <casa/BasicSL/String.h>
#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>

#include "STSubTable.h"

namespace asap {

/**
The Focus subtable of the Scantable

@author Malte Marquarding
*/
class STFocus : public STSubTable {
public:
  STFocus( casa::Table::TableType tt = casa::Table::Memory);

  virtual ~STFocus();

  casa::uInt addEntry( casa::Float rotation, casa::Float angle,
                       casa::Float ftan);

  void getEntry( casa::Float& rotation, casa::Float& angle,
                       casa::Float& ftan, casa::uInt id);

private:
  void setup();
  static const casa::String name_;
  //casa::Table table_;
  //casa::ScalarColumn<casa::uInt> freqidCol_;
  casa::ScalarColumn<casa::Float> rotationCol_, angleCol_,
                                  tanCol_;
};

}

#endif
