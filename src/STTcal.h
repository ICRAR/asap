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
  STTcal() {;}
  explicit STTcal(casacore::Table tab);
  explicit STTcal( const Scantable& parent);

  virtual ~STTcal();

  STTcal& operator=(const STTcal& other);

  casacore::uInt addEntry( const casacore::String& time,
                       const casacore::Vector<casacore::Float>& tcal);
  void getEntry( casacore::String& time, casacore::Vector<casacore::Float>& tcal,
                 casacore::uInt id );

  const casacore::String& name() const { return name_; }

private:
  void setup();
  static const casacore::String name_;
  //casacore::Table table_;
  casacore::ArrayColumn<casacore::Float> tcalCol_;
  casacore::ScalarColumn<casacore::String> timeCol_;
};

}

#endif
