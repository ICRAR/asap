//
// C++ Interface: STFit
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSTFIT_H
#define ASAPSTFIT_H

#include <casa/aips.h>
#include <casa/BasicSL/String.h>
#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>

#include "STSubTable.h"
namespace asap {

class STFitEntry;
/**
The Fit subtable of the Scantable

@author Malte Marquarding
*/
class STFit : public STSubTable {
public:
  STFit() {;}
  explicit STFit(casacore::Table tab);
  explicit STFit( const Scantable& parent);

  virtual ~STFit();

  STFit& operator=(const STFit& other);

  casacore::uInt addEntry( const STFitEntry& fit, casacore::Int id=-1 );
  void getEntry( STFitEntry& fit, casacore::uInt id ) const;

  const casacore::String& name() const { return name_; }

private:
  void setup();
  static const casacore::String name_;
  casacore::ArrayColumn<casacore::String> funcCol_;
  casacore::ArrayColumn<casacore::Int> compCol_;
  casacore::ArrayColumn<casacore::Double> parCol_;
  //  casacore::ArrayColumn<casacore::Double> errCol_;
  casacore::ArrayColumn<casacore::Bool> maskCol_;
  casacore::ArrayColumn<casacore::String> frameCol_;
};

}

#endif
