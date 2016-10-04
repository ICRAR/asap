//
// C++ Interface: STMolecules
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSTMOLECULES_H
#define ASAPSTMOLECULES_H

#include <casa/aips.h>
#include <casa/BasicSL/String.h>
#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>
#include <casa/Arrays/Array.h>

#include "STSubTable.h"

namespace asap {

/**
The Molecules subtable of the Scantable

@author Malte Marquarding
*/
class STMolecules : public STSubTable {
public:
  STMolecules() {;}
  explicit STMolecules(casacore::Table tab);
  explicit STMolecules( const Scantable& parent);

  virtual ~STMolecules();

  STMolecules& operator=(const STMolecules& other);

/***
  casacore::uInt addEntry( casacore::Double restfreq, const casacore::String& name="",
                       const casacore::String& formattedname="");
***/

  casacore::uInt addEntry( casacore::Vector<casacore::Double> restfreq,
		       const casacore::Vector<casacore::String>& name=casacore::Vector<casacore::String>(0),
                       const casacore::Vector<casacore::String>& formattedname=casacore::Vector<casacore::String>(0));

/***
  void getEntry( casacore::Double& restfreq, casacore::String& name,
                 casacore::String& formattedname, casacore::uInt id) const;
***/
  void getEntry( casacore::Vector<casacore::Double>& restfreq,
		 casacore::Vector<casacore::String>& name,
                 casacore::Vector<casacore::String>& formattedname,
		 casacore::uInt id) const;

  std::vector<double> getRestFrequencies() const;
  std::vector<double> getRestFrequency( casacore::uInt id ) const;
  const casacore::String& name() const { return name_; }
  int nrow() const;

private:
  void setup();
  static const casacore::String name_;
  //casacore::Table table_;
  //casacore::ScalarColumn<casacore::uInt> freqidCol_;
  //casacore::ScalarColumn<casacore::Double> restfreqCol_;
  casacore::ArrayColumn<casacore::Double> restfreqCol_;
  //casacore::ScalarColumn<casacore::String> nameCol_;
  casacore::ArrayColumn<casacore::String> nameCol_;
  //casacore::ScalarColumn<casacore::String> formattednameCol_; // e.g. latex
  casacore::ArrayColumn<casacore::String> formattednameCol_; // e.g. latex

};

}

#endif
