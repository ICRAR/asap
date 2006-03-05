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

#include "STSubTable.h"

namespace asap {

/**
The Molecules subtable of the Scantable

@author Malte Marquarding
*/
class STMolecules : public STSubTable {
public:
  STMolecules() {;}
  STMolecules(casa::Table tab);
  STMolecules( const Scantable& parent);

  virtual ~STMolecules();

  STMolecules& operator=(const STMolecules& other);

  casa::uInt addEntry( casa::Double restfreq, const casa::String& name="",
                       const casa::String& formattedname="");

  void getEntry( casa::Double& restfreq, casa::String& name,
                 casa::String& formattedname, casa::uInt id) const;

  std::vector<double> getRestFrequencies() const;
  double getRestFrequency( casa::uInt id ) const;
  const casa::String& name() const { return name_; }

private:
  void setup();
  static const casa::String name_;
  //casa::Table table_;
  //casa::ScalarColumn<casa::uInt> freqidCol_;
  casa::ScalarColumn<casa::Double> restfreqCol_;
  casa::ScalarColumn<casa::String> nameCol_;
  casa::ScalarColumn<casa::String> formattednameCol_; // e.g. latex

};

}

#endif
