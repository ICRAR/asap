//
// C++ Interface: STWeather
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSTWeaTHER_H
#define ASAPSTWeaTHER_H

#include <casa/aips.h>
#include <casa/BasicSL/String.h>
#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>

#include "STSubTable.h"

namespace asap {

/**
The Weather subtable of the Scantable

@author Malte Marquarding
*/
class STWeather : public STSubTable {
public:
  STWeather( casa::Table::TableType tt = casa::Table::Memory);

  virtual ~STWeather();

  casa::uInt addEntry( casa::Float temperature, casa::Float pressure,
                       casa::Float humidity,
                       casa::Float windspeed, casa::Float windaz);

  void getEntry( casa::Float& temperature, casa::Float& pressure,
                       casa::Float& humidity,
                       casa::Float& windspeed, casa::Float& windaz,
                       casa::uInt id);

private:
  void setup();
  static const casa::String name_;
  //casa::Table table_;
  //casa::ScalarColumn<casa::uInt> freqidCol_;
  casa::ScalarColumn<casa::Float> pressureCol_, temperatureCol_,
                                  humidityCol_,
                                  windspeedCol_, windazCol_;
};

}

#endif
