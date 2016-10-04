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
#ifndef ASAPSTWEATHER_H
#define ASAPSTWEATHER_H

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
  STWeather() {;}
  explicit STWeather(casacore::Table tab);
  explicit STWeather(const Scantable& parent);

  virtual ~STWeather();

  STWeather& operator=(const STWeather& other);

  casacore::uInt addEntry( casacore::Float temperature, casacore::Float pressure,
                       casacore::Float humidity,
                       casacore::Float windspeed, casacore::Float windaz);

  void getEntry( casacore::Float& temperature, casacore::Float& pressure,
                       casacore::Float& humidity,
                       casacore::Float& windspeed, casacore::Float& windaz,
                       casacore::uInt id) const;

  const casacore::String& name() const { return name_; }

private:
  void setup();
  static const casacore::String name_;
  //casacore::Table table_;
  //casacore::ScalarColumn<casacore::uInt> freqidCol_;
  casacore::ScalarColumn<casacore::Float> pressureCol_, temperatureCol_,
                                  humidityCol_,
                                  windspeedCol_, windazCol_;
};

}

#endif
