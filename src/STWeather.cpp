//
// C++ Implementation: STWeather
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

#include "STWeather.h"


using namespace casa;

namespace asap {

const casa::String STWeather::name_ = "WEATHER";

STWeather::STWeather(casa::Table::TableType tt) :
  STSubTable( name_, tt )
{
  setup();
}


STWeather::~STWeather()
{
}

void asap::STWeather::setup( )
{
  // add to base class table
  table_.addColumn(ScalarColumnDesc<Float>("TEMPERATURE"));
  table_.addColumn(ScalarColumnDesc<Float>("PRESSURE"));
  table_.addColumn(ScalarColumnDesc<Float>("HUMIDITY"));
  table_.addColumn(ScalarColumnDesc<Float>("WINDSPEED"));
  table_.addColumn(ScalarColumnDesc<Float>("WINDAZ"));

  // new cached columns
  temperatureCol_.attach(table_,"TEMPERATURE");
  pressureCol_.attach(table_,"PRESSURE");
  humidityCol_.attach(table_,"HUMIDITY");
  windspeedCol_.attach(table_,"WINDSPEED");
  windazCol_.attach(table_,"WINDAZ");
}

uInt STWeather::addEntry( Float temp, Float pressure, Float humidity,
                          Float wspeed, Float waz )
{
  /// @todo this is a zero implementation as none of the telescopes
  /// fills in this information (yet)
  uInt nrow = table_.nrow();
  if ( nrow == 0 ) {
    table_.addRow();
    TableRow row(table_);
    TableRecord& rec = row.record();
    RecordFieldPtr< Float > rfp;
    rfp.attachToRecord(rec,"TEMPERATURE");
    *rfp = temp;
    rfp.attachToRecord(rec,"PRESSURE");
    *rfp = pressure;
    rfp.attachToRecord(rec,"HUMIDITY");
    *rfp = humidity;
    rfp.attachToRecord(rec,"WINDSPEED");
    *rfp = wspeed;
    rfp.attachToRecord(rec,"WINDAZ");
    *rfp = waz;
    row.put(table_.nrow()-1, rec);
  }
  return 0;
}

}
