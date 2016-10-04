//
// C++ Interface: STHistory
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSTHISTORY_H
#define ASAPSTHISTORY_H

#include <casa/aips.h>
#include <casa/BasicSL/String.h>
#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>

#include "STSubTable.h"

namespace asap {

/**
The History subtable of the Scantable

@author Malte Marquarding
*/
class STHistory : public STSubTable {
public:
  STHistory() {;}
  explicit STHistory(casacore::Table tab);
  explicit STHistory( const Scantable& parent );

  virtual ~STHistory();

  STHistory& operator=(const STHistory& other);

  /**
   * add another row to this table
   * @param item the histor string
   * @return an ID (irrelevant)
   */
  casacore::uInt addEntry( const casacore::String& item);

  /**
   * Just for consistency sake, doesn't do anything
   * as IDs are irrelevant in this class
   * @param  item a string
   * @param id an uInt ID
   */
  void getEntry( casacore::String& item, casacore::uInt id);

  /**
   * add all rows from another STHistory table
   * @param other a STHistory reference
   */
  void append(const STHistory& other);

  /**
   * Get all rows in the table as a vector
   * @return a vector of strings
   */
  std::vector<std::string> getHistory(int nrow=-1, int start=0 ) const;

  const casacore::String& name() const { return name_; }

  int nrow() const { return table_.nrow(); }


  void drop();

private:
  void setup();
  static const casacore::String name_;
  casacore::ScalarColumn<casacore::String> itemCol_;
};

}

#endif
