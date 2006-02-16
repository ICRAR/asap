//
// C++ Interface: STSubTable
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSTSUBTABLE_H
#define ASAPSTSUBTABLE_H

#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>

#include "SDLog.h"

namespace asap {

/**
Abstract base class for all subtables in the Scantable class.

@author Malte Marquarding
@date $Date:$
@version $Revision:$
*/
class STSubTable : public SDLog {
public:
  STSubTable( const casa::String& name,
              casa::Table::TableType tt = casa::Table::Memory );

  virtual ~STSubTable();


  /**
   * Add extra columns. To be implemented in derived class
   */
  virtual void setup() = 0;
  // -> virtual bool conformant(const STSubTable& other) = 0;

  /**
   * Recalculate IDs to be 0-based and incremented by 1 i.e.
   * rowno == ID
   * @return the 'old' IDs
   */
  casa::Vector<casa::uInt> repopulate();

  const casa::Table& table() const { return table_; }
  casa::Table table() { return table_; }

  casa::Table::TableType tableType() const { return type_; }

protected:
  casa::Table table_;
  casa::ScalarColumn<casa::uInt> idCol_;

private:
  casa::Table::TableType type_;

};

}

#endif
