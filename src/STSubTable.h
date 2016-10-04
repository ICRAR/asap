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



namespace asap {

class Scantable;
/**
Abstract base class for all subtables in the Scantable class.

@author Malte Marquarding
@date $Date:$
@version $Revision:$
*/
class STSubTable  {
public:
  STSubTable() {;}
  STSubTable( casacore::Table tab, const casacore::String& name);
  STSubTable( const Scantable& parent, const casacore::String& name );

  virtual ~STSubTable();

  STSubTable& operator=(const STSubTable& other);



  /**
   * Add extra columns. To be implemented in derived class
   */
  virtual void setup() = 0;
  // -> virtual bool conformant(const STSubTable& other) = 0;

  virtual const casacore::String& name() const = 0;

  /**
   * Recalculate IDs to be 0-based and incremented by 1 i.e.
   * rowno == ID
   * @return the 'old' IDs
   */
  casacore::Vector<casacore::uInt> repopulate();

  const casacore::Table& table() const { return table_; }
  casacore::Table table() { return table_; }

protected:
  casacore::Table table_;
  casacore::ScalarColumn<casacore::uInt> idCol_;

private:

};

}

#endif
