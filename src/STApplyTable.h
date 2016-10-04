//
// C++ Interface: STApplyTable
//
// Description:
//
// Base class for application tables.
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp> (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAP_APPLY_TABLE_H
#define ASAP_APPLY_TABLE_H

#include <casa/Arrays/Vector.h>
#include <casa/Logging/LogIO.h>
#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>

#include "Scantable.h"
#include "STSelector.h"
#include "STCalEnum.h"

namespace asap {

/**
Abstract base class for all application tables.

@author Takeshi Nakazato
@date $Date:$
@version $Revision:$
*/
class STApplyTable  {
public:
  STApplyTable() {;}
  STApplyTable(const Scantable& parent, const casacore::String& name);
  STApplyTable(const casacore::String &name);

  virtual ~STApplyTable();

  /**
   * Add extra columns. To be implemented in derived class
   */
  virtual void setup() = 0;

  /***
   * Name of the table
   ***/
  virtual const casacore::String& name() const = 0;

  const casacore::Table& table() const { return table_; }
  casacore::Table table() { return table_; }
  void attach();
  void attachBaseColumns();
  virtual void attachOptionalColumns() = 0;

  casacore::uInt nrow() {return table_.nrow();}

  casacore::Vector<casacore::uInt> getScan() const {return scanCol_.getColumn();}
  casacore::Vector<casacore::uInt> getCycle() const {return cycleCol_.getColumn();}
  casacore::Vector<casacore::uInt> getBeam() const {return beamCol_.getColumn();}
  casacore::Vector<casacore::uInt> getIF() const {return ifCol_.getColumn();}
  casacore::Vector<casacore::uInt> getPol() const {return polCol_.getColumn();}
  casacore::Vector<casacore::Double> getTime() const {return timeCol_.getColumn();}

  void setSelection(STSelector &sel, bool sortByTime=false);
  void unsetSelection();
  casacore::String caltype();

  void save(const casacore::String &name);

  virtual casacore::uInt nchan(casacore::uInt ifno) = 0;

  // static methods
  static STCalEnum::CalType getCalType(const casacore::String &name);
  static STCalEnum::CalType getCalType(casacore::CountedPtr<STApplyTable> tab);
  static STCalEnum::CalType getCalType(STApplyTable *tab);

protected:
  void setbasedata(casacore::uInt irow, casacore::uInt scanno, casacore::uInt cycleno,
                   casacore::uInt beamno, casacore::uInt ifno, casacore::uInt polno, 
                   casacore::uInt freqid, casacore::Double time);
  casacore::Block<casacore::Double> getFrequenciesRow(casacore::uInt id);

  casacore::Table table_, originaltable_;
  casacore::ScalarColumn<casacore::uInt> scanCol_, cycleCol_, beamCol_, ifCol_, polCol_, freqidCol_;
  casacore::ScalarColumn<casacore::Double> timeCol_;
  casacore::MEpoch::ScalarColumn timeMeasCol_;
  STSelector sel_;
  casacore::LogIO os_;

private:
  static STCalEnum::CalType stringToType(const casacore::String &caltype);
};

}

#endif
