//
// C++ Interface: STCalTsysTable
//
// Description:
//
// ApplyTable for Tsys calibration.
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp> (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAP_CALTSYS_TABLE_H
#define ASAP_CALTSYS_TABLE_H

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include "Scantable.h"
#include "STApplyTable.h"

namespace asap {

/**
ApplyTable for Tsys calibration

@author Takeshi Nakazato
@date $Date:$
@version $Revision:$
*/
class STCalTsysTable : public STApplyTable {
public:
  STCalTsysTable() {;}
  STCalTsysTable(const Scantable& parent);
  STCalTsysTable(const casacore::String &name);

  virtual ~STCalTsysTable();

  void setup();
  const casacore::String& name() const {return name_;};

  void attachOptionalColumns();

  void setdata(casacore::uInt irow, casacore::uInt scanno, casacore::uInt cycleno, 
               casacore::uInt beamno, casacore::uInt ifno, casacore::uInt polno, 
               casacore::uInt freqid, casacore::Double time, casacore::Float elevation, 
               const casacore::Vector<casacore::Float> &tsys,
	       const casacore::Vector<casacore::uChar> &flagtra);
  void appenddata(casacore::uInt scanno, casacore::uInt cycleno, 
                  casacore::uInt beamno, casacore::uInt ifno, casacore::uInt polno, 
                  casacore::uInt freqid, casacore::Double time, casacore::Float elevation, 
                  const casacore::Vector<casacore::Float> &tsys,
		  const casacore::Vector<casacore::uChar> &flagtra);
  
  casacore::Vector<casacore::Float> getElevation() const {return elCol_.getColumn();}
  casacore::Matrix<casacore::Float> getTsys() const {return tsysCol_.getColumn();}
  casacore::Matrix<casacore::uChar> getFlagtra() const {return flagtraCol_.getColumn();}
  casacore::uInt nchan(casacore::uInt ifno);

  casacore::Vector<casacore::Double> getBaseFrequency(casacore::uInt whichrow);

private:
  static const casacore::String name_ ;
  casacore::ArrayColumn<casacore::Float> tsysCol_;
  casacore::ArrayColumn<casacore::uChar> flagtraCol_;
  casacore::ScalarColumn<casacore::Float> elCol_;
};

}

#endif
