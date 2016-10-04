//
// C++ Interface: STCalSkyTable
//
// Description:
//
// ApplyTable for sky calibration.
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp> (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAP_CALSKY_TABLE_H
#define ASAP_CALSKY_TABLE_H

#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include "Scantable.h"
#include "STApplyTable.h"

namespace asap {

/**
ApplyTable for sky calibration

@author Takeshi Nakazato
@date $Date:$
@version $Revision:$
*/
class STCalSkyTable : public STApplyTable {
public:
  STCalSkyTable() {;}
  STCalSkyTable(const Scantable& parent, const casacore::String &caltype);
  STCalSkyTable(const casacore::String &name);

  virtual ~STCalSkyTable();

  void setup();
  void attachOptionalColumns();

  const casacore::String& name() const {return name_;}

  void setdata(casacore::uInt irow, casacore::uInt scannos, casacore::uInt cycleno, 
               casacore::uInt beamno, casacore::uInt ifno, casacore::uInt polno, 
               casacore::uInt freqid, casacore::Double time, casacore::Float elevation, 
               const casacore::Vector<casacore::Float> &spectra,
	       const casacore::Vector<casacore::uChar> &flagtra);
  void appenddata(casacore::uInt scanno, casacore::uInt cycleno, casacore::uInt beamno, 
                  casacore::uInt ifno, casacore::uInt polno, casacore::uInt freqid,  
                  casacore::Double time, casacore::Float elevation, 
                  const casacore::Vector<casacore::Float> &spectra,
		  const casacore::Vector<casacore::uChar> &flagtra);
  
  casacore::Vector<casacore::Float> getElevation() const {return elCol_.getColumn();}
  casacore::Matrix<casacore::Float> getSpectra() const {return spectraCol_.getColumn();}
  casacore::Matrix<casacore::uChar> getFlagtra() const {return flagtraCol_.getColumn();}
  casacore::uInt nchan(casacore::uInt ifno);

  //casacore::Vector<casacore::Double> getBaseFrequency(casacore::uInt whichrow);

private:
  static const casacore::String name_;
  const casacore::String caltype_;
  casacore::ArrayColumn<casacore::Float> spectraCol_;
  casacore::ArrayColumn<casacore::uChar> flagtraCol_;
  casacore::ScalarColumn<casacore::Float> elCol_;
};

}

#endif
