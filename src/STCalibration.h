//
// C++ Interface: STCalibration
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSTCALIBRATION_H
#define ASAPSTCALIBRATION_H

#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <casa/BasicSL/String.h>
#include <casa/Utilities/CountedPtr.h>
#include <casa/Logging/LogIO.h>
#include <casa/Containers/Record.h>

#include <scimath/Mathematics/InterpolateArray1D.h>

#include "Scantable.h"
#include "STDefs.h"
#include "STApplyTable.h"
#include "STSelector.h"

namespace asap {

/**
 * Calibration operations on Scantable objects
 * @author TakeshiNakazato
 */
class STCalibration {
public:
  STCalibration(casacore::CountedPtr<Scantable> &s, const casacore::String target_column);

  void calibrate();

  virtual ~STCalibration() {;}

  void save(casacore::String name) {applytable_->save(name);}
  //const STApplyTable &applytable() {return *applytable_;}
  const casacore::CountedPtr<STApplyTable> applytable() {return applytable_;}

  void setOption(casacore::Record &rec) {options_ = rec;}

protected:
  virtual void setupSelector(const STSelector &sel) = 0;
  virtual void fillCalTable();
  virtual void appenddata(casacore::uInt scanno, casacore::uInt cycleno, 
			  casacore::uInt beamno, casacore::uInt ifno, casacore::uInt polno, 
			  casacore::uInt freqid, casacore::Double time, casacore::Float elevation, 
			  const casacore::Vector<casacore::Float> &any_data,
			  const casacore::Vector<casacore::uChar> &channel_flag) = 0;

  STSelector sel_;
  casacore::CountedPtr<Scantable> scantable_;
  casacore::CountedPtr<STApplyTable> applytable_; 
  casacore::LogIO os_;
  casacore::Record options_;
  const casacore::String target_column_;
};
 
}
#endif
