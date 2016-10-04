//
// C++ Interface: STApplyCal
//
// Description:
//
// Apply any apply tables to target data.
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp> (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAP_APPLY_CAL_H
#define ASAP_APPLY_CAL_H

#include <map>
#include <vector>

#include <casa/Utilities/CountedPtr.h>
#include <casa/Arrays/Vector.h>
#include <casa/Logging/LogIO.h>
#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>

#include "Scantable.h"
#include "STSelector.h"
#include "STApplyTable.h"
#include "STCalEnum.h"
//#include "Calibrator.h"
//#include "Interpolator1D.h"
#include "STCalSkyTable.h"
#include "STCalTsysTable.h"

namespace asap {

template<class T, class U> class Interpolator1D;
class Calibrator;

/**
Apply any apply tables to target data

@author Takeshi Nakazato
@date $Date:$
@version $Revision:$
*/
class STApplyCal  {
public:
  STApplyCal();
  STApplyCal(casacore::CountedPtr<Scantable> target);

  ~STApplyCal();

  // set data
  void setTarget(casacore::CountedPtr<Scantable> target);
  void setTarget(const casacore::String &name);

  // push new caltable 
  void push(STCalSkyTable *table);
  void push(STCalTsysTable *table);

  // set interpolation method
  //void setInterpolation(STCalEnum::InterpolationAxis axis, STCalEnum::InterpolationType itype, casacore::Int order=-1);
  void setTimeInterpolation(STCalEnum::InterpolationType itype, casacore::Int order=-1);
  void setFrequencyInterpolation(STCalEnum::InterpolationType itype, casacore::Int order=-1);

  // set IF (spw) mapping for Tsys transfer
  void setTsysTransfer(casacore::uInt from, casacore::Vector<casacore::uInt> to);

  // apply tables
  void apply(casacore::Bool insitu=false, casacore::Bool filltsys=true);

  // split target data and store it to disk
  void save(const casacore::String &name);

  // reset all settings except target scantable
  void reset();

  // reset all settings
  void completeReset();
  
private:
  // initialization
  void init();

  // setup interpolator
  void initInterpolator();

  // single loop element in apply()
  void doapply(casacore::uInt beamno, casacore::uInt ifno, casacore::uInt polno, 
               casacore::Vector<casacore::uInt> &rows,
               casacore::Vector<casacore::uInt> &skylist,
               casacore::Bool filltsys=true);

  // get frequency information from FREQUENCIES subtable
  casacore::Vector<casacore::Double> getBaseFrequency(casacore::uInt whichrow);

  // search spwmap_ to get IFNO for Tsys
  casacore::uInt getIFForTsys(casacore::uInt to);

  // target data
  casacore::CountedPtr<Scantable> target_;

  // working data
  casacore::CountedPtr<Scantable> work_;

  // calibrator
  casacore::CountedPtr<Calibrator> calibrator_;

  // interpolation method
  STCalEnum::InterpolationType iTime_;
  STCalEnum::InterpolationType iFreq_;
  casacore::Int order_;
  casacore::CountedPtr<Interpolator1D<casacore::Double, casacore::Float> > interpolatorT_;
  casacore::CountedPtr<Interpolator1D<casacore::Double, casacore::Float> > interpolatorF_;
  casacore::CountedPtr<Interpolator1D<casacore::Double, casacore::Float> > interpolatorS_;

  // IF (spw) mapping for Tsys transfer
  map<casacore::uInt, casacore::Vector<casacore::uInt> > spwmap_;

  // list of apply tables
  std::vector<STCalSkyTable*> skytable_;
  std::vector<STCalTsysTable*> tsystable_;

  // calibration type
  STCalEnum::CalType caltype_;
  casacore::Bool doTsys_;

  // selector
  STSelector sel_;

  // logger
  casacore::LogIO os_;
};

}

#endif
