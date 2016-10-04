//
// C++ Interface: CalibrationManager
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAP_CALIBRATION_MANAGER_H
#define ASAP_CALIBRATION_MANAGER_H

// ASAP
// ScantableWrapper.h must be included first to avoid compiler warnings
// related with _XOPEN_SOURCE
#include "ScantableWrapper.h"
#include "STApplyCal.h"
#include "STCalTsys.h"
#include "STCalibration.h"
#include "STCalEnum.h"

#include <string>
#include <vector>

//#include <boost/scoped_ptr>

#include <casa/BasicSL/String.h>
#include <casa/Utilities/CountedPtr.h>
#include <casa/Logging/LogIO.h>

namespace asap {

/**
 * Class for calibration management. 
 * It also intends to be an interface for Python layer.
 * @author TakeshiNakazato
 */
class CalibrationManager {
public:
  CalibrationManager();

  virtual ~CalibrationManager();

  void setScantable(ScantableWrapper &s);
  void setScantableByName(const std::string &s);
  void addApplyTable(const std::string &c);
  void addSkyTable(const std::string &c);
  void addTsysTable(const std::string &c);
  void setMode(const std::string &mode);
  void setTimeInterpolation(const std::string &interp, int order=-1);
  void setFrequencyInterpolation(const std::string &interp, int order=-1);
  void setTsysSpw(const std::vector<int> &spwlist);
  void setTsysSpwWithRange(const casacore::Record &spwlist, bool average=false);
  void setTsysTransfer(unsigned int from, 
                       const std::vector<unsigned int> &to);
  void setCalibrationOptions(const casacore::Record &options) {options_ = options;}
  void resetCalSetup();
  void reset();
  
  void calibrate();
  void apply(bool insitu=false, bool filltsys=true);
  void saveCaltable(const std::string &name);
  void split(const std::string &name);
private:
  STCalEnum::InterpolationType stringToInterpolationEnum(const std::string &s);

  casacore::Bool isAlmaAntenna();

  casacore::CountedPtr<STApplyCal> applicator_;

  std::vector<casacore::CountedPtr<STApplyTable> > skytables_;
  std::vector<casacore::CountedPtr<STApplyTable> > tsystables_;

  casacore::CountedPtr<Scantable> target_;

  casacore::String calmode_;
  std::vector<int> spwlist_;
  casacore::Record spwlist_withrange_;
  bool do_average_;

  casacore::LogIO os_;

  casacore::Record options_;
};

}
#endif
