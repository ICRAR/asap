//
// C++ Implementation: CalibrationManager
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <assert.h>

#include <casa/Arrays/Vector.h>
#include <casa/Exceptions/Error.h>
#include <casa/Utilities/Regex.h>
#include <casa/BasicSL/String.h>
#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>


#include "CalibrationManager.h"
#include "Scantable.h"
#include "STCalTsys.h"
#include "STCalSkyPSAlma.h"

using namespace casa;
using namespace std;

namespace asap {

CalibrationManager::CalibrationManager()
  : target_(0),
    calmode_(""),
    spwlist_(0)
{
  applicator_ = new STApplyCal();
}

CalibrationManager::~CalibrationManager()
{
}

void CalibrationManager::setScantable(ScantableWrapper &s)
{
  os_.origin(LogOrigin("CalibrationManager","setScantable",WHERE));
  os_ << LogIO::DEBUGGING << "set scantable object." << LogIO::POST;
  target_ = s.getCP();
}

void CalibrationManager::setScantableByName(const string &s)
{
  os_.origin(LogOrigin("CalibrationManager","setScantableAsName",WHERE));
  os_ << LogIO::DEBUGGING << "set scantable " << s << "." << LogIO::POST;
  // always plain table
  target_ = new Scantable(s, Table::Plain);
}

void CalibrationManager::addSkyTable(const string &c)
{
  os_.origin(LogOrigin("CalibrationManager","addSkyTable",WHERE));
  os_ << LogIO::DEBUGGING << "add STCalSkyTable " << c << "." << LogIO::POST;
  skytables_.push_back(new STCalSkyTable(c));
}

void CalibrationManager::addTsysTable(const string &c)
{
  os_.origin(LogOrigin("CalibrationManager","addTsysTable",WHERE));
  os_ << LogIO::DEBUGGING << "add STCalTsysTable " << c << "." << LogIO::POST;
  tsystables_.push_back(new STCalTsysTable(c));
}

void CalibrationManager::setMode(const string &mode)
{
  os_.origin(LogOrigin("CalibrationManager","setMode",WHERE));
  os_ << LogIO::DEBUGGING << "set calibration mode to " << mode << "." << LogIO::POST;
  calmode_ = mode;
  calmode_.upcase();
}

void CalibrationManager::setTimeInterpolation(const string &interp, int order)
{
  os_.origin(LogOrigin("CalibrationManager","setTimeInterpolation",WHERE));
  os_ << LogIO::DEBUGGING << "set interpolation method for time axis to " << interp << "." <<  LogIO::POST;
  applicator_->setTimeInterpolation(stringToInterpolationEnum(interp),
                                    order);
}

void CalibrationManager::setFrequencyInterpolation(const string &interp, int order)
{
  os_.origin(LogOrigin("CalibrationManager","setFrequencyInterpolation",WHERE));
  os_ << LogIO::DEBUGGING << "set interpolation method for frequency axis to " << interp << "." << LogIO::POST;
  applicator_->setFrequencyInterpolation(stringToInterpolationEnum(interp),
                                         order);
}

void CalibrationManager::setTsysTransfer(unsigned int from,
                                         const vector<unsigned int> &to)
{
  os_.origin(LogOrigin("CalibrationManager","setTsysTransfer",WHERE));
  os_ << LogIO::DEBUGGING << "associate Tsys IFNO " << from << " with science IFNO [";
  for (size_t i = 0; i < to.size() ; i++) {
    os_ << to[i];
    if (i == to.size() - 1)
      os_ << "].";
    else
      os_ << ", ";
  }
  os_ << LogIO::POST;
  Vector<uInt> v(to);
  applicator_->setTsysTransfer(from, v); 
}

void CalibrationManager::setTsysSpw(const vector<int> &spwlist)
{
  os_.origin(LogOrigin("CalibrationManager","setTsysSpw",WHERE));
  os_ << LogIO::DEBUGGING << "set IFNO for Tsys calibration to [";
  for (size_t i = 0; i < spwlist.size() ; i++) {
    os_ << spwlist[i];
    if (i == spwlist.size() - 1)
      os_ << "].";
    else
      os_ << ", ";
  }
  os_ << LogIO::POST;
  spwlist_ = spwlist;
}

void CalibrationManager::resetCalSetup()
{
  os_.origin(LogOrigin("CalibrationManager","resetCalSetup",WHERE));
  os_ << LogIO::DEBUGGING << "reset all calibration settings except target data ." << LogIO::POST;
  applicator_->reset();
  calmode_ = "";
  spwlist_.clear();
}

void CalibrationManager::reset()
{
  os_.origin(LogOrigin("CalibrationManager","reset",WHERE));
  os_ << LogIO::DEBUGGING << "reset all calibration settings." << LogIO::POST;
  applicator_->completeReset();
  calmode_ = "";
  spwlist_.clear();
}

void CalibrationManager::calibrate()
{
  os_.origin(LogOrigin("CalibrationManager","calibrate",WHERE));
  os_ << LogIO::DEBUGGING << "start calibration with mode " << calmode_ << "." << LogIO::POST;
  assert(!target_.null());
  if (calmode_ == "TSYS") {
    assert(spwlist_.size() > 0);
    STCalTsys cal(target_, spwlist_);
    cal.calibrate();
    tsystables_.push_back(cal.applytable_);
  }
  else if (calmode_ == "PS") {
    // will match DV01-25, DA41-65, PM01-04, CM01-12
    Regex reant("^(DV(0[1-9]|1[0-9]|2[0-5])|DA(4[1-9]|5[0-9]|6[0-5])|PM0[1-4]|CM(0[1-9]|1[1,2]))$");
    const String antname = target_->getAntennaName();
    if (reant.match(antname.c_str(), antname.size()) != String::npos) {
      os_ << LogIO::DEBUGGING << "ALMA specific position-switch calibration." << LogIO::POST; 
      STCalSkyPSAlma cal(target_);
      cal.calibrate();
      skytables_.push_back(cal.applytable_);
    }
    else {
      String msg = "Calibration type " + calmode_ + " for antenna " + antname + " is not supported.";
      os_.origin(LogOrigin("CalibrationManager","calibrate",WHERE));
      os_ << LogIO::SEVERE << msg << LogIO::POST;
      throw AipsError(msg);
    }      
  }
  else {
    String msg = "Calibration type " + calmode_ + " is not supported.";
    os_.origin(LogOrigin("CalibrationManager","calibrate",WHERE));
    os_ << LogIO::SEVERE << msg << LogIO::POST;
    throw AipsError(msg);
  }
}

void CalibrationManager::apply()
{
  os_.origin(LogOrigin("CalibrationManager","apply",WHERE));
  os_ << LogIO::DEBUGGING << "apply calibration to the data." << LogIO::POST;
  applicator_->setTarget(target_);
  for (size_t i = 0; i < tsystables_.size() ; i++)
    applicator_->push(dynamic_cast<STCalTsysTable*>(&(*tsystables_[i])));
  for (size_t i = 0; i < skytables_.size(); i++) 
    applicator_->push(dynamic_cast<STCalSkyTable*>(&(*skytables_[i])));
  applicator_->apply(false);
}

void CalibrationManager::saveCaltable(const string &name)
{
  os_.origin(LogOrigin("CalibrationManager","saveCaltable",WHERE));
  if (calmode_ == "TSYS") {
    os_ << LogIO::DEBUGGING << "save latest STCalTsysTable as " << name << "." << LogIO::POST;
    tsystables_[tsystables_.size()-1]->save(name);
  }
  else {
    os_ << LogIO::DEBUGGING << "save latest STCalSkyTable as " << name << "." << LogIO::POST;
    skytables_[skytables_.size()-1]->save(name);
  }
}

void CalibrationManager::split(const string &name)
{
  os_.origin(LogOrigin("CalibrationManager","split",WHERE));
  os_ << LogIO::DEBUGGING << "split science data and save them to " << name << "." << LogIO::POST;
  applicator_->save(name);
}

STCalEnum::InterpolationType CalibrationManager::stringToInterpolationEnum(const string &s)
{
  String itype(s);
  itype.upcase();
  const Char *c = itype.c_str();
  String::size_type len = itype.size();
  Regex nearest("^NEAREST(NEIGHBOR)?$");
  Regex linear("^LINEAR$");
  Regex spline("^(C(UBIC)?)?SPLINE$");
  Regex poly("^POLY(NOMIAL)?$");
  if (nearest.match(c, len) != String::npos) {
    return STCalEnum::NearestInterpolation;
  }
  else if (linear.match(c, len) != String::npos) {
    return STCalEnum::LinearInterpolation;
  }
  else if (spline.match(c, len) != String::npos) {
    return STCalEnum::CubicSplineInterpolation;
  }
  else if (poly.match(c, len) != String::npos) {
    return STCalEnum::PolynomialInterpolation;
  }

  os_.origin(LogOrigin("CalibrationManager","stringToInterpolationEnum",WHERE));
  os_ << LogIO::WARN << "Interpolation type " << s << " is not available. Use default interpolation method." << LogIO::POST;
  return STCalEnum::DefaultInterpolation;
}

}
