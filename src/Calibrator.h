//
// C++ Interface: Calibrator
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAP_CALIBRATOR_H
#define ASAP_CALIBRATOR_H

#include <memory>

#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <casa/BasicSL/String.h>
#include <casa/Utilities/CountedPtr.h>

namespace asap {

/**
 * Base class for calibration operation 
 * @author TakeshiNakazato
 */
class Calibrator {
public:
  Calibrator();
  Calibrator(unsigned int nchan);

  virtual ~Calibrator();

  void setSource(casacore::Vector<casacore::Float> &v);
  void setReference(casacore::Vector<casacore::Float> &v);
  void setReference2(casacore::Vector<casacore::Float> &v);
  void setScaler(casacore::Vector<casacore::Float> &v);

  const casacore::Vector<casacore::Float> getCalibrated();

  virtual void calibrate() = 0;

protected:
  void initStorage();
  void freeStorage();
  void set(casacore::Float *p, casacore::Vector<casacore::Float> &v);

  unsigned int nchan_;
  unsigned int nchanS_;

  casacore::Float *source_;
  casacore::Float *ref_;
  casacore::Float *ref2_;
  casacore::Float *scaler_;
  casacore::Float *calibrated_;
};

}
#endif
