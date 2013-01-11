//
// C++ Interface: Interpolator1D
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAP_INTERPOLATOR_1D_H
#define ASAP_INTERPOLATOR_1D_H

#include "Locator.h"

namespace asap {

/**
 * Base class for interpolation operation 
 * @author TakeshiNakazato
 */
class Interpolator1D {
public:
  Interpolator1D();

  virtual ~Interpolator1D();

  void setData(double *x, float *y, unsigned int n);
  void setX(double *x, unsigned int n);
  void setY(float *y, unsigned int n);
  void reset();

  // currently only effective for polynomial interpolation
  void setOrder(unsigned int order) {order_ = order;}

  virtual float interpolate(double x) = 0;

protected:
  int locate(double x);
  bool isready();
  void createLocator();

  unsigned int order_;
  unsigned int n_;
  double *x_;
  float *y_;
  Locator *locator_;
};

}
#endif
