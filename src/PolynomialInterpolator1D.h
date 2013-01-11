//
// C++ Interface: PolynomialInterpolator1D
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAP_POLYNOMIAL_INTERPOLATOR_1D_H
#define ASAP_POLYNOMIAL_INTERPOLATOR_1D_H

#include "Interpolator1D.h"

namespace asap {

/**
 * Polynomial interpolation
 * @author TakeshiNakazato
 */
class PolynomialInterpolator1D : public Interpolator1D {
public:
  PolynomialInterpolator1D();

  virtual ~PolynomialInterpolator1D();

  float interpolate(double x);
private:
  float polint(double x, unsigned int loc, unsigned int left, unsigned int n);
};

}
#endif
