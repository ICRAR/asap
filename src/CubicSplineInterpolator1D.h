//
// C++ Interface: CubicSplineInterpolator1D
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAP_CUBIC_SPLINE_INTERPOLATOR_1D_H
#define ASAP_CUBIC_SPLINE_INTERPOLATOR_1D_H

#include "Interpolator1D.h"

namespace asap {

/**
 * Implementation of (natural) cubic spline interpolation.
 * @author TakeshiNakazato
 */
class CubicSplineInterpolator1D : public Interpolator1D {
public:
  // Default constructor.
  CubicSplineInterpolator1D();

  // Destructor.
  virtual ~CubicSplineInterpolator1D();

  // Override Interpolator1D::setData.
  // @see Interpolator1D::setData
  void setData(double *x, float *y, unsigned int n);

  // Override Interpolator1D::setY.
  // @see Interpolator1D::setY()
  void setY(float *y, unsigned int n);

  // Perform interpolation.
  // @param[in] x horizontal location where the value is evaluated 
  //              by interpolation.
  // @return interpolated value at x.
  float interpolate(double x);
private:
  // Determine second derivatives of each point based on 
  // natural cubic spline condition (second derivative at each 
  // end is zero).
  void evaly2();

  // Do interpolation using second derivatives determined by evaly2().
  // @param[in] x horizontal location where the value is evaluated 
  //              by interpolation.
  // @param[in] i location index for x.
  // @return interpolated value at x.
  float dospline(double x, unsigned int i);
  
  // Array to store second derivatives on the data points.
  float *y2_;

  // number of data points for second derivatives
  unsigned int ny2_;

  // Boolean parameter whether buffered values are effective or not.
  bool reusable_;
};

}
#endif
