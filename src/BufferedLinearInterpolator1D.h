//
// C++ Interface: BufferedLinearInterpolator1D
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAP_BUFFERED_LINEAR_INTERPOLATOR_1D_H
#define ASAP_BUFFERED_LINEAR_INTERPOLATOR_1D_H

#include "Interpolator1D.h"

namespace asap {

/**
 * Linear interpolation with some buffers for acceleration.
 * @author TakeshiNakazato
 */
class BufferedLinearInterpolator1D : public Interpolator1D {
public:
  // Default constructor.
  BufferedLinearInterpolator1D();

  // Destructor.
  virtual ~BufferedLinearInterpolator1D();

  // Set horizontal (x) and vertical (y) data.
  // @param[in] x pointer to horizontal data.
  // @param[in] y pointer to vertical data.
  // @param[in] n number of data.
  // @see Interpolator1D::setData()
  void setData(double *x, float *y, unsigned int n);

  // Set horizontal data (x).
  // @param[in] x pointer to horizontal data.
  // @param[in] n number of data.
  // @see Interpolator1D::setX()
  void setX(double *x, unsigned int n);

  // Perform interpolation.
  // @param[in] x horizontal location where the value is evaluated 
  //              by interpolation.
  // @return interpolated value at x.
  // @see Interpolator1D::interpolate()
  float interpolate(double x);

private:
  // Numerical factor for linear interpolation.
  double factor_;

  // Previous location.
  double xold_;

  // Previous location as an index 
  unsigned int prev_;

  // Boolean parameter whether buffered values are effective or not.
  bool reusable_;
};

}
#endif
