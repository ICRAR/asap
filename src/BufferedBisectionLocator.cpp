//
// C++ Implementation: BufferedBisectionLocator
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

#include "BufferedBisectionLocator.h"

namespace asap {
BufferedBisectionLocator::BufferedBisectionLocator()
  : Locator(),
    prev_(0)
{}

BufferedBisectionLocator::BufferedBisectionLocator(double *v, unsigned int n, 
                                                   bool copystorage)
  : Locator(v, n, copystorage),
    prev_(0)
{}

BufferedBisectionLocator::~BufferedBisectionLocator()
{}

unsigned int BufferedBisectionLocator::locate(double x)
{
  if (n_ == 1)
    return 0;

  unsigned int jl = 0;
  unsigned int ju = n_;
  if (ascending_) {
    // ascending order
    if (x <= x_[0])
      return 0;
    else if (x > x_[n_-1])
      return n_;

    if (x < x_[prev_]) {
      ju = bisection(x, jl, prev_);
    }
    else {
      ju = bisection(x, prev_, ju);
    }

  }
  else {
    // descending order
    if (x >= x_[0])
      return 0;
    else if (x < x_[n_-1])
      return n_;

    if (x > x_[prev_]) {
      ju = bisection(x, jl, prev_);
    }
    else {
      ju = bisection(x, prev_, ju);
    }
  }

  prev_ = (ju > 0) ? ju - 1 : 0;
  return ju;    
}

}
