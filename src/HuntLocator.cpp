//
// C++ Implementation: HuntLocator
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

#include "HuntLocator.h"

namespace asap {

HuntLocator::HuntLocator()
  : Locator(),
    prev_(0)
{}

HuntLocator::HuntLocator(double *v, unsigned int n, bool copystorage)
  : Locator(v, n, copystorage),
    prev_(0)
{}

HuntLocator::~HuntLocator()
{}

unsigned int HuntLocator::locate(double x)
{
  if (n_ == 1)
    return 0;

  if (ascending_) {
    if (x <= x_[0])
      return 0;
    else if (x > x_[n_-1])
      return n_;
  }
  else {
    if (x > x_[0])
      return 0;
    else if (x <= x_[n_-1])
      return n_;
  }

  unsigned int jl = 0;
  unsigned int ju = n_;

  // hunt phase
  if (prev_ > 0 && prev_ < n_) {
    jl = prev_;
    hunt(x, jl, ju);
  }

  // final bisection phase
  unsigned int j = bisection(x, jl, ju);
  prev_ = (j > 0) ? j - 1 : 0;
  return j;
}

}
