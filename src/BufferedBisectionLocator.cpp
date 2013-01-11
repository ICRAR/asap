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

BufferedBisectionLocator::BufferedBisectionLocator(double *v, unsigned int n)
  : Locator(v, n),
    prev_(0)
{}

BufferedBisectionLocator::~BufferedBisectionLocator()
{}

unsigned int BufferedBisectionLocator::locate(double x)
{
  if (n_ == 1)
    return 0;
  bool ascending = (x_[n_-1] >= x_[0]);
  if (ascending) {
    if (x <= x_[0])
      return 0;
    else if (x > x_[n_-1])
      return n_;

    unsigned int jl = 0;
    unsigned int ju = n_;
    unsigned int jm;

    if (x < x_[prev_]) {
      ju = prev_;
      prev_ = 0;
    }
    else 
      jl = prev_;

    while (ju - jl > 1) {
      jm = (ju + jl) >> 1;
      if (x > x_[jm])
        jl = jm;
      else
        ju = jm;
    }
    prev_ = jl;
    return ju;
  }
  else {
    if (x >= x_[0])
      return 0;
    else if (x < x_[n_-1])
      return n_;

    unsigned int jl = 0;
    unsigned int ju = n_;
    unsigned int jm;

    if (x > x_[prev_]) {
      ju = prev_;
      prev_ = 0;
    }
    else 
      jl = prev_;

    while (ju - jl > 1) {
      jm = (ju + jl) >> 1;
      if (x < x_[jm])
        jl = jm;
      else
        ju = jm;
    }
    prev_ = jl;
    return ju;
  }
}

}
