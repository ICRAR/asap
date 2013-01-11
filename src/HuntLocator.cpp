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

HuntLocator::HuntLocator(double *v, unsigned int n)
  : Locator(v, n),
    prev_(0)
{}

HuntLocator::~HuntLocator()
{}

unsigned int HuntLocator::locate(double x)
{
  if (n_ == 1)
    return 0;
  bool ascending = (x_[n_-1] >= x_[0]);
  if (ascending) {
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
  if (prev_ > 0 && prev_ < n_) {
    // do hunt
    jl = prev_;
    unsigned int inc = 1;
    if ((x >= x_[jl]) == ascending) {
      // hunt up 
      if (jl >= n_ - 1)
        return jl;
      ju = jl + inc;
      while ((x >= x_[ju]) == ascending) {
        jl = ju;
        inc <<= 1;
        ju = jl + inc;
        if (ju > n_ - 1) {
          ju = n_;
          break;
        }
      }
    }
    else {
      // hunt down
      if (jl == 0) 
        return jl;
      ju = jl;
      jl -= inc;
      while ((x < x_[jl]) == ascending) {
        ju = jl;
        inc <<= 1;
        if (inc >= ju) {
          jl = 0;
          break;
        }
        else
          jl = ju - inc;
      }
    }
  }

  // final bisection phase
  unsigned int jm;
  if (ascending) {
    while (ju - jl > 1) {
      jm = (ju + jl) >> 1;
      if (x > x_[jm])
        jl = jm;
      else
        ju = jm;
    }
  }
  else {
    while (ju - jl > 1) {
      jm = (ju + jl) >> 1;
      if (x < x_[jm])
        jl = jm;
      else
        ju = jm;
    }
  }
  prev_ = jl;
  return ju;
}

}
