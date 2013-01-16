//
// C++ Implementation: BisectionLocator
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

#include "BisectionLocator.h"

namespace asap {
BisectionLocator::BisectionLocator()
  : Locator()
{
}


BisectionLocator::BisectionLocator(double *v, unsigned int n, bool copystorage)
  : Locator(v, n, copystorage)
{}

BisectionLocator::~BisectionLocator()
{}

unsigned int BisectionLocator::locate(double x)
{
  if (n_ == 1)
    return 0;

  return bisection(x, 0, n_);
}

}
