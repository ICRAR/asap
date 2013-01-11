//
// C++ Implementation: Locator
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

#include "Locator.h"

namespace asap {

Locator::Locator(double *v, unsigned int n)
{
  set(v, n);
}

Locator::~Locator()
{}

void Locator::set(double *v, unsigned int n)
{
  x_ = v;
  n_ = n;
}

}
