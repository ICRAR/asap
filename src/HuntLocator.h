//
// C++ Interface: HuntLocator
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAP_HUNT_LOCATOR_H
#define ASAP_HUNT_LOCATOR_H

#include "Locator.h"

namespace asap {

/**
 * Implementation of locate operation by bisection search 
 * @author TakeshiNakazato
 */
class HuntLocator : public Locator {
public:
  HuntLocator() {;}
  HuntLocator(double *v, unsigned int n);

  virtual ~HuntLocator();

  unsigned int locate(double x);
private:
  unsigned int prev_;
};

}
#endif
