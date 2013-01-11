//
// C++ Interface: BufferedBisectionLocator
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAP_BUFFERED_BISECTION_LOCATOR_H
#define ASAP_BUFFERED_BISECTION_LOCATOR_H

#include "Locator.h"

namespace asap {

/**
 * Implementation of locate operation by bisection search 
 * @author TakeshiNakazato
 */
class BufferedBisectionLocator : public Locator {
public:
  BufferedBisectionLocator() {;}
  BufferedBisectionLocator(double *v, unsigned int n);

  virtual ~BufferedBisectionLocator();

  unsigned int locate(double x);
private:
  unsigned int prev_;
};

}
#endif
