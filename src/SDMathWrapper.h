//#---------------------------------------------------------------------------
//# SDMathWrapper.h: Wrapper classes to use CountedPtr
//#---------------------------------------------------------------------------
//# Copyright (C) 2004
//# ATNF
//#
//# This program is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by the Free
//# Software Foundation; either version 2 of the License, or (at your option)
//# any later version.
//#
//# This program is distributed in the hope that it will be useful, but
//# WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
//# Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with this program; if not, write to the Free Software Foundation, Inc.,
//# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning this software should be addressed as follows:
//#        Internet email: Malte.Marquarding@csiro.au
//#        Postal address: Malte Marquarding,
//#                        Australia Telescope National Facility,
//#                        P.O. Box 76,
//#                        Epping, NSW, 2121,
//#                        AUSTRALIA
//#
//# $Id:
//#---------------------------------------------------------------------------
#ifndef SDMATHWRAPPER_H
#define SDMATHWRAPPER_H

#include <vector>
#include <string>

#include <boost/python/tuple.hpp>

#include "SDMemTableWrapper.h"
#include "SDMath.h"

namespace asap {

namespace SDMathWrapper {

// Quotient

  SDMemTableWrapper quotient(const SDMemTableWrapper& on,
                             const SDMemTableWrapper& off) {
    return SDMemTableWrapper(SDMath::quotient(on.getCP(),
                                             off.getCP()));
  }

// Multiply

  void scaleInSitu(SDMemTableWrapper& in, casa::Float factor, casa::Bool all)
  {
    SDMemTable* sdmt = in.getPtr();
    SDMath::multiplyInSitu(in.getPtr(),factor, all);
  }
  SDMemTableWrapper scale(const SDMemTableWrapper& in,
                          casa::Float factor, casa::Bool all)
  {
    return SDMemTableWrapper(SDMath::multiply(in.getCP(), factor, all));
  }

// Add

  void addInSitu(SDMemTableWrapper& in, casa::Float offset, casa::Bool all)
  {
    SDMemTable* sdmt = in.getPtr();
    SDMath::addInSitu(in.getPtr(), offset, all);
  }
  SDMemTableWrapper add(const SDMemTableWrapper& in, casa::Float offset, casa::Bool all) 
  {
    return SDMemTableWrapper(SDMath::add(in.getCP(), offset, all));
  }

// Hanning

  SDMemTableWrapper hanning(const SDMemTableWrapper& in) {
    return SDMemTableWrapper(SDMath::hanning(in.getCP()));
  }

// Bin up

  SDMemTableWrapper bin(const SDMemTableWrapper& in,
                        int width) {
    return SDMath::bin(in.getCP(), width);
  }

// Average in time

  SDMemTableWrapper average (boost::python::tuple tpl,
                             const std::vector<bool>& mask,
                             bool scanAv, const std::string& wt);

// Average polarizations

  SDMemTableWrapper averagePol(const SDMemTableWrapper& in,
                               const std::vector<bool>& mask) {
    return SDMath::averagePol(in.getCP(), mask);
  }

// Statistics

  std::vector<float> statistic(const SDMemTableWrapper& in,
            const std::vector<bool>& mask, const std::string& which) {
    return SDMath::statistic(in.getCP(), mask, which);
  }

};

} // namespace
#endif
