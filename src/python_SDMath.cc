//#---------------------------------------------------------------------------
//# python_SDMath.cc: python exposure of c++ SDMath class
//#---------------------------------------------------------------------------
//# Copyright (C) 2004
//# Malte Marquarding, ATNF
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
#include <boost/python.hpp>

#include "SDMathWrapper.h"

using namespace boost::python;

namespace atnf_sd {
  namespace python {

SDMemTableWrapper SDMathWrapper::averages(boost::python::tuple tp,
					  const std::vector<bool>& mask) {
  int n;
  n = extract<int>(tp.attr("__len__")());
  Block<CountedPtr<atnf_sd::SDMemTable> > b(n);
  for (int i=0;i< n;++i) {
    SDMemTableWrapper sdmw = 
      extract<SDMemTableWrapper>( tp.attr("__getitem__")(i) );
    b[i] = sdmw.getCP();
  }
  Vector<Bool> msk(mask);
  return SDMemTableWrapper(SDMath::averages(b,msk));
};

void python_SDMath() {
  class_<SDMathWrapper>("sdmath")
    .def("average", &SDMathWrapper::average)
    .def("quotient", &SDMathWrapper::quotient)
    .def("multiply", &SDMathWrapper::multiply)
    .def("baseline", &SDMathWrapper::baseline)
    .def("hanning", &SDMathWrapper::hanning)
    .def("averages", &SDMathWrapper::averages)
    .def("averagepol", &SDMathWrapper::averagePol)
    .def("bin", &SDMathWrapper::bin)
    .def("rms", &SDMathWrapper::rms)
    ;
};

  } // python
} // atnf_sd

