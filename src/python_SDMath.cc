//#---------------------------------------------------------------------------
//# python_SDMath.cc: python exposure of c++ SDMath class
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
#include <boost/python.hpp>

#include <casa/aips.h>
#include <casa/Containers/Block.h>
#include <casa/Utilities/CountedPtr.cc>

#include "SDMathWrapper.h"

using namespace casa;
using namespace boost::python;

namespace asap {
  namespace SDMathWrapper {
    SDMemTableWrapper SDMathWrapper::average(boost::python::tuple tp,
                                             const std::vector<bool>& mask,
                                             bool scanAv,
                                             const std::string& weightStr) {
      int n;
      n = extract<int>(tp.attr("__len__")());
      Block<CountedPtr<asap::SDMemTable> > b(n);
      for (int i=0;i< n;++i) {
        SDMemTableWrapper sdmw =
          extract<SDMemTableWrapper>( tp.attr("__getitem__")(i) );
        b[i] = sdmw.getCP();
      }
      Vector<Bool> msk(mask);
      return SDMemTableWrapper(SDMath::average(b,msk,scanAv,weightStr));
    };
  } // namespace SDMathWrapper

  namespace python {
    void python_SDMath() {
      def("quotient", &SDMathWrapper::quotient);
      def("scale", &SDMathWrapper::scale);
      def("scale_insitu", &SDMathWrapper::scaleInSitu);
      def("add", &SDMathWrapper::add);
      def("add_insitu", &SDMathWrapper::addInSitu);
      def("hanning", &SDMathWrapper::hanning);
      def("average", &SDMathWrapper::average);
      def("averagepol", &SDMathWrapper::averagePol);
      def("averagepol_insitu", &SDMathWrapper::averagePolInSitu);
      def("bin", &SDMathWrapper::bin);
      def("stats", &SDMathWrapper::statistic);
    };

  } // python
} // asap

