//#---------------------------------------------------------------------------
//# python_SDFitter.cc: python exposure of c++ SDFitter class
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

#include "SDFitter.h"

using namespace boost::python;

namespace asap {
  namespace python {

    void python_SDFitter() {
      class_<SDFitter>("fitter")
        .def( init <> () )
        .def("setexpression", &SDFitter::setExpression)
        .def("setdata", &SDFitter::setData)
        .def("getresidual", &SDFitter::getResidual)
        .def("getfit", &SDFitter::getFit)
        .def("getfixedparameters", &SDFitter::getFixedParameters)
        .def("setfixedparameters", &SDFitter::setFixedParameters)
        .def("getparameters", &SDFitter::getParameters)
        .def("setparameters", &SDFitter::setParameters)
        .def("getestimate", &SDFitter::getEstimate)
        .def("estimate", &SDFitter::computeEstimate)
        .def("geterrors", &SDFitter::getErrors)
        .def("getchi2", &SDFitter::getChisquared)
        .def("reset", &SDFitter::reset)
        .def("fit", &SDFitter::fit)
      ;
    };

  } //namespace python
} // namespace asap
