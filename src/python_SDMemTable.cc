//#---------------------------------------------------------------------------
//# python_SDMemTable.cc: python exposure of c++ SDMemTable class
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
#include <vector>

#include <boost/python.hpp>

#include "SDMemTableWrapper.h"

using namespace boost::python;

namespace atnf_sd {
  namespace python {

void python_SDMemTable() {
  class_<SDMemTableWrapper>("sdtable")
    .def( init < std::string > () )
    .def( init < SDMemTableWrapper, int > () )
    .def("getscan", &SDMemTableWrapper::getScan)
    .def("getspectrum", &SDMemTableWrapper::getSpectrum)
    //.def("getmask", &SDMemTableWrapper::getMask)
    .def("gettsys", &SDMemTableWrapper::getTsys)
    .def("getsourcename", &SDMemTableWrapper::getSourceName)
    .def("gettime", &SDMemTableWrapper::getTime)
    .def("getif", &SDMemTableWrapper::getIF)
    .def("getbeam", &SDMemTableWrapper::getBeam)
    .def("getpol", &SDMemTableWrapper::getPol)
    .def("setif", &SDMemTableWrapper::setIF)
    .def("setbeam", &SDMemTableWrapper::setBeam)
    .def("setpol", &SDMemTableWrapper::setPol)
    .def("makepersistent",  &SDMemTableWrapper::makePersistent)
    .def("summary",  &SDMemTableWrapper::summary)
    .def("name",  &SDMemTableWrapper::name)
  ;
};

  } // python
} // atnf_sd
