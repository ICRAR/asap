//
// C++ Implementation: python_STSelector
//
// Description: This file is the boost::python wrapper for asap::STSelector
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <boost/python.hpp>

#include "STSelector.h"

using namespace boost::python;

namespace asap {
  namespace python {

    void python_STSelector() {
      class_<STSelector>("selector")
        .def( init <> () )
        .def("_getbeams", &STSelector::getBeams)
        .def("_getifs", &STSelector::getIFs)
        .def("_getpols", &STSelector::getPols)
        .def("_getscans", &STSelector::getScans)
        .def("_getcycles", &STSelector::getCycles)
        .def("_reset", &STSelector::reset)
        .def("_setbeams", &STSelector::setBeams)
        .def("_setifs", &STSelector::setIFs)
        .def("_setpols", &STSelector::setPolarizations)
        .def("_setscans", &STSelector::setScans)
        .def("_setcycles", &STSelector::setCycles)
        .def("_setname", &STSelector::setName)
        .def("_settaql", &STSelector::setTaQL)
      ;
    };
  }
}