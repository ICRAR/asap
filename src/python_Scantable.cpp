//#---------------------------------------------------------------------------
//# python_Scantable.cc: python exposure of c++ Scantable class
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
#include <vector>

#include <boost/python.hpp>
#include <boost/python/args.hpp>

#include "ScantableWrapper.h"

using namespace boost::python;

namespace asap {
  namespace python {

void python_Scantable() {
  class_<ScantableWrapper>("Scantable")
    //.def( init <> () )
    .def( init < const std::string&, const std::string& > () )
    .def( init < const ScantableWrapper& > () )
    .def("_copy", &ScantableWrapper::copy)
    .def("_assign", &ScantableWrapper::assign)
    .def("getif", &ScantableWrapper::getIF)
    .def("getbeam", &ScantableWrapper::getBeam)
    .def("getpol", &ScantableWrapper::getPol)
    .def("getscan", &ScantableWrapper::getScan)
    .def("getcycle", &ScantableWrapper::getCycle)
    .def("nif", &ScantableWrapper::nif)
    .def("nbeam", &ScantableWrapper::nbeam)
    .def("npol", &ScantableWrapper::npol)
    .def("nchan", &ScantableWrapper::nchan)
    .def("nscan", &ScantableWrapper::nscan)
    .def("nrow", &ScantableWrapper::nrow)
    .def("get_fluxunit", &ScantableWrapper::getFluxUnit)
    .def("set_fluxunit", &ScantableWrapper::setFluxUnit)
    .def("_setInstrument", &ScantableWrapper::setInstrument)
    .def("_getspectrum", &ScantableWrapper::getSpectrum,
         (boost::python::arg("whichRow")=0))
    /*
    .def("nstokes", &ScantableWrapper::nStokes)
    .def("_getstokesspectrum", &ScantableWrapper::getStokesSpectrum,
         (boost::python::arg("whichRow")=0),
         (boost::python::arg("linpol")=false) )
    .def("_stokestopolspectrum", &ScantableWrapper::stokesToPolSpectrum,
         (boost::python::arg("whichRow")=0),
         (boost::python::arg("linear")=false),
         (boost::python::arg("thepol")=-1) )
    */
    .def("_getpolarizationlabel", &ScantableWrapper::getPolarizationLabel,
         (boost::python::arg("linear")=false),
         (boost::python::arg("stokes")=false),
         (boost::python::arg("linpol")=false) )
//         (boost::python::arg("thepol")=0) )        // Boost fails with 4 arguments
/*    .def("_setspectrum",&ScantableWrapper::setSpectrum,
         (boost::python::arg("whichRow")=0) )*/
    .def("_getabcissa", &ScantableWrapper::getAbcissa,
         (boost::python::arg("whichRow")=0) )
    .def("_getabcissalabel", &ScantableWrapper::getAbcissaLabel,
         (boost::python::arg("whichRow")=0) )
    .def("_getmask", &ScantableWrapper::getMask,
         (boost::python::arg("whichRow")=0) )
    .def("_gettsys", &ScantableWrapper::getTsys)
    .def("_getsourcename", &ScantableWrapper::getSourceName,
         (boost::python::arg("whichRow")=0) )
    .def("_getelevation", &ScantableWrapper::getElevation,
         (boost::python::arg("whichRow")=0) )
    .def("_getazimuth", &ScantableWrapper::getAzimuth,
         (boost::python::arg("whichRow")=0) )
    .def("_getparangle", &ScantableWrapper::getParAngle,
         (boost::python::arg("whichRow")=0) )
    .def("_gettime", &ScantableWrapper::getTime,
         (boost::python::arg("whichRow")=0) )
    .def("_flag", &ScantableWrapper::flag)
    .def("_save",  &ScantableWrapper::makePersistent)
    .def("_summary",  &ScantableWrapper::summary,
         (boost::python::arg("verbose")=true) )
    .def("_getrestfreqs",  &ScantableWrapper::getRestFrequencies)
    //.def("_setrestfreqs",  &ScantableWrapper::setRestFrequencies)
    .def("_setcoordinfo", &ScantableWrapper::setCoordInfo)
    .def("_getcoordinfo", &ScantableWrapper::getCoordInfo)
    .def("_gethistory", &ScantableWrapper::getHistory)
    .def("_addhistory", &ScantableWrapper::addHistory)
    .def("_getselection", &ScantableWrapper::getSelection)
    .def("_setselection", &ScantableWrapper::setSelection)
    /*
    .def("_addfit", &ScantableWrapper::addFit)
    .def("_getfit", &ScantableWrapper::getSDFitTable)
    */
    .def("_recalcazel", &ScantableWrapper::calculateAZEL)
  ;
};

  } // python
} // asap
