//#---------------------------------------------------------------------------
//# python_SDMemTable.cc: python exposure of c++ SDMemTable class
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

#include "SDMemTableWrapper.h"

using namespace boost::python;

namespace asap {
  namespace python {

void python_SDMemTable() {
  class_<SDMemTableWrapper>("sdtable")
    //.def( init <> () )
    .def( init < std::string > () )
    .def( init < const SDMemTableWrapper& > () )
    .def( init < const SDMemTableWrapper&, std::string > () )
    .def("_copy", &SDMemTableWrapper::copy)
    .def("getif", &SDMemTableWrapper::getIF)
    .def("getbeam", &SDMemTableWrapper::getBeam)
    .def("getpol", &SDMemTableWrapper::getPol)
    .def("_lines", &SDMemTableWrapper::spectralLines)
    .def("nif", &SDMemTableWrapper::nIF)
    .def("nbeam", &SDMemTableWrapper::nBeam)
    .def("npol", &SDMemTableWrapper::nPol)
    .def("nchan", &SDMemTableWrapper::nChan)
    .def("nscan", &SDMemTableWrapper::nScan)
    .def("nrow", &SDMemTableWrapper::nRow)
    .def("setif", &SDMemTableWrapper::setIF,
         (boost::python::arg("whichIF")=0) )
    .def("setbeam", &SDMemTableWrapper::setBeam)
    .def("setpol", &SDMemTableWrapper::setPol)
    .def("rotate_xyphase", &SDMemTableWrapper::rotateXYPhase,
         (boost::python::arg("angle")=0.0) )
    .def("_setmask", &SDMemTableWrapper::setMask)
    .def("get_fluxunit", &SDMemTableWrapper::getFluxUnit)
    .def("set_fluxunit", &SDMemTableWrapper::setFluxUnit)
    .def("_setInstrument", &SDMemTableWrapper::setInstrument)
    .def("_getscan", &SDMemTableWrapper::getScan)
    .def("_getsource", &SDMemTableWrapper::getSource)
    .def("_getspectrum", &SDMemTableWrapper::getSpectrum,
         (boost::python::arg("whichRow")=0))
    .def("_getstokesspectrum", &SDMemTableWrapper::getStokesSpectrum,
         (boost::python::arg("whichRow")=0),
         (boost::python::arg("pol")=false),
         (boost::python::arg("pa")=0.0) )
    .def("_getcircularspectrum", &SDMemTableWrapper::getCircularSpectrum,
         (boost::python::arg("whichRow")=0),
         (boost::python::arg("rr")=true))
    .def("_setspectrum",&SDMemTableWrapper::setSpectrum,
         (boost::python::arg("whichRow")=0) )
    .def("_getabcissa", &SDMemTableWrapper::getAbcissa,
         (boost::python::arg("whichRow")=0) )
    .def("_getabcissalabel", &SDMemTableWrapper::getAbcissaString,
         (boost::python::arg("whichRow")=0) )
    .def("_getmask", &SDMemTableWrapper::getMask,
         (boost::python::arg("whichRow")=0) )
    .def("_gettsys", &SDMemTableWrapper::getTsys)
    .def("_getsourcename", &SDMemTableWrapper::getSourceName,
         (boost::python::arg("whichRow")=0) )
    .def("_gettime", &SDMemTableWrapper::getTime,
         (boost::python::arg("whichRow")=0) )
    .def("_flag", &SDMemTableWrapper::flag,
         (boost::python::arg("whichRow")=-1) )
    .def("_save",  &SDMemTableWrapper::makePersistent)
    .def("_summary",  &SDMemTableWrapper::summary,
	 (boost::python::arg("verbose")=true) )
    .def("_getrestfreqs",  &SDMemTableWrapper::getRestFreqs)
    .def("_setrestfreqs",  &SDMemTableWrapper::setRestFreqs)
    .def("_setcoordinfo", &SDMemTableWrapper::setCoordInfo)
    .def("_getcoordinfo", &SDMemTableWrapper::getCoordInfo)
    .def("_history", &SDMemTableWrapper::history,
         (boost::python::arg("whichRow")=0) )
  ;
};

  } // python
} // asap
