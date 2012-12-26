//#---------------------------------------------------------------------------
//# python_STSideBandSep.cpp: python exposure of c++ STSideBandSep class
//#---------------------------------------------------------------------------
//# Author: Kanako Sugimoto, (C) 2012
//#
//# Copyright: See COPYING file that comes with this distribution
//#
//#---------------------------------------------------------------------------
#include <boost/python.hpp>
#include <boost/python/args.hpp>


#include "STSideBandSep.h"

using namespace boost::python;

namespace asap {
  namespace python {

void python_STSideBandSep() {
  class_<STSideBandSep>("SBSeparator")
    .def( init <> () )
    .def( "set_freq", &STSideBandSep::setFrequency,
	  (boost::python::arg("frame")="") )
    .def( "set_lo1", &STSideBandSep::setLO1,
	  (boost::python::arg("frame")="TOPO",
	   boost::python::arg("reftime")=-1,
	   boost::python::arg("refdir")="") )
    //.def( "set_asdm", &STSideBandSep::setLO1Asdm )
    // temporal methods
    .def( "set_imgtable", &STSideBandSep::setImageTable )
    .def( "solve_imgfreq", &STSideBandSep::solveImageFreqency )
  ;
};

  } // python
} // asap
