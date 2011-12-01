//#---------------------------------------------------------------------------
//# python_STGrid.cc: python exposure of c++ STGrid class
//#---------------------------------------------------------------------------
#include <string>

#include <boost/python.hpp>
#include <boost/python/args.hpp>

#include "STGrid.h"
//#include "STGridWrapper.h"

using namespace boost::python;

namespace asap {
  namespace python {

void python_STGrid() {
  //class_<STGridWrapper>("stgrid")
  class_<STGrid>("stgrid")
    .def( init <> () )
    .def( init < const std::string > () )
    //.def( init < const STGrid& > () )
//     .def("_defineimage", &STGridWrapper::defineImage)
//     .def("_setoption", &STGridWrapper::setOption)
//     .def("_grid", &STGridWrapper::grid)
    .def("_defineimage", &STGrid::defineImage)
    .def("_setoption", &STGrid::setOption)
    .def("_grid", &STGrid::grid)
    .def("_setin", &STGrid::setFileIn)
//     .def("_setout", &STGrid::setFileOut)
    .def("_save", &STGrid::saveData)
    ;
};

  } // python
} // asap
