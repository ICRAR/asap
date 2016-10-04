//
// C++ Interface: PlotHelper
//
// Description:
//    A small helper class to handle direction coordinate in asapplotter
//
// Author: Kana Sugimoto <kana.sugi@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPPLOTHELPER_H
#define ASAPPLOTHELPER_H

// ASAP
// ScantableWrapper.h must be included first to avoid compiler warnings
// related with _XOPEN_SOURCE
#include "ScantableWrapper.h"
#include "Scantable.h"

// STL
#include <iostream>
#include <string>
#include <vector>
// casacore
#include <casa/aips.h>
#include <casa/Utilities/CountedPtr.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>

namespace asap {

class PlotHelper
{
public:
  /**
   * constructors and a destructor
   **/
  PlotHelper();
  explicit PlotHelper( const ScantableWrapper &s);

  virtual ~PlotHelper();
  /**
   * Set scantable
   **/
  void setScantable( const ScantableWrapper &s ) ;

  /**
   * Set grid parameters for plot panel (by values in radian)
   **/
  void setGridParamVal(const int nx, const int ny, 
		    const double cellx, const double celly,
		    const double centx, const double centy,
		    string epoch="J2000", const string projname="SIN") ;
  /// TODO
  /**
   * Set grid parameters for plot panel (by quantity)
   **/
  void setGridParam(const int nx, const int ny,
		    const string cellx="", const string celly="",
		    string center="", const string projname="SIN") ;

  /**
   * Get Pixel position of a row in the scantable
   **/
  vector<double> getGridPixel(const int whichrow=0);

  /**
   * Get the reference direction of grid coordinate (grid center)
   **/
  string getGridRef();

  /**
   * Get the cell size (>0) of the grid coordinate (in radian)
   **/
  vector<double> getGridCellVal();


private:
  /** Generate temporal coordinate from the DIRECTION column of a scantable**/
  casacore::DirectionCoordinate getSTCoord(const int nx, const int ny,
				       const casacore::Projection::Type ptype);

  /** Generation of direction coordinate **/
  void setupCoord(const casacore::MDirection::Types mdt,
		  const casacore::Projection::Type pjt,
		  const casacore::Double centx, const casacore::Double centy,
		  const casacore::Double incx, const casacore::Double incy,
		  const casacore::Double refx, const casacore::Double refy);

  casacore::DirectionCoordinate *dircoord_p;
  casacore::CountedPtr<Scantable> data_p;

};

} // namespace

#endif
