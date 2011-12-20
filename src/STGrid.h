//
// C++ Interface: STGrid
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp>, (C) 2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSTGRID_H
#define ASAPSTGRID_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <casa/BasicSL/String.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Cube.h>
// #include <casa/Arrays/ArrayMath.h>
// #include <casa/Inputs/Input.h>
// #include <casa/Quanta/Quantum.h>
// #include <casa/Quanta/QuantumHolder.h>
// #include <casa/Utilities/CountedPtr.h>

#include <tables/Tables/Table.h>
// #include <tables/Tables/ScalarColumn.h>
// #include <tables/Tables/ArrayColumn.h>

// #include <measures/Measures/MDirection.h>

// #include "Scantable.h"

using namespace std ;
using namespace casa ;

namespace asap {
class STGrid
{
public:
  STGrid() ;
  STGrid( const string infile ) ;
  virtual ~STGrid() {} ;

  void setFileIn( const string infile ) ;

  void setIF( unsigned int ifno ) { ifno_ = ifno ; } ;

  void setPolList( vector<unsigned int> pols ) ;

  void setScanList( vector<unsigned int> scans ) ;

  void defineImage( int nx=-1,
                    int ny=-1,
                    string scellx="",
                    string scelly="",
                    string scenter="" ) ;
  void setFunc( string convtype="box",
                int convsupport=-1 ) ;

  void setWeight( const string wType="uniform" ) ;

  void grid() ;
  void gridPerRow() ;
  
  string saveData( string outfile="" ) ;

private:
  void init() ;

  void setupGrid( Int &nx, 
                  Int &ny, 
                  String &cellx, 
                  String &celly, 
                  Double &xmin,
                  Double &xmax,
                  Double &ymin,
                  Double &ymax,
                  String &center ) ;

  void setData( Array<Complex> &gdata,
                Array<Float> &gwgt ) ;
  
  void getData( Array<Complex> &spectra,
                Array<Double> &direction,
                Array<Int> &flagtra,
                Array<Int> &rflag,
                Array<Float> &weight ) ;
  void getData( Array<Complex> &spectra,
                Array<Double> &direction,
                Array<uChar> &flagtra,
                Array<uInt> &rflag,
                Array<Float> &weight ) ;
  void getDataChunk( Array<Complex> &spectra,
                     Array<Double> &direction,
                     Array<Int> &flagtra,
                     Array<Int> &rflag,
                     Array<Float> &weight ) ;

  void getWeight( Array<Float> &w,
                  Array<Float> &tsys,
                  Array<Double> &tint ) ;

  void toInt( Array<uChar> &u, Array<Int> &v ) ;
  void toInt( Array<uInt> &u, Array<Int> &v ) ;

  void toPixel( Array<Double> &world, Array<Double> &pixel ) ;
  
  void boxFunc( Vector<Float> &convFunc, Int &convSize ) ;
  void spheroidalFunc( Vector<Float> &convFunc ) ;
  void gaussFunc( Vector<Float> &convFunc ) ;
  void pbFunc( Vector<Float> &convFunc ) ;
  void setConvFunc( Vector<Float> &convFunc ) ;
  void selectData( Table &tab ) ;

  Float polMean( const Float *p ) ;
  Double polMean( const Double *p ) ;

  void setupArray( Table &tab ) ;

  void prepareTable( Table &tab, String &name ) ;

  Bool pastEnd() ;


  String infile_ ;
  Int ifno_ ;
  Int nx_ ;
  Int ny_ ;
  Int npol_ ;
  Int nchan_ ;
  Int nrow_ ;
  Double cellx_ ;
  Double celly_ ;
  Vector<Double> center_ ;
  String convType_ ;
  Int convSupport_ ;
  Int userSupport_ ;
  Int convSampling_ ;
  Array<Float> data_ ;
  Vector<uInt> pollist_ ;
  Vector<uInt> scanlist_ ;
  String wtype_ ;

  Table tab_ ;
  Int irow_ ;
};
}
#endif
