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
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>
#include <casa/Containers/RecordField.h>

#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>
//#include <tables/Tables/TableRow.h>

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
  STGrid( const vector<string> infile ) ;
  virtual ~STGrid() {} ;

  void setFileIn( const string infile ) ;
  void setFileList( const vector<string> infile ) ;

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
  
  string saveData( string outfile="" ) ;

private:
  void init() ;

  // actual gridding
  void gridPerRow() ;
  void gridPerPol() ;

  void setupGrid() ;
  void setupGrid( Int &nx, 
                  Int &ny, 
                  String &cellx, 
                  String &celly, 
                  Double &xmin,
                  Double &xmax,
                  Double &ymin,
                  Double &ymax,
                  String &center ) ;
  void mapExtent( Double &xmin, Double &xmax,
                  Double &ymin, Double &ymax ) ;

  void setData( Array<Complex> &gdata,
                Array<Float> &gwgt ) ;
  
  void getDataPerPol( Array<Complex> &spectra,
                      Array<Double> &direction,
                      Array<Int> &flagtra,
                      Array<Int> &rflag,
                      Array<Float> &weight ) ;
  void getDataPerPol( Array<Float> &spectra,
                      Array<Double> &direction,
                      Array<uChar> &flagtra,
                      Array<uInt> &rflag,
                      Array<Float> &weight ) ;
  Int getDataChunk( Array<Complex> &spectra,
                    Array<Double> &direction,
                    Array<Int> &flagtra,
                    Array<Int> &rflag,
                    Array<Float> &weight ) ;
  Int getDataChunk( Array<Float> &spectra,
                    Array<Double> &direction,
                    Array<uChar> &flagtra,
                    Array<uInt> &rflag,
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

  void prepareTable( Table &tab, String &name ) ;

  Bool pastEnd() ;

  void selectData() ;
  void setupArray() ;

  Bool examine() ;
  void attach( Table &tab ) ;

  void call_ggridsd( Array<Double> &xy,
                     Array<Complex> &values,
                     Int &nvispol,
                     Int &nvischan,
                     Array<Int> &flag,
                     Array<Int> &rflag,
                     Array<Float> &weight,
                     Int &nrow,
                     Int &irow,
                     Array<Complex> &grid,
                     Array<Float> &wgrid,
                     Int &nx,
                     Int &ny,
                     Int &npol,
                     Int &nchan,
                     Int &support,
                     Int &sampling,
                     Vector<Float> &convFunc,
                     Int *chanMap,
                     Int *polMap ) ;

  void initPol( Int ipol ) ;
  void initTable( uInt idx ) ;
  Bool isMultiIF( Table &tab ) ;


  // user input
  Int nxUI_ ;
  Int nyUI_ ;
  String cellxUI_ ;
  String cellyUI_ ;
  String centerUI_ ;

  Block<String> infileList_ ;
  uInt nfile_ ;
  Int ifno_ ;
  Int nx_ ;
  Int ny_ ;
  Int npol_ ;
  Int npolOrg_ ;
  Int nchan_ ;
  Int nrow_ ;
  Block<Int> rows_ ;
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
  Block<Table> tableList_ ;
  Table ptab_ ;
  ROArrayColumn<Float> spectraCol_ ;
  ROArrayColumn<uChar> flagtraCol_ ;
  ROArrayColumn<Double> directionCol_ ;
  ROScalarColumn<uInt> flagRowCol_ ;
  ROArrayColumn<Float> tsysCol_ ;
  ROScalarColumn<Double> intervalCol_ ;
  Int nprocessed_ ;
  Int nchunk_ ;

  Array<Float> spectraF_ ;
  Array<uChar> flagtraUC_ ;
  Array<uInt> rflagUI_ ;

};
}
#endif
