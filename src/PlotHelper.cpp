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

// casacore
#include <casa/Arrays/Vector.h>
#include <casa/Quanta/Quantum.h>
#include <casa/Quanta/QuantumHolder.h>
#include <casa/Logging/LogIO.h>

#include <measures/Measures/MDirection.h>
#include <tables/Tables/TableRecord.h>


#include "PlotHelper.h"

//#ifndef KS_DEBUG
//#define KS_DEBUG
//#endif

using namespace std ;
using namespace casa ;
using namespace asap ;

namespace asap {

PlotHelper::PlotHelper() : dircoord_(0) 
{
#ifdef KS_DEBUG
  cout << "Default constructor nrow = " << data_p->nrow() << endl;
#endif
};

PlotHelper::PlotHelper( const ScantableWrapper &s) : dircoord_(0)
{
#ifdef KS_DEBUG
  cout << "Constructing PlotHelper with scantable wrapper" << endl;
#endif
  setScantable(s);
};

PlotHelper::~PlotHelper(){
#ifdef KS_DEBUG
  cout << "Called PlotHelper destructor" << endl;
#endif
  if (dircoord_){
#ifdef KS_DEBUG
    cout << "Destructing dircoord_" << endl;
#endif
    delete dircoord_;
    dircoord_ = 0;
  }
};

void PlotHelper::setScantable( const ScantableWrapper &s )
{
#ifdef KS_DEBUG
  cout << "Setting scantable" << endl;
#endif
  data_p = s.getCP();
};


DirectionCoordinate PlotHelper::getSTCoord(const int nx, const int ny,
					   const Projection::Type ptype)
{
  LogIO os(LogOrigin("PlotHelper","getSTCoord()", WHERE));
  os << "Start getSTCoord()" << LogIO::POST;
  if (data_p->nrow() < 1)
    throw AipsError("Scantable is not set. Please set a scantable first.");

  // First, generate rough direction coordinate.
  DirectionCoordinate coord;
  Double incx, incy;
  MDirection::Types mdt;
  ROArrayColumn<Double> dircol;
  Double xmax, xmin, ymax, ymin;
  Double centx, centy;
  Matrix<Double> xform(2,2);
  xform = 0.0;
  xform.diagonal() = 1.0;
  // Rough estimates of center and cell from scantable DIRECTIONs.
  dircol.attach( data_p->table(), "DIRECTION" );
  const Vector<String> udir = dircol.keywordSet().asArrayString("QuantumUnits");
  const Matrix<Double> direction = dircol.getColumn();
  minMax(xmin, xmax, direction.row(0));
  minMax(ymin, ymax, direction.row(1));
  if (!MDirection::getType(mdt, data_p->getDirectionRefString()))
    throw AipsError("Failed to get direction reference from scantable.");
  centx = 0.5 * (xmin + xmax);
  centy = 0.5 * (ymin + ymax);
  incx = abs(xmax - xmin) / (double) nx * cos(centy);
  incy = abs(ymax - ymin) / (double) ny;
  // Generate a temporal direction coordinte
  coord = DirectionCoordinate(mdt, ptype,
				centx, centy, incx, incy, xform,
				0.5*Double(nx), 0.5*Double(ny));
  coord.setWorldAxisUnits(udir);
#ifdef KS_DEBUG
  {//Debug outputs
    cout << "Generating a coordinate from DIRECTIONs in a scantable: " << endl;
    Vector<String> units = coord.worldAxisUnits();
    Vector<Double> refv = coord.referenceValue();
    cout <<"- Reference: " << MDirection::showType(coord.directionType()) << " " << refv[0] << units[0] << " " << refv[1] << units[1]  << endl;
    Vector<Double> refpix = coord.referencePixel();
    cout <<"- Reference Pixel: [" << refpix[0] << ", " << refpix[1] << "]" << endl;
    Vector<Double> inc = coord.increment();
    cout <<"- Increments: [" << inc[0] << ", " << inc[1] << "]" << endl;
    cout <<"- Projection: " << coord.projection().name() << endl;
  }
#endif
  return coord;
}

void PlotHelper::setGridParam(const int nx, const int ny,
			      const string cellx, const string celly,
			      string center, const string projname)
{
  LogIO os(LogOrigin("PlotHelper","setGridParam()", WHERE));
  os << "Start setGridParam()" << LogIO::POST;
  // Value check of nx and ny
  if (nx < 1)
    throw(AipsError("nx should be > 0"));
  if (ny < 1)
    throw(AipsError("ny should be > 0"));
  // Destroy old coord
  if (dircoord_){
#ifdef KS_DEBUG
    cout << "Destructing dircoord_" << endl;
#endif
    delete dircoord_;
    dircoord_ = 0;
  }

  // Check for availability of data_p
  const bool stset = ((data_p->nrow() > 0) ? true : false);
  const bool needst = center.empty() || (cellx.empty() && celly.empty());
  if (needst && !stset)
    throw AipsError("Could not resolve grid parameters. Please set a scantable first.");

  // projection
  Projection::Type projtype(Projection::type(String(projname)));

  // Calculate projected map center (in world coordinate)
  // and extent (in pixel coordinate) from scantable DIRECIONs (if necessary).
  MDirection stcent;
  Double stxmax, stxmin, stymax, stymin;
  MDirection::Types stdt;
  DirectionCoordinate stcoord;
  Quantum<Double> stincx, stincy;
  if (needst && stset) {
    stcoord = getSTCoord(nx, ny, projtype);
    Vector<Double> inc = stcoord.increment();
    Vector<String> units = stcoord.worldAxisUnits();
    stincx = Quantum<Double>(inc[0], units[0]);
    stincy = Quantum<Double>(inc[1], units[0]);
    stdt = stcoord.directionType();
    // Get the extent of directions in Pixel coordinate
    ROArrayColumn<Double> dircol;
    dircol.attach( data_p->table(), "DIRECTION" );
    const Matrix<Double> direction = dircol.getColumn();
    Matrix<Double> dirpix;
    Vector<Bool> failures;
    if (!stcoord.toPixelMany(dirpix, direction, failures))
      throw AipsError("Failed to get directions in pixel coordinate.");
    minMax(stxmin, stxmax, dirpix.row(0));
    minMax(stymin, stymax, dirpix.row(1));
    // Get the direction center in World coordinate.
    Vector<Double> centpix(2);
    centpix[0] = 0.5 * (stxmin + stxmax);
    centpix[1] = 0.5 * (stymin + stymax);
    stcoord.toWorld(stcent, centpix);
#ifdef KS_DEBUG
    {//Debug output
      Quantum< Vector<Double> > qcent;
      cout << "Got the center and map extent of scantable." << endl;
      qcent = stcent.getAngle();
      cout << "- Center of DIRECTIONs: " << stcent.getRefString() << " "
	   << qcent.getValue()[0] << qcent.getUnit() << " "
	   << qcent.getValue()[1] << qcent.getUnit() << endl;
      cout << "- Map extent: [" << (stxmax-stxmin)*stincx.getValue() << ", "
	   << (stymax-stymin)*stincy.getValue() << "] ( " << stincx.getUnit() << ")" << endl;
    }
#endif
  }

  // Now, define direction coordinate from input parameters (in radian)
  Double incx, incy;
  Double centx, centy;
  MDirection::Types mdt;
  // center
  if (center.empty()){
    if (!stset)
      throw AipsError("Scantable is not set. Could not resolve map center.");
#ifdef KS_DEBUG
    cout << "Using pointing center from DIRECTION column" << endl;
#endif
    centx = stcent.getAngle("rad").getValue()[0];
    centy = stcent.getAngle("rad").getValue()[1];
    mdt = stdt;
  } else {
#ifdef KS_DEBUG
    cout << "Using user defined grid center" << endl;
#endif
    // Parse center string
    string::size_type pos0 = center.find(" ");
    
    if (pos0 == string::npos)
      throw AipsError("bad string format in center direction");

    string::size_type pos1 = center.find(" ", pos0+1);
    String sepoch, sra, sdec;
    if (pos1 != string::npos) {
      sepoch = center.substr(0, pos0);
      sra = center.substr(pos0+1, pos1-pos0);
      sdec = center.substr(pos1+1);
    } else {
      sepoch = "J2000";
      sra = center.substr(0, pos0);
      sdec = center.substr(pos0+1);
    }
    if (!MDirection::getType(mdt,sepoch))
      throw AipsError("Invalid direction reference in center");
    if (stset && mdt != stdt)
      throw AipsError("Direction reference of center should be the same as input scantable");
    QuantumHolder qh ;
    String err ;
    qh.fromString(err, sra);
    Quantum<Double> ra  = qh.asQuantumDouble();
    qh.fromString(err, sdec) ;
    Quantum<Double> dec = qh.asQuantumDouble();
    centx = ra.getValue("rad");
    centy = dec.getValue("rad");
    // rotaion
    if (stset) {
      Double stcentx = stcent.getAngle("rad").getValue()[0];
      Double rotnum = round( abs(centx - stcentx) / (C::_2pi) );
      if (centx < stcentx) rotnum *= -1;
      centx -= (rotnum * C::_2pi);
    }
  }
#ifdef KS_DEBUG
  cout << "The center direction of plotting grids: [" << centx  << ", " << centy << "] (rad)" <<endl;
#endif

  // cell
  if (cellx.empty() && celly.empty()){
    if (!stset)
      throw AipsError("Scantable is not set. Could not resolve cell size.");
#ifdef KS_DEBUG
    cout << "Using cell size defined from DIRECTION column" << endl;
#endif
    Vector<Double> centpix;
    MDirection centmd = MDirection(Quantum<Double>(centx, "rad"),
				   Quantum<Double>(centy, "rad"), mdt);
    stcoord.toPixel(centpix, centmd);
#ifdef KS_DEBUG
    cout << "- got centpix [" << centpix[0] << ", " << centpix[1] << "]" <<endl;
#endif
    Double wx = max( abs(stxmax-centpix[0]), abs(stxmin-centpix[0]) )
      * 2 * stincx.getValue("rad");
    Double wy = max( abs(stymax-centpix[1]), abs(stymin-centpix[1]) )
      * 2 * stincy.getValue("rad");
    incx = wx / max(nx - 1., 1.);
    incy = wy / max(ny - 1., 1.);
  } else {
#ifdef KS_DEBUG
    cout << "Using user defined cell size" << endl;
#endif
    Quantum<Double> qcellx, qcelly;
    if (!cellx.empty() && !celly.empty()){
      readQuantity(qcellx, String(cellx));
      readQuantity(qcelly, String(celly));
    } else if (celly.empty()) {
      readQuantity(qcellx, String(cellx));
      qcelly = qcellx;
    } else { //only celly is defined
      readQuantity(qcelly, String(celly));
      qcellx = qcelly;
    }
    incx = qcellx.getValue("rad");
    incy = qcelly.getValue("rad");
  }
#ifdef KS_DEBUG
  cout << "The cell size of plotting grids: [" << incx << ", " << incy << "] (rad)" <<endl;
#endif

  Matrix<Double> xform(2,2) ;
  xform = 0.0 ;
  xform.diagonal() = 1.0 ;
  dircoord_ = new DirectionCoordinate(mdt, projtype,
				      centx, centy, incx, incy,
				      xform,
				      0.5*Double(nx), 
				      0.5*Double(ny)) ; // pixel at center
  os << "Successfully finished generation of Direction Coordinate" << LogIO::POST;
};

void PlotHelper::setGridParamVal(const int nx, const int ny,
				 const double cellx, const double celly,
				 const double centx, const double centy,
				 const string epoch, const string projname){
  // Value check of nx and ny
  if (nx < 1)
    throw(AipsError("nx should be > 0"));
  if (ny < 1)
    throw(AipsError("ny should be > 0"));
  // Destroy old coord
  if (dircoord_){
#ifdef KS_DEBUG
    cout << "Destructing dircoord_" << endl;
#endif
    delete dircoord_;
    dircoord_ = 0;
  }

  // center (in rad)
  Double centX(centx), centY(centy);
  // cell size (in rad)
  Double incX(cellx), incY(celly);
  // epoch
  MDirection::Types mdt;
  MDirection::getType(mdt, String(epoch));
  // projection
  Projection::Type projType(Projection::type(String(projname)));

  Matrix<Double> xform(2,2) ;
  xform = 0.0 ;
  xform.diagonal() = 1.0 ;
  dircoord_ = new DirectionCoordinate(mdt, projType,
				      centX, centY, incX, incY,
				      xform,
				      0.5*Double(nx), 
				      0.5*Double(ny)) ; // pixel at center
// 				      0.5*Double(nx-1), 
// 				      0.5*Double(ny-1)) ; // pixel at grid
#ifdef KS_DEBUG
  {//Debug outputs
  cout << "Direction coordinate is set: " << endl;
  Vector<String> units = dircoord_->worldAxisUnits();
  Vector<Double> refv = dircoord_->referenceValue();
  cout <<"Reference: " << MDirection::showType(dircoord_->directionType()) << " " << refv[0] << units[0] << " " << refv[1] << units[1]  << endl;
  Vector<Double> refpix = dircoord_->referencePixel();
  cout <<"Reference Pixel: [" << refpix[0] << ", " << refpix[1] << "]" << endl;
  Vector<Double> inc = dircoord_->increment();
  cout <<"Increments: [" << inc[0] << ", " << inc[1] << "]" << endl;
  cout <<"Projection: " << dircoord_->projection().name() << endl;
  }
#endif
};

vector<double>  PlotHelper::getGridPixel(const int whichrow){
  if (data_p->nrow() < 1)
    throw AipsError("Scantable is not set. Could not get direction.");
  else if (whichrow > int(data_p->nrow()) - 1)
    throw AipsError("Row index out of range.");
  if (!dircoord_)
    throw AipsError("Direction coordinate is not defined.");

  Vector<Double> pixel;
  MDirection world;
  vector<double> outvec;
  world = data_p->getDirection(whichrow);
#ifdef KS_DEBUG
  cout << "searching pixel position (world = " << data_p->getDirectionString(whichrow) << " = [" << world.getAngle("rad").getValue()[0] << ", " << world.getAngle("rad").getValue()[1] << "])" << endl;
#endif
  dircoord_->toPixel(pixel, world);
#ifdef KS_DEBUG
  cout << "got pixel = [" << pixel[0] << ", " << pixel[1] << "]" << endl;
#endif
  // convert pixel to std vector
  pixel.tovector(outvec);
  return outvec;
};

} //namespace asap
