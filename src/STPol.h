//
// C++ Interface: STPol
//
// Description:
//
//
// Author: Malte Marquarding <Malte.Marquarding@csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSTPOL_H
#define ASAPSTPOL_H

#include <map>
#include <utility>

#include <casa/aips.h>
#include <casa/Exceptions/Error.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Vector.h>
#include "Factory.h"

namespace asap {

/**
Convert betweeen the possible polarisations (linear, circular, stokes, stokes2)

@author Malte Marquarding
@date $Date:$

*/
class STPol {
public:

  typedef  void (STPol::*polOperation)( casa::Float phase );
  STPol() {}
  virtual ~STPol() {}

  typedef FactoryBase<STPol> STPolFactory;

  static STPol* getPolClass( std::map<std::string,STPol::STPolFactory *> factories,
                             const std::string& type )
    { return factories[type]->create(); }

  casa::Vector<casa::Float> getSpectrum( casa::uInt index, const std::string& mode )
    {
      if (mode == "linear")
        return getLinear(index);
      else if ( mode == "stokes" )
        return getStokes(index);
      else if ( mode == "linpol" )
        return getLinPol(index);
      else if ( mode == "cicular" )
        return getLinPol(index);
      else
        throw(casa::AipsError("Polarisation type unknown"));
    }

  virtual casa::Vector<casa::Float> getCircular( casa::uInt index ) = 0;

  virtual casa::Vector<casa::Float> getStokes( casa::uInt index ) = 0;

  virtual casa::Vector<casa::Float> getLinPol( casa::uInt index ) = 0;

  virtual casa::Vector<casa::Float> getLinear( casa::uInt index ) = 0;

  virtual void rotatePhase( casa::Float phase ) {}
  virtual void rotateLinPolPhase( casa::Float phase ) {}

  virtual void invertPhase( casa::Float phase ) {}

  casa::uInt nspec() const { return basespectra_.ncolumn(); }

  const casa::Vector<casa::Float> getSpectrum(casa::uInt index) const
    { return basespectra_.column(index); }

  casa::Matrix<casa::Float>& getSpectra()
    { return basespectra_; }

  void setSpectra(const casa::Matrix<casa::Float>& spec)
    { basespectra_.resize(); basespectra_ = spec; }

  void setPhaseCorrections(casa::Float, casa::Float, casa::Float) {}

  static std::pair<int, std::string> polFromString(const std::string& key);
  static std::string getPolLabel(int index, const std::string& type = "linear");

private:
  static void initPolMap();
  static void initLabelMap();
  static std::map<std::string, std::pair<int, std::string> > polmap_;
  static std::map<std::string, std::map<int, std::string> > labelmap_;

  casa::Vector<casa::Float> phaseCorrections_;
  std::string mode_;
  casa::Matrix<casa::Float> basespectra_;

};

}

#endif
