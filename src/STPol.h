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

@author Malte Marquarding @date $Date:$

*/
class STPol {
public:

  typedef  void (STPol::*polOperation)( casacore::Float phase );
  STPol(): totalangle_(0.0),feedhand_(1.0) {}
  virtual ~STPol() {}

  typedef FactoryBase<STPol> STPolFactory;

  static STPol* getPolClass( std::map<std::string,STPol::STPolFactory *> factories,
                             const std::string& type )
    { return factories[type]->create(); }

  casacore::Vector<casacore::Float> getSpectrum( casacore::uInt index, const std::string& mode )
    {
      if (mode == "linear")
        return getLinear(index);
      else if ( mode == "stokes" )
        return getStokes(index);
      else if ( mode == "linpol" )
        return getLinPol(index);
      else if ( mode == "circular" )
        return getCircular(index);
      else
        throw(casacore::AipsError("Polarisation type unknown"));
    }

  virtual casacore::Vector<casacore::Float> getCircular( casacore::uInt index ) = 0;

  virtual casacore::Vector<casacore::Float> getStokes( casacore::uInt index ) = 0;

  virtual casacore::Vector<casacore::Float> getLinPol( casacore::uInt index ) = 0;

  virtual casacore::Vector<casacore::Float> getLinear( casacore::uInt index ) = 0;

  virtual void rotatePhase( casacore::Float ) {}
  virtual void rotateLinPolPhase( casacore::Float) {}

  virtual void invertPhase( casacore::Float ) {}

  casacore::uInt nspec() const { return basespectra_.ncolumn(); }

  const casacore::Vector<casacore::Float> getSpectrum(casacore::uInt index) const
    { return basespectra_.column(index); }

  casacore::Matrix<casacore::Float>& getSpectra()
    { return basespectra_; }

  void setSpectra(const casacore::Matrix<casacore::Float>& spec)
    { basespectra_.resize(); basespectra_ = spec; }


  void setPhaseCorrections(casacore::Float totalang=0.0, casacore::Float feedhand=1.0)
    { totalangle_=totalang;feedhand_=feedhand;}

  casacore::Float getTotalPhase() const { return totalangle_; }
  casacore::Float getFeedHand() const { return feedhand_; }

  static std::pair<int, std::string> polFromString(const std::string& key);
  static std::string getPolLabel(int index, const std::string& type = "linear");

private:
  static void initPolMap();
  static void initLabelMap();
  static std::map<std::string, std::pair<int, std::string> > polmap_;
  static std::map<std::string, std::map<int, std::string> > labelmap_;

  casacore::Float totalangle_,feedhand_;
  std::string mode_;
  casacore::Matrix<casacore::Float> basespectra_;

};

}

#endif
