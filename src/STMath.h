//
// C++ Interface: STMath
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSTMATH_H
#define ASAPSTMATH_H

#include <string>
#include <map>

#include <casa/aips.h>
#include <casa/Utilities/CountedPtr.h>
#include <casa/BasicSL/String.h>
#include <casa/Arrays/Vector.h>
#include <scimath/Mathematics/InterpolateArray1D.h>

#include "Scantable.h"
#include "SDDefs.h"
#include "SDLog.h"

namespace asap {

/**
Mathmatical operations on Scantable objects

@author Malte Marquarding
*/
class STMath : private SDLog {
public:

  typedef casa::InterpolateArray1D<casa::Double,
                                   casa::Float>::InterpolationMethod imethod;

  typedef std::map<std::string, imethod> imap;

  STMath(bool insitu=true);

  ~STMath();

  /**
   * set the @attr insitu attribute
   * @param b
   */
  bool insitu() const { return insitu_;};
  void setInsitu(bool b) { insitu_ = b; };

  casa::CountedPtr<Scantable>
    average( const std::vector<casa::CountedPtr<Scantable> >& in,
             const casa::Vector<casa::Bool>& mask,
             const std::string& weight,
             const std::string& avmode = "",
             bool alignfreq = false );

  casa::CountedPtr<Scantable>
    unaryOperate( const casa::CountedPtr<Scantable>& in, float val,
                  const std::string& mode, bool tsys=false );

  casa::CountedPtr<Scantable> quotient( const casa::CountedPtr<Scantable>& in,
                                        const std::string& mode = "NEAREST",
                                        bool preserve = true );

  casa::CountedPtr<Scantable>
    freqSwitch( const casa::CountedPtr<Scantable>& in );

  std::vector<float> statistic(const casa::CountedPtr<Scantable>& in,
                               const std::vector<bool>& mask,
                               const std::string& which);

  casa::CountedPtr<Scantable> bin( const casa::CountedPtr<Scantable>& in,
                                   int width=5);
  casa::CountedPtr<Scantable>
    resample(const casa::CountedPtr<Scantable>& in,
             const std::string& method, float width);

  casa::CountedPtr<Scantable>
    smooth(const casa::CountedPtr<Scantable>& in, const std::string& kernel,
		      float width);

  casa::CountedPtr<Scantable>
    gainElevation(const casa::CountedPtr<Scantable>& in,
                  const casa::Vector<casa::Float>& coeffs,
                  const std::string& fileName,
		  const std::string& method);
  casa::CountedPtr<Scantable>
    convertFlux(const casa::CountedPtr<Scantable>& in, float d,
                float etaap, float jyperk);

  casa::CountedPtr<Scantable> opacity(const casa::CountedPtr<Scantable>& in,
                                      float tau);

  /// @todo frequency alignment

private:
  static imethod stringToIMethod(const std::string& in);
  static WeightType stringToWeight(const std::string& in);

  void scaleByVector(casa::Table& in,
                     const casa::Vector<casa::Float>& factor,
                     bool dotsys);

  void scaleFromAsciiTable(casa::Table& in, const std::string& filename,
                           const std::string& method,
                           const casa::Vector<casa::Float>& xout,
                           bool dotsys);

  void scaleFromTable(casa::Table& in, const casa::Table& table,
                      const std::string& method,
                      const casa::Vector<casa::Float>& xout, bool dotsys);

  void convertBrightnessUnits(casa::CountedPtr<Scantable>& in,
                              bool tokelvin, float cfac);

  casa::CountedPtr< Scantable >
    getScantable(const casa::CountedPtr< Scantable >& in, bool droprows);

  casa::MaskedArray<casa::Float>
    maskedArray( const casa::Vector<casa::Float>& s,
                 const casa::Vector<casa::uChar>& f );
  casa::Vector<casa::uChar>
    flagsFromMA(const casa::MaskedArray<casa::Float>& ma);

  bool insitu_;
};

}
#endif