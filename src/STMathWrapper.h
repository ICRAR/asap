//
// C++ Interface: STMathWrapper
//
// Description:
//
//
// Author: Malte Marquarding <Malte.Marquarding@csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSTMATHWRAPPER_H
#define ASAPSTMATHWRAPPER_H

#include <vector>
#include <string>

#include <casa/Utilities/CountedPtr.h>

#include "STMath.h"
#include "Scantable.h"
#include "ScantableWrapper.h"

namespace asap {

/**
Wrapper class to handle ScantableWrapper

@author Malte Marquarding
*/
class STMathWrapper : public STMath {
public:
  STMathWrapper() {;}
  STMathWrapper(bool insitu) : STMath(insitu) {;}

  virtual ~STMathWrapper() {;}

  ScantableWrapper
    average( const std::vector<ScantableWrapper>& in,
             const std::vector<bool>& mask,
             const std::string& weight,
             const std::string& avmode )
  {
    std::vector<casa::CountedPtr<Scantable> > sts;
    for (unsigned int i=0; i<in.size(); ++i) sts.push_back(in[i].getCP());
    return ScantableWrapper(STMath::average(sts, mask, weight, avmode));
  }

  ScantableWrapper
    averageChannel( const ScantableWrapper& in,
                    const std::string& mode = "MEDIAN",
                    const std::string& avmode = "NONE")
  {
    return ScantableWrapper(STMath::averageChannel(in.getCP(), mode, avmode));
  }

  ScantableWrapper
    averagePolarisations( const ScantableWrapper& in,
                          const std::vector<bool>& mask,
                          const std::string& weight)
  { return ScantableWrapper(STMath::averagePolarisations(in.getCP(),mask, weight));}

  ScantableWrapper
    averageBeams( const ScantableWrapper& in,
                  const std::vector<bool>& mask,
                  const std::string& weight)

  { return ScantableWrapper(STMath::averageBeams(in.getCP(),mask, weight));}

  ScantableWrapper
    unaryOperate( const ScantableWrapper& in, float val,
                  const std::string& mode, bool tsys=false )
  { return ScantableWrapper(STMath::unaryOperate(in.getCP(), val, mode, tsys)); }

  ScantableWrapper autoQuotient( const ScantableWrapper& in,
                                 const std::string& mode = "NEAREST",
                                 bool preserve = true )
  { return ScantableWrapper(STMath::autoQuotient(in.getCP(), mode, preserve)); }

  ScantableWrapper quotient( const ScantableWrapper& on,
                             const ScantableWrapper& off,
                             bool preserve = true )
  { return ScantableWrapper( STMath::quotient( on.getCP(), off.getCP(),
                                               preserve ) ); }

  ScantableWrapper
    freqSwitch( const ScantableWrapper& in )
  { return ScantableWrapper(STMath::freqSwitch(in.getCP())); }

  std::vector<float> statistic(const ScantableWrapper& in,
                               const std::vector<bool>& mask,
                               const std::string& which)
  { return STMath::statistic(in.getCP(), mask, which); }

  ScantableWrapper bin( const ScantableWrapper& in, int width=5)
  { return ScantableWrapper(STMath::bin(in.getCP(), width)); }

  ScantableWrapper
    resample(const ScantableWrapper& in,
             const std::string& method, float width)
  { return ScantableWrapper(STMath::resample(in.getCP(), method, width)); }

  ScantableWrapper
    smooth(const ScantableWrapper& in, const std::string& kernel, float width)
  { return ScantableWrapper(STMath::smooth(in.getCP(), kernel, width)); }

  ScantableWrapper
    gainElevation(const ScantableWrapper& in,
                  const std::vector<float>& coeff,
                  const std::string& filename,
                  const std::string& method)

  { return
      ScantableWrapper(STMath::gainElevation(in.getCP(), coeff, filename, method)); }

  ScantableWrapper
    convertFlux(const ScantableWrapper& in, float d,
                float etaap, float jyperk)
  { return ScantableWrapper(STMath::convertFlux(in.getCP(), d, etaap, jyperk)); }

  ScantableWrapper opacity(const ScantableWrapper& in,
                                      float tau)
  { return ScantableWrapper(STMath::opacity(in.getCP(), tau)); }

  ScantableWrapper
    merge(const std::vector<ScantableWrapper >& in)

  {
    std::vector<casa::CountedPtr<Scantable> > sts;
    for (unsigned int i=0; i<in.size(); ++i) sts.push_back(in[i].getCP());
    return ScantableWrapper(STMath::merge(sts)); }

  ScantableWrapper rotateXYPhase( const ScantableWrapper& in, float angle)
  { return ScantableWrapper(STMath::rotateXYPhase(in.getCP(), angle)); }

  ScantableWrapper rotateLinPolPhase( const ScantableWrapper& in, float angle)
  { return ScantableWrapper(STMath::rotateLinPolPhase(in.getCP(), angle)); }

  ScantableWrapper invertPhase( const ScantableWrapper& in )
  { return ScantableWrapper(STMath::invertPhase(in.getCP())); }

  ScantableWrapper swapPolarisations( const ScantableWrapper& in )
  { return ScantableWrapper(STMath::swapPolarisations(in.getCP())); }

  ScantableWrapper frequencyAlign( const ScantableWrapper& in,
                                   const std::string& refTime,
                                   const std::string& method  )
  { return ScantableWrapper(STMath::frequencyAlign(in.getCP())); }

  ScantableWrapper convertPolarisation( const ScantableWrapper& in,
                                        const std::string& newtype )
  { return ScantableWrapper(STMath::convertPolarisation(in.getCP(),newtype)); }

  ScantableWrapper mxExtract( const ScantableWrapper& in,
                              const std::string& scantype="on" )
  { return ScantableWrapper(STMath::mxExtract(in.getCP(),scantype)); }

  ScantableWrapper lagFlag( const ScantableWrapper& in,
                            double frequency, int width )
  { return ScantableWrapper(STMath::lagFlag(in.getCP(), frequency, width)); }

};

}

#endif
