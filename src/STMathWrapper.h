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
             const std::string& avmode,
             bool alignfreq)
  {
    std::vector<casa::CountedPtr<Scantable> > sts;
    for (int i=0; i<in.size(); ++i) sts.push_back(in[i].getCP());
    return ScantableWrapper(STMath::average(sts, mask, weight, avmode, alignfreq));
  }

  ScantableWrapper
    unaryOperate( const ScantableWrapper& in, float val,
                  const std::string& mode, bool tsys=false )
  { return ScantableWrapper(STMath::unaryOperate(in.getCP(), val, mode, tsys)); }

  ScantableWrapper quotient( const ScantableWrapper& in,
                             const std::string& mode = "NEAREST",
                             bool preserve = true )
  { return ScantableWrapper(STMath::quotient(in.getCP(), mode, preserve)); }

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
    for (int i=0; i<in.size(); ++i) sts.push_back(in[i].getCP());
    return ScantableWrapper(STMath::merge(sts)); }

};

}

#endif