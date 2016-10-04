//
// C++ Interface: RowAccumulator
//
// Description:
//
//
// Author: Malte Marquarding <Malte.Marquarding@csiro.au>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPROWACCUMULATOR_H
#define ASAPROWACCUMULATOR_H

#include <math.h>
#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/MaskedArray.h>
#include "STDefs.h"

namespace asap {
/**
  * This class accumulates spectra and weights and returns the averaged data
  * @brief Class for averaging of spectra
  * @author Malte Marquarding
  * @date $Date:$
  * @version
  */
class RowAccumulator {

public:

  /**
   * Constructor taking a weight type as defined in @ref STDefs
   */
  explicit RowAccumulator(WeightType wt = asap::W_NONE);

 ~RowAccumulator();

  /**
    * add a new "row" to the accumulator
    * @param v the spectrum
    * @param m the mask for the spectrum
    * @param tsys the Tsys corresponing to the spectrum
    * @param interval the intergration time
    * @param time the time of the observation
    */
  void add(const casacore::Vector<casacore::Float>& v,
           const casacore::Vector<casacore::Bool>& m,
           const casacore::Vector<casacore::Float>& tsys,
           const casacore::Double interval,
           const casacore::Double time);
  /**
    * Also set a user mask which get combined with the individual masks
    * from the spectra
    * @param m a boolean mask of teh same length as the spectrum
    */
  void setUserMask(const casacore::Vector<casacore::Bool>& m);
  /**
    * Get the spectrum. Applies the normalisation (averaging)
    * @return the spectrum vector
    */
  casacore::Vector<casacore::Float> getSpectrum() const;
  /**
    * Get the Tsys. Applies the normalisation (averaging)
    * @return the Tsys vector
    */
  casacore::Vector<casacore::Float> getTsys() const;
  /**
    * Get the spectrum's mask. Applies the normalisation (averaging)
    * @return the mask vector
    */
  casacore::Vector<casacore::Bool> getMask() const;
  /**
    * Get the total interval.
    * @return the integration time
    */
  casacore::Double getInterval() const;
  /**
    * Get the time of the observation. Retrieves the "mean" time.
    * @return the integration time
    */
  casacore::Double getTime() const;
  /**
    * Reset the acummulator to the state at construction.
    */
  void reset(const casacore::uInt size=0, const casacore::uInt tsysSize=0);
  void initialize(const casacore::uInt size, const casacore::uInt tsysSize);
  /**
    * check the initialization state 
    */ 
  casacore::Bool state() const;
  /**
    * replace NaN values with (normal) values at the same channels in the given spetrum.
    * (CAS-2776; 2011/04/07 by Wataru Kawasaki)
    */
  void replaceNaN();

private:
  void addSpectrum(const casacore::Vector<casacore::Float>& v,
		   const casacore::Vector<casacore::Bool>& m,
		   const casacore::Vector<casacore::Float>& tsys,
		   const casacore::Double interval,
		   const casacore::Double time);
  void doAddSpectrum(const casacore::Vector<casacore::Float>& v,
		     const casacore::Vector<casacore::Bool>& m,
		     const casacore::Vector<casacore::Float>& tsys,
		     const casacore::Double interval,
		     const casacore::Double time,
		     const casacore::Bool inverseMask);
  void doAddSpectrum2(const casacore::Vector<casacore::Float>& v,
                      const casacore::Vector<casacore::Bool>& m,
                      const casacore::Vector<casacore::Float>& tsys,
                      const casacore::Double interval,
                      const casacore::Double time);
  casacore::Float getTotalWeight(const casacore::MaskedArray<casacore::Float>& data,
			     const casacore::Vector<casacore::Float>& tsys,
			     const casacore::Double interval,
			     const casacore::Double time,
			     const casacore::Bool inverseMask);
  casacore::Float addTsys(const casacore::Vector<casacore::Float>& v, casacore::Bool inverseMask);
  casacore::Float addInterval(casacore::Double inter, casacore::Bool inverseMask);
  void addTime(casacore::Double t, casacore::Bool inverseMask);

  WeightType weightType_;
  casacore::Bool initialized_;
  //these are Vectors
  casacore::MaskedArray<casacore::Float> spectrum_;
  casacore::MaskedArray<casacore::Float> weightSum_;
  casacore::MaskedArray<casacore::uInt> n_;

  //these three are used for normalise() (CAS-2776; 2011/04/07 by WK)
  casacore::MaskedArray<casacore::Float> spectrumNoMask_;
  casacore::MaskedArray<casacore::Float> weightSumNoMask_;
  casacore::MaskedArray<casacore::uInt> nNoMask_;

  casacore::Vector<casacore::Bool> userMask_;

  casacore::Vector<casacore::Float> tsysSum_, tsysSumNoMask_;
  casacore::Double timeSum_, timeSumNoMask_;
  casacore::Double intervalSum_, intervalSumNoMask_;
};

}
#endif
