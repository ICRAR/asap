//
// C++ Implementation: RowAccumulator
//
// Description:
//
//
// Author: Malte Marquarding <Malte.Marquarding@csiro.au>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <iostream>

#include <casa/iomanip.h>
#include <casa/Arrays/MaskArrMath.h>
#include <casa/Arrays/MaskArrLogi.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayLogical.h>
#include "RowAccumulator.h"

using namespace casa;
using namespace asap;

RowAccumulator::RowAccumulator(WeightType wt) : weightType_(wt), initialized_(False)
{
  reset();
}

RowAccumulator::~RowAccumulator()
{
}


void RowAccumulator::reset(const uInt size, const uInt tsysSize)
{
  Vector<Bool> m(size, True);

  spectrum_.setData(Vector<Float>(size, 0.0), Vector<Bool>(size, True));
  spectrumNoMask_.setData(Vector<Float>(size, 0.0), Vector<Bool>(size, True));

  n_.setData(Vector<uInt>(size, 0), Vector<Bool>(size, True));
  nNoMask_.setData(Vector<uInt>(size, 0), Vector<Bool>(size, True));

  weightSum_.setData(Vector<Float>(size, 0.0), Vector<Bool>(size, True));
  weightSumNoMask_.setData(Vector<Float>(size, 0.0), Vector<Bool>(size, True));

  tsysSum_.resize(tsysSize); tsysSum_=0.0;
  tsysSumNoMask_.resize(tsysSize); tsysSumNoMask_=0.0;

  intervalSum_ = 0.0;
  intervalSumNoMask_ = 0.0;

  timeSum_ = 0.0;
  timeSumNoMask_ = 0.0;
}

void RowAccumulator::initialize(const uInt size, const uInt tsysSize)
{
    reset(size, tsysSize);
    initialized_ = True;
}

void RowAccumulator::add(const Vector<Float>& v,
                         const Vector<Bool>& m,
                         const Vector<Float>& tsys,
                         const Double interval,
                         const Double time)
{
  uInt size = v.nelements();
  //if (size != m.nelements()) raiseError;
  if (!initialized_) initialize(size, tsys.nelements());

  addSpectrum(v, m, tsys, interval, time);
}

void RowAccumulator::addSpectrum(const Vector<Float>& v,
				 const Vector<Bool>& m,
				 const Vector<Float>& tsys,
				 const Double interval,
				 const Double time)
{
  doAddSpectrum(v, m, tsys, interval, time, False);
  doAddSpectrum(v, m, tsys, interval, time, True);  // CAS-2776
}

void RowAccumulator::doAddSpectrum(const Vector<Float>& v,
				   const Vector<Bool>& m,
				   const Vector<Float>& tsys,
				   const Double interval,
				   const Double time,
				   const Bool ignoreMask)
{
  Vector<Float> vUse = v.copy();
  Vector<Bool> mUse = m.copy();
  if (ignoreMask) mUse = !mUse;

  MaskedArray<Float> vadd(vUse, mUse);
  Float totalWeight = getTotalWeight(vadd, tsys, interval, time, ignoreMask);
  vadd *= totalWeight;
  MaskedArray<Float> wadd(Vector<Float>(mUse.nelements(), totalWeight), mUse);
  MaskedArray<uInt> inc(Vector<uInt>(mUse.nelements(), 1), mUse);

  if (ignoreMask) {
    spectrumNoMask_ += vadd;
    weightSumNoMask_ += wadd;
    nNoMask_ += inc;
  } else {
    spectrum_ += vadd;
    weightSum_ += wadd;
    n_ += inc;
  }

  //
  cout << "***" << endl;
  Vector<Float> spe = spectrum_.getArray();
  Vector<Float> spe0 = spectrumNoMask_.getArray();
  Vector<Float> wei = weightSum_.getArray();
  Vector<Float> wei0 = weightSumNoMask_.getArray();
  Vector<uInt>  n = n_.getArray();
  Vector<uInt>  n0 = nNoMask_.getArray();
  cout << "S__" << "[" << spe[0] << "][" << spe[1] << "][" << spe[2] << "][" << spe[3] << "][" << spe[4] << "][" << spe[5] << "][" << spe[6] << "][" << spe[7] << "][" << spe[8] << "][" << spe[9] << "][" << spe[10] << "][" << spe[11] << "][" << spe[12] << "]" << endl;
  cout << "S0_" << "[" << spe0[0] << "][" << spe0[1] << "][" << spe0[2] << "][" << spe0[3] << "][" << spe0[4] << "][" << spe0[5] << "][" << spe0[6] << "][" << spe0[7] << "][" << spe0[8] << "][" << spe0[9] << "][" << spe0[10] << "][" << spe0[11] << "][" << spe0[12] << "]" << endl;
  cout << "W__" << "[" << wei[0] << "][" << wei[1] << "][" << wei[2] << "][" << wei[3] << "][" << wei[4] << "][" << wei[5] << "][" << wei[6] << "][" << wei[7] << "][" << wei[8] << "][" << wei[9] << "][" << wei[10] << "][" << wei[11] << "][" << wei[12] << "]" << endl;
  cout << "W0_" << "[" << wei0[0] << "][" << wei0[1] << "][" << wei0[2] << "][" << wei0[3] << "][" << wei0[4] << "][" << wei0[5] << "][" << wei0[6] << "][" << wei0[7] << "][" << wei0[8] << "][" << wei0[9] << "][" << wei0[10] << "][" << wei0[11] << "][" << wei0[12] << "]" << endl;
  cout << "N__" << "[" << n[0] << "][" << n[1] << "][" << n[2] << "][" << n[3] << "][" << n[4] << "][" << n[5] << "][" << n[6] << "][" << n[7] << "][" << n[8] << "][" << n[9] << "][" << n[10] << "][" << n[11] << "][" << n[12] << "]" << endl;
  cout << "N0_" << "[" << n0[0] << "][" << n0[1] << "][" << n0[2] << "][" << n0[3] << "][" << n0[4] << "][" << n0[5] << "][" << n0[6] << "][" << n0[7] << "][" << n0[8] << "][" << n0[9] << "][" << n0[10] << "][" << n0[11] << "][" << n0[12] << "]" << endl;
  cout << "***" << endl;
  //
}

Float RowAccumulator::getTotalWeight(const MaskedArray<Float>& data,
				     const Vector<Float>& tsys,
				     const Double interval,
				     const Double time,
				     const Bool ignoreMask)
{
  Float totalWeight = 1.0;

  Vector<Bool> m = data.getMask();
  if (!allEQ(m, False)) {  // only add these if not everything masked
    totalWeight *= addTsys(tsys, ignoreMask);
    totalWeight *= addInterval(interval, ignoreMask);
    addTime(time, ignoreMask);
  }

  if (weightType_ == W_VAR) {
    Float fac = 1.0/variance(data);
    if (!ignoreMask && (m.nelements() == userMask_.nelements()))
      fac = 1.0/variance(data(userMask_));

    totalWeight *= fac;
  }

  return totalWeight;
}

Float RowAccumulator::addTsys(const Vector<Float>& v, Bool ignoreMask)
{
  // @fixme this assume tsys is the same for all channels

  Float w = 1.0;
  if (ignoreMask) {
    tsysSumNoMask_ += v[0];
  } else {
    tsysSum_ += v[0];
  }
  if ( weightType_ == W_TSYS  || weightType_ == W_TINTSYS ) {
    w /= (v[0]*v[0]);
  }
  return w;
}

void RowAccumulator::addTime(Double t, Bool ignoreMask)
{
  if (ignoreMask) {
    timeSumNoMask_ += t;
  } else {
    timeSum_ += t;
  }
}

Float RowAccumulator::addInterval(Double inter, Bool ignoreMask)
{
  Float w = 1.0;
  if (ignoreMask) {
    intervalSumNoMask_ += inter;
  } else {
    intervalSum_ += inter;
  }
  if (weightType_ == W_TINT || weightType_ == W_TINTSYS) {
    w /= Float(inter);
  }
  return w;
}

Vector<Float> RowAccumulator::getSpectrum() const
{
  return (spectrum_/weightSum_).getArray();
}

Double RowAccumulator::getTime() const
{
  return timeSum_/Float(max(n_));
}

Double RowAccumulator::getInterval() const
{
  return intervalSum_;
}

Vector<Bool> RowAccumulator::getMask() const
{
  // Return the "total" mask - False where no points have been accumulated.
  return (n_.getArray() > uInt(0));
}

Vector<Float> RowAccumulator::getTsys() const
{
  // @fixme this assumes tsys.nelements() == 1
  return tsysSum_/Float(max(n_));
}

void RowAccumulator::setUserMask(const Vector<Bool>& m)
{
  userMask_.resize();
  userMask_ = m;
}

// Added by TT  check the state of RowAccumulator
Bool RowAccumulator::state() const
{
  return initialized_;
}

void RowAccumulator::replaceNaN()
{
  Vector<Float> v = spectrum_.getArray();
  Vector<Float> w = weightSum_.getArray();
  Vector<Float> vRef = spectrumNoMask_.getArray();
  Vector<Float> wRef = weightSumNoMask_.getArray();

  //-------
  cout << "SB-";
  for (uInt i=0; i<13; ++i) {
    cout << "[" << v[i] << "]";
  }
  cout << endl;
  cout << "WB-";
  for (uInt i=0; i<13; ++i) {
    cout << "[" << w[i] << "]";
  }
  cout << endl;
  //-------

  for (uInt i = 0; i < v.size(); ++i) {
    if (w[i] == 0.0) {
      v[i] = vRef[i];
      w[i] = wRef[i];
    }
  }

  spectrum_.setData(v, Vector<Bool>(v.nelements(), True));
  weightSum_.setData(w, Vector<Bool>(w.nelements(), True));

  //-------
  cout << "SA-";
  for (uInt i=0; i<13; ++i) {
    cout << "[" << v[i] << "]";
  }
  cout << endl;
  cout << "WA-";
  for (uInt i=0; i<13; ++i) {
    cout << "[" << w[i] << "]";
  }
  cout << endl;
  //-------
}
