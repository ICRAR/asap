//#---------------------------------------------------------------------------
//# SDMathWrapper.cc: Wrapper classes to use CountedPtr
//#---------------------------------------------------------------------------
//# Copyright (C) 2004
//# ATNF
//#
//# This program is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by the Free
//# Software Foundation; either version 2 of the License, or (at your option)
//# any later version.
//#
//# This program is distributed in the hope that it will be useful, but
//# WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
//# Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with this program; if not, write to the Free Software Foundation, Inc.,
//# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning this software should be addressed as follows:
//#        Internet email: Malte.Marquarding@csiro.au
//#        Postal address: Malte Marquarding,
//#                        Australia Telescope National Facility,
//#                        P.O. Box 76,
//#                        Epping, NSW, 2121,
//#                        AUSTRALIA
//#
//# $Id:
//#---------------------------------------------------------------------------


#include "SDMathWrapper.h"
#include "SDMath.h"


using namespace asap;
using namespace casa;

SDMemTableWrapper SDMathWrapper::quotient(const SDMemTableWrapper& on,
                                          const SDMemTableWrapper& off,
                                          Bool preserveContinuum)
{
    SDMath sdm;
    Bool doTSys = True;
    return SDMemTableWrapper(sdm.binaryOperate(on.getCP(), off.getCP(), 
                             String("QUOTIENT"), preserveContinuum, doTSys));
}


SDMemTableWrapper SDMathWrapper::binaryOperate(const SDMemTableWrapper& left,
                                               const SDMemTableWrapper& right,
                                               const std::string& op, bool doTSys)
{
    SDMath sdm;
    return SDMemTableWrapper(sdm.binaryOperate(left.getCP(), right.getCP(), 
                                               String(op), False, Bool(doTSys)));
}


void SDMathWrapper::scaleInSitu(SDMemTableWrapper& in, float factor, bool doAll, bool doTSys)
{
  SDMemTable* pIn = in.getPtr();
  const uInt what = 0;
//
  SDMath sdm;
  SDMemTable* pOut = sdm.unaryOperate (*pIn, Float(factor), 
					Bool(doAll), what, Bool(doTSys));
  *pIn = *pOut;
   delete pOut;
}

SDMemTableWrapper SDMathWrapper::scale(const SDMemTableWrapper& in,
                                       float factor, bool doAll, bool doTSys)
{
  const CountedPtr<SDMemTable>& pIn = in.getCP();
  const uInt what = 0;
  SDMath sdm;
  return CountedPtr<SDMemTable>(sdm.unaryOperate(*pIn, Float(factor), Bool(doAll), 
                                                 what, Bool (doTSys)));
}



void SDMathWrapper::addInSitu(SDMemTableWrapper& in, float offset, bool doAll)
{
  SDMemTable* pIn = in.getPtr();
  const uInt what = 1;
//
  SDMath sdm;
  Bool doTSys = False;
  SDMemTable* pOut = sdm.unaryOperate (*pIn, Float(offset), 
					Bool(doAll), what, doTSys);
  *pIn = *pOut;
   delete pOut;
}

SDMemTableWrapper SDMathWrapper::add(const SDMemTableWrapper& in,
                                     float offset, bool doAll)
{
  const CountedPtr<SDMemTable>& pIn = in.getCP();
  const uInt what = 1;
  SDMath sdm;
  Bool doTSys = False;
  return CountedPtr<SDMemTable>(sdm.unaryOperate(*pIn, Float(offset),
						  Bool(doAll), what, doTSys));
}


void SDMathWrapper::smoothInSitu(SDMemTableWrapper& in, 
				 const std::string& kernel, float width,
				 bool doAll)
{
  SDMemTable* pIn = in.getPtr();
  SDMath sdm;
  SDMemTable* pOut = sdm.smooth(*pIn, String(kernel), 
				Float(width), Bool(doAll));
  *pIn = *pOut;
   delete pOut;
}


SDMemTableWrapper SDMathWrapper::smooth (const SDMemTableWrapper& in, 
					 const std::string& kernel, 
                                         float width, bool doAll)
{
  const CountedPtr<SDMemTable>& pIn = in.getCP();
  SDMath sdm;
  return CountedPtr<SDMemTable>(sdm.smooth(*pIn, String(kernel), 
					   Float(width), Bool(doAll)));
}



void SDMathWrapper::binInSitu(SDMemTableWrapper& in, int width)
{
  SDMemTable* pIn = in.getPtr();
  SDMath sdm;
  SDMemTable* pOut = sdm.bin (*pIn, Int(width));
  *pIn = *pOut;
   delete pOut;
}

SDMemTableWrapper SDMathWrapper::bin (const SDMemTableWrapper& in,
                                      int width)
{
  const CountedPtr<SDMemTable>& pIn = in.getCP();
  SDMath sdm;
  return CountedPtr<SDMemTable>(sdm.bin(*pIn, Int(width)));
}

void SDMathWrapper::resampleInSitu(SDMemTableWrapper& in, const std::string& method,
                                   float width)
{
  SDMemTable* pIn = in.getPtr();
  SDMath sdm;
  SDMemTable* pOut = sdm.resample(*pIn, String(method), Float(width));
  *pIn = *pOut;
   delete pOut;
}

SDMemTableWrapper SDMathWrapper::resample (const SDMemTableWrapper& in,
                                           const std::string& method, float width)
{
  const CountedPtr<SDMemTable>& pIn = in.getCP();
  SDMath sdm;
  return CountedPtr<SDMemTable>(sdm.resample(*pIn, String(method), Float(width)));
}


void SDMathWrapper::averagePolInSitu(SDMemTableWrapper& in, 
				     const std::vector<bool>& mask,
                                     const std::string& weightStr)
{
  SDMemTable* pIn = in.getPtr();
  SDMath sdm;
  Vector<Bool> tMask(mask);
  SDMemTable* pOut = sdm.averagePol (*pIn, tMask, String(weightStr));
  *pIn = *pOut;
   delete pOut;
}

SDMemTableWrapper SDMathWrapper::averagePol (const SDMemTableWrapper& in,
					     const std::vector<bool>& mask,
                                             const std::string& weightStr)
{
  const CountedPtr<SDMemTable>& pIn = in.getCP();
  SDMath sdm;
  Vector<Bool> tMask(mask);
  return CountedPtr<SDMemTable>(sdm.averagePol(*pIn, tMask, String(weightStr)));
}


std::vector<float> SDMathWrapper::statistic(const SDMemTableWrapper& in,
                                            const std::vector<bool>& mask, 
                                            const std::string& which, int row) 
{
  SDMath sdm;
  Vector<Bool> tMask(mask);
  return sdm.statistic(in.getCP(), tMask, String(which), Int(row));
}


void SDMathWrapper::convertFluxInSitu(SDMemTableWrapper& in, 
                                      float area, float eta, bool doAll)
{
  SDMemTable* pIn = in.getPtr();
  SDMath sdm;
  SDMemTable* pOut = sdm.convertFlux (*pIn, Float(area), Float(eta), Bool(doAll));
  *pIn = *pOut;
  delete pOut;
}


SDMemTableWrapper SDMathWrapper::convertFlux(const SDMemTableWrapper& in, 
                                             float area, float eta, bool doAll)
{
  const CountedPtr<SDMemTable>& pIn = in.getCP();
  SDMath sdm;
  return CountedPtr<SDMemTable>(sdm.convertFlux(*pIn, Float(area), Float(eta), Bool(doAll)));
}


void SDMathWrapper::gainElevationInSitu(SDMemTableWrapper& in, 
                                        const std::vector<float>& coeffs,
                                        const string& fileName,
                                        const string& method, bool doAll)
{
  SDMemTable* pIn = in.getPtr();
  Vector<Float> tCoeffs(coeffs);
  SDMath sdm;
  SDMemTable* pOut = sdm.gainElevation(*pIn, tCoeffs, String(fileName), 
                                       String(method), Bool(doAll));
  *pIn = *pOut;
  delete pOut;
}


SDMemTableWrapper SDMathWrapper::gainElevation(const SDMemTableWrapper& in, 
                                               const std::vector<float>& coeffs,
                                               const string& fileName, 
                                               const string& method, bool doAll)
{
  const CountedPtr<SDMemTable>& pIn = in.getCP();
  Vector<Float> tCoeffs(coeffs);
  SDMath sdm;
  return CountedPtr<SDMemTable>(sdm.gainElevation(*pIn, tCoeffs, String(fileName), 
                                                  String(method), Bool(doAll)));
}

void SDMathWrapper::frequencyAlignmentInSitu (SDMemTableWrapper& in, const std::string& refTime)
{
  SDMemTable* pIn = in.getPtr();
  SDMath sdm;
  SDMemTable* pOut = sdm.frequencyAlignment(*pIn, String(refTime));
  *pIn = *pOut;
  delete pOut;
}


SDMemTableWrapper SDMathWrapper::frequencyAlignment (const SDMemTableWrapper& in,
                                                     const std::string& refTime)
{
  const CountedPtr<SDMemTable>& pIn = in.getCP();
  SDMath sdm;
  return CountedPtr<SDMemTable>(sdm.frequencyAlignment(*pIn, String(refTime)));
}

void SDMathWrapper::opacityInSitu(SDMemTableWrapper& in, float tau, bool doAll)
{
  SDMemTable* pIn = in.getPtr();
  SDMath sdm;
  SDMemTable* pOut = sdm.opacity(*pIn, Float(tau), Bool(doAll));
  *pIn = *pOut;
  delete pOut;
}


SDMemTableWrapper SDMathWrapper::opacity(const SDMemTableWrapper& in, 
                                         float tau, bool doAll)
{
  const CountedPtr<SDMemTable>& pIn = in.getCP();
  SDMath sdm;
  return CountedPtr<SDMemTable>(sdm.opacity(*pIn, Float(tau), Bool(doAll)));
}

