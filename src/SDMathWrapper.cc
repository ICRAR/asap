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

using namespace asap;
using namespace casa;

SDMemTableWrapper SDMathWrapper::quotient(const SDMemTableWrapper& on,
                                          const SDMemTableWrapper& off)
{
    return SDMemTableWrapper(SDMath::quotient(on.getCP(),
                                             off.getCP()));
}


void SDMathWrapper::scaleInSitu(SDMemTableWrapper& in, casa::Float factor, casa::Bool doAll)
{
  SDMemTable* pIn = in.getPtr();
  const uInt what = 0;
  SDMemTable* pOut = SDMath::simpleOperate (*pIn, factor, doAll, what);
  *pIn = *pOut;
   delete pOut;
}

SDMemTableWrapper SDMathWrapper::scale(const SDMemTableWrapper& in,
                          casa::Float factor, casa::Bool doAll)
{
  const CountedPtr<SDMemTable>& pIn = in.getCP();
  const uInt what = 0;
  return CountedPtr<SDMemTable>(SDMath::simpleOperate(*pIn, factor, doAll, what));
}



void SDMathWrapper::addInSitu(SDMemTableWrapper& in, casa::Float offset, casa::Bool doAll)
{
  SDMemTable* pIn = in.getPtr();
  const uInt what = 1;
  SDMemTable* pOut = SDMath::simpleOperate (*pIn, offset, doAll, what);
  *pIn = *pOut;
   delete pOut;
}

SDMemTableWrapper SDMathWrapper::add(const SDMemTableWrapper& in,
                                     casa::Float offset, casa::Bool doAll)
{
  const CountedPtr<SDMemTable>& pIn = in.getCP();
  const uInt what = 1;
  return CountedPtr<SDMemTable>(SDMath::simpleOperate(*pIn, offset, doAll, what));
}


SDMemTableWrapper SDMathWrapper::hanning(const SDMemTableWrapper& in) 
{
  return SDMemTableWrapper(SDMath::hanning(in.getCP()));
}


void SDMathWrapper::binInSitu(SDMemTableWrapper& in, int width)
{
  SDMemTable* pIn = in.getPtr();
  SDMemTable* pOut = SDMath::bin (*pIn, Int(width));
  *pIn = *pOut;
   delete pOut;
}

SDMemTableWrapper SDMathWrapper::bin (const SDMemTableWrapper& in,
                                      int width)
{
  const CountedPtr<SDMemTable>& pIn = in.getCP();
  return CountedPtr<SDMemTable>(SDMath::bin(*pIn, Int(width)));
}


void SDMathWrapper::averagePolInSitu(SDMemTableWrapper& in, const std::vector<bool>& mask)
{
  SDMemTable* pIn = in.getPtr();
  SDMemTable* pOut = SDMath::averagePol (*pIn, mask);
  *pIn = *pOut;
   delete pOut;
}

SDMemTableWrapper SDMathWrapper::averagePol (const SDMemTableWrapper& in, const std::vector<bool>& mask)

{
  const CountedPtr<SDMemTable>& pIn = in.getCP();
  return CountedPtr<SDMemTable>(SDMath::averagePol(*pIn, mask));
}


std::vector<float> SDMathWrapper::statistic(const SDMemTableWrapper& in,
                                            const std::vector<bool>& mask, 
                                            const std::string& which) 
{
  return SDMath::statistic(in.getCP(), mask, which);
}
