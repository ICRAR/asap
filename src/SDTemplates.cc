//#---------------------------------------------------------------------------
//# SDTemplates.cc: explicit templates for aips++
//#---------------------------------------------------------------------------
//# Copyright (C) 2004
//# Malte Marquarding, ATNF
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
#include <casa/Containers/Block.h>
#include <casa/Exceptions/Error.cc>
#include <casa/Exceptions/Error2.cc>
#include <casa/Utilities/CountedPtr.cc>
#include <casa/Utilities/CountedPtr2.cc>
#include <casa/Utilities/Copy.cc>
#include "SDMemTable.h"

namespace asap {
  template class CountedConstPtr<SDMemTable>;
  template class CountedPtr<SDMemTable>;
  template class PtrRep<SDMemTable>;
  template class SimpleCountedConstPtr<SDMemTable>;
  template class SimpleCountedPtr<SDMemTable>;
  template class Block<CountedPtr<SDMemTable> >;
}
template void objcopy(CountedPtr<asap::SDMemTable> *, CountedPtr<asap::SDMemTable> const *, uInt);
template void objmove(CountedPtr<asap::SDMemTable> *, CountedPtr<asap::SDMemTable> const *, uInt);
template void objset(CountedPtr<asap::SDMemTable> *, CountedPtr<asap::SDMemTable>, uInt);

#include <casa/Arrays/ArrayLogical.cc>
#include <casa/Arrays/ArrayMath.cc>
#include <casa/Arrays/MaskArrMath.cc>
#include <casa/Arrays/Array.h>
#include <scimath/Functionals/CompiledFunction.cc>
#include <scimath/Functionals/CompiledParam.cc>
#include <scimath/Mathematics/AutoDiff.h>
#include <scimath/Mathematics/AutoDiffMath.h>
#include <casa/Arrays/Vector2.cc>
#include <images/Images/ImageUtilities2.cc>
#include <casa/Utilities/PtrHolder.cc>
#include <lattices/Lattices/Lattice.h>
#include "MathUtils.cc"

template void convertArray(Array<Bool> &, Array<uChar> const &);
template void convertArray(Array<uChar> &, Array<Bool> const &);
template LogicalArray operator!=(Array<Float> const &, Float const &);
template LogicalArray operator==(Array<Float> const &, Float const &);
template LogicalArray operator>(Array<Float> const &, Float const &);
template LogicalArray operator>=(Array<Float> const &, Float const &);
template Array<Float>& operator/=(Array<Float>&, MaskedArray<Float> const&);
template MaskedArray<Float> const& operator*=(MaskedArray<Float> const&, Float const&);
template MaskedArray<Float> operator-(MaskedArray<Float> const&, MaskedArray<Float> const&);
template MaskedArray<Float> const& operator/=(MaskedArray<Float> const&, Float const&);
template Float stddev(MaskedArray<Float> const&);
template class CompiledFunction<AutoDiff<Float> >;
template class CompiledParam<AutoDiff<Float> >;
template Vector<Bool>::Vector(const vector<bool> &);
template void Array<float>::tovector(vector<float> &) const;
template void hanning(Vector<Float>&, Vector<Bool>&, 
		      const Vector<Float>&, const Vector<Bool>&, Bool, Bool);
template void ImageUtilities::bin(MaskedArray<float>&, Coordinate&, MaskedArray<float> const&, Coordinate const&, uInt, uInt);
template class PtrHolder<Lattice<Float> >;

