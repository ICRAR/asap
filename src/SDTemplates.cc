//#---------------------------------------------------------------------------
//# SDTemplates.cc: explicit templates for aips++
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
#include "SDMemTable.h"

#include <casa/aips.h>
#include <casa/namespace.h>
#include <casa/Containers/Block.h>
#include <casa/Exceptions/Error.cc>
#include <casa/Exceptions/Error2.cc>
#include <casa/Utilities/CountedPtr.cc>
#include <casa/Utilities/CountedPtr2.cc>
#include <casa/Utilities/Copy.cc>

namespace asap {
  template class casa::CountedConstPtr<SDMemTable>;
  template class casa::CountedPtr<SDMemTable>;
  template class casa::PtrRep<SDMemTable>;
  template class casa::SimpleCountedConstPtr<SDMemTable>;
  template class casa::SimpleCountedPtr<SDMemTable>;
  template class casa::Block<CountedPtr<SDMemTable> >;
}

template void objcopy<CountedPtr<asap::SDMemTable> >(CountedPtr<asap::SDMemTable> *, CountedPtr<asap::SDMemTable> const *, uInt);
template void objmove<CountedPtr<asap::SDMemTable> >(CountedPtr<asap::SDMemTable> *, CountedPtr<asap::SDMemTable> const *, uInt);
template void objset<CountedPtr<asap::SDMemTable> >(CountedPtr<asap::SDMemTable> *, CountedPtr<asap::SDMemTable>, uInt);

#include <casa/Arrays/ArrayLogical.cc>
#include <casa/Arrays/ArrayMath.cc>
#include <casa/Arrays/MaskArrMath.cc>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Vector2.cc>
#include <casa/Utilities/PtrHolder.cc>
#include <casa/Utilities/BinarySearch.cc>
#include <coordinates/Coordinates/VelocityAligner.cc>
#include <coordinates/Coordinates/FrequencyAligner.cc>
#include <lattices/Lattices/Lattice.h>
#include <lattices/Lattices/LatticeUtilities.cc>
#include <scimath/Functionals/CompiledFunction.cc>
#include <scimath/Functionals/CompiledParam.cc>
#include <scimath/Mathematics/AutoDiff.h>
#include <scimath/Mathematics/AutoDiffMath.h>
#include <scimath/Mathematics/InterpolateArray1D.cc>

template void convertArray<Bool, uChar>(Array<Bool> &, Array<uChar> const &);
template void convertArray<uChar, Bool>(Array<uChar> &, Array<Bool> const &);
template LogicalArray operator!=<Float>(Array<Float> const &, Float const &);
template LogicalArray operator==<Float>(Array<Float> const &, Float const &);
template LogicalArray operator><Float>(Array<Float> const &, Float const &);
template LogicalArray operator>=<Float>(Array<Float> const &, Float const &);
template Array<Float>& operator/=<Float>(Array<Float>&, MaskedArray<Float> const&);
template MaskedArray<Float> const& operator*=<Float>(MaskedArray<Float> const&, Float const&);
template MaskedArray<Float> operator+<Float>(MaskedArray<Float> const&, MaskedArray<Float> const&);
template MaskedArray<Float> operator-<Float>(MaskedArray<Float> const&, MaskedArray<Float> const&);
template MaskedArray<Float> operator*<Float>(MaskedArray<Float> const&, MaskedArray<Float> const&);
template MaskedArray<Float> operator/<Float>(MaskedArray<Float> const&, MaskedArray<Float> const&);
template MaskedArray<Float> const& operator/=<Float>(MaskedArray<Float> const&, Float const&);
template MaskedArray<Float> operator-<Float>(MaskedArray<Float> const&, Array<Float> const&);
template MaskedArray<Float> operator*<Float>(Array<Float> const&, MaskedArray<Float> const&);
template Float stddev<Float>(MaskedArray<Float> const&);
template Float median<Float>(MaskedArray<Float> const&);
template Float sumsquares<Float>(MaskedArray<Float> const&);
template Float avdev<Float>(MaskedArray<Float> const&);
template class CompiledFunction<AutoDiff<Float> >;
template class CompiledParam<AutoDiff<Float> >;
template Vector<Bool>::Vector(const vector<bool> &);
template Vector<Float>::Vector(const vector<float> &);
template Vector<Double>::Vector(const vector<double> &);
template void Array<float>::tovector(vector<float> &) const;
template void Array<Bool>::tovector(vector<bool> &) const;
template void Array<double>::tovector(vector<double> &) const;
template void LatticeUtilities::bin(MaskedArray<float>&, MaskedArray<float> const&, uInt, uInt);
template class PtrHolder<Lattice<Float> >;
template class VelocityAligner<Float>;
template class FrequencyAligner<Float>;
template class InterpolateArray1D<Double, Float>;
template Int binarySearchBrackets<Vector<Double>,Double>(Bool &, Vector<Double> const &, Double const &, uInt, Int);
//
#include "MathUtils2.cc"
namespace mathutil {
  template void hanning(Vector<Float>&, Vector<Bool>&,
			const Vector<Float>&, 
			const Vector<Bool>&, 
			Bool, Bool);
  template uInt addEntry(Vector<uInt>&, uInt);
}
