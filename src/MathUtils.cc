//#---------------------------------------------------------------------------
//# MathUtilities.cc: General math operations
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

#include <aips/Arrays/Vector.h>
#include <aips/Arrays/VectorSTLIterator.h>

template <class T>
void hanning(Vector<T>& out, Vector<Bool>& outmask, 
	     const Vector<T>& in, const Vector<Bool>& mask, 
	     Bool relaxed, Bool ignoreOther) {

  Vector< Vector<T> > weights(8);
  Vector<Float> vals(3);
  vals = 0.0;weights[0] = vals;// FFF
  vals[0] = 1.0; vals[1] = 0.0; vals[2] = 0.0; weights[1] = vals;// TFF
  vals[0] = 0.0; vals[1] = 1.0; vals[2] = 0.0; weights[2] = vals;// FTF
  vals[0] = 1.0/3.0; vals[1] = 2.0/3.0; vals[2] = 0.0; weights[3] = vals;// TTF
  vals[0] = 0.0; vals[1] = 0.0; vals[2] = 1.0;weights[4] = vals;// FFT
  vals[0] = 0.5; vals[1] = 0.0; vals[2] = 0.5; weights[5] = vals;// TFT
  vals[0] = 0.0; vals[1] = 2.0/3.0; vals[2] = 1.0/3.0; weights[6] = vals;// FTT
  vals[0] = 0.25; vals[1] = 0.5; vals[2] = 0.25; weights[7] = vals;// TTT  
  // Chris' case
  Vector<Bool> weighted(8); 
  if (relaxed) {
    weighted = False;
    weighted[7] = True;

  } else {
    weighted = True;
    weighted[0] = False;
  }
  
  out.resize(in.nelements());
  outmask.resize(mask.nelements());

  // make special case for first and last
  /// ...here
  // loop from 1..n-2
  uInt i = 1;
  VectorSTLIterator<T> outit(out);
  outit++;
  VectorSTLIterator<Bool> outmit(outmask);outmit++;
  uInt m;Vector<T>* w;
  for (VectorSTLIterator<T> it = outit;it != out.end()-1;++it) {

    m = mask[i-1] + 2*mask[i] + 4*mask[i+1];
    w = &(weights[m]);
    if (weighted[m]) {
      (*it) = (*w)[0]*in[i-1] + (*w)[1]*in[i] + (*w)[2]*in[i+1];
      (*outmit) = True;
    } else { // mask it
      (*outmit) = False;
    }
    ++i;
    ++outmit;
  }
}
