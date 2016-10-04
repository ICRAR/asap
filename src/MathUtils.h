//#---------------------------------------------------------------------------
//# MathUtilities.h: General math operations
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
//# $Id$
//#---------------------------------------------------------------------------
#ifndef MATHUTILS_H
#define MATHUTILS_H

#include <string>
#include <vector>
#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <casa/BasicSL/String.h>
#include <casa/Arrays/IPosition.h>

namespace mathutil {

// Hanning smoothing
/**
 * Hanning smooth a masked vector
 * @param out the smoothed vector
 * @param outmask  the smoothed mask
 * @param in the input vector
 * @param mask the input mask
 * @param relaxed a weighting scheme
 * @param ignoreOther drop every second channel (NYI)
 */
void hanning(casacore::Vector<casacore::Float>& out,
	     casacore::Vector<casacore::Bool>& outmask,
             const casacore::Vector<casacore::Float>& in,
	     const casacore::Vector<casacore::Bool>& mask,
             casacore::Bool relaxed=casacore::False,
             casacore::Bool ignoreOther=casacore::False);

/**
 * Apply a running median to  a masked vector.
 * Edge solution:  The first and last hwidth channels will be replicated
 * from the first/last value from a full window.
 * @param out the smoothed vector
 * @param outmask  the smoothed mask
 * @param in the input vector
 * @param mask the input mask
 * @param hwidth half-width of the smoothing window
 */
void runningMedian(casacore::Vector<casacore::Float>& out,
                   casacore::Vector<casacore::Bool>& outflag,
                   const casacore::Vector<casacore::Float>& in,
                   const casacore::Vector<casacore::Bool>& flag,
                   float hwidth);

void polyfit(casacore::Vector<casacore::Float>& out,
             casacore::Vector<casacore::Bool>& outmask,
             const casacore::Vector<casacore::Float>& in,
             const casacore::Vector<casacore::Bool>& mask,
             float hwidth, int order);

// Generate specified statistic
float statistics(const casacore::String& which,
                 const casacore::MaskedArray<casacore::Float>& data);

// Return a position of min or max value
casacore::IPosition minMaxPos(const casacore::String& which,
                 const casacore::MaskedArray<casacore::Float>& data);

// Replace masked value by zero
void replaceMaskByZero(casacore::Vector<casacore::Float>& data,
                       const casacore::Vector<casacore::Bool>& mask);

/**
 * Convert casa implementations to stl
 * @param in casa string
 * @return a std vector of std strings
 */
std::vector<std::string> tovectorstring(const casacore::Vector<casacore::String>& in);

/**
 * convert stl implementations to casa versions
 * @param in
 * @return
 */
casacore::Vector<casacore::String> toVectorString(const std::vector<std::string>& in);

void doZeroOrderInterpolation(casacore::Vector<casacore::Float>& data,
			      std::vector<bool>& mask);

/**
 * RA nomalization: n*2pi rotation if necessary
 **/
void rotateRA( const casacore::Vector<casacore::Double> &in,
               casacore::Vector<casacore::Double> &out ) ;
void rotateRA( casacore::Vector<casacore::Double> &v ) ;

/**
 * tool to record current time stamp
 **/
double gettimeofday_sec() ;

}

#endif
