//#---------------------------------------------------------------------------
//# SDAttr.h: Return known attributes about telescopes
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
#ifndef SDATTR_H
#define SDATTR_H

#include "SDDefs.h"
#include <casa/aips.h>
#include <measures/Measures/MEpoch.h>

template<class T> class casa::Vector;
//class casa::MEpoch;



namespace asap {

class SDAttr {

 public:

// Constructor
   SDAttr();

// Destructor
   ~SDAttr();

// Copy Constructor (copy semantics)
   SDAttr (const SDAttr& other);

// Assignment  (copy semantics)
   SDAttr &operator=(const SDAttr& other);

// Telescope diameter (m). Throws exception if unknown.
   casa::Float diameter (Instrument inst) const;

// Beam efficiency.  Frequency in Hz.  Returns 1 if unknown.
   casa::Vector<casa::Float> beamEfficiency (Instrument, const casa::MEpoch& dateObs, 
                                             const casa::Vector<casa::Float>& freqs) const;
 
// Aperture efficiency. Frequency in Hz.  Returns 1 if unknown.
   casa::Vector<casa::Float> apertureEfficiency (Instrument, const casa::MEpoch& dateObs, 
                                                 const casa::Vector<casa::Float>& freqs) const;

// Find factor to convert Jy -> K for this telescope, date of observation and frequency (Hz)
   casa::Vector<casa::Float> JyPerK (Instrument, const casa::MEpoch& dateObs, 
                                     const casa::Vector<casa::Float>& freqs) const;

// Factor to convert K -> Jy. Provide aperture efficiency and dish geometric diameter (m)
   static casa::Float findJyPerKFac (casa::Float etaAp, casa::Float D);

 private:

// Static data

   casa::Vector<casa::Float> MopEtaBeamX_;
   casa::Vector<casa::Float> MopEtaBeam2003Y_;
   casa::Vector<casa::Float> MopEtaBeam2004Y_;
   casa::Vector<casa::Float> MopEtaApX_;
   casa::Vector<casa::Float> MopEtaAp2004Y_;

// Init private data
   void init();

// Linear interpolation
   casa::Vector<casa::Float> interp (const casa::Vector<casa::Float>& xOut, const casa::Vector<casa::Float>& xIn, 
                                     const casa::Vector<casa::Float>& yIn) const;
};

} // namespace

#endif
