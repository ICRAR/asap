//#---------------------------------------------------------------------------
//# STAttr.h: Return known attributes about telescopes
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
#ifndef STATTR_H
#define STATTR_H

#include <string>

#include <casa/aips.h>
#include <measures/Measures/MEpoch.h>
#include <casa/Utilities/CountedPtr.h>

#include "STDefs.h"

namespace casacore {
  template<class T> class Vector;
}

namespace asap {

class STAttr {

 public:

  // Constructor
  STAttr();

  // Destructor
  virtual ~STAttr();

  // Copy Constructor (copy semantics)
  STAttr(const STAttr& other);

  // Assignment  (copy semantics)
  STAttr &operator=(const STAttr& other);

  // Telescope diameter(m). Throws exception if unknown.
  casacore::Float diameter(Instrument inst) const;

  // Beam efficiency.  Frequency in Hz.  Returns 1 if unknown.
  casacore::Vector<casacore::Float>
  beamEfficiency(Instrument instr, const casacore::MEpoch& dateObs,
                 const casacore::Vector<casacore::Float>& freqs) const;

  // Aperture efficiency. Frequency in Hz.  Returns 1 if unknown.
  casacore::Vector<casacore::Float>
  apertureEfficiency(Instrument instr,
                     const casacore::MEpoch& dateObs,
                     const casacore::Vector<casacore::Float>& freqs) const;

  // Find factor to convert Jy -> K for this telescope, date of
  // observation and frequency (Hz)
  casacore::Vector<casacore::Float> JyPerK(Instrument instr,
                                   const casacore::MEpoch& dateObs,
                                   const casacore::Vector<casacore::Float>& freqs)
    const;

  // Gain Elevation polynomial correction coefficients (elevation in
  // degrees) Returns length 0 if not known.
   casacore::Vector<casacore::Float> gainElevationPoly(Instrument instr) const;

  // Find feed polarization type of feeds. In future this needs to come
  // from the data themselves
  std::string feedPolType(Instrument) const;

  // Helper function to check instrument (antenna) name and give enum
  static Instrument convertInstrument(const casacore::String& instrument,
                                      casacore::Bool throwIt);

  // Helper function.  Finds factor to convert K -> Jy. Provide
  // aperture efficiency and dish geometric diameter (m)
  static casacore::Float findJyPerK(casacore::Float etaAp, casacore::Float D);

private:

  // Static data
  casacore::Vector<casacore::Float> MopEtaBeamX_;     // Beam efficiency
  casacore::Vector<casacore::Float> MopEtaBeam2003Y_;
  casacore::Vector<casacore::Float> MopEtaBeam2004Y_;

  casacore::Vector<casacore::Float> MopEtaApX_;       // Aperture efficiency
  casacore::Vector<casacore::Float> MopEtaAp2004Y_;

  casacore::Vector<casacore::Float> TidEtaApX_;       // Aperture efficiency
  casacore::Vector<casacore::Float> TidEtaApY_;

  casacore::Vector<casacore::Float> TidGainElPoly_;   // Gain-el correction poly coeffs
  casacore::Vector<casacore::Float> ParkesGainElPoly_;// K-band Gain-el correction
                                              //poly coeffs

  // Init private data
  void initData();

  // Linear interpolation
  casacore::Vector<casacore::Float> interp(const casacore::Vector<casacore::Float>& xOut,
                                   const casacore::Vector<casacore::Float>& xIn,
                                   const casacore::Vector<casacore::Float>& yIn) const;


};

} // namespace

#endif
