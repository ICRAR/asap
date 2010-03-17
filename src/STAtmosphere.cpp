//#---------------------------------------------------------------------------
//# STAtmosphere.h: Model of atmospheric opacity
//#---------------------------------------------------------------------------
//# Copyright (C) 2004
//# ATNF
//#
//# The code is based on the Fortran code written by Bob Sault for MIRIAD.
//# Converted to C++ by Max Voronkov. This code uses a simple model of the
//# atmosphere and Liebe's model (1985) of the complex refractive index of
//# air.
//# 
//# The model of the atmosphere is one with an exponential fall-off in
//# the water vapour content (scale height of 1540 m) and a temperature lapse
//# rate of 6.5 mK/m. Otherwise the atmosphere obeys the ideal gas equation
//# and hydrostatic equilibrium.
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
//# $Id: STAtmosphere.h 1346 2007-04-26 03:24:41Z mar637 $
//#---------------------------------------------------------------------------

// own includes
#include "STAtmosphere.h"

// casa includes
#include <casa/Utilities/Assert.h>

using namespace casa;
using namespace asap;

/**
 * Default Constructor (apart from optional parameters).
 * The class set up this way will assume International Standard Atmosphere (ISA) conditions,
 * except for humidity. The latter is assumed to be 50%, which seems more realistic for 
 * Australian telescopes than 0%. 
 * @param[in] wvScale water vapour scale height (m), default is 1540m to match MIRIAD's model
 * @param[in] maxAlt maximum altitude of the model atmosphere (m), plane parallel layers are spread linearly up to
 *            this height, default is 10000m to match MIRIAD.
 * @param[in] nLayers number of plane parallel layers in the model (essentially for a numberical integration),
 *            default is 50 to match MIRIAD.
 **/
STAtmosphere::STAtmosphere(double wvScale, double maxAlt, size_t nLayers) :
   itsGndTemperature(288.), itsPressure(101325.), itsHumidity(0.5), 
   itsLapseRate(0.0065), itsWVScale(wvScale), itsMaxAlt(maxAlt),
   itsHeights(nLayers), itsTemperatures(nLayers), 
   itsDryPressures(nLayers), itsVapourPressures(nLayers) 
{
  recomputeAtmosphereModel();
}

/**
 * Constructor with explicitly given parameters of the atmosphere 
 * @param[in] temperature air temperature at the observatory (K)
 * @param[in] pressure air pressure at the observatory (Pascals)
 * @param[in] humidity air humidity at the observatory (fraction)
 * @param[in] lapseRate temperature lapse rate (K/m), default is 0.0065 K/m to match MIRIAD and ISA 
 * @param[in] wvScale water vapour scale height (m), default is 1540m to match MIRIAD's model
 * @param[in] maxAlt maximum altitude of the model atmosphere (m), plane parallel layers are spread linearly up to
 *            this height, default is 10000m to match MIRIAD.
 * @param[in] nLayers number of plane parallel layers in the model (essentially for a numberical integration),
 *            default is 50 to match MIRIAD.
 **/
STAtmosphere::STAtmosphere(double temperature, double pressure, double humidity, double lapseRate, 
               double wvScale, double maxAlt, size_t nLayers) :
   itsGndTemperature(temperature), itsPressure(pressure), itsHumidity(humidity), 
   itsLapseRate(lapseRate), itsWVScale(wvScale), itsMaxAlt(maxAlt),
   itsHeights(nLayers), itsTemperatures(nLayers), 
   itsDryPressures(nLayers), itsVapourPressures(nLayers) 
{
  recomputeAtmosphereModel();
}
               
/**
 * Set the new weather station data, recompute the model 
 * @param[in] temperature air temperature at the observatory (K)
 * @param[in] pressure air pressure at the observatory (Pascals)
 * @param[in] humidity air humidity at the observatory (fraction)
 **/
void STAtmosphere::setWeather(double temperature, double pressure, double humidity)
{
  itsGndTemperature = temperature;
  itsPressure = pressure;
  itsHumidity = humidity;
  recomputeAtmosphereModel();
}

/**
 * Build the atmosphere model based on exponential fall-off, ideal gas and hydrostatic
 * equilibrium. The model parameters are taken from the data members of this class.
 **/
void STAtmosphere::recomputeAtmosphereModel()
{
}
  
/**
 * Obtain the number of model layers, do consistency check that everything is
 * resized accordingly
 * @retrun number of model layers
 **/
size_t STAtmosphere::nLayers() const
{
  const size_t result = itsHeights.size();
  DebugAssert(itsTemperatures.size() == result, AipsError);
  DebugAssert(itsDryPressures.size() == result, AipsError);
  DebugAssert(itsVapourPressures.size() == result, AipsError);  
  return result;
}

