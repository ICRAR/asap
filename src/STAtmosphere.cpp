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
#include <casa/Quanta.h>

// std includes
#include <cmath>

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
   itsGndTemperature(288.), itsGndPressure(101325.), itsGndHumidity(0.5), 
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
   itsGndTemperature(temperature), itsGndPressure(pressure), itsGndHumidity(humidity), 
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
  itsGndPressure = pressure;
  itsGndHumidity = humidity;
  recomputeAtmosphereModel();
}

/**
 * Build the atmosphere model based on exponential fall-off, ideal gas and hydrostatic
 * equilibrium. The model parameters are taken from the data members of this class.
 **/
void STAtmosphere::recomputeAtmosphereModel()
{
  AlwaysAssert(itsGndTemperature > 0, AipsError);
  AlwaysAssert(itsGndPressure > 0., AipsError);
  AlwaysAssert((itsGndHumidity >= 0.) && (itsGndHumidity<=1.), AipsError);
  AlwaysAssert(itsMaxAlt > 0., AipsError);
  AlwaysAssert(itsWVScale > 0., AipsError);
  
  const double heightStep = itsMaxAlt/double(nLayers());
  // molar mass of the air
  const double M = 28.96e-3;
  // free-fall acceleration
  const double g = 9.81;
  const double wvGndSaturationPressure = wvSaturationPressure(itsGndTemperature);
  for (size_t layer = 0; layer < nLayers(); ++layer) {
       const double height = double(layer)*heightStep;
       itsHeights[layer] = height;
       itsTemperatures[layer] = itsGndTemperature/(1.+itsLapseRate*height/itsGndTemperature);
       const double pressure = itsGndPressure * exp(-M*g/(QC::R.get().getValue()*itsGndTemperature)*
                   (height+0.5*itsLapseRate*height*height/itsGndTemperature));
       itsVapourPressures[layer] = casa::min(itsGndHumidity*exp(-height/itsWVScale)*wvGndSaturationPressure,
                                             wvSaturationPressure(itsTemperatures[layer]));
       itsDryPressures[layer] = pressure - itsVapourPressures[layer];                                      
  }
}
  
/**
 * Obtain the number of model layers, do consistency check that everything is
 * resized accordingly
 * @retrun number of model layers
 **/
size_t STAtmosphere::nLayers() const
{
  const size_t result = itsHeights.size();
  DebugAssert(result > 0, AipsError);
  DebugAssert(itsTemperatures.size() == result, AipsError);
  DebugAssert(itsDryPressures.size() == result, AipsError);
  DebugAssert(itsVapourPressures.size() == result, AipsError);  
  return result;
}

/**
 * Determine the saturation pressure of water vapour for the given temperature.
 *
 * Reference:
 * Waters, Refraction effects in the neutral atmosphere. Methods of
 * Experimental Physics, vol 12B, p 186-200 (1976).
 *   
 * @param[in] temperature temperature in K
 * @return vapour saturation pressure (Pascals) 
 **/
double STAtmosphere::wvSaturationPressure(double temperature)
{
  if (temperature > 215.) {
      return 0.;
  }
  const double theta = 300.0/temperature;
  return 1e5/(41.51/std::pow(theta,5)*std::pow(10.,9.834*theta-10.0));
}

/**
 * Compute the complex refractivity of the dry components of the atmosphere
 * (oxygen lines) at the given frequency.
 * @param[in] freq frequency (Hz)
 * @param[in] temperature air temperature (K)
 * @param[in] pDry partial pressure of dry components (Pascals)
 * @param[in] pVapour partial pressure of water vapour (Pascals)
 * @return complex refractivity
 **/
std::complex<double> STAtmosphere::dryRefractivity(double freq, double temperature, 
                     double pDry, double pVapour)
{
  // the number of parameters per atmospheric line and the number of lines taken into account
  const size_t nLineParams = 7;
  const size_t nLines = 48;
  // actual tabulated values
  const double lines[nLines][nLineParams] = 
    {{49.452379,    0.12E-6, 11.830,  8.40E-3, 0.0,  5.60E-3,  1.7},
     {49.962257,    0.34E-6, 10.720,  8.50E-3, 0.0,  5.60E-3,  1.7},
     {50.474238,    0.94E-6,  9.690,  8.60E-3, 0.0,  5.60E-3,  1.7},
     {50.987748,    2.46E-6,  8.690,  8.70E-3, 0.0,  5.50E-3,  1.7},
     {51.503350,    6.08E-6,  7.740,  8.90E-3, 0.0,  5.60E-3,  1.8},
     {52.021409,   14.14E-6,  6.840,  9.20E-3, 0.0,  5.50E-3,  1.8},
     {52.542393,   31.02E-6,  6.000,  9.40E-3, 0.0,  5.70E-3,  1.8},
     {53.066906,   64.10E-6,  5.220,  9.70E-3, 0.0,  5.30E-3,  1.9},
     {53.595748,  124.70E-6,  4.480, 10.00E-3, 0.0,  5.40E-3,  1.8},
     {54.129999,  228.00E-6,  3.810, 10.20E-3, 0.0,  4.80E-3,  2.0},
     {54.671157,  391.80E-6,  3.190, 10.50E-3, 0.0,  4.80E-3,  1.9},
     {55.221365,  631.60E-6,  2.620, 10.79E-3, 0.0,  4.17E-3,  2.1},
     {55.783800,  953.50E-6,  2.115, 11.10E-3, 0.0,  3.75E-3,  2.1},
     {56.264777,  548.90E-6,  0.010, 16.46E-3, 0.0,  7.74E-3,  0.9},
     {56.363387, 1344.00E-6,  1.655, 11.44E-3, 0.0,  2.97E-3,  2.3},
     {56.968180, 1763.00E-6,  1.255, 11.81E-3, 0.0,  2.12E-3,  2.5},
     {57.612481, 2141.00E-6,  0.910, 12.21E-3, 0.0,  0.94E-3,  3.7},
     {58.323874, 2386.00E-6,  0.621, 12.66E-3, 0.0, -0.55E-3, -3.1},
     {58.446589, 1457.00E-6,  0.079, 14.49E-3, 0.0,  5.97E-3,  0.8},
     {59.164204, 2404.00E-6,  0.386, 13.19E-3, 0.0, -2.44E-3,  0.1},
     {59.590982, 2112.00E-6,  0.207, 13.60E-3, 0.0,  3.44E-3,  0.5},
     {60.306057, 2124.00E-6,  0.207, 13.82E-3, 0.0, -4.13E-3,  0.7},
     {60.434775, 2461.00E-6,  0.386, 12.97E-3, 0.0,  1.32E-3, -1.0},
     {61.150558, 2504.00E-6,  0.621, 12.48E-3, 0.0, -0.36E-3,  5.8},
     {61.800152, 2298.00E-6,  0.910, 12.07E-3, 0.0, -1.59E-3,  2.9},
     {62.411212, 1933.00E-6,  1.255, 11.71E-3, 0.0, -2.66E-3,  2.3},
     {62.486253, 1517.00E-6,  0.078, 14.68E-3, 0.0, -4.77E-3,  0.9},
     {62.997974, 1503.00E-6,  1.660, 11.39E-3, 0.0, -3.34E-3,  2.2},
     {63.568515, 1087.00E-6,  2.110, 11.08E-3, 0.0, -4.17E-3,  2.0},
     {64.127764,  733.50E-6,  2.620, 10.78E-3, 0.0, -4.48E-3,  2.0},
     {64.678900,  463.50E-6,  3.190, 10.50E-3, 0.0, -5.10E-3,  1.8},
     {65.224067,  274.80E-6,  3.810, 10.20E-3, 0.0, -5.10E-3,  1.9},
     {65.764769,  153.00E-6,  4.480, 10.00E-3, 0.0, -5.70E-3,  1.8},
     {66.302088,   80.09E-6,  5.220,  9.70E-3, 0.0, -5.50E-3,  1.8},
     {66.836827,   39.46E-6,  6.000,  9.40E-3, 0.0, -5.90E-3,  1.7},
     {67.369595,   18.32E-6,  6.840,  9.20E-3, 0.0, -5.60E-3,  1.8},
     {67.900862,    8.01E-6,  7.740,  8.90E-3, 0.0, -5.80E-3,  1.7},
     {68.431001,    3.30E-6,  8.690,  8.70E-3, 0.0, -5.70E-3,  1.7},
     {68.960306,    1.28E-6,  9.690,  8.60E-3, 0.0, -5.60E-3,  1.7},
     {69.489021,    0.47E-6, 10.720,  8.50E-3, 0.0, -5.60E-3,  1.7},
     {70.017342,    0.16E-6, 11.830,  8.40E-3, 0.0, -5.60E-3,  1.7},
     {118.750341,  945.00E-6,  0.000, 15.92E-3, 0.0, -0.44E-3,  0.9},
     {368.498350,   67.90E-6,  0.020, 19.20E-3, 0.6,  0.00E00,  1.0},
     {424.763120,  638.00E-6,  0.011, 19.16E-3, 0.6,  0.00E00,  1.0},
     {487.249370,  235.00E-6,  0.011, 19.20E-3, 0.6,  0.00E00,  1.0},
     {715.393150,   99.60E-6,  0.089, 18.10E-3, 0.6,  0.00E00,  1.0},
     {773.838730,  671.00E-6,  0.079, 18.10E-3, 0.6,  0.00E00,  1.0},
     {834.145330,  180.00E-6,  0.079, 18.10E-3, 0.6,  0.00E00,  1.0}};
     
  // convert to the units of Liebe
  const double theta = 300./temperature;
  const double kPaPVap = 0.001*pVapour;
  const double kPaPDry = 0.001*pDry;
  const double fGHz = freq * 1e-9;
  
  // some coefficients
  const double ap = 1.4e-10*(1-1.2e-5*std::pow(fGHz,1.5));
  const double gamma0 = 5.6e-3*(kPaPDry + 1.1*kPaPVap)*std::pow(theta,0.8);
  // initial refractivity
  std::complex<double> result(2.588*kPaPDry*theta +
         3.07e-4*(1.0/(1.0+std::pow(fGHz/gamma0,2))-1)*kPaPDry*theta*theta,
         (2*3.07e-4/(gamma0*(1+std::pow(fGHz/gamma0,2))*(1+std::pow(fGHz/60,2))) + 
          ap*kPaPDry*std::pow(theta,2.5))*fGHz*kPaPDry*theta*theta);
  // sum the contributions of all the lines
  for (size_t l = 0; l < nLines; ++l) {
       const double S = lines[l][1]*kPaPDry*std::pow(theta,3)*exp(lines[l][2]*(1.-theta));
       const double gamma = lines[l][3]*(kPaPDry*std::pow(theta,0.8-lines[l][4]) + 1.1*kPaPVap*theta);
       const double delta = lines[l][5]*kPaPDry*std::pow(theta,lines[l][6]);
       const double x = (lines[l][0]-fGHz)*(lines[l][0]-fGHz) + gamma*gamma;
       const double y = (lines[l][0]+fGHz)*(lines[l][0]+fGHz) + gamma*gamma;
       const double z = (lines[l][0]+gamma*gamma/lines[l][0]);
       result += std::complex<double> (S*( (z-fGHz)/x + (z+fGHz)/y - 2./lines[l][0] + 
                                  delta*(1/x-1/y)*gamma*fGHz/lines[l][0]),
               S*( (1/x+1/y)*gamma*fGHz/lines[l][0] -
               delta*((lines[l][0]-fGHz)/x + (lines[l][0]+fGHz)/y)*fGHz/lines[l][0]));       
  }
  
  return result;
}
