//
// C++ Interface: STMath
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSTMATH_H
#define ASAPSTMATH_H

#include <map>
#include <string>
#include <iostream>

#include <casa/aips.h>
#include <casa/Arrays/Vector.h>
#include <casa/BasicSL/String.h>
#include <casa/Utilities/CountedPtr.h>

#include <scimath/Mathematics/InterpolateArray1D.h>

#include "Scantable.h"
#include "STDefs.h"
#include "STPol.h"

namespace asap {

/**
	* Mathmatical operations on Scantable objects
	* @author Malte Marquarding
*/
class STMath {
public:
	// typedef for long method name
  typedef casacore::InterpolateArray1D<casacore::Double,
                                   casacore::Float>::InterpolationMethod imethod;

  // typedef for std::map
  typedef std::map<std::string, imethod> imap;

/**
  * whether to operate on the given Scantable or return a new one
  * @param insitu the toggle for this behaviour
  */
  explicit STMath(bool insitu=true);

  virtual ~STMath();

  /**
   * get the currnt @attr inistu_ state
   */
  bool insitu() const { return insitu_;};

  /**
   * set the currnt @attr inistu state
   * @param b the new state
   */
  void setInsitu(bool b) { insitu_ = b; };


  /**
    * average a vector of Scantables
    * @param in the vector of Scantables to average
    * @param mask an optional mask to apply on specific weights
    * @param weight weighting scheme
    * @param avmode the mode ov averaging. Per "SCAN" or "ALL".
    * @return a casacore::CountedPtr<Scantable> which either holds a new Scantable
    * or returns the imput pointer.
    */
  casacore::CountedPtr<Scantable>
    average( const std::vector<casacore::CountedPtr<Scantable> >& in,
             const std::vector<bool>& mask = std::vector<bool>(),
             const std::string& weight = "NONE",
             const std::string& avmode = "SCAN");

  /**
    * median average a vector of Scantables. See also STMath::average
    * @param in the Scantable to average
    * @param mode the averaging mode. Currently only "MEDIAN"
    * @param avmode the mode ov averaging. Per "SCAN" or "ALL".
    * @return a casacore::CountedPtr<Scantable> which either holds a new Scantable
    * or returns the imput pointer.
    */
  casacore::CountedPtr<Scantable>
    averageChannel( const casacore::CountedPtr<Scantable> & in,
                    const std::string& mode = "MEDIAN",
                    const std::string& avmode = "SCAN");

  /**
    * Average polarisations together. really only useful if only linears are
    * available.
    * @param in the input Scantable
    * @param mask an optional mask if weight allows one
    * @param weight weighting scheme
    * @return
    */
  casacore::CountedPtr< Scantable >
    averagePolarisations( const casacore::CountedPtr< Scantable > & in,
                          const std::vector<bool>& mask,
                          const std::string& weight );

  /**
    * Average beams together.
    * @param in the input Scantable
    * @param mask an optional mask if weight allows one
    * @param weight weighting scheme
    * @return
    */
  casacore::CountedPtr< Scantable >
    averageBeams( const casacore::CountedPtr< Scantable > & in,
                   const std::vector<bool>& mask,
                   const std::string& weight );

  casacore::CountedPtr<Scantable>
    unaryOperate( const casacore::CountedPtr<Scantable>& in, float val,
                  const std::string& mode, bool tsys=false,
                  bool skip_flaggedrow=false );

  // array operation
  casacore::CountedPtr<Scantable>
    arrayOperate( const casacore::CountedPtr<Scantable>& in,
                  const std::vector<float> val,
                  const std::string& mode,
                  const std::string& opmode="channel",  
                  bool tsys=false,
                  bool skip_flaggedrow=false );

  // channel operation
  casacore::CountedPtr<Scantable>
    arrayOperateChannel( const casacore::CountedPtr<Scantable>& in,
                         const std::vector<float> val,
                         const std::string& mode, bool tsys=false,
                         bool skip_flaggedrow=false );

  // row operation
  casacore::CountedPtr<Scantable>
    arrayOperateRow( const casacore::CountedPtr<Scantable>& in,
                     const std::vector<float> val,
                     const std::string& mode, bool tsys=false,
                     bool skip_flaggedrow=false );

  // 2d array operation
  casacore::CountedPtr<Scantable>
    array2dOperate( const casacore::CountedPtr<Scantable>& in,
                  const std::vector< std::vector<float> > val,
                  const std::string& mode, bool tsys=false );

  casacore::CountedPtr<Scantable>
    binaryOperate( const casacore::CountedPtr<Scantable>& left,
		   const casacore::CountedPtr<Scantable>& right,
		   const std::string& mode);

  casacore::CountedPtr<Scantable> autoQuotient(const casacore::CountedPtr<Scantable>& in,
                                           const std::string& mode = "NEAREST",
                                           bool preserve = true);

  casacore::CountedPtr<Scantable> quotient( const casacore::CountedPtr<Scantable>& on,
                                        const casacore::CountedPtr<Scantable>& off,
                                        bool preserve = true );

  /**
    * Calibrate total power scans (translated from GBTIDL)
    * @param calon uncalibrated Scantable with CAL noise signal 
    * @param caloff uncalibrated Scantable with no CAL signal
    * @param tcal optional scalar Tcal, CAL temperature (K)
    * @return casacore::CountedPtr<Scantable> which holds a calibrated Scantable
    * (spectrum - average of the two CAL on and off spectra;
    * tsys - mean Tsys = <caloff>*Tcal/<calon-caloff> + Tcal/2)
    */  	    
  casacore::CountedPtr<Scantable> dototalpower( const casacore::CountedPtr<Scantable>& calon,
                                            const casacore::CountedPtr<Scantable>& caloff,
                                            casacore::Float tcal=1.0 );

  /**
    * Combine signal and reference scans (translated from GBTIDL)
    * @param sig Scantable which contains signal scans
    * @param ref Scantable which contains reference scans
    * @param smoothref optional Boxcar smooth width of the reference scans
    * default: no smoothing (=1)
    * @param tsysv optional scalar Tsys value at the zenith, required to 
    * set tau, as well 
    * @param tau optional scalar Tau value
    * @return casacore::CountedPtr<Scantable> which holds combined scans
    * (spectrum = (sig-ref)/ref * Tsys )
    */
  casacore::CountedPtr<Scantable> dosigref( const casacore::CountedPtr<Scantable>& sig,
                                        const casacore::CountedPtr<Scantable>& ref,
                                        int smoothref=1,
                                        casacore::Float tsysv=0.0,
                                        casacore::Float tau=0.0 );

  /**
    * Calibrate GBT Nod scan pairs (translated from GBTIDL)
    * @param s Scantable which contains Nod scans
    * @param scans Vector of scan numbers
    * @param smoothref optional Boxcar smooth width of the reference scans
    * @param tsysv optional scalar Tsys value at the zenith, required to
    * set tau, as well
    * @param tau optional scalar Tau value 
    * @param tcal optional scalar Tcal, CAL temperature (K)
    * @return casacore::CountedPtr<Scantable> which holds calibrated scans
    */
  casacore::CountedPtr<Scantable> donod( const casacore::CountedPtr<Scantable>& s,
                                     const std::vector<int>& scans,
                                     int smoothref=1,
                                     casacore::Float tsysv=0.0,
                                     casacore::Float tau=0.0,
                                     casacore::Float tcal=0.0 );

  /**
    * Calibrate frequency switched scans (translated from GBTIDL)
    * @param s Scantable which contains frequency switched  scans
    * @param scans Vector of scan numbers
    * @param smoothref optional Boxcar smooth width of the reference scans
    * @param tsysv optional scalar Tsys value at the zenith, required to
    * set tau, as well
    * @param tau optional scalar Tau value
    * @param tcal optional scalar Tcal, CAL temperature (K)
    * @return casacore::CountedPtr<Scantable> which holds calibrated scans
    */
  casacore::CountedPtr<Scantable> dofs( const casacore::CountedPtr<Scantable>& s,
                                    const std::vector<int>& scans,
                                    int smoothref=1,
                                    casacore::Float tsysv=0.0,
                                    casacore::Float tau=0.0,
                                    casacore::Float tcal=0.0 );

  /**
   * Calibrate data with Chopper-Wheel like calibration method 
   * which adopts position switching by antenna motion, 
   * wobbler (nutator) switching and On-The-Fly observation.
   * 
   * The method is applicable to APEX, and other telescopes other than GBT.
   *
   * @param a Scantable which contains ON and OFF scans
   * @param a string that indicates calibration mode 
   * @param a string that indicates antenna name
   **/
  casacore::CountedPtr<Scantable> cwcal( const casacore::CountedPtr<Scantable>& s,
                                       const casacore::String calmode,
                                       const casacore::String antname );

  /**
   * Calibrate frequency switched scans with Chopper-Wheel like 
   * calibration method.
   *
   * The method is applicable to APEX, and other telescopes other than GBT.
   * 
   * @param a Scantable which contains ON and OFF scans
   * @param a string that indicates antenna name
   **/
  casacore::CountedPtr<Scantable> cwcalfs( const casacore::CountedPtr<Scantable>& s,
                                       const casacore::String antname );


  /**
   * Folding frequency-switch data
   * @param sig
   * @param ref
   * @param choffset
   **/
  casacore::CountedPtr<Scantable> dofold( const casacore::CountedPtr<Scantable> &sig,
                                      const casacore::CountedPtr<Scantable> &ref,
                                      casacore::Double choffset,
                                      casacore::Double choffset2 = 0.0 );

  /**
   * ALMA calibration
   **/
  casacore::CountedPtr<Scantable> almacal( const casacore::CountedPtr<Scantable>& s,
                                       const casacore::String calmode ) ;
  casacore::CountedPtr<Scantable> almacalfs( const casacore::CountedPtr<Scantable>& s ) ;

  casacore::CountedPtr<Scantable>
    freqSwitch( const casacore::CountedPtr<Scantable>& in );

  std::vector<float> statistic(const casacore::CountedPtr<Scantable>& in,
                               const std::vector<bool>& mask,
                               const std::string& which);

  std::vector<float> statisticRow(const casacore::CountedPtr<Scantable>& in,
                               const std::vector<bool>& mask,
			       const std::string& which,
			       int row);

  std::vector< int > minMaxChan(const casacore::CountedPtr<Scantable>& in,
                                const std::vector<bool>& mask,
                                const std::string& which);

  casacore::CountedPtr<Scantable> bin( const casacore::CountedPtr<Scantable>& in,
                                   int width=5);
  casacore::CountedPtr<Scantable>
    resample(const casacore::CountedPtr<Scantable>& in,
             const std::string& method, float width);

  casacore::CountedPtr<Scantable>
    smooth(const casacore::CountedPtr<Scantable>& in, const std::string& kernel,
                      float width, int order=2);

  casacore::CountedPtr<Scantable>
    gainElevation(const casacore::CountedPtr<Scantable>& in,
                  const std::vector<float>& coeff,
                  const std::string& fileName,
                  const std::string& method);
  casacore::CountedPtr<Scantable>
    convertFlux(const casacore::CountedPtr<Scantable>& in, float d,
                float etaap, float jyperk);

  casacore::CountedPtr<Scantable> opacity(const casacore::CountedPtr<Scantable>& in,
                                      const std::vector<float>& tau);

  casacore::CountedPtr<Scantable>
    merge(const std::vector<casacore::CountedPtr<Scantable> >& in,
	  const std::string &freqTol = "");

  casacore::CountedPtr<Scantable>
    invertPhase( const casacore::CountedPtr<Scantable>& in);

  casacore::CountedPtr<Scantable>
    rotateXYPhase( const casacore::CountedPtr<Scantable>& in, float phase);

  casacore::CountedPtr<Scantable>
    rotateLinPolPhase( const casacore::CountedPtr<Scantable>& in, float phase);

  casacore::CountedPtr<Scantable>
    swapPolarisations(const casacore::CountedPtr<Scantable>& in);

  casacore::CountedPtr<Scantable>
    frequencyAlign( const casacore::CountedPtr<Scantable>& in,
                    const std::string& refTime = "",
                    const std::string& method = "cubic" );

  casacore::CountedPtr<Scantable>
    convertPolarisation( const casacore::CountedPtr<Scantable>& in,
                         const std::string& newtype);

  casacore::CountedPtr<Scantable>
    mxExtract( const casacore::CountedPtr<Scantable>& in,
               const std::string& srctype = "on");

  /**
   * "hard" flag the data, this flags everything selected in setSelection()
   * @param frequency the frequency to remove
   * @param width the number of lags to flag left to the side of the frequency
   */
  casacore::CountedPtr<Scantable>
    lagFlag( const casacore::CountedPtr<Scantable>& in, double start,
             double end, const std::string& mode="frequency");

  std::vector<float>
    fft( const casacore::CountedPtr<Scantable>& in,
	 const std::vector<int>& whichrow, 
	 bool getRealImag=false );

  // test for average spectra with different channel/resolution
  casacore::CountedPtr<Scantable>
    new_average( const std::vector<casacore::CountedPtr<Scantable> >& in,
		 const bool& compel, 
		 const std::vector<bool>& mask = std::vector<bool>(),
		 const std::string& weight = "NONE",
		 const std::string& avmode = "SCAN" )
    throw (casacore::AipsError) ;

private:
  casacore::CountedPtr<Scantable>  applyToPol( const casacore::CountedPtr<Scantable>& in,
                                           STPol::polOperation fptr,
                                           casacore::Float phase);

  static imethod stringToIMethod(const std::string& in);
  static WeightType stringToWeight(const std::string& in);

  void scaleByVector(casacore::Table& in,
                     const casacore::Vector<casacore::Float>& factor,
                     bool dotsys);

  void scaleFromAsciiTable(casacore::Table& in, const std::string& filename,
                           const std::string& method,
                           const casacore::Vector<casacore::Float>& xout,
                           bool dotsys);

  void scaleFromTable(casacore::Table& in, const casacore::Table& table,
                      const std::string& method,
                      const casacore::Vector<casacore::Float>& xout, bool dotsys);

  void convertBrightnessUnits(casacore::CountedPtr<Scantable>& in,
                              bool tokelvin, float cfac);

  casacore::CountedPtr< Scantable >
    smoothOther( const casacore::CountedPtr< Scantable >& in,
                 const std::string& kernel,
                 float width, int order=2 );

  casacore::CountedPtr< Scantable >
    getScantable(const casacore::CountedPtr< Scantable >& in, bool droprows);

  casacore::MaskedArray<casacore::Float>
    maskedArray( const casacore::Vector<casacore::Float>& s,
                 const casacore::Vector<casacore::uChar>& f );
  casacore::MaskedArray<casacore::Double>
    maskedArray( const casacore::Vector<casacore::Double>& s,
                 const casacore::Vector<casacore::uChar>& f );
  casacore::Vector<casacore::uChar>
    flagsFromMA(const casacore::MaskedArray<casacore::Float>& ma);

  // Frequency switching
  void calibrateFS( casacore::CountedPtr<Scantable> &sig,
                    casacore::CountedPtr<Scantable> &ref,
                    const casacore::CountedPtr<Scantable> &rsig,
                    const casacore::CountedPtr<Scantable> &rref,
                    const casacore::CountedPtr<Scantable> &sky,
                    const casacore::CountedPtr<Scantable> &hot,
                    const casacore::CountedPtr<Scantable> &cold,
                    const casacore::Vector<casacore::uInt> &rows ) ;
  void calibrateAPEXFS( casacore::CountedPtr<Scantable> &sig,
                        casacore::CountedPtr<Scantable> &ref,
                        const vector< casacore::CountedPtr<Scantable> > &on,
                        const vector< casacore::CountedPtr<Scantable> > &sky,
                        const vector< casacore::CountedPtr<Scantable> > &hot,
                        const vector< casacore::CountedPtr<Scantable> > &cold,
                        const casacore::Vector<casacore::uInt> &rows ) ;
  void copyRows( casacore::Table &out,
                 const casacore::Table &in,
                 casacore::uInt startout,
                 casacore::uInt startin,
                 casacore::uInt nrow,
                 casacore::Bool copySpectra=true,
                 casacore::Bool copyFlagtra=true,
                 casacore::Bool copyTsys=true ) ;
  casacore::CountedPtr<Scantable> averageWithinSession( casacore::CountedPtr<Scantable> &s,
                                                    vector<bool> &mask,
                                                    string weight ) ;

  bool insitu_;
};

}
#endif
