//
// C++ Interface: STFrequencies
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSTFREQUENCIES_H
#define ASAPSTFREQUENCIES_H

#include <casa/aips.h>
#include <casa/BasicSL/String.h>
#include <coordinates/Coordinates/SpectralCoordinate.h>
#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>

#include "STSubTable.h"

namespace asap {

/**
The Frequencies subtable of the Scantable

@author Malte Marquarding
@brief The frequency subtable of the Scantable
*/
class STFrequencies : public STSubTable {
public:
  STFrequencies() {;}
  explicit STFrequencies(casacore::Table tab);
  explicit STFrequencies(const Scantable& parent);

  virtual ~STFrequencies();

  STFrequencies& operator=(const STFrequencies& other);

  /**
   * Add a new Entry to the Frequency subtable. This checks for duplicates.
   * @param[in] refpix the reference pixel
   * @param[in] refval the reference value
   * @param[in] inc the increment
   * @return an index into the frequency table
   */
  casacore::uInt addEntry( casacore::Double refpix, casacore::Double refval,
                       casacore::Double inc );

  /**
   * Retrieve the frequency values for a specific id via references
   * @param refpix the reference pixel
   * @param refval the reference value
   * @param inc the increment
   * @param id the identifier
   */
  void getEntry( casacore::Double& refpix, casacore::Double& refval,
                 casacore::Double& inc, casacore::uInt id );

  /***
   * Set the frequency values for a specific id via references
   * @param refpix the reference pixel
   * @param refval the reference value
   * @param inc    the increment
   * @param id     the identifier
   *
   * 17/09/2008 Takeshi Nakazato
   ***/
  void setEntry( casacore::Double refpix, casacore::Double refval,
		 casacore::Double inc, casacore::uInt id ) ;


  bool conformant(const STFrequencies& other) const;

  /**
   * Retrieve  the frequency values as a casacore::SpectralCoordinate
   * @param freqID
   * @return casacore::SpectralCoordinate
   */
  casacore::SpectralCoordinate getSpectralCoordinate( casacore::uInt freqID ) const;

  /**
  casacore::SpectralCoordinate getSpectralCoordinate( const casacore::MDirection& md,
                                                  const casacore::MPosition& mp,
                                                  const casacore::MEpoch& me,
                                                  casacore::Double restfreq,
                                                  casacore::uInt freqID
                                                  ) const;
  **/
  casacore::SpectralCoordinate getSpectralCoordinate( const casacore::MDirection& md,
                                                  const casacore::MPosition& mp,
                                                  const casacore::MEpoch& me,
                                                  casacore::Vector<casacore::Double> restfreq,
                                                  casacore::uInt freqID
                                                  ) const;

  /**
   * Return the unit of the frequency values
   * @return casacore::Unit
   */
  casacore::Unit getUnit() const;
  std::string getUnitString() const;

  /**
   * Return the doppler type of the values
   * @return casacore::MDoppler::Types
   */
  casacore::MDoppler::Types getDoppler() const;
  std::string getDopplerString() const;


  /**
   * Return the frame type, e.g MFrequency::TOPO
   * @param base return the base frame or the user frame
   * @return casacore::MFrequency::Types
   */
  casacore::MFrequency::Types getFrame(bool base=false) const;

  /**
   * Return a string representation of the frame type, e.g TOPO
   * @param base return the base frame or the user frame
   * @return the string representation of the frame
   */
  std::string getFrameString(bool base=false) const;

  /**
   * set the frequency frame from a string value
   * @param frame a string identifier
   */
  void setFrame(const std::string& frame, bool base=false);
  /**
   * set the frequency frame from a casacore::MFrequency::Types
   * @param frame casacore::MFrequency::Types
   */
  void setFrame(casacore::MFrequency::Types frame, bool base=false);
  void setUnit( const std::string & unit );
  void setDoppler( const std::string & doppler );
  /**
   * rescale the whole table by a given factor
   * @param factor the factor to bin or resample by
   * @param mode the rescaling mode
   * @li "BIN"
   * @li "RESAMPLE"
   */
  void rescale(casacore::Float factor, const std::string& mode);

  /**
   * get the reference frequency at a given channel for a specidif identifier
   * @param id the identifier
   * @param channel the channel number
   * @return teh reference frequency
   */
  float getRefFreq(casacore::uInt id, casacore::uInt channel);

  /**
    * shift the reference pixel by an integer amount
    * @param npix the shift in pixels
    * @param id the coordinate id
    */
  void shiftRefPix(int npix, casacore::uInt id);
  /**
   * Return this table or s specific row as a string representation
   * @param id the identifier. If id<0 all rows are returned
   * @return a string
   */
  std::string print(int id=-1, casacore::Bool strip=casacore::False) const;

  std::vector<std::string> getInfo() const;
  void setInfo( const std::vector<std::string>& theinfo );

  const casacore::String& name() const { return name_; }

  /**
   * Examine given set of refpix, refval, and increment matches
   * any of the rows within a tolerance of freqTolInHz. If match,
   * return true and id is filled properly. Otherwise, return false
   * and id may have invalid value.
   *
   * @param[in] refpix
   * @param[in] refval
   * @param[in] inc
   * @param[in] freqTolInHz
   * @param[out] id
   * @return boolean indicating match with any rows or not
   */
  bool match( casacore::Double refpix, casacore::Double refval, casacore::Double inc,
	      casacore::Double freqTolInHz, casacore::uInt &id);

private:

  /**
   * setup the the column structure of the casacore::table
   */
  void setup();
  /**
   * the actual binning of the SpectralCoordinate as called by rescale
   * @param sc
   * @param factor the bin factor
   * @return casacore::SpectralCoordinate
   */
  casacore::SpectralCoordinate binCsys(const casacore::SpectralCoordinate& sc, casacore::Int factor);
  /**
   * the actual resampling of the SpectralCoordinate as called by rescale
   * @param sc
   * @param width the resacle width. Can be decimal.
   * @return
   */
  casacore::SpectralCoordinate resampleCsys(const casacore::SpectralCoordinate& sc, casacore::Float width);

  static const casacore::String name_;
  casacore::ScalarColumn<casacore::Double> refvalCol_, refpixCol_, incrCol_;
};

}

#endif
