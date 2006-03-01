//
// C++ Interface: Scantable
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSCANTABLE_H
#define ASAPSCANTABLE_H

// STL
#include <string>
#include <vector>
// AIPS++
#include <casa/aips.h>
#include <casa/Arrays/MaskedArray.h>
#include <casa/BasicSL/String.h>
#include <casa/Utilities/CountedPtr.h>

#include <coordinates/Coordinates/SpectralCoordinate.h>

#include <tables/Tables/Table.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <measures/TableMeasures/ScalarMeasColumn.h>

#include "SDLog.h"
#include "SDContainer.h"
#include "STFrequencies.h"
#include "STWeather.h"
#include "STFocus.h"
#include "STTcal.h"
#include "STMolecules.h"
#include "STSelector.h"
#include "STHistory.h"



namespace asap {

/**
  * This class contains and wraps a casa::Table, which is used to store
  * all the information. This can be either a MemoryTable or a
  * disk based Table.
  * It provides access functions to the underlying table
  * It contains n subtables:
  * @li weather
  * @li frequencies
  * @li molecules
  * @li tcal
  * @li focus
  * @li fits
  *
  * @brief The main ASAP data container
  * @author Malte Marquarding
  * @date
  * @version
*/
class Scantable : private SDLog
{
public:
  /**
   * Default constructor
   */
  Scantable(casa::Table::TableType ttype = casa::Table::Memory);

  /**
   * Create a Scantable object form an existing table on disk
   * @param[in] name the name of the existing Scantable
   */
  Scantable(const std::string& name, casa::Table::TableType ttype = casa::Table::Memory);

  /// @fixme this is only sensible for MemoryTables....
  Scantable(const Scantable& other, bool clear=true);

  /**
   * Destructor
   */
  virtual ~Scantable();

  /**
   * get a const reference to the underlying casa::Table
   * @return consantcasa::Table reference
   */
  const casa::Table& table() const;

  /**
   * get a reference to the underlying casa::Table with the Selection
   * object applied if set
   * @return casa::Table reference
   */
  casa::Table& table();


  /**
   * Get a handle to the selection object
   * @return constant STSelector reference
   */
  const STSelector& getSelection() const { return selector_; }

  /**
   * Set the data to be a subset as defined by the STSelector
   * @param selection a STSelector object
   */
  void setSelection(const STSelector& selection);

  /**
   * unset the selection of the data
   */
  void unsetSelection();
  /**
   * set the header
   * @param[in] sdh an SDHeader object
   */
  void setHeader( const SDHeader& sdh );

  /**
   * get the header information
   * @return an SDHeader object
   */
  SDHeader getHeader( ) const;
  /**
   * Checks if the "other" Scantable is conformant with this,
   * i.e. if  header values are the same.
   * @param[in] other another Scantable
   * @return true or false
   */
  bool conformant( const Scantable& other);

  /**
   * return the number of scans in the table
   * @return number of scans in the table
   */
  int nscan() const;

  //casa::MDirection::Types getDirectionReference() const;
  //casa::MEpoch::Types getTimeReference() const;

  /**
   * Get global antenna position
   * @return casa::MPosition
   */
  casa::MPosition getAntennaPosition() const;

  /**
   *  Return the Flux unit of the data, e.g. "Jy" or "K"
   * @return string
   */
  std::string getFluxUnit() const;

  /**
   * Set the Flux unit of the data
   * @param unit a string representing the unit, e.g "Jy" or "K"
   */
  void setFluxUnit( const std::string& unit );

  /**
   *
   * @param instrument a string representing an insturment. see xxx
   */
  void setInstrument( const std::string& instrument );

  /**
   * (Re)calculate the azimuth and elevationnfor all rows
   */
  void calculateAZEL();

  /**
   * "hard" flag the data, this flags everything selected in setSelection()
   */
  void flag();

  /**
   * Return a list of row numbers with respect to the original table.
   * @return a lsi of rownumbers with respect to the original table
   */
  std::vector<unsigned int> rownumbers() const;


  /**
   * Get the number of beams in the data or a specific scan
   * @param scanno the scan number to get the number of beams for.
   * If scanno<0 the number is retrieved from the header.
   * @return an integer number
   */
  int nbeam(int scanno=-1) const;
  /**
   * Get the number of IFs in the data or a specific scan
   * @param scanno the scan number to get the number of IFs for.
   * If scanno<0 the number is retrieved from the header.
   * @return an integer number
   */
  int nif(int scanno=-1) const;
  /**
   * Get the number of polarizations in the data or a specific scan
   * @param scanno the scan number to get the number of polarizations for.
   * If scanno<0 the number is retrieved from the header.
   * @return an integer number
   */
  int npol(int scanno=-1) const;

  /**
   * Get the number of channels in the data or a specific IF. This currently
   * varies only with IF number
   * @param ifno the IF number to get the number of channels for.
   * If ifno<0 the number is retrieved from the header.
   * @return an integer number
   */
  int nchan(int ifno=-1) const;

  /**
   * Get the number of integartion cycles
   * @param scanno the scan number to get the number of rows for.
   * If scanno<0 the number is retrieved from the header.
   * @return
   */
  int nrow(int scanno=-1) const;

  int ncycle(int scanno=-1) const;

  int getBeam(int whichrow) const;
  int getIF(int whichrow) const;
  int getPol(int whichrow) const;

  double getInterval(int whichrow) const
    { return integrCol_(whichrow); }

  float getTsys(int whichrow) const;
  float getElevation(int whichrow) const
    { return elCol_(whichrow); }
  float getAzimuth(int whichrow) const
    { return azCol_(whichrow); }
  float getParangle(int whichrow) const
    { return paraCol_(whichrow); }

  std::vector<bool> getMask(int whichrow) const;
  std::vector<float> getSpectrum(int whichrow) const;

  std::vector<float> getStokesSpectrum( int whichrow=0,
                                        bool dopol=false) const;
  std::string getPolarizationLabel(bool linear, bool stokes,
                                   bool linpol,
                                   int polidx=-1) const;

  /**
   * Write the Scantable to disk
   * @param filename the output file name
   */
  void makePersistent(const std::string& filename);

  std::vector<std::string> getHistory() const
    { return historyTable_.getHistory(); };

  void addHistory(const std::string& hist) { historyTable_.addEntry(hist); }

  void appendToHistoryTable(const STHistory& otherhist)
    { historyTable_.append(otherhist); }

  std::string summary(bool verbose=false);
  std::string getTime(int whichrow=-1, bool showdate=true) const;

  // returns unit, conversion frame, doppler, base-frame

  /**
   * Get the frequency set up
   * This is forwarded to the STFrequencies subtable
   * @return unit, frame, doppler
   */
  std::vector<std::string> getCoordInfo() const
    { return freqTable_.getInfo(); };

  void setCoordInfo(std::vector<string> theinfo)
    { return freqTable_.setInfo(theinfo); };

  std::string getAbcissaLabel(int whichrow) const;
  std::vector<double> getRestFrequencies() const
    { return moleculeTable_.getRestFrequencies(); }

  void setRestFrequencies(double rf, const std::string& = "Hz");
  void setRestFrequencies(const std::string& name);

  STFrequencies& frequencies() { return freqTable_; }
  STWeather& weather() { return weatherTable_; }
  STFocus& focus() { return focusTable_; }
  STTcal& tcal() { return tcalTable_; }
  STMolecules& molecules() { return moleculeTable_; }
  STHistory& history() { return historyTable_; }

private:
  /**
   * Turns a time vale into a formatted string
   * @param x
   * @return
   */
  std::string formatSec(casa::Double x) const;

  std::string formatTime(const casa::MEpoch& me, bool showdate)const;

  /**
   *  Turns a casa::MDirection into a nicely formatted string
   * @param md an casa::MDirection
   * @return
   */
  std::string formatDirection(const casa::MDirection& md) const;


  /**
   * Create a unique file name for the paged (temporary) table
   * @return just the name
   */
  static casa::String generateName();

  /**
   * attach to cached columns
   */
  void attach();

  /**
   * Set up the main casa::Table
   */
  void setupMainTable();

  void setupFitTable();

  void attachSubtables();
  /**
   * Convert an "old" asap1 style row index into a new index
   * @param[in] therow
   * @return and index into @table_
   */
  int rowToScanIndex(int therow);

  static const unsigned int version_ = 2;

  STSelector selector_;

  casa::Table::TableType type_;

  // the actual data
  casa::Table table_;
  casa::Table originalTable_;

  STTcal tcalTable_;
  STFrequencies freqTable_;
  STWeather weatherTable_;
  STFocus focusTable_;
  STMolecules moleculeTable_;
  STHistory historyTable_;

  casa::Table fitTable_;

  // Cached Columns to avoid reconstructing them for each row get/put
  casa::ScalarColumn<casa::Double> integrCol_;
  casa::MDirection::ScalarColumn dirCol_;
  casa::MEpoch::ScalarColumn timeCol_;
  casa::ScalarColumn<casa::Double> azCol_;
  casa::ScalarColumn<casa::Double> elCol_;
  casa::ScalarColumn<casa::Float> paraCol_;
  casa::ScalarColumn<casa::String> srcnCol_, fldnCol_;
  casa::ScalarColumn<casa::uInt> scanCol_, beamCol_, ifCol_, polCol_, cycleCol_;
  casa::ScalarColumn<casa::Int> rbeamCol_;
  casa::ArrayColumn<casa::Float> specCol_, tsCol_;
  casa::ArrayColumn<casa::uChar> flagsCol_;

  // id in frequencies table
  casa::ScalarColumn<casa::uInt> mfreqidCol_;
  // id in tcal table
  casa::ScalarColumn<casa::uInt> mtcalidCol_;

  casa::ArrayColumn<casa::String> histitemCol_;
  casa::ScalarColumn<casa::uInt> mfitidCol_, fitidCol_;
  // id in weather table and main table
  casa::ScalarColumn<casa::uInt> mweatheridCol_;

  casa::ScalarColumn<casa::uInt> mfocusidCol_;

  casa::ScalarColumn<casa::uInt> mmolidCol_;

};


} // namespace

#endif
