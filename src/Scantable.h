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
   * @return const casa::Table reference
   */
  const casa::Table& table() const;

  /**
   * get a reference to the underlying casa::Table with the Selection
   * object applied if set
   * @return casa::Table reference
   */
  casa::Table& table();

  void setSelection(const STSelector& selection);
  void unsetSelection();
  /**
   * set the header
   * @param[in] sdh an SDHeader object
   */
  void putSDHeader( const SDHeader& sdh );

  /**
   * get the header information
   * @return an SDHeader object
   */
  SDHeader getSDHeader( ) const;


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
  int nScan() const;

  //casa::MDirection::Types getDirectionReference() const;
  //casa::MEpoch::Types getTimeReference() const;

  /**
   * Get global antenna position
   * @return
   */
  casa::MPosition getAntennaPosition() const;

  /**
   *
   * @return
   */

  std::string getFluxUnit() const;

  /**
   *
   * @param unit
   */
  void setFluxUnit( const std::string& unit );

  /**
   *
   * @param instrument
   */
  void setInstrument( const std::string& instrument );

  void calculateAZEL();

  /**
   * "hard" flags
   * @param[in] whichrow
   */
  void flag();

  int nbeam(int scanno=-1) const;
  int nif(int scanno=-1) const;
  int npol(int scanno=-1) const;
  int nchan(int scanno=-1, int ifno=-1) const;

  int nrow(int scanno=-1) const;

  double getInterval(int whichrow=0) const;

  float getTsys(int whichrow=0) const;

  std::vector<bool> getMask(int whichrow=0) const;
  std::vector<float> getSpectrum(int whichrow=0) const;

  std::vector<float> getStokesSpectrum( int whichrow=0,
                                        bool dopol=false) const;
  std::string getPolarizationLabel(bool linear, bool stokes,
                                   bool linpol,
                                   int polidx=-1) const;

  void makePersistent(const std::string& filename);


  void select(const STSelector& sel);

  const STSelector& selection() const { return selector_; }

  std::vector<std::string> getHistory() const;
  void addHistory(const std::string& hist);

  casa::Table getHistoryTable() const;
  void appendToHistoryTable(const casa::Table& otherHist);

  std::string summary(bool verbose=false);
  std::string getTime(int whichrow=-1, bool showdate=true) const;


  STFrequencies& frequencies() { return freqTable_; }
  STWeather& weather() { return weatherTable_; }
  STFocus& focus() { return focusTable_; }
  STTcal& tcal() { return tcalTable_; }
  STMolecules& molecules() { return moleculeTable_; }

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

  void setupHistoryTable();
  void setupMoleculeTable();
  void setupFitTable();

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
  casa::Table fitTable_;
  casa::Table historyTable_;

  // Cached Columns to avoid reconstructing them for each row get/put
  casa::ScalarColumn<casa::Double> timeCol_, integrCol_;
  casa::MDirection::ScalarColumn dirCol_;
  casa::ScalarColumn<casa::Double> azCol_;
  casa::ScalarColumn<casa::Double> elCol_;
  casa::ScalarColumn<casa::Float> paraCol_;
  casa::ScalarColumn<casa::String> srcnCol_, fldnCol_;
  casa::ScalarColumn<casa::uInt> scanCol_, beamCol_, cycleCol_;
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
