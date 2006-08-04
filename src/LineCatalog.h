//
// C++ Interface: LineCatalog
//
// Description:
//
//
// Author: Malte Marquarding <Malte.Marquarding@csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef LINECATALOG_H
#define LINECATALOG_H

#include <string>

#include <casa/aips.h>
#include <tables/Tables/Table.h>

namespace asap {
/**
  * A represenation of a line catalog, which can be ASCII, or an aips++ table
  *
  * ASCII catalogs have to be formatted like JPL.
  * Name frequency error log(I)
  * Only  "name", "frequency" and log(I) are used at this stage
  *
  * @author Malte Marquarding
  * @date $Date:$
  */
class LineCatalog {
public:
  /**
    *
    * @param name the name of the ASCII file or aips++ table
    */
  LineCatalog(const std::string& name = "jpl");
  /**
    * select a subset of the table by frequency range
    * @param fmin the lower frequency bound
    * @param fmin the upper frequency bound
    */

  virtual ~LineCatalog() {}

  /**
   * select a subset of the data by frequency range
   * @param fmin the lower frequency bound
   * @param fmax the upper frequency bound
   */
  void setFrequencyLimits(float fmin, float fmax);
  /**
    * select a subset of the table by line strength range
    * @param smin the lower strength bound
    * @param smin the upper strength bound
    */
  void setStrengthLimits(float smin, float smax);
  /**
    * select a subset of the data by name pattern match (unix-style)
    * @param name the string pattern e.g. "*CS*"
    * @param ptype pattern type e.g.
    * @li "pattern"
    * @li "regex"
    */
  void setPattern(const std::string& name, const std::string& ptype="pattern");
  /**
    * save the table  with current limits to disk (as an aips++ table)
    * @param name the filename
    */
  void save(const std::string& name);
  /**
    * Return a string representation of this table
    * @param an integer descriing the row number to show
    * default -1 is all rows
    * @return std::string
    */
  std::string summary(int row=-1) const;

  double getFrequency(uint row) const;

  std::string getName(uint row) const;

private:
  /**
   * utility function to hadle range limits
   * @param lmin the lower limit
   * @param lmax the upper limit
   * @param colname the columd to apply the limits to
   * @return a new casa::Table
   */
  casa::Table setLimits(float lmin, float lmax, const std::string& colname);

  // the table with seelection
  casa::Table table_;
  // the pristine table
  casa::Table baseTable_;
};

} // namespace

#endif //LINECATALOG_H
