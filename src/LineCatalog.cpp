
//
// C++ Implementation: LineCatalog
//
// Description: A class representing a line catalog
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//

// std includes

// casa includes
#include <casa/Exceptions/Error.h>
#include <casa/OS/Path.h>
#include <casa/OS/File.h>
#include <tables/Tables/ReadAsciiTable.h>
#include <tables/Tables/TableParse.h>

#include "LineCatalog.h"

using namespace casa;

namespace asap
{

LineCatalog::LineCatalog(const std::string& name)
{
  Path path(name);
  std::string inname = path.expandedName();
  File f(inname);
  if (f.isDirectory()) { //assume its a table
    table_ = Table(inname);
  } else {
    String formatString;
    // formatSring , TableType, ascii file name, TableDesc name, table name, autoheader
  table_ = readAsciiTable(formatString, Table::Plain, inname, "", "", True);
  }
  baseTable_ = table_;
}

void LineCatalog::setStrengthLimits(float smin, float smax)
{
  table_ = setLimits(smin, smax, "Column4");
}

void LineCatalog::setFrequencyLimits(float fmin, float fmax)
{
  table_ = setLimits(fmin, fmax, "Column2");
}

void LineCatalog::setPattern(const std::string& name, const std::string& stype)
{
  std::string mode = stype+"('"+name+"')";
  std::string taql = "SELECT from $1 WHERE Column1 == "+mode;
  Table tmp = tableCommand(taql, table_);
  if (tmp.nrow() > 0) table_ = tmp;
  else throw(AipsError("No match."));
}

Table LineCatalog::setLimits(float lmin, float lmax, const std::string& colname)
{
  Table tmp = table_(table_.col(colname) > lmin && table_.col(colname) > lmax);
  if (tmp.nrow() > 0) return tmp;
  else throw(AipsError("No match."));
}

void LineCatalog::save(const std::string& name)
{
  Path path(name);
  std::string inname = path.expandedName();
  table_.deepCopy(inname, Table::New);
}

} // namespace
