//
// C++ Implementation: STSelector
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <tables/Tables/ExprNode.h>
#include <tables/Tables/TableParse.h>
#include <casa/BasicSL/String.h>
#include <casa/iostream.h>
#include <casa/iomanip.h>

#include "STSelector.h"

using namespace asap;
using namespace casa;

STSelector::STSelector() :
  taql_("")
{
}

STSelector::STSelector( const STSelector & other ) :
  selections_(other.selections_),
  taql_(other.taql_) {
}

STSelector& STSelector::operator=( const STSelector& other )
{
  if (&other != this) {
    this->selections_ = other.selections_;
    this->taql_ = other.taql_;
  }
  return *this;
}

STSelector::~STSelector()
{
}

void STSelector::setScans( const std::vector< int >& scans )
{
  set("SCANNO", scans);
}

void STSelector::setBeams( const std::vector< int >& beams )
{
  set("BEAMNO", beams);
}

void STSelector::setIFs( const std::vector< int >& ifs )
{
  set("IFNO", ifs);
}

void STSelector::setPolarizations( const std::vector< int >& pols )
{
  set("POLNO", pols);
}

void STSelector::set(const std::string& key, const std::vector< int >& val)
{
  if ( val.size() > 0 ) {
    selections_[key] = val;
  }
}

void STSelector::setTaQL( const std::string& taql )
{
  taql_ = taql;
}

Table STSelector::apply( const Table& tab )
{
  if ( empty() ) {
    return tab;
  }
  TableExprNode query;
  idmap::const_iterator it = selections_.begin();
  for (it; it != selections_.end(); ++it) {
    TableExprNode theset(Vector<Int>( (*it).second ));
    if ( query.isNull() ) {
      query = tab.col((*it).first).in(theset);
    } else {
      query = tab.col((*it).first).in(theset) && query;
    }
  }
  // add taql query
  if ( taql_.size() > 0 ) {
    Table tmpt = tab;

    if ( !query.isNull() ) { // taql and selection
      tmpt = tableCommand(taql_, tab(query));
    } else { // taql only
      tmpt = tableCommand(taql_, tab);
    }
    return tmpt.copyToMemoryTable("dummy");
  } else {
    return tab(query).copyToMemoryTable("dummy");
  }
}

std::vector< int > STSelector::get( const std::string& key)
{
  if (selections_.count(key) > 0) {
    return  std::vector<int>();//selections_[key];
  }
}

std::vector< int > STSelector::getScans( )
{
  return get("SCANNO");
}

std::vector< int > STSelector::getBeams( )
{
  return get("BEAMNO");
}

std::vector< int > STSelector::getIFs( )
{
  return get("IFNO");
}

std::vector< int > STSelector::getPols( )
{
  return get("POLNO");
}

std::string asap::STSelector::print( )
{
  ostringstream oss;
  oss.flags(std::ios_base::left);
  oss << setw(15) << "Selection:";
  if ( empty() ) {
    oss << "none";
    return String(oss);
  }

  idmap::const_iterator it = selections_.begin();
  while (it != selections_.end()) {
    if ( it != selections_.begin() )
      oss << setw(15) << " ";
    oss << it->first << ": " << Vector<Int>(it->second);
    ++it;
    if ( it != selections_.end() ) oss << endl;
  }
  if ( taql_.size() > 0 ) {
    oss << endl << setw(15) << "" << taql_;
  }
  return String(oss);
}

bool asap::STSelector::empty( ) const
{
  return (selections_.empty() && taql_.size() == 0 );
}
