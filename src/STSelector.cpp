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
#include <tables/Tables/ExprNode.h>
#include <casa/BasicSL/String.h>
#include <casa/iostream.h>
#include <casa/iomanip.h>
#include <casa/Exceptions/Error.h>

#include "MathUtils.h"
#include "STPol.h"
#include "STSelector.h"

using namespace asap;
using namespace casa;

STSelector::STSelector() :
  taql_("")
{
}

STSelector::STSelector( const STSelector&  other ) :
  intselections_(other.intselections_),
  stringselections_(other.stringselections_),
  taql_(other.taql_) {
}

STSelector& STSelector::operator=( const STSelector& other )
{
  if (&other != this) {
    this->intselections_ = other.intselections_;
    this->stringselections_ = other.stringselections_;
    this->taql_ = other.taql_;
  }
  return *this;
}

STSelector::~STSelector()
{
}

void STSelector::setScans( const std::vector< int >& scans )
{
  setint("SCANNO", scans);
}

void STSelector::setBeams( const std::vector< int >& beams )
{
  setint("BEAMNO", beams);
}

void STSelector::setIFs( const std::vector< int >& ifs )
{
  setint("IFNO", ifs);
}

void STSelector::setPolarizations( const std::vector< int >& pols )
{
  setint("POLNO", pols);
}

void asap::STSelector::setCycles( const std::vector< int >& cycs )
{
  setint("CYCLENO", cycs);
}

void asap::STSelector::setName( const std::string& sname )
{
  std::string sql = "SELECT from $1 WHERE SRCNAME == pattern('"+sname+"')";
  setTaQL(sql);
}

void STSelector::setint(const std::string& key, const std::vector< int >& val)
{
  if ( val.size() > 0 ) {
    intselections_[key] = val;
  }
}

void STSelector::setstring( const std::string& key,
                            const std::vector<std::string>& val )
{
  if ( val.size() > 0 ) {
    stringselections_[key] = val;
  }
}

void STSelector::setTaQL( const std::string& taql )
{
  taql_ = taql;
}


void asap::STSelector::setSortOrder( const std::vector< std::string > & order )
{
  order_.resize(order.size(), True);
  for (int i=0;i<order.size();++i) {
    order_[i] = order[i];
  }
}

Table STSelector::apply( const Table& tab )
{
  if ( empty() ) {
    return sort(tab);
  }
  TableExprNode query;
  intidmap::const_iterator it = intselections_.begin();
  for (it; it != intselections_.end(); ++it) {
    TableExprNode theset(Vector<Int>( (*it).second ));
    if ( query.isNull() ) {
      query = tab.col((*it).first).in(theset);
    } else {
      query = tab.col((*it).first).in(theset) && query;
    }
  }
  stringidmap::const_iterator it1 = stringselections_.begin();
  for (it1; it1 != stringselections_.end(); ++it1) {
    TableExprNode theset(mathutil::toVectorString( (*it1).second ));
    if ( query.isNull() ) {
      query = tab.col((*it1).first).in(theset);
    } else {
      query = tab.col((*it1).first).in(theset) && query;
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
    return sort(tmpt);
  } else {
    if ( query.isNull() ) {
      return sort(tab);
    } else {
      return sort(tab(query));
    }
  }
}

std::vector< int > STSelector::getint( const std::string& key )
{
  if (intselections_.count(key) > 0) {
    return  intselections_[key];
  }
  return  std::vector<int>();
}

std::vector< int > STSelector::getScans( )
{
  return getint("SCANNO");
}

std::vector< int > STSelector::getBeams( )
{
  return getint("BEAMNO");
}

std::vector< int > STSelector::getIFs( )
{
  return getint("IFNO");
}

std::vector< int > STSelector::getPols( )
{
  return getint("POLNO");
}

std::vector< int > asap::STSelector::getCycles( )
{
  return getint("CYCLENO");
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

  intidmap::const_iterator it = intselections_.begin();
  while (it != intselections_.end()) {
    if ( it != intselections_.begin() )
      oss << setw(15) << " ";
    oss << it->first << ": " << Vector<Int>(it->second);
    ++it;
    if ( it != intselections_.end() ) oss << endl;
  }
  stringidmap::const_iterator it1 = stringselections_.begin();
  while (it1 != stringselections_.end()) {
    if ( it1 != stringselections_.begin() )
      oss << setw(15) << " ";
    oss << it1->first << ": " << mathutil::toVectorString(it1->second);
    ++it1;
    if ( it1 != stringselections_.end() ) oss << endl;
  }
  if ( taql_.size() > 0 ) {
    oss << endl << setw(15) << "" << taql_;
  }
  return String(oss);
}

bool asap::STSelector::empty( ) const
{
  return (intselections_.empty() && taql_.size() == 0 );
}

casa::Table asap::STSelector::sort( const casa::Table & tab )
{
  if (order_.nelements() > 0) {
    cout << "sorting" << endl;
    Table t = tab.sort(order_);
    return t;
  } else {
    return tab;
  }
}


void asap::STSelector::setPolFromStrings( const std::vector< std::string >& pols )
{
  poltypes_.clear();
  std::vector<int> polints;
  std::vector<std::string>::const_iterator strit = pols.begin();
  for (strit; strit != pols.end(); ++strit) {
    std::pair<int, std::string> val;
    try {
       val = STPol::polFromString(*strit);
    } catch (AipsError& e) {
      poltypes_.clear();
      throw(e);
    }
    polints.push_back(val.first);
    poltypes_.push_back(val.second);
  }
  setint("POLNO", polints);
}

std::vector< std::string > asap::STSelector::getPolTypes( )
{
  return poltypes_;
}
