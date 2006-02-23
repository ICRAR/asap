//
// C++ Interface: STSelector
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSTSELECTOR_H
#define ASAPSTSELECTOR_H

#include <string>
#include <vector>
#include <map>

#include <tables/Tables/Table.h>

namespace asap {

/**
A class to set a subselection of a Scantable

@author Malte Marquarding
*/
class STSelector {

public:
  STSelector();
  STSelector(const STSelector& other);

  STSelector& operator=(const STSelector& other);

  virtual ~STSelector();

  void setScans(const std::vector<int>& scans);
  void setBeams(const std::vector<int>& beams);
  void setIFs(const std::vector<int>& ifs);
  void setPolarizations(const std::vector<int>& pols);
  void setCycles(const std::vector<int>& cycs);
  void setTaQL(const std::string& taql);

  std::vector<int> getScans();
  std::vector<int> getBeams();
  std::vector<int> getIFs();
  std::vector<int> getPols();
  std::vector<int> getCycles();

  casa::Table apply(const casa::Table& tab);
  casa::Table operator()(const casa::Table& tab) { return apply(tab); };

  void reset() { selections_.clear(); taql_ = "";};

  bool empty() const;

  std::string print();

private:
  //
  std::vector< int > get( const std::string& key);
  void set(const std::string& key, const std::vector< int >& val);

  typedef std::map<std::string, std::vector<int> > idmap;
  // has to be mutable, as to stl limitations
  mutable idmap selections_;
  std::string taql_;
};

}

#endif
