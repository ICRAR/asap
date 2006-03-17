//
// C++ Interface: Scantable
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPSCANTABLEWRAPPER_H
#define ASAPSCANTABLEWRAPPER_H

#include <vector>
#include <string>
#include <casa/Arrays/Vector.h>

#include "MathUtils.h"
#include "STFit.h"
#include "Scantable.h"

namespace asap {
/**
  * This class contains and wraps a CountedPtr<Scantable>, as the CountedPtr
  * class does not provide the dor operator which is need for references
  * in boost::python
  * see Scantable for interfce description
  *
  * @brief The main ASAP data container wrapper
  * @author Malte Marquarding
  * @date 2006-02-23
  * @version 2.0a
*/
class ScantableWrapper {

public:
  ScantableWrapper( const std::string& name,
                    const std::string& type="")
  {
    casa::Table::TableType tp = casa::Table::Memory;
    if ( type == "disk") tp = casa::Table::Plain;
    table_ = new Scantable(name, tp);
  }

  ScantableWrapper(const std::string& type="")
  {
    casa::Table::TableType tp = casa::Table::Memory;
    if ( type == "disk") tp = casa::Table::Plain;
    table_= new Scantable(tp);
  }

  ScantableWrapper(casa::CountedPtr<Scantable> cp) : table_(cp) {;}

  ScantableWrapper(const ScantableWrapper& mt) :
    table_(mt.getCP()) {;}

  void assign(const ScantableWrapper& mt)
    { table_= mt.getCP(); }

  ScantableWrapper copy() {
    return ScantableWrapper(new Scantable(*(this->getCP()), false));
  }

  std::vector<float> getSpectrum( int whichrow=0,
                                  const std::string& poltype="" ) const {
    return table_->getSpectrum(whichrow, poltype);
  }
  //  std::string getPolarizationLabel(bool linear, bool stokes, bool linPol, int polIdx) const {
  // Boost fails with 4 arguments.
  std::string getPolarizationLabel(int index, const std::string& ptype) const {
    return table_->getPolarizationLabel(index, ptype);
  }

  std::vector<double> getAbcissa(int whichrow=0) const
    { return table_->getAbcissa(whichrow); }

  std::string getAbcissaLabel(int whichrow=0) const
    { return table_->getAbcissaLabel(whichrow); }

  float getTsys(int whichrow=0) const
    { return table_->getTsys(whichrow); }

  std::string getTime(int whichrow=0) const
    { return table_->getTime(whichrow); }

  std::string getFluxUnit() const { return table_->getFluxUnit(); }

  void setFluxUnit(const std::string& unit) { table_->setFluxUnit(unit); }

  void setInstrument(const std::string& name) {table_->setInstrument(name);}

  std::vector<bool> getMask(int whichrow=0) const
    { return table_->getMask(whichrow); }

  void flag() { table_->flag(); }

  std::string getSourceName(int whichrow=0) const
    { return table_->getSourceName(whichrow); }

  float getElevation(int whichrow=0) const
    { return table_->getElevation(whichrow); }

  float getAzimuth(int whichrow=0) const
    { return table_->getAzimuth(whichrow); }

  float getParAngle(int whichrow=0) const
    { return table_->getParAngle(whichrow); }


  void setSpectrum(std::vector<float> spectrum, int whichrow=0)
    { table_->setSpectrum(spectrum, whichrow); }

  int getIF(int whichrow) const {return table_->getIF(whichrow);}
  int getBeam(int whichrow) const {return table_->getBeam(whichrow);}
  int getPol(int whichrow) const {return table_->getPol(whichrow);}
  int getCycle(int whichrow) const {return table_->getCycle(whichrow);}
  int getScan(int whichrow) const {return table_->getScan(whichrow);}

  STSelector getSelection() const { return table_->getSelection(); }
  void setSelection(const STSelector& sts)
    { return table_->setSelection(sts);}

  std::string getPolType() const { return table_->getPolType(); }

  int nif(int scanno=-1) const {return table_->nif(scanno);}
  int nbeam(int scanno=-1) const {return table_->nbeam(scanno);}
  int npol(int scanno=-1) const {return table_->npol(scanno);}
  int nchan(int ifno=-1) const {return table_->nchan(ifno);}
  int nscan() const {return table_->nscan();}
  int nrow() const {return table_->nrow();}
  int ncycle(int scanno) const {return table_->ncycle(scanno);}
  ///@todo int nstokes() {return table_->nStokes();}

  void makePersistent(const std::string& fname)
    { table_->makePersistent(fname); }

  void setRestFrequencies(double rf, const std::string& unit)
    { table_->setRestFrequencies(rf, unit); }

  void setRestFrequencies(const std::string& name) {
    table_->setRestFrequencies(name);
  }

  std::vector<double> getRestFrequencies() const
    { return table_->getRestFrequencies(); }

  void setCoordInfo(std::vector<string> theinfo) {
    table_->setCoordInfo(theinfo);
  }
  std::vector<string> getCoordInfo() const {
    return table_->getCoordInfo();
  }

  casa::CountedPtr<Scantable> getCP() const {return table_;}
  Scantable* getPtr() {return &(*table_);}

  std::string summary(bool verbose=false) const {
    return table_->summary(verbose);
  }

  std::vector<std::string> getHistory()const
    { return table_->getHistory(); }

  void addHistory(const std::string& hist)
    { table_->addHistory(hist); }
  /*
  void addFit(int whichrow, const std::vector<double>& p,
              const std::vector<bool>& m, const std::vector<string>& f,
              const std::vector<int>& c) {

    casa::Vector<casa::Double> p2(p);
    casa::Vector<casa::Bool> m2(m);
    casa::Vector<casa::String> f2 = mathutil::toVectorString(f);
    casa::Vector<casa::Int> c2(c);
    table_->addFit(casa::uInt(whichrow), p2,m2,f2,c2);
  }
  SDFitTable getSDFitTable(int whichrow) {
    return table_->getSDFitTable(casa::uInt(whichrow));
  }
  */
  void calculateAZEL() { table_->calculateAZEL(); };

private:
  casa::CountedPtr<Scantable> table_;
};

} // namespace
#endif

