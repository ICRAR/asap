//
// C++ Interface: FillerBase
//
// Description:
//
// This class is the Base class for all data fillers.
// The derived filler needs to implement
// open()
// close()
// fill()
//
// The fill() method usually iterates over the source data and calls
// the setXYZ() methods for. After all the data for a row has been set via
// these methods, the fill() method needs to call commitRow() to write the
// data to the scantable.
// All arguments which are defaulted in the setXYZ() methods are optional. All
// others should be set explicitly.
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ASAPFILLERBASE_H
#define ASAPFILLERBASE_H

// STL
#include <string>
#include <vector>
// AIPS++
#include <casa/aips.h>
#include <casa/Utilities/CountedPtr.h>
#include <casa/Arrays/Vector.h>
#include <tables/Tables/TableRow.h>
#include "Scantable.h"

namespace asap
{

class FillerBase
{
  public:
    explicit FillerBase(casacore::CountedPtr<Scantable> stable);
    virtual ~FillerBase() {;}

    virtual bool open(const std::string& filename, const casacore::Record& rec) = 0;
    //    virtual bool open(const std::string& filename) = 0;
    virtual void fill() = 0;
    virtual void close() = 0;

    void setReferenceRegex(const std::string& rx) { referenceRx_ = rx; }
    std::string getReferenceRegex() { return referenceRx_;  }

  protected:

    void commitRow();
    void setHeader(const STHeader& header);
    void setSpectrum(const casacore::Vector<casacore::Float>& spectrum,
                             const casacore::Vector<casacore::uChar>& flags,
                             const casacore::Vector<casacore::Float>& tsys);
    void setFlagrow(casacore::uInt flag);
    void setOpacity(casacore::Float opacity=0.0f);
    void setIndex(casacore::uInt scanno, casacore::uInt cycleno,
                          casacore::uInt ifno, casacore::uInt polno,
                          casacore::uInt beamno=0);
    void setFrequency(casacore::Double refpix, casacore::Double refval,
                              casacore::Double incr);
    void setMolecule(const casacore::Vector<casacore::Double>& restfreq);
    void setDirection(const casacore::Vector<casacore::Double>& dir,
                              casacore::Float az=0.0f, casacore::Float el=0.0f);

    void setFocus(casacore::Float pa=0.0f, casacore::Float faxis=0.0f,
                          casacore::Float ftan=0.0f, casacore::Float frot=0.0f);
    void setTime(casacore::Double mjd, casacore::Double integration);
    void setWeather(casacore::Float temperature=0.0f,
                            casacore::Float pressure=0.0f,
                            casacore::Float humidity=0.0f,
                            casacore::Float windspeed=0.0f,
                            casacore::Float windaz=0.0f);
    void setWeather2(casacore::Float temperature=0.0f,
                            casacore::Float pressure=0.0f,
                            casacore::Float humidity=0.0f,
                            casacore::Float windspeed=0.0f,
                            casacore::Float windaz=0.0f);
    void setTcal(const casacore::String& caltime="",
                         const casacore::Vector<casacore::Float>& tcal=casacore::Vector<casacore::Float>());
    void setTcal2(const casacore::String& caltime="",
                         const casacore::Vector<casacore::Float>& tcal=casacore::Vector<casacore::Float>());
    void setScanRate(const casacore::Vector<casacore::Double>& srate=casacore::Vector<casacore::Double>());
    void setReferenceBeam(casacore::Int beamno=-1);
    void setSource(const std::string& name, casacore::Int type,
                           const std::string& fieldname="",
                           const casacore::Vector<casacore::Double>& dir=casacore::Vector<casacore::Double>(),
                           const casacore::Vector<casacore::Double>& propermot=casacore::Vector<casacore::Double>(),
                           casacore::Double velocity=0.0);

    casacore::CountedPtr< Scantable > table_;

  private:

    FillerBase();
    FillerBase(const FillerBase&);
    FillerBase& operator=(const FillerBase&);
    casacore::String referenceRx_;
    casacore::TableRow row_;
 
    std::vector< casacore::Vector<casacore::Double> > mEntry_ ;
    std::vector<casacore::uInt> mIdx_ ;
    std::vector< casacore::Vector<casacore::Double> > fEntry_ ;
    std::vector<casacore::uInt> fIdx_ ;
    std::vector< casacore::Vector<casacore::Float> > wEntry_ ;
    std::vector<casacore::uInt> wIdx_ ;
};


};
#endif
