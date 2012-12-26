// C++ Interface: STSideBandSep
//
// Description:
//    A class to invoke sideband separation of Scantable
//
// Author: Kana Sugimoto <kana.sugi@nao.ac.jp>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//

// STL
#include <ctype.h>

// cascore
#include <casa/OS/File.h>
#include <casa/Logging/LogIO.h>
#include <casa/Quanta/QuantumHolder.h>

#include <measures/Measures/MFrequency.h>
#include <measures/Measures/MCFrequency.h>

#include <tables/Tables/TableRecord.h>
#include <tables/Tables/TableVector.h>

// asap
#include "STSideBandSep.h"

using namespace std ;
using namespace casa ;
using namespace asap ;

//#ifndef KS_DEBUG
//#define KS_DEBUG
//#endif

namespace asap {

// constructors
STSideBandSep::STSideBandSep()
  : sigIfno_(0), ftol_(-1), 
    loTime_(-1), lo1Freq_(-1), loDir_("")
{
#ifdef KS_DEBUG
  cout << "Default constructor STSideBandSep()" << endl;
#endif
  // Default LO frame is TOPO
  loFrame_ = MFrequency::TOPO;
};

STSideBandSep::~STSideBandSep()
{
#ifdef KS_DEBUG
  cout << "Destructor ~STSideBandSep()" << endl;
#endif
};


void STSideBandSep::setFrequency(const unsigned int ifno, const double freqtol, string frame)
{
#ifdef KS_DEBUG
  cout << "STSideBandSep::setFrequency" << endl;
  cout << "IFNO = " << ifno << endl;
  cout << "freq tol = " << freqtol << " [Hz]" << endl;
  cout << "frame = " << frame << endl;
#endif
};

// Temporal function to set scantable of image sideband
void STSideBandSep::setImageTable(const ScantableWrapper &s)
{
#ifdef KS_DEBUG
  cout << "STSideBandSep::setImageTable" << endl;
  cout << "got image table nrow = " << s.nrow() << endl;
#endif
  imgTab_p = s.getCP();
  AlwaysAssert(!imgTab_p.null(),AipsError);

};

//
void STSideBandSep::setLO1(const double lo1, string frame, double reftime, string refdir)
{
  lo1Freq_ = lo1;
  MFrequency::getType(loFrame_, frame);
  loTime_ = reftime;
  loDir_ = refdir;

#ifdef KS_DEBUG
  cout << "STSideBandSep::setLO1" << endl;
  if (lo1Freq_ > 0.)
    cout << "lo1 = " << lo1Freq_ << " [Hz] (" << frame << ")" << endl;
  if (loTime_ > 0.)
    cout << "ref time = " << loTime_ << " [day]" << endl;
  if (!loDir_.empty())
    cout << "ref direction = " << loDir_ << " [day]" << endl;
#endif
};

void STSideBandSep::setLO1Asdm(string asdmname)
{
  // Check for existance of the asdm.
  if (!checkFile(asdmname)) {
     throw(AipsError("File does not exist"));
  }
  if (asdmname[(asdmname.size()-1)] == '/')
    asdmname = asdmname.substr(0,(asdmname.size()-2));

  asdmName_ = asdmname;

#ifdef KS_DEBUG
  cout << "STSideBandSep::setLO1Asdm" << endl;
  cout << "asdm name = " << asdmName_ << endl;
#endif
};

void STSideBandSep::solveImageFreqency()
{
#ifdef KS_DEBUG
  cout << "STSideBandSep::solveImageFrequency" << endl;
#endif
  LogIO os(LogOrigin("STSideBandSep","solveImageFreqency()", WHERE));
  os << "Start calculating frequencies of image side band" << LogIO::POST;

  if (imgTab_p.null())
    throw AipsError("STSideBandSep::solveImageFreqency - an image side band scantable should be set first");

  // Check for the availability of LO1
  os << "Looking for LO1. lo1Freq_ = " << lo1Freq_ << LogIO::POST; 
  if (lo1Freq_ > 0.) {
    os << "Using user defined LO1 frequency." << LogIO::POST;
  } else if (!asdmName_.empty() > 0) {
      // ASDM name is set.
    os << "Using user defined ASDM: " << asdmName_ << LogIO::POST;
    if (!getLo1FromAsdm(asdmName_)) {
      throw AipsError("Failed to get LO1 frequency from ASDM");
    }
  } else {
    // Try getting ASDM name from scantable header
    os << "Try getting information from scantable header" << LogIO::POST;
    if (!getLo1FromTab(imgTab_p)) {
      throw AipsError("Failed to get LO1 frequency from asis table");
    }
  }
  // LO1 should now be ready.
  if (lo1Freq_ < 0.)
    throw(AipsError("Got negative LO1 Frequency"));

  // The code assumes that imgTab_p has only an IF and only a FREQ_ID
  // is associated to an IFNO
  // TODO: More complete Procedure would be
  // 1. Get freq IDs associated to sigIfno_
  // 2. Get freq information of the freq IDs
  // 3. For each freqIDs, get freq infromation in TOPO and an LO1 
  //    frequency and calculate image band frequencies.
  STFrequencies freqTab_ = imgTab_p->frequencies();
  // get the base frame of table
  const MFrequency::Types tabframe = freqTab_.getFrame(true);
  TableVector<uInt> freqIdVec( imgTab_p->table(), "FREQ_ID" );
  // assuming single freqID per IFNO
  uInt freqid = freqIdVec(0);
  double refpix, refval, increment ;
  freqTab_.getEntry(refpix, refval, increment, freqid);
  //MFrequency sigrefval = MFrequency(MVFrequency(refval),tabframe);
  // get freq infromation of sigIfno_ in loFrame_
  const MPosition mp = imgTab_p->getAntennaPosition();
  MEpoch me;
  MDirection md;
  if (loTime_ < 0.)
    me = imgTab_p->getEpoch(-1);
  else
    me = MEpoch(MVEpoch(loTime_));
  if (loDir_.empty()) {
    ArrayColumn<Double> srcdirCol_;
    srcdirCol_.attach(imgTab_p->table(), "SRCDIRECTION");
    // Assuming J2000 and SRCDIRECTION in unit of rad
    Quantum<Double> srcra = Quantum<Double>(srcdirCol_(0)(IPosition(1,0)), "rad");
    Quantum<Double> srcdec = Quantum<Double>(srcdirCol_(0)(IPosition(1,1)), "rad");
    md = MDirection(srcra, srcdec, MDirection::J2000);
    //imgTab_p->getDirection(0);
  } else {
    // parse direction string
    string::size_type pos0 = loDir_.find(" ");
    
    if (pos0 == string::npos) {
      throw AipsError("bad string format in LO1 direction");
    }
    string::size_type pos1 = loDir_.find(" ", pos0+1);
    String sepoch, sra, sdec;
    if (pos1 != string::npos) {
      sepoch = loDir_.substr(0, pos0);
      sra = loDir_.substr(pos0+1, pos1-pos0);
      sdec = loDir_.substr(pos1+1);
    }
    MDirection::Types epoch;
    MDirection::getType(epoch, sepoch);
    QuantumHolder qh ;
    String err ;
    qh.fromString( err, sra);
    Quantum<Double> ra = qh.asQuantumDouble() ;
    qh.fromString( err, sdec ) ;
    Quantum<Double> dec = qh.asQuantumDouble() ;
    md = MDirection(ra.getValue("rad"), dec.getValue("rad"),epoch);
  }
  // Summary (LO1)
  Vector<Double> dirvec = md.getAngle(Unit(String("rad"))).getValue();
  os << "[LO1 settings]" << LogIO::POST;
  os << "- Frequency: " << lo1Freq_ << " [Hz] ("
     << MFrequency::showType(loFrame_) << ")" << LogIO::POST;
  os << "- Reference time: " << me.get(Unit(String("d"))).getValue()
     << " [day]" << LogIO::POST;
  os << "- Reference direction: [" << dirvec[0] << ", " << dirvec[1]
     << "] (" << md.getRefString() << ") " << LogIO::POST;
  //
  MeasFrame mframe( me, mp, md );
  MFrequency::Convert toloframe(loFrame_, MFrequency::Ref(tabframe, mframe));
  MFrequency::Convert tobframe(tabframe, MFrequency::Ref(loFrame_, mframe));
  // convert refval to loFrame_
  double sigrefval;
  if (tabframe == loFrame_)
    sigrefval = refval;
  else
    sigrefval = toloframe(Quantum<Double>(refval, "Hz")).get("Hz").getValue();

  //Summary (signal)
  os << "[Signal side band]" << LogIO::POST;
  os << "- IFNO: " << imgTab_p->getIF(0) << " (FREQ_ID = " << freqid << ")"
     << LogIO::POST;
  os << "- Reference value: " << refval << " [Hz] ("
     << MFrequency::showType(tabframe) << ") = "
     << sigrefval << " [Hz] (" <<  MFrequency::showType(loFrame_)
     << ")" << LogIO::POST;
  os << "- Reference pixel: " << refpix  << LogIO::POST;
  os << "- Increment: " << increment << " [Hz]" << LogIO::POST;
  // calculate image band incr and refval in loFrame_
  Double imgincr = -increment;
  Double imgrefval = 2 * lo1Freq_ - sigrefval;
  Double imgrefval_tab = imgrefval;
  // convert imgrefval back to table base frame
  if (tabframe != loFrame_)
    imgrefval = tobframe(Quantum<Double>(imgrefval, "Hz")).get("Hz").getValue();
  // set new frequencies table
  uInt fIDnew = freqTab_.addEntry(refpix, imgrefval, imgincr);
  // update freq_ID in table.
  freqIdVec = fIDnew;
  // Summary (Image side band)
  os << "[Image side band]" << LogIO::POST;
  os << "- IFNO: " << imgTab_p->getIF(0) << " (FREQ_ID = " << freqIdVec(0)
     << ")" << LogIO::POST;
  os << "- Reference value: " << imgrefval << " [Hz] ("
     << MFrequency::showType(tabframe) << ") = "
     << imgrefval_tab << " [Hz] (" <<  MFrequency::showType(loFrame_)
     << ")" << LogIO::POST;
  os << "- Reference pixel: " << refpix  << LogIO::POST;
  os << "- Increment: " << imgincr << " [Hz]" << LogIO::POST;
};

Bool STSideBandSep::checkFile(const string name, string type)
{
  File file(name);
  if (!file.exists()){
    return false;
  } else if (type.empty()) {
    return true;
  } else {
    // Check for file type
    switch (tolower(type[0])) {
    case 'f':
      return file.isRegular(True);
    case 'd':
      return file.isDirectory(True);
    case 's':
      return file.isSymLink();
    default:
      throw AipsError("Invalid file type. Available types are 'file', 'directory', and 'symlink'.");
    }
  }
};

bool STSideBandSep::getLo1FromAsdm(const string asdmname)
{
  // Check for relevant tables.
  string spwname = asdmname + "/SpectralWindow.xml";
  string recname = asdmname + "/Receiver.xml";
  if (!checkFile(spwname) || !checkFile(recname)) {
    throw(AipsError("Could not find subtables in ASDM"));
  }

  return false;

};


bool STSideBandSep::getLo1FromTab(CountedPtr< Scantable > &scantab)
{
  Table& tab = scantab->table();
  // Check for relevant tables.
  const String spwname = tab.keywordSet().asString("ASDM_SPECTRALWINDOW");
  const String recname = tab.keywordSet().asString("ASDM_RECEIVER");
  if (!checkFile(spwname,"directory") || !checkFile(recname,"directory")) {
    throw(AipsError("Could not find subtables in MS"));
  }
  // get

  return false;

};


} //namespace asap
