//
// C++ Implementation: STCalTsys
//
// Description:
//
//
// Author: Takeshi Nakazato <takeshi.nakazato@nao.ac.jp> (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <vector>

#include <casa/Arrays/ArrayMath.h>
#include <casa/Logging/LogIO.h>

#include "STSelector.h"
#include "STCalTsys.h"
#include "STDefs.h"
#include <atnf/PKSIO/SrcType.h>

using namespace std;
using namespace casa;

namespace asap {
  STCalTsys::STCalTsys(CountedPtr<Scantable> &s, vector<int> &iflist)
    : STCalibration(s, "TSYS"),
      iflist_(iflist)
{
  applytable_ = new STCalTsysTable(*s);
}

void STCalTsys::setupSelector(const STSelector &sel)
{
  sel_ = sel;
  vector<int> ifnos = sel_.getIFs();
  if (ifnos.size() > 0) {
    int nif = 0;
    int nifOrg = iflist_.size();
    vector<int> iflistNew(iflist_);
    for (int i = 0; i < nifOrg; i++) {
      if (find(ifnos.begin(), ifnos.end(), iflist_[i]) != ifnos.end()) {
        iflistNew[nif] = iflist_[i];
        ++nif;
      }
    }
    if (nif == 0) {
      LogIO os(LogOrigin("STCalTsys", "setupSelector", WHERE));
      os << LogIO::SEVERE << "Selection contains no data." << LogIO::EXCEPTION;
    }
    sel_.setIFs(iflistNew);
  }
  else {
    sel_.setIFs(iflist_);
  }
}

void STCalTsys::appenddata(uInt scanno, uInt cycleno, 
			   uInt beamno, uInt ifno, uInt polno, 
			   uInt freqid, Double time, Float elevation, 
			   Vector<Float> any_data)
{
  STCalTsysTable *p = dynamic_cast<STCalTsysTable *>(&(*applytable_));
  p->appenddata(scanno, cycleno, beamno, ifno, polno,
		freqid, time, elevation, any_data);
}

}
