#include <vector>

#include <casa/Arrays/Vector.h>
#include <casa/BasicSL/String.h>
#include <casa/Utilities/CountedPtr.h>

#include "Scantable.h"
#include "STTcal.h"
#include "STIdxIter.h"
#include "STSelector.h"

using namespace casa;
using namespace asap;

namespace {
// Iteration Helper

// template<class T>
// class IterationHelper
// {
// public:
//   static void Iterate(T &calibrator)
//   {
//     vector<string> cols( 3 ) ;
//     cols[0] = "BEAMNO" ;
//     cols[1] = "POLNO" ;
//     cols[2] = "IFNO" ;
//     STIdxIter2 calibrator.GetIterator(cols);//iter( out, cols ) ;
//     STSelector sel ;
//     while ( !iter.pastEnd() ) {
//       Record ids = iter.currentValue() ;
//       stringstream ss ;
//       ss << "SELECT FROM $1 WHERE "
//          << "BEAMNO==" << ids.asuInt(cols[0]) << "&&"
//          << "POLNO==" << ids.asuInt(cols[1]) << "&&"
//          << "IFNO==" << ids.asuInt(cols[2]) ;
//       //cout << "TaQL string: " << ss.str() << endl ;
//       sel.setTaQL( ss.str() ) ;
//       Vector<uInt> rows = iter.getRows( SHARE ) ;
//       calibrator.Calibrate(sel, rows);
//       aoff->setSelection( sel ) ;
//       // out should be an exact copy of s except that SPECTRA column is empty
//       calibrateALMA( out, s, aoff, rows ) ;
//       aoff->unsetSelection() ;
//       sel.reset() ;
//       iter.next() ;
//     }    
//   }
// };

// Interpolation Helper
class TcalData
{
public:
  TcalData(CountedPtr<Scantable> s)
    : table_(s)
  {}
  ~TcalData() {}
  const String method_name() const {return "getTcalFromTime";}
  uInt nrow() const {return table_->nrow();}
  Vector<Float> GetEntry(int idx) const
  {
    String time;
    uInt tcalid = table_->getTcalId(idx);
    Vector<Float> return_value;
    table_->tcal().getEntry(time, return_value, tcalid);
    return return_value;
  }  
private:
  CountedPtr<Scantable> table_;
};
  
class TsysData
{
public:
  TsysData(CountedPtr<Scantable> s)
    : tsyscolumn_(s->table(), "TSYS")
  {}
  ~TsysData() {}
  const String method_name() const {return "getTsysFromTime";}
  uInt nrow() const {return tsyscolumn_.nrow();}
  Vector<Float> GetEntry(int idx) const
  {
    return tsyscolumn_(idx);
  }  
private:
  ROArrayColumn<Float> tsyscolumn_;
};

class SpectralData
{
public:
  SpectralData(Matrix<Float> s)
    : data_(s)
  {}
  ~SpectralData() {}
  const String method_name() const {return "getSpectraFromTime";}
  uInt nrow() const {return data_.ncolumn();}
  Vector<Float> GetEntry(int idx) const
  {
    return data_.column(idx);
  }  
private:
  Matrix<Float> data_;
};


vector<int> getRowIdFromTime(double reftime, const Vector<Double> &t)
{
  //   double reft = reftime ;
  double dtmin = 1.0e100 ;
  double dtmax = -1.0e100 ;
  //   vector<double> dt ;
  int just_before = -1 ;
  int just_after = -1 ;
  Vector<Double> dt = t - reftime ;
  for ( unsigned int i = 0 ; i < dt.size() ; i++ ) {
    if ( dt[i] > 0.0 ) {
      // after reftime
      if ( dt[i] < dtmin ) {
	just_after = i ;
	dtmin = dt[i] ;
      }
    }
    else if ( dt[i] < 0.0 ) {
      // before reftime
      if ( dt[i] > dtmax ) {
	just_before = i ;
	dtmax = dt[i] ;
      }
    }
    else {
      // just a reftime
      just_before = i ;
      just_after = i ;
      dtmax = 0 ;
      dtmin = 0 ;
      break ;
    }
  }
  
  vector<int> v(2) ;
  v[0] = just_before ;
  v[1] = just_after ;
  
  return v ;
}
  
template<class T>
class SimpleInterpolationHelper
{
 public:
  static Vector<Float> GetFromTime(double reftime,
				   const Vector<Double> &timeVec,
				   const vector<int> &idx,
				   const T &data,
				   const string mode)
  {
    Vector<Float> return_value;
    LogIO os_;
    LogIO os( LogOrigin( "STMath", data.method_name(), WHERE ) ) ;
    if ( data.nrow() == 0 ) {
      os << LogIO::SEVERE << "No row in the input scantable. Return empty tcal." << LogIO::POST ;
    }
    else if ( data.nrow() == 1 ) {
      return_value = data.GetEntry(0);
    }
    else {
      if ( mode == "before" ) {
	int id = -1 ;
	if ( idx[0] != -1 ) {
	  id = idx[0] ;
	}
	else if ( idx[1] != -1 ) {
	  os << LogIO::WARN << "Failed to find a scan before reftime. return a spectrum just after the reftime." << LogIO::POST ;
	  id = idx[1] ;
	}
	
	return_value = data.GetEntry(id);
      }
      else if ( mode == "after" ) {
	int id = -1 ;
	if ( idx[1] != -1 ) {
	  id = idx[1] ;
	}
	else if ( idx[0] != -1 ) {
	  os << LogIO::WARN << "Failed to find a scan after reftime. return a spectrum just before the reftime." << LogIO::POST ;
	  id = idx[1] ;
	}
	
	return_value = data.GetEntry(id);
      }
      else if ( mode == "nearest" ) {
	int id = -1 ;
	if ( idx[0] == -1 ) {
	  id = idx[1] ;
	}
	else if ( idx[1] == -1 ) {
	  id = idx[0] ;
	}
	else if ( idx[0] == idx[1] ) {
	  id = idx[0] ;
	}
	else {
	  double t0 = timeVec[idx[0]] ;
	  double t1 = timeVec[idx[1]] ;
	  if ( abs( t0 - reftime ) > abs( t1 - reftime ) ) {
	    id = idx[1] ;
	  }
	  else {
	    id = idx[0] ;
	  }
	}
	return_value = data.GetEntry(id);
      }
      else if ( mode == "linear" ) {
	if ( idx[0] == -1 ) {
	  // use after
	  os << LogIO::WARN << "Failed to interpolate. return a spectrum just after the reftime." << LogIO::POST ;
	  int id = idx[1] ;
	  return_value = data.GetEntry(id);
	}
	else if ( idx[1] == -1 ) {
	  // use before
	  os << LogIO::WARN << "Failed to interpolate. return a spectrum just before the reftime." << LogIO::POST ;
	  int id = idx[0] ;
	  return_value = data.GetEntry(id);
	}
	else if ( idx[0] == idx[1] ) {
	  // use before
	  //os << "No need to interporate." << LogIO::POST ;
	  int id = idx[0] ;
	  return_value = data.GetEntry(id);
	}
	else {
	  // do interpolation

	  double t0 = timeVec[idx[0]] ;
	  double t1 = timeVec[idx[1]] ;
	  Vector<Float> value0 = data.GetEntry(idx[0]);
	  Vector<Float> value1 = data.GetEntry(idx[1]);
	  double tfactor = (reftime - t0) / (t1 - t0) ;
	  for ( unsigned int i = 0 ; i < value0.size() ; i++ ) {
	    value1[i] = ( value1[i] - value0[i] ) * tfactor + value0[i] ;
	  }
	  return_value = value1;
	}
      }
      else {
	os << LogIO::SEVERE << "Unknown mode" << LogIO::POST ;
      }
    }
    return return_value ;
  }
};
  
  
} // anonymous namespace
