//
// C++ Interface: FillerBase
//
// Description:
//
//
// Author: Malte Marquarding <asap@atnf.csiro.au>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <casa/Containers/RecordField.h>
#include <tables/Tables/ExprNode.h>

#include "FillerBase.h"

using namespace casa;

namespace asap {

FillerBase::FillerBase(casa::CountedPtr<Scantable> stable) :
  table_(stable)
{
    row_ = TableRow(table_->table());

    // FIT_ID is -1 by default
    RecordFieldPtr<Int> fitIdCol( row_.record(), "FIT_ID" ) ;
    *fitIdCol = -1 ;
}

void FillerBase::setHeader(const STHeader& header)
{
  table_->setHeader(header);
}

void FillerBase::setSpectrum(const Vector<Float>& spectrum,
                             const Vector<uChar>& flags,
                             const Vector<Float>& tsys)
{
  RecordFieldPtr< Array<Float> > specCol(row_.record(), "SPECTRA");
  RecordFieldPtr< Array<uChar> > flagCol(row_.record(), "FLAGTRA");
  RecordFieldPtr< Array<Float> > tsysCol(row_.record(), "TSYS");

  //*specCol = spectrum;
  //*flagCol = flags;
  //*tsysCol = tsys;
  specCol.define(spectrum);
  flagCol.define(flags);
  tsysCol.define(tsys);
}

void FillerBase::setFlagrow(uInt flag)
{
  RecordFieldPtr<uInt> flagrowCol(row_.record(), "FLAGROW");
  *flagrowCol = flag;
}

void FillerBase::setOpacity(Float opacity)
{
  RecordFieldPtr<Float> tauCol(row_.record(), "OPACITY") ;
  *tauCol = opacity ;
}

void FillerBase::setIndex(uInt scanno, uInt cycleno, uInt ifno, uInt polno,
                          uInt beamno)
{
  RecordFieldPtr<uInt> beamCol(row_.record(), "BEAMNO");
  RecordFieldPtr<uInt> ifCol(row_.record(), "IFNO");
  RecordFieldPtr<uInt> polCol(row_.record(), "POLNO");
  RecordFieldPtr<uInt> cycleCol(row_.record(), "CYCLENO");
  RecordFieldPtr<uInt> scanCol(row_.record(), "SCANNO");
  *beamCol = beamno;
  *cycleCol = cycleno;
  *ifCol = ifno;
  *polCol = polno;
  *scanCol = scanno;
}

void FillerBase::setFrequency(Double refpix, Double refval,
                              Double incr)
{
  /// @todo this has to change when nchan isn't global anymore
  uInt id= table_->frequencies().addEntry(refpix, refval, incr);
  RecordFieldPtr<uInt> mfreqidCol(row_.record(), "FREQ_ID");
  *mfreqidCol = id;

}


void FillerBase::setMolecule(const Vector<Double>& restfreq)
{
  Vector<String> tmp;
  uInt id = table_->molecules().addEntry(restfreq, tmp, tmp);
  RecordFieldPtr<uInt> molidCol(row_.record(), "MOLECULE_ID");
  *molidCol = id;
}

void FillerBase::setDirection(const Vector<Double>& dir,
                              Float az, Float el)
{
  RecordFieldPtr<Array<Double> > dirCol(row_.record(), "DIRECTION");
  *dirCol = dir;
  RecordFieldPtr<Float> azCol(row_.record(), "AZIMUTH");
  *azCol = az;
  RecordFieldPtr<Float> elCol(row_.record(), "ELEVATION");
  *elCol = el;
}

void FillerBase::setFocus(Float pa, Float faxis,
                      Float ftan, Float frot)
{
  RecordFieldPtr<uInt> mfocusidCol(row_.record(), "FOCUS_ID");
  uInt id = table_->focus().addEntry(pa, faxis, ftan, frot);
  *mfocusidCol = id;
}

void FillerBase::setTime(Double mjd, Double interval)
{
    RecordFieldPtr<Double> mjdCol(row_.record(), "TIME");
    *mjdCol = mjd;
    RecordFieldPtr<Double> intCol(row_.record(), "INTERVAL");
    *intCol = interval;

}

void FillerBase::setWeather(Float temperature, Float pressure,
                        Float humidity,
                        Float windspeed, Float windaz)
{
    uInt id = table_->weather().addEntry(temperature, pressure,
                                         humidity, windspeed, windaz);
    RecordFieldPtr<uInt> mweatheridCol(row_.record(), "WEATHER_ID");
    *mweatheridCol = id;
}

void FillerBase::setTcal(const String& tcaltime,
                     const Vector<Float>& tcal)
{
    uInt id = table_->tcal().addEntry(tcaltime, tcal);
    RecordFieldPtr<uInt> mcalidCol(row_.record(), "TCAL_ID");
    *mcalidCol = id;
}

void FillerBase::setScanRate(const Vector<Double>& srate)
{
    RecordFieldPtr<Array<Double> > srateCol(row_.record(), "SCANRATE");
    *srateCol = srate;
}

void FillerBase::setReferenceBeam(Int beamno)
{
  RecordFieldPtr<Int> rbCol(row_.record(), "REFBEAMNO");
  *rbCol = beamno;
}

void FillerBase::setSource(const std::string& name, Int type,
                           const std::string& fieldname,
                           const Vector<Double>& dir,
                           const Vector<Double>& propermot,
                           Double velocity)
{
    RecordFieldPtr<String> srcnCol(row_.record(), "SRCNAME");
    *srcnCol = name;
    RecordFieldPtr<Int> srctCol(row_.record(), "SRCTYPE");
    *srctCol = type;
    RecordFieldPtr<String> fieldnCol(row_.record(), "FIELDNAME");
    *fieldnCol = fieldname;
    RecordFieldPtr<Array<Double> > spmCol(row_.record(), "SRCPROPERMOTION");
    *spmCol = propermot;
    RecordFieldPtr<Array<Double> > sdirCol(row_.record(), "SRCDIRECTION");
    *sdirCol = dir;
    RecordFieldPtr<Double> svelCol(row_.record(), "SRCVELOCITY");
    *svelCol = velocity;
}

void FillerBase::commitRow()
{
  table_->table().addRow();
  row_.put(table_->table().nrow()-1);
}

void FillerBase::setWeather2(Float temperature, 
                             Float pressure,
                             Float humidity,
                             Float windspeed, 
                             Float windaz)
{
  uInt id ;
  Table tab = table_->weather().table() ;
  Table subt = tab( tab.col("TEMPERATURE") == temperature \
                    && tab.col("PRESSURE") == pressure \
                    && tab.col("HUMIDITY") == humidity \
                    && tab.col("WINDSPEED") == windspeed \
                    && tab.col("WINDAZ") == windaz ) ;
  Int nrow = tab.nrow() ;
  Int nrowSel = subt.nrow() ;
  if ( nrowSel == 0 ) {
    tab.addRow( 1, True ) ;
    TableRow row( tab ) ;
    TableRecord &rec = row.record() ;
    RecordFieldPtr<casa::uInt> rfpi ;
    rfpi.attachToRecord( rec, "ID" ) ;
    *rfpi = (uInt)nrow ;
    RecordFieldPtr<casa::Float> rfp ;
    rfp.attachToRecord( rec, "TEMPERATURE" ) ;
    *rfp = temperature ;
    rfp.attachToRecord( rec, "PRESSURE" ) ;
    *rfp = pressure ;
    rfp.attachToRecord( rec, "HUMIDITY" ) ;
    *rfp = humidity ;
    rfp.attachToRecord( rec, "WINDSPEED" ) ;
    *rfp = windspeed ;
    rfp.attachToRecord( rec, "WINDAZ" ) ;
    *rfp = windaz ;
    row.put( nrow, rec ) ;
    id = (uInt)nrow ;
  }
  else {
    ROTableColumn tc( subt, "ID" ) ;
    id = tc.asuInt( 0 ) ;
  }
  RecordFieldPtr<uInt> mweatheridCol(row_.record(), "WEATHER_ID");
  *mweatheridCol = id;
}

void FillerBase::setTcal2(const String& tcaltime,
                          const Vector<Float>& tcal)
{
  uInt id ;
  Table tab = table_->tcal().table() ;
  Int nrow = tab.nrow() ;
  Vector<uInt> rowList( 0 ) ;
  ArrayColumn<Float> tcalCol( tab, "TCAL" ) ;
  TableColumn timeCol( tab, "TIME" ) ;
  TableColumn idCol( tab, "ID" ) ;
  uInt nelem = tcal.nelements() ;
  for ( Int irow = 0 ; irow < nrow ; irow++ ) {
    if ( tcalCol.shape(irow)[0] == nelem ) {
      rowList.resize( rowList.nelements()+1, True ) ;
      rowList[rowList.nelements()-1] = irow ;
    }
  }
  
  //cout << "rowList = " << rowList << endl ;

  if ( rowList.nelements() == 0 ) {
    // add new row
    tab.addRow( 1 ) ;
    //cout << "tab.nrow() = " << tab.nrow() << endl ;
    tcalCol.put( nrow, tcal ) ;
    timeCol.putScalar( nrow, tcaltime ) ; 
    id = (uInt)nrow ;
    idCol.putScalar( nrow, id ) ;
  }
  else {
    uInt ichan = 0 ;
    while ( rowList.nelements() > 1 && ichan < nelem ) {
      Vector<uInt> tmp = rowList.copy() ;
      rowList.resize( 0 ) ;
      for ( uInt irow = 0 ; irow < tmp.nelements() ; irow++ ) {
        Vector<Float> t = tcalCol( tmp[irow] ) ;
        if ( t[ichan] == tcal[ichan] ) {
          rowList.resize( rowList.nelements()+1, True ) ;
          rowList[rowList.nelements()-1] = irow ;
        }
      }
      ichan++ ;
    }
    
    //cout << "updated rowList = " << rowList << endl ;

    if ( rowList.nelements() == 0 ) {
      // add new row
      tab.addRow( 1, True ) ;
      //cout << "tab.nrow() = " << tab.nrow() << endl ;
      tcalCol.put( nrow, tcal ) ;
      timeCol.putScalar( nrow, tcaltime ) ; 
      id = (uInt)nrow ;
      idCol.putScalar( nrow, id ) ;
    }
    else {
      Vector<Float> t = tcalCol( rowList[0] ) ;
      if ( allEQ( t, tcal ) ) {
        ROTableColumn tc( tab, "ID" ) ;
        id = tc.asuInt( rowList[0] ) ;
      }
      else {
        // add new row
        tab.addRow( 1, True ) ;
        //cout << "tab.nrow() = " << tab.nrow() << endl ;
        tcalCol.put( nrow, tcal ) ;
        timeCol.putScalar( nrow, tcaltime ) ; 
        id = (uInt)nrow ;
        idCol.putScalar( nrow, id ) ;
      }
    }
  }
  
  RecordFieldPtr<uInt> mcalidCol(row_.record(), "TCAL_ID");
  *mcalidCol = id;
}

};
