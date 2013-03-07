//#---------------------------------------------------------------------------
//# ASTEDataset.cc: Class for ASTE dataset.
//#---------------------------------------------------------------------------
//# Copyright (C) 2000-2006
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id$
//#---------------------------------------------------------------------------
//# Original: 2008/10/30, Takeshi Nakazato, NAOJ
//#---------------------------------------------------------------------------

#include <atnf/PKSIO/ASTEDataset.h>

#include <iostream>
#include <cstring>

using namespace std ;

// constructor 
ASTEDataset::ASTEDataset( string name )
  : NRODataset( name ) 
{}

// destructor 
ASTEDataset::~ASTEDataset() 
{}

// data initialization
void ASTEDataset::initialize()
{
  int arymax = arrayMax() ;

  // it must be called
  NRODataset::initialize() ;

  // additional initialization
  datasize_ += sizeof( char ) * arymax * 16 // RX
    + sizeof( double ) * arymax * 6         // HPBW, EFFA, EFFB, EFFL, EFSS GAIN
    + sizeof( char ) * arymax * 4           // HORN
    + sizeof( char ) * arymax * 4           // POLTP
    + sizeof( double ) * arymax * 3         // POLDR, POLAN, DFRQ
    + sizeof( char ) * arymax * 4           // SIDBID
    + sizeof( int ) * arymax * 3            // REFN, IPINT, MULTN
    + sizeof( double ) * arymax             // MLTSCF
    + sizeof( char ) * arymax * 8           // LAGWIND
    + sizeof( double ) * arymax * 3         // BEBW, BERES, CHWID
    + sizeof( int ) * arymax * 2            // ARRY, NFCAL
    + sizeof( double ) * arymax             // F0CAL
    + sizeof( double ) * arymax * 10 * 3    // FQCAL, CHCAL, CWCAL
    + sizeof( char ) * 116                       // CDMY1
    + sizeof( double ) * arymax ;           // DSBFC
}

int ASTEDataset::fillHeader( int sameEndian )
{
  LogIO os( LogOrigin( "ASTEDataset", "fillHeader()", WHERE ) ) ;

  int arymax = arrayMax();

  string str4( 4, ' ' ) ;
  string str8( 8, ' ' ) ;
  string str16( 16, ' ' ) ;
  string str18( 18, ' ' ) ;
  string str24( 24, ' ' ) ;
  string str40( 40, ' ' ) ;
  string str120( 120, ' ' ) ;
  char c4[4] ;
  char c8[8] ;
  char c16[16] ;
  char c18[18] ;
  char c24[24] ;
  char c40[40] ;
  char c120[120] ;
  // make sure file pointer points a beginning of the file
  fseek( fp_, 0, SEEK_SET ) ;

  // read data header
  strcpy( c8, str8.c_str() ) ;
  if ( readHeader( c8, 8 ) == -1 ) {
    os << LogIO::WARN << "Error while reading data LOFIL." << LogIO::POST ;
    return -1 ;
  }
  LOFIL = string( c8 ) ;
  // DEBUG
  //cout << "LOFIL = " << LOFIL << endl ;
  //
  strcpy( c8, str8.c_str() ) ;
  if ( readHeader( c8, 8 ) == -1 ) {
    os << LogIO::WARN << "Error while reading data VER." << LogIO::POST ;
    return -1 ;
  }
  VER = string( c8 ) ;
  // DEBUG
  //cout << "VER = " << VER << endl ;
  //
  strcpy( c16, str16.c_str() ) ;
  if ( readHeader( c16, 16 ) == -1 ) {
    os << LogIO::WARN << "Error while reading data GROUP." << LogIO::POST ;
    return -1 ;
  }
  GROUP = string( c16 ) ;
  // DEBUG
  //cout << "GROUP = " << GROUP << endl ;
  //
  strcpy( c16, str16.c_str() ) ;
  if ( readHeader( c16, 16 ) == -1 ) {
    os << LogIO::WARN << "Error while reading data PROJ." << LogIO::POST ;
    return -1 ;
  }
  PROJ = string( c16 ) ;
  // DEBUG
  //cout << "PROJ = " << PROJ << endl ;
  //
  strcpy( c24, str24.c_str() ) ;
  if ( readHeader( c24, 24 ) == -1 ) {
    os << LogIO::WARN << "Error while reading data SCHED." << LogIO::POST ;
    return -1 ;
  }
  SCHED = string( c24 ) ;
  // DEBUG
  //cout << "SCHED = " << SCHED << endl ;
  //
  strcpy( c40, str40.c_str() ) ;
  if ( readHeader( c40, 40 ) == -1 ) {
    os << LogIO::WARN << "Error while reading data OBSVR." << LogIO::POST ;
    return -1 ;
  }  
  OBSVR = string( c40 ) ;
  // DEBUG
  //cout << "OBSVR = " << OBSVR << endl ;
  //
  strcpy( c16, str16.c_str() ) ;
  if ( readHeader( c16, 16 ) == -1 ) {
    os << LogIO::WARN << "Error while reading data LOSTM." << LogIO::POST ;
    return -1 ;
  }
  LOSTM = string( c16 ) ;
  // DEBUG
  //cout << "LOSTM = " << LOSTM << endl ;
  //
  strcpy( c16, str16.c_str() ) ;
  if ( readHeader( c16, 16 ) == -1 ) {
    os << LogIO::WARN << "Error while reading data LOETM." << LogIO::POST ;
    return -1 ;
  }
  LOETM = string( c16 ) ;
  // DEBUG
  //cout << "LOETM = " << LOETM << endl ;
  //
  if ( readHeader( ARYNM, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data ARYNM." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "ARYNM = " << ARYNM << endl ;
  //
  if ( readHeader( NSCAN, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data NSCAN." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "NSCAN = " << NSCAN << endl ;
  //
  strcpy( c120, str120.c_str() ) ;
  if ( readHeader( c120, 120 ) == -1 ) {
    os << LogIO::WARN << "Error while reading data TITLE." << LogIO::POST ;
    return -1 ;
  }
  TITLE = string( c120 ) ;
  // DEBUG
  //cout << "TITLE = " << TITLE << endl ;
  //
  strcpy( c16, str16.c_str() ) ;
  if ( readHeader( c16, 16 ) == -1 ) {
    os << LogIO::WARN << "Error while reading data OBJ." << LogIO::POST ;
    return -1 ;
  }
  OBJ = string( c16 ) ;
  // DEBUG
  //cout << "OBJ = " << OBJ << endl ;
  //
  strcpy( c8, str8.c_str() ) ;
  if ( readHeader( c8, 8 ) == -1 ) {
    os << LogIO::WARN << "Error while reading data EPOCH." << LogIO::POST ;
    return -1 ;
  }
  EPOCH = string( c8 ) ;
  // DEBUG
  //cout << "EPOCH = " << EPOCH << endl ;
  //
  if ( readHeader( RA0, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data RA0." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "RA0 = " << RA0 << endl ;
  //
  if ( readHeader( DEC0, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data DEC0." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "DEC0 = " << DEC0 << endl ;
  //
  if ( readHeader( GLNG0, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data GLNG0." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "GLNG0 = " << GLNG0 << endl ;
  //
  if ( readHeader( GLAT0, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data GLAT0." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "GLAT0 = " << GLAT0 << endl ;
  //
  if ( readHeader( NCALB, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data NCALB." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "NCALB = " << NCALB << endl ;
  //
  if ( readHeader( SCNCD, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data SCNCD." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "SCNCD = " << SCNCD << endl ;
  //
  strcpy( c120, str120.c_str() ) ;
  if ( readHeader( c120, 120 ) == -1 ) {
    os << LogIO::WARN << "Error while reading data SCMOD." << LogIO::POST ;
    return -1 ;
  }
  SCMOD = string( c120 ) ;
  // DEBUG
  //cout << "SCMOD = " << SCMOD << endl ;
  //
  if ( readHeader( URVEL, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data URVEL." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "URVEL = " << URVEL << endl ;
  //
  strcpy( c4, str4.c_str() ) ;
  if ( readHeader( c4, 4 ) == -1 ) {
    os << LogIO::WARN << "Error while reading data VREF." << LogIO::POST ;
    return -1 ;
  }
  VREF = string( c4 ) ;
  // DEBUG
  //cout << "VREF = " << VREF << endl ;
  //
  strcpy( c4, str4.c_str() ) ;
  if ( readHeader( c4, 4 ) == -1 ) {
    os << LogIO::WARN << "Error while reading data VDEF." << LogIO::POST ;
    return -1 ;
  }
  VDEF = string( c4 ) ;
  // DEBUG
  //cout << "VDEF = " << VDEF << endl ;
  //
  strcpy( c8, str8.c_str() ) ;
  if ( readHeader( c8, 8 ) == -1 ) {
    os << LogIO::WARN << "Error while reading data SWMOD." << LogIO::POST ;
    return -1 ;
  }
  SWMOD = string( c8 ) + "::OTF" ;
  // DEBUG
  //cout << "SWMOD = " << SWMOD << endl ;
  //
  if ( readHeader( FRQSW, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data FRQSW." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "FRQSW = " << FRQSW << endl ;
  //
  if ( readHeader( DBEAM, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data DBEAM." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "DBEAM = " << DBEAM << endl ;
  //
  if ( readHeader( MLTOF, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data MLTOF." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "MLTOF = " << MLTOF << endl ;
  //
  if ( readHeader( CMTQ, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data CMTQ." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "CMTQ = " << CMTQ << endl ;
  //
  if ( readHeader( CMTE, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data CMTE." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "CMTE = " << CMTE << endl ;
  //
  if ( readHeader( CMTSOM, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data CMTSOM." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "CMTSOM = " << CMTSOM << endl ;
  //
  if ( readHeader( CMTNODE, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data CMTNODE." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "CMTNODE = " << CMTNODE << endl ;
  //
  if ( readHeader( CMTI, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data CMTI." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "CMTI = " << CMTI << endl ;
  //
  strcpy( c24, str24.c_str() ) ;
  if ( readHeader( c24, 24 ) == -1 ) {
    os << LogIO::WARN << "Error while reading data CMTTM." << LogIO::POST ;
    return -1 ;
  }
  CMTTM = string( c24 ) ;
  // DEBUG
  //cout << "CMTTM = " << CMTTM << endl ;
  //
  if ( readHeader( SBDX, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data SBDX." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "SBDX = " << SBDX << endl ;
  //
  if ( readHeader( SBDY, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data SBDY." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "SBDY = " << SBDY << endl ;
  //
  if ( readHeader( SBDZ1, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data SBDZ1." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "SBDZ1 = " << SBDZ1 << endl ;
  //
  if ( readHeader( SBDZ2, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data SBDZ2." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "SBDZ2 = " << SBDZ2 << endl ;
  //
  if ( readHeader( DAZP, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data DAZP." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "DAZP = " << DAZP << endl ;
  //
  if ( readHeader( DELP, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data DELP." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "DELP = " << DELP << endl ;
  //
  if ( readHeader( CHBIND, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data CHBIND." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "CHBIND = " << CHBIND << endl ;
  //
  if ( readHeader( NUMCH, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data NUMCH." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "NUMCH = " << NUMCH << endl ;
  //
  if ( readHeader( CHMIN, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data CHMIN." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "CHMIN = " << CHMIN << endl ;
  //
  if ( readHeader( CHMAX, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data CHMAX." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "CHMAX = " << CHMAX << endl ;
  //
  if ( readHeader( ALCTM, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data ALCTM." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "ALCTM = " << ALCTM << endl ;
  //
  if ( readHeader( IPTIM, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data IPTIM." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "IPTIM = " << IPTIM << endl ;
  //
  if ( readHeader( PA, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data PA." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "PA = " << PA << endl ;
  //
  for ( int i = 0 ; i < arymax ; i++ ) {
    strcpy( c18, str18.c_str() ) ;
    if ( readHeader( c18, 16 ) == -1 ) {
      os << LogIO::WARN << "Error while reading data RX[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
    c18[16] = '\0' ;
    RX[i] = string( c18 ) ;
  }
  // DEBUG
//   nro_debug_output( "RX", arymax, RX ) ;
  // 
  for ( int i = 0 ; i < arymax ; i++ ) {
    if ( readHeader( HPBW[i], sameEndian ) == -1 ) {
      os << LogIO::WARN << "Error while reading data HPBW[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
  }
  // DEBUG
//   nro_debug_output( "HPBW", arymax, HPBW ) ;
  // 
  for ( int i = 0 ; i < arymax ; i++ ) {
    if ( readHeader( EFFA[i], sameEndian ) == -1 ) {
      os << LogIO::WARN << "Error while reading data EFFA[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
  }
  // DEBUG
//   nro_debug_output( "EFFA", arymax, EFFA ) ;
  // 
  for ( int i = 0 ; i < arymax ; i++ ) {
    if ( readHeader( EFFB[i], sameEndian ) == -1 ) {
      os << LogIO::WARN << "Error while reading data EFFB[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
  }
  // DEBUG
//   nro_debug_output( "EFFB", arymax, EFFB ) ;
  // 
  for ( int i = 0 ; i < arymax ; i++ ) {
    if ( readHeader( EFFL[i], sameEndian ) == -1 ) {
      os << LogIO::WARN << "Error while reading data EFFL[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
  }
  // DEBUG
//   nro_debug_output( "EFFL", arymax, EFFL ) ;
  // 
  for ( int i = 0 ; i < arymax ; i++ ) {
    if ( readHeader( EFSS[i], sameEndian ) == -1 ) {
      os << LogIO::WARN << "Error while reading data EFSS[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
  }
  // DEBUG
//   nro_debug_output( "EFSS", arymax, EFSS ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    if ( readHeader( GAIN[i], sameEndian ) == -1 ) {
      os << LogIO::WARN << "Error while reading data GAIN[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
  }
  // DEBUG
//   nro_debug_output( "GAIN", arymax, GAIN ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    strcpy( c4, str4.c_str() ) ;
    if ( readHeader( c4, 4 ) == -1 ) {
      os << LogIO::WARN << "Error while reading data HORN[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
    HORN[i] = string( c4 ) ;
  }
  // DEBUG
//   nro_debug_output( "HORN", arymax, HORN ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    strcpy( c4, str4.c_str() ) ;
    if ( readHeader( c4, 4 ) == -1 ) {
      os << LogIO::WARN << "Error while reading data POLTP[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
    POLTP[i] = string( c4 ) ;
  }
  // DEBUG
//   nro_debug_output( "POLTP", arymax, POLTP ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    if ( readHeader( POLDR[i], sameEndian ) == -1 ) {
      os << LogIO::WARN << "Error while reading data POLDR[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
  }
  // DEBUG
//   nro_debug_output( "POLDR", arymax, POLDR ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    if ( readHeader( POLAN[i], sameEndian ) == -1 ) {
      os << LogIO::WARN << "Error while reading data POLAN[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
  }
  // DEBUG
//   nro_debug_output( "POLAN", arymax, POLAN ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    if ( readHeader( DFRQ[i], sameEndian ) == -1 ) {
      os << LogIO::WARN << "Error while reading data DFRQ[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
  }
  // DEBUG
//   nro_debug_output( "DFRQ", arymax, DFRQ ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    strcpy( c4, str4.c_str() ) ;
    if ( readHeader( c4, 4 ) == -1 ) {
      os << LogIO::WARN << "Error while reading data SIDBD[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
    SIDBD[i] = string( c4 ) ;
  }
  // DEBUG
//   nro_debug_output( "SIDBD", arymax, SIDBD ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    if ( readHeader( REFN[i], sameEndian ) == -1 ) {
      os << LogIO::WARN << "Error while reading data REFN[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
  }
  // DEBUG
//   nro_debug_output( "REFN", arymax, REFN ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    if ( readHeader( IPINT[i], sameEndian ) == -1 ) {
      os << LogIO::WARN << "Error while reading data IPINT[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
  }
  // DEBUG
//   nro_debug_output( "IPINT", arymax, IPINT ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    if ( readHeader( MULTN[i], sameEndian ) == -1 ) {
      os << LogIO::WARN << "Error while reading data MULTN[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
  }
  // DEBUG
//   nro_debug_output( "MULTN", arymax, MULTN ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    if ( readHeader( MLTSCF[i], sameEndian ) == -1 ) {
      os << LogIO::WARN << "Error while reading data MLTSCF[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
  }
  // DEBUG
//   nro_debug_output( "MLTSCF", arymax, MLTSCF ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    strcpy( c8, str8.c_str() ) ;
    if ( readHeader( c8, 8 ) == -1 ) {
      os << LogIO::WARN << "Error while reading data LAGWIND[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
    LAGWIND[i] = string( c8 ) ;
  }
  // DEBUG
//   nro_debug_output( "LAGWIND", arymax, LAGWIND ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    if ( readHeader( BEBW[i], sameEndian ) == -1 ) {
      os << LogIO::WARN << "Error while reading data BEBW[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
  }
  // DEBUG
//   nro_debug_output( "BEBW", arymax, BEBW ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    if ( readHeader( BERES[i], sameEndian ) == -1 ) {
      os << LogIO::WARN << "Error while reading data BERES[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
  }
  // DEBUG
//   nro_debug_output( "BERES", arymax, BERES ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    if ( readHeader( CHWID[i], sameEndian ) == -1 ) {
      os << LogIO::WARN << "Error while reading data CHWID[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
  }
  // DEBUG
//   nro_debug_output( "CHWID", arymax, CHWID ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    if ( readHeader( ARRY[i], sameEndian ) == -1 ) {
      os << LogIO::WARN << "Error while reading data ARRY[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
  }
  // DEBUG
//   nro_debug_output( "ARRY", arymax, ARRY ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    if ( readHeader( NFCAL[i], sameEndian ) == -1 ) {
      os << LogIO::WARN << "Error while reading data NFCAL[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
  }
  // DEBUG
//   nro_debug_output( "NFCAL", arymax, NFCAL ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    if ( readHeader( F0CAL[i], sameEndian ) == -1 ) {
      os << LogIO::WARN << "Error while reading data F0CAL[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
  }
  // DEBUG
//   nro_debug_output( "F0CAL", arymax, F0CAL ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    for ( int j = 0 ; j < 10 ; j++ ) {
      if ( readHeader( FQCAL[i][j], sameEndian ) == -1 ) {
        os << LogIO::WARN << "Error while reading data FQCAL[" << i << "][" << j << "]." << LogIO::POST ;
        return -1 ;
      }
    }
  }
  // DEBUG
//   nro_debug_output( "FQCAL", arymax, 10,  FQCAL ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    for ( int j = 0 ; j < 10 ; j++ ) {
      if ( readHeader( CHCAL[i][j], sameEndian ) == -1 ) {
        os << LogIO::WARN << "Error while reading data CHCAL[" << i << "][" << j << "]." << LogIO::POST ;
        return -1 ;
      }
    }
  }
  // DEBUG
//   nro_debug_output( "CHCAL", arymax, 10, CHCAL ) ;
  // 
  for ( int i= 0 ; i < arymax ; i++) {
    for ( int j = 0 ; j < 10 ; j++ ) {
      if ( readHeader( CWCAL[i][j], sameEndian ) == -1 ) {
        os << LogIO::WARN << "Error while reading data CWCAL[" << i << "][" << j << "]." << LogIO::POST ;
        return -1 ;
      }
    }
  }
  // DEBUG
//   nro_debug_output( "CWCAL", arymax, 10, CWCAL ) ;
  // 
  if ( readHeader( SCNLEN, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data SCNLEN." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "SCNLEN = " << SCNLEN << endl ;
  //
  if ( readHeader( SBIND, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data SBIND." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "SBIND = " << SBIND << endl ;
  //
  if ( readHeader( IBIT, sameEndian ) == -1 ) {
    os << LogIO::WARN << "Error while reading data IBIT." << LogIO::POST ;
    return -1 ;
  }
  // DEBUG
  //cout << "IBIT = " << IBIT << endl ;
  //
  strcpy( c8, str8.c_str() ) ;
  if ( readHeader( c8, 8 ) == -1 ) {
    os << LogIO::WARN << "Error while reading data SITE." << LogIO::POST ;
    return -1 ;
  }
  SITE = string( c8 ) ;
  // DEBUG
  //cout << "SITE = " << SITE << endl ;
  //
  strcpy( c120, str120.c_str() ) ;
  if ( readHeader( c120, 116 ) == -1 ) {
    os << LogIO::WARN << "Error while reading data CDMY1." << LogIO::POST ;
    return -1 ;
  }
  c120[116] = '\0' ;
  CDMY1 = string( c120 ) ;
  // DEBUG
  //cout << "CDMY1 = " << CDMY1 << endl ;
  //
  for ( int i = 0 ; i < arymax ; i++ ) {
    if ( readHeader( DSBFC[i], sameEndian ) == -1 ) {
      os << LogIO::WARN << "Error while reading data DSBFC[" << i << "]." << LogIO::POST ;
      return -1 ;
    }
  }
  // DEBUG
//   nro_debug_output( "DSBFC", arymax, DSBFC ) ;
  // 

  return 0 ;
}
