//#---------------------------------------------------------------------------
//# SDLineFinder.cc: A class for automated spectral line search
//#--------------------------------------------------------------------------
//# Copyright (C) 2004
//# ATNF
//#
//# This program is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by the Free
//# Software Foundation; either version 2 of the License, or (at your option)
//# any later version.
//#
//# This program is distributed in the hope that it will be useful, but
//# WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
//# Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with this program; if not, write to the Free Software Foundation, Inc.,
//# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning this software should be addressed as follows:
//#        Internet email: Malte.Marquarding@csiro.au
//#        Postal address: Malte Marquarding,
//#                        Australia Telescope National Facility,
//#                        P.O. Box 76,
//#                        Epping, NSW, 2121,
//#                        AUSTRALIA
//#
//# $Id:
//#---------------------------------------------------------------------------


// ASAP
#include "SDLineFinder.h"

// STL
#include <functional>
#include <algorithm>
#include <iostream>

using namespace asap;
using namespace casa;
using namespace std;
using namespace boost::python;

///////////////////////////////////////////////////////////////////////////////
//
// LFRunningMean - a running mean algorithm for line detection
//
//

namespace asap {

// An auxiliary class implementing one pass of the line search algorithm,
// which uses a running mean. We define this class here because it is
// used in SDLineFinder only. The incapsulation of this code into a separate
// class will provide a possibility to add new algorithms with minor changes
class LFRunningMean {
   // The input data to work with. Use reference symantics to avoid
   // an unnecessary copying   
   const casa::Vector<casa::Float>  &spectrum; // a buffer for the spectrum
   const casa::Vector<casa::Bool>   &mask; // associated mask
   const std::pair<int,int>         &edge; // start and stop+1 channels
                                           // to work with
   
   // statistics for running mean filtering
   casa::Float sum;       // sum of fluxes
   casa::Float sumsq;     // sum of squares of fluxes
   int box_chan_cntr;     // actual number of channels in the box
   int max_box_nchan;     // maximum allowed number of channels in the box
                          // (calculated from boxsize and actual spectrum size)

   // temporary line edge channels and flag, which is True if the line
   // was detected in the previous channels.
   std::pair<int,int> cur_line;
   casa::Bool is_detected_before;
   int  min_nchan;                         // A minimum number of consequtive
                                           // channels, which should satisfy
					   // the detection criterion, to be
					   // a detection
   casa::Float threshold;                  // detection threshold - the 
                                           // minimal signal to noise ratio
public:
   // set up the object with the references to actual data
   // as well as the detection criterion (min_nchan and threshold, see above)
   // and the number of channels in the running box
   LFRunningMean(const casa::Vector<casa::Float>  &in_spectrum,
                 const casa::Vector<casa::Bool>   &in_mask,
		 const std::pair<int,int>         &in_edge,
		 int in_max_box_nchan, 
		 int in_min_nchan = 3,
		 casa::Float in_threshold = 5);
		 
   // replace the detection criterion
   void setCriterion(int in_min_nchan, casa::Float in_threshold) throw();

   // find spectral lines and add them into list
   void findLines(std::list<pair<int,int> > &lines) throw(casa::AipsError);
   
protected:
   // supplementary function to control running mean calculations.
   // It adds a specified channel to the running mean box and
   // removes (ch-maxboxnchan+1)'th channel from there
   // Channels, for which the mask is false or index is beyond the
   // allowed range, are ignored
   void advanceRunningBox(int ch) throw(casa::AipsError);
   

   // test a channel against current running mean & rms
   // if channel specified is masked out or beyond the allowed indexes,
   // false is returned
   casa::Bool testChannel(int ch) const
                     throw(std::exception, casa::AipsError);

   // process a channel: update curline and is_detected before and
   // add a new line to the list, if necessary using processCurLine()
   void processChannel(std::list<pair<int,int> > &lines,
                       int ch) throw(casa::AipsError);

   // process the interval of channels stored in curline
   // if it satisfies the criterion, add this interval as a new line
   void processCurLine(std::list<pair<int,int> > &lines) throw(casa::AipsError);

};
} // namespace asap

//
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
//
// LFRunningMean - a running mean algorithm for line detection
//
//

// set up the object with the references to actual data
// as well as the detection criterion (min_nchan and threshold, see above)
// and the number of channels in the running box
LFRunningMean::LFRunningMean(const casa::Vector<casa::Float>  &in_spectrum,
                             const casa::Vector<casa::Bool>   &in_mask,
     	                     const std::pair<int,int>         &in_edge,
			     int in_max_box_nchan, 
     	                     int in_min_nchan, casa::Float in_threshold) :
        spectrum(in_spectrum), mask(in_mask), edge(in_edge),
	max_box_nchan(in_max_box_nchan),
	min_nchan(in_min_nchan),threshold(in_threshold) {}

// replace the detection criterion
void LFRunningMean::setCriterion(int in_min_nchan, casa::Float in_threshold)
                                 throw()
{
  min_nchan=in_min_nchan;
  threshold=in_threshold;
}


// supplementary function to control running mean calculations.
// It adds a specified channel to the running mean box and
// removes (ch-max_box_nchan+1)'th channel from there
// Channels, for which the mask is false or index is beyond the
// allowed range, are ignored
void LFRunningMean::advanceRunningBox(int ch) throw(AipsError)
{
  if (ch>=edge.first && ch<edge.second)
      if (mask[ch]) { // ch is a valid channel
          ++box_chan_cntr;
          sum+=spectrum[ch];
          sumsq+=square(spectrum[ch]);          
      }
  int ch2remove=ch-max_box_nchan;
  if (ch2remove>=edge.first && ch2remove<edge.second)
      if (mask[ch2remove]) { // ch2remove is a valid channel
          --box_chan_cntr;
          sum-=spectrum[ch2remove];
          sumsq-=square(spectrum[ch2remove]);  
      }
}

// test a channel against current running mean & rms
// if channel specified is masked out or beyond the allowed indexes,
// false is returned
Bool LFRunningMean::testChannel(int ch) const throw(exception, AipsError)
{
  if (ch<edge.first || ch>=edge.second) return False;
  if (!mask[ch]) return False;
  DebugAssert(box_chan_cntr, AipsError);
  Float mean=sum/Float(box_chan_cntr);
  Float variance=sqrt(sumsq/Float(box_chan_cntr)-square(mean));
  /*
  if (ch>3900 && ch<4100)
    cout<<"Tested "<<ch<<" mean="<<mean<<" variance="<<variance<<" sp-mean="<<spectrum[ch]-mean<<endl;
  */
  return fabs(spectrum[ch]-mean)>=threshold*variance;
}

// process a channel: update cur_line and is_detected before and
// add a new line to the list, if necessary
void LFRunningMean::processChannel(std::list<pair<int,int> > &lines,
                                   int ch) throw(casa::AipsError)
{
  try {
       if (testChannel(ch)) {
           if (is_detected_before)
               cur_line.second=ch+1;
	   else {
	       is_detected_before=True;
	       cur_line.first=ch;
	       cur_line.second=ch+1;
	   }
       } else processCurLine(lines);   
  }
  catch (const AipsError &ae) {
      throw;
  }  
  catch (const exception &ex) {
      throw AipsError(String("SDLineFinder::processChannel - STL error: ")+ex.what());
  }
}

// process the interval of channels stored in cur_line
// if it satisfies the criterion, add this interval as a new line
void LFRunningMean::processCurLine(std::list<pair<int,int> > &lines)
                                   throw(casa::AipsError)
{
  try {
       if (is_detected_before) {	              
           if (cur_line.second-cur_line.first>min_nchan) {
	       // it was a detection. We need to change the list
	       Bool add_new_line=False;
	       if (lines.size()) { 
	           for (int i=lines.back().second;i<cur_line.first;++i)
		        if (mask[i]) { // one valid channel in between
			        //  means that we deal with a separate line
			    add_new_line=True;
			    break;
			}
	       } else add_new_line=True; 
	       if (add_new_line) 
	           lines.push_back(cur_line);
               else lines.back().second=cur_line.second;		   
	   }
	   is_detected_before=False;
       }      
  }
  catch (const AipsError &ae) {
      throw;
  }  
  catch (const exception &ex) {
      throw AipsError(String("SDLineFinder::processCurLine - STL error: ")+ex.what());
  }
}

// find spectral lines and add them into list
void LFRunningMean::findLines(std::list<pair<int,int> > &lines)
                        throw(casa::AipsError)
{
  const int minboxnchan=4;

  // fill statistics for initial box
  box_chan_cntr=0; // no channels are currently in the box
  sum=0;           // initialize statistics
  sumsq=0;
  int initial_box_ch=edge.first;
  for (;initial_box_ch<edge.second && box_chan_cntr<max_box_nchan;
        ++initial_box_ch)
       advanceRunningBox(initial_box_ch);
    
  if (initial_box_ch==edge.second)       
      throw AipsError("LFRunningMean::findLines - too much channels are masked");

  // actual search algorithm
  is_detected_before=False;

  if (box_chan_cntr>=minboxnchan) 
      // there is a minimum amount of data. We can search in the
      // half of the initial box   
      for (int n=edge.first;n<initial_box_ch-max_box_nchan/2;++n)
           processChannel(lines,n);           	   
  
  // now the box can be moved. n+max_box_nchan/2 is a new index which haven't
  // yet been included in the running mean.
  for (int n=initial_box_ch-max_box_nchan/2;n<edge.second;++n) {
      advanceRunningBox(n+max_box_nchan/2); // update running mean & variance
      if (box_chan_cntr>=minboxnchan) // have enough data to process
          processChannel(lines,n);
      else processCurLine(lines); // just finish what was accumulated before
  }  
}

//
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//
// SDLineFinder::IntersectsWith  -  An auxiliary object function to test
// whether two lines have a non-void intersection
//


// line1 - range of the first line: start channel and stop+1
SDLineFinder::IntersectsWith::IntersectsWith(const std::pair<int,int> &in_line1) :
                          line1(in_line1) {}


// return true if line2 intersects with line1 with at least one
// common channel, and false otherwise
// line2 - range of the second line: start channel and stop+1
bool SDLineFinder::IntersectsWith::operator()(const std::pair<int,int> &line2)
                          const throw()
{
  if (line2.second<line1.first) return false; // line2 is at lower channels
  if (line2.first>line1.second) return false; // line2 is at upper channels
  return true; // line2 has an intersection or is adjacent to line1
}

//
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//
// SDLineFinder::BuildUnion - An auxiliary object function to build a union
// of several lines to account for a possibility of merging the nearby lines
//

// set an initial line (can be a first line in the sequence)
SDLineFinder::BuildUnion::BuildUnion(const std::pair<int,int> &line1) :
                             temp_line(line1) {}

// update temp_line with a union of temp_line and new_line
// provided there is no gap between the lines
void SDLineFinder::BuildUnion::operator()(const std::pair<int,int> &new_line)
                                   throw()
{
  if (new_line.first<temp_line.first) temp_line.first=new_line.first;
  if (new_line.second>temp_line.second) temp_line.second=new_line.second;
}

// return the result (temp_line)
const std::pair<int,int>& SDLineFinder::BuildUnion::result() const throw()
{
  return temp_line;
}

//
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//
// SDLineFinder::LaterThan - An auxiliary object function to test whether a
// specified line is at lower spectral channels (to preserve the order in
// the line list)
//

// setup the line to compare with
SDLineFinder::LaterThan::LaterThan(const std::pair<int,int> &in_line1) :
                         line1(in_line1) {}

// return true if line2 should be placed later than line1
// in the ordered list (so, it is at greater channel numbers)
bool SDLineFinder::LaterThan::operator()(const std::pair<int,int> &line2)
                          const throw()
{
  if (line2.second<line1.first) return false; // line2 is at lower channels
  if (line2.first>line1.second) return true; // line2 is at upper channels
  
  // line2 intersects with line1. We should have no such situation in
  // practice
  return line2.first>line1.first;
}

//
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
//
// SDLineFinder  -  a class for automated spectral line search
//
//

SDLineFinder::SDLineFinder() throw() : edge(0,0)
{
  // detection threshold - the minimal signal to noise ratio
  threshold=3.; // 3 sigma is a default
  box_size=1./16.; // default box size for running mean calculations is
                  // 1/16 of the whole spectrum
  // A minimum number of consequtive channels, which should satisfy
  // the detection criterion, to be a detection
  min_nchan=3;     // default is 3 channels
}

SDLineFinder::~SDLineFinder() throw(AipsError) {}

// set scan to work with (in_scan parameter), associated mask (in_mask
// parameter) and the edge channel rejection (in_edge parameter)
//   if in_edge has zero length, all channels chosen by mask will be used
//   if in_edge has one element only, it represents the number of
//      channels to drop from both sides of the spectrum
//   in_edge is introduced for convinience, although all functionality
//   can be achieved using a spectrum mask only   
void SDLineFinder::setScan(const SDMemTableWrapper &in_scan,
               const std::vector<bool> &in_mask,
	       const boost::python::tuple &in_edge) throw(AipsError)
{
  try {
       scan=in_scan.getCP();
       AlwaysAssert(!scan.null(),AipsError);
       if (scan->nRow()!=1)
           throw AipsError("SDLineFinder::setScan - in_scan contains more than 1 row."
	                   "Choose one first.");       
       mask=in_mask;
       if (mask.nelements()!=scan->nChan())
           throw AipsError("SDLineFinder::setScan - in_scan and in_mask have different"
	                   "number of spectral channels.");

       // number of elements in the in_edge tuple
       int n=extract<int>(in_edge.attr("__len__")());
       if (n>2 || n<0)
           throw AipsError("SDLineFinder::setScan - the length of the in_edge parameter"
	                   "should not exceed 2");
       if (!n) {
           // all spectrum, no rejection
           edge.first=0;
	   edge.second=scan->nChan();
       } else {
           edge.first=extract<int>(in_edge.attr("__getitem__")(0));
	   if (edge.first<0)
	       throw AipsError("SDLineFinder::setScan - the in_edge parameter has a negative"
	                        "number of channels to drop");
           if (edge.first>=scan->nChan())
	       throw AipsError("SDLineFinder::setScan - all channels are rejected by the in_edge parameter");
           if (n==2) {
	       edge.second=extract<int>(in_edge.attr("__getitem__")(1));
	       if (edge.second<0)
  	           throw AipsError("SDLineFinder::setScan - the in_edge parameter has a negative"
	                           "number of channels to drop");
               edge.second=scan->nChan()-edge.second;
	   } else edge.second=scan->nChan()-edge.first;
           if (edge.second<0 || (edge.second+edge.first)>scan->nChan())
	       throw AipsError("SDLineFinder::setScan - all channels are rejected by the in_edge parameter");
       }       
  }
  catch(const AipsError &ae) {
       // setScan is unsuccessfull, reset scan/mask/edge
       scan=CountedConstPtr<SDMemTable>(); // null pointer
       mask.resize(0);
       edge=pair<int,int>(0,0);
       throw;
  }
}

// search for spectral lines. Number of lines found is returned
int SDLineFinder::findLines() throw(casa::AipsError)
{
  const int minboxnchan=4;
  if (scan.null())
      throw AipsError("SDLineFinder::findLines - a scan should be set first,"
                      " use set_scan");
  DebugAssert(mask.nelements()==scan->nChan(), AipsError);
  int max_box_nchan=int(scan->nChan()*box_size); // number of channels in running
                                                 // box
  if (max_box_nchan<2)
      throw AipsError("SDLineFinder::findLines - box_size is too small");

  scan->getSpectrum(spectrum);

  lines.resize(0); // search from the scratch
  Vector<Bool> temp_mask(mask);
  Bool need2iterate=True;
  while (need2iterate) {
     // line find algorithm
     LFRunningMean lfalg(spectrum,temp_mask,edge,max_box_nchan,min_nchan,
                         threshold);
     std::list<pair<int,int> > new_lines;
     lfalg.findLines(new_lines);
     if (!new_lines.size()) need2iterate=False;
     temp_mask=getMask();
     // update the list (lines) merging intervals, if necessary
     addNewSearchResult(new_lines);
  }
  return int(lines.size());
}


// get the mask to mask out all lines that have been found (default)
// if invert=true, only channels belong to lines will be unmasked
// Note: all channels originally masked by the input mask (in_mask
//       in setScan) or dropped out by the edge parameter (in_edge
//       in setScan) are still excluded regardless on the invert option
std::vector<bool> SDLineFinder::getMask(bool invert)
                                        const throw(casa::AipsError)
{
  try {
       if (scan.null())
           throw AipsError("SDLineFinder::getMask - a scan should be set first,"
                      " use set_scan followed by find_lines");
       DebugAssert(mask.nelements()==scan->nChan(), AipsError);
       /*
       if (!lines.size())
           throw AipsError("SDLineFinder::getMask - one have to search for "
	                   "lines first, use find_lines");
       */			   
       std::vector<bool> res_mask(mask.nelements());
       // iterator through lines
       std::list<std::pair<int,int> >::const_iterator cli=lines.begin();
       for (int ch=0;ch<res_mask.size();++ch) 
            if (ch<edge.first || ch>=edge.second) res_mask[ch]=false;
	    else if (!mask[ch]) res_mask[ch]=false;
	    else {            
	            res_mask[ch]=!invert; // no line by default
		    if (cli==lines.end()) continue;
	            if (ch>=cli->first && ch<cli->second)
		        res_mask[ch]=invert; // this is a line
                    if (ch>=cli->second)
		        ++cli; // next line in the list
	         }
       
       return res_mask;
  }
  catch (const AipsError &ae) {
      throw;
  }  
  catch (const exception &ex) {
      throw AipsError(String("SDLineFinder::getMask - STL error: ")+ex.what());
  }
}

// get range for all lines found. If defunits is true (default), the
// same units as used in the scan will be returned (e.g. velocity
// instead of channels). If defunits is false, channels will be returned
std::vector<int> SDLineFinder::getLineRanges(bool defunits)
                             const throw(casa::AipsError)
{
  try {
       if (scan.null())
           throw AipsError("SDLineFinder::getLineRanges - a scan should be set first,"
                      " use set_scan followed by find_lines");
       DebugAssert(mask.nelements()==scan->nChan(), AipsError);
       
       if (!lines.size())
           throw AipsError("SDLineFinder::getLineRanges - one have to search for "
	                   "lines first, use find_lines");
       			   
       // temporary
       if (defunits)
           throw AipsError("SDLineFinder::getLineRanges - sorry, defunits=true have not "
	                   "yet been implemented");
       //
       std::vector<int> res(2*lines.size());
       // iterator through lines & result
       std::list<std::pair<int,int> >::const_iterator cli=lines.begin();
       std::vector<int>::iterator ri=res.begin();
       for (;cli!=lines.end() && ri!=res.end();++cli,++ri) {            
	    *ri=cli->first;
	    if (++ri!=res.end()) 
	        *ri=cli->second-1;	    
       }
       return res;
  }
  catch (const AipsError &ae) {
      throw;
  }  
  catch (const exception &ex) {
      throw AipsError(String("SDLineFinder::getLineRanges - STL error: ")+ex.what());
  }
}

// concatenate two lists preserving the order. If two lines appear to
// be adjacent, they are joined into the new one
void SDLineFinder::addNewSearchResult(const std::list<pair<int, int> > &newlines)
                        throw(AipsError)
{
  try {
      for (std::list<pair<int,int> >::const_iterator cli=newlines.begin();
           cli!=newlines.end();++cli) {
	   
	   // the first item, which has a non-void intersection or touches
	   // the new line
	   std::list<pair<int,int> >::iterator pos_beg=find_if(lines.begin(),
	                  lines.end(), IntersectsWith(*cli));           
	   // the last such item	  
	   std::list<pair<int,int> >::iterator pos_end=find_if(pos_beg,
	                  lines.end(), not1(IntersectsWith(*cli)));

           // extract all lines which intersect or touch a new one into
	   // a temporary buffer. This may invalidate the iterators
	   // line_buffer may be empty, if no lines intersects with a new
	   // one.
	   std::list<pair<int,int> > lines_buffer;
	   lines_buffer.splice(lines_buffer.end(),lines, pos_beg, pos_end);

	   // build a union of all intersecting lines 
	   pair<int,int> union_line=for_each(lines_buffer.begin(),
	           lines_buffer.end(),BuildUnion(*cli)).result();
           
	   // search for a right place for the new line (union_line) and add
	   std::list<pair<int,int> >::iterator pos2insert=find_if(lines.begin(),
	                  lines.end(), LaterThan(union_line));
	   lines.insert(pos2insert,union_line);
      }
  }
  catch (const AipsError &ae) {
      throw;
  }  
  catch (const exception &ex) {
      throw AipsError(String("SDLineFinder::addNewSearchResult - STL error: ")+ex.what());
  }
}
