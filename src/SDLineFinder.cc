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
#include <fstream>

using namespace asap;
using namespace casa;
using namespace std;
using namespace boost::python;

namespace asap {

///////////////////////////////////////////////////////////////////////////////
//
// IStatHolder - an abstract class to collect statistics from the running
//               mean calculator, if necessary.
//               We define it here, because it is used in LFRunningMean and
//               SDLineFinder only
//

struct IStatHolder  {
   // This function is called for each spectral channel processed by
   // the running mean calculator. The order of channel numbers may be
   // arbitrary
   // ch - a number of channel, corresponding to (approximately) the centre
   //      of the running box
   // box_nchan - number of channels in the box
   //
   virtual void accumulate(int ch, Float sum, Float sum2, int box_nchan)
                      throw(AipsError) = 0;		      
};

//
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//
// SignAccumulator - a simple class to deal with running mean statistics:
//                   it stores the sign of the value-mean only
//
class SignAccumulator : public IStatHolder {
   Vector<Int>  sign;                // either +1, -1 or 0
   const Vector<Float>   &spectrum;  // a reference to the spectrum
                                     // to calculate the sign   
public:
   // all channels >=nchan are ignored
   SignAccumulator(uInt nchan, const Vector<Float> &in_spectrum);
   
   // This function is called for each spectral channel processed by
   // the running mean calculator. The order of channel numbers may be
   // arbitrary
   // ch - a number of channel, corresponding to (approximately) the centre
   //      of the running box
   // box_nchan - number of channels in the box
   //
   virtual void accumulate(int ch, Float sum, Float, int box_nchan)
                      throw(AipsError);

   // access to the sign
   const Vector<Int>& getSigns() const throw();
};

//
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//
// RunningBox -    a running box calculator. This class implements 
//                 interations over the specified spectrum and calculates
//                 running box filter statistics.
//

class RunningBox {
   // The input data to work with. Use reference symantics to avoid
   // an unnecessary copying   
   const casa::Vector<casa::Float>  &spectrum; // a buffer for the spectrum
   const casa::Vector<casa::Bool>   &mask; // associated mask
   const std::pair<int,int>         &edge; // start and stop+1 channels
                                           // to work with
   
   // statistics for running box filtering
   casa::Float sumf;       // sum of fluxes
   casa::Float sumf2;     // sum of squares of fluxes
   casa::Float sumch;       // sum of channel numbers (for linear fit)
   casa::Float sumch2;     // sum of squares of channel numbers (for linear fit)
   casa::Float sumfch;     // sum of flux*(channel number) (for linear fit)
   
   int box_chan_cntr;     // actual number of channels in the box
   int max_box_nchan;     // maximum allowed number of channels in the box
                          // (calculated from boxsize and actual spectrum size)
   // cache for derivative statistics
   mutable casa::Bool need2recalculate; // if true, values of the statistics
                                       // below are invalid
   mutable casa::Float linmean;  // a value of the linear fit to the
                                 // points in the running box
   mutable casa::Float linvariance; // the same for variance
   int cur_channel;       // the number of the current channel
   int start_advance;     // number of channel from which the box can
                          // be moved (the middle of the box, if there is no
			  // masking)
public:
   // set up the object with the references to actual data
   // as well as the number of channels in the running box
   RunningBox(const casa::Vector<casa::Float>  &in_spectrum,
                 const casa::Vector<casa::Bool>   &in_mask,
		 const std::pair<int,int>         &in_edge,
		 int in_max_box_nchan) throw(AipsError);
		 
   // access to the statistics
   const casa::Float& getLinMean() const throw(AipsError);
   const casa::Float& getLinVariance() const throw(AipsError);
   const casa::Float aboveMean() const throw(AipsError);
   int getChannel() const throw();
   
   // actual number of channels in the box (max_box_nchan, if no channels
   // are masked)
   int getNumberOfBoxPoints() const throw();

   // next channel
   void next() throw(AipsError);

   // checking whether there are still elements
   casa::Bool haveMore() const throw();

   // go to start
   void rewind() throw(AipsError);
 
protected:
   // supplementary function to control running mean calculations.
   // It adds a specified channel to the running mean box and
   // removes (ch-maxboxnchan+1)'th channel from there
   // Channels, for which the mask is false or index is beyond the
   // allowed range, are ignored
   void advanceRunningBox(int ch) throw(casa::AipsError);

   // calculate derivative statistics. This function is const, because
   // it updates the cache only
   void updateDerivativeStatistics() const throw(AipsError);
};

//
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//
// LFAboveThreshold   An algorithm for line detection using running box
//                    statistics.  Line is detected if it is above the
//                    specified threshold at the specified number of
//                    consequtive channels. Prefix LF stands for Line Finder
//
class LFAboveThreshold {
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
   std::list<pair<int,int> > &lines;       // list where detections are saved
                                           // (pair: start and stop+1 channel)
   RunningBox *running_box;                // running box filter
public:

   // set up the detection criterion
   LFAboveThreshold(std::list<pair<int,int> > &in_lines,
                    int in_min_nchan = 3,
		    casa::Float in_threshold = 5) throw();
   virtual ~LFAboveThreshold() throw();
   
   // replace the detection criterion
   void setCriterion(int in_min_nchan, casa::Float in_threshold) throw();

   // find spectral lines and add them into list
   // if statholder is not NULL, the accumulate function of it will be
   // called for each channel to save statistics
   //    spectrum, mask and edge - reference to the data
   //    max_box_nchan  - number of channels in the running box
   void findLines(const casa::Vector<casa::Float> &spectrum,
		  const casa::Vector<casa::Bool> &mask,
		  const std::pair<int,int> &edge,
		  int max_box_nchan,
                  IStatHolder* statholder = NULL) throw(casa::AipsError);

protected:

   // process a channel: update curline and is_detected before and
   // add a new line to the list, if necessary using processCurLine()
   // detect=true indicates that the current channel satisfies the criterion
   void processChannel(Bool detect, const casa::Vector<casa::Bool> &mask)
                                                  throw(casa::AipsError);

   // process the interval of channels stored in curline
   // if it satisfies the criterion, add this interval as a new line
   void processCurLine(const casa::Vector<casa::Bool> &mask)
                                                 throw(casa::AipsError);
};

//
///////////////////////////////////////////////////////////////////////////////

} // namespace asap

///////////////////////////////////////////////////////////////////////////////
//
// SignAccumulator - a simple class to deal with running mean statistics:
//                   it stores the sign of the value-mean only
//

// all channels >=nchan are ignored
SignAccumulator::SignAccumulator(uInt nchan,
                   const Vector<Float> &in_spectrum) : sign(nchan,0),
		   spectrum(in_spectrum) {}


// This function is called for each spectral channel processed by
// the running mean calculator. The order of channel numbers may be
// arbitrary
// ch - a number of channel, corresponding to (approximately) the centre
//      of the running box
// box_nchan - number of channels in the box
//
void SignAccumulator::accumulate(int ch, Float sum, Float sum2, int box_nchan)
                   throw(AipsError)
{
   if (ch>=sign.nelements()) return;
   DebugAssert(ch>=0,AipsError);
   DebugAssert(ch<=spectrum.nelements(), AipsError);
   if (box_nchan) {
       Float buf=spectrum[ch]-sum; ///Float(box_nchan);
       if (buf>0) sign[ch]=1;
       else if (buf<0) sign[ch]=-1;
       else sign[ch]=0;
   } else sign[ch]=0;    
}

// access to the sign
const Vector<Int>& SignAccumulator::getSigns() const throw()
{
   return sign;
}


//
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//
// RunningBox -    a running box calculator. This class implements 
//                 interations over the specified spectrum and calculates
//                 running box filter statistics.
//

// set up the object with the references to actual data
// and the number of channels in the running box
RunningBox::RunningBox(const casa::Vector<casa::Float>  &in_spectrum,
                       const casa::Vector<casa::Bool>   &in_mask,
     	               const std::pair<int,int>         &in_edge,
		       int in_max_box_nchan) throw(AipsError) :
        spectrum(in_spectrum), mask(in_mask), edge(in_edge),
	max_box_nchan(in_max_box_nchan)
{
  rewind();
}

void RunningBox::rewind() throw(AipsError) {
  // fill statistics for initial box
  box_chan_cntr=0; // no channels are currently in the box
  sumf=0.;           // initialize statistics
  sumf2=0.;
  sumch=0.;
  sumch2=0.;
  sumfch=0.;
  int initial_box_ch=edge.first;
  for (;initial_box_ch<edge.second && box_chan_cntr<max_box_nchan;
        ++initial_box_ch)
       advanceRunningBox(initial_box_ch);
  
  if (initial_box_ch==edge.second)       
      throw AipsError("RunningBox::rewind - too much channels are masked");

  cur_channel=edge.first;
  start_advance=initial_box_ch-max_box_nchan/2;  
}

// access to the statistics
const casa::Float& RunningBox::getLinMean() const throw(AipsError)
{
  DebugAssert(cur_channel<edge.second, AipsError);
  if (need2recalculate) updateDerivativeStatistics();
  return linmean;
}

const casa::Float& RunningBox::getLinVariance() const throw(AipsError)
{
  DebugAssert(cur_channel<edge.second, AipsError);
  if (need2recalculate) updateDerivativeStatistics();
  return linvariance;
}

const casa::Float RunningBox::aboveMean() const throw(AipsError)
{
  DebugAssert(cur_channel<edge.second, AipsError);
  if (need2recalculate) updateDerivativeStatistics();
  return spectrum[cur_channel]-linmean;
}

int RunningBox::getChannel() const throw()
{
  return cur_channel;
}

// actual number of channels in the box (max_box_nchan, if no channels
// are masked)
int RunningBox::getNumberOfBoxPoints() const throw()
{
  return box_chan_cntr;
}

// supplementary function to control running mean calculations.
// It adds a specified channel to the running mean box and
// removes (ch-max_box_nchan+1)'th channel from there
// Channels, for which the mask is false or index is beyond the
// allowed range, are ignored
void RunningBox::advanceRunningBox(int ch) throw(AipsError)
{
  if (ch>=edge.first && ch<edge.second)
      if (mask[ch]) { // ch is a valid channel
          ++box_chan_cntr;
          sumf+=spectrum[ch];
          sumf2+=square(spectrum[ch]);
	  sumch+=Float(ch);
	  sumch2+=square(Float(ch));
	  sumfch+=spectrum[ch]*Float(ch);
	  need2recalculate=True;
      }
  int ch2remove=ch-max_box_nchan;
  if (ch2remove>=edge.first && ch2remove<edge.second)
      if (mask[ch2remove]) { // ch2remove is a valid channel
          --box_chan_cntr;
          sumf-=spectrum[ch2remove];
          sumf2-=square(spectrum[ch2remove]);  
	  sumch-=Float(ch2remove);
	  sumch2-=square(Float(ch2remove));
	  sumfch-=spectrum[ch2remove]*Float(ch2remove);
	  need2recalculate=True;
      }
}

// next channel
void RunningBox::next() throw(AipsError)
{
   AlwaysAssert(cur_channel<edge.second,AipsError);
   ++cur_channel;
   if (cur_channel+max_box_nchan/2<edge.second && cur_channel>=start_advance)
       advanceRunningBox(cur_channel+max_box_nchan/2); // update statistics
}

// checking whether there are still elements
casa::Bool RunningBox::haveMore() const throw()
{
   return cur_channel<edge.second;
}

// calculate derivative statistics. This function is const, because
// it updates the cache only
void RunningBox::updateDerivativeStatistics() const throw(AipsError)
{
  AlwaysAssert(box_chan_cntr, AipsError);
  
  Float mean=sumf/Float(box_chan_cntr);

  // linear LSF formulae
  Float meanch=sumch/Float(box_chan_cntr);
  Float meanch2=sumch2/Float(box_chan_cntr);
  if (meanch==meanch2 || box_chan_cntr<3) {
      // vertical line in the spectrum, can't calculate linmean and linvariance
      linmean=0.;
      linvariance=0.;
  } else {
      Float coeff=(sumfch/Float(box_chan_cntr)-meanch*mean)/
                (meanch2-square(meanch));
      linmean=coeff*(Float(cur_channel)-meanch)+mean;
      linvariance=sqrt(sumf2/Float(box_chan_cntr)-square(mean)-
                    square(coeff)*(meanch2-square(meanch)));
  }
  need2recalculate=False;
}


//
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//
// LFAboveThreshold - a running mean algorithm for line detection
//
//


// set up the detection criterion
LFAboveThreshold::LFAboveThreshold(std::list<pair<int,int> > &in_lines,
                                   int in_min_nchan,
                                   casa::Float in_threshold) throw() :
             min_nchan(in_min_nchan), threshold(in_threshold),
	     lines(in_lines), running_box(NULL) {}

LFAboveThreshold::~LFAboveThreshold() throw()
{
  if (running_box!=NULL) delete running_box;
}

// replace the detection criterion
void LFAboveThreshold::setCriterion(int in_min_nchan, casa::Float in_threshold)
                                 throw()
{
  min_nchan=in_min_nchan;
  threshold=in_threshold;
}


// process a channel: update cur_line and is_detected before and
// add a new line to the list, if necessary
void LFAboveThreshold::processChannel(Bool detect,
                 const casa::Vector<casa::Bool> &mask) throw(casa::AipsError)
{
  try {
       if (detect) {
           if (is_detected_before)
               cur_line.second=running_box->getChannel()+1;
	   else {
	       is_detected_before=True;
	       cur_line.first=running_box->getChannel();
	       cur_line.second=running_box->getChannel()+1;
	   }
       } else processCurLine(mask);   
  }
  catch (const AipsError &ae) {
      throw;
  }  
  catch (const exception &ex) {
      throw AipsError(String("LFAboveThreshold::processChannel - STL error: ")+ex.what());
  }
}

// process the interval of channels stored in cur_line
// if it satisfies the criterion, add this interval as a new line
void LFAboveThreshold::processCurLine(const casa::Vector<casa::Bool> &mask)
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
      throw AipsError(String("LFAboveThreshold::processCurLine - STL error: ")+ex.what());
  }
}

// find spectral lines and add them into list
void LFAboveThreshold::findLines(const casa::Vector<casa::Float> &spectrum,
		              const casa::Vector<casa::Bool> &mask,
		              const std::pair<int,int> &edge,
		              int max_box_nchan,
                              IStatHolder* statholder)
                        throw(casa::AipsError)
{
  const int minboxnchan=4;
  try {

      if (running_box!=NULL) delete running_box;
      running_box=new RunningBox(spectrum,mask,edge,max_box_nchan);
    
      // actual search algorithm
      is_detected_before=False;
    
      for (;running_box->haveMore();running_box->next()) {
           const int ch=running_box->getChannel();
           if (running_box->getNumberOfBoxPoints()>=minboxnchan)
	       processChannel(mask[ch] && (fabs(spectrum[ch]-
	          running_box->getLinMean()) >=
		  threshold*running_box->getLinVariance()), mask);
	   else processCurLine(mask); // just finish what was accumulated before
	   if (statholder!=NULL)
              statholder->accumulate(running_box->getChannel(),
	                             running_box->getLinMean(),
				     running_box->getLinVariance(),
				     running_box->getNumberOfBoxPoints());
      }
  }
  catch (const AipsError &ae) {
      throw;
  }  
  catch (const exception &ex) {
      throw AipsError(String("LFAboveThreshold::findLines - STL error: ")+ex.what());
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
  box_size=1./5.; // default box size for running mean calculations is
                  // 1/5 of the whole spectrum
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

  Bool first_pass=True;
  while (true) {
     // a buffer for new lines found at this iteration
     std::list<pair<int,int> > new_lines;

     // line find algorithm
     LFAboveThreshold lfalg(new_lines,min_nchan, threshold);

     SignAccumulator sacc(spectrum.nelements(),spectrum);


     try {
         lfalg.findLines(spectrum,temp_mask,edge,max_box_nchan,&sacc);
     }
     catch(const AipsError &ae) {
         if (first_pass) throw;
	 break; // nothing new
     }
     first_pass=False;
     if (!new_lines.size()) break; // nothing new

     searchForWings(new_lines, sacc.getSigns());

     // update the list (lines) merging intervals, if necessary
     addNewSearchResult(new_lines,lines);
     // get a new mask
     temp_mask=getMask();     
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
void SDLineFinder::addNewSearchResult(const std::list<pair<int, int> > &newlines,
                         std::list<std::pair<int, int> > &lines_list) 
                        throw(AipsError)
{
  try {
      for (std::list<pair<int,int> >::const_iterator cli=newlines.begin();
           cli!=newlines.end();++cli) {
	   
	   // the first item, which has a non-void intersection or touches
	   // the new line
	   std::list<pair<int,int> >::iterator pos_beg=find_if(lines_list.begin(),
	                  lines_list.end(), IntersectsWith(*cli));           
	   // the last such item	  
	   std::list<pair<int,int> >::iterator pos_end=find_if(pos_beg,
	                  lines_list.end(), not1(IntersectsWith(*cli)));

           // extract all lines which intersect or touch a new one into
	   // a temporary buffer. This may invalidate the iterators
	   // line_buffer may be empty, if no lines intersects with a new
	   // one.
	   std::list<pair<int,int> > lines_buffer;
	   lines_buffer.splice(lines_buffer.end(),lines_list, pos_beg, pos_end);

	   // build a union of all intersecting lines 
	   pair<int,int> union_line=for_each(lines_buffer.begin(),
	           lines_buffer.end(),BuildUnion(*cli)).result();
           
	   // search for a right place for the new line (union_line) and add
	   std::list<pair<int,int> >::iterator pos2insert=find_if(lines_list.begin(),
	                  lines_list.end(), LaterThan(union_line));
	   lines_list.insert(pos2insert,union_line);
      }
  }
  catch (const AipsError &ae) {
      throw;
  }  
  catch (const exception &ex) {
      throw AipsError(String("SDLineFinder::addNewSearchResult - STL error: ")+ex.what());
  }
}

// extend all line ranges to the point where a value stored in the
// specified vector changes (e.g. value-mean change its sign)
// This operation is necessary to include line wings, which are below
// the detection threshold. If lines becomes adjacent, they are
// merged together. Any masked channel stops the extension
void SDLineFinder::searchForWings(std::list<std::pair<int, int> > &newlines,
           const casa::Vector<casa::Int> &signs) throw(casa::AipsError)
{
  try {
      for (std::list<pair<int,int> >::iterator li=newlines.begin();
           li!=newlines.end();++li) {
	   // update the left hand side
	   for (int n=li->first-1;n>=edge.first;--n) {
	        if (!mask[n]) break;
	        if (signs[n]==signs[li->first] && signs[li->first])
		    li->first=n;
		else break;    
	   }
	   // update the right hand side
	   for (int n=li->second;n<edge.second;++n) {
	        if (!mask[n]) break;
		if (signs[n]==signs[li->second-1] && signs[li->second-1])
		    li->second=n;
		else break;    
	   }
      }
      // need to search for possible mergers.
      std::list<std::pair<int, int> >  result_buffer;
      addNewSearchResult(newlines,result_buffer);
      newlines.clear();
      newlines.splice(newlines.end(),result_buffer);
  }
  catch (const AipsError &ae) {
      throw;
  }  
  catch (const exception &ex) {
      throw AipsError(String("SDLineFinder::extendLines - STL error: ")+ex.what());
  }
}
