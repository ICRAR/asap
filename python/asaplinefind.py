import _asap

class linefinder:
    """
    The class for automated spectral line search in ASAP.

    Example:
       fl=linefinder()
       fl.set_scan(sc,edge=(50,))
       fl.set_options(threshold=3)
       nlines=fl.find_lines()
       if nlines!=0:
          print "Found ",nlines," spectral lines"
          print fl.get_ranges(False)
       else:
          print "No lines found!"
       sc2=poly_baseline(sc,fl.get_mask(),7)
    
    The algorithm involves a simple threshold criterion. The line is
    considered to be detected if a specified number of consequtive
    channels (default is 3) is brighter (with respect to the current baseline
    estimate) than the threshold times the noise level. This criterion is
    applied in the iterative procedure updating baseline estimate and trying
    reduced spectral resolutions to detect broad lines as well. The off-line
    noise level is determined at each iteration as an average of 80% of the
    lowest variances across the spectrum (i.e. histogram equalization is
    used to avoid missing weak lines if strong ones are present). For
    bad baseline shapes it is reccommended to increase the threshold and 
    possibly switch the averaging option off (see set_options) to
    detect strong lines only, fit a high order baseline and repeat the line
    search. 

    """

    def __init__(self):
	"""
	Create a line finder object.
	"""
	self.finder = _asap.linefinder()
	return

    def set_options(self,threshold=1.7320508075688772,min_nchan=3,
        avg_limit=8,box_size=0.2):
	"""
	Set the parameters of the algorithm
	Parameters:
	     threshold    a single channel S/N ratio above which the 
	                  channel is considered to be a detection
			  Default is sqrt(3), which together with
			  min_nchan=3 gives a 3-sigma criterion
	     min_nchan    a minimal number of consequtive channels, 
                          which should satisfy a threshold criterion to
			  be a detection. Default is 3.
	     avg_limit    A number of consequtive channels not greater than
	                  this parameter can be averaged to search for
			  broad lines. Default is 8.
	     box_size     A running mean box size specified as a fraction
                          of the total spectrum length. Default is 1/5
	Note:  For bad baselines threshold should be increased, 
	       and avg_limit decreased (or even switched off completely by
	       setting this parameter to 1) to avoid detecting baseline
	       undulations instead of real lines.  
        """
        self.finder.setoptions(threshold,min_nchan,avg_limit,box_size)
	return
	     
    def set_scan(self,scan,mask=None,edge=(0,0)):
	"""
	Set the 'data' (scantable) to work with.
	Parameters:
	     scan:    a scantable
	     mask:       an optional mask retreived from scantable
	     edge:       an optional number of channel to drop at
			 the edge of spectrum. If only one value is
			 specified, the same number will be dropped from
			 both sides of the spectrum. Default is to keep
			 all channels
        """
	if not scan:
	   raise RuntimeError, 'Please give a correct scan'
	if len(edge)>2:
	   raise RuntimeError, "The edge parameter should have two \
           or less elements"
	if mask is None:
	    from numarray import ones
	    self.finder.setscan(scan,ones(scan.nchan()),edge)
	else:    
	    self.finder.setscan(scan,mask,edge)
	return 
    def find_lines(self,nRow=0):
	"""
	Search for spectral lines in the scan assigned in set_scan.
	Current Beam/IF/Pol is used, Row is specified by parameter
	A number of lines found will be returned
	"""
	return self.finder.findlines(nRow)
    def get_mask(self,invert=False):
	"""
	Get the mask to mask out all lines that have been found (default)

	Parameters:
	      invert  if True, only channels belong to lines will be unmasked

	Note: all channels originally masked by the input mask or
	      dropped out by the edge parameter will still be excluded
	      regardless on the invert option
        """
	return self.finder.getmask(invert)
    def get_ranges(self,defunits=True):
	"""
	Get ranges (start and end channels or velocities) for all spectral
	lines found. 

	Parameters:
	      defunits  if True (default), the range will use the same units
		        as set for the scan (e.g. LSR velocity)
			if False, the range will be expressed in channels 
	"""
	if (defunits):
	    return self.finder.getlineranges()
	else:
	    return self.finder.getlinerangesinchannels()

def auto_poly_baseline(scan, mask=None, edge=(0,0), order=0,
    threshold=3,insitu=None):
    """
    Return a scan which has been baselined (all rows) by a polynomial.
    Spectral lines are detected first using linefinder and masked out
    to avoid them affecting the baseline solution.

    Parameters:
        scan:    a scantable
        mask:       an optional mask retreived from scantable
        edge:       an optional number of channel to drop at
                    the edge of spectrum. If only one value is
                    specified, the same number will be dropped from
                    both sides of the spectrum. Default is to keep
                    all channels
        order:      the order of the polynomial (default is 0)
	threshold:  the threshold used by line finder. It is better to
                    keep it large as only strong lines affect the
		    baseline solution.
        insitu:     if False a new scantable is returned.
                    Otherwise, the scaling is done in-situ
                    The default is taken from .asaprc (False)

    Example:
        sc2=auto_poly_baseline(sc,order=7)
    """
    from asap.asapfitter import fitter
    from asap import scantable

    # setup fitter

    f = fitter()
    f._verbose(True)
    f.set_function(poly=order)

    # setup line finder

    fl=linefinder()
    fl.set_options(threshold=threshold)

    if not insitu:
        workscan=scan.copy()
    else:
        workscan=scan

    vb=workscan._vb
    # remember the verbose parameter and selection
    workscan._vb=False
    sel=workscan.get_cursor()
    rows=range(workscan.nrow()) 
    for i in range(workscan.nbeam()):
        workscan.setbeam(i)
	for j in range(workscan.nif()):
	    workscan.setif(j)
	    for k in range(workscan.npol()):
	        scan.setpol(k)
		if f._vb:
		   print "Processing:"
		   print 'Beam[%d], IF[%d], Pol[%d]' % (i,j,k)
		for iRow in rows:
                   fl.set_scan(workscan,mask,edge)
                   fl.find_lines(iRow)
                   f.set_scan(workscan, fl.get_mask())
		   f.x=workscan._getabcissa(iRow)
		   f.y=workscan._getspectrum(iRow)
                   f.data=None
		   f.fit()
		   workscan._setspectrum(f.getresidual(),iRow)
    workscan.set_cursor(sel[0],sel[1],sel[2])
    workscan._vb = vb
    if not insitu:
       return scan
