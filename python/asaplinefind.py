import _asap

class linefinder:
    """
    The class for automated spectral line search in ASAP.
    """

    def __init__(self):
	"""
	Create a line finder object.
	"""
	self.finder = _asap.linefinder()
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
    def find_lines(self):
	"""
	Search for spectral lines in the scan assigned in set_scan.
	A number of lines found will be returned
	"""
	return self.finder.findlines()
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
	return self.finder.getlineranges(defunits)
