from scantable import scantable

def average_time(*args, **kwargs):
    """
    Return the (time) average of a scan or list of scans. [in channels only]
    Parameters:
        one scan or comma separated  scans
        mask:     an optional mask (only used for 'var' and 'tsys' weighting)
        scanav:   False (default) averages all scans together,
                  True averages each scan separately
        weight:   Weighting scheme. 'none' (default), 'var' (variance
                  weighted), 'tsys'
    Example:
        # return a time averaged scan from scana and scanb
        # without using a mask
        scanav = average_time(scana,scanb)
        # return the (time) averaged scan, i.e. the average of
        # all correlator cycles
        scanav = average_time(scan)

    """
    scanAv = False
    if kwargs.has_key('scanav'):
       scanAv = kwargs.get('scanav')
#
    weight = 'none'
    if kwargs.has_key('weight'):
       weight = kwargs.get('weight')
#
    mask = ()
    if kwargs.has_key('mask'):
        mask = kwargs.get('mask')
#
    lst = tuple(args)
    from asap._asap import average as _av
    for s in lst:
        if not isinstance(s,scantable):
            print "Please give a list of scantables"
            return
    return scantable(_av(lst, mask, scanAv, weight))

def quotient(source, reference):
    """
    Return the quotient of a 'source' scan and a 'reference' scan
    Parameters:
        source:        the 'on' scan
        reference:     the 'off' scan
    """
    from asap._asap import quotient as _quot
    return scantable(_quot(source, reference))

def scale(scan, factor, insitu=False, all=True):
    """
    Return a scan where all spectra are scaled by the give 'factor'
    Parameters:
        scan:        a scantable
        factor:      the scaling factor
        insitu:      if False (default) a new scantable is returned.
                     Otherwise, the scaling is done in-situ
        all:         if True (default) apply to all spectra. Otherwise
                     apply only to the selected (beam/pol/if)spectra only
    """
    if not insitu:
        from asap._asap import scale as _scale
        return scantable(_scale(scan, factor, all))
    else:
        from asap._asap import scale_insitu as _scale
        _scale(scan, factor, all)
        return
        

def add(scan, offset, insitu=False, all=True):
    """
    Return a scan where all spectra have the offset added
    Parameters:
        scan:        a scantable
        offset:      the offset
        insitu:      if False (default) a new scantable is returned.
                     Otherwise, the addition is done in-situ
        all:         if True (default) apply to all spectra. Otherwise
                     apply only to the selected (beam/pol/if)spectra only
    """
    if not insitu:
        from asap._asap import add as _add
        return scantable(_add(scan, offset, all))
    else:
        from asap._asap import add_insitu as _add
        _add(scan, offset, all)
        return
        
def bin(scan, binwidth=5):
    """
    """
    from asap._asap import bin as _bin
    return scantable(_bin(scan, binwidth))

def average_pol(scan, mask=None):
    """
    Average the Polarisations together.
    Parameters:
        scan   - a scantable
        mask   - an optional mask defining the region, where
                 the averaging will be applied. The output
                 will have all specified points masked.
                 The dimension won't be reduced and
                 all polarisations will contain the
                 averaged spectrum.
    Example:
        polav = average_pols(myscan)
    """
    from asap._asap import averagepol as _avpol
    from numarray import ones
    if mask is None:
        mask = tuple(ones(scan.nchan()))
    return scantable(_avpol(scan, mask))
    
def hanning(scan):
    """
    Hanning smooth the channels.
    Parameters:
         scan    - the input scan
    Example:
         none
    """
    from asap._asap import hanning as _han
    return scantable(_han(scan))

    
def poly_baseline(scan, mask=None, order=0):
    """
    Return a scan which has been baselined by a polynomial.
    Parameters:
        scan:    a scantable
        mask:    an optional mask
        order:   the order of the polynomial (default is 0)
    Example:
        # return a scan baselined by a third order polynomial,
        # not using a mask
        bscan = poly_baseline(scan, order=3)
    """
    from asap.asapfitter import fitter
    if mask is None:
        from numarray import ones
        mask = tuple(ones(scan.nchan()))
    f = fitter()
    f._verbose(True)
    f.set_scan(scan, mask)
    f.set_function(poly=order)    
    sf = f.auto_fit()
    return sf
