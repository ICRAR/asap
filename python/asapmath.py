from scantable import scantable

def average_time(*args, **kwargs):
    """
    Return the (time) average of a scan or list of scans. [in channels only]
    Parameters:
        one scan or comma separated  scans
        mask:     an optional mask
    Example:
        # return a time averaged scan from scana and scanb
        # without using a mask
        scanav = average_scans(scana,scanb)
        # return the (time) averaged scan, i.e. the average of
        # all correlator cycles
        scanav = average_time(scan)
        
    """
    lst = args
    if len(args) < 2:
        if type(args[0]) is list:
            if  len(args[0]) < 2:
                print "Please give at least two scantables"
                return
        else:
            s = args[0]
            if s.nrow() > 1:
                from asap._asap import average as _av
                return scantable(_av(s))
            else:
                print "Given scantable is already time averaged"
                return
        lst = tuple(args[0])
    else:
        lst = tuple(args)
    from asap._asap import averages as _avs
    d = [lst[0].nbeam(),lst[0].nif(),lst[0].npol(),lst[0].nchan()]
    for s in lst:
        if not isinstance(s,scantable):
            print "Please give a list of scantables"
            return
        dim = [s.nbeam(),s.nif(),s.npol(),s.nchan()]
        if (dim != d):
            print "All scans have to have the same numer of Beams/IFs/Pols/Chans"
            return
    if kwargs.has_key('mask'):
        return scantable(_avs(lst, kwargs.get('mask')))
    else:
        from numarray import ones
        mask = tuple(ones(d[3]))
        return scantable(_avs(lst, mask))

def quotient(source, reference):
    """
    Return the quotient of a 'source' scan and a 'reference' scan
    Parameters:
        source:        the 'on' scan
        reference:     the 'off' scan
    """
    from asap._asap import quotient as _quot
    return scantable(_quot(source, reference))

def scale(scan, factor):
    """
    Return a scan where all spectra are scaled by the give 'factor'
    Parameters:
        scan:        a scantable
        factor:      the scaling factor
    Note:
        This currently applies the all beams/IFs/pols
    """
    from asap._asap import scale as _scale
    return scantable(_scale(scan, factor))

def add(scan, offset):
    """
    Return a scan where the offset is added.
    Parameters:
        scan:        a scantable
        offset:      the value to add
    Note:
        This currently applies the all beams/IFs/pols
    """
    from asap._asap import add as _add
    return scantable(_add(scan, offset))


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
