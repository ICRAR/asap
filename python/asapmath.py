from scantable import scantable
def average_scan(scan):
    """
        Return a (time) averaged a scan, i.e. all correlator cycles
        are averaged into one "scan".
    """
    from asap._asap import average as _av
    return scantable(_av(scan))


def average_scans(*args, **kwargs):
    """
        Return the (time) average of a list of scans. [in channels only]
        Parameters:
            a comma separated list of scans
            mask:     an optional mask
        Example:
            scanav = average_scans(scana,scanb)
            # return a time averaged scan from scan and scanb
            # without using a mask
    """

    if len(args) < 2:
        print "Please give at least two scantables"
        return

    from asap._asap import averages as _av
    d = [args[0].nbeam(),args[0].nif(),args[0].npol(),args[0].nchan()]
    for s in args:
        if not isinstance(s,scantable):
            print "Please give a list of scantables"
            return
        dim = [s.nbeam(),s.nif(),s.npol(),s.nchan()]
        if (dim != d):
            print "All scans have to have the same numer of Beams/IFs/Pols/Chans"
            return
    if kwargs.has_key('mask'):
        return scantable(_av(args, kwargs.get('mask')))
    else:
        from numarray import ones
        mask = list(ones(d[3]))
        return scantable(_av((args), mask))

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
        factor:      the sclaing factor
    Note:
        This currently applies the all beams/IFs/pols
    """
    from asap._asap import scale as _scale
    return scantable(_scale(scan, factor))


def bin(scan, binwidth=5):
    """
    """
    from asap._asap import bin as _bin
    return scantable(_bin(scan, binwidth))
