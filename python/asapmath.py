from scantable import scantable
from asap import rcParams

def average_time(*args, **kwargs):
    """
    Return the (time) average of a scan or list of scans. [in channels only]
    The cursor of the output scan is set to 0
    Parameters:
        one scan or comma separated  scans
        mask:     an optional mask (only used for 'var' and 'tsys' weighting)
        scanav:   True (default) averages each scan separately
                  False averages all scans together,
        weight:   Weighting scheme. 'none' (default), 'var' (1/var(spec) 
                  weighted), 'tsys' (1/Tsys**2 weighted), or 'tint'
                  (integration time weighted)
    Example:
        # return a time averaged scan from scana and scanb
        # without using a mask
        scanav = average_time(scana,scanb)
        # return the (time) averaged scan, i.e. the average of
        # all correlator cycles
        scanav = average_time(scan)

    """
    scanAv = True
    if kwargs.has_key('scanav'):
       scanAv = kwargs.get('scanav')
    weight = 'none'
    if kwargs.has_key('weight'):
       weight = kwargs.get('weight')
    mask = ()
    if kwargs.has_key('mask'):
        mask = kwargs.get('mask')
    varlist = vars()
    lst = tuple(args)
    del varlist["kwargs"]
    varlist["args"] = "%d scantables" % len(lst)
    # need special formatting her for history...
    
    from asap._asap import average as _av
    for s in lst:
        if not isinstance(s,scantable):
            print "Please give a list of scantables"
            return
    s = scantable(_av(lst, mask, scanAv, weight))
    s._add_history("average_time",varlist)
    return s

def quotient(source, reference, preserve=True):
    """
    Return the quotient of a 'source' (signal) scan and a 'reference' scan.
    The reference can have just one row, even if the signal has many. Otherwise
    they must have the same number of rows.
    The cursor of the output scan is set to 0
    Parameters:
        source:        the 'on' scan
        reference:     the 'off' scan
        preserve:      you can preserve (default) the continuum or 
                       remove it.  The equations used are 
                          preserve - Output = Tref * (sig/ref) - Tref
                          remove   - Output = Tref * (sig/ref) - Tsig
    """
    varlist = vars()
    from asap._asap import quotient as _quot
    s = scantable(_quot(source, reference, preserve))
    s._add_history("quotient",varlist)
    return s

def simple_math(left, right, op='add', tsys=True):
    """
    Apply simple mathematical binary operations to two 
    scan tables,  returning the result in a new scan table.
    The operation is applied to both the correlations and the TSys data
    The cursor of the output scan is set to 0
    Parameters:
        left:          the 'left' scan
        right:         the 'right' scan
        op:            the operation: 'add' (default), 'sub', 'mul', 'div'
        tsys:          if True (default) then apply the operation to Tsys
                       as well as the data
    """
    varlist = vars()
    if not isinstance(left,scantable) and not isinstance(right,scantable):
        print "Please provide two scantables as input"
        return
    from asap._asap import b_operate as _bop
    s = scantable(_bop(left, right, op, tsys))
    s._add_history("simple_math", varlist)
    return s
