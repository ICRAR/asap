from asap.scantable import scantable
from asap import rcParams
from asap import print_log

def average_time(*args, **kwargs):
    """
    Return the (time) average of a scan or list of scans. [in channels only]
    The cursor of the output scan is set to 0
    Parameters:
        one scan or comma separated  scans or a list of scans
        mask:     an optional mask (only used for 'var' and 'tsys' weighting)
        scanav:   True averages each scan separately.
                  False (default) averages all scans together,
        weight:   Weighting scheme.
                    'none'     (mean no weight)
                    'var'      (1/var(spec) weighted)
                    'tsys'     (1/Tsys**2 weighted)
                    'tint'     (integration time weighted)
                    'tintsys'  (Tint/Tsys**2)
                    'median'   ( median averaging)
        align:    align the spectra in velocity before averaging. It takes
                  the time of the first spectrum in the first scantable
                  as reference time.
    Example:
        # return a time averaged scan from scana and scanb
        # without using a mask
        scanav = average_time(scana,scanb)
	# or equivalent
	# scanav = average_time([scana, scanb])
        # return the (time) averaged scan, i.e. the average of
        # all correlator cycles
        scanav = average_time(scan, scanav=True)
    """
    scanav = False
    if kwargs.has_key('scanav'):
       scanav = kwargs.get('scanav')
    weight = 'tint'
    if kwargs.has_key('weight'):
       weight = kwargs.get('weight')
    mask = ()
    if kwargs.has_key('mask'):
        mask = kwargs.get('mask')
    align = False
    if kwargs.has_key('align'):
        align = kwargs.get('align')
    varlist = vars()
    if isinstance(args[0],list):
        lst = args[0]
    elif isinstance(args[0],tuple):
        lst = list(args[0])
    else:
        lst = list(args)

    del varlist["kwargs"]
    varlist["args"] = "%d scantables" % len(lst)
    # need special formatting here for history...

    from asap._asap import stmath
    stm = stmath()
    for s in lst:
        if not isinstance(s,scantable):
            msg = "Please give a list of scantables"
            if rcParams['verbose']:
                print msg
                return
            else:
                raise TypeError(msg)
    if scanav: scanav = "SCAN"
    else: scanav = "NONE"
    alignedlst = []
    if align:
        refepoch = lst[0].get_time(0)
        for scan in lst:
            alignedlst.append(scan.freq_align(refepoch,insitu=False))
    else:
        alignedlst = lst
    if weight.upper() == 'MEDIAN':
        # median doesn't support list of scantables - merge first
        merged = None
        if len(alignedlst) > 1:
            merged = merge(alignedlst)
        else:
            merged = alignedlst[0]
        s = scantable(stm._averagechannel(merged, 'MEDIAN', scanav))
        del merged
    else:
        s = scantable(stm._average(alignedlst, mask, weight.upper(), scanav))
    s._add_history("average_time",varlist)
    print_log()
    return s

def quotient(source, reference, preserve=True):
    """
    Return the quotient of a 'source' (signal) scan and a 'reference' scan.
    The reference can have just one scan, even if the signal has many. Otherwise
    they must have the same number of scans.
    The cursor of the output scan is set to 0
    Parameters:
        source:        the 'on' scan
        reference:     the 'off' scan
        preserve:      you can preserve (default) the continuum or
                       remove it.  The equations used are
                       preserve:  Output = Toff * (on/off) - Toff
                       remove:    Output = Toff * (on/off) - Ton
    """
    varlist = vars()
    from asap._asap import stmath
    stm = stmath()
    stm._setinsitu(False)
    s = scantable(stm._quotient(source, reference, preserve))
    s._add_history("quotient",varlist)
    print_log()
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
    print "simple_math is deprecated use +=/* instead."

def merge(*args):
    """
    Merge a list of scanatables, or comma-sperated scantables into one
    scnatble.
    Parameters:
        A list [scan1, scan2] or scan1, scan2.
    Example:
        myscans = [scan1, scan2]
	allscans = merge(myscans)
	# or equivalent
	sameallscans = merge(scan1, scan2)
    """
    varlist = vars()
    if isinstance(args[0],list):
        lst = tuple(args[0])
    elif isinstance(args[0],tuple):
        lst = args[0]
    else:
        lst = tuple(args)
    varlist["args"] = "%d scantables" % len(lst)
    # need special formatting her for history...
    from asap._asap import stmath
    stm = stmath()
    for s in lst:
        if not isinstance(s,scantable):
            msg = "Please give a list of scantables"
            if rcParams['verbose']:
                print msg
                return
            else:
                raise TypeError(msg)
    s = scantable(stm._merge(lst))
    s._add_history("merge", varlist)
    print_log()
    return s

