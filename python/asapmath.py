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
    scanAv = True
    if kwargs.has_key('scanav'):
       scanAv = kwargs.get('scanav')
    weight = 'none'
    if kwargs.has_key('weight'):
       weight = kwargs.get('weight')
    mask = ()
    if kwargs.has_key('mask'):
        mask = kwargs.get('mask')
    lst = tuple(args)
    from asap._asap import average as _av
    for s in lst:
        if not isinstance(s,scantable):
            print "Please give a list of scantables"
            return
    return scantable(_av(lst, mask, scanAv, weight))

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
    from asap._asap import quotient as _quot
    return scantable(_quot(source, reference, preserve))

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
    if not isinstance(left,scantable) and not isinstance(right,scantable):
        print "Please provide two scantables as input"
        return
    from asap._asap import b_operate as _bop
    return scantable(_bop(left, right, op, tsys))

def scale(scan, factor, insitu=None, allaxes=None, tsys=True):
    """
    Return a scan where all spectra are scaled by the give 'factor'
    Parameters:
        scan:        a scantable
        factor:      the scaling factor
        insitu:      if False a new scantable is returned.
                     Otherwise, the scaling is done in-situ
                     The default is taken from .asaprc (False)
        allaxes:     if True apply to all spectra. Otherwise
                     apply only to the selected (beam/pol/if)spectra only.
                     The default is taken from .asaprc (True if none)
        tsys:        if True (default) then apply the operation to Tsys
                     as well as the data
    """
    if allaxes is None: allaxes = rcParams['scantable.allaxes']
    if insitu is None: insitu = rcParams['insitu']
    if not insitu:
        from asap._asap import scale as _scale
        return scantable(_scale(scan, factor, allaxes, tsys))
    else:
        from asap._asap import scale_insitu as _scale
        _scale(scan, factor, allaxes, tsys)
        return
        

def add(scan, offset, insitu=None, allaxes=None):
    """
    Return a scan where all spectra have the offset added
    Parameters:
        scan:        a scantable
        offset:      the offset
        insitu:      if False a new scantable is returned.
                     Otherwise, the scaling is done in-situ
                     The default is taken from .asaprc (False)
        allaxes:     if True apply to all spectra. Otherwise
                     apply only to the selected (beam/pol/if)spectra only
                     The default is taken from .asaprc (True if none)
    """
    if allaxes is None: allaxes = rcParams['scantable.allaxes']
    if insitu is None: insitu = rcParams['insitu']
    if not insitu:
        from asap._asap import add as _add
        return scantable(_add(scan, offset, allaxes))
    else:
        from asap._asap import add_insitu as _add
        _add(scan, offset, allaxes)
        return
        
def convert_flux(scan, jyperk=None, eta=None, d=None, insitu=None, allaxes=None):
    """
    Return a scan where all spectra are converted to either Jansky or Kelvin
        depending upon the flux units of the scan table.  By default the
        function tries to look the values up internally. If it can't find
        them (or if you want to over-ride), you must specify EITHER jyperk
        OR eta (and D which it will try to look up also if you don't
        set it).  jyperk takes precedence if you set both.
    Parameters:
        scan:        a scantable
        jyperk:      the Jy / K conversion factor
        eta:         the aperture efficiency  
        d:           the geomtric diameter (metres)
        insitu:      if False a new scantable is returned.
                     Otherwise, the scaling is done in-situ
                     The default is taken from .asaprc (False)
        allaxes:         if True apply to all spectra. Otherwise
                     apply only to the selected (beam/pol/if)spectra only
                     The default is taken from .asaprc (True if none)
    """
    if allaxes is None: allaxes = rcParams['scantable.allaxes']
    if insitu is None: insitu = rcParams['insitu']
    if jyperk is None: jyperk = -1.0
    if d is None: d = -1.0
    if eta is None: eta = -1.0
    if not insitu:
        from asap._asap import convertflux as _convert
        return scantable(_convert(scan, d, eta, jyperk, allaxes))
    else:
        from asap._asap import convertflux_insitu as _convert
        _convert(scan, d, eta, jyperk, allaxes)
        return

def gain_el(scan, poly=None, filename="", method="linear",
            insitu=None, allaxes=None):
    """
    Return a scan after applying a gain-elevation correction. The correction
    can be made via either a polynomial or a table-based interpolation 
    (and extrapolation if necessary).
    You specify polynomial coefficients, an ascii table or neither.
    If you specify neither, then a polynomial correction will be made
    with built in coefficients known for certain telescopes (an error will
    occur if the instrument is not known).
    Parameters:
        scan:        a scantable
        poly:        Polynomial coefficients (default None) to compute a
                     gain-elevation correction as a function of
                     elevation (in degrees).
        filename:    The name of an ascii file holding correction factors.
                     The first row of the ascii file must give the column 
                     names and these MUST include columns
                     "ELEVATION" (degrees) and "FACTOR" (multiply data by this)
                     somewhere.
                     The second row must give the data type of the column. Use
                     'R' for Real and 'I' for Integer.  An example file
                     would be:

                     TIME ELEVATION FACTOR
                     R R R
                     0.1 0 1.5
                     0.2 20 1.4
                     0.3 40 1.3
                     0.4 60 1.2
                     0.5 80 1.1
                     0.6 90 1.0
        method:      Interpolation method when correcting from a table. Values 
                     are  "nearest", "linear" (default), "cubic" and "spline"
        insitu:      if False a new scantable is returned.
                     Otherwise, the scaling is done in-situ
                     The default is taken from .asaprc (False)
        allaxes:         if True apply to all spectra. Otherwise
                     apply only to the selected (beam/pol/if) spectra only
                     The default is taken from .asaprc (True if none)
    """
    if allaxes is None: allaxes = rcParams['scantable.allaxes']
    if poly is None:
       poly = ()
    if insitu is None: insitu = rcParams['insitu']
    from os.path import expandvars
    filename = expandvars(filename)
    if not insitu:
        from asap._asap import gainel as _gainEl
        return scantable(_gainEl(scan, poly, filename, method, allaxes))
    else:
        from asap._asap import gainel_insitu as _gainEl
        _gainEl(scan, poly, filename, method, allaxes)
        return
        
def freq_align(scan, reftime=None, method='cubic', perif=False, insitu=None):
    """
        Return a scan where all rows have been aligned in frequency. The
        alignment frequency frame (e.g. LSRK) is that set by function
        set_freqframe.  
        scan:        a scantable
        reftime:     reference time to align at. By default, the time of
                     the first row of data is used.  
        method:      Interpolation method for regridding the spectra. Choose
                     from "nearest", "linear", "cubic" (default) and "spline"
        perif:       Generate aligners per freqID (no doppler tracking) or
                     per IF (scan-based doppler tracking)
        insitu:      if False a new scantable is returned.
                     Otherwise, the scaling is done in-situ
                     The default is taken from .asaprc (False)
    """
    if reftime is None: reftime = ''
    if insitu is None: insitu = rcParams['insitu']
    perfreqid = not perif
    if not insitu:
        from asap._asap import freq_align as _align
        return scantable(_align(scan, reftime, method, perfreqid))
    else:
        from asap._asap import freq_align_insitu as _align
        _align(scan, reftime, method, perfreqid)
        return
        
def opacity(scan, tau, insitu=None, allaxes=None):
    """
    Return a scan after applying an opacity correction.
    Parameters:
        scan:        a scantable
        tau:         Opacity from which the correction factor is exp(tau*ZD)
                     where ZD is the zenith-distance
        insitu:      if False a new scantable is returned.
                     Otherwise, the scaling is done in-situ
                     The default is taken from .asaprc (False)
        allaxes:     if True apply to all spectra. Otherwise
                     apply only to the selected (beam/pol/if)spectra only
                     The default is taken from .asaprc (True if none)
    """
    if allaxes is None: allaxes = rcParams['scantable.allaxes']
    if insitu is None: insitu = rcParams['insitu']
    if not insitu:
        from asap._asap import opacity as _opacity
        return scantable(_opacity(scan, tau, allaxes))
    else:
        from asap._asap import opacity_insitu as _opacity
        _opacity(scan, tau, allaxes)
        return
        
def bin(scan, width=5, insitu=None):
    """
    Return a scan where all spectra have been binned up
        width:       The bin width (default=5) in pixels
        insitu:      if False a new scantable is returned.
                     Otherwise, the scaling is done in-situ
                     The default is taken from .asaprc (False)
    """
    if insitu is None: insitu = rcParams['insitu']
    if not insitu:
        from asap._asap import bin as _bin
        return scantable(_bin(scan, width))
    else:
        from asap._asap import bin_insitu as _bin
        _bin(scan, width)
        return

def resample(scan, width=5, method='cubic', insitu=None):
    """
    Return a scan where all spectra have been binned up
        width:       The bin width (default=5) in pixels
        method:      Interpolation method when correcting from a table. Values 
                     are  "nearest", "linear", "cubic" (default) and "spline"
        insitu:      if False a new scantable is returned.
                     Otherwise, the scaling is done in-situ
                     The default is taken from .asaprc (False)
    """
    if insitu is None: insitu = rcParams['insitu']
    if not insitu:
        from asap._asap import resample as _resample
        return scantable(_resample(scan, method, width))
    else:
        from asap._asap import resample_insitu as _resample
        _resample(scan, method, width)
        return

def average_pol(scan, mask=None, weight='none', insitu=None):
    """
    Average the Polarisations together.
    The polarisation cursor of the output scan is set to 0
    Parameters:
        scan:        The scantable
        mask:        An optional mask defining the region, where the
                     averaging will be applied. The output will have all 
                     specified points masked. 
        weight:      Weighting scheme. 'none' (default), or 'var' (variance
                     weighted)
        insitu:      if False a new scantable is returned.
                     Otherwise, the scaling is done in-situ
                     The default is taken from .asaprc (False)
    Example:
        polav = average_pols(myscan)
    """
    if mask is None:
        mask = ()
    if insitu is None: insitu = rcParams['insitu']
    if not insitu:
        from asap._asap import averagepol as _avpol
        return scantable(_avpol(scan, mask, weight))
    else:
        from asap._asap import averagepol_insitu as _avpol
        _avpol(scan, mask, weight)
        return
    
def smooth(scan, kernel="hanning", width=5.0, insitu=None, allaxes=None):
    """
    Smooth the spectrum by the specified kernel (conserving flux).
    Parameters:
        scan:       The input scan
        kernel:     The type of smoothing kernel. Select from
                    'hanning' (default), 'gaussian' and 'boxcar'.
                    The first three characters are sufficient.
        width:      The width of the kernel in pixels. For hanning this is
                    ignored otherwise it defauls to 5 pixels.
                    For 'gaussian' it is the Full Width Half
                    Maximum. For 'boxcar' it is the full width.
        insitu:     if False a new scantable is returned.
                    Otherwise, the scaling is done in-situ
                    The default is taken from .asaprc (False)
        allaxes:    If True (default) apply to all spectra. Otherwise
                    apply only to the selected (beam/pol/if)spectra only
                    The default is taken from .asaprc (True if none)
    Example:
         none
    """
    if allaxes is None: allaxes = rcParams['scantable.allaxes']
    if insitu is None: insitu = rcParams['insitu']
    if not insitu:
        from asap._asap import smooth as _smooth
        return scantable(_smooth(scan,kernel,width,allaxes))
    else:
        from asap._asap import smooth_insitu as _smooth
        _smooth(scan,kernel,width,allaxes)
        return
    
def poly_baseline(scan, mask=None, order=0, insitu=None):
    """
    Return a scan which has been baselined (all rows) by a polynomial. 
    Parameters:
        scan:    a scantable
        mask:    an optional mask
        order:   the order of the polynomial (default is 0)
        insitu:      if False a new scantable is returned.
                     Otherwise, the scaling is done in-situ
                     The default is taken from .asaprc (False)
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
    sf = f.auto_fit(insitu)
    return sf
