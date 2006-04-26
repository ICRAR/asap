from asap._asap import Scantable
from asap import rcParams
from asap import print_log, asaplog
from asap import selector
from numarray import ones,zeros
import sys

class scantable(Scantable):
    """
        The ASAP container for scans
    """

    def __init__(self, filename, average=None, unit=None):
        """
        Create a scantable from a saved one or make a reference
        Parameters:
            filename:    the name of an asap table on disk
                         or
                         the name of a rpfits/sdfits/ms file
                         (integrations within scans are auto averaged
                         and the whole file is read)
                         or
                         [advanced] a reference to an existing
                         scantable
            average:     average all integrations withinb a scan on read.
                         The default (True) is taken from .asaprc.
            unit:         brightness unit; must be consistent with K or Jy.
                         Over-rides the default selected by the reader
                         (input rpfits/sdfits/ms) or replaces the value
                         in existing scantables
        """
        if average is None:
            average = rcParams['scantable.autoaverage']
        varlist = vars()
        from asap._asap import stmath
        self._math = stmath()
        if isinstance(filename, Scantable):
            Scantable.__init__(self, filename)
        else:
            if isinstance(filename,str):
                import os.path
                filename = os.path.expandvars(filename)
                filename = os.path.expanduser(filename)
                if not os.path.exists(filename):
                    s = "File '%s' not found." % (filename)
                    if rcParams['verbose']:
                        asaplog.push(s)
                        print asaplog.pop().strip()
                        return
                    raise IOError(s)
                if os.path.isdir(filename):
                    # crude check if asap table
                    if os.path.exists(filename+'/table.info'):
                        Scantable.__init__(self, filename, "memory")
                        if unit is not None:
                            self.set_fluxunit(unit)
                        self.set_freqframe(rcParams['scantable.freqframe'])
                    else:
                        msg = "The given file '%s'is not a valid asap table." % (filename)
                        if rcParams['verbose']:
                            print msg
                            return
                        else:
                            raise IOError(msg)
                else:
                    self._fill([filename],unit, average)
            elif (isinstance(filename,list) or isinstance(filename,tuple)) \
                  and isinstance(filename[-1], str):
                self._fill(filename, unit, average)
        print_log()

    def save(self, name=None, format=None, overwrite=False):
        """
        Store the scantable on disk. This can be an asap (aips++) Table, SDFITS,
        Image FITS or MS2 format.
        Parameters:
            name:        the name of the outputfile. For format="FITS" this
                         is the directory file name into which all the files
                         will be written (default is 'asap_FITS'). For format
                         "ASCII" this is the root file name (data in 'name'.txt
                         and header in 'name'_header.txt)
            format:      an optional file format. Default is ASAP.
                         Allowed are - 'ASAP' (save as ASAP [aips++] Table),
                                       'SDFITS' (save as SDFITS file)
                                       'FITS' (saves each row as a FITS Image)
                                       'ASCII' (saves as ascii text file)
                                       'MS2' (saves as an aips++
                                              MeasurementSet V2)
            overwrite:   If the file should be overwritten if it exists.
                         The default False is to return with warning
                         without writing the output. USE WITH CARE.
        Example:
            scan.save('myscan.asap')
            scan.save('myscan.sdfits','SDFITS')
        """
        from os import path
        if format is None: format = rcParams['scantable.save']
        suffix = '.'+format.lower()
        if name is None or name =="":
            name = 'scantable'+suffix
            from asap import asaplog
            msg = "No filename given. Using default name %s..." % name
            asaplog.push(msg)
        name = path.expandvars(name)
        if path.isfile(name) or path.isdir(name):
            if not overwrite:
                msg = "File %s exists." % name
                if rcParams['verbose']:
                    print msg
                    return
                else:
                    raise IOError(msg)
        format2 = format.upper()
        if format2 == 'ASAP':
            self._save(name)
        else:
            from asap._asap import stwriter as stw
            w = stw(format2)
            w.write(self, name)
        print_log()
        return

    def copy(self):
        """
        Return a copy of this scantable.
        Parameters:
            none
        Example:
            copiedscan = scan.copy()
        """
        sd = scantable(Scantable._copy(self))
        return sd

    def get_scan(self, scanid=None):
        """
        Return a specific scan (by scanno) or collection of scans (by
        source name) in a new scantable.
        Parameters:
            scanid:    a (list of) scanno or a source name, unix-style
                       patterns are accepted for source name matching, e.g.
                       '*_R' gets all 'ref scans
        Example:
            # get all scans containing the source '323p459'
            newscan = scan.get_scan('323p459')
            # get all 'off' scans
            refscans = scan.get_scan('*_R')
            # get a susbset of scans by scanno (as listed in scan.summary())
            newscan = scan.get_scan([0,2,7,10])
        """
        if scanid is None:
            if rcParams['verbose']:
                print "Please specify a scan no or name to retrieve from the scantable"
                return
            else:
                raise RuntimeError("No scan given")

        try:
            bsel = self.get_selection()
            sel = selector()
            if type(scanid) is str:
                sel.set_name(scanid)
                self.set_selection(bsel+sel)
                scopy = self._copy()
                self.set_selection(bsel)
                return scantable(scopy)
            elif type(scanid) is int:
                sel.set_scans([scanid])
                self.set_selection(bsel+sel)
                scopy = self._copy()
                self.set_selection(bsel)
                return scantable(scopy)
            elif type(scanid) is list:
                sel.set_scans(scanid)
                self.set_selection(sel)
                scopy = self._copy()
                self.set_selection(bsel)
                return scantable(scopy)
            else:
                msg = "Illegal scanid type, use 'int' or 'list' if ints."
                if rcParams['verbose']:
                    print msg
                else:
                    raise TypeError(msg)
        except RuntimeError:
            if rcParams['verbose']: print "Couldn't find any match."
            else: raise

    def __str__(self):
        return Scantable._summary(self,True)

    def summary(self, filename=None):
        """
        Print a summary of the contents of this scantable.
        Parameters:
            filename:    the name of a file to write the putput to
                         Default - no file output
            verbose:     print extra info such as the frequency table
                         The default (False) is taken from .asaprc
        """
        info = Scantable._summary(self, True)
        #if verbose is None: verbose = rcParams['scantable.verbosesummary']
        if filename is not None:
            if filename is "":
                filename = 'scantable_summary.txt'
            from os.path import expandvars, isdir
            filename = expandvars(filename)
            if not isdir(filename):
                data = open(filename, 'w')
                data.write(info)
                data.close()
            else:
                msg = "Illegal file name '%s'." % (filename)
                if rcParams['verbose']:
                    print msg
                else:
                    raise IOError(msg)
        if rcParams['verbose']:
            try:
                from IPython.genutils import page as pager
            except ImportError:
                from pydoc import pager
            pager(info)
        else:
            return info


    def get_selection(self):
        """
        Get the selection object currently set on this scantable.
        Parameters:
            none
        Example:
            sel = scan.get_selection()
            sel.set_ifs(0)              # select IF 0
            scan.set_selection(sel)     # apply modified selection
        """
        return selector(self._getselection())

    def set_selection(self, selection=selector()):
        """
        Select a subset of the data. All following operations on this scantable
        are only applied to thi selection.
        Parameters:
            selection:    a selector object (default unset the selection)
        Examples:
            sel = selector()         # create a selection object
            self.set_scans([0,3])    # select SCANNO 0 and 3
            scan.set_selection(sel)  # set the selection
            scan.summary()           # will only print summary of scanno 0 an 3
            scan.set_selection()     # unset the selection
        """
        self._setselection(selection)

    def set_cursor(self, beam=0, IF=0, pol=0):
        print "DEPRECATED - use set_selection"

    def get_cursor(self):
        print "DEPRECATED - use get_selection"

    def stats(self, stat='stddev', mask=None):
        """
        Determine the specified statistic of the current beam/if/pol
        Takes a 'mask' as an optional parameter to specify which
        channels should be excluded.
        Parameters:
            stat:    'min', 'max', 'sumsq', 'sum', 'mean'
                     'var', 'stddev', 'avdev', 'rms', 'median'
            mask:    an optional mask specifying where the statistic
                     should be determined.
        Example:
            scan.set_unit('channel')
            msk = scan.create_mask([100,200],[500,600])
            scan.stats(stat='mean', mask=m)
        """
        from numarray import array,zeros,Float
        if mask == None:
            mask = []
        axes = ['Beam','IF','Pol','Time']
        if not self._check_ifs():
             raise ValueError("Cannot apply mask as the IFs have different number of channels"
                              "Please use setselection() to select individual IFs")

        statvals = self._math._stats(self, mask, stat)
        out = ''
        axes = []
        for i in range(self.nrow()):
            axis = []
            axis.append(self.getscan(i))
            axis.append(self.getbeam(i))
            axis.append(self.getif(i))
            axis.append(self.getpol(i))
            axis.append(self.getcycle(i))
            axes.append(axis)
            tm = self._gettime(i)
            src = self._getsourcename(i)
            out += 'Scan[%d] (%s) ' % (axis[0], src)
            out += 'Time[%s]:\n' % (tm)
            if self.nbeam(-1) > 1: out +=  ' Beam[%d] ' % (axis[1])
            if self.nif(-1) > 1: out +=  ' IF[%d] ' % (axis[2])
            if self.npol(-1) > 1: out +=  ' Pol[%d] ' % (axis[3])
            out += '= %3.3f\n' % (statvals[i])
            out +=  "--------------------------------------------------\n"

        if rcParams['verbose']:
            print "--------------------------------------------------"
            print " ",stat
            print "--------------------------------------------------"
            print out
        retval = { 'axesnames': ['scanno','beamno','ifno','polno','cycleno'],
                   'axes' : axes,
                   'data': statvals}
        return retval

    def stddev(self,mask=None):
        """
        Determine the standard deviation of the current beam/if/pol
        Takes a 'mask' as an optional parameter to specify which
        channels should be excluded.
        Parameters:
            mask:    an optional mask specifying where the standard
                     deviation should be determined.

        Example:
            scan.set_unit('channel')
            msk = scan.create_mask([100,200],[500,600])
            scan.stddev(mask=m)
        """
        return self.stats(stat='stddev',mask=mask);


    def column_names(self):
        """
        Return a  list of column names, which can be used for selection.
        """
        return list(Scantable.column_names(self))

    def get_tsys(self):
        """
        Return the System temperatures.
        Parameters:

        Returns:
            a list of Tsys values for the current selection
        """

        return self._row_callback(self._gettsys, "Tsys")

    def _row_callback(self, callback, label):
        axes = []
        axesnames = ['scanno','beamno','ifno','polno','cycleno']
        out = ""
        outvec =[]
        for i in range(self.nrow()):
            axis = []
            axis.append(self.getscan(i))
            axis.append(self.getbeam(i))
            axis.append(self.getif(i))
            axis.append(self.getpol(i))
            axis.append(self.getcycle(i))
            axes.append(axis)
            tm = self._gettime(i)
            src = self._getsourcename(i)
            out += 'Scan[%d] (%s) ' % (axis[0], src)
            out += 'Time[%s]:\n' % (tm)
            if self.nbeam(-1) > 1: out +=  ' Beam[%d] ' % (axis[1])
            if self.nif(-1) > 1: out +=  ' IF[%d] ' % (axis[2])
            if self.npol(-1) > 1: out +=  ' Pol[%d] ' % (axis[3])
            outvec.append(callback(i))
            out += '= %3.3f\n' % (outvec[i])
            out +=  "--------------------------------------------------\n"
        if rcParams['verbose']:
            print "--------------------------------------------------"
            print " %s" % (label)
            print "--------------------------------------------------"
            print out
        retval = {'axesnames': axesnames, 'axes': axes, 'data': outvec}
        return retval


    def get_time(self, row=-1):
        """
        Get a list of time stamps for the observations.
        Return a string for each integration in the scantable.
        Parameters:
            row:    row no of integration. Default -1 return all rows
        Example:
            none
        """
        out = []
        if row == -1:
            for i in range(self.nrow()):
                out.append(self._gettime(i))
            return out
        else:
            if row < self.nrow():
                return self._gettime(row)

    def get_sourcename(self, row=-1):
        """
        Get a list source names for the observations.
        Return a string for each integration in the scantable.
        Parameters:
            row:    row no of integration. Default -1 return all rows
        Example:
            none
        """
        out = []
        if row == -1:
            return [self._getsourcename(i) for i in range(self.nrow())]
        else:
            if  0 <= row < self.nrow():
                return self._getsourcename(row)

    def get_elevation(self, row=-1):
        """
        Get a list of elevations for the observations.
        Return a float for each integration in the scantable.
        Parameters:
            row:    row no of integration. Default -1 return all rows
        Example:
            none
        """
        out = []
        if row == -1:
            return [self._getelevation(i) for i in range(self.nrow())]
        else:
            if  0 <= row < self.nrow():
                return self._getelevation(row)

    def get_azimuth(self, row=-1):
        """
        Get a list of azimuths for the observations.
        Return a float for each integration in the scantable.
        Parameters:
            row:    row no of integration. Default -1 return all rows
        Example:
            none
        """
        out = []
        if row == -1:
            return [self._getazimuth(i) for i in range(self.nrow())]
        else:
            if  0 <= row < self.nrow():
                return self._getazimuth(row)

    def get_parangle(self, row=-1):
        """
        Get a list of parallactic angles for the observations.
        Return a float for each integration in the scantable.
        Parameters:
            row:    row no of integration. Default -1 return all rows
        Example:
            none
        """
        out = []
        if row == -1:
            return [self._getparangle(i) for i in range(self.nrow())]
        else:
            if  0 <= row < self.nrow():
                return self._getparangle(row)

    def set_unit(self, unit='channel'):
        """
        Set the unit for all following operations on this scantable
        Parameters:
            unit:    optional unit, default is 'channel'
                     one of '*Hz','km/s','channel', ''
        """
        varlist = vars()
        if unit in ['','pixel', 'channel']:
            unit = ''
        inf = list(self._getcoordinfo())
        inf[0] = unit
        self._setcoordinfo(inf)
        self._add_history("set_unit",varlist)

    def set_instrument(self, instr):
        """
        Set the instrument for subsequent processing
        Parameters:
            instr:    Select from 'ATPKSMB', 'ATPKSHOH', 'ATMOPRA',
                      'DSS-43' (Tid), 'CEDUNA', and 'HOBART'
        """
        self._setInstrument(instr)
        self._add_history("set_instument",vars())
        print_log()

    def set_doppler(self, doppler='RADIO'):
        """
        Set the doppler for all following operations on this scantable.
        Parameters:
            doppler:    One of 'RADIO', 'OPTICAL', 'Z', 'BETA', 'GAMMA'
        """
        varlist = vars()
        inf = list(self._getcoordinfo())
        inf[2] = doppler
        self._setcoordinfo(inf)
        self._add_history("set_doppler",vars())
        print_log()

    def set_freqframe(self, frame=None):
        """
        Set the frame type of the Spectral Axis.
        Parameters:
            frame:   an optional frame type, default 'LSRK'. Valid frames are:
                     'REST','TOPO','LSRD','LSRK','BARY',
                     'GEO','GALACTO','LGROUP','CMB'
        Examples:
            scan.set_freqframe('BARY')
        """
        if frame is None: frame = rcParams['scantable.freqframe']
        varlist = vars()
        valid = ['REST','TOPO','LSRD','LSRK','BARY', \
                   'GEO','GALACTO','LGROUP','CMB']

        if frame in valid:
            inf = list(self._getcoordinfo())
            inf[1] = frame
            self._setcoordinfo(inf)
            self._add_history("set_freqframe",varlist)
        else:
            msg  = "Please specify a valid freq type. Valid types are:\n",valid
            if rcParams['verbose']:
                print msg
            else:
                raise TypeError(msg)
        print_log()

    def set_dirframe(self, frame=""):
        """
        Set the frame type of the Direction on the sky.
        Parameters:
            frame:   an optional frame type, default ''. Valid frames are:
                     'J2000', 'B1950', 'GALACTIC'
        Examples:
            scan.set_dirframe('GALACTIC')
        """
        varlist = vars()
        try:
            Scantable.set_dirframe(self, frame)
        except RuntimeError,msg:
            if rcParams['verbose']:
                print msg
            else:
                raise
        self._add_history("set_dirframe",varlist)

    def get_unit(self):
        """
        Get the default unit set in this scantable
        Parameters:
        Returns:
            A unit string
        """
        inf = self._getcoordinfo()
        unit = inf[0]
        if unit == '': unit = 'channel'
        return unit

    def get_abcissa(self, rowno=0):
        """
        Get the abcissa in the current coordinate setup for the currently
        selected Beam/IF/Pol
        Parameters:
            rowno:    an optional row number in the scantable. Default is the
                      first row, i.e. rowno=0
        Returns:
            The abcissa values and it's format string (as a dictionary)
        """
        abc = self._getabcissa(rowno)
        lbl = self._getabcissalabel(rowno)
        print_log()
        return abc, lbl

    def flag(self, mask=[]):
        """
        Flag the selected data using an optional channel mask.
        Parameters:
            mask:   an optional channel mask, created with create_mask. Default
                    (no mask) is all channels.
        """
        varlist = vars()
        try:
            self._flag(mask)
        except RuntimeError,msg:
            if rcParams['verbose']:
                print msg
                return
            else: raise
        self._add_history("flag", varlist)


    def create_mask(self, *args, **kwargs):
        """
        Compute and return a mask based on [min,max] windows.
        The specified windows are to be INCLUDED, when the mask is
        applied.
        Parameters:
            [min,max],[min2,max2],...
                Pairs of start/end points (inclusive)specifying the regions
                to be masked
            invert:     optional argument. If specified as True,
                        return an inverted mask, i.e. the regions
                        specified are EXCLUDED
            row:        create the mask using the specified row for
                        unit conversions, default is row=0
                        only necessary if frequency varies over rows.
        Example:
            scan.set_unit('channel')
            a)
            msk = scan.create_mask([400,500],[800,900])
            # masks everything outside 400 and 500
            # and 800 and 900 in the unit 'channel'

            b)
            msk = scan.create_mask([400,500],[800,900], invert=True)
            # masks the regions between 400 and 500
            # and 800 and 900 in the unit 'channel'
            c)
            mask only channel 400
            msk =  scan.create_mask([400,400])
        """
        row = 0
        if kwargs.has_key("row"):
            row = kwargs.get("row")
        data = self._getabcissa(row)
        u = self._getcoordinfo()[0]
        if rcParams['verbose']:
            if u == "": u = "channel"
            from asap import asaplog
            msg = "The current mask window unit is %s" % u
            if not self._check_ifs():
                msg += "\nThis mask is only valid for IF=%d" % (self.getif(i))
            asaplog.push(msg)
        n = self.nchan()
        msk = zeros(n)
        # test if args is a 'list' or a 'normal *args - UGLY!!!

        ws = (isinstance(args[-1][-1],int) or isinstance(args[-1][-1],float)) and args or args[0]
        for window in ws:
            if (len(window) != 2 or window[0] > window[1] ):
                raise TypeError("A window needs to be defined as [min,max]")
            for i in range(n):
                if data[i] >= window[0] and data[i] <= window[1]:
                    msk[i] = 1
        if kwargs.has_key('invert'):
            if kwargs.get('invert'):
                from numarray import logical_not
                msk = logical_not(msk)
        print_log()
        return msk

    def get_restfreqs(self):
        """
        Get the restfrequency(s) stored in this scantable.
        The return value(s) are always of unit 'Hz'
        Parameters:
            none
        Returns:
            a list of doubles
        """
        return list(self._getrestfreqs())


    def set_restfreqs(self, freqs=None, unit='Hz'):
        """
        Set or replace the restfrequency specified and
        If the 'freqs' argument holds a scalar,
        then that rest frequency will be applied to all the selected
        data.  If the 'freqs' argument holds
        a vector, then it MUST be of equal or smaller length than
        the number of IFs (and the available restfrequencies will be
        replaced by this vector).  In this case, *all* data have
        the restfrequency set per IF according
        to the corresponding value you give in the 'freqs' vector.
        E.g. 'freqs=[1e9,2e9]'  would mean IF 0 gets restfreq 1e9 and
        IF 1 gets restfreq 2e9.
        You can also specify the frequencies via known line names
        from the built-in Lovas table.
        Parameters:
            freqs:   list of rest frequency values or string idenitfiers
            unit:    unit for rest frequency (default 'Hz')

        Example:
            # set the given restfrequency for the whole table
            scan.set_restfreqs(freqs=1.4e9)
            # If thee number of IFs in the data is >= 2 the IF0 gets the first
            # value IF1 the second...
            scan.set_restfreqs(freqs=[1.4e9,1.67e9])
            #set the given restfrequency for the whole table (by name)
            scan.set_restfreqs(freqs="OH1667")

        Note:
            To do more sophisticate Restfrequency setting, e.g. on a
            source and IF basis, use scantable.set_selection() before using
            this function.
            # provide your scantable is call scan
            selection = selector()
            selection.set_name("ORION*")
            selection.set_ifs([1])
            scan.set_selection(selection)
            scan.set_restfreqs(freqs=86.6e9)

        """
        varlist = vars()

        t = type(freqs)
        if isinstance(freqs, int) or isinstance(freqs,float):
           self._setrestfreqs(freqs, unit)
        elif isinstance(freqs, list) or isinstance(freqs,tuple):
            if isinstance(freqs[-1], int) or isinstance(freqs[-1],float):
                sel = selector()
                savesel = self._getselection()
                for i in xrange(len(freqs)):
                    sel.set_ifs([i])
                    self._setselection(sel)
                    self._setrestfreqs(freqs[i], unit)
                self._setselection(savesel)
            elif isinstance(freqs[-1], str):
                # not yet implemented
                pass
        else:
            return
        self._add_history("set_restfreqs", varlist)



    def history(self):
        hist = list(self._gethistory())
        out = "-"*80
        for h in hist:
            if h.startswith("---"):
                out += "\n"+h
            else:
                items = h.split("##")
                date = items[0]
                func = items[1]
                items = items[2:]
                out += "\n"+date+"\n"
                out += "Function: %s\n  Parameters:" % (func)
                for i in items:
                    s = i.split("=")
                    out += "\n   %s = %s" % (s[0],s[1])
                out += "\n"+"-"*80
        try:
            from IPython.genutils import page as pager
        except ImportError:
            from pydoc import pager
        pager(out)
        return

    #
    # Maths business
    #

    def average_time(self, mask=None, scanav=False, weight='tint', align=False):
        """
        Return the (time) average of a scan, or apply it 'insitu'.
        Note:
            in channels only
            The cursor of the output scan is set to 0.
        Parameters:
            one scan or comma separated  scans
            mask:     an optional mask (only used for 'var' and 'tsys'
                      weighting)
            scanav:   True averages each scan separately
                      False (default) averages all scans together,
            weight:   Weighting scheme. 'none', 'var' (1/var(spec)
                      weighted), 'tsys' (1/Tsys**2 weighted), 'tint'
                      (integration time weighted) or 'tintsys' (Tint/Tsys**2).
                      The default is 'tint'
            align:    align the spectra in velocity before averaging. It takes
                      the time of the first spectrum as reference time.
        Example:
            # time average the scantable without using a mask
            newscan = scan.average_time()
        """
        varlist = vars()
        if weight is None: weight = 'TINT'
        if mask is None: mask = ()
        if scanav:
            scanav = "SCAN"
        else:
            scanav = "NONE"
        scan = (self,)
        try:
          if align:
              scan = (self.freq_align(insitu=False),)
          s = scantable(self._math._average(scan, mask, weight.upper(),
                        scanav))
        except RuntimeError,msg:
            if rcParams['verbose']:
                print msg
                return
            else: raise
        s._add_history("average_time",varlist)
        print_log()
        return s

    def convert_flux(self, jyperk=None, eta=None, d=None, insitu=None):
        """
        Return a scan where all spectra are converted to either
        Jansky or Kelvin depending upon the flux units of the scan table.
        By default the function tries to look the values up internally.
        If it can't find them (or if you want to over-ride), you must
        specify EITHER jyperk OR eta (and D which it will try to look up
        also if you don't set it). jyperk takes precedence if you set both.
        Parameters:
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
        if insitu is None: insitu = rcParams['insitu']
        self._math._setinsitu(insitu)
        varlist = vars()
        if jyperk is None: jyperk = -1.0
        if d is None: d = -1.0
        if eta is None: eta = -1.0
        s = scantable(self._math._convertflux(self, d, eta, jyperk))
        s._add_history("convert_flux", varlist)
        print_log()
        if insitu: self._assign(s)
        else: return s

    def gain_el(self, poly=None, filename="", method="linear", insitu=None):
        """
        Return a scan after applying a gain-elevation correction.
        The correction can be made via either a polynomial or a
        table-based interpolation (and extrapolation if necessary).
        You specify polynomial coefficients, an ascii table or neither.
        If you specify neither, then a polynomial correction will be made
        with built in coefficients known for certain telescopes (an error
        will occur if the instrument is not known).
        The data and Tsys are *divided* by the scaling factors.
        Parameters:
            poly:        Polynomial coefficients (default None) to compute a
                         gain-elevation correction as a function of
                         elevation (in degrees).
            filename:    The name of an ascii file holding correction factors.
                         The first row of the ascii file must give the column
                         names and these MUST include columns
                         "ELEVATION" (degrees) and "FACTOR" (multiply data
                         by this) somewhere.
                         The second row must give the data type of the
                         column. Use 'R' for Real and 'I' for Integer.
                         An example file would be
                         (actual factors are arbitrary) :

                         TIME ELEVATION FACTOR
                         R R R
                         0.1 0 0.8
                         0.2 20 0.85
                         0.3 40 0.9
                         0.4 60 0.85
                         0.5 80 0.8
                         0.6 90 0.75
            method:      Interpolation method when correcting from a table.
                         Values are  "nearest", "linear" (default), "cubic"
                         and "spline"
            insitu:      if False a new scantable is returned.
                         Otherwise, the scaling is done in-situ
                         The default is taken from .asaprc (False)
        """

        if insitu is None: insitu = rcParams['insitu']
        self._math._setinsitu(insitu)
        varlist = vars()
        if poly is None:
           poly = ()
        from os.path import expandvars
        filename = expandvars(filename)
        s = scantable(self._math._gainel(self, poly, filename, method))
        s._add_history("gain_el", varlist)
        print_log()
        if insitu: self._assign(s)
        else: return s

    def freq_align(self, reftime=None, method='cubic', insitu=None):
        """
        Return a scan where all rows have been aligned in frequency/velocity.
        The alignment frequency frame (e.g. LSRK) is that set by function
        set_freqframe.
        Parameters:
            reftime:     reference time to align at. By default, the time of
                         the first row of data is used.
            method:      Interpolation method for regridding the spectra.
                         Choose from "nearest", "linear", "cubic" (default)
                         and "spline"
            insitu:      if False a new scantable is returned.
                         Otherwise, the scaling is done in-situ
                         The default is taken from .asaprc (False)
        """
        if insitu is None: insitu = rcParams["insitu"]
        self._math._setinsitu(insitu)
        varlist = vars()
        if reftime is None: reftime = ""
        s = scantable(self._math._freq_align(self, reftime, method))
        s._add_history("freq_align", varlist)
        print_log()
        if insitu: self._assign(s)
        else: return s

    def opacity(self, tau, insitu=None):
        """
        Apply an opacity correction. The data
        and Tsys are multiplied by the correction factor.
        Parameters:
            tau:         Opacity from which the correction factor is
                         exp(tau*ZD)
                         where ZD is the zenith-distance
            insitu:      if False a new scantable is returned.
                         Otherwise, the scaling is done in-situ
                         The default is taken from .asaprc (False)
        """
        if insitu is None: insitu = rcParams['insitu']
        self._math._setinsitu(insitu)
        varlist = vars()
        s = scantable(self._math._opacity(self, tau))
        s._add_history("opacity", varlist)
        print_log()
        if insitu: self._assign(s)
        else: return s

    def bin(self, width=5, insitu=None):
        """
        Return a scan where all spectra have been binned up.
            width:       The bin width (default=5) in pixels
            insitu:      if False a new scantable is returned.
                         Otherwise, the scaling is done in-situ
                         The default is taken from .asaprc (False)
        """
        if insitu is None: insitu = rcParams['insitu']
        self._math._setinsitu(insitu)
        varlist = vars()
        s = scantable(self._math._bin(self, width))
        s._add_history("bin",varlist)
        print_log()
        if insitu: self._assign(s)
        else: return s


    def resample(self, width=5, method='cubic', insitu=None):
        """
        Return a scan where all spectra have been binned up
            width:       The bin width (default=5) in pixels
            method:      Interpolation method when correcting from a table.
                         Values are  "nearest", "linear", "cubic" (default)
                         and "spline"
            insitu:      if False a new scantable is returned.
                         Otherwise, the scaling is done in-situ
                         The default is taken from .asaprc (False)
        """
        if insitu is None: insitu = rcParams['insitu']
        self._math._setinsitu(insitu)
        varlist = vars()
        s = scantable(self._math._resample(self, method, width))
        s._add_history("resample",varlist)
        print_log()
        if insitu: self._assign(s)
        else: return s


    def average_pol(self, mask=None, weight='none'):
        """
        Average the Polarisations together.
        Parameters:
            mask:        An optional mask defining the region, where the
                         averaging will be applied. The output will have all
                         specified points masked.
            weight:      Weighting scheme. 'none' (default), 'var' (1/var(spec)
                         weighted), or 'tsys' (1/Tsys**2 weighted)
        """
        varlist = vars()
        if mask is None:
            mask = ()
        s = scantable(self._math._averagepol(self, mask, weight.upper()))
        s._add_history("average_pol",varlist)
        print_log()
        return s

    def convert_pol(self, poltype=None):
        """
        Convert the data to a different polarisation type.
        Parameters:
            poltype:    The new polarisation type. Valid types are:
                        "linear", "stokes" and "circular"
        """
        varlist = vars()
        try:
            s = scantable(self._math._convertpol(self, poltype))
        except RuntimeError,msg:
            if rcParams['verbose']:
              print msg
              return
            else:
                raise
        s._add_history("convert_pol",varlist)
        print_log()
        return s

    def smooth(self, kernel="hanning", width=5.0, insitu=None):
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
        Example:
             none
        """
        if insitu is None: insitu = rcParams['insitu']
        self._math._setinsitu(insitu)
        varlist = vars()
        s = scantable(self._math._smooth(self,kernel.lower(),width))
        s._add_history("smooth", varlist)
        print_log()
        if insitu: self._assign(s)
        else: return s


    def poly_baseline(self, mask=None, order=0, insitu=None):
        """
        Return a scan which has been baselined (all rows) by a polynomial.
        Parameters:
            scan:       a scantable
            mask:       an optional mask
            order:      the order of the polynomial (default is 0)
            insitu:     if False a new scantable is returned.
                        Otherwise, the scaling is done in-situ
                        The default is taken from .asaprc (False)
            allaxes:    If True (default) apply to all spectra. Otherwise
                        apply only to the selected (beam/pol/if)spectra only
                        The default is taken from .asaprc (True if none)
        Example:
            # return a scan baselined by a third order polynomial,
            # not using a mask
            bscan = scan.poly_baseline(order=3)
        """
        if insitu is None: insitu = rcParams['insitu']
        varlist = vars()
        if mask is None:
            from numarray import ones
            mask = list(ones(self.nchan(-1)))
        from asap.asapfitter import fitter
        f = fitter()
        f.set_scan(self, mask)
        f.set_function(poly=order)
        s = f.auto_fit(insitu)
        s._add_history("poly_baseline", varlist)
        print_log()
        if insitu: self._assign(s)
        else: return s

    def auto_poly_baseline(self, mask=[], edge=(0,0), order=0,
                           threshold=3, insitu=None):
        """
        Return a scan which has been baselined (all rows) by a polynomial.
        Spectral lines are detected first using linefinder and masked out
        to avoid them affecting the baseline solution.

        Parameters:
            mask:       an optional mask retreived from scantable
            edge:       an optional number of channel to drop at
                        the edge of spectrum. If only one value is
                        specified, the same number will be dropped from
                        both sides of the spectrum. Default is to keep
                        all channels. Nested tuples represent individual
                        edge selection for different IFs (a number of spectral
                        channels can be different)
            order:      the order of the polynomial (default is 0)
            threshold:  the threshold used by line finder. It is better to
                        keep it large as only strong lines affect the
                        baseline solution.
            insitu:     if False a new scantable is returned.
                        Otherwise, the scaling is done in-situ
                        The default is taken from .asaprc (False)

        Example:
            scan2=scan.auto_poly_baseline(order=7)
        """
        if insitu is None: insitu = rcParams['insitu']
        varlist = vars()
        from asap.asapfitter import fitter
        from asap.asaplinefind import linefinder
        from asap import _is_sequence_or_number as _is_valid

        # check whether edge is set up for each IF individually
        individualEdge = False;
        if len(edge)>1:
           if isinstance(edge[0],list) or isinstance(edge[0],tuple):
               individualEdge = True;

        if not _is_valid(edge, int) and not individualEdge:
            raise ValueError, "Parameter 'edge' has to be an integer or a \
            pair of integers specified as a tuple. Nested tuples are allowed \
            to make individual selection for different IFs."

        curedge = (0,0)
        if individualEdge:
           for edge_par in edge:
               if not _is_valid(edge,int):
                  raise ValueError, "Each element of the 'edge' tuple has \
                  to be a pair of integers or an integer."
        else:
           curedge = edge;

        # setup fitter
        f = fitter()
        f.set_function(poly=order)

        # setup line finder
        fl=linefinder()
        fl.set_options(threshold=threshold)

        if not insitu:
            workscan=self.copy()
        else:
            workscan=self

        fl.set_scan(workscan)

        rows=range(workscan.nrow())
        from asap import asaplog
        asaplog.push("Processing:")
        for r in rows:
            msg = " Scan[%d] Beam[%d] IF[%d] Pol[%d] Cycle[%d]" %        (workscan.getscan(r),workscan.getbeam(r),workscan.getif(r),workscan.getpol(r), workscan.getcycle(r))
            asaplog.push(msg, False)

            # figure out edge parameter
            if individualEdge:
               if len(edge)>=workscan.getif(r):
                  raise RuntimeError, "Number of edge elements appear to be less than the number of IFs"
                  curedge = edge[workscan.getif(r)]

            # setup line finder
            fl.find_lines(r,mask,curedge)
            f.set_scan(workscan, fl.get_mask())
            f.x = workscan._getabcissa(r)
            f.y = workscan._getspectrum(r)
            f.data = None
            f.fit()
            x = f.get_parameters()
            workscan._setspectrum(f.fitter.getresidual(), r)
        workscan._add_history("poly_baseline", varlist)
        if insitu:
            self._assign(workscan)
        else:
            return workscan

    def rotate_linpolphase(self, angle):
        """
        Rotate the phase of the complex polarization O=Q+iU correlation.
        This is always done in situ in the raw data.  So if you call this
        function more than once then each call rotates the phase further.
        Parameters:
            angle:   The angle (degrees) to rotate (add) by.
        Examples:
            scan.rotate_linpolphase(2.3)
        """
        varlist = vars()
        self._math._rotate_linpolphase(self, angle)
        self._add_history("rotate_linpolphase", varlist)
        print_log()
        return


    def rotate_xyphase(self, angle):
        """
        Rotate the phase of the XY correlation.  This is always done in situ
        in the data.  So if you call this function more than once
        then each call rotates the phase further.
        Parameters:
            angle:   The angle (degrees) to rotate (add) by.
        Examples:
            scan.rotate_xyphase(2.3)
        """
        varlist = vars()
        self._math._rotate_xyphase(self, angle)
        self._add_history("rotate_xyphase", varlist)
        print_log()
        return

    def swap_linears(self):
        """
        Swap the linear polarisations XX and YY
        """
        varlist = vars()
        self._math._swap_linears(self)
        self._add_history("swap_linears", varlist)
        print_log()
        return

    def invert_phase(self):
        """
        Invert the phase of the complex polarisation
        """
        varlist = vars()
        self._math._invert_phase(self)
        self._add_history("invert_phase", varlist)
        print_log()
        return

    def add(self, offset, insitu=None):
        """
        Return a scan where all spectra have the offset added
        Parameters:
            offset:      the offset
            insitu:      if False a new scantable is returned.
                         Otherwise, the scaling is done in-situ
                         The default is taken from .asaprc (False)
        """
        if insitu is None: insitu = rcParams['insitu']
        self._math._setinsitu(insitu)
        varlist = vars()
        s = scantable(self._math._unaryop(self, offset, "ADD", False))
        s._add_history("add",varlist)
        print_log()
        if insitu:
            self._assign(s)
        else:
            return s

    def scale(self, factor, tsys=True, insitu=None,):
        """
        Return a scan where all spectra are scaled by the give 'factor'
        Parameters:
            factor:      the scaling factor
            insitu:      if False a new scantable is returned.
                         Otherwise, the scaling is done in-situ
                         The default is taken from .asaprc (False)
            tsys:        if True (default) then apply the operation to Tsys
                         as well as the data
        """
        if insitu is None: insitu = rcParams['insitu']
        self._math._setinsitu(insitu)
        varlist = vars()
        s = scantable(self._math._unaryop(self, factor, "MUL", tsys))
        s._add_history("scale",varlist)
        print_log()
        if insitu:
            self._assign(s)
        else:
            return s

    def auto_quotient(self, mode='time', preserve=True):
        """
        This function allows to build quotients automatically.
        It assumes the observation to have the same numer of
        "ons" and "offs"
        It will support "closest off in time" in the future
        Parameters:
            mode:           the on/off detection mode; 'suffix' (default)
                            'suffix' identifies 'off' scans by the
                            trailing '_R' (Mopra/Parkes) or
                            '_e'/'_w' (Tid)
            preserve:       you can preserve (default) the continuum or
                            remove it.  The equations used are
                            preserve: Output = Toff * (on/off) - Toff
                            remove:   Output = Tref * (on/off) - Ton
        """
        modes = ["time"]
        if not mode in modes:
            msg = "please provide valid mode. Valid modes are %s" % (modes)
            raise ValueError(msg)
        varlist = vars()
        s = scantable(self._math._quotient(self, mode, preserve))
        s._add_history("auto_quotient",varlist)
        print_log()
        return s




    def freq_switch(self, insitu=None):
        """
        Apply frequency switching to the data.
        Parameters:
            insitu:      if False a new scantable is returned.
                         Otherwise, the swictching is done in-situ
                         The default is taken from .asaprc (False)
        Example:
            none
        """
        if insitu is None: insitu = rcParams['insitu']
        self._math._setinsitu(insitu)
        varlist = vars()
        s = scantable(self._math._freqswitch(self))
        s._add_history("freq_switch",varlist)
        print_log()
        if insitu: self._assign(s)
        else: return s

    def recalc_azel(self):
        """
        Recalculate the azimuth and elevation for each position.
        Parameters:
            none
        Example:
        """
        varlist = vars()
        self._recalcazel()
        self._add_history("recalc_azel", varlist)
        print_log()
        return

    def __add__(self, other):
        varlist = vars()
        s = None
        if isinstance(other, scantable):
            print "scantable + scantable NYI"
            return
        elif isinstance(other, float):
            s = scantable(self._math._unaryop(self, other, "ADD", False))
        else:
            raise TypeError("Other input is not a scantable or float value")
        s._add_history("operator +", varlist)
        print_log()
        return s

    def __sub__(self, other):
        """
        implicit on all axes and on Tsys
        """
        varlist = vars()
        s = None
        if isinstance(other, scantable):
            print "scantable - scantable NYI"
            return
        elif isinstance(other, float):
            s = scantable(self._math._unaryop(self, other, "SUB", False))
        else:
            raise TypeError("Other input is not a scantable or float value")
        s._add_history("operator -", varlist)
        print_log()
        return s

    def __mul__(self, other):
        """
        implicit on all axes and on Tsys
        """
        varlist = vars()
        s = None
        if isinstance(other, scantable):
            print "scantable * scantable NYI"
            return
        elif isinstance(other, float):
            s = scantable(self._math._unaryop(self, other, "MUL", False))
        else:
            raise TypeError("Other input is not a scantable or float value")
        s._add_history("operator *", varlist)
        print_log()
        return s


    def __div__(self, other):
        """
        implicit on all axes and on Tsys
        """
        varlist = vars()
        s = None
        if isinstance(other, scantable):
            print "scantable / scantable NYI"
            return
        elif isinstance(other, float):
            if other == 0.0:
                raise ZeroDivisionError("Dividing by zero is not recommended")
            s = scantable(self._math._unaryop(self, other, "DIV", False))
        else:
            raise TypeError("Other input is not a scantable or float value")
        s._add_history("operator /", varlist)
        print_log()
        return s

    def get_fit(self, row=0):
        """
        Print or return the stored fits for a row in the scantable
        Parameters:
            row:    the row which the fit has been applied to.
        """
        if row > self.nrow():
            return
        from asap.asapfit import asapfit
        fit = asapfit(self._getfit(row))
        if rcParams['verbose']:
            print fit
            return
        else:
            return fit.as_dict()

    def _add_history(self, funcname, parameters):
        # create date
        sep = "##"
        from datetime import datetime
        dstr = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
        hist = dstr+sep
        hist += funcname+sep#cdate+sep
        if parameters.has_key('self'): del parameters['self']
        for k,v in parameters.iteritems():
            if type(v) is dict:
                for k2,v2 in v.iteritems():
                    hist += k2
                    hist += "="
                    if isinstance(v2,scantable):
                        hist += 'scantable'
                    elif k2 == 'mask':
                        if isinstance(v2,list) or isinstance(v2,tuple):
                            hist += str(self._zip_mask(v2))
                        else:
                            hist += str(v2)
                    else:
                        hist += str(v2)
            else:
                hist += k
                hist += "="
                if isinstance(v,scantable):
                    hist += 'scantable'
                elif k == 'mask':
                    if isinstance(v,list) or isinstance(v,tuple):
                        hist += str(self._zip_mask(v))
                    else:
                        hist += str(v)
                else:
                    hist += str(v)
            hist += sep
        hist = hist[:-2] # remove trailing '##'
        self._addhistory(hist)


    def _zip_mask(self, mask):
        mask = list(mask)
        i = 0
        segments = []
        while mask[i:].count(1):
            i += mask[i:].index(1)
            if mask[i:].count(0):
                j = i + mask[i:].index(0)
            else:
                j = len(mask)
            segments.append([i,j])
            i = j
        return segments

    def _get_ordinate_label(self):
        fu = "("+self.get_fluxunit()+")"
        import re
        lbl = "Intensity"
        if re.match(".K.",fu):
            lbl = "Brightness Temperature "+ fu
        elif re.match(".Jy.",fu):
            lbl = "Flux density "+ fu
        return lbl

    def _check_ifs(self):
        nchans = [self.nchan(i) for i in range(self.nif(-1))]
        nchans = filter(lambda t: t > 0, nchans)
        return (sum(nchans)/len(nchans) == nchans[0])

    def _fill(self, names, unit, average):
        import os
        varlist = vars()
        from asap._asap import stfiller
        first = True
        fullnames = []
        for name in names:
            name = os.path.expandvars(name)
            name = os.path.expanduser(name)
            if not os.path.exists(name):
                msg = "File '%s' does not exists" % (name)
                if rcParams['verbose']:
                    asaplog.push(msg)
                    print asaplog.pop().strip()
                    return
                raise IOError(msg)
            fullnames.append(name)
        if average:
            asaplog.push('Auto averaging integrations')
        for name in fullnames:
            r = stfiller()
            msg = "Importing %s..." % (name)
            asaplog.push(msg,False)
            print_log()
            r._open(name,-1,-1)
            r._read()
            tbl = r._getdata()
            if average:
                tbl = self._math._average((tbl,),(),'NONE','SCAN')
                #tbl = tbl2
            if not first:
                tbl = self._math._merge([self, tbl])
                #tbl = tbl2
            Scantable.__init__(self, tbl)
            r._close()
            del r,tbl
            first = False
        if unit is not None:
            self.set_fluxunit(unit)
        self.set_freqframe(rcParams['scantable.freqframe'])
        #self._add_history("scantable", varlist)

