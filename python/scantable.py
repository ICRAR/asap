from asap._asap import sdtable
from asap import rcParams
from numarray import ones,zeros
import sys

class scantable(sdtable):
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
        if average is None or type(average) is not bool:
            average = rcParams['scantable.autoaverage']            

        varlist = vars()
        self._vb = rcParams['verbose']
        self._p = None

        if isinstance(filename,sdtable):
            sdtable.__init__(self, filename)            
            if unit is not None:
                self.set_fluxunit(unit)                       
        else:
            import os.path
            if not os.path.exists(filename):
                print "File '%s' not found." % (filename)
                return
            filename = os.path.expandvars(filename)
            if os.path.isdir(filename):
                # crude check if asap table
                if os.path.exists(filename+'/table.info'):
                    sdtable.__init__(self, filename)
                    if unit is not None:
                        self.set_fluxunit(unit)                       
                else:
                    print "The given file '%s'is not a valid asap table." % (filename)
                    return
            else:            
                from asap._asap import sdreader
                ifSel = -1
                beamSel = -1
                r = sdreader(filename,ifSel,beamSel)
                print 'Importing data...'
                r._read([-1])
                tbl = r._getdata()
                if unit is not None:
                    tbl.set_fluxunit(unit)
                if average:
                    from asap._asap import average as _av
                    print 'Auto averaging integrations...'
                    tbl2 = _av((tbl,),(),True,'none')
                    sdtable.__init__(self,tbl2)
                    del tbl2
                else:
                    sdtable.__init__(self,tbl)
                del r,tbl
                self._add_history("scantable", varlist)

    def save(self, name=None, format=None, stokes=False, overwrite=False):
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
            stokes:      Convert to Stokes parameters (only available
                         currently with FITS and ASCII formats.
                         Default is False.
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
            print "No filename given. Using default name %s..." % name
        name = path.expandvars(name)
        if path.isfile(name) or path.isdir(name):
            if not overwrite:
                print "File %s already exists." % name
                return
        format2 = format.upper()
        if format2 == 'ASAP':
            self._save(name)
        else:
            from asap._asap import sdwriter as _sw
            w = _sw(format2)
            w.write(self, name, stokes)
        return

    def copy(self):
        """
        Return a copy of this scantable.
        Parameters:
            none
        Example:
            copiedscan = scan.copy()
        """
        sd = scantable(sdtable._copy(self))
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
            print "Please specify a scan no or name to retrieve from the scantable"
        try:
            if type(scanid) is str:
                s = sdtable._getsource(self,scanid)
                return scantable(s)
            elif type(scanid) is int:
                s = sdtable._getscan(self,[scanid])
                return scantable(s)
            elif type(scanid) is list:
                s = sdtable._getscan(self,scanid)
                return scantable(s)
            else:
                print "Illegal scanid type, use 'int' or 'list' if ints."
        except RuntimeError:
            print "Couldn't find any match."

    def __str__(self):
        return sdtable._summary(self,True)

    def summary(self,filename=None, verbose=None):
        """
        Print a summary of the contents of this scantable.
        Parameters:
            filename:    the name of a file to write the putput to
                         Default - no file output
            verbose:     print extra info such as the frequency table
                         The default (False) is taken from .asaprc
        """
        info = sdtable._summary(self, verbose)
        if verbose is None: verbose = rcParams['scantable.verbosesummary']
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
                print "Illegal file name '%s'." % (filename)
        print info
            
    def set_cursor(self, beam=0, IF=0, pol=0):
        """
        Set the spectrum for individual operations.
        Parameters:
            beam, IF, pol:    a number
        Example:
            scan.set_cursor(0,0,1)
            pol1sig = scan.stats(all=False) # returns std dev for beam=0
                                            # if=0, pol=1
        """
        varlist = vars()
        self.setbeam(beam)
        self.setpol(pol)
        self.setif(IF)
        self._add_history("set_cursor",varlist)
        return

    def get_cursor(self):
        """
        Return/print a the current 'cursor' into the Beam/IF/Pol cube.
        Parameters:
            none
        Returns:
            a list of values (currentBeam,currentIF,currentPol)
        Example:
            none
        """
        i = self.getbeam()
        j = self.getif()
        k = self.getpol()
        if self._vb:
            print "--------------------------------------------------"
            print " Cursor position"
            print "--------------------------------------------------"
            out = 'Beam=%d IF=%d Pol=%d ' % (i,j,k)
            print out
        return i,j,k

    def stats(self, stat='stddev', mask=None, allaxes=None):
        """
        Determine the specified statistic of the current beam/if/pol
        Takes a 'mask' as an optional parameter to specify which
        channels should be excluded.
        Parameters:
            stat:    'min', 'max', 'sumsq', 'sum', 'mean'
                     'var', 'stddev', 'avdev', 'rms', 'median'
            mask:    an optional mask specifying where the statistic
                     should be determined.
            allaxes: if True apply to all spectra. Otherwise
                     apply only to the selected (beam/pol/if)spectra only.
                     The default is taken from .asaprc (True if none)
        Example:
            scan.set_unit('channel')
            msk = scan.create_mask([100,200],[500,600])
            scan.stats(stat='mean', mask=m)
        """
        if allaxes is None: allaxes = rcParams['scantable.allaxes']
        from asap._asap import stats as _stats
        from numarray import array,zeros,Float
        if mask == None:
            mask = ones(self.nchan())
        axes = ['Beam','IF','Pol','Time']

        beamSel,IFSel,polSel = (self.getbeam(),self.getif(),self.getpol())
        if allaxes:
            n = self.nbeam()*self.nif()*self.npol()*self.nrow()
            shp = [self.nbeam(),self.nif(),self.npol(),self.nrow()]
            arr = array(zeros(n),shape=shp,type=Float)

            for i in range(self.nbeam()):
                self.setbeam(i)
                for j in range(self.nif()):
                    self.setif(j)
                    for k in range(self.npol()):
                        self.setpol(k)
                        arr[i,j,k,:] = _stats(self,mask,stat,-1)
            retval = {'axes': axes, 'data': arr, 'cursor':None}
            tm = [self._gettime(val) for val in range(self.nrow())]
            if self._vb:
                self._print_values(retval,stat,tm)
            self.setbeam(beamSel)
            self.setif(IFSel)
            self.setpol(polSel)
            return retval

        else:
            statval = _stats(self,mask,stat,-1)
            out = ''
            for l in range(self.nrow()):
                tm = self._gettime(l)
                out += 'Time[%s]:\n' % (tm)
                if self.nbeam() > 1: out +=  ' Beam[%d] ' % (beamSel)
                if self.nif() > 1: out +=  ' IF[%d] ' % (IFSel)
                if self.npol() > 1: out +=  ' Pol[%d] ' % (polSel)
                out += '= %3.3f\n' % (statval[l])
                out +=  "--------------------------------------------------\n"

            if self._vb:
                print "--------------------------------------------------"
                print " ",stat
                print "--------------------------------------------------"
                print out
            retval = {'axes': axes, 'data': array(statval), 'cursor':(i,j,k)}
            return retval

    def stddev(self,mask=None, allaxes=None):
        """
        Determine the standard deviation of the current beam/if/pol
        Takes a 'mask' as an optional parameter to specify which
        channels should be excluded.
        Parameters:
            mask:    an optional mask specifying where the standard
                     deviation should be determined.
            allaxes: optional flag to show all or a cursor selected
                     spectrum of Beam/IF/Pol. Default is all or taken
                     from .asaprc

        Example:
            scan.set_unit('channel')
            msk = scan.create_mask([100,200],[500,600])
            scan.stddev(mask=m)
        """
        if allaxes is None: allaxes = rcParams['scantable.allaxes']
        return self.stats(stat='stddev',mask=mask, allaxes=allaxes);

    def get_tsys(self, allaxes=None):
        """
        Return the System temperatures.
        Parameters:
           allaxes:     if True apply to all spectra. Otherwise
                        apply only to the selected (beam/pol/if)spectra only.
                        The default is taken from .asaprc (True if none)
        Returns:
            a list of Tsys values.
        """
        if allaxes is None: allaxes = rcParams['scantable.allaxes']
        from numarray import array,zeros,Float
        axes = ['Beam','IF','Pol','Time']

        if allaxes:
            n = self.nbeam()*self.nif()*self.npol()*self.nrow()
            shp = [self.nbeam(),self.nif(),self.npol(),self.nrow()]
            arr = array(zeros(n),shape=shp,type=Float)

            for i in range(self.nbeam()):
                self.setbeam(i)
                for j in range(self.nif()):
                    self.setif(j)
                    for k in range(self.npol()):
                        self.setpol(k)
                        arr[i,j,k,:] = self._gettsys()
            retval = {'axes': axes, 'data': arr, 'cursor':None}
            tm = [self._gettime(val) for val in range(self.nrow())]
            if self._vb:
                self._print_values(retval,'Tsys',tm)
            return retval

        else:
            i,j,k = (self.getbeam(),self.getif(),self.getpol())
            statval = self._gettsys()
            out = ''
            for l in range(self.nrow()):
                tm = self._gettime(l)
                out += 'Time[%s]:\n' % (tm)
                if self.nbeam() > 1: out +=  ' Beam[%d] ' % (i)
                if self.nif() > 1: out +=  ' IF[%d] ' % (j)
                if self.npol() > 1: out +=  ' Pol[%d] ' % (k)
                out += '= %3.3f\n' % (statval[l])
                out +=  "--------------------------------------------------\n"

            if self._vb:
                print "--------------------------------------------------"
                print " TSys"
                print "--------------------------------------------------"
                print out
            retval = {'axes': axes, 'data': array(statval), 'cursor':(i,j,k)}
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
        if self._p: self.plot()
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
        if self._p: self.plot()
        self._add_history("set_doppler",vars())
 
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
            print "Please specify a valid freq type. Valid types are:\n",valid
            
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
        return abc, lbl
        #return {'abcissa':abc,'label':lbl}

    def create_mask(self, *args, **kwargs):
        """
        Compute and return a mask based on [min,max] windows.
        The specified windows are to be INCLUDED, when the mask is
        applied.
        Parameters:
            [min,max],[min2,max2],...
                Pairs of start/end points specifying the regions
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
            msk = scan.set_mask([400,500],[800,900])
            # masks everything outside 400 and 500
            # and 800 and 900 in the unit 'channel'

            b)
            msk = scan.set_mask([400,500],[800,900], invert=True)
            # masks the regions between 400 and 500
            # and 800 and 900 in the unit 'channel'
           
        """
        row = 0
        if kwargs.has_key("row"):
            row = kwargs.get("row")
        data = self._getabcissa(row)
        u = self._getcoordinfo()[0]
        if self._vb:
            if u == "": u = "channel"
            print "The current mask window unit is", u
        n = self.nchan()
        msk = zeros(n)
        for window in args:
            if (len(window) != 2 or window[0] > window[1] ):
                print "A window needs to be defined as [min,max]"
                return
            for i in range(n):
                if data[i] >= window[0] and data[i] < window[1]:
                    msk[i] = 1
        if kwargs.has_key('invert'):
            if kwargs.get('invert'):
                from numarray import logical_not
                msk = logical_not(msk)            
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

    def lines(self):
        """
        Print the list of known spectral lines
        """
        sdtable._lines(self)

    def set_restfreqs(self, freqs=None, unit='Hz', lines=None, source=None,
                      theif=None):
        """
        Select the restfrequency for the specified source and IF OR
        replace for all IFs.  If the 'freqs' argument holds a scalar,
        then that rest frequency will be applied to the selected
        data (and added to the list of available rest frequencies).
        In this way, you can set a rest frequency for each
        source and IF combination.   If the 'freqs' argument holds
        a vector, then it MUST be of length the number of IFs
        (and the available restfrequencies will be replaced by
        this vector).  In this case, *all* data ('source' and 
        'theif' are ignored) have the restfrequency set per IF according 
        to the corresponding value you give in the 'freqs' vector.  
        E.g. 'freqs=[1e9,2e9]'  would mean IF 0 gets restfreq 1e9 and 
        IF 1 gets restfreq 2e9.

        You can also specify the frequencies via known line names
        in the argument 'lines'.  Use 'freqs' or 'lines'.  'freqs'
        takes precedence. See the list of known names via function
        scantable.lines()
        Parameters:
            freqs:   list of rest frequencies
            unit:    unit for rest frequency (default 'Hz')
            lines:   list of known spectral lines names (alternative to freqs).
                     See possible list via scantable.lines()
            source:  Source name (blank means all)
            theif:   IF (-1 means all)
        Example:
            scan.set_restfreqs(freqs=1.4e9, source='NGC253', theif=2)
            scan.set_restfreqs(freqs=[1.4e9,1.67e9])
        """
        varlist = vars()
        if source is None:
            source = ""
        if theif is None:
            theif = -1
        t = type(freqs)
        if t is int or t is float:
           freqs = [freqs]
        if freqs is None:
           freqs = []
        t = type(lines)
        if t is str:
           lines = [lines]
        if lines is None:
           lines = []
        sdtable._setrestfreqs(self, freqs, unit, lines, source, theif)
        self._add_history("set_restfreqs", varlist)
        


    def flag_spectrum(self, thebeam, theif, thepol):
        """
        This flags a selected spectrum in the scan 'for good'.
        USE WITH CARE - not reversible.
        Use masks for non-permanent exclusion of channels.
        Parameters:
            thebeam,theif,thepol:    all have to be explicitly
                                     specified
        Example:
            scan.flag_spectrum(0,0,1)
            flags the spectrum for Beam=0, IF=0, Pol=1
        """
        if (thebeam < self.nbeam() and
            theif < self.nif() and
            thepol < self.npol()):
            sdtable.setbeam(self, thebeam)
            sdtable.setif(self, theif)
            sdtable.setpol(self, thepol)
            sdtable._flag(self)
            self._add_history("flag_spectrum", vars())
        else:
            print "Please specify a valid (Beam/IF/Pol)"
        return

    def plot(self, what='spectrum',col='Pol', panel=None):
        """
        Plot the spectra contained in the scan. Alternatively you can also
        Plot Tsys vs Time
        Parameters:
            what:     a choice of 'spectrum' (default) or 'tsys'
            col:      which out of Beams/IFs/Pols should be colour stacked
            panel:    set up multiple panels, currently not working.
        """
        print "Warning! Not fully functional. Use plotter.plot() instead"
        
        validcol = {'Beam':self.nbeam(),'IF':self.nif(),'Pol':self.npol()}

        validyax = ['spectrum','tsys']
        from asap.asaplot import ASAPlot
        if not self._p:
            self._p = ASAPlot()
            #print "Plotting not enabled"
            #return
        if self._p.is_dead:
            del self._p
            self._p = ASAPlot()
        npan = 1
        x = None
        if what == 'tsys':
            n = self.nrow()
            if n < 2:
                print "Only one integration. Can't plot."
                return
        self._p.hold()
        self._p.clear()
        if panel == 'Time':
            npan = self.nrow()
            self._p.set_panels(rows=npan)
        xlab,ylab,tlab = None,None,None
        self._vb = False
        sel = self.get_cursor()        
        for i in range(npan):
            if npan > 1:
                self._p.subplot(i)
            for j in range(validcol[col]):
                x = None
                y = None
                m = None
                tlab = self._getsourcename(i)
                import re
                tlab = re.sub('_S','',tlab)
                if col == 'Beam':
                    self.setbeam(j)
                elif col == 'IF':
                    self.setif(j)
                elif col == 'Pol':
                    self.setpol(j)
                if what == 'tsys':
                    x = range(self.nrow())
                    xlab = 'Time [pixel]'
                    m = list(ones(len(x)))
                    y = []
                    ylab = r'$T_{sys}$'
                    for k in range(len(x)):
                        y.append(self._gettsys(k))
                else:
                    x,xlab = self.get_abcissa(i)
                    y = self._getspectrum(i)
                    ylab = r'Flux'
                    m = self._getmask(i)
                llab = col+' '+str(j)
                self._p.set_line(label=llab)
                self._p.plot(x,y,m)
            self._p.set_axes('xlabel',xlab)
            self._p.set_axes('ylabel',ylab)
            self._p.set_axes('title',tlab)
        self._p.release()
        self.set_cursor(sel[0],sel[1],sel[2])
        self._vb = rcParams['verbose']
        return

        print out 

    def _print_values(self, dat, label='', timestamps=[]):
        d = dat['data']
        a = dat['axes']
        shp = d.getshape()
        out = ''
        for i in range(shp[3]):
            out += '%s [%s]:\n' % (a[3],timestamps[i])
            t = d[:,:,:,i]
            for j in range(shp[0]):
                if shp[0] > 1: out +=  ' %s[%d] ' % (a[0],j)
                for k in range(shp[1]):
                    if shp[1] > 1: out +=  ' %s[%d] ' % (a[1],k)
                    for l in range(shp[2]):
                        if shp[2] > 1: out +=  ' %s[%d] ' % (a[2],l)
                        out += '= %3.3f\n' % (t[j,k,l])
            out += "-"*80
            out += "\n"
        print "-"*80
        print " ", label
        print "-"*80
        print out 

    def history(self):
        hist = list(self._gethistory())
        print "-"*80
        for h in hist:
            if h.startswith("---"):
                print h
            else:
                items = h.split("##")
                date = items[0]
                func = items[1]
                items = items[2:]
                print date            
                print "Function: %s\n  Parameters:" % (func)
                for i in items:
                    s = i.split("=")
                    print "   %s = %s" % (s[0],s[1])
                print "-"*80
        return

    #
    # Maths business
    #

    def average_time(self, mask=None, scanav=False, weight='tint'):
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
        Example:
            # time average the scantable without using a mask
            newscan = scan.average_time()            
        """
        varlist = vars()
        if weight is None: weight = 'tint'
        if mask is None: mask = ()
        from asap._asap import average as _av        
        s = scantable(_av((self,), mask, scanav, weight))
        s._add_history("average_time",varlist)
        return s
        
    def convert_flux(self, jyperk=None, eta=None, d=None, insitu=None,
                     allaxes=None):
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
        if allaxes is None: allaxes = rcParams['scantable.allaxes']
        if insitu is None: insitu = rcParams['insitu']
        varlist = vars()
        if jyperk is None: jyperk = -1.0
        if d is None: d = -1.0
        if eta is None: eta = -1.0
        if not insitu:
            from asap._asap import convertflux as _convert
            s = scantable(_convert(self, d, eta, jyperk, allaxes))
            s._add_history("convert_flux", varlist)
            return s
        else:
            from asap._asap import convertflux_insitu as _convert
            _convert(self, d, eta, jyperk, allaxes)
            self._add_history("convert_flux", varlist)
            return

    def gain_el(self, poly=None, filename="", method="linear",
                insitu=None, allaxes=None):
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
            allaxes:     If True apply to all spectra. Otherwise
                         apply only to the selected (beam/pol/if) spectra only
                         The default is taken from .asaprc (True if none)
        """

        if allaxes is None: allaxes = rcParams['scantable.allaxes']
        if insitu is None: insitu = rcParams['insitu']
        varlist = vars()
        if poly is None:
           poly = ()
        from os.path import expandvars
        filename = expandvars(filename)
        if not insitu:
            from asap._asap import gainel as _gainEl
            s = scantable(_gainEl(self, poly, filename, method, allaxes))
            s._add_history("gain_el", varlist)
            return s
        else:
            from asap._asap import gainel_insitu as _gainEl
            _gainEl(self, poly, filename, method, allaxes)
            self._add_history("gain_el", varlist)
            return
    
    def freq_align(self, reftime=None, method='cubic', perif=False,
                   insitu=None):
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
            perif:       Generate aligners per freqID (no doppler tracking) or
                         per IF (scan-based doppler tracking)
            insitu:      if False a new scantable is returned.
                         Otherwise, the scaling is done in-situ
                         The default is taken from .asaprc (False)
        """
        if insitu is None: insitu = rcParams['insitu']
        varlist = vars()
        if reftime is None: reftime = ''
        perfreqid = not perif
        if not insitu:
            from asap._asap import freq_align as _align
            s = scantable(_align(self, reftime, method, perfreqid))
            s._add_history("freq_align", varlist)
            return s
        else:
            from asap._asap import freq_align_insitu as _align
            _align(self, reftime, method, perfreqid)
            self._add_history("freq_align", varlist)
            return

    def opacity(self, tau, insitu=None, allaxes=None):
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
            allaxes:     if True apply to all spectra. Otherwise
                         apply only to the selected (beam/pol/if)spectra only
                         The default is taken from .asaprc (True if none)
        """
        if allaxes is None: allaxes = rcParams['scantable.allaxes']
        if insitu is None: insitu = rcParams['insitu']
        varlist = vars()
        if not insitu:
            from asap._asap import opacity as _opacity
            s = scantable(_opacity(self, tau, allaxes))
            s._add_history("opacity", varlist)
            return s
        else:
            from asap._asap import opacity_insitu as _opacity
            _opacity(self, tau, allaxes)
            self._add_history("opacity", varlist)
            return

    def bin(self, width=5, insitu=None):
        """
        Return a scan where all spectra have been binned up.
            width:       The bin width (default=5) in pixels
            insitu:      if False a new scantable is returned.
                         Otherwise, the scaling is done in-situ
                         The default is taken from .asaprc (False)
        """
        if insitu is None: insitu = rcParams['insitu']
        varlist = vars()
        if not insitu:
            from asap._asap import bin as _bin
            s = scantable(_bin(self, width))
            s._add_history("bin",varlist)
            return s
        else:
            from asap._asap import bin_insitu as _bin
            _bin(self, width)
            self._add_history("bin",varlist)
            return

    
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
        varlist = vars()
        if not insitu:
            from asap._asap import resample as _resample
            s = scantable(_resample(self, method, width))
            s._add_history("resample",varlist)
            return s
        else:
            from asap._asap import resample_insitu as _resample
            _resample(self, method, width)
            self._add_history("resample",varlist)
            return

    def average_pol(self, mask=None, weight='none', insitu=None):
        """
        Average the Polarisations together.
        The polarisation cursor of the output scan is set to 0
        Parameters:
            scan:        The scantable
            mask:        An optional mask defining the region, where the
                         averaging will be applied. The output will have all
                         specified points masked.
            weight:      Weighting scheme. 'none' (default), 'var' (1/var(spec)
                         weighted), or 'tsys' (1/Tsys**2 weighted)
            insitu:      if False a new scantable is returned.
                         Otherwise, the scaling is done in-situ
                         The default is taken from .asaprc (False)
        Example:
            polav = average_pols(myscan)
        """
        if insitu is None: insitu = rcParams['insitu']
        varlist = vars()

        if mask is None:
            mask = ()
        if not insitu:
            from asap._asap import averagepol as _avpol
            s = scantable(_avpol(self, mask, weight))
            s._add_history("average_pol",varlist)
            return s
        else:
            from asap._asap import averagepol_insitu as _avpol
            _avpol(self, mask, weight)
            self._add_history("average_pol",varlist)
            return

    def smooth(self, kernel="hanning", width=5.0, insitu=None, allaxes=None):
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
        varlist = vars()
        if not insitu:
            from asap._asap import smooth as _smooth
            s = scantable(_smooth(self,kernel,width,allaxes))
            s._add_history("smooth", varlist)
            return s
        else:
            from asap._asap import smooth_insitu as _smooth
            _smooth(self,kernel,width,allaxes)
            self._add_history("smooth", varlist)
            return

    def poly_baseline(self, mask=None, order=0, insitu=None):
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
            bscan = scan.poly_baseline(order=3)
        """
        if insitu is None: insitu = rcParams['insitu']
        varlist = vars()
        if mask is None:
            from numarray import ones
            mask = list(ones(scan.nchan()))
        from asap.asapfitter import fitter
        f = fitter()
        f._verbose(True)
        f.set_scan(self, mask)
        f.set_function(poly=order)
        sf = f.auto_fit(insitu)
        if insitu:
            self._add_history("poly_baseline", varlist)
            return
        else:
            sf._add_history("poly_baseline", varlist)
            return sf

    def auto_poly_baseline(self, mask=None, edge=(0,0), order=0,
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
            scan2=scan.auto_poly_baseline(order=7)
        """
        if insitu is None: insitu = rcParams['insitu']
        varlist = vars()
        from asap.asapfitter import fitter
        from asap.asaplinefind import linefinder
        from asap import _is_sequence_or_number as _is_valid
        
        if not _is_valid(edge, int):
            raise RuntimeError, "Parameter 'edge' has to be an integer or a \
            pair of integers specified as a tuple"
        
        # setup fitter
        f = fitter()
        f._verbose(True)
        f.set_function(poly=order)

        # setup line finder
        fl=linefinder()
        fl.set_options(threshold=threshold)

        if not insitu:
            workscan=self.copy()
        else:
            workscan=self

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
                    workscan.setpol(k)
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
                       x=f.get_parameters()
                       workscan._setspectrum(f.fitter.getresidual(),iRow)
        workscan.set_cursor(sel[0],sel[1],sel[2])
        workscan._vb = vb
        if not insitu:
           return workscan

    def rotate_linpolphase(self, angle, allaxes=None):
        """
        Rotate the phase of the complex polarization O=Q+iU correlation.  
        This is always done in situ in the raw data.  So if you call this 
        function more than once then each call rotates the phase further.
        Parameters:
            angle:   The angle (degrees) to rotate (add) by.
            allaxes: If True apply to all spectra. Otherwise
                     apply only to the selected (beam/pol/if)spectra only.
                     The default is taken from .asaprc (True if none)
        Examples:
            scan.rotate_linpolphase(2.3)
    """
        if allaxes is None: allaxes = rcParams['scantable.allaxes']
        varlist = vars()
        from asap._asap import _rotate_linpolphase as _rotate
        _rotate(self, angle, allaxes)
        self._add_history("rotate_linpolphase", varlist)
        return
    

    def rotate_xyphase(self, angle, allaxes=None):
        """
        Rotate the phase of the XY correlation.  This is always done in situ
        in the data.  So if you call this function more than once
        then each call rotates the phase further.
        Parameters:
            angle:   The angle (degrees) to rotate (add) by.
            allaxes: If True apply to all spectra. Otherwise
                     apply only to the selected (beam/pol/if)spectra only.
                     The default is taken from .asaprc (True if none)
        Examples:
            scan.rotate_xyphase(2.3)
        """
        if allaxes is None: allaxes = rcParams['scantable.allaxes']
        varlist = vars()
        from asap._asap import rotate_xyphase as _rotate_xyphase
        _rotate_xyphase(self, angle, allaxes)
        self._add_history("rotate_xyphase", varlist)
        return


    def add(self, offset, insitu=None, allaxes=None):
        """
        Return a scan where all spectra have the offset added
        Parameters:
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
        varlist = vars()
        if not insitu:
            from asap._asap import add as _add
            s = scantable(_add(self, offset, allaxes))
            s._add_history("add",varlist)
            return s
        else:
            from asap._asap import add_insitu as _add
            _add(self, offset, allaxes)
            self._add_history("add",varlist)
            return

    def scale(self, factor, insitu=None, allaxes=None, tsys=True):
        """
        Return a scan where all spectra are scaled by the give 'factor'
        Parameters:
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
        varlist = vars()
        if not insitu:
            from asap._asap import scale as _scale
            s = scantable(_scale(self, factor, allaxes, tsys))
            s._add_history("scale",varlist)
            return s
        else:
            from asap._asap import scale_insitu as _scale
            _scale(self, factor, allaxes, tsys)
            self._add_history("scale",varlist)
            return


    def quotient(self, other, isreference=True, preserve=True):
        """
        Return the quotient of a 'source' (on) scan and a 'reference' (off)
        scan.
        The reference can have just one row, even if the signal has many.
        Otherwise they must have the same number of rows.
        The cursor of the output scan is set to 0
        Parameters:
            other:          the 'other' scan
            isreference:    if the 'other' scan is the reference (default)
                            or source
            preserve:       you can preserve (default) the continuum or 
                            remove it.  The equations used are 
                            preserve: Output = Toff * (on/off) - Toff
                            remove:   Output = Tref * (on/off) - Ton
        Example:
            # src is a scantable for an 'on' scan, ref for an 'off' scantable
            q1 = src.quotient(ref)
            q2 = ref.quotient(src, isreference=False)
            # gives the same result
        """
        from asap._asap import quotient as _quot
        if isreference:
            return scantable(_quot(self, other, preserve))
        else:
            return scantable(_quot(other, self, preserve))

    def __add__(self, other):
        varlist = vars()
        s = None
        if isinstance(other, scantable):
            from asap._asap import b_operate as _bop            
            s = scantable(_bop(self, other, 'add', True))
        elif isinstance(other, float):
            from asap._asap import add as _add
            s = scantable(_add(self, other, True))
        else:
            print "Other input is not a scantable or float value"
            return
        s._add_history("operator +", varlist)
        return s


    def __sub__(self, other):
        """
        implicit on all axes and on Tsys
        """
        varlist = vars()
        s = None
        if isinstance(other, scantable):
            from asap._asap import b_operate as _bop            
            s = scantable(_bop(self, other, 'sub', True))
        elif isinstance(other, float):
            from asap._asap import add as _add
            s = scantable(_add(self, -other, True))
        else:
            print "Other input is not a scantable or float value"
            return
        s._add_history("operator -", varlist)
        return s
    
    def __mul__(self, other):
        """
        implicit on all axes and on Tsys
        """
        varlist = vars()
        s = None
        if isinstance(other, scantable):
            from asap._asap import b_operate as _bop            
            s = scantable(_bop(self, other, 'mul', True))
        elif isinstance(other, float):
            if other == 0.0:
                print "Multiplying by zero is not recommended."
                return
            from asap._asap import scale as _sca
            s = scantable(_sca(self, other, True))
        else:
            print "Other input is not a scantable or float value"
            return
        s._add_history("operator *", varlist)
        return s
    

    def __div__(self, other):
        """
        implicit on all axes and on Tsys
        """
        varlist = vars()
        s = None
        if isinstance(other, scantable):
            from asap._asap import b_operate as _bop            
            s = scantable(_bop(self, other, 'div', True))
        elif isinstance(other, float):
            if other == 0.0:
                print "Dividing by zero is not recommended"
                return
            from asap._asap import scale as _sca
            s = scantable(_sca(self, 1.0/other, True))
        else:
            print "Other input is not a scantable or float value"
            return
        s._add_history("operator /", varlist)
        return s

    def get_fit(self, row=0):
        """
        Print or return the stored fits for a row in the scantable
        Parameters:
            row:    the row which the fit has been applied to.
        """
        if row > self.nrow():
            return
        from asap import asapfit
        fit = asapfit(self._getfit(row))
        if self._vb:
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
