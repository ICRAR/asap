from asap._asap import sdtable
from asap import rcParams
from numarray import ones,zeros
import sys

class scantable(sdtable):
    """
        The ASAP container for scans
    """
    
    def __init__(self, filename, unit=None):
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
           unit:         brightness unit; must be consistent with K or Jy.
                         Over-rides the default selected by the reader
                         (input rpfits/sdfits/ms) or replaces the value
                         in existing scantables
        """
        self._vb = rcParams['verbose']
        self._p = None
        from os import stat as st
        import stat
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
                autoav = rcParams['scantable.autoaverage']

                from asap._asap import sdreader
                ifSel = -1
                beamSel = -1
                r = sdreader(filename,ifSel,beamSel)
                print 'Importing data...'
                r._read([-1])
                tbl = r._getdata()
                if unit is not None:
                    tbl.set_fluxunit(unit)
                if autoav:
                    from asap._asap import average
                    tmp = tuple([tbl])
                    print 'Auto averaging integrations...'
                    tbl2 = average(tmp,(),True,'none')
                    sdtable.__init__(self,tbl2)
                    del r,tbl
                else:
                    sdtable.__init__(self,tbl)

    def save(self, name=None, format=None, overwrite=False):
        """
        Store the scantable on disk. This can be an asap (aips++) Table, SDFITS, 
        Image FITS or MS2 format.
        Parameters:
            name:        the name of the outputfile. For format="FITS" this
                         is the directory file name into which all the files will
                         be written (default is 'asap_FITS')
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
            print "No filename given. Using default name %s..." % name
        name = path.expandvars(name)
        if path.isfile(name) or path.isdir(name):
            if not overwrite:
                print "File %s already exists." % name
                return
        if format == 'ASAP':
            self._save(name)
        else:
            from asap._asap import sdwriter as _sw
            w = _sw(format)
            w.write(self, name)
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
            scanid:    a scanno or a source name
        Example:
            scan.get_scan('323p459')
            # gets all scans containing the source '323p459'
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
            from os.path import expandvars
            filename = expandvars(filename)
            data = open(filename, 'w')
            data.write(info)
            data.close()
        print info

    def set_cursor(self, thebeam=0,theif=0,thepol=0):
        """
        Set the spectrum for individual operations.
        Parameters:
            thebeam,theif,thepol:    a number
        Example:
            scan.set_cursor(0,0,1)
            pol1sig = scan.stats(all=False) # returns std dev for beam=0
                                            # if=0, pol=1
        """
        self.setbeam(thebeam)
        self.setpol(thepol)
        self.setif(theif)
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

    def stats(self, stat='stddev', mask=None, all=None):
        """
        Determine the specified statistic of the current beam/if/pol
        Takes a 'mask' as an optional parameter to specify which
        channels should be excluded.
        Parameters:
            stat:    'min', 'max', 'sumsq', 'sum', 'mean'
                     'var', 'stddev', 'avdev', 'rms', 'median'
            mask:    an optional mask specifying where the statistic
                     should be determined.
            all:     if true show all (default or .asaprc) rather
                     that the cursor selected spectrum of Beam/IF/Pol

        Example:
            scan.set_unit('channel')
            msk = scan.create_mask([100,200],[500,600])
            scan.stats(stat='mean', mask=m)
        """
        if all is None: all = rcParams['scantable.allaxes']
        from asap._asap import stats as _stats
        from numarray import array,zeros,Float
        if mask == None:
            mask = ones(self.nchan())
        axes = ['Beam','IF','Pol','Time']

        if all:
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
            return retval

        else:
            i,j,k = (self.getbeam(),self.getif(),self.getpol())
            statval = _stats(self,mask,stat,-1)
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
                print " ",stat
                print "--------------------------------------------------"
                print out
            retval = {'axes': axes, 'data': array(statval), 'cursor':(i,j,k)}
            return retval

    def stddev(self,mask=None, all=None):
        """
        Determine the standard deviation of the current beam/if/pol
        Takes a 'mask' as an optional parameter to specify which
        channels should be excluded.
        Parameters:
            mask:    an optional mask specifying where the standard
                     deviation should be determined.
            all:     optional flag to show all or a cursor selected
                     spectrum of Beam/IF/Pol. Default is all or taken
                     from .asaprc

        Example:
            scan.set_unit('channel')
            msk = scan.create_mask([100,200],[500,600])
            scan.stddev(mask=m)
        """
        if all is None: all = rcParams['scantable.allaxes']
        return self.stats(stat='stddev',mask=mask, all=all);

    def get_tsys(self, all=None):
        """
        Return the System temperatures.
        Parameters:
            all:    optional parameter to get the Tsys values for all
                    Beams/IFs/Pols (default) or just the one selected
                    with scantable.set_cursor()
                    [True or False]
        Returns:
            a list of Tsys values.
        """
        if all is None: all = rcParams['scantable.allaxes']
        from numarray import array,zeros,Float
        axes = ['Beam','IF','Pol','Time']

        if all:
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

        if unit in ['','pixel', 'channel']:
            unit = ''
        inf = list(self._getcoordinfo())
        inf[0] = unit
        self._setcoordinfo(inf)
        if self._p: self.plot()

    def set_instrument (self, instr):
        """
        Set the instrument for subsequent processing
        Parameters:
            instr:    Select from 'ATPKSMB', 'ATPKSHOH', 'ATMOPRA', 
                      'DSS-43' (Tid), 'CEDUNA', and 'HOBART'
        """
        self._setInstrument(instr)

    def set_doppler(self, doppler='RADIO'):
        """
        Set the doppler for all following operations on this scantable.
        Parameters:
            doppler:    One of 'RADIO', 'OPTICAL', 'Z', 'BETA', 'GAMMA'
        """

        inf = list(self._getcoordinfo())
        inf[2] = doppler
        self._setcoordinfo(inf)
        if self._p: self.plot()

    def set_freqframe(self, frame=None):
        """
        Set the frame type of the Spectral Axis.
        Parameters:
            frame:   an optional frame type, default 'LSRK'.
        Examples:
            scan.set_freqframe('BARY')
        """
        if not frame: frame = rcParams['scantable.freqframe']
        valid = ['REST','TOPO','LSRD','LSRK','BARY', \
                   'GEO','GALACTO','LGROUP','CMB']
        if frame in valid:
            inf = list(self._getcoordinfo())
            inf[1] = frame
            self._setcoordinfo(inf)
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
        u = self._getcoordinfo()[0]
        if self._vb:
            if u == "": u = "channel"
            print "The current mask window unit is", u
        n = self.nchan()
        data = self._getabcissa()
        msk = zeros(n)
        for  window in args:
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

    def set_restfreqs(self, freqs=None, unit='Hz', lines=None, source=None, theif=None):
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
        return


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
            out += "--------------------------------------------------\n"
        print "--------------------------------------------------"
        print " ", label
        print "--------------------------------------------------"
        print out 
