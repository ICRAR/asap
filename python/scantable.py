from asap._asap import sdtable
from numarray import ones
import sys

class scantable(sdtable):
    """
        The ASAP container for scans
    """
    
    def __init__(self, filename):
        """
        Create a scantable from a saved one or make a reference
        Parameters:
            filename:    the name of an asap table on disk, or
                         [advanced] a refernce to an existing
                         scantable
        """
        self._vb = True
        self._p = None
        sdtable.__init__(self, filename)

    def save(self, name, format='ASAP'):
        """
        Store the scantable on disk. This can be a asap file or SDFITS/MS2.
        Parameters:
            name:        the name of the outputfile
            format:      an optional file format. Default is ASAP.
                         Alllowed are 'ASAP', 'SDFITS' and 'MS2'
        Example:
            scan.save('myscan.asap')
            scan.save('myscan.sdfits','SDFITS')
        """
        if format == 'ASAP':
            self._save(name)
        else:
            from asap._asap import sdwriter as _sw
            w = _sw()
            if format == 'SDFITS':
                w.write(self, name)
        return

    def _verbose(self, *args):
        """
        Set the verbose level True or False, to indicate if output
        should be printed as well as returned.
        """
        if type(args[0]) is bool:
            self._vb = args[0]
            return
        elif len(args) == 0:
            return self._vb


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
                s = sdtable._getscan(self,scanid)
                return scantable(s)
        except RuntimeError:
            print "Couldn't find any match."

    def __str__(self):
        return sdtable.summary(self)

    def summary(self,filename=None):
        """
        Print a summary of the contents of this scantable.
        Parameters:
            filename:    the name of a file to write the putput to
                         Default - no file output
        """
        info = sdtable.summary(self)
        if filename is not None:
            data = open(filename, 'w')
            data.write(info)
            data.close()
        print info

    def set_selection(self, thebeam=0,theif=0,thepol=0):
        """
        Set the spectrum for individual operations.
        Parameters:
            thebeam,theif,thepol:    a number
        Example:
            scan.set_selection(0,0,1)
            pol1rms = scan.rms(all=False) # returns rms for beam=0
                                         # if=0, pol=1
        """
        self.setbeam(thebeam)
        self.setpol(thepol)
        self.setif(theif)
        return

    def get_selection(self):
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
            out = 'Beam=%d IF=%d Pol=%d '% (i,j,k)
            print out
        return i,j,k

    def rms(self,mask=None, all=True):
        """
        Determine the root mean square of the current beam/if/pol
        Takes a 'mask' as an optional parameter to specify which
        channels should be excluded.
        Parameters:
            mask:    an optional mask specifying where the rms
                     should be determined.
            all:     optional flag to show all or a selected
                     spectrum of Beam/IF/Pol

        Example:
            scan.set_unit('channel')
            msk = scan.create_mask([100,200],[500,600])
            scan.rms(mask=m)
        """
        from asap._asap import rms as _rms
        if mask == None:
            mask = ones(self.nchan())
        if all:
            out = ''
            tmp = []
            for i in range(self.nbeam()):
                self.setbeam(i)
                for j in range(self.nif()):
                    self.setif(j)
                    for k in range(self.npol()):
                        self.setpol(k)
                        rmsval = _rms(self,mask)
                        tmp.append(rmsval)
                        out += 'Beam[%d], IF[%d], Pol[%d] = %3.3f\n' % (i,j,k,rmsval)
            if self._vb:
                print out
            return tmp

        else:
            i = self.getbeam()
            j = self.getif()
            k = self.getpol()
            rmsval = _rms(self,mask)
            out = 'Beam[%d], IF[%d], Pol[%d] = %3.3f' % (i,j,k,rmsval)
            if self._vb:
                print out
            return rmsval

    def get_tsys(self, all=True):
        """
        Return the System temperatures.
        Parameters:
            all:    optional parameter to get the Tsys values for all
                    Beams/IFs/Pols (default) or just the one selected
                    with scantable.set_selection()
                    [True or False]
        Returns:
            a list of Tsys values.
        """
        if all:
            tmp = []
            out = ''
            for i in range(self.nbeam()):
                self.setbeam(i)
                for j in range(self.nif()):
                    self.setif(j)
                    for k in range(self.npol()):
                        self.setpol(k)
                        ts = self._gettsys()
                        tmp.append(ts)
                        out += 'TSys: Beam[%d], IF[%d], Pol[%d] = %3.3f\n' % (i,j,k,ts)
            if self._vb:
                print out
            return tmp
        else:
            i = self.getbeam()
            j = self.getif()
            k = self.getpol()
            ts = self._gettsys()
            out = 'TSys: Beam[%d], IF[%d], Pol[%d] = %3.3f' % (i,j,k,ts)
            if self._vb:
                print out
            return ts
        
    def get_time(self):
        """
        Get a list of time stamps for the observations.
        Return a string for each intergration in the scantable.
        Parameters:
            none
        Example:
            none
        """
        out = []
        for i in range(self.nrow()):
            out.append(self._gettime(i))
        return out

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

    def set_freqframe(self, frame='LSRK'):
        """
        Set the frame type of the Spectral Axis.
        Parameters:
            frame:   an optional frame type, default 'LSRK'.
        Examples:
            scan.set_freqframe('BARY')
        """
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

    def get_abscissa(self, rowno=0):
        """
        Get the abscissa in the current coordinate setup for the currently
        selected Beam/IF/Pol
        Parameters:
            none
        Returns:
            The abscissa values and it's format string.
        """
        absc = self.getabscissa(rowno)
        lbl = self.getabscissalabel(rowno)
        return absc, lbl

    def create_mask(self, *args, **kwargs):
        """
        Compute and return a mask based on [min,max] windows.
        The specified windows are to be EXCLUDED, when the mask is
        applied.
        Parameters:
            [min,max],[min2,max2],...
                Pairs of start/end points specifying the regions
                to be masked
            invert:     return an inverted mask, i.e. the regions
                        specified are not masked (INCLUDED)
        Example:
            scan.set_unit('channel')

            a)
            msk = scan.set_mask([400,500],[800,900])
            # masks the regions between 400 and 500
            # and 800 and 900 in the unit 'channel'

            b)
            msk = scan.set_mask([400,500],[800,900], invert=True)
            # masks the regions outside 400 and 500
            # and 800 and 900 in the unit 'channel'
           
        """
        u = self._getcoordinfo()[0]
        if self._vb:
            if u == "": u = "channel"
            print "The current mask window unit is", u
        n = self.nchan()
        data = self.getabscissa()
        msk = ones(n)
        for  window in args:
            if (len(window) != 2 or window[0] > window[1] ):
                print "A window needs to be defined as [min,max]"
                return
            for i in range(n):
                if data[i] >= window[0] and data[i] < window[1]:
                    msk[i] = 0
        if kwargs.has_key('invert'):
            if kwargs.get('invert'):
                from numarray import logical_not
                msk = logical_not(msk)
        return msk
    
    def set_restfreqs(self, freqs, unit='Hz'):
        """
        Set the restfrequency(s) for this scantable.
        Parameters:
            freqs:    one or more frequencies
            unit:     optional 'unit', default 'Hz'
        Example:
            scan.set_restfreqs([1000000000.0])
        """
        if type(freqs) is float or int:
            freqs = (freqs)
        sdtable._setrestfreqs(self,freqs, unit)
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
        if (thebeam < self.nbeam() and  theif < self.nif() and thepol < self.npol()):
            stable.setbeam(thebeam)
            stable.setif(theif)
            stable.setpol(thepol)
            stable._flag(self)
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
        validcol = {'Beam':self.nbeam(),'IF':self.nif(),'Pol':self.npol()}

        validyax = ['spectrum','tsys']
        if not self._p:
            from asap.asaplot import ASAPlot
            self._p = ASAPlot()
#            print "Plotting not enabled"
#            return
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
        vb = self._verbose()
        self._verbose(False)
        sel = self.get_selection()
        for i in range(npan):
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
                    x,xlab = self.get_abscissa(i)
                    y = self.getspectrum(i)
                    ylab = r'Flux'
                    m = self.getmask(i)
                llab = col+' '+str(j)
                self._p.set_line(label=llab)
                self._p.plot(x,y,m)
            self._p.set_axes('xlabel',xlab)
            self._p.set_axes('ylabel',ylab)
            self._p.set_axes('title',tlab)
        self._p.release()
        self.set_selection(sel[0],sel[1],sel[2])
        self._verbose(vb)
        return
