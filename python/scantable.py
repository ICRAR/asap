from asap._asap import sdtable
from numarray import ones
import sys

class scantable(sdtable):
    """
        The ASAP container for scans
    """
    def functions(self):
        return ['copy','get_scan','summary','rms','set_unit','create_mask']

    def __init__(self, filename):
        """
        Create a scantable from a saved one or make a reference
        Parameters:
            filename:    the name of an asap table on disk, or
                         [advanced] a refernce to an existing
                         scantable
        """
        sdtable.__init__(self, filename)
        self._unit = 'channel'

    def copy(self):
        """
        Return a copy of this scantable.
        Parameters:

        Example:
            copiedscan = scan.copy()
        """
        sd = sdtable._copy(self)
        return scantable(sd)
    def get_scan(self, scanid=None):
        """
        Return a specific scan (by scanno) or collection of scans (by
        source name) in a new scantable.
        Parameters:
            scanid:    a scanno or a source name
        Example:
            scan.getscan('323p459')
            # gets all scans containing the source '323p459'
        """
        if scanid is None:
            print "Please specify a scan no or name to retrieve from the scantable"
        try:
            if type(scanid) is str:
                s = sdtable.getsource(self,scanid)
                return scantable(s)
            elif type(scanid) is int:
                s = sdtable.getscan(self,scanid)
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
        i = self.getbeam()
        j = self.getif()
        k = self.getpol()
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
            scan.setuint('channel')
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
            print out
            return tmp

        else:
            i = self.getbeam()
            j = self.getif()
            k = self.getpol()
            rmsval = _rms(self,mask)
            out = 'Beam[%d], IF[%d], Pol[%d] = %3.3f' % (i,j,k,rmsval)
            print out
            return rmsval

    def get_tsys(self, all=True):
        if all:
            tmp = []
            out = ''
            for i in range(self.nbeam()):
                self.setbeam(i)
                for j in range(self.nif()):
                    self.setif(j)
                    for k in range(self.npol()):
                        self.setpol(k)
                        ts = self.gettsys()
                        tmp.append(ts)
                        out += 'TSys: Beam[%d], IF[%d], Pol[%d] = %3.3f\n' % (i,j,k,ts)
            print out
            return tmp
        else:
            i = self.getbeam()
            j = self.getif()
            k = self.getpol()
            ts = self.gettsys()
            out = 'TSys: Beam[%d], IF[%d], Pol[%d] = %3.3f' % (i,j,k,ts)
            print out
            return ts

    def set_unit(self, unit='channel'):
        """
        Set the unit for all following operations on this scantable
        Parameters:
            unit:    optional unit, default is 'channel'
                     one of 'GHz','MHz','km/s','channel'
        """
        units = ['GHz','MHz','km/s','channel','pixel','']
        if unit in units:
            if unit in ['','pixel']:
                unit = 'channel'
            self._unit = unit
        else:
            print "Invalid unit given, please use one of:"
            print units

    def create_mask(self, *args):
        """
        Compute and return a mask based on [min,max] windows.
        Parameters:
            [min,max],[min2,max2],...
                Pairs of start/end points specifying the regions
                to be masked
        Example:
            scan.setunit('channel')
            msk = scan.set_mask([400,500],[800,900])
            masks the regions between 400 and 500
            and 800 and 900 in the unit 'channel'
        """
        print "The current mask window unit is", self._unit
        n = self.nchan()
        if self._unit == 'channel':
            data = range(n)
        else:
            data = self.getabscissa(unit=self._unit)
        msk = ones(n)
        for  window in args:
            if (len(window) != 2 or window[0] > window[1] ):
                print "A window needs to be defined as [min,max]"
                return
            for i in range(n):
                if data[i] > window[0] and data[i] < window[1]:
                    msk[i] = 0
        return msk
    def set_restfrequencies(self, freqs, unit='Hz'):
        """
        Set the restfrequency(s) for this scantable.
        Parameters:
            freqs:    one or more frequencies
            unit:     optional 'unit', default 'Hz'
        Example:
            scan.set_restfrequencies([1000000000.0])
        """
        self.setrestfreqs(freqs, unit)
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
