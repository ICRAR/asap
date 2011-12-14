import numpy
from asap import rcParams
from asap.scantable import scantable
from asap.selector import selector
from asap._asap import stgrid
import pylab as pl
from logging import asaplog

class asapgrid:
    def __init__( self, infile ):
        self.infile = infile
        self.outfile = None
        self.gridder = stgrid( self.infile )
        self.ifno = None

    def setData( self, infile ):
        self.gridder._setin( infile )

    def setIF( self, ifno ):
        self.ifno = ifno
        self.gridder._setif( self.ifno )

    def setPolList( self, pollist ):
        self.gridder._setpollist( pollist )

    def setScanList( self, scanlist ):
        self.gridder._setscanlist( scanlist )

    def defineImage( self, nx=-1, ny=-1, cellx='', celly='', center='' ):
        self.gridder._defineimage( nx, ny, cellx, celly, center )

    def setFunc( self, func='box', width=-1 ):
        self.gridder._setfunc( func, width )

    def setWeight( self, weightType='uniform' ):
        self.gridder._setweight( weightType ) 

    def grid( self ):
        self.gridder._grid()

    def save( self, outfile='' ):
        self.outfile = self.gridder._save( outfile ) 

    def plot( self, plotchan=-1, plotpol=-1 ):
        import time
        t0=time.time()
        # to load scantable on disk
        storg = rcParams['scantable.storage']
        rcParams['scantable.storage'] = 'disk'
        plotter = _SDGridPlotter( self.infile, self.outfile, self.ifno )
        plotter.plot( chan=plotchan, pol=plotpol )
        # back to original setup
        rcParams['scantable.storage'] = storg
        t1=time.time()
        asaplog.push('plot: elapsed time %s sec'%(t1-t0))
        asaplog.post('DEBUG','asapgrid.plot')
        
class _SDGridPlotter:
    def __init__( self, infile, outfile=None, ifno=0 ):
        self.infile = infile
        self.outfile = outfile
        if self.outfile is None:
            self.outfile = self.infile.rstrip('/')+'.grid'
        self.grid = None
        self.pointing = None
        self.nx = -1
        self.ny = -1
        self.nchan = 0
        self.npol = 0
        self.pollist = []
        self.cellx = 0.0
        self.celly = 0.0
        self.center = [0.0,0.0]
        self.nonzero = [[0.0],[0.0]]
        self.ifno = ifno
        self.get()

    def get( self ):
        s = scantable( self.infile, average=False )
        sel = selector()
        sel.set_ifs( self.ifno )
        s.set_selection( sel ) 
        self.pointing = numpy.array( s.get_directionval() ).transpose()
        self.nchan = len(s._getspectrum(0))
        s.set_selection()
        del s
        del sel

        s = scantable( self.outfile, average=False )
        nrow = s.nrow()
        pols = numpy.ones( nrow, dtype=int )
        for i in xrange(nrow):
            pols[i] = s.getpol(i)
        self.pollist, indices = numpy.unique( pols, return_inverse=True )
        self.npol = len(self.pollist)
        self.pollist = self.pollist[indices[:self.npol]]
        #print 'pollist=',self.pollist
        #print 'npol=',self.npol
        #print 'nrow=',nrow
        dirstring = numpy.array(s.get_direction()).take(range(0,nrow,self.npol))
        self.grid = numpy.array( s.get_directionval() ).take(range(0,nrow,self.npol),axis=0).transpose()

        idx = 0
        d0 = dirstring[0].split()[-1]
        while ( dirstring[idx].split()[-1] == d0 ):  
            idx += 1
        
        self.ny = idx
        self.nx = nrow / (self.npol * idx )
        #print 'nx,ny=',self.nx,self.ny
        
        self.cellx = abs( self.grid[0][0] - self.grid[0][1] )
        self.celly = abs( self.grid[1][0] - self.grid[1][self.ny] )
        #print 'cellx,celly=',self.cellx,self.celly

    def plot( self, chan=-1, pol=-1 ):
        if pol < 0:
            opt = 'averaged over pol'
        else:
            opt = 'pol %s'%(pol)
        if chan < 0:
            opt += ', averaged over channel'
        else:
            opt += ', channel %s'%(chan)
        data = self.getData( chan, pol ) 
        title = 'Gridded Image (%s)'%(opt)
        pl.figure(10)
        pl.clf()
        pl.plot(self.grid[0],self.grid[1],',',color='blue')
        pl.plot(self.pointing[0],self.pointing[1],',',color='green')
        extent=[self.grid[0].min()-0.5*self.cellx,
                self.grid[0].max()+0.5*self.cellx,
                self.grid[1].min()-0.5*self.celly,
                self.grid[1].max()+0.5*self.celly]
        pl.imshow(data,extent=extent,origin='lower',interpolation='nearest')
        pl.colorbar()
        pl.xlabel('R.A. [rad]')
        pl.ylabel('Dec. [rad]')
        pl.title( title )

    def getData( self, chan=-1, pol=-1 ):
        if chan == -1:
            spectra = self.__chanAverage()
        else:
            spectra = self.__chanIndex( chan )
        data = spectra.reshape( (self.npol,self.nx,self.ny) )
        if pol == -1:
            retval = data.mean(axis=0)
        else:
            retval = data[pol]
        return retval

    def __chanAverage( self ):
        s = scantable( self.outfile, average=False )
        nrow = self.nx * self.ny
        spectra = numpy.zeros( (self.npol,nrow/self.npol), dtype=float )
        irow = 0
        sp = [0 for i in xrange(self.nchan)]
        for i in xrange(nrow/self.npol):
            for ip in xrange(self.npol):
                sp = s._getspectrum( irow )
                spectra[ip,i] = numpy.mean( sp )
                irow += 1
        return spectra

    def __chanIndex( self, idx ):
        s = scantable( self.outfile, average=False )
        nrow = self.nx * self.ny
        spectra = numpy.zeros( (self.npol,nrow/self.npol), dtype=float )
        irow = 0
        sp = [0 for i in xrange(self.nchan)]
        for i in xrange(nrow/self.npol):
            for ip in xrange(self.npol):
                sp = s._getspectrum( irow )
                spectra[ip,i] = sp[idx]
                irow += 1
        return spectra
        
            
