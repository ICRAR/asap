import numpy
from asap.scantable import scantable
from asap._asap import stgrid
import pylab as pl

class asapgrid:
    def __init__( self, infile ):
        self.infile = infile
        self.outfile = None
        self.gridder = stgrid( self.infile )

    def setData( self, infile ):
        self.gridder._setin( infile )

    def defineImage( self, nx=-1, ny=-1, cellx='', celly='', center='' ):
        self.gridder._defineimage( nx, ny, cellx, celly, center )

    def setOption( self, convType='box', convSupport=-1 ):
        self.gridder._setoption( convType, convSupport )

    def grid( self ):
        self.gridder._grid()

    def save( self, outfile='' ):
        self.outfile = self.gridder._save( outfile ) 

    def plot( self, plotchan=-1 ):
        plotter = _SDGridPlotter( self.infile, self.outfile )
        plotter.plot( chan=plotchan )
        
class _SDGridPlotter:
    def __init__( self, infile, outfile=None ):
        self.infile = infile
        self.outfile = outfile
        if self.outfile is None:
            self.outfile = self.infile.rstrip('/')+'.grid'
        self.grid = None
        self.pointing = None
        self.data = None
        self.nx = -1
        self.ny = -1
        self.nchan = 0
        self.cellx = 0.0
        self.celly = 0.0
        self.center = [0.0,0.0]
        self.nonzero = [[0.0],[0.0]]
        self.get()

    def get( self ):
        s = scantable( self.infile, average=False )
        self.pointing = numpy.array( s.get_directionval() ).transpose()
        spectra = []
        for i in xrange(s.nrow()):
            spectra.append( s._getspectrum( i ) )
        spectra = numpy.array( spectra ).transpose()
        self.nchan = spectra.shape[0]
        del s

        idx = spectra.nonzero()[1]
        #self.nonzero = self.pointing.take( idx, axis=1 )
        
        s = scantable( self.outfile, average=False )
        self.grid = numpy.array( s.get_directionval() ).transpose()
        dirstring = s.get_direction()
        nrow = s.nrow()
        spectra = []
        for i in xrange(nrow):
            spectra.append( s._getspectrum( i ) )
        spectra = numpy.array( spectra ).transpose()

        idx = 0
        d0 = dirstring[0].split()[-1]
        while ( dirstring[idx].split()[-1] == d0 ):  
            idx += 1

        self.ny = idx
        self.nx = nrow / idx
        
        self.cellx = abs( self.grid[0][0] - self.grid[0][1] )
        self.celly = abs( self.grid[1][0] - self.grid[1][self.ny] )

        self.data = spectra.reshape( (self.nchan,self.nx,self.ny) )

    def plot( self, chan=-1 ):
        if chan < 0:
            data = self.data.mean(axis=0)
            title = 'Gridded Image (averaged over channel)'
        else:
            data = self.data[chan]
            title = 'Gridded Image (channel %s)'%(chan)
        pl.figure(10)
        pl.clf()
        pl.plot(self.grid[0],self.grid[1],'.',color='blue')
        pl.plot(self.pointing[0],self.pointing[1],'.',color='red')
        #pl.plot(self.nonzero[0],self.nonzero[1],'o',color='green')
        extent=[self.grid[0].min()-0.5*self.cellx,
                self.grid[0].max()+0.5*self.cellx,
                self.grid[1].min()-0.5*self.celly,
                self.grid[1].max()+0.5*self.celly]
        pl.imshow(data,extent=extent,origin='lower',interpolation='nearest')
        pl.colorbar()
        pl.xlabel('R.A. [rad]')
        pl.ylabel('Dec. [rad]')
        pl.title( title )
