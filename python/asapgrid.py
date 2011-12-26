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

    def plot( self, plotchan=-1, plotpol=-1, plotobs=False, plotgrid=False ):
        import time
        t0=time.time()
        # to load scantable on disk
        storg = rcParams['scantable.storage']
        rcParams['scantable.storage'] = 'disk'
        plotter = _SDGridPlotter( self.infile, self.outfile, self.ifno )
        plotter.plot( chan=plotchan, pol=plotpol, plotobs=plotobs, plotgrid=plotgrid )
        # back to original setup
        rcParams['scantable.storage'] = storg
        t1=time.time()
        asaplog.push('plot: elapsed time %s sec'%(t1-t0))
        asaplog.post('DEBUG','asapgrid.plot')
        
class _SDGridPlotter:
    def __init__( self, infile, outfile=None, ifno=-1 ):
        self.infile = infile
        self.outfile = outfile
        if self.outfile is None:
            self.outfile = self.infile.rstrip('/')+'.grid'
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
        self.tablein = None
        self.nrow = 0
        self.blc = None
        self.trc = None
        self.get()

    def get( self ):
        s = scantable( self.outfile, average=False )
        self.nchan = len(s._getspectrum(0))
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

        idx = 0
        d0 = s.get_direction( 0 ).split()[-1]
        while ( s.get_direction(self.npol*idx).split()[-1] == d0 ):  
            idx += 1
        
        self.nx = idx
        self.ny = nrow / (self.npol * idx )
        #print 'nx,ny=',self.nx,self.ny

        self.blc = s.get_directionval( 0 )
        self.trc = s.get_directionval( nrow-self.npol )
        #print self.blc
        #print self.trc
        incrx = s.get_directionval( self.npol )
        incry = s.get_directionval( self.nx*self.npol ) 
        self.cellx = abs( self.blc[0] - incrx[0] )
        self.celly = abs( self.blc[1] - incry[1] )
        #print 'cellx,celly=',self.cellx,self.celly

    def plot( self, chan=-1, pol=-1, plotobs=False, plotgrid=False ):
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
        # plot grid position
        if plotgrid:
            x = numpy.arange(self.blc[0],self.trc[0]+0.5*self.cellx,self.cellx,dtype=float)
            #print 'len(x)=',len(x)
            #print 'x=',x
            ybase = numpy.ones(self.nx,dtype=float)*self.blc[1]
            #print 'len(ybase)=',len(ybase)
            incr = self.celly 
            for iy in xrange(self.ny):
                y = ybase + iy * incr
                #print y
                pl.plot(x,y,',',color='blue')
        # plot observed position
        if plotobs:
            self.createTableIn()
            irow = 0 
            while ( irow < self.nrow ):
                chunk = self.getPointingChunk( irow )
                #print chunk
                pl.plot(chunk[0],chunk[1],',',color='green')
                irow += chunk.shape[1]
                #print irow
        # show image
        extent=[self.blc[0]-0.5*self.cellx,
                self.trc[0]+0.5*self.cellx,
                self.blc[1]-0.5*self.celly,
                self.trc[1]+0.5*self.celly]
        pl.imshow(data,extent=extent,origin='lower',interpolation='nearest')
        pl.colorbar()
        pl.xlabel('R.A. [rad]')
        pl.ylabel('Dec. [rad]')
        pl.title( title )

    def createTableIn( self ):
        self.tablein = scantable( self.infile, average=False )
        if self.ifno < 0:
            ifno = self.tablein.getif(0)
            print 'ifno=',ifno
        else:
            ifno = self.ifno
        sel = selector()
        sel.set_ifs( ifno )
        self.tablein.set_selection( sel ) 
        self.nchan = len(self.tablein._getspectrum(0))
        self.nrow = self.tablein.nrow() 
        del sel
        

    def getPointingChunk( self, irow ):
        numchunk = 1000
        nrow = min( self.nrow-irow, numchunk )
        #print 'nrow=',nrow
        v = numpy.zeros( (2,nrow), dtype=float )
        idx = 0
        for i in xrange(irow,irow+nrow):
            d = self.tablein.get_directionval( i )
            v[0,idx] = d[0]
            v[1,idx] = d[1]
            idx += 1
        return v

    def getData( self, chan=-1, pol=-1 ):
        if chan == -1:
            spectra = self.__chanAverage()
        else:
            spectra = self.__chanIndex( chan )
        data = spectra.reshape( (self.npol,self.ny,self.nx) )
        if pol == -1:
            retval = data.mean(axis=0)
        else:
            retval = data[pol]
        #retval[0][self.nx-1] = -1.0
        return retval

    def __chanAverage( self ):
        s = scantable( self.outfile, average=False )
        nrow = s.nrow() 
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
        nrow = s.nrow()
        spectra = numpy.zeros( (self.npol,nrow/self.npol), dtype=float )
        irow = 0
        sp = [0 for i in xrange(self.nchan)]
        for i in xrange(nrow/self.npol):
            for ip in xrange(self.npol):
                sp = s._getspectrum( irow )
                spectra[ip,i] = sp[idx]
                irow += 1
        return spectra
        
            
