from asap.asaplot import ASAPlot
from asap import rcParams

class asapplotter:
    """
    The ASAP plotter.
    By default the plotter is set up to plot polarisations
    'colour stacked' and scantables across panels.
    The defaul plotter is called 'plotter'.
    Note:
        Currenly it only plots 'spectra' not Tsys or
        other variables.
    """
    def __init__(self):
        self._plotter = ASAPlot()

        self._tdict = {'Time':'t','time':'t','t':'t','T':'t'}
        self._bdict = {'Beam':'b','beam':'b','b':'b','B':'b'}
        self._idict = {'IF':'i','if':'i','i':'i','I':'i'}
        self._pdict = {'Pol':'p','pol':'p','p':'p'}
        self._sdict = {'scan':'s','Scan':'s','s':'s','S':'s'}
        self._cdict = {'t':'scan.nrow()',
                       'b':'scan.nbeam()',
                       'i':'scan.nif()',
                       'p':'scan.npol()',
                       's':'len(scans)'}
        self._ldict = {'b':'Beam',
                       'i':'IF',
                       'p':'Pol',
                       's':'Scan'}
        self._dicts = [self._tdict,self._bdict,
                       self._idict,self._pdict,
                       self._sdict]
        self._panels = 's'
        self._stacking = rcParams['plotter.stacking']
        self._autoplot = False
        self._minmax = None
        self._data = None
        self._lmap = []
        self._title = None

    def _translate(self, name):
        for d in self._dicts:
            if d.has_key(name):
                return d[name]
        return None
        
    def plot(self,*args):
        """
        Plot a (list of) scantables.
        Parameters:
            one or more comma separated scantables 
        Note:
            If a (list) of scantables was specified in a previous call
            to plot, no argument has to be given to 'replot'
            NO checking is done that the abscissas of the scantables
            are consistent e.g. all 'channel' or all 'velocity' etc.
        """
        if self._plotter.is_dead:
            self._plotter = ASAPlot()
        self._plotter.clear()
        self._plotter.hold()
        if len(args) > 0:
            self._data = tuple(args)            
        if self._panels == 't':
            if self._data[0].nrow() > 25:
                print "Scan to be plotted contains more than 25 rows.\nCan't plot that many panels..."
                return
            self._plot_time(self._data[0], self._stacking)
        elif self._panels == 's':
            self._plot_scans(self._data, self._stacking)
        else:
            self._plot_other(self._data, self._stacking)
        if self._minmax is not None:
            self._plotter.set_limits(xlim=self._minmax)
        self._plotter.release()
        return

    def _plot_time(self, scan, colmode):
        if colmode == 't':
            return
        n = scan.nrow()
        cdict = {'b':'scan.setbeam(j)',
                 'i':'scan.setif(j)',
                 'p':'scan.setpol(j)'}
        if self._stacking is not None:
            ncol = eval(self._cdict.get(colmode))
        self._plotter.set_panels()
        if n > 1:
            self._plotter.set_panels(rows=n)
        for i in range(n):
            if n > 1:
                self._plotter.palette(0)
                self._plotter.subplot(i)
            for j in range(ncol):
                eval(cdict.get(colmode))
                x = None
                y = None
                m = None
                if not self._title:
                    tlab = scan._getsourcename(i)                    
                else:
                    if len(self._title) == n:
                        tlab = self._title[i]
                    else:
                        tlab = scan._getsourcename(i)                   
                x,xlab = scan.get_abcissa(i)
                y = scan.getspectrum(i)
                ylab = 'Flux ('+scan.get_fluxunit()+')'
                m = scan.getmask(i)
                if self._lmap and len(self._lmap) > 0:
                    llab = self._lmap[j]
                else:
                    llab = self._ldict.get(colmode)+' '+str(j)
                self._plotter.set_line(label=llab)
                self._plotter.plot(x,y,m)
                xlim=[min(x),max(x)]
                self._plotter.axes.set_xlim(xlim)
            self._plotter.set_axes('xlabel',xlab)
            self._plotter.set_axes('ylabel',ylab)
            self._plotter.set_axes('title',tlab)            
        return

    def _plot_scans(self, scans, colmode):        
        if colmode == 's':
            return
        cdict = {'b':'scan.setbeam(j)',
                 'i':'scan.setif(j)',
                 'p':'scan.setpol(j)'}
        n = len(scans)
        if self._stacking is not None:
            scan = scans[0]
            ncol = eval(self._cdict.get(colmode))
        self._plotter.set_panels()
        if n > 1:
            self._plotter.set_panels(rows=n)
        i = 0
        for scan in scans:
            if n > 1:
                self._plotter.subplot(i)
                self._plotter.palette(0)
            for j in range(ncol):
                eval(cdict.get(colmode))
                x = None
                y = None
                m = None
                tlab = self._title
                if not self._title:
                    tlab = scan._getsourcename()
                x,xlab = scan.get_abcissa()
                y = scan.getspectrum()
                ylab = 'Flux ('+scan.get_fluxunit()+')'
                m = scan.getmask()
                if len(self._lmap) > 0:
                    llab = self._lmap[j]
                else:
                    llab = self._ldict.get(colmode)+' '+str(j)
                self._plotter.set_line(label=llab)
                self._plotter.plot(x,y,m)
                xlim=[min(x),max(x)]
                self._plotter.axes.set_xlim(xlim)

            self._plotter.set_axes('xlabel',xlab)
            self._plotter.set_axes('ylabel',ylab)
            self._plotter.set_axes('title',tlab)
            i += 1
        return
    
    def _plot_other(self,scans,colmode):
        if colmode == self._panels:
            return
        cdict = {'b':'scan.setbeam(j)',
                 'i':'scan.setif(j)',
                 'p':'scan.setpol(j)',
                 's':'scans[j]'}
        scan = scans[0]
        n = eval(self._cdict.get(self._panels))
        if self._stacking is not None:            
            ncol = eval(self._cdict.get(colmode))
        self._plotter.set_panels()
        if n > 1:
            self._plotter.set_panels(rows=n)
        for i in range(n):
            if n>1:
                self._plotter.subplot(i)
                self._plotter.palette(0)
            k=0
            j=i
            eval(cdict.get(self._panels))
            for j in range(ncol):
                if colmode == 's':
                    scan = eval(cdict.get(colmode))
                elif colmode == 't':
                    k = j
                else:
                    eval(cdict.get(colmode))
                x = None
                y = None
                m = None
                x,xlab = scan.get_abcissa(k)
                y = scan.getspectrum(k)
                ylab = 'Flux ('+scan.get_fluxunit()+')'
                m = scan.getmask(k)
                if colmode == 's' or colmode == 't':
                    if not self._title:
                        tlab = self._ldict.get(self._panels)+' '+str(i)
                    else:
                        if len(self.title) == n:
                            tlab = self._title[i]
                        else:
                            tlab = self._ldict.get(self._panels)+' '+str(i)
                    llab = scan._getsourcename(k)
                else:
                    if self._title and len(self._title) > 0:
                        tlab = self._title[k]
                    else:
                        tlab = scan._getsourcename(k)
                    if self._lmap and len(self._lmap) > 0:
                        llab = self._lmap[j]
                    else:
                        llab = self._ldict.get(colmode)+' '+str(j)
                self._plotter.set_line(label=llab)
                self._plotter.plot(x,y,m)
                xlim=[min(x),max(x)]
                self._plotter.axes.set_xlim(xlim)

            self._plotter.set_axes('xlabel',xlab)
            self._plotter.set_axes('ylabel',ylab)
            self._plotter.set_axes('title',tlab)
        
        return


    def set_mode(self, stacking=None, panelling=None):
        """
        Parameters:
            stacking:     tell the plotter which variable to plot
                          as line colour overlays (default 'pol')
            panelling:    tell the plotter which variable to plot
                          across multiple panels (default 'scan'
        Note:
            Valid modes are:
                 'beam' 'Beam' 'b':     Beams
                 'if' 'IF' 'i':         IFs
                 'pol' 'Pol' 'p':       Polarisations
                 'scan' 'Scan' 's':     Scans
                 'time' 'Time' 't':     Times
        """
        if not self.set_panels(panelling):
            print "Invalid mode"
            return
        if not self.set_stacking(stacking):
            print "Invalid mode"
            return
        if self._data: self.plot()
        return

    def set_panels(self, what=None):        
        if not what:
             what = rcParams['plotter.panelling']
        md = self._translate(what)
        if md:
            self._panels = md
            self._title = None
            return True
        return False

    def set_stacking(self, what=None):  
        if not what:
             what = rcParams['plotter.stacking']        
        md = self._translate(what)
        if md:
            self._stacking = md
            self._lmap = None
            return True
        return False

    def set_range(self,start=None,end=None):
        """
        Set the range of interest on the abcissa of the plot
        Parameters:
            start,end:    The start an end point of the 'zoom' window
        Note:
            These become non-sensical when the unit changes.
            use plotter.set_range() without parameters to reset

        """
        if start is None and end is None:
            self._minmax = None
            if self._data: self.plot()
        else:
            self._minmax = [start,end]
            if self._data: self.plot()
        return
    
    def set_legend_map(self, mp=[]):
        """
        Specify a mapping for the legend instead of using the default
        indices:
        Parameters:
             mp:    a list of 'strings'. This should have the same length
                    as the number of elements on the legend and then maps
                    to the indeces in order

        Example:
             If the data has to IFs/rest frequencies with index 0 and 1
             for CO and SiO:
             plotter.set_stacking('i')
             plotter.set_legend_map(['CO','SiO'])
             plotter.plot()
        """
        self._lmap = mp
        if self._data: self.plot()
        return

    def set_title(self, title=None):
        self._title = title
        if self._data: self.plot()
        return

if __name__ == '__main__':
    plotter = asapplotter()
