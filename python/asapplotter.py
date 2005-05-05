from asap.asaplot import ASAPlot
from asap import rcParams

class asapplotter:
    """
    The ASAP plotter.
    By default the plotter is set up to plot polarisations
    'colour stacked' and scantables across panels.
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
        self._cdict = {'t':'len(self._cursor["t"])',
                       'b':'len(self._cursor["b"])',
                       'i':'len(self._cursor["i"])',
                       'p':'len(self._cursor["p"])',
                       's':'len(scans)'}
        self._ldict = {'b':'Beam',
                       'i':'IF',
                       'p':'Pol',
                       's':'Scan'}
        self._dicts = [self._tdict,self._bdict,
                       self._idict,self._pdict,
                       self._sdict]
        self._panelling = None
        self._stacking = None
        self.set_panelling()
        self.set_stacking()
        self._rows = None
        self._cols = None
        self._autoplot = False
        self._minmaxx = None
        self._minmaxy = None
        self._data = None
        self._lmap = None
        self._title = None
        self._ordinate = None
        self._abcissa = None
        self._cursor = {'t':None, 'b':None,
                        'i':None, 'p':None
                        }

    def _translate(self, name):
        for d in self._dicts:
            if d.has_key(name):
                return d[name]
        return None
        
    def plot(self, *args):
        """
        Plot a (list of) scantables.
        Parameters:
            one or more comma separated scantables 
        Note:
            If a (list) of scantables was specified in a previous call
            to plot, no argument has to be given to 'replot'
            NO checking is done that the abcissas of the scantables
            are consistent e.g. all 'channel' or all 'velocity' etc.
        """
        if self._plotter.is_dead:
            self._plotter = ASAPlot()
        self._plotter.hold()
        self._plotter.clear()
        if len(args) > 0:
            if self._data is not None:                
                if list(args) != self._data:
                    self._data = list(args)
                    # reset cursor
                    self.set_cursor(refresh=False)
            else:
                self._data = list(args)
                self.set_cursor(refresh=False)
        if self._panelling == 't':
            maxrows = 25
            if self._data[0].nrow() > maxrows:
                if self._cursor["t"] is None or \
                       (isinstance(self._cursor["t"],list) and \
                        len(self._cursor["t"]) > maxrows ):
                    print "Scan to be plotted contains more than %d rows.\n" \
                          "Selecting first %d rows..." % (maxrows,maxrows)
                    self._cursor["t"] = range(maxrows)
            self._plot_time(self._data[0], self._stacking)
        elif self._panelling == 's':
            self._plot_scans(self._data, self._stacking)
        else:
            self._plot_other(self._data, self._stacking)
        if self._minmaxx is not None or self._minmaxy is not None:
            self._plotter.set_limits(xlim=self._minmaxx,ylim=self._minmaxy)
        self._plotter.release()
        return

    def _plot_time(self, scan, colmode):
        if colmode == 't':
            return
        n = len(self._cursor["t"])
        cdict = {'b':'scan.setbeam(j)',
                 'i':'scan.setif(j)',
                 'p':'scan.setpol(j)'}
        cdict2 = {'b':'self._cursor["b"]',
                  'i':'self._cursor["i"]',
                  'p':'self._cursor["p"]'}
        ncol = 1
        if self._stacking is not None:
            ncol = eval(self._cdict.get(colmode))
        if n > 1:
            if self._rows and self._cols:
                n = min(n,self._rows*self._cols)
                self._plotter.set_panels(rows=self._rows,cols=self._cols,
                                         nplots=n)
            else:
                self._plotter.set_panels(rows=n,cols=0,nplots=n)
        else:
            self._plotter.set_panels()
        rows = self._cursor["t"]
        self._plotter.palette(1)
        for rowsel in rows:
            i = self._cursor["t"].index(rowsel)
            if n > 1:
                self._plotter.palette(1)
                self._plotter.subplot(i)
            colvals = eval(cdict2.get(colmode))
            for j in colvals:
                polmode = "raw"
                jj = colvals.index(j)
                savej = j
                for k in cdict.keys():
                    sel = eval(cdict2.get(k))                    
                    j = sel[0]
                    if k == "p":
                        which = self._cursor["p"].index(j)
                        polmode = self._polmode[which]
                        j = which
                    eval(cdict.get(k))
                j = savej
                if colmode == "p":
                    polmode = self._polmode[self._cursor["p"].index(j)]
                    j = jj
                eval(cdict.get(colmode))
                x = None
                y = None
                m = None
                if self._title is None:
                    tlab = scan._getsourcename(rowsel)                    
                else:
                    if len(self._title) >= n:
                        tlab = self._title[rowsel]
                    else:
                        tlab = scan._getsourcename(rowsel)
                x,xlab = scan.get_abcissa(rowsel)
                if self._abcissa: xlab = self._abcissa
                y = None
                if polmode == "stokes":
                    y = scan._getstokesspectrum(rowsel)
                elif polmode == "stokes2":
                    y = scan._getstokesspectrum(rowsel,True)
                elif polmode == "circular":
                    y = scan._stokestopolspectrum(rowsel,False,-1)
                else:
                    y = scan._getspectrum(rowsel)
                if self._ordinate:
                    ylab = self._ordinate
                else:
                    ylab = scan._get_ordinate_label()
                m = scan._getmask(rowsel)
                if self._lmap and len(self._lmap) > 0:
                    llab = self._lmap[jj]
                else:
                    if colmode == 'p':
                        llab = self._get_pollabel(scan, polmode)
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
        print "Can only plot one row per scan."
        if colmode == 's':
            return
        cdict = {'b':'scan.setbeam(j)',
                 'i':'scan.setif(j)',
                 'p':'scan.setpol(j)'}
        cdict2 = {'b':'self._cursor["b"]',
                  'i':'self._cursor["i"]',
                  'p':'self._cursor["p"]'}
        
        n = len(scans)
        ncol = 1
        if self._stacking is not None:
            scan = scans[0]
            ncol = eval(self._cdict.get(colmode))
        if n > 1:
            if self._rows and self._cols:
                n = min(n,self._rows*self._cols)
                self._plotter.set_panels(rows=self._rows,cols=self._cols,
                                         nplots=n)
            else:
                self._plotter.set_panels(rows=n,cols=0,nplots=n)
        else:
            self._plotter.set_panels()

        for scan in scans:
            self._plotter.palette(1)
            if n > 1:
                self._plotter.subplot(scans.index(scan))
                self._plotter.palette(1)
            colvals = eval(cdict2.get(colmode))
            rowsel = self._cursor["t"][0]
            for j in colvals:
                polmode = "raw"
                jj = colvals.index(j)
                savej = j
                for k in cdict.keys():
                    sel = eval(cdict2.get(k))                    
                    j = sel[0]
                    eval(cdict.get(k))
                    if k == "p":
                        which = self._cursor["p"].index(j)
                        polmode = self._polmode[which]
                        j = which
                j = savej
                if colmode == "p":
                    polmode = self._polmode[self._cursor["p"].index(j)]
                    j = jj
                eval(cdict.get(colmode))
                x = None
                y = None
                m = None
                tlab = self._title
                if not self._title:
                    tlab = scan._getsourcename(rowsel)
                x,xlab = scan.get_abcissa(rowsel)
                if self._abcissa: xlab = self._abcissa
                if polmode == "stokes":
                    y = scan._getstokesspectrum(rowsel)
                elif polmode == "stokes2":
                    y = scan._getstokesspectrum(rowsel,True)
                elif polmode == "circular":
                    y = scan._stokestopolspectrum(rowsel,False,-1)
                else:
                    y = scan._getspectrum(rowsel)
                if self._ordinate:
                    ylab = self._ordinate
                else:
                    ylab = scan._get_ordinate_label()
                m = scan._getmask(rowsel)
                if self._lmap and len(self._lmap) > 0:
                    llab = self._lmap[jj]
                else:
                    if colmode == 'p':
                        llab = self._get_pollabel(scan, polmode)
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
    
    def _plot_other(self,scans,colmode):
        if colmode == self._panelling:
            return
        cdict = {'b':'scan.setbeam(i)',
                 'i':'scan.setif(i)',
                 'p':'scan.setpol(i)'}
        cdict2 = {'b':'self._cursor["b"]',
                  'i':'self._cursor["i"]',
                  'p':'self._cursor["p"]',
                  's': 'scans',
                  't': 'self._cursor["t"]'}
        scan = scans[0]
        n = eval(self._cdict.get(self._panelling))
        ncol=1
        if self._stacking is not None:            
            ncol = eval(self._cdict.get(colmode))
        if n > 1:
            if self._rows and self._cols:
                n = min(n,self._rows*self._cols)
                self._plotter.set_panels(rows=self._rows,cols=self._cols,
                                         nplots=n)
            else:
                self._plotter.set_panels(rows=n,cols=0,nplots=n)
        else:
            self._plotter.set_panels()            
        panels = self._cursor[self._panelling]        
        for i in panels:
            self._plotter.palette(1)
            polmode = "raw"
            ii = self._cursor[self._panelling].index(i)
            if n>1:
                self._plotter.subplot(ii)
            if self._panelling == "p":
                polmode = self._polmode[ii]
                eval(cdict.get(self._panelling))
            else:
                eval(cdict.get(self._panelling))
            colvals = eval(cdict2.get(colmode))
            for j in colvals:
                rowsel = self._cursor["t"][0]
                jj = colvals.index(j)
                savei = i
                for k in cdict.keys():
                    if k != self._panelling:
                        sel = eval(cdict2.get(k))
                        i = sel[0]
                        if k == "p":
                            which = self._cursor["p"].index(i)
                            polmode = self._polmode[which]
                            i = which                        
                        eval(cdict.get(k))
                i = savei
                if colmode == 's':
                    scan = j
                elif colmode == 't':
                    rowsel = j                    
                else:
                    savei = i
                    if colmode == 'p':
                        polmode = self._polmode[self._cursor["p"].index(j)]
                    i = j
                    eval(cdict.get(colmode))
                    i = savei
                x = None
                y = None
                m = None
                x,xlab = scan.get_abcissa(rowsel)
                if self._abcissa: xlab = self._abcissa
                if polmode == "stokes":
                    y = scan._getstokesspectrum(rowsel)
                elif polmode == "stokes2":
                    y = scan._getstokesspectrum(rowsel,True)
                elif polmode == "circular":
                    y = scan._stokestopolspectrum(rowsel,False,-1)
                else:
                    y = scan._getspectrum(rowsel)

                if self._ordinate:
                    ylab = self._ordinate
                else:
                    ylab = scan._get_ordinate_label()
                m = scan._getmask(rowsel)
                if colmode == 's' or colmode == 't':
                    if self._title and len(self._title) > 0:
                        tlab = self._title[ii]
                    else:                        
                        tlab = self._ldict.get(self._panelling)+' '+str(i)
                    if self._lmap and len(self._lmap) > 0:
                        llab = self._lmap[jj]
                    else:
                        llab = scan._getsourcename(rowsel)
                else:
                    if self._title and len(self._title) > 0:
                        tlab = self._title[ii]
                    else:
                        if self._panelling == 'p':
                            tlab = self._get_pollabel(scan, polmode)
                        else:
                            tlab = self._ldict.get(self._panelling)+' '+str(i)
                    if self._lmap and len(self._lmap) > 0:
                        llab = self._lmap[jj]
                    else:
                        if colmode == 'p':
                            llab = self._get_pollabel(scan, polmode)
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
        Set the plots look and feel, i.e. what you want to see on the plot.
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
        if not self.set_panelling(panelling):
            print "Invalid mode"
            return
        if not self.set_stacking(stacking):
            print "Invalid mode"
            return
        if self._data: self.plot()
        return

    def set_panelling(self, what=None):
        mode = what
        if mode is None:
             mode = rcParams['plotter.panelling']
        md = self._translate(mode)
        if md:
            self._panelling = md
            self._title = None
            return True
        return False

    def set_layout(self,rows=None,cols=None):
        """
        Set the multi-panel layout, i.e. how many rows and columns plots
        are visible.
        Parameters:
             rows:   The number of rows of plots
             cols:   The number of columns of plots
        Note:
             If no argument is given, the potter reverts to its auto-plot
             behaviour.
        """
        self._rows = rows
        self._cols = cols
        if self._data: self.plot()
        return

    def set_stacking(self, what=None):  
        mode = what
        if mode is None:            
             mode = rcParams['plotter.stacking']        
        md = self._translate(mode)
        if md:
            self._stacking = md
            self._lmap = None
            return True
        return False

    def set_range(self,xstart=None,xend=None,ystart=None,yend=None):
        """
        Set the range of interest on the abcissa of the plot
        Parameters:
            [x,y]start,[x,y]end:  The start and end points of the 'zoom' window
        Note:
            These become non-sensical when the unit changes.
            use plotter.set_range() without parameters to reset

        """
        if xstart is None and xend is None:
            self._minmaxx = None
        else:
            self._minmaxx = [xstart,xend]
        if ystart is None and yend is None:
            self._minmaxy = None
        else:
            self._minmaxy = [ystart,yend]            
        if self._data: self.plot()
        return
    
    def set_legend(self, mp=None):
        """
        Specify a mapping for the legend instead of using the default
        indices:
        Parameters:
             mp:    a list of 'strings'. This should have the same length
                    as the number of elements on the legend and then maps
                    to the indeces in order

        Example:
             If the data has two IFs/rest frequencies with index 0 and 1
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

    def set_ordinate(self, ordinate=None):
        self._ordinate = ordinate
        if self._data: self.plot()
        return

    def set_abcissa(self, abcissa=None):
        self._abcissa = abcissa
        if self._data: self.plot()
        return

    def save(self, filename=None):
        """
        Save the plot to a file. The know formats are 'png', 'ps', 'eps'.
        Parameters:
             filename:    The name of the output file. This is optional
                          and autodetects the image format from the file
                          suffix. If non filename is specified a file
                          called 'yyyymmdd_hhmmss.png' is created in the
                          current directory.
        """
        self._plotter.save(filename)
        return
    
    def set_cursor(self, row=None,beam=None,IF=None,pol=None, refresh=True):
        """
        Specify a 'cursor' for plotting selected spectra. Time (rows),
        Beam, IF, Polarisation ranges can be specified.
        Parameters:
            Default for all paramaters is to select all available
            row:    selects the rows (time stamps) to be plotted, this has
                    to be a vector of row indices, e.g. row=[0,2,5] or row=[2]
            beam:   select a range of beams
            IF:     select a range of IFs
            pol:    select Polarisations for plotting these can be by index
                    (raw polarisations (default)) or by names any of:
                    ["I", "Q", "U", "V"] or
                    ["I", "Plinear", "Pangle", "V"] or
                    ["XX", "YY", "Real(XY)", "Imag(XY)"] or
                    ["RR", "LL"]
        Example:
            plotter.set_mode('pol','time')
            plotter.plot(myscan) # plots all raw polarisations colour stacked
            plotter.set_cursor(pol=["I"]) # plot "I" only for all rows
            # plot "I" only for two time stamps row=0 and row=2
            plotter.set_cursor(row=[0,2],pol=["I"])

        Note:
            Be careful to select only exisiting polarisations.            
        """
        if not self._data:
            print "Can only set cursor after a first call to plot()"
            return
        
        n = self._data[0].nrow()
        if row is None:
            self._cursor["t"] = range(n)
        else:
            for i in row:
                if i < 0 or i >= n:
                    print "Row index '%d' out of range" % i
                    return
            self._cursor["t"] = row

        n = self._data[0].nbeam()
        if beam is None:
            self._cursor["b"] = range(n)
        else:
            for i in beam:
                if i < 0 or  i >= n:
                    print "Beam index '%d' out of range" % i
                    return            
            self._cursor["b"] = beam

        n = self._data[0].nif()
        if IF is None:
            self._cursor["i"] = range(n)
        else:
            for i in IF:
                if i < 0 or i >= n:
                    print "IF index '%d' out of range" %i
                    return            
            self._cursor["i"] = IF            

        n = self._data[0].npol()
        dstokes = {"I":0,"Q":1,"U":2,"V":3}
        dstokes2 = {"I":0,"Plinear":1,"Pangle":2,"V":3}
        draw = {"XX":0, "YY":1,"Real(XY)":2, "Imag(XY)":3}
        dcirc = { "RR":0,"LL":1}#,"Real(RL)":2,"Image(RL)":3}
        
        if pol is None:
            self._cursor["p"] = range(n)
            self._polmode = ["raw" for i in range(n)]
        else:
            if isinstance(pol,str):
                pol = pol.split()
            polmode = []
            pols = []
            for i in pol:
                if isinstance(i,str):
                    if draw.has_key(i):
                        pols.append(draw.get(i))
                        polmode.append("raw")
                    elif dstokes.has_key(i):
                        pols.append(dstokes.get(i))
                        polmode.append("stokes")
                    elif dstokes2.has_key(i):
                        pols.append(dstokes2.get(i))
                        polmode.append("stokes2")
                    elif dcirc.has_key(i):
                        pols.append(dcirc.get(i))
                        polmode.append("circular")
                    else:
                        "Pol type '%s' not valid" %i
                        return
                elif 0 > i >= n:
                    print "Pol index '%d' out of range" %i
                    return
                else:
                    pols.append(i)
                    polmode.append("raw")
            self._cursor["p"] = pols
            self._polmode = polmode
        if self._data and refresh: self.plot()

    def _get_pollabel(self, scan, polmode):
        tlab = ""
        if polmode == "stokes":
            tlab = scan._getpolarizationlabel(0,1,0)
        elif polmode == "stokes2":
            tlab = scan._getpolarizationlabel(0,1,1)
        elif polmode == "circular":
            tlab = scan._getpolarizationlabel(0,0,0)
        else:
            tlab = scan._getpolarizationlabel(1,0,0)
        return tlab
            
if __name__ == '__main__':
    plotter = asapplotter()
