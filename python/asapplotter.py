from asap import rcParams, print_log, selector
from numarray import logical_and

class asapplotter:
    """
    The ASAP plotter.
    By default the plotter is set up to plot polarisations
    'colour stacked' and scantables across panels.
    Note:
        Currenly it only plots 'spectra' not Tsys or
        other variables.
    """
    def __init__(self, visible=None):
        self._visible = rcParams['plotter.gui']
        if visible is not None:
            self._visible = visible
        self._plotter = self._newplotter()


        self._panelling = None
        self._stacking = None
        self.set_panelling()
        self.set_stacking()
        self._rows = None
        self._cols = None
        self._autoplot = False
        self._minmaxx = None
        self._minmaxy = None
        self._datamask = None
        self._data = None
        self._lmap = None
        self._title = None
        self._ordinate = None
        self._abcissa = None
        self._abcunit = None
        self._usermask = []
        self._maskselection = None
        self._selection = selector()
        self._hist = rcParams['plotter.histogram']

    def _translate(self, instr):
        keys = "s b i p t".split()
        if isinstance(instr, str):
            for key in keys:
                if instr.lower().startswith(key):
                    return key
        return None

    def _newplotter(self):
        if self._visible:
            from asap.asaplotgui import asaplotgui as asaplot
        else:
            from asap.asaplot import asaplot
        return asaplot()


    def plot(self, scan=None):
        """
        Plot a scantable.
        Parameters:
            scan:   a scantable
        Note:
            If a scantable was specified in a previous call
            to plot, no argument has to be given to 'replot'
            NO checking is done that the abcissas of the scantable
            are consistent e.g. all 'channel' or all 'velocity' etc.
        """
        if self._plotter.is_dead:
            self._plotter = self._newplotter()
        self._plotter.hold()
        self._plotter.clear()
        from asap import scantable
        if not self._data and not scan:
            print "please provide a scantable to plot"
        if isinstance(scan, scantable):
            if self._data is not None:
                if scan != self._data:
                    self._data = scan
                    # reset
                    self._reset()
            else:
                self._data = scan
                self._reset()
        # ranges become invalid when unit changes
        if self._abcunit and self._abcunit != self._data.get_unit():
            self._minmaxx = None
            self._minmaxy = None
            self._abcunit = self._data.get_unit()
            self._datamask = None
        self._plot(self._data)
        if self._minmaxy is not None:
            self._plotter.set_limits(ylim=self._minmaxy)
        self._plotter.release()
        print_log()
        return

    def set_mode(self, stacking=None, panelling=None):
        """
        Set the plots look and feel, i.e. what you want to see on the plot.
        Parameters:
            stacking:     tell the plotter which variable to plot
                          as line color overlays (default 'pol')
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
        msg = "Invalid mode"
        if not self.set_panelling(panelling) or \
               not self.set_stacking(stacking):
            if rcParams['verbose']:
                print msg
                return
            else:
                raise TypeError(msg)
        if self._data: self.plot(self._data)
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
        if self._data: self.plot(self._data)
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
        if self._data: self.plot(self._data)
        return

    def set_legend(self, mp=None):
        """
        Specify a mapping for the legend instead of using the default
        indices:
        Parameters:
             mp:    a list of 'strings'. This should have the same length
                    as the number of elements on the legend and then maps
                    to the indeces in order. It is possible to uses latex
                    math expression. These have to be enclosed in r'', e.g. r'$x^{2}$'

        Example:
             If the data has two IFs/rest frequencies with index 0 and 1
             for CO and SiO:
             plotter.set_stacking('i')
             plotter.set_legend(['CO','SiO'])
             plotter.plot()
             plotter.set_legend([r'$^{12}CO$', r'SiO'])
        """
        self._lmap = mp
        if self._data: self.plot(self._data)
        return

    def set_title(self, title=None):
        """
        Set the title of the plot. If multiple panels are plotted,
        multiple titles have to be specified.
        Example:
             # two panels are visible on the plotter
             plotter.set_title(["First Panel","Second Panel"])
        """
        self._title = title
        if self._data: self.plot(self._data)
        return

    def set_ordinate(self, ordinate=None):
        """
        Set the y-axis label of the plot. If multiple panels are plotted,
        multiple labels have to be specified.
        Parameters:
            ordinate:    a list of ordinate labels. None (default) let
                         data determine the labels
        Example:
             # two panels are visible on the plotter
             plotter.set_ordinate(["First Y-Axis","Second Y-Axis"])
        """
        self._ordinate = ordinate
        if self._data: self.plot(self._data)
        return

    def set_abcissa(self, abcissa=None):
        """
        Set the x-axis label of the plot. If multiple panels are plotted,
        multiple labels have to be specified.
        Parameters:
            abcissa:     a list of abcissa labels. None (default) let
                         data determine the labels
        Example:
             # two panels are visible on the plotter
             plotter.set_ordinate(["First X-Axis","Second X-Axis"])
        """
        self._abcissa = abcissa
        if self._data: self.plot(self._data)
        return

    def set_colors(self, colormap):
        """
        Set the colors to be used. The plotter will cycle through
        these colors when lines are overlaid (stacking mode).
        Parameters:
            colormap:     a list of colour names
        Example:
             plotter.set_colors("red green blue")
             # If for example four lines are overlaid e.g I Q U V
             # 'I' will be 'red', 'Q' will be 'green', U will be 'blue'
             # and 'V' will be 'red' again.
        """
        if isinstance(colormap,str):
            colormap = colormap.split()
        self._plotter.palette(0,colormap=colormap)
        if self._data: self.plot(self._data)

    def set_histogram(self, hist=True):
        """
        Enable/Disable histogram-like plotting.
        Parameters:
            hist:        True (default) or False. The fisrt default
                         is taken from the .asaprc setting
                         plotter.histogram
        """
        self._hist = hist
        if self._data: self.plot(self._data)

    def set_linestyles(self, linestyles):
        """
        Set the linestyles to be used. The plotter will cycle through
        these linestyles when lines are overlaid (stacking mode) AND
        only one color has been set.
        Parameters:
             linestyles:     a list of linestyles to use.
                             'line', 'dashed', 'dotted', 'dashdot',
                             'dashdotdot' and 'dashdashdot' are
                             possible

        Example:
             plotter.set_colors("black")
             plotter.set_linestyles("line dashed dotted dashdot")
             # If for example four lines are overlaid e.g I Q U V
             # 'I' will be 'solid', 'Q' will be 'dashed',
             # U will be 'dotted' and 'V' will be 'dashdot'.
        """
        if isinstance(linestyles,str):
            linestyles = linestyles.split()
        self._plotter.palette(color=0,linestyle=0,linestyles=linestyles)
        if self._data: self.plot(self._data)

    def save(self, filename=None, orientation=None, dpi=None):
        """
        Save the plot to a file. The know formats are 'png', 'ps', 'eps'.
        Parameters:
             filename:    The name of the output file. This is optional
                          and autodetects the image format from the file
                          suffix. If non filename is specified a file
                          called 'yyyymmdd_hhmmss.png' is created in the
                          current directory.
             orientation: optional parameter for postscript only (not eps).
                          'landscape', 'portrait' or None (default) are valid.
                          If None is choosen for 'ps' output, the plot is
                          automatically oriented to fill the page.
             dpi:         The dpi of the output non-ps plot
        """
        self._plotter.save(filename,orientation,dpi)
        return


    def set_mask(self, mask=None, selection=None):
        """
        Set a plotting mask for a specific polarization.
        This is useful for masking out "noise" Pangle outside a source.
        Parameters:
             mask:           a mask from scantable.create_mask
             selection:      the spectra to apply the mask to.
        Example:
             select = selector()
             select.setpolstrings("Pangle")
             plotter.set_mask(mymask, select)
        """
        if not self._data:
            msg = "Can only set mask after a first call to plot()"
            if rcParams['verbose']:
                print msg
                return
            else:
                raise RuntimeError(msg)
        if len(mask):
            if isinstance(mask, list) or isinstance(mask, tuple):
                self._usermask = array(mask)
            else:
                self._usermask = mask
        if mask is None and selection is None:
            self._usermask = []
            self._maskselection = None
        if isinstance(selection, selector):
            self._maskselection = {'b': selection.get_beams(),
                                   's': selection.get_scans(),
                                   'i': selection.get_ifs(),
                                   'p': selection.get_pols(),
                                   't': [] }
        else:
            self._maskselection = None
        self.plot(self._data)

    def _slice_indeces(self, data):
        mn = self._minmaxx[0]
        mx = self._minmaxx[1]
        asc = data[0] < data[-1]
        start=0
        end = len(data)-1
        inc = 1
        if not asc:
            start = len(data)-1
            end = 0
            inc = -1
        # find min index
        while data[start] < mn:
            start+= inc
        # find max index
        while data[end] > mx:
            end-=inc
        end +=1
        if start > end:
            return end,start
        return start,end

    def _reset(self):
        self._usermask = []
        self._usermaskspectra = None
        self.set_selection(None, False)

    def _plot(self, scan):
        savesel = scan.get_selection()
        sel = savesel +  self._selection
        d0 = {'s': 'SCANNO', 'b': 'BEAMNO', 'i':'IFNO',
              'p': 'POLNO', 'c': 'CYCLENO', 't' : 'TIME' }
        order = [d0[self._panelling],d0[self._stacking]]
        sel.set_order(order)
        scan.set_selection(sel)
        d = {'b': scan.getbeam, 's': scan.getscan,
             'i': scan.getif, 'p': scan.getpol, 't': scan._gettime }

        polmodes = dict(zip(self._selection.get_pols(),self._selection.get_poltypes()))
        n,nstack = self._get_selected_n(scan)
        maxpanel, maxstack = 16,8
        if n > maxpanel or nstack > maxstack:
            from asap import asaplog
            msg ="Scan to be plotted contains more than %d selections.\n" \
                  "Selecting first %d selections..." % (maxpanel,maxpanel)
            asaplog.push(msg)
            print_log()
            n = min(n,maxpanel)
            nstack = min(nstack,maxstack)

        if n > 1:
            ganged = rcParams['plotter.ganged']
            if self._rows and self._cols:
                n = min(n,self._rows*self._cols)
                self._plotter.set_panels(rows=self._rows,cols=self._cols,
                                         nplots=n,ganged=ganged)
            else:
                self._plotter.set_panels(rows=n,cols=0,nplots=n,ganged=ganged)
        else:
            self._plotter.set_panels()
        r=0
        nr = scan.nrow()
        a0,b0 = -1,-1
        allxlim = []
        allylim = []
        newpanel=True
        panelcount,stackcount = 0,0
        while r < nr:
            a = d[self._panelling](r)
            b = d[self._stacking](r)
            if a > a0 and panelcount < n:
                if n > 1:
                    self._plotter.subplot(panelcount)
                self._plotter.palette(0)
                #title
                xlab = self._abcissa and self._abcissa[panelcount] \
                       or scan._getabcissalabel()
                ylab = self._ordinate and self._ordinate[panelcount] \
                       or scan._get_ordinate_label()
                self._plotter.set_axes('xlabel',xlab)
                self._plotter.set_axes('ylabel',ylab)
                lbl = self._get_label(scan, r, self._panelling, self._title)
                if isinstance(lbl, list) or isinstance(lbl, tuple):
                    if 0 <= panelcount < len(lbl):
                        lbl = lbl[panelcount]
                    else:
                        # get default label
                        lbl = self._get_label(scan, r, self._panelling, None)
                self._plotter.set_axes('title',lbl)
                newpanel = True
                stackcount =0
                panelcount += 1
            if (b > b0 or newpanel) and stackcount < nstack:
                y = []
                if len(polmodes):
                    y = scan._getspectrum(r, polmodes[scan.getpol(r)])
                else:
                    y = scan._getspectrum(r)
                m = scan._getmask(r)
                if self._maskselection and len(self._usermask) == len(m):
                    if d[self._stacking](r) in self._maskselection[self._stacking]:
                        m = logical_and(m, self._usermask)
                x = scan._getabcissa(r)
                if self._minmaxx is not None:
                    s,e = self._slice_indeces(x)
                    x = x[s:e]
                    y = y[s:e]
                    m = m[s:e]
                if len(x) > 2048 and rcParams['plotter.decimate']:
                    fac = len(x)/2048
                    x = x[::fac]
                    m = m[::fac]
                    y = y[::fac]
                llbl = self._get_label(scan, r, self._stacking, self._lmap)
                if isinstance(llbl, list) or isinstance(llbl, tuple):
                    if 0 <= stackcount < len(llbl):
                        # use user label
                        llbl = llbl[stackcount]
                    else:
                        # get default label
                        llbl = self._get_label(scan, r, self._stacking, None)
                self._plotter.set_line(label=llbl)
                plotit = self._plotter.plot
                if self._hist: plotit = self._plotter.hist
                plotit(x,y,m)
                xlim= self._minmaxx or [min(x),max(x)]
                allxlim += xlim
                ylim= self._minmaxy or [min(y),max(y)]
                allylim += ylim
                stackcount += 1
                # last in colour stack -> autoscale x
                if stackcount == nstack:
                    allxlim.sort()
                    self._plotter.axes.set_xlim([allxlim[0],allxlim[-1]])
                    # clear
                    allxlim =[]

            newpanel = False
            a0=a
            b0=b
            # ignore following rows
            if (panelcount == n) and (stackcount == nstack):
                # last panel -> autoscale y if ganged
                if rcParams['plotter.ganged']:
                    allylim.sort()
                    self._plotter.set_limits(ylim=[allylim[0],allylim[-1]])
                break
            r+=1 # next row
        #reset the selector to the scantable's original
        scan.set_selection(savesel)

    def set_selection(self, selection=None, refresh=True):
        self._selection = isinstance(selection,selector) and selection or selector()
        d0 = {'s': 'SCANNO', 'b': 'BEAMNO', 'i':'IFNO',
              'p': 'POLNO', 'c': 'CYCLENO', 't' : 'TIME' }
        order = [d0[self._panelling],d0[self._stacking]]
        self._selection.set_order(order)
        if self._data and refresh: self.plot(self._data)

    def _get_selected_n(self, scan):
        d1 = {'b': scan.nbeam, 's': scan.nscan,
             'i': scan.nif, 'p': scan.npol, 't': scan.ncycle }
        d2 = { 'b': len(self._selection.get_beams()),
               's': len(self._selection.get_scans()),
               'i': len(self._selection.get_ifs()),
               'p': len(self._selection.get_pols()),
               't': len(self._selection.get_cycles()) }
        n =  d2[self._panelling] or d1[self._panelling]()
        nstack = d2[self._stacking] or d1[self._stacking]()
        return n,nstack

    def _get_label(self, scan, row, mode, userlabel=None):
        pms = dict(zip(self._selection.get_pols(),self._selection.get_poltypes()))
        if len(pms):
            poleval = scan._getpollabel(scan.getpol(row),pms[scan.getpol(row)])
        else:
            poleval = scan._getpollabel(scan.getpol(row),scan.poltype())
        d = {'b': "Beam "+str(scan.getbeam(row)),
             's': scan._getsourcename(row),
             'i': "IF"+str(scan.getif(row)),
             'p': poleval,
             't': scan._gettime(row) }
        return userlabel or d[mode]
