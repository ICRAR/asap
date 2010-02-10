from asap import rcParams, print_log, print_log_dec
from asap import selector, scantable
import matplotlib.axes
from matplotlib.font_manager import FontProperties
from matplotlib.text import Text

import re

class asapplotter:
    """
    The ASAP plotter.
    By default the plotter is set up to plot polarisations
    'colour stacked' and scantables across panels.
    Note:
        Currenly it only plots 'spectra' not Tsys or
        other variables.
    """
    def __init__(self, visible=None , **kwargs):
        self._visible = rcParams['plotter.gui']
        if visible is not None:
            self._visible = visible
        self._plotter = self._newplotter(**kwargs)

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
        self._fp = FontProperties()

    def _translate(self, instr):
        keys = "s b i p t".split()
        if isinstance(instr, str):
            for key in keys:
                if instr.lower().startswith(key):
                    return key
        return None

    def _newplotter(self, **kwargs):
        if self._visible:
            from asap.asaplotgui import asaplotgui as asaplot
        else:
            from asap.asaplot import asaplot
        return asaplot(**kwargs)

    @print_log_dec
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
            msg = "Input is not a scantable"
            if rcParams['verbose']:
                print msg
                return
            raise TypeError(msg)
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
        self._plotter.tidy()
        self._plotter.show(hardrefresh=False)
        return

    def gca(self):
        return self._plotter.figure.gca()

    def refresh(self):
        """Do a soft refresh"""
        self._plotter.figure.show()

    def create_mask(self, nwin=1, panel=0, color=None):
        """
        Interactively define a mask.It retruns a mask that is equivalent to
        the one created manually with scantable.create_mask.
        Parameters:
            nwin:       The number of mask windows to create interactively
                        default is 1.
            panel:      Which panel to use for mask selection. This is useful
                        if different IFs are spread over panels (default 0)
        """
        if self._data is None:
            return []
        outmask = []
        self._plotter.subplot(panel)
        xmin, xmax = self._plotter.axes.get_xlim()
        marg = 0.05*(xmax-xmin)
        self._plotter.axes.set_xlim(xmin-marg, xmax+marg)
        self.refresh()

        def cleanup(lines=False, texts=False, refresh=False):
            if lines:
                del self._plotter.axes.lines[-1]
            if texts:
                del self._plotter.axes.texts[-1]
            if refresh:
                self.refresh()

        for w in xrange(nwin):
            wpos = []
            self.text(0.05,1.0, "Add start boundary",
                      coords="relative", fontsize=10)
            point = self._plotter.get_point()
            cleanup(texts=True)
            if point is None:
                continue
            wpos.append(point[0])
            self.axvline(wpos[0], color=color)
            self.text(0.05,1.0, "Add end boundary", coords="relative", fontsize=10)
            point = self._plotter.get_point()
            cleanup(texts=True, lines=True)
            if point is None:
                self.refresh()
                continue
            wpos.append(point[0])
            self.axvspan(wpos[0], wpos[1], alpha=0.1,
                         edgecolor=color, facecolor=color)
            ymin, ymax = self._plotter.axes.get_ylim()
            outmask.append(wpos)

        self._plotter.axes.set_xlim(xmin, xmax)
        self.refresh()
        if len(outmask) > 0:
            return self._data.create_mask(*outmask)
        return []

    # forwards to matplotlib axes
    def text(self, *args, **kwargs):
        if kwargs.has_key("interactive"):
            if kwargs.pop("interactive"):
                pos = self._plotter.get_point()
                args = tuple(pos)+args
        self._axes_callback("text", *args, **kwargs)

    text.__doc__ = matplotlib.axes.Axes.text.__doc__

    def arrow(self, *args, **kwargs):
        if kwargs.has_key("interactive"):
            if kwargs.pop("interactive"):
                pos = self._plotter.get_region()
                dpos = (pos[0][0], pos[0][1],
                        pos[1][0]-pos[0][0],
                        pos[1][1] - pos[0][1])
                args = dpos + args
        self._axes_callback("arrow", *args, **kwargs)

    arrow.__doc__ = matplotlib.axes.Axes.arrow.__doc__

    def annotate(self, text, xy=None, xytext=None, **kwargs):
        if kwargs.has_key("interactive"):
            if kwargs.pop("interactive"):
                xy = self._plotter.get_point()
                xytext = self._plotter.get_point()
        if not kwargs.has_key("arrowprops"):
            kwargs["arrowprops"] = dict(arrowstyle="->")
        self._axes_callback("annotate", text, xy, xytext, **kwargs)

    annotate.__doc__ = matplotlib.axes.Axes.annotate.__doc__

    def axvline(self, *args, **kwargs):
        if kwargs.has_key("interactive"):
            if kwargs.pop("interactive"):
                pos = self._plotter.get_point()
                args = (pos[0],)+args
        self._axes_callback("axvline", *args, **kwargs)

    axvline.__doc__ = matplotlib.axes.Axes.axvline.__doc__

    def axhline(self, *args, **kwargs):
        if kwargs.has_key("interactive"):
            if kwargs.pop("interactive"):
                pos = self._plotter.get_point()
                args = (pos[1],)+args
        self._axes_callback("axhline", *args, **kwargs)

    axhline.__doc__ = matplotlib.axes.Axes.axhline.__doc__

    def axvspan(self, *args, **kwargs):
        if kwargs.has_key("interactive"):
            if kwargs.pop("interactive"):
                pos = self._plotter.get_region()
                dpos = (pos[0][0], pos[1][0])
                args = dpos + args
        self._axes_callback("axvspan", *args, **kwargs)
        # hack to preventy mpl from redrawing the patch
        # it seem to convert the patch into lines on every draw.
        # This doesn't happen in a test script???
        #del self._plotter.axes.patches[-1]

    axvspan.__doc__ = matplotlib.axes.Axes.axvspan.__doc__

    def axhspan(self, *args, **kwargs):
        if kwargs.has_key("interactive"):
            if kwargs.pop("interactive"):
                pos = self._plotter.get_region()
                dpos = (pos[0][1], pos[1][1])
                args = dpos + args

        self._axes_callback("axhspan", *args, **kwargs)
        # hack to preventy mpl from redrawing the patch
        # it seem to convert the patch into lines on every draw.
        # This doesn't happen in a test script???
        #del self._plotter.axes.patches[-1]

    axhspan.__doc__ = matplotlib.axes.Axes.axhspan.__doc__

    def _axes_callback(self, axesfunc, *args, **kwargs):
        panel = 0
        if kwargs.has_key("panel"):
            panel = kwargs.pop("panel")
        coords = None
        if kwargs.has_key("coords"):
            coords = kwargs.pop("coords")
            if coords.lower() == 'world':
                kwargs["transform"] = self._plotter.axes.transData
            elif coords.lower() == 'relative':
                kwargs["transform"] = self._plotter.axes.transAxes
        self._plotter.subplot(panel)
        self._plotter.axes.set_autoscale_on(False)
        getattr(self._plotter.axes, axesfunc)(*args, **kwargs)
        self._plotter.show(False)
        self._plotter.axes.set_autoscale_on(True)
    # end matplotlib.axes fowarding functions


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

    def set_legend(self, mp=None, fontsize = None, mode = 0):
        """
        Specify a mapping for the legend instead of using the default
        indices:
        Parameters:
            mp:        a list of 'strings'. This should have the same length
                       as the number of elements on the legend and then maps
                       to the indeces in order. It is possible to uses latex
                       math expression. These have to be enclosed in r'',
                       e.g. r'$x^{2}$'
            fontsize:  The font size of the label (default None)
            mode:      where to display the legend
                       Any other value for loc else disables the legend:
                        0: auto
                        1: upper right
                        2: upper left
                        3: lower left
                        4: lower right
                        5: right
                        6: center left
                        7: center right
                        8: lower center
                        9: upper center
                        10: center

        Example:
             If the data has two IFs/rest frequencies with index 0 and 1
             for CO and SiO:
             plotter.set_stacking('i')
             plotter.set_legend(['CO','SiO'])
             plotter.plot()
             plotter.set_legend([r'$^{12}CO$', r'SiO'])
        """
        self._lmap = mp
        self._plotter.legend(mode)
        if isinstance(fontsize, int):
            from matplotlib import rc as rcp
            rcp('legend', fontsize=fontsize)
        if self._data:
            self.plot(self._data)
        return

    def set_title(self, title=None, fontsize=None):
        """
        Set the title of the plot. If multiple panels are plotted,
        multiple titles have to be specified.
        Example:
             # two panels are visible on the plotter
             plotter.set_title(["First Panel","Second Panel"])
        """
        self._title = title
        if isinstance(fontsize, int):
            from matplotlib import rc as rcp
            rcp('axes', titlesize=fontsize)
        if self._data: self.plot(self._data)
        return

    def set_ordinate(self, ordinate=None, fontsize=None):
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
        if isinstance(fontsize, int):
            from matplotlib import rc as rcp
            rcp('axes', labelsize=fontsize)
            rcp('ytick', labelsize=fontsize)
        if self._data: self.plot(self._data)
        return

    def set_abcissa(self, abcissa=None, fontsize=None):
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
        if isinstance(fontsize, int):
            from matplotlib import rc as rcp
            rcp('axes', labelsize=fontsize)
            rcp('xtick', labelsize=fontsize)
        if self._data: self.plot(self._data)
        return

    def set_colors(self, colmap):
        """
        Set the colours to be used. The plotter will cycle through
        these colours when lines are overlaid (stacking mode).
        Parameters:
            colmap:     a list of colour names
        Example:
             plotter.set_colors("red green blue")
             # If for example four lines are overlaid e.g I Q U V
             # 'I' will be 'red', 'Q' will be 'green', U will be 'blue'
             # and 'V' will be 'red' again.
        """
        if isinstance(colmap,str):
            colmap = colmap.split()
        self._plotter.palette(0, colormap=colmap)
        if self._data: self.plot(self._data)

    # alias for english speakers
    set_colours = set_colors

    def set_histogram(self, hist=True, linewidth=None):
        """
        Enable/Disable histogram-like plotting.
        Parameters:
            hist:        True (default) or False. The fisrt default
                         is taken from the .asaprc setting
                         plotter.histogram
        """
        self._hist = hist
        if isinstance(linewidth, float) or isinstance(linewidth, int):
            from matplotlib import rc as rcp
            rcp('lines', linewidth=linewidth)
        if self._data: self.plot(self._data)

    def set_linestyles(self, linestyles=None, linewidth=None):
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
        if isinstance(linewidth, float) or isinstance(linewidth, int):
            from matplotlib import rc as rcp
            rcp('lines', linewidth=linewidth)
        if self._data: self.plot(self._data)

    def set_font(self, **kwargs):
        """
        Set font properties.
        Parameters:
            family:    one of 'sans-serif', 'serif', 'cursive', 'fantasy', 'monospace'
            style:     one of 'normal' (or 'roman'), 'italic'  or 'oblique'
            weight:    one of 'normal or 'bold'
            size:      the 'general' font size, individual elements can be adjusted
                       seperately
        """
        from matplotlib import rc as rcp
        fdict = {}
        for k,v in kwargs.iteritems():
            if v:
                fdict[k] = v
        self._fp = FontProperties(**fdict)
        if self._data:
            self.plot()

    def plot_lines(self, linecat=None, doppler=0.0, deltachan=10, rotate=90.0,
                   location=None):
        """
        Plot a line catalog.
        Parameters:
            linecat:      the linecatalog to plot
            doppler:      the velocity shift to apply to the frequencies
            deltachan:    the number of channels to include each side of the
                          line to determine a local maximum/minimum
            rotate:       the rotation (in degrees) )for the text label (default 90.0)
            location:     the location of the line annotation from the 'top',
                          'bottom' or alternate (None - the default)
        Notes:
        If the spectrum is flagged no line will be drawn in that location.
        """
        if not self._data:
            raise RuntimeError("No scantable has been plotted yet.")
        from asap._asap import linecatalog
        if not isinstance(linecat, linecatalog):
            raise ValueError("'linecat' isn't of type linecatalog.")
        if not self._data.get_unit().endswith("Hz"):
            raise RuntimeError("Can only overlay linecatalogs when data is in frequency.")
        from matplotlib.numerix import ma
        for j in range(len(self._plotter.subplots)):
            self._plotter.subplot(j)
            lims = self._plotter.axes.get_xlim()
            for row in range(linecat.nrow()):
                # get_frequency returns MHz
                base = { "GHz": 1000.0, "MHz": 1.0, "Hz": 1.0e-6 }
                restf = linecat.get_frequency(row)/base[self._data.get_unit()]
                c = 299792.458
                freq = restf*(1.0-doppler/c)
                if lims[0] < freq < lims[1]:
                    if location is None:
                        loc = 'bottom'
                        if row%2: loc='top'
                    else: loc = location
                    maxys = []
                    for line in self._plotter.axes.lines:
                        v = line._x
                        asc = v[0] < v[-1]

                        idx = None
                        if not asc:
                            if v[len(v)-1] <= freq <= v[0]:
                                i = len(v)-1
                                while i>=0 and v[i] < freq:
                                    idx = i
                                    i-=1
                        else:
                           if v[0] <= freq <= v[len(v)-1]:
                                i = 0
                                while  i<len(v) and v[i] < freq:
                                    idx = i
                                    i+=1
                        if idx is not None:
                            lower = idx - deltachan
                            upper = idx + deltachan
                            if lower < 0: lower = 0
                            if upper > len(v): upper = len(v)
                            s = slice(lower, upper)
                            y = line._y[s]
                            maxy = ma.maximum(y)
                            if isinstance( maxy, float):
                                maxys.append(maxy)
                    if len(maxys):
                        peak = max(maxys)
                        if peak > self._plotter.axes.get_ylim()[1]:
                            loc = 'bottom'
                    else:
                        continue
                    self._plotter.vline_with_label(freq, peak,
                                                   linecat.get_name(row),
                                                   location=loc, rotate=rotate)
        self._plotter.show(hardrefresh=False)


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
        while start > 0 and data[start] < mn:
            start+= inc
        # find max index
        while end > 0 and data[end] > mx:
            end-=inc
        if end > 0: end +=1
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

        polmodes = dict(zip(self._selection.get_pols(),
                            self._selection.get_poltypes()))
        # this returns either a tuple of numbers or a length  (ncycles)
        # convert this into lengths
        n0,nstack0 = self._get_selected_n(scan)
        if isinstance(n0, int): n = n0
        else: n = len(n0)
        if isinstance(nstack0, int): nstack = nstack0
        else: nstack = len(nstack0)
        maxpanel, maxstack = 16,16
        if n > maxpanel or nstack > maxstack:
            from asap import asaplog
            maxn = 0
            if nstack > maxstack: maxn = maxstack
            if n > maxpanel: maxn = maxpanel
            msg ="Scan to be plotted contains more than %d selections.\n" \
                  "Selecting first %d selections..." % (maxn, maxn)
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
                self._plotter.set_axes('xlabel', xlab)
                self._plotter.set_axes('ylabel', ylab)
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
                from matplotlib.numerix import logical_not, logical_and
                if self._maskselection and len(self._usermask) == len(m):
                    if d[self._stacking](r) in self._maskselection[self._stacking]:
                        m = logical_and(m, self._usermask)
                x = scan._getabcissa(r)
                from matplotlib.numerix import ma, array
                y = ma.masked_array(y,mask=logical_not(array(m,copy=False)))
                if self._minmaxx is not None:
                    s,e = self._slice_indeces(x)
                    x = x[s:e]
                    y = y[s:e]
                if len(x) > 1024 and rcParams['plotter.decimate']:
                    fac = len(x)/1024
                    x = x[::fac]
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
                if len(x) > 0:
                    plotit(x,y)
                    xlim= self._minmaxx or [min(x),max(x)]
                    allxlim += xlim
                    ylim= self._minmaxy or [ma.minimum(y),ma.maximum(y)]
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
        if self._fp is not None:
            for o in self._plotter.figure.findobj(Text):
                o.set_fontproperties(self._fp)


    def set_selection(self, selection=None, refresh=True, **kw):
        if selection is None:
            # reset
            if len(kw) == 0:
                self._selection = selector()
            else:
                # try keywords
                for k in kw:
                    if k not in selector.fields:
                        raise KeyError("Invalid selection key '%s', valid keys are %s" % (k, selector.fields))
                self._selection = selector(**kw)
        elif isinstance(selection, selector):
            self._selection = selection
        else:
            raise TypeError("'selection' is not of type selector")

        d0 = {'s': 'SCANNO', 'b': 'BEAMNO', 'i':'IFNO',
              'p': 'POLNO', 'c': 'CYCLENO', 't' : 'TIME' }
        order = [d0[self._panelling],d0[self._stacking]]
        self._selection.set_order(order)
        if self._data and refresh: self.plot(self._data)

    def _get_selected_n(self, scan):
        d1 = {'b': scan.getbeamnos, 's': scan.getscannos,
             'i': scan.getifnos, 'p': scan.getpolnos, 't': scan.ncycle }
        d2 = { 'b': self._selection.get_beams(),
               's': self._selection.get_scans(),
               'i': self._selection.get_ifs(),
               'p': self._selection.get_pols(),
               't': self._selection.get_cycles() }
        n =  d2[self._panelling] or d1[self._panelling]()
        nstack = d2[self._stacking] or d1[self._stacking]()
        return n,nstack

    def _get_label(self, scan, row, mode, userlabel=None):
        if isinstance(userlabel, list) and len(userlabel) == 0:
            userlabel = " "
        pms = dict(zip(self._selection.get_pols(),self._selection.get_poltypes()))
        if len(pms):
            poleval = scan._getpollabel(scan.getpol(row),pms[scan.getpol(row)])
        else:
            poleval = scan._getpollabel(scan.getpol(row),scan.poltype())
        d = {'b': "Beam "+str(scan.getbeam(row)),
             's': scan._getsourcename(row),
             'i': "IF"+str(scan.getif(row)),
             'p': poleval,
             't': str(scan.get_time(row)) }
        return userlabel or d[mode]

    def plotazel(self):
        """
        plot azimuth and elevation  versus time of a scantable
        """
        import pylab as PL
        from matplotlib.dates import DateFormatter, timezone, HourLocator, MinuteLocator, DayLocator
        from matplotlib.ticker import MultipleLocator
        from matplotlib.numerix import array, pi
        dates = self._data.get_time(asdatetime=True)
        t = PL.date2num(dates)
        tz = timezone('UTC')
        PL.cla()
        PL.ioff()
        PL.clf()
        tdel = max(t) - min(t)
        ax = PL.subplot(2,1,1)
        el = array(self._data.get_elevation())*180./pi
        PL.ylabel('El [deg.]')
        dstr = dates[0].strftime('%Y/%m/%d')
        if tdel > 1.0:
            dstr2 = dates[len(dates)-1].strftime('%Y/%m/%d')
            dstr = dstr + " - " + dstr2
            majloc = DayLocator()
            minloc = HourLocator(range(0,23,12))
            timefmt = DateFormatter("%b%d")
        else:
            timefmt = DateFormatter('%H')
            majloc = HourLocator()
            minloc = MinuteLocator(20)
        PL.title(dstr)
        PL.plot_date(t,el,'b,', tz=tz)
        #ax.grid(True)
        ax.yaxis.grid(True)
        yloc = MultipleLocator(30)
        ax.set_ylim(0,90)
        ax.xaxis.set_major_formatter(timefmt)
        ax.xaxis.set_major_locator(majloc)
        ax.xaxis.set_minor_locator(minloc)
        ax.yaxis.set_major_locator(yloc)
        if tdel > 1.0:
            labels = ax.get_xticklabels()
        #    PL.setp(labels, fontsize=10, rotation=45)
            PL.setp(labels, fontsize=10)
        # Az plot
        az = array(self._data.get_azimuth())*180./pi
        if min(az) < 0:
            for irow in range(len(az)):
                if az[irow] < 0: az[irow] += 360.0

        ax = PL.subplot(2,1,2)
        PL.xlabel('Time (UT)')
        PL.ylabel('Az [deg.]')
        PL.plot_date(t,az,'b,', tz=tz)
        ax.set_ylim(0,360)
        #ax.grid(True)
        ax.yaxis.grid(True)
        #hfmt = DateFormatter('%H')
        #hloc = HourLocator()
        yloc = MultipleLocator(60)
        ax.xaxis.set_major_formatter(timefmt)
        ax.xaxis.set_major_locator(majloc)
        ax.xaxis.set_minor_locator(minloc)
        ax.yaxis.set_major_locator(yloc)
        if tdel > 1.0:
            labels = ax.get_xticklabels()
            PL.setp(labels, fontsize=10)
        PL.ion()
        PL.draw()

    def plotpointing(self):
        """
        plot telescope pointings
        """
        import pylab as PL
        from matplotlib.dates import DateFormatter, timezone
        from matplotlib.ticker import MultipleLocator
        from matplotlib.numerix import array, pi, zeros
        dir = array(self._data.get_directionval()).transpose()
        ra = dir[0]*180./pi
        dec = dir[1]*180./pi
        PL.cla()
        PL.ioff()
        PL.clf()
        ax = PL.axes([0.1,0.1,0.8,0.8])
        ax = PL.axes([0.1,0.1,0.8,0.8])
        ax.set_aspect('equal')
        PL.plot(ra,dec, 'b,')
        PL.xlabel('RA [deg.]')
        PL.ylabel('Declination [deg.]')
        PL.title('Telescope pointings')
        [xmin,xmax,ymin,ymax] = PL.axis()
        PL.axis([xmax,xmin,ymin,ymax])
        PL.ion()
        PL.draw()
