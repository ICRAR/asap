"""
ASAP plotting class based on matplotlib.
"""

import sys
from re import match
import Tkinter as Tk

import matplotlib
matplotlib.use("TkAgg")

from matplotlib.backends import new_figure_manager, show
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, \
	FigureManagerTkAgg
from matplotlib.figure import Figure, Text
from matplotlib.font_manager import FontProperties
from matplotlib.numerix import sqrt
from matplotlib import rc, rcParams

# Force use of the newfangled toolbar.
matplotlib.rcParams['toolbar'] = 'toolbar2'

class ASAPlot:
    """
    ASAP plotting class based on matplotlib.
    """

    def __init__(self, rows=1, cols=0, title='', size=(8,6), buffering=False):
	"""
	Create a new instance of the ASAPlot plotting class.

	If rows < 1 then a separate call to set_panels() is required to define
	the panel layout; refer to the doctext for set_panels().
	"""
	self.window = Tk.Tk()
	self.is_dead = False
	def dest_callback():
	    self.is_dead = True
	    self.window.destroy()

	self.window.protocol("WM_DELETE_WINDOW", dest_callback)

	self.figure = Figure(figsize=size, facecolor='#ddddee')
	self.canvas = FigureCanvasTkAgg(self.figure, master=self.window)
	self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

	# Simply instantiating this is enough to get a working toolbar.
	self.figmgr = FigureManagerTkAgg(self.canvas, 1, self.window)
	self.window.wm_title('ASAPlot graphics window')

	self.events = {'button_press':None,
		       'button_release':None,
		       'motion_notify':None}

	self.set_title(title)
	self.subplots = []
	if rows > 0:
	    self.set_panels(rows, cols)


	# Set matplotlib default colour sequence.
	self.colormap = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'purple', 'orange', 'pink']
        self.color = 0;
	self.attributes = {}
	self.loc = 0

	matplotlib.interactive = True
	self.buffering = buffering

	self.canvas.show()


    def clear(self):
	"""
	Delete all lines from the plot.  Line numbering will restart from 1.
	"""

	for i in range(len(self.lines)):
	   self.delete(i)
	self.axes.clear()
	self.color = 0
	self.lines = []


    def palette(self, color, colormap=None):
        if colormap:
            self.colormap = colormap
        if 0 <= color < len(self.colormap):
            self.color = color

    def delete(self, numbers=None):
	"""
	Delete the 0-relative line number, default is to delete the last.
	The remaining lines are NOT renumbered.
	"""

	if numbers is None: numbers = [len(self.lines)-1]

	if not hasattr(numbers, '__iter__'):
	    numbers = [numbers]

	for number in numbers:
	    if 0 <= number < len(self.lines):
		if self.lines[number] is not None:
		    for line in self.lines[number]:
			line.set_linestyle('None')
			self.lines[number] = None
        self.show()        

    def get_line(self):
	"""
	Get the current default line attributes.
	"""
	return self.attributes


    def hist(self, x=None, y=None, fmt=None):
	"""
	Plot a histogram.  N.B. the x values refer to the start of the
	histogram bin.

	fmt is the line style as in plot().
	"""

	if x is None:
	    if y is None: return
	    x = range(0,len(y))

	if len(x) != len(y):
	    return

	l2 = 2*len(x)
	x2 = range(0,l2)
	y2 = range(0,l2)

	for i in range(0,l2):
	    x2[i] = x[i/2]

	y2[0] = 0
	for i in range(1,l2):
	    y2[i] = y[(i-1)/2]

	self.plot(x2, y2, fmt)


    def hold(self, hold=True):
	"""
	Buffer graphics until subsequently released.
	"""
	self.buffering = hold


    def legend(self, loc=None):
	"""
	Add a legend to the plot.

	Any other value for loc else disables the legend:
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

	"""
        if isinstance(loc,int):
            if 0 > loc > 10: loc = 0
            self.loc = loc
	self.show()


    def map(self):
	"""
	Reveal the ASAPlot graphics window and bring it to the top of the
	window stack.
	"""
	self.window.wm_deiconify()
	self.window.lift()



    def plot(self, x=None, y=None, mask=None, fmt=None, add=None):
	"""
	Plot the next line in the current frame using the current line
	attributes.  The ASAPlot graphics window will be mapped and raised.

	The argument list works a bit like the matlab plot() function.
	"""

	if x is None:
	    if y is None: return
	    x = range(len(y))

	elif y is None:
	    y = x
	    x = range(len(y))

	if mask is None:
	    if fmt is None:
		line = self.axes.plot(x, y)
	    else:
		line = self.axes.plot(x, y, fmt)
	else:
	    segments = []

	    mask = list(mask)
	    i = 0
	    while mask[i:].count(1):
		i += mask[i:].index(1)
		if mask[i:].count(0):
		    j = i + mask[i:].index(0)
		else:
		    j = len(mask)

		segments.append(x[i:j])
		segments.append(y[i:j])

		i = j

	    line = self.axes.plot(*segments)

	# Add to an existing line?
	if add is None or len(self.lines) < add < 0:
	    # Don't add.
	    self.lines.append(line)
	    i = len(self.lines) - 1
	else:
	    if add == 0: add = len(self.lines)
	    i = add - 1
	    self.lines[i].extend(line)

	# Set/reset attributes for the line.
	gotcolour = False
	for k, v in self.attributes.iteritems():
	    if k == 'color': gotcolour = True
	    for segment in self.lines[i]:
		getattr(segment, "set_%s"%k)(v)

	if not gotcolour and len(self.colormap):
	    for segment in self.lines[i]:
		getattr(segment, "set_color")(self.colormap[self.color])

	    self.color += 1
	    if self.color >= len(self.colormap):
		self.color = 0

	self.show()


    def position(self):
	"""
	Use the mouse to get a position from a graph.
	"""

	def position_disable(event):
	    self.register('button_press', None)
	    print '%.4f, %.4f' % (event.xdata, event.ydata)

	print 'Press any mouse button...'
	self.register('button_press', position_disable)


    def quit(self):
	"""
	Destroy the ASAPlot graphics window.
	"""
	self.window.destroy()


    def region(self):
	"""
	Use the mouse to get a rectangular region from a plot.

	The return value is [x0, y0, x1, y1] in world coordinates.
	"""

	def region_start(event):
	    height = self.canvas.figure.bbox.height()
	    self.rect = {'fig': None, 'height': height,
			 'x': event.x, 'y': height - event.y,
			 'world': [event.xdata, event.ydata,
				   event.xdata, event.ydata]}
	    self.register('button_press', None)
	    self.register('motion_notify', region_draw)
	    self.register('button_release', region_disable)

	def region_draw(event):
	    self.canvas._tkcanvas.delete(self.rect['fig'])
	    self.rect['fig'] = self.canvas._tkcanvas.create_rectangle(
				self.rect['x'], self.rect['y'],
				event.x, self.rect['height'] - event.y)

	def region_disable(event):
	    self.register('motion_notify', None)
	    self.register('button_release', None)

	    self.canvas._tkcanvas.delete(self.rect['fig'])

	    self.rect['world'][2:4] = [event.xdata, event.ydata]
	    print '(%.2f, %.2f)  (%.2f, %.2f)' % (self.rect['world'][0],
		self.rect['world'][1], self.rect['world'][2],
		self.rect['world'][3])

	self.register('button_press', region_start)

	# This has to be modified to block and return the result (currently
	# printed by region_disable) when that becomes possible in matplotlib.

	return [0.0, 0.0, 0.0, 0.0]


    def register(self, type=None, func=None):
	"""
	Register, reregister, or deregister events of type 'button_press',
	'button_release', or 'motion_notify'.

	The specified callback function should have the following signature:

	    def func(event)

	where event is an MplEvent instance containing the following data:

	    name		# Event name.
	    canvas		# FigureCanvas instance generating the event.
	    x      = None	# x position - pixels from left of canvas.
	    y      = None	# y position - pixels from bottom of canvas.
	    button = None	# Button pressed: None, 1, 2, 3.
	    key    = None	# Key pressed: None, chr(range(255)), shift,
				  win, or control
	    inaxes = None	# Axes instance if cursor within axes.
	    xdata  = None	# x world coordinate.
	    ydata  = None	# y world coordinate.

	For example:

	    def mouse_move(event):
		print event.xdata, event.ydata

	    a = asaplot()
	    a.register('motion_notify', mouse_move)

	If func is None, the event is deregistered.

	Note that in TkAgg keyboard button presses don't generate an event.
	"""

	if not self.events.has_key(type): return

	if func is None:
	    if self.events[type] is not None:
		# It's not clear that this does anything.
		self.canvas.mpl_disconnect(self.events[type])
		self.events[type] = None

		# It seems to be necessary to return events to the toolbar.
		if type == 'motion_notify':
		    self.canvas.mpl_connect(type + '_event',
			self.figmgr.toolbar.mouse_move)
		elif type == 'button_press':
		    self.canvas.mpl_connect(type + '_event',
			self.figmgr.toolbar.press)
		elif type == 'button_release':
		    self.canvas.mpl_connect(type + '_event',
			self.figmgr.toolbar.release)

	else:
	    self.events[type] = self.canvas.mpl_connect(type + '_event', func)


    def release(self):
	"""
	Release buffered graphics.
	"""
	self.buffering = False
	self.show()


    def save(self, fname=None):
	"""
	Save the plot to a file.

	fname is the name of the output file.  The image format is determined
	from the file suffix; 'png', 'ps', and 'eps' are recognized.  If no
	file name is specified 'yyyymmdd_hhmmss.png' is created in the current
	directory.
	"""
	if fname is None:
	    from datetime import datetime
	    dstr = datetime.now().strftime('%Y%m%d_%H%M%S')
	    fname = 'asap'+dstr+'.png'

	d = ['png','.ps','eps']

	from os.path import expandvars
	fname = expandvars(fname)

	if fname[-3:].lower() in d:
	    try:
		self.canvas.print_figure(fname)
		print 'Written file %s' % (fname)
	    except IOError, msg:
		print 'Failed to save %s: Error msg was\n\n%s' % (fname, err)
		return
	else:
	    print "Invalid image type. Valid types are:"
	    print "ps, eps, png"


    def set_axes(self, what=None, *args, **kwargs):
	"""
	Set attributes for the axes by calling the relevant Axes.set_*()
	method.  Colour translation is done as described in the doctext
	for palette().
	"""

	if what is None: return
	if what[-6:] == 'colour': what = what[:-6] + 'color'

	newargs = {}
        
	for k, v in kwargs.iteritems():
	    k = k.lower()
	    if k == 'colour': k = 'color'
	    newargs[k] = v

	getattr(self.axes, "set_%s"%what)(*args, **newargs)
	self.show()


    def set_figure(self, what=None, *args, **kwargs):
	"""
	Set attributes for the figure by calling the relevant Figure.set_*()
	method.  Colour translation is done as described in the doctext
	for palette().
	"""

	if what is None: return
	if what[-6:] == 'colour': what = what[:-6] + 'color'
	#if what[-5:] == 'color' and len(args):
	#    args = (get_colour(args[0]),)

	newargs = {}
	for k, v in kwargs.iteritems():
	    k = k.lower()
	    if k == 'colour': k = 'color'
	    newargs[k] = v

	getattr(self.figure, "set_%s"%what)(*args, **newargs)
	self.show()


    def set_limits(self, xlim=None, ylim=None):
	"""
	Set x-, and y-limits for each subplot.

	xlim = [xmin, xmax] as in axes.set_xlim().
	ylim = [ymin, ymax] as in axes.set_ylim().
	"""
	for s in self.subplots:
	    self.axes  = s['axes']
	    self.lines = s['lines']
            oldxlim =  list(self.axes.get_xlim())
            oldylim =  list(self.axes.get_ylim())
            if xlim is not None:
                for i in range(len(xlim)):
                    if xlim[i] is not None:
                        oldxlim[i] = xlim[i]
            if ylim is not None:                        
                for i in range(len(ylim)):
                    if ylim[i] is not None:
                        oldylim[i] = ylim[i]
            self.axes.set_xlim(oldxlim)
            self.axes.set_ylim(oldylim)
        return


    def set_line(self, number=None, **kwargs):
	"""
	Set attributes for the specified line, or else the next line(s)
	to be plotted.

	number is the 0-relative number of a line that has already been
	plotted.  If no such line exists, attributes are recorded and used
	for the next line(s) to be plotted.

	Keyword arguments specify Line2D attributes, e.g. color='r'.  Do

	    import matplotlib
	    help(matplotlib.lines)

	The set_* methods of class Line2D define the attribute names and
	values.  For non-US usage, "colour" is recognized as synonymous with
	"color".

	Set the value to None to delete an attribute.

	Colour translation is done as described in the doctext for palette().
	"""

	redraw = False
	for k, v in kwargs.iteritems():
	    k = k.lower()
	    if k == 'colour': k = 'color'

	    if 0 <= number < len(self.lines):
		if self.lines[number] is not None:
		    for line in self.lines[number]:
			getattr(line, "set_%s"%k)(v)
		    redraw = True
	    else:
		if v is None:
		    del self.attributes[k]
		else:
		    self.attributes[k] = v

	if redraw: self.show()


    def set_panels(self, rows=1, cols=0, n=-1, nplots=-1, ganged=True):
	"""
	Set the panel layout.

	rows and cols, if cols != 0, specify the number of rows and columns in
	a regular layout.   (Indexing of these panels in matplotlib is row-
	major, i.e. column varies fastest.)

	cols == 0 is interpreted as a retangular layout that accomodates
	'rows' panels, e.g. rows == 6, cols == 0 is equivalent to
	rows == 2, cols == 3.

	0 <= n < rows*cols is interpreted as the 0-relative panel number in
	the configuration specified by rows and cols to be added to the
	current figure as its next 0-relative panel number (i).  This allows
	non-regular panel layouts to be constructed via multiple calls.  Any
	other value of n clears the plot and produces a rectangular array of
	empty panels.  The number of these may be limited by nplots.
	"""
	if n < 0 and len(self.subplots):
	    self.figure.clear()
	    self.set_title()

	if rows < 1: rows = 1

	if cols <= 0:
	    i = int(sqrt(rows))
	    if i*i < rows: i += 1
	    cols = i

	    if i*(i-1) >= rows: i -= 1
	    rows = i

	if 0 <= n < rows*cols:
	    i = len(self.subplots)
	    self.subplots.append({})

	    self.subplots[i]['axes']  = self.figure.add_subplot(rows,
					    cols, n+1)
	    self.subplots[i]['lines'] = []

	    if i == 0: self.subplot(0)

	    self.rows = 0
	    self.cols = 0

	else:
	    self.subplots = []

	    if nplots < 1 or rows*cols < nplots:
		nplots = rows*cols

	    for i in range(nplots):
		self.subplots.append({})

		self.subplots[i]['axes']  = self.figure.add_subplot(rows,
						cols, i+1)
		self.subplots[i]['lines'] = []
                xfsize = self.subplots[i]['axes'].xaxis.label.get_size()-cols/2
                yfsize = self.subplots[i]['axes'].yaxis.label.get_size()-rows/2
                self.subplots[i]['axes'].xaxis.label.set_size(xfsize)
                self.subplots[i]['axes'].yaxis.label.set_size(yfsize)
                
                if ganged:
                    if rows > 1 or cols > 1:
                        # Squeeze the plots together.
                        pos = self.subplots[i]['axes'].get_position()
                        if cols > 1: pos[2] *= 1.2
                        if rows > 1: pos[3] *= 1.2
                        self.subplots[i]['axes'].set_position(pos)

                    # Suppress tick labelling for interior subplots.
                    if i <= (rows-1)*cols - 1:
                        if i+cols < nplots:
                            # Suppress x-labels for frames width
                            # adjacent frames
                            for tick in \
                                    self.subplots[i]['axes'].xaxis.majorTicks:
                                tick.label1On = False
                                self.subplots[i]['axes'].xaxis.label.set_visible(False)
                    if i%cols:
                        # Suppress y-labels for frames not in the left column.
                        for tick in self.subplots[i]['axes'].yaxis.majorTicks:
                            tick.label1On = False
                        self.subplots[i]['axes'].yaxis.label.set_visible(False)
                        

		self.rows = rows
		self.cols = cols

	    self.subplot(0)

    def set_title(self, title=None):
	"""
	Set the title of the plot window.  Use the previous title if title is
	omitted.
	"""
	if title is not None:
	    self.title = title

	self.figure.text(0.5, 0.95, self.title, horizontalalignment='center')


    def show(self):
	"""
	Show graphics dependent on the current buffering state.
	"""
	if not self.buffering:
	    if self.loc is not None:
		for j in range(len(self.subplots)):
		    lines  = []
		    labels = []
		    i = 0
		    for line in self.subplots[j]['lines']:
			i += 1
			if line is not None:
			    lines.append(line[0])
			    lbl = line[0].get_label()
			    if lbl == '':
				lbl = str(i)
			    labels.append(lbl)

		    if len(lines):
			self.subplots[j]['axes'].legend(tuple(lines),
							tuple(labels),
							self.loc)
		    else:
			self.subplots[j]['axes'].legend((' '))

	    self.window.wm_deiconify()
	    self.canvas.show()

    def subplot(self, i=None, inc=None):
	"""
	Set the subplot to the 0-relative panel number as defined by one or
	more invokations of set_panels().
	"""
	l = len(self.subplots)
	if l:
	    if i is not None:
		self.i = i

	    if inc is not None:
		self.i += inc

	    self.i %= l
	    self.axes  = self.subplots[self.i]['axes']
	    self.lines = self.subplots[self.i]['lines']


    def terminate(self):
	"""
	Clear the figure.
	"""
	self.window.destroy()


    def text(self, *args, **kwargs):
	"""
	Add text to the figure.
	"""
	self.figure.text(*args, **kwargs)
	self.show()


    def unmap(self):
	"""
	Hide the ASAPlot graphics window.
	"""
	self.window.wm_withdraw()
