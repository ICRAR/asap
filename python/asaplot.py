"""
ASAP plotting class based on matplotlib.
"""

import sys
from re import match
import Tkinter as Tk

print "Importing matplotlib with TkAgg backend."
import matplotlib
matplotlib.use("TkAgg")

from matplotlib.backends import new_figure_manager, show
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, \
	FigureManagerTkAgg
from matplotlib.figure import Figure, Text
from matplotlib.numerix import sqrt

# Force use of the newfangled toolbar.
matplotlib.rcParams['toolbar'] = 'toolbar2'

# Colour dictionary.
colours = {}

class ASAPlot:
    """
    ASAP plotting class based on matplotlib.
    """

    def __init__(self, rows=1, cols=0, title='', size=(7,5), buffering=False):
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
	#self.frame1 = Tk.Frame(self.window, relief=Tk.RIDGE, borderwidth=4)
	#self.frame1.pack(fill=Tk.BOTH)

	self.figure = Figure(figsize=size,facecolor='#ddddee')
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
	self.colours = [1, 'b', 'g', 'r', 'c', 'm', 'y', 'k']
	self.attributes = {}
	self.loc = 1

	matplotlib.interactive = True
	self.buffering = buffering

	self.canvas.show()


    def clear(self):
	"""
	Delete all lines from the plot.  Line numbering will restart from 1.
	"""

	for i in range(1,len(self.lines)+1):
	   self.delete(i)

	self.axes.clear()
	self.colours[0] = 1
	self.lines = []
	

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


    def legend(self, loc=1):
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
	if 1 > loc > 10: loc = 0
	self.loc = loc
	self.show()


    def map(self):
	"""
	Reveal the ASAPlot graphics window and bring it to the top of the
	window stack.
	"""
	self.window.wm_deiconify()
	self.window.lift()


    def palette(self, pen=None, colours=None):
	"""
	Redefine the colour sequence.

	pen is the pen number to use for the next plot; this will be auto-
	incremented.

	colours is the list of pen colours.  Colour may be specified via
	the single letter values understood by matplotlib:

	    b: blue
	    g: green
	    r: red
	    c: cyan
	    m: magenta
	    y: yellow
	    k: black
	    w: white

	or via the full name as listed in the colour dictionary which is
	loaded by default by load_colours() from rgb.txt and listed by
	list_colours().
	"""

	if pen is None and colours is None:
	    self.colours = []
	    return

	if pen is None:
	    if not len(self.colours):
		self.colours = [1]
	else:
	    self.colours[0] = pen

	if colours is None:
	    return

	cols = []
	for col in colours:
	    cols.append(get_colour(col))

	self.colours[1:] = cols

	if 0 > self.colours[0] > len(self.colours):
	    self.colours[0] = 1


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

	if not gotcolour and len(self.colours):
	    for segment in self.lines[i]:
		getattr(segment, "set_color")(self.colours[self.colours[0]])

	    self.colours[0] += 1
	    if self.colours[0] >= len(self.colours):
		self.colours[0] = 1

	self.show()


    def quit(self):
	"""
	Destroy the ASAPlot graphics window.
	"""
	self.window.destroy()


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

	    if k == 'color':
		v = get_colour(v)

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
	if what[-5:] == 'color' and len(args):
	    args = (get_colour(args[0]),)

	newargs = {}
	for k, v in kwargs.iteritems():
	    k = k.lower()
	    if k == 'colour': k = 'color'

	    if k == 'color':
		v = get_colour(v)

	    newargs[k] = v

	getattr(self.figure, "set_%s"%what)(*args, **newargs)
	self.show()


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

	    if k == 'color':
		v = get_colour(v)

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


    def set_panels(self, rows=1, cols=0, n=-1):
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
	empty panels.
	"""
	if n < 0 and len(self.subplots):
	    self.figure.clear()
	    self.set_title()

	if rows < 1:
	    rows = 1
        nel = 1
	if cols == 0:
            nel = rows
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

	else:
	    self.subplots = []
	    for i in range(0,nel):
		self.subplots.append({})
		self.subplots[i]['axes']  = self.figure.add_subplot(rows,
						cols, i+1)
		self.subplots[i]['lines'] = []

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
	    if self.loc:
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
                        self.subplots[j]['axes'].legend(tuple(lines), tuple(labels), self.loc)
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

    def set_limits(self,xlim=None,ylim=None):
        for s in self.subplots:
	    self.axes  = s['axes']
	    self.lines = s['lines']
            if xlim is not None:
                self.axes.set_xlim(xlim)
            if ylim is not None:
                self.axes.set_ylim(ylim)
        return

    def save(self, fname=None):
        if fname is None:
            from datetime import datetime
            dstr = datetime.now().strftime('%Y%m%d_%H%M%S')
            fname = 'asap'+dstr+'.png'
            
        d = ['png','.ps','eps']
        if fname[-3:].lower() in d:
            try:
                self.canvas.print_figure(fname)
            except IOError, msg:
                print 'Failed to save %s: Error msg was\n\n%s' % (fname, err)
                return
        else:
            print "Invalid image type. Valid types are:"
            print d


def get_colour(colour='black'):
    """
    Look up a colour by name in the colour dictionary.  Matches are
    case-insensitive, insensitive to blanks, and 'gray' matches 'grey'.
    """

    if colour is None: return None

    if match('[rgbcmykw]$', colour): return colour
    if match('#[\da-fA-F]{6}$', colour): return colour

    if len(colours) == 0: load_colours()

    # Try a quick match.
    if colours.has_key(colour): return colours[colour]

    colour = colour.replace(' ','').lower()
    colour = colour.replace('gray','grey')
    for name in colours.keys():
	if name.lower() == colour:
	    return colours[name]

    return '#000000'


def list_colours():
    """
    List the contents of the colour dictionary sorted by name.
    """

    if len(colours) == 0: load_colours()

    names = colours.keys()
    names.sort()
    for name in names:
	print colours[name], name


def load_colours(file='/usr/local/lib/rgb.txt'):
    """
    Load the colour dictionary from the specified file.
    """
    print 'Loading colour dictionary from', file
    rgb = open(file, 'r')

    while True:
	line = rgb.readline()
	if line == '': break
	tmp = line.split()

	if len(tmp) == 4:
	    if tmp[3][:4] == 'gray': continue
	    if tmp[3].lower().find('gray') != -1: continue

	    name = tmp[3][0].upper() + tmp[3][1:]
	    r, g, b = int(tmp[0]), int(tmp[1]), int(tmp[2])
	    colours[name] = '#%2.2x%2.2x%2.2x' % (r, g, b)
