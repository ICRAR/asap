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

# Force use of the newfangled toolbar.
matplotlib.rcParams['toolbar'] = 'toolbar2'

# Colour dictionary.
colours = {}

class ASAPlot:
    """
    ASAP plotting class based on matplotlib.
    """

    def __init__(self, rowcol='11', title='', size=(8,6), buffering=False):
	"""
	Create a new instance of the ASAPlot plotting class.
	"""
	self.window = Tk.Tk()
	
	self.frame1 = Tk.Frame(self.window, relief=Tk.RIDGE, borderwidth=4)
	self.frame1.pack(fill=Tk.BOTH)

	self.figure = Figure(figsize=size, facecolor='#ddddee')
	self.canvas = FigureCanvasTkAgg(self.figure, master=self.frame1)
	self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

	# Simply instantiating this is enough to get a working toolbar.
	self.figmgr = FigureManagerTkAgg(self.canvas, 1, self.window)
	self.window.wm_title('ASAPlot graphics window')

	self.figure.text(0.5, 0.95, title, horizontalalignment='center')

	self.rows = int(rowcol[0])
	self.cols = int(rowcol[1])
	self.subplots = []
	for i in range(0,self.rows*self.cols):
	    self.subplots.append({})
	    self.subplots[i]['axes']  = self.figure.add_subplot(self.rows,
					    self.cols, i+1)
	    self.subplots[i]['lines'] = []

	self.figmgr.toolbar.set_active([0])

	self.axes  = self.subplots[0]['axes']
	self.lines = self.subplots[0]['lines']


	# Set matplotlib default colour sequence.
	self.colours = [1, 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
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


    def show(self):
	"""
	Show graphics dependent on the current buffering state.
	"""
	if not self.buffering:
	    if self.loc:
		lines  = []
		labels = []
		i = 0
		for line in self.lines:
		    i += 1
		    if line is not None:
			lines.append(line[0])
			lbl = line[0].get_label()
			if lbl == '':
			    lbl = str(i)
			labels.append(lbl)

		if len(lines):
		    self.axes.legend(tuple(lines), tuple(labels), self.loc)
		else:
		    self.axes.legend((' '))

	    self.window.wm_deiconify()
	    self.canvas.show()


    def subplot(self, col=0, row=0):
	"""
	Set the subplot; 0-relative column and row numbers.  Overrange column
	numbers map onto successive rows.
	"""
	i = (self.cols*row + col)%(self.rows*self.cols)
	self.axes  = self.subplots[i]['axes']
	self.lines = self.subplots[i]['lines']


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
