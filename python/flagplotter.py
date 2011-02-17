from asap.asapplotter import asapplotter
from asap.logging import asaplog, asaplog_post_dec

from asap.parameters import rcParams
from asap.selector import selector
from asap.scantable import scantable
import matplotlib.axes
from matplotlib.font_manager import FontProperties
from matplotlib.text import Text

class flagplotter(asapplotter):
    """
    The flag plotter
    Only row based panneling is allowed.

    Example:
       scan = asap.scantable(filename='your_filename',average=False)
       guiflagger = asap.flagplotter(visible=True)
       guiflagger.plot(scan)
       ### flag/Unflag data graphically.
       scan.save(name='flagged_file.asap',format='ASAP')
    
    NOTICE: 
       The flagged data is not saved until you explicitly run scantable.save
    """
    def __init__(self, visible=None, **kwargs):
        self._scan=None
        asapplotter.__init__(self,visible=visible, **kwargs)
        self._plotter.window.title('Flag Plotter')
        self._panelling = 'r'
        self._stacking = 's'

    def _newcasabar(self):
        backend=matplotlib.get_backend()
        if self._visible and backend == "TkAgg":
            #from asap.casatoolbar import CustomToolbarTkAgg
            #return CustomToolbarTkAgg(self)
            from asap.flagtoolbar import CustomFlagToolbarTkAgg
            return CustomFlagToolbarTkAgg(self)
        return None

    def set_mode(self, stacking=None, panelling=None, refresh=True):
        """ This function is not available for the class flagplotter """
        self._invalid_func(name='set_mode')

    def set_panelling(self, what=None):
        """ This function is not available for the class flagplotter """
        self._invalid_func(name='set_panelling')

    def set_stacking(self, what=None):
        """ This function is not available for the class flagplotter """
        self._invalid_func(name='set_stacking')

    @asaplog_post_dec
    def _invalid_func(self, name):
        msg = "Invalid function 'flagplotter."+name+"'"
        asaplog.push(msg)
        asaplog.post('ERROR')

    def save_data(self, name=None, format=None, overwrite=False):
        # simply calls scantable.save
        self._data.save(name,format,overwrite)
