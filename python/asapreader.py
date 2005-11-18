from asap._asap import sdreader
from asap import print_log

class reader(sdreader):
    """
    This class allows the user to import single dish files
    (rpfits,sdfits,ms).
    The reader reads in integrations from the file and reamins at
    the fileposition afterwards.
    Available functions are:

    read(integrations)
    summary    CURRENTLY DISABLED

    Example:
        r = reader('/tmp/P389.rpf')
        scans r.read() # reads in the complete file into 'scans'
        print scans    # summarises the contents
        del r          # destroys the reader

    IMPORTANT: Due to limitations in the rpfits library, only one reader
               can be created at a time.
               r = reader('XYZ.rpf')
               r2 = reader('ABC.rpf')
               is NOT possible. This is a limitation affecting
               rpfits ONLY.
    """

    def __init__(self, filename, unit=None, theif=None, thebeam=None):
        self.unit = unit
        """
        Parameters:
            filename:    the name of an rpfits/sdfits/ms file on disk
            unit:        brightness unit; must be consistent with K or Jy.
                         The default is that a unit is set depending on
                         the telescope.  Setting this over-rides that choice.
            theif:       select a specific IF (default is all)
            thebeam:     select a specific beam (default is all)
        Example:
            r = reader('/tmp/2001-09-01_0332_P363.rpf', theif=2)
        """
        if theif is None:
            theif = -1
        if thebeam is None:
            thebeam = -1
        from os.path import expandvars
        filename = expandvars(filename)
        sdreader.__init__(self, filename, theif, thebeam)
        print_log()

    def read(self,integrations=None):
        """
        Reads in an returns a specified sequence of integrations.
        If no list is given all integrations a read in.
        Parameters:
            integrations:    a 'range' of integration numbers, e.g.
                             range(100) or [0,1,2,3,4,10,11,100]
                             If not given (default) all integrations
                             are read in
        Example:
            r.read([0,1,2,3,4])    # reads in the first 5 integatrions
                                   # NOT scans
            r.read(range(100))     # read in the first 100 integrations
        """
        from asap import scantable
        from asap import asaplog
        if integrations is None:
            integrations = [-1]
        asaplog.push("Reading integrations from disk...")
        sdreader._read(self,integrations)
        tbl = sdreader._getdata(self)
        sdreader._reset(self) # reset to the beginning of the file
        if self.unit is not None:
            tbl.set_fluxunit(self.unit)
        print_log()
        return scantable(tbl)

    def summary(self, name=None):
        """
        Print a summary of all scans/integrations. This reads through the
        whole file once.
        Parameters:
             None
        Example:
             r.summary()
        """
        sdreader._reset(self)
        sdreader._read(self,[-1])
        from asap import scantable
        tbl = scantable(sdreader._getdata(self))
        sdreader._reset(self)
        tbl.summary(name)
        return
