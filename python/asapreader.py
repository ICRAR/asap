from asap._asap import sdreader

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
        sdreader.__init__(self, filename, theif, thebeam)

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
        if integrations is None:
            integrations = [-1]
        print "Reading integrations from disk..."
        sdreader._read(self,integrations)
        tbl = sdreader._getdata(self)
        if self.unit is not None:
            tbl.set_fluxunit(self.unit)
        return scantable(tbl)

    def summary(self):
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
        tbl = sdreader._getdata(self)
        sdreader._reset(self)
        print tbl._summary()
        return
    def reset(self):
        """
        Reset to the beginning of the file.
        Parameters:
             none
        Examples:
             r = reader('xyz.sdfits')
             scans0 = r.read(range(10))
             r.reset()
             scans1 = r.read(range(5,15))
             # produces 2 scatables  one containing integrations(rows) 0-9
             # and a second one with rows 5-14
        """
        sdreader._reset(self)
        
