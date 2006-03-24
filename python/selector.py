from asap._asap import selector as _selector
class selector(_selector):
    def __init(self):
        _selector.__init__(self)

    def reset(self):
        """
        Unset all selections.
        """
        self._reset()

        
    def set_polarisations(self, pols=[]):
        """
        Set the polarisations to be selected in the scantable.
        Parameters:
             pols:     a list of integers of 0-3, or strings, e.g ["I","Q"].
                       Default [] is no selection
        Example:
             sel = selector()
             # These are equivalent if data is 'linear'
             sel.set_polarisations(["XX","Re(XY)"])
             sel.set_polarisations([0,2])
             # reset the polarisation selection
             sel.set_polarisations()
                          
        """
        vec = _to_list(pols, str) and _to_list(pols, int)
        if vec: # is an empty and/or valid vector
            if len(vec) and isinstance(vec[-1],str):
                self._setpolstrings(self, vec)
                return
            self._setpols(vec)
        else:
            raise TypeError('Unknown pol type. Please use [0,1...] or ["XX","YY"...]')

    def set_ifs(self, ifs=[]):
        """
        Set a sequence of IF numbers (0-based).
        Parameters:
            ifs:    a list of integers. Default [] is to unset the selection.
        """
        vec = _to_list(ifs, int)
        if vec:
            self._setifs(vec)
        else:
            raise TypeError('Unknown IFno type. Use lists of integers.')

    def set_scans(self, scans=[]):
        """
        Set a sequence of Scan numbers (0-based).
        Parameters:
            scans:    a list of integers. Default [] is to unset the selection.
        """
        vec = _to_list(scans, int)
        if vec:
            self._setscans(vec)
        else:
            raise TypeError('Unknown Scan number type. Use lists of integers.')

    def set_beams(self, beams=[]):
        """
        Set a sequence of Beam numbers (0-based).
        Parameters:
            beams:    a list of integers. Default [] is to unset the selection.
        """
        vec = _to_list(beams, int)
        if vec:
            self._setbeams(vec)
        else:
            raise TypeError('Unknown Beam number type. Use lists of integers.')

    def set_cycles(self, cycles=[]):
        """
        Set a sequence of IF numbers (0-based).
        Parameters:
            cycless:    a list of integers. Default [] is to unset the selection.
        """
        vec = _to_list(cycles, int)
        if vec:
            self._setcycles(vec)
        else:
            raise TypeError('Unknown Cycle number type. Use lists of integers.')
                                    

    def set_name(self, name):
        """
        Set a selection based on a name. This can be a unix pattern , e.g. "*_R"
        Parameters:
            name:    a string containing a source name or pattern
        Examples:
            # select all reference scans which start with "Orion"
            selection.set_name("Orion*_R")
        """
        if isinstance(name, str):
            self._setame(name)
        else:
            raise TypeError('name must be a string')
        
    def set_tsys(self, tsysmin=0.0, tsysmax=None):
        """
        Select by Tsys range.
        Parameters:
            tsysmin:     the lower threshold. Default 0.0
            tsysmax:     the upper threshold. Default None.
        Examples:
            # select all spectra with Tsys <= 500.0
            selection.set_tsys(tsysmax=500.0)
            
        """
        taql =  "SELECT FROM $1 WHERE TSYS >= %f" % (tsysmin)
        if isinstance(tsysmax, float):
            taql = taql + " AND TSYS <= %f" % ( tsysmax)
        self._settaql(taql)

    def set_query(self, query):
        """
        Select by Column query. Power users only!
        Example:
            # select all off scans with integration times over 60 seconds.
            selection.set_query("SRCTYPE == 1 AND INTERVAL > 60.0")
        """
        taql = "SELECT FROM $1 WHERE " + query
        self._settaql(taql)
        
