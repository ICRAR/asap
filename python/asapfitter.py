import _asap

class fitter:
    """
    The fitting class for ASAP.
    """
    def _verbose(self, *args):
        """
        Set stdout output.
        """
        if type(args[0]) is bool:
            self._vb = args[0]
            return
        elif len(args) == 0:
            return self._vb
        
    def __init__(self):
        """
        Create a fitter object. No state is set.
        """
        self.fitter = _asap.fitter()
        self.x = None
        self.y = None
        self.mask = None
        self.fitfunc = None
        self.fitted = False
        self.data = None
        self._p = None
        self._vb = True

    def set_data(self, xdat, ydat, mask=None):
        """
        Set the absissa and ordinate for the fit. Also set the mask
        indicationg valid points.
        This can be used for data vectors retrieved from a scantable.
        For scantable fitting use 'fitter.set_scan(scan, mask)'.
        Parameters:
            xdat:    the abcissa values
            ydat:    the ordinate values
            mask:    an optional mask
        
        """
        self.fitted = False
        self.x = xdat
        self.y = ydat
        if mask == None:
            from numarray import ones
            self.mask = ones(len(xdat))
        else:
            self.mask = mask
        return

    def set_scan(self, thescan=None, mask=None):
        """
        Set the 'data' (a scantable) of the fitter.
        Parameters:
            thescan:     a scantable
            mask:        a msk retireved from the scantable
        """
        if not thescan:
            print "Please give a correct scan"
        self.fitted = False
        self.data = thescan
        if mask is None:
            from numarray import ones
            self.mask = ones(self.data.nchan())
        else:
            self.mask = mask
        return

    def set_function(self, **kwargs):
        """
        Set the function to be fit.
        Parameters:
            poly:    use a polynomial of the order given
            gauss:   fit the number of gaussian specified
        Example:
            fitter.set_function(gauss=2) # will fit two gaussians
            fitter.set_function(poly=3)  # will fit a 3rd order polynomial
        """
        #default poly order 0
        self.fitfunc = 'poly'
        n=0
        if kwargs.has_key('poly'):
            self.fitfunc = 'poly'
            n = kwargs.get('poly')
        elif kwargs.has_key('gauss'):
            n = kwargs.get('gauss')
            self.fitfunc = 'gauss'
        
        self.fitter.setexpression(self.fitfunc,n)
        return
            
    def fit(self):
        """
        Execute the actual fitting process. All the state has to be set.
        Parameters:
            none
        Example:
            s= scantable('myscan.asap')
            f = fitter()
            f.set_scan(s)
            f.set_function(poly=0)
            f.fit()
        """
        if ((self.x is None or self.y is None) and self.data is None) \
               or self.fitfunc is None:
            print "Fitter not yet initialised. Please set data & fit function"
            return
        else:
            if self.data is not None:
                self.x = self.data.getabcissa()
                self.y = self.data.getspectrum()
                print "Fitting:"
                vb = self.data._verbose
                self.data._verbose(True)
                s = self.data.get_selection()
                self.data._verbose(vb)
        
        self.fitter.setdata(self.x,self.y,self.mask)
        if self.fitfunc == 'gauss':
            ps = self.fitter.getparameters()
            if len(ps) == 0:
                self.fitter.estimate()
        self.fitter.fit()
        self.fitted = True
        return

    def set_parameters(self, params, fixed=None):
        self.fitter.setparameters(params)
        if fixed is not None:
            self.fitter.setfixedparameters(fixed)
        return
    
    def get_parameters(self):
        """
        Return the fit paramters.
        
        """
        if not self.fitted:
            print "Not yet fitted."
        pars = list(self.fitter.getparameters())
        fixed = list(self.fitter.getfixedparameters())
        if self._vb:
            print self._format_pars(pars)
        return pars,fixed
    
    def _format_pars(self, pars):
        out = ''
        if self.fitfunc == 'poly':
            c = 0
            for i in pars:
                out += '  p%d = %3.3f, ' % (c,i)
                c+=1
        elif self.fitfunc == 'gauss':
            i = 0
            c = 0
            unit = ''
            if self.data:
                unit = self.data.get_unit()
            while i < len(pars):
                out += '  %d: peak = %3.3f , centre = %3.3f %s, FWHM = %3.3f %s \n' % (c,pars[i],pars[i+1],unit,pars[i+2],unit)
                c+=1
                i+=3
        return out
        
    def get_estimate(self):
        """
        Return the paramter estimates (for non-linear functions).
        """
        pars = self.fitter.getestimate()
        if self._vb:
            print self._format_pars(pars)
        return pars
       

    def get_residual(self):
        """
        Return the residual of the fit.
        """
        if not self.fitted:
            print "Not yet fitted."
        return self.fitter.getresidual()

    def get_chi2(self):
        """
        Return chi^2.
        """
        
        if not self.fitted:
            print "Not yet fitted."
        ch2 = self.fitter.getchi2()
        if self._vb:
            print 'Chi^2 = %3.3f' % (ch2)
        return ch2 

    def get_fit(self):
        """
        Return the fitted ordinate values.
        """
        if not self.fitted:
            print "Not yet fitted."
        return self.fitter.getfit()

    def commit(self):
        """
        Return a new scan where teh fits have been commited.
        """
        if not self.fitted:
            print "Not yet fitted."
        if self.data is not scantable:
            print "Only works with scantables"
            return
        scan = self.data.copy()
        scan.setspectrum(self.fitter.getresidual())

    def plot(self, residual=False):
        """
        Plot the last fit.
        Parameters:
            residual:    an optional parameter indicating if the residual
                         should be plotted (default 'False')
        """
        if not self.fitted:
            return
        if not self._p:
            from asap.asaplot import ASAPlot
            self._p = ASAPlot()
        self._p.clear()
        tlab = 'Spectrum'
        xlab = 'Abcissa'
        if self.data:
            tlab = self.data._getsourcename(0)
            xlab = self.data.getabcissalabel(0)
        ylab = r'Flux'
        m = self.data.getmask(0)
        self._p.set_line(colour='blue',label='Spectrum')
        self._p.plot(self.x, self.y, m)
        if residual:
            self._p.set_line(colour='green',label='Residual')
            self._p.plot(self.x, self.get_residual(), m)
        self._p.set_line(colour='red',label='Fit')
        self._p.plot(self.x, self.get_fit(), m)
        
        self._p.set_axes('xlabel',xlab)
        self._p.set_axes('ylabel',ylab)
        self._p.set_axes('title',tlab)
        self._p.release()


    def auto_fit(self):
        """
        Return a scan where the function is applied to all rows for all Beams/IFs/Pols.
        
        """
        from asap import scantable
        if not isinstance(self.data,scantable) :
            print "Only works with scantables"
            return
        scan = self.data.copy()
        vb = scan._verbose
        scan._verbose(False)
        sel = scan.get_selection()
        rows = range(scan.nrow())
        for i in range(scan.nbeam()):
            scan.setbeam(i)
            for j in range(scan.nif()):
                scan.setif(j)
                for k in range(scan.npol()):
                    scan.setpol(k)
                    if self._vb:
                        print "Fitting:"
                        print 'Beam[%d], IF[%d], Pol[%d]' % (i,j,k)
                    for iRow in rows:
                        self.x = scan.getabcissa(iRow)
                        self.y = scan.getspectrum(iRow)
                        self.data = None
                        self.fit()                    
                        x = self.get_parameters()
                        scan.setspectrum(self.fitter.getresidual(),iRow)
        scan.set_selection(sel[0],sel[1],sel[2])
        scan._verbose(vb)
        return scan
