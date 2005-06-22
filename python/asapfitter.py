import _asap
from asap import rcParams

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
        self.fitfuncs = None
        self.fitted = False
        self.data = None
        self.components = 0
        self._fittedrow = 0
        self._p = None
        self._vb = True
        self._selection = None

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
        n=0
        if kwargs.has_key('poly'):
            self.fitfunc = 'poly'
            n = kwargs.get('poly')
            self.components = [n]
        elif kwargs.has_key('gauss'):
            n = kwargs.get('gauss')
            self.fitfunc = 'gauss'
            self.fitfuncs = [ 'gauss' for i in range(n) ]
            self.components = [ 3 for i in range(n) ]
        else:
            print "Invalid function type."
            return
        self.fitter.setexpression(self.fitfunc,n)
        return
            
    def fit(self, row=0):
        """
        Execute the actual fitting process. All the state has to be set.
        Parameters:
            row:    specify the row in the scantable
        Example:
            s = scantable('myscan.asap')
            s.set_cursor(thepol=1)        # select second pol
            f = fitter()
            f.set_scan(s)
            f.set_function(poly=0)
            f.fit(row=0)                  # fit first row 
        """
        if ((self.x is None or self.y is None) and self.data is None) \
               or self.fitfunc is None:
            print "Fitter not yet initialised. Please set data & fit function"
            return
        else:
            if self.data is not None:
                self.x = self.data._getabcissa(row)
                self.y = self.data._getspectrum(row)
                print "Fitting:"
                vb = self.data._vb
                self.data._vb = True
                self.selection = self.data.get_cursor()
                self.data._vb = vb
        self.fitter.setdata(self.x, self.y, self.mask)
        if self.fitfunc == 'gauss':
            ps = self.fitter.getparameters()
            if len(ps) == 0:
                self.fitter.estimate()
        try:
            self.fitter.fit()
        except RuntimeError, msg:
            print msg            
        self._fittedrow = row
        self.fitted = True
        return

    def store_fit(self):
        """
        Store the fit parameters in the scantable.
        """
        if self.fitted and self.data is not None:
            pars = list(self.fitter.getparameters())
            fixed = list(self.fitter.getfixedparameters())
            self.data._addfit(self._fittedrow, pars, fixed,
                              self.fitfuncs, self.components)

    def set_parameters(self, params, fixed=None, component=None):
        """
        Set the parameters to be fitted.
        Parameters:
              params:    a vector of parameters
              fixed:     a vector of which parameters are to be held fixed
                         (default is none)
              component: in case of multiple gaussians, the index of the
                         component
             """
        if self.fitfunc is None:
            print "Please specify a fitting function first."
            return
        if self.fitfunc == "gauss" and component is not None:
            if not self.fitted:
                from numarray import zeros
                pars = list(zeros(len(self.components)*3))
                fxd = list(zeros(len(pars)))
            else:
                pars = list(self.fitter.getparameters())              
                fxd = list(self.fitter.getfixedparameters())
            i = 3*component
            pars[i:i+3] = params
            fxd[i:i+3] = fixed
            params = pars
            fixed = fxd          
        self.fitter.setparameters(params)
        if fixed is not None:
            self.fitter.setfixedparameters(fixed)
        return

    def set_gauss_parameters(self, peak, centre, fhwm,
                             peakfixed=False, centerfixed=False,
                             fhwmfixed=False,
                             component=0):
        """
        Set the Parameters of a 'Gaussian' component, set with set_function.
        Parameters:
            peak, centre, fhwm:  The gaussian parameters
            peakfixed,
            centerfixed,
            fhwmfixed:           Optional parameters to indicate if
                                 the paramters should be held fixed during
                                 the fitting process. The default is to keep
                                 all parameters flexible.
            component:           The number of the component (Default is the
                                 component 0)
        """
        if self.fitfunc != "gauss":
            print "Function only operates on Gaussian components."
            return
        if 0 <= component < len(self.components):
            self.set_parameters([peak, centre, fhwm],
                                [peakfixed, centerfixed, fhwmfixed],
                                component)
        else:
            print "Please select a valid  component."
            return
        
    def get_parameters(self, component=None):
        """
        Return the fit paramters.
        Parameters:
             component:    get the parameters for the specified component
                           only, default is all components
        """
        if not self.fitted:
            print "Not yet fitted."
        pars = list(self.fitter.getparameters())
        fixed = list(self.fitter.getfixedparameters())
        if component is not None:            
            if self.fitfunc == "gauss":
                i = 3*component
                cpars = pars[i:i+3]
                cfixed = fixed[i:i+3]
            else:
                cpars = pars
                cfixed = fixed                
        else:
            cpars = pars
            cfixed = fixed
        fpars = self._format_pars(cpars, cfixed)
        if self._vb:
            print fpars
        return cpars, cfixed, fpars
    
    def _format_pars(self, pars, fixed):
        out = ''
        if self.fitfunc == 'poly':
            c = 0
            for i in range(len(pars)):
                fix = ""
                if fixed[i]: fix = "(fixed)"
                out += '  p%d%s= %3.3f,' % (c,fix,pars[i])
                c+=1
            out = out[:-1]  # remove trailing ','
        elif self.fitfunc == 'gauss':
            i = 0
            c = 0
            aunit = ''
            ounit = ''
            if self.data:
                aunit = self.data.get_unit()
                ounit = self.data.get_fluxunit()
            while i < len(pars):
                out += '  %d: peak = %3.3f %s , centre = %3.3f %s, FWHM = %3.3f %s \n' % (c,pars[i],ounit,pars[i+1],aunit,pars[i+2],aunit)
                c+=1
                i+=3
        return out
        
    def get_estimate(self):
        """
        Return the parameter estimates (for non-linear functions).
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
        Return a new scan where the fits have been commited (subtracted)
        """
        if not self.fitted:
            print "Not yet fitted."
        if self.data is not scantable:
            print "Only works with scantables"
            return
        scan = self.data.copy()
        scan._setspectrum(self.fitter.getresidual())

    def plot(self, residual=False, components=None, plotparms=False):
        """
        Plot the last fit.
        Parameters:
            residual:    an optional parameter indicating if the residual
                         should be plotted (default 'False')
            components:  a list of components to plot, e.g [0,1],
                         -1 plots the total fit. Default is to only
                         plot the total fit.
            plotparms:   Inidicates if the parameter values should be present
                         on the plot
        """
        if not self.fitted:
            return
        if not self._p:
            from asap.asaplot import ASAPlot
            self._p = ASAPlot()
        if self._p.is_dead:
            from asap.asaplot import ASAPlot
            self._p = ASAPlot()
        self._p.clear()
        self._p.set_panels()
        self._p.palette(0)
        tlab = 'Spectrum'
        xlab = 'Abcissa'        
        m = ()
        if self.data:
            tlab = self.data._getsourcename(self._fittedrow)
            xlab = self.data._getabcissalabel(self._fittedrow)
            m = self.data._getmask(self._fittedrow)
            ylab = self.data._get_ordinate_label()

        colours = ["grey60","grey80","red","orange","purple","green","magenta", "cyan"]
        self._p.palette(0,colours)
        self._p.set_line(label='Spectrum')
        self._p.plot(self.x, self.y, m)
        if residual:
            self._p.palette(1)
            self._p.set_line(label='Residual')
            self._p.plot(self.x, self.get_residual(), m)
        self._p.palette(2)
        if components is not None:
            cs = components
            if isinstance(components,int): cs = [components]
            if plotparms:
                self._p.text(0.15,0.15,str(self.get_parameters()[2]),size=8)
            n = len(self.components)
            self._p.palette(3)
            for c in cs:
                if 0 <= c < n:
                    lab = self.fitfuncs[c]+str(c)
                    self._p.set_line(label=lab)
                    self._p.plot(self.x, self.fitter.evaluate(c), m)
                elif c == -1:
                    self._p.palette(2)
                    self._p.set_line(label="Total Fit")
                    self._p.plot(self.x, self.get_fit(), m)                    
        else:
            self._p.palette(2)
            self._p.set_line(label='Fit')
            self._p.plot(self.x, self.get_fit(), m)
        self._p.set_axes('xlabel',xlab)
        self._p.set_axes('ylabel',ylab)
        self._p.set_axes('title',tlab)
        self._p.release()

    def auto_fit(self, insitu=None):
        """
        Return a scan where the function is applied to all rows for
        all Beams/IFs/Pols.
        
        """
        from asap import scantable
        if not isinstance(self.data, scantable) :
            print "Only works with scantables"
            return
        if insitu is None: insitu = rcParams['insitu']
        if not insitu:
            scan = self.data.copy()
        else:
            scan = self.data
        vb = scan._vb
        scan._vb = False
        sel = scan.get_cursor()
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
                        self.x = scan._getabcissa(iRow)
                        self.y = scan._getspectrum(iRow)
                        self.data = None
                        self.fit()                    
                        x = self.get_parameters()
                        scan._setspectrum(self.fitter.getresidual(),iRow)
        scan.set_cursor(sel[0],sel[1],sel[2])
        scan._vb = vb
        if not insitu:
            return scan
