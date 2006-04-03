import _asap
from asap import rcParams
from asap import print_log

class fitter:
    """
    The fitting class for ASAP.
    """

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
            msg = "Please give a correct scan"
            if rcParams['verbose']:
                print msg
                return
            else:
                raise TypeError(msg)
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
            msg = "Invalid function type."
            if rcParams['verbose']:
                print msg
                return
            else:
                raise TypeError(msg)

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
            msg = "Fitter not yet initialised. Please set data & fit function"
            if rcParams['verbose']:
                print msg
                return
            else:
                raise RuntimeError(msg)

        else:
            if self.data is not None:
                self.x = self.data._getabcissa(row)
                self.y = self.data._getspectrum(row)
                from asap import asaplog
                asaplog.push("Fitting:")
                i = row
                out = "Scan[%d] Beam[%d] IF[%d] Pol[%d] Cycle[%d]" % (self.data.getscan(i),self.data.getbeam(i),self.data.getif(i),self.data.getpol(i), self.data.getcycle(i))
                asaplog.push(out)
        self.fitter.setdata(self.x, self.y, self.mask)
        if self.fitfunc == 'gauss':
            ps = self.fitter.getparameters()
            if len(ps) == 0:
                self.fitter.estimate()
        try:
            self.fitter.fit()
        except RuntimeError, msg:
            if rcParams['verbose']:
                print msg
            else:
                raise
        self._fittedrow = row
        self.fitted = True
        print_log()
        return

    def store_fit(self):
        """
        Store the fit parameters in the scantable.
        """
        if self.fitted and self.data is not None:
            pars = list(self.fitter.getparameters())
            fixed = list(self.fitter.getfixedparameters())
            from asap.asapfit import asapfit
            fit = asapfit()
            fit.setparameters(pars)
            fit.setfixedparameters(fixed)
            fit.setfunctions(self.fitfuncs)
            fit.setcomponents(self.components)
            fit.setframeinfo(self.data._getcoordinfo())
            self.data._addfit(fit,self._fittedrow)

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
            msg = "Please specify a fitting function first."
            if rcParams['verbose']:
                print msg
                return
            else:
                raise RuntimeError(msg)
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
        print_log()
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
            msg = "Function only operates on Gaussian components."
            if rcParams['verbose']:
                print msg
                return
            else:
                raise ValueError(msg)
        if 0 <= component < len(self.components):
            self.set_parameters([peak, centre, fhwm],
                                [peakfixed, centerfixed, fhwmfixed],
                                component)
        else:
            msg = "Please select a valid  component."
            if rcParams['verbose']:
                print msg
                return
            else:
                raise ValueError(msg)

    def get_area(self, component=None):
        """
        Return the area under the fitted gaussian component.
        Parameters:
              component:   the gaussian component selection,
                           default (None) is the sum of all components
        Note:
              This will only work for gaussian fits.
        """
        if not self.fitted: return
        if self.fitfunc == "gauss":
            pars = list(self.fitter.getparameters())
            from math import log,pi,sqrt
            fac = sqrt(pi/log(16.0))
            areas = []
            for i in range(len(self.components)):
                j = i*3
                cpars = pars[j:j+3]
                areas.append(fac * cpars[0] * cpars[2])
        else:
            return None
        if component is not None:
            return areas[component]
        else:
            return sum(areas)

    def get_parameters(self, component=None):
        """
        Return the fit paramters.
        Parameters:
             component:    get the parameters for the specified component
                           only, default is all components
        """
        if not self.fitted:
            msg = "Not yet fitted."
            if rcParams['verbose']:
                print msg
                return
            else:
                raise RuntimeError(msg)
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
        fpars = self._format_pars(cpars, cfixed, self.get_area(component))
        if rcParams['verbose']:
            print fpars
        return cpars, cfixed, fpars

    def _format_pars(self, pars, fixed, area):
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
                out += '  %2d: peak = %3.3f %s , centre = %3.3f %s, FWHM = %3.3f %s\n      area = %3.3f %s %s\n' % (c,pars[i],ounit,pars[i+1],aunit,pars[i+2],aunit, area,ounit,aunit)
                c+=1
                i+=3
        return out

    def get_estimate(self):
        """
        Return the parameter estimates (for non-linear functions).
        """
        pars = self.fitter.getestimate()
        fixed = self.fitter.getfixedparameters()
        if rcParams['verbose']:
            print self._format_pars(pars,fixed)
        return pars

    def get_residual(self):
        """
        Return the residual of the fit.
        """
        if not self.fitted:
            msg = "Not yet fitted."
            if rcParams['verbose']:
                print msg
                return
            else:
                raise RuntimeError(msg)
        return self.fitter.getresidual()

    def get_chi2(self):
        """
        Return chi^2.
        """
        if not self.fitted:
            msg = "Not yet fitted."
            if rcParams['verbose']:
                print msg
                return
            else:
                raise RuntimeError(msg)
        ch2 = self.fitter.getchi2()
        if rcParams['verbose']:
            print 'Chi^2 = %3.3f' % (ch2)
        return ch2

    def get_fit(self):
        """
        Return the fitted ordinate values.
        """
        if not self.fitted:
            msg = "Not yet fitted."
            if rcParams['verbose']:
                print msg
                return
            else:
                raise RuntimeError(msg)
        return self.fitter.getfit()

    def commit(self):
        """
        Return a new scan where the fits have been commited (subtracted)
        """
        if not self.fitted:
            print "Not yet fitted."
            msg = "Not yet fitted."
            if rcParams['verbose']:
                print msg
                return
            else:
                raise RuntimeError(msg)
        from asap import scantable
        if not isinstance(self.data, scantable):
            msg = "Not a scantable"
            if rcParams['verbose']:
                print msg
                return
            else:
                raise TypeError(msg)
        scan = self.data.copy()
        scan._setspectrum(self.fitter.getresidual())
        print_log()

    def plot(self, residual=False, components=None, plotparms=False, filename=None):
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
        if not self._p or self._p.is_dead:
            if rcParams['plotter.gui']:
                from asap.asaplotgui import asaplotgui as asaplot
            else:
                from asap.asaplot import asaplot
            self._p = asaplot()
        self._p.hold()
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

        colours = ["#777777","#bbbbbb","red","orange","purple","green","magenta", "cyan"]
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
        xlim=[min(self.x),max(self.x)]
        self._p.axes.set_xlim(xlim)
        self._p.set_axes('xlabel',xlab)
        self._p.set_axes('ylabel',ylab)
        self._p.set_axes('title',tlab)
        self._p.release()
        if (not rcParams['plotter.gui']):
            self._p.save(filename)
        print_log()

    def auto_fit(self, insitu=None):
        """
        Return a scan where the function is applied to all rows for
        all Beams/IFs/Pols.

        """
        from asap import scantable
        if not isinstance(self.data, scantable) :
            msg = "Data is not a scantable"
            if rcParams['verbose']:
                print msg
                return
            else:
                raise TypeError(msg)
        if insitu is None: insitu = rcParams['insitu']
        if not insitu:
            scan = self.data.copy()
        else:
            scan = self.data
        rows = xrange(scan.nrow())
        from asap import asaplog
        asaplog.push("Fitting:")
        for r in rows:
            out = " Scan[%d] Beam[%d] IF[%d] Pol[%d] Cycle[%d]" %        (scan.getscan(r),scan.getbeam(r),scan.getif(r),scan.getpol(r), scan.getcycle(r))
            asaplog.push(out, False)
            self.x = scan._getabcissa(r)
            self.y = scan._getspectrum(r)
            self.data = None
            self.fit()
            x = self.get_parameters()
            scan._setspectrum(self.fitter.getresidual(), r)
        print_log()
        return scan

