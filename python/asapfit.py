from _asap import fitentry
from asap import rcParams

class asapfit(fitentry):

    def __init__(self, other=None):
        if isinstance(other, fitentry):
            fitentry.__init__(self,other)
        else:
            fitentry.__init__(self)

    def __str__(self):
        out = ""
        out += "Fit:"
        pars = self.getparameters()
        mask = self.getfixedparameters()
        funcs = self.getfunctions()
        comps = self.getcomponents()
        finfo = self.getframeinfo()
        pos=0
        k = 0
        for f in funcs:
            out += "\n Type:         "
            out += f
            s = pos
            pos += comps[k]
            ps = pars[s:pos]
            out += "\n  Parameters:  "
            out += self._format_pars(pars[s:pos],f, finfo[0])
            out += "\n  Fixed Parms: "
            out += str(mask[s:pos])
            out += "\n  Frame:       "
            out += str(finfo)
            out += "\n"
        out += "\n"
        return out

    def as_dict(self):
        out = []
        pars = self.getparameters()
        mask = self.getfixedparameters()
        funcs = self.getfunctions()
        comps = self.getcomponents()
        pos=0
        k=0
        comp = []
        for f in funcs:
            s = pos
            pos += comps[k]
            ps = pars[s:pos]
            d = {'function' : f,
                  'parameters' : pars[s:pos],
                  'fixed' : mask[s:pos],
                  'frame' : self.getframeinfo()
                  }
            comp.append(d)
        out.append(comp)
        return out

    def _format_pars(self, pars, ftype, unit):
        out = ''
        if ftype == 'poly':
            i = 0
            for p in pars:
                out += ' p%d = %3.3f %s,' % (i,p,unit)
                i+=1
            out = out[1:-1]
        elif ftype == 'gauss':
            out += 'peak = %3.3f , centre = %3.3f %s, FWHM = %3.3f %s' % (pars[0],pars[1],unit,pars[2],unit)
        return out
