from _asap import sdfit
from asap import rcParams

class asapfit(sdfit):
    
    def __init__(self, other):        
        sdfit.__init__(self,other)

    def __str__(self):
        if self.__len__() == 0:
            return "No fits"
        out = ""
        for i in range(self.__len__()):
            out += "Fit No %d:" % (i)
            pars = self.getparameters(i)
            mask = self.getfixedparameters(i)
            funcs = self.getfunctions(i)
            comps = self.getcomponents(i)
            finfo = self.getframeinfo(i)
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
        for i in range(self.__len__()):
            pars = self.getparameters(i)
            mask = self.getfixedparameters(i)
            funcs = self.getfunctions(i)
            comps = self.getcomponents(i)
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
                     'frame' : self.getframeinfo(i)
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
