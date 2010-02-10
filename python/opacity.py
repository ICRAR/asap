import os
import math
from asap import scantable
from asap import merge
from asap import fitter
from asap import selector
from asap import rcParams

def _import_data(data):
    if not hasattr(data, "__len__"):
        if isinstance(data, scantable):
            return data
        elif isinstance(data, str):
            return scantable(data)
    tables = []
    for d in data:
        if isinstance(d, scantable):
            tables.append(d)
        elif isinstance(d, str):
            if os.path.exists(d):
                tables.append(scantable(d))
            else:
                raise IOError("Data file doesn't exists")
        else:
            raise TypeError("data is not a scantable or valid file")
    return merge(tables)

def skydip(data, averagepol=True, tsky=300., plot=False):
    """Determine the opacity from a set of 'skydip' obervations.
    This can be any set of observations over a range of elevations,
    but will ususally be a dedicated (set of) scan(s).
    Return a list of 'n' opacities for 'n' IFs. In case of averagepol
    being 'False' a list of 'n*m' elements where 'm' is the number of
    polarisations, e.g.
    nIF = 3, nPol = 2 => [if0pol0, if0pol1, if1pol0, if1pol1, if2pol0, if2pol1]

    The opacity is determined by fitting a first order polynomial to:


        Tsys(airmass) = p0 + airmass*p1

    where

        airmass = 1/sin(elevation)

        tau =  p1/Tsky

    Parameters:
        data:       a list of file names or scantables or a single
                    file name or scantable.
        averagepol: Return the average of the opacities for the polarisations
                    This might be useful to set to 'False' if one polarisation
                    is corrupted (Mopra). If set to 'False', an opacity value
                    per polarisation is returned.
        tksy:       The sky temperature (default 300.0K). This might
                    be read from the data in the future.
        plot:       Plot each fit (airmass vs. Tsys). Default is 'False'
    """
    scan = _import_data(data)
    f = fitter()
    f.set_function(poly=1)
    sel = selector()
    basesel = scan.get_selection()
    inos = scan.getifnos()
    pnos = scan.getpolnos()
    opacities = []
    for ino in inos:
        sel.set_ifs(ino)
        opacity = []
        for pno in pnos:
            sel.set_polarisations(pno)
            scan.set_selection(basesel+sel)
            tsys = scan.get_tsys()
            elev = scan.get_elevation()
            airmass = [ 1./math.sin(i) for i in elev ]
            f.set_data(airmass, tsys)
            f.fit()
            if plot:
                f.plot(residual=True, plotparms=True)
                raw_input("Hit <return> for next fit...")
            params = f.get_parameters()["params"]
            opacity.append(params[1]/tsky)
        if averagepol:
            opacities.append(sum(opacity)/len(opacity))
        else:
            opacities += opacity
        sel.reset()
    scan.set_selection(basesel)
    return opacities
