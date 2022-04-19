"""
decam_reduce.plotting
=====================

Plotting utilities, particularly for quality assurance purposes.

"""

import matplotlib.pyplot as plt
import numpy as np # assuming this will get used ...
import decam_reduce.common as common

def exp_sky_locations(nightsum, coordsys='equ'):
    """
    Make a full-sky plot of exposure locations.

    Parameters
    ----------
        nightsum : pandas.core.frame.DataFrame
            Night summary data frame as might be obtained from e.g., the
            util.query_night function. For coordsys='equ', needs columns
            named 'ra_min' and 'dec_min'. For coordsys='gal', needs columns
            named 'lgal', 'bgal'.
        coordsys : str, optional
            Should be either 'equ' for equatorial coordinates or 'gal'
            for Galactic coordinates.

    Notes 
    -----
        Do I need more input parameters to control things like the desired
            output directory and output image file name?
        Doesn't return anything but saves an image to disk...
        For now, do this in simple-minded rectangular projection -- use
            better projection scheme in the future?
        Different color points for different filters?
        Add overlays showing Galactic and/or ecliptic planes?

    """

    assert(len(nightsum) > 0)

    lon = nightsum['lgal'] if coordsys == 'gal' else nightsum['ra_min']
    lat = nightsum['bgal'] if coordsys == 'gal' else nightsum['dec_min']
    
    xlabel = 'Galactic longitude ' if coordsys == 'gal' else 'RA'
    ylabel = 'Galacticl latitude ' if coordsys == 'gal' else 'Dec'

    xlabel += ' (degrees)'
    ylabel += ' (degrees)'

    plt.scatter(lon, lat, s=10, edgecolor='none', c='k')

    par = common.decam_params()

    plt.xlim((362.5, -2.5))
    plt.ylim((-92.5, par['northern_limit_deg'] + 2.5))

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    caldat = nightsum['caldat'].iloc[0]

    title = caldat

    plt.title(title)

    plt.show()
