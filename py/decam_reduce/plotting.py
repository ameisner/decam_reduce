"""
decam_reduce.plotting
=====================

Plotting utilities, particularly for quality assurance purposes.

"""

import matplotlib.pyplot as plt
import numpy as np
import decam_reduce.common as common
import decam_reduce.util as util
from astropy.coordinates import SkyCoord
import astropy.units as u

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
            Need to finish this for filters other than grz.
        Add overlays showing Galactic and/or ecliptic planes?
        Options for saved image file format e.g., PNG versus EPS?
        Option to either show or save file?

    """

    assert(len(nightsum) > 0)

    lon = nightsum['lgal'] if coordsys == 'gal' else nightsum['ra_min']
    lat = nightsum['bgal'] if coordsys == 'gal' else nightsum['dec_min']
    
    xlabel = 'Galactic longitude ' if coordsys == 'gal' else 'RA'
    ylabel = 'Galacticl latitude ' if coordsys == 'gal' else 'Dec'

    xlabel += ' (degrees)'
    ylabel += ' (degrees)'

    color_dict = {util.full_filter_name('g') : 'g',
                  util.full_filter_name('r') : 'r',
                  util.full_filter_name('z') : 'm'}

    for band in color_dict.keys():
        plt.scatter(lon[nightsum['ifilter'] == band],
                    lat[nightsum['ifilter'] == band],
                    s=10, edgecolor='none', c=color_dict[band],
                    label=util.abbrev_filter_name(band))

    par = common.decam_params()

    # now plot the Galactic plane, if coordsys is 'equ'
    if coordsys == 'equ':

        for bgal in [-10, 0, 10]:
            lgal_samp = np.arange(361)
            bgal_samp = lgal_samp*0 + bgal

            skycoords = SkyCoord(lgal_samp*u.deg, bgal_samp*u.deg,
                                 frame='galactic')

            linestyle = '-' if bgal == 0 else '--'
            plt.plot(skycoords.icrs.ra, skycoords.icrs.dec, c='gray',
                     linestyle=linestyle)

    pad = 2.5 # deg
    plt.xlim((360 + pad, -1*pad))
    plt.ylim((-90 - pad, par['northern_limit_deg'] + pad))

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    caldat = nightsum['caldat'].iloc[0]

    title = caldat

    plt.title(title)

    plt.legend(loc='lower left')

    plt.show()
