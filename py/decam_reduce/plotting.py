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
from pkg_resources import resource_filename
import os
import glob
import astropy.io.fits as fits

def exp_sky_locations(nightsum, coordsys='equ', save=False):
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
        save : bool, optional
            If true, save to file instead of popping up the plot.

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
    lat_upper = par['northern_limit_deg'] if coordsys is 'equ' else 90

    plt.xlim((360 + pad, -1*pad))
    plt.ylim((-90 - pad, lat_upper + pad))

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    caldat = nightsum['caldat'].iloc[0]

    title = caldat

    plt.title(title)

    plt.legend(loc='lower left')

    if not save:
        plt.show()
    else:
        outname = 'exp_sky_locations_' + caldat + '-' + coordsys + '.png'
        plt.savefig(outname)

def outputs_fp_map(fname_raw, rerun_dir, save=False):
    """
    Focal plane map of which CCDs were reduced successfully.
    """

    assert(os.path.exists(fname_raw))

    ccds, ccdnums = util.get_raw_ccds_list(fname_raw)

    assert(os.path.exists(rerun_dir))

    expid = util.get_expid(fname_raw)

    expid_str = str(expid).zfill(7)
    dir = os.path.join(rerun_dir, expid_str, 'calexp')

    flist = glob.glob(dir + '/calexp-' + expid_str + '_??.fits')

    ccdnum_calexp = [int(f[-7:-5]) for f in flist]

    # read in the fits file with the CCD boundaries
    fname_corners = resource_filename('decam_reduce',
        os.path.join('data', 'cornerCoords_SN-C1-reordered-TAN.fits'))

    corners = fits.getdata(fname_corners)

    # plot the edges ; red for missing, black otherwise

    for _ccdnum in np.unique(corners['CCD']):
        _corners = corners[corners['CCD'] == _ccdnum]
        missing = (_ccdnum in ccdnums) and (_ccdnum not in ccdnum_calexp)
        color = 'r' if missing else 'k'
        ra_center = np.mean(_corners['RA'])
        dec_center = np.mean(_corners['DEC'])
        _corners = _corners[[0, 1, 2, 3, 0]]
        plt.plot(_corners['RA'], _corners['DEC'], c=color)
        plt.text(ra_center, dec_center, str(_ccdnum), color=color)

    plt.text(np.max(corners['RA']) - 0.2,
             np.min(corners['DEC']) - 0.08, 
             'black = reductions succeeded')

    plt.text(np.max(corners['RA']) - 1.4,
             np.min(corners['DEC']) - 0.08, 
             'red = reductions missing', color='r')

    plt.xlim((np.max(corners['RA']) + 0.05, np.min(corners['RA']) - 0.05))
    plt.xticks([])
    plt.yticks([])
    plt.title('EXPNUM = ' + str(expid))

    if not save:
        plt.show()
    else:
        outname = 'fp_output_summary-' + str(expid).zfill(7) + '.png'
        plt.savefig(outname, dpi=200, bbox_inches='tight')
 
    
