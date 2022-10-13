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
from astropy.table import Table

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

def outputs_fp_map(fname_raw, rerun_dir, save=False, outname_extra='',
                   title_extra=''):
    """
    Focal plane map of which CCDs were/weren't reduced successfully.

    Parameters
    ----------
        fname_raw : str
            File name of raw DECam image.
        rerun_dir : str
            Base rerun directory within which to look for calexp outputs.
        save : bool, optional
            If set True, save the plot to a PNG file. Default is False, in
            which case the plot is popped up rather than saved.
        outname_extra : str, optional
            Add-on to the output PNG filename that will appear immediately
            before the .png extension. Empty by default.
        title_extra : str, optional
            Add-on to plot title. Empty by default.

    Notes
    -----
        Might be interesting to add a boolean keyword argument that dictates
        whether the CCD labels are CCDNUM or CCDNAME.

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
        _annot = util.ccdnum_to_ccdname(_ccdnum) + ' / ' + str(_ccdnum)
        plt.text(ra_center + 0.15, dec_center - 0.03, _annot, color=color)

    plt.text(np.max(corners['RA']) - 0.2,
             np.min(corners['DEC']) - 0.08, 
             'black = reductions succeeded')

    plt.text(np.max(corners['RA']) - 1.4,
             np.min(corners['DEC']) - 0.08, 
             'red = reductions missing', color='r')

    plt.xlim((np.max(corners['RA']) + 0.05, np.min(corners['RA']) - 0.05))
    plt.xticks([])
    plt.yticks([])

    title = 'EXPNUM = ' + str(expid)
    if title_extra != '':
        title += ' ; ' + title_extra
    plt.title(title)

    if not save:
        plt.show()
    else:
        outname = 'fp_output_summary-' + str(expid).zfill(7) + \
            outname_extra + '.png'
        plt.savefig(outname, dpi=200, bbox_inches='tight')
 
    
def zeropoint_trend(caldat, rerun_dir, _filter, save=False,
                    return_table=False):
    """
    Nightly strip chart of the LSST pipeline measured photometric zeropoint.

    Parameters
    ----------
        caldat : str
            Observing night as a string in the format YYYY-MM-DD.
        rerun_dir : str
            Relative or full path of the specific LSST pipeline rerun to use.
        _filter : str
            Filter name.
        save: bool, optional
            Save the resulting plot to PNG? Default is to not save the plot.
        return_table: bool, optional
            Return a summary table of the zeropoints and corresponding time
            stamps? Default is to not return such a table even though it gets
            generated internally as part of the plot-making process.

    Returns
    -------
        t : astropy.table.table.Table
            Table of zeropoint plus time stamp information. Only returned if
            return_table=True. By default, return_table=False and nothing is
            returned.

    """

    assert(os.path.exists(rerun_dir))

    flist = glob.glob(rerun_dir + '/*/metadata/metadata*.yaml')

    flist.sort()

    zeropoints = []
    expids = []

    for f in flist:
        _f = open(f, 'r')
        lines = _f.readlines()
        line = lines[2575]
        line = line.replace('\n', '')
        line = line.replace(' - ', '')
        line = line.replace(' ', '')
        zeropoints.append(float(line))
        _f.close()

        expid = int(os.path.basename(f)[9:16])
        expids.append(expid)


    expids = np.array(expids)
    plt.scatter(expids - np.min(expids), zeropoints, s=10)
    plt.plot(expids - np.min(expids), zeropoints)
    plt.ylabel('MAGZERO')
    plt.xlabel('EXPNUM - ' + str(np.min(expids)))
    plt.title(caldat + ' ; filter = ' + _filter + ' ; CCDNUM = 10') # HACK !

    if not save:
        plt.show()
    else:
        outname = 'MAGZERO_' + caldat + '_' + _filter + '.png'
        plt.savefig(outname, dpi=200, bbox_inches='tight')

    if return_table:
        t = Table()
        t['expid'] = expids
        t['zeropoints'] = zeropoints
        return t

def ref_star_locations(visit, ccdnum, repo_path, icSrc=False):
    """
    Overplot locations of crossmatched reference stars on top of all detections

    Parameters
    ----------
        visit : int
            Exposure number.
        ccdnum : int
            CCDNUM CCD identifier.
        repo_path : str
            Top-level output directory for the relevant rerun output.
        icSrc : bool
            Set True to use src/ catalog output rather than icSrc/ catalog
            output for the 'detection list'. Default value is False.

    Notes
    -----
        Intended to show one CCD at a time...

        Regarding src/ versus icSrc/, see:
            https://pipelines.lsst.io/v/v19_0_0/modules/lsst.pipe.tasks/tasks/lsst.pipe.tasks.processCcd.ProcessCcdTask.html?highlight=icsrc#output-datasets

        Probably this should all be done via Butler rather than constructing
        various file/directory names.

    """

    assert(os.path.exists(repo_path))

    expnum_str = str(visit).zfill(7)
    ccdnum_str = str(ccdnum).zfill(2)

    det_subdir = 'icSrc' if icSrc else 'src'
    ref_subdir = 'srcMatch'

    det_path = os.path.join(repo_path, expnum_str, det_subdir)
    
    print(det_path)
    assert(os.path.exists(det_path))

    fname_det = det_subdir + '-' + expnum_str + '_' + ccdnum_str + '.fits'

    fname_det = os.path.join(det_path, fname_det)

    print(fname_det)
    assert(os.path.exists(fname_det))

    det = fits.getdata(fname_det) # ex = 1

    ref_path = os.path.join(repo_path, expnum_str, ref_subdir)

    assert(os.path.exists(ref_path))

    fname_ref = ref_subdir + '-' + expnum_str + '_' + ccdnum_str + '.fits'
    fname_ref = os.path.join(ref_path, fname_ref)

    assert(os.path.exists(fname_ref))
    
    ref = fits.getdata(fname_ref)

    plt.scatter(det['BASE_SDSSCENTROID_X'],
                det['BASE_SDSSCENTROID_Y'],
                s=1, edgecolor='none', c='k')

    intersect, ind_ref, ind_det = np.intersect1d(ref['SECOND'], det['ID'],
        return_indices=True)

    det_row_matched = det[ind_det]

    print(ref.dtype)
    plt.scatter(det_row_matched['BASE_SDSSCENTROID_X'],
                det_row_matched['BASE_SDSSCENTROID_Y'],
                s=10, marker='o', facecolor='none', edgecolor='red')

    par = common.decam_params()

    plt.xlim([0, par['nx_active']])
    plt.ylim([0, par['ny_active']])

    plt.axes().set_aspect('equal')
    plt.show()

    # construct path of the detection catalog
    # check that detection catalog exists
    # construct path of the srcMatch matched detection-reference catalog
    # check that the srcMatch catalog exists

    # read in the two catalogs

    # plot the detection catalog as small black dots overplot the detection 
    # overplot the detectionr-reference cross-matches as larger red squares
