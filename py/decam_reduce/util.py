"""
decam_reduce.util
=================

A collection of DECam reduction related utility functions.

"""

import decam_reduce.common as common
from astropy.table import Table
import requests
import json
import pandas as pd
from astropy.table import vstack
import os
import numpy as np
import stat

def print_hostname():
    """
    Print name of host computer.

    Notes
    -----
        Meant to be a minimal utility for logging purposes...

    """

    try:
        print('Running on host: ' + str(os.environ.get('HOSTNAME')))
    except:
        print('Could not retrieve hostname!')

def print_process_id():
    """
    Print procesc ID.

    Notes
    -----
        For logging purposes...

    """

    print('PID = ' + str(os.getpid()))

def full_filter_name(filter):
    '''
    Retrieve full DECam filter name used by Astro Data Archive.

    Parameters
    ----------
        filter : str
            Abbreviated filter name; of 'u', 'g', 'r', 'i', 'z'
            'Y', 'VR'

    Returns
    -------
        str
            Full filter name.

    '''

    full_names = {'u' : 'u DECam c0006 3500.0 1000.0',
                  'g' : 'g DECam SDSS c0001 4720.0 1520.0',
                  'r' : 'r DECam SDSS c0002 6415.0 1480.0',
                  'i' : 'i DECam SDSS c0003 7835.0 1470.0',
                  'z' : 'z DECam SDSS c0004 9260.0 1520.0',
                  'Y' : 'Y DECam c0005 10095.0 1130.0',
                  'VR' : 'VR DECam c0007 6300.0 2600.0'}

    return full_names[filter]

def abbrev_filter_name(filter):
    '''
    Retrieve abbreviated DECam filter name based on full DECam filter name.

    Parameters
    ----------
        filter : str
            Full DECam filter name.

    Returns
    -------
        str
            Abbreviated filter name.

    '''

    abbrev_names = {'u DECam c0006 3500.0 1000.0' : 'u',
                    'g DECam SDSS c0001 4720.0 1520.0' : 'g',
                    'r DECam SDSS c0002 6415.0 1480.0' : 'r',
                    'i DECam SDSS c0003 7835.0 1470.0' : 'i',
                    'z DECam SDSS c0004 9260.0 1520.0' : 'z',
                    'Y DECam c0005 10095.0 1130.0' : 'Y',
                    'VR DECam c0007 6300.0 2600.0' : 'VR'}

    return abbrev_names[filter]

def getShards(ra, dec, radius, depth=7):
        '''
        Get HTM trixel numbers within a given radius of a sky location.

        Parameters
        ----------
            ra : float
                Right Ascension in decimal degrees.
            dec : float
                Declination in decimal degrees.
            radius : float
                Radius within which to return nearby HTM trixel numbers.
            depth : int (optional)
                HTM pixelization depth parameter.

        Returns
        -------
            t : astropy.table.table.Table
                Table with the list of HTM trixel numbers and corresponding
                boolean flag indicating whether each trixel is on the boundary
                of the search region.

        Notes
        -----
            adapted from:
            https://tigress-web.princeton.edu/~pprice/ps1_pv3_3pi_20170110/README.txt

            takes ~0.001-0.002 seconds for a ~1 deg radius

            takes ~7 seconds for a 180 deg radius (all-sky), almost all of
            which is accounted for by my reformatting to Astropy table. actual
            calculations take ~5e-4 seconds for radius = 1 deg, 3e-2 seconds
            for radius = 180 deg. almost all of the run time is accounted for
            by the conversion of the onBoundary generator to a list.

        '''

        from lsst.meas.algorithms.htmIndexer import HtmIndexer
        from lsst.geom import degrees
        from lsst.geom import SpherePoint

        htm = HtmIndexer(depth=depth)

        coord = SpherePoint(ra, dec, degrees)

        shards, onBoundary = htm.getShardIds(coord, radius*degrees)

        t = Table()
        t['shard'] = shards
        t['is_on_boundary'] = list(onBoundary)

        return t

def getShards_decam_pointing(ra, dec, depth=7, margin=0.0):
    '''
    Get a list of HTM trixels possibly overlapping with a DECam pointing

    Parameters
    ----------
        ra : float
            DECam pointing center RA.
        dec : float
            DECam pointing center Dec.
        depth : int
            HTM pixelization depth parameter.
        margin : float
            Amount of padding for DECam field radius in degrees.

    Returns
    -------
        result : astropy.table.table.Table
            Table with the list of HTM trixel numbers and corresponding
            boolean flag indicating whether each trixel is on the boundary
            of the search region.

    '''

    # maximum angular radius of any DECam detector pixel from the field center
    # eventually should be factored out into some sort of repository for
    # for 'special numbers', rather than hardcoded here (and elsewhere)...

    par = common.decam_params()

    radius = par['max_pixel_radius_deg']

    radius += margin

    result = getShards(ra, dec, radius, depth=7)

    return result

def query_night(night):
    """
    Issue Astro Data Archive query to get a nightly list of DECam images.

    Parameters
    ----------
        night : str
            Observing night in format YYYY-MM-DD.

    Notes
    -----
        Queries the Astro Data Archive /short API.
        Handling of failed queries? Retry in that case?

    """

    par = common.decam_params()
    url_short =  par['base_url_short'] + night + '/'
    nightsum = pd.DataFrame(requests.get(url_short).json()[1:])

    return nightsum

def select_raw_science(nightsum, min_exptime_s=None, _filter=None, propid=None):
    """
    Select raw science frames.

    Parameters
    ----------
        nightsum : pandas.core.frame.DataFrame
            Night summary DataFrame of the sort that would be returned by
            the query_night function. Needs to have columns 'proc_type',
            'prod_type', 'obs_type', 'exposure', 'original_filename'.
        min_exptime_s : float
            Minimum exposure time in seconds to consider reducing. Default
            value is obtained from the set of parameters in common.py.
        _filter : str, optional
            Should be one of ['u', 'g', 'r', 'i', 'z', 'Y', 'VR']. Underscore 
            in keyword argument name because filter is apparently a Python
            built-in (?). Default is None, in which case no cuts will be
            applied to the set of raw science images to reduce.
        propid : str, optional
            Proposal ID as a string. DECaLS example is 2014B-0404. Still need
            to investigate the set of possible/allowed proposal ID values so
            that I can eventually add some checks. Default of None means no
            filtering on proposal ID.

    Returns
    -------
        result : pandas.core.frame.DataFrame
            Filtered version of input DataFrame.

    Notes
    -----
        Restrict to u, g, r, i, z, Y ?
        What is the right quantity (if any) to sort the output by?

    """

    par = common.decam_params()

    # if user really wants no min exptime, they should specify either
    # a value <= 0 rather than None...
    if min_exptime_s is None:
        min_exptime_s = par['min_exptime_s']

    keep = (nightsum['proc_type'] == 'raw') & \
           (nightsum['prod_type'] == 'image') & \
           (nightsum['obs_type'] == 'object') & \
           (nightsum['exposure'] >= 1)

    if _filter is not None:
        _full_filter_name = full_filter_name(_filter)
        keep = np.logical_and(keep, nightsum['ifilter'] == _full_filter_name)

    if propid is not None:
        keep = np.logical_and(keep, nightsum['proposal'] == propid)

    # what to do for edge case in which nothing is retained?
    result = nightsum[keep]

    result = result.sort_values('original_filename')

    return result

def select_mastercal(nightsum, raw=None):
    """
    Select the master calibration records for an observing night.

    Parameters
    ----------
        nightsum : pandas.core.frame.DataFrame
            Observing night summary DataFrame from /short API query.
        raw : pandas.core.frame.DataFrame, optional
            Optional dataframe listing the raw science data to be processed.
            If specified, will be used to limit the calibrations to only
            those in the relevant band(s).

    Returns
    -------
        result : pandas.core.frame.DataFrame
            Subset of rows of input DataFrame corresponding to mastercal
            dome flats and zeros.

    """

    keep = (nightsum['proc_type'] == 'mastercal') & \
           (nightsum['prod_type'] == 'image') & \
           ((nightsum['obs_type'] == 'dome flat') | \
            (nightsum['obs_type'] == 'zero'))

    result = nightsum[keep]

    if raw is not None:
        filters = np.unique(raw['ifilter'])

    keep = np.ones(len(result), dtype=bool)
    for i in range(len(result)):
        if result['obs_type'].iloc[i] == 'dome flat':
            if result['ifilter'].iloc[i] not in filters:
                keep[i] = False

    result = result[keep]

    return result

def download_images(df, outdir):
    """
    Download a list of images based on their URL's.

    Parameters
    ----------
        df : pandas.core.frame.DataFrame
            pandas DataFrame. Should have a nonzero number of rows and
            columns 'url', 'archive_filename'
        outdir : str
            Output directory into which to write the downloaded images.

    Notes
    -----
        df could be, for instance, either a set of raw science exposures or a
        set of calibration images to download.

        Could parallelize for speed-up? At the risk of getting blacklisted...

    """

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for i in range(len(df)):
        print(i+1, ' of ', len(df))
        url = df['url'].iloc[i]

        outname = os.path.join(outdir,
                               os.path.basename(df['archive_filename'].iloc[i]))

        print(url, outname)

        assert(not os.path.exists(outname))

        r = requests.get(url, allow_redirects=True)
        open(outname, 'wb').write(r.content)

def download_raw_science(df):
    """
    Download raw science images.

    Parameters
    ----------
        df : pandas.core.frame.DataFrame
            pandas DataFrame with one row per raw science image file to
            download. Needs to include columns 'archive_filename' and 'url'.

    Notes
    -----
        No outputs, but attempts to download files to disk.
        Basically just a minimal wrapper for download_images function.

    """

    print('DOWNLOADING RAW SCIENCE FRAMES')

    outdir = 'raw' # make this not be hardcoded?
    download_images(df, outdir)

def download_calibs(df):
    """
    Download DECam master calibration images.

    Parameters
    ----------
        df : pandas.core.frame.DataFrame
            pandas DataFrame with one row per master calibration file to
            download. Needs to include columns 'archive_filename' and 'url'.

    Notes
    -----
        No outputs, but attempts to download files to disk.
        Basically just a minimal wrapper for download_images function.

    """

    print('DOWNLOADING NIGHTLY MASTER CALIBRATIONS')

    outdir = 'flats_biases' # make this not be hardcoded
    download_images(df, outdir)

def download_1shard(url, outname):
    """
    Download one reference catalog file using requests.

    Parameters
    ----------
        url : str
            Download URL.
        outname : str
            Output file name to be written.

    Notes
    -----
        No outputs, but attempts to download files to disk.
        Perhaps an alternative way of doing this would be to issue a wget
        command. Not sure about trade-offs with that (would wget be faster?).
        Is it assumed that the output directory corresponding to outname
        already exists?

    """

    r = requests.get(url, allow_redirects=True)
    open(outname, 'wb').write(r.content)

def download_ps1_shards(ras, decs, nmp=None):
    """
    Download PS1 catalogs for a list of DECam pointing (RA, Dec) field centers.

    Parameters
    ----------
        ras : list
            List of DECam field center RA values in decimal degrees. Must
            have same number of elements as decs input.
        decs : list
            List of DECam field center Dec values in decimal degrees. Must
            have same number of elements as ras input.
        nmp : int, optional
            Number of multiprocessing processes to use. Currently not yet
            implemented, and may be a bad idea to implement from the
            perspective of avoiding getting blacklisted.

    Notes
    -----
        Add a check that ras, decs have same number of elements?
        Gets shards list for each ra, dec then does downloads for the unique
        set.

    """

    par = common.decam_params()
    margin = par['shard_cone_margin_deg']
    tables = []
    for i in range(len(ras)):
        shards = getShards_decam_pointing(ras[i], decs[i], depth=7,
                                          margin=margin)
        tables.append(shards)

    table = vstack(tables)

    shards = np.unique(table['shard'])

    base_url = par['base_url_pton']

    outdir = 'ps1_pv3_3pi_20170110'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for i, shard in enumerate(shards):
        print(i+1, ' of ', len(shards))
        _name = str(shard) + '.fits'
        url = os.path.join(base_url, _name)
        print(url)
        outname = os.path.join(outdir, _name)
        download_1shard(url, outname)

def add_exec_permission(fname):
    """
    Add executable permission to a file.

    Parameters
    ----------
        fname : str
            Name of the file for which to add executable permission.

    Notes
    -----
        Should this be made to work with an input consisting of a list of
        file names?

    """

    st = os.stat(fname)
    os.chmod(fname, st.st_mode | stat.S_IXUSR)
