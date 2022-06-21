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
from astropy.coordinates import SkyCoord
import astropy.units as u
import astropy.io.fits as fits
from astropy.coordinates import Angle
import copy
from functools import lru_cache

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

    # augment with lgal, bgal columns
    # this may not be the best / proper place to do this -- may need 
    # clean-up eventually
    lgal, bgal = galactic_coords(nightsum['ra_min'], nightsum['dec_min'])

    nightsum['lgal'] = lgal
    nightsum['bgal'] = bgal

    return nightsum

def select_raw_science(nightsum, min_exptime_s=None, _filter=None, 
                       propid=None, bgal_min=None, expnum=None):
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
        bgal_min : float, optional
            Minimum (field center) absolute Galactic latitude in degrees.
            The idea is to avoid running jobs near the Galactic plane/center
            that use excessive amounts of memory. Default value of None means
            no |b_gal| but will get enforced.
        expnum : int, optional
            Specific single EXPNUM to restrict to. Defaul is None, in which
            case no cut on EXPNUM gets made. The thought is that this optional
            argument could be useful for testing/debugging purposes.

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

    if bgal_min is not None:
        # any checking of bgal_min (value only makes sense if >= 0 and < 90)
        keep = np.logical_and(keep, np.abs(nightsum['bgal']) > bgal_min)

    if expnum is not None:
        # clean this up...
        expnums = np.array([int(os.path.basename(origname)[6:14]) for origname in nightsum['original_filename']])
        keep = np.logical_and(keep, expnums == expnum)
        assert(np.sum(keep) == 1)

    # what to do for edge case in which nothing is retained?
    result = nightsum[keep]

    result = result.sort_values('original_filename')

    return result

def downselect_zeros(tab):
    """
    Choose best master calibration zero among many options for a single night.

    Parameters
    ----------
        tab : pandas.core.frame.DataFrame
            pandas DataFrame with one row per mastercal zero product.

    Returns
    -------
        tab : pandas.core.frame.DataFrame
            Possibly modified version of input pandas DataFrame.

    Notes
    -----
        Assume input table 'tab' includes only valid mastercal images.
        Can this function be deleted?

    """

    n_zeros = np.sum(tab['exposure'] == 0)

    # don't do anything in this case
    if n_zeros < 2:
        return tab
    else:
        keep = np.ones(len(tab), dtype=bool)
        for i in range(len(tab)):
            if tab['exposure'].iloc[i] != 0:
                continue
            if tab['ifilter'].iloc[i].strip() != 'solid plate 0.0 0.0':
                keep[i] = False
            if (n_zeros-np.sum(np.logical_not(keep))) == 1:
                break

    tab = tab[keep]

    print('TRIMMED OUT ', np.sum(np.logical_not(keep)), ' EXTRA ZEROS')
    return tab

def _remove_morning_master(tab):
    """
    Remove morning master calibrations, which appear to cause problems.

    Parameters
    ----------
        tab : pandas.core.frame.DataFrame
            pandas DataFrame with one row per mastercal product.

    Returns
    -------
        tab : pandas.core.frame.DataFrame
            Possibly modified version of input pandas DataFrame.

    Notes
    -----
        This is basically a hack; should eventually figure out a better
        solution for this situation.

    """

    _char = np.zeros(len(tab), dtype='U1')
    for i in range(len(tab)):
        # dateobs_min, dateobs_max are the same usually?
        date = tab['dateobs_min'].iloc[i]
        pos = date.find('T') + 1
        _char[i] = date[pos]

    _char = np.array(_char)

    keep = (_char != '0')

    return tab[keep]

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
        if (result['obs_type'].iloc[i] == 'dome flat') or \
           (result['exposure'].iloc[i] != 0):
            if result['ifilter'].iloc[i] not in filters:
                keep[i] = False
        # for now, throw out dome flats marked as zeros
        # but eventually will want to try to salvage these
        # and use them as the dome flats that they are
        if (result['obs_type'].iloc[i] == 'zero') and \
           (result['exposure'].iloc[i] != 0):
                keep[i] = False

    result = result[keep]

    # now remove duplicates by requiring unique EXPNUM
    # apparently we don't directly have EXPNUM from the /short API
    # but we can (hopefully) use original_filename in an equivalent way
    _, unique_indices = np.unique(result['original_filename'],
                                  return_index=True)

    result = result.iloc[unique_indices]

    result = _remove_morning_master(result)

    result = downselect_zeros(result)

    return result

def download_images(df, outdir, return_outnames=False):
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

    outnames = []
    for i in range(len(df)):
        print(i+1, ' of ', len(df))
        url = df['url'].iloc[i]

        outname = os.path.join(outdir,
                               os.path.basename(df['archive_filename'].iloc[i]))

        print(url, outname)

        assert(not os.path.exists(outname))

        r = requests.get(url, allow_redirects=True)
        open(outname, 'wb').write(r.content)
        outnames.append(outname)

    if return_outnames:
        return outnames

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
    fnames = download_images(df, outdir, return_outnames=True)

    for fname in fnames:
        dummy_mastercal_hdu(fname)

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

def trixel_number_to_depth(trixel):
    """
    Infer the HTM depth from the trixel number.

    Parameters
    ----------
        trixel : int
            HTM trixel number.

    Returns
    -------
        depth : int
            HTM depth parameter.
        is_valid : bool
            Boolean indicating whether the integer supplied is in fact a
            valid HTM trixel number.

    Notes
    -----
        Currently only works for scalar input 'trixel' value.

    """

    def log4(x):
        """
        Log base 4.

        Parameters
        ----------
            x : float
                Value for which to compute log base 4.

        Returns
        -------
            y : float
                Log base 4 of x.

        """

        return np.log(x) / np.log(4)

    depth = int(np.floor(log4(trixel / 2.0) - 1))

    is_valid = trixel < (4**(depth+1)*4)

    return depth, is_valid

def galactic_coords(ras, decs):
    """
    Compute galactic coordinates for a list of equatorial coordinates.

    Parameters
    ----------
        ras : numpy.ndarray
            List of RA coordinates in decimal degrees. Must have same size as
            decs input
        decs : numpy.ndarray
            List of Dec coordinates in decimal degrees. Must have same size as
            ras input.

    Returns
    -------
        lgal : numpy.ndarray
            Galactic longitude values in decimal degrees.
        bgal : numpy.ndarray
            Galactic latitude values in decimal degrees.

    """

    skycoords = SkyCoord(ras*u.deg, decs*u.deg,
                         frame='icrs')

    lgal = np.array(skycoords.galactic.l, dtype=float)
    bgal = np.array(skycoords.galactic.b, dtype=float)

    return lgal, bgal

def get_raw_ccds_list(fname):
    """
    Get list of CCD names/numbers present in a raw DECam image file.

    Parameters
    ----------
        fname : str
            Full file name of a raw DECam exposure.

    Returns
    -------
        ccds : list
            List of string CCD names (e.g., 'N10')
        ccdnums : list
            Corresponding list of integer CCD numbers (e.g., 41)

    """

    hdul = fits.open(fname)

    ccds = []
    ccdnums = []
    for i in range(len(hdul)):
        h = hdul[i].header
        if 'EXTNAME' in h:
            if h['EXTNAME'][0] in ['N', 'S']:
                ccds.append(h['EXTNAME'])
                ccdnums.append(h['CCDNUM'])

    return ccds, ccdnums

def get_expid(fname_raw):
    """
    Retrieve a raw DECam file's exposure number from its header.

    Parameters
    ----------
        fname_raw : str
            Full file name of a raw DECam exposure.

    Returns
    -------
         : int
        Exposure number from the EXPNUM header card.

    Notes
    -----
        Is it safe to just go based on the raw file name, assuming that
        always has the form DECam_????????.fits.fz ?

    """

    hdul = fits.open(fname_raw)

    return hdul[0].header['EXPNUM']

def header_radec_to_decimal(h):
    """
    Convert DECam raw header format (RA, Dec) to decimal degrees.

    Parameters
    ----------
        h : astropy.io.fits.header.Header
            Primary extension header of raw DECam image file. Expected to
            have 'RA', 'DEC' cards with values that are sexagesimal/string.

    Returns
    -------
        ra_deg : float
            RA of field center in decimal degrees.
        dec_deg : float
            Dec of field center in decimal degrees.

    """

    ra = Angle(h['RA'] + ' hours')
    dec = Angle(h['DEC'] + ' degrees')

    ra_deg = ra.deg
    dec_deg = dec.deg

    return ra_deg, dec_deg

def dummy_mastercal_hdu(fname):
    """
    Add dummy extensions to keep CCDNUM and extension number in alignment.

    Parameters
    ----------
        fname : str
            Full file name of master calibration file.

    Notes
    -----
        File gets overwritten if nonzero number of dummy extensions are needed.
        Needs major rewrite/cleanup.

    """

    hdul = fits.open(fname)

    if len(hdul) == 61:
        # assume that this means that CCDNUM=2 is missing AND
        # CCDNUM = 61 is also missing
        print('REPAIRING MISSING MASTERCAL EXTENSION : ', fname)
        _hdul = hdul[0:2]
        _hdul.append(_hdul[1]) #dummy, choice of 1 here shouldn't matter
        
        #_hdul.append(hdul[2:])
        for i in range(2, len(hdul)):
            _hdul.append(hdul[i])
        del hdul
        hdul = copy.deepcopy(_hdul)

    if len(hdul) == 62:
        # assume this means that CCDNUM=61 is missing
        print('REPAIRING MISSING MASTERCAL EXTENSION : ', fname)

        _hdul = hdul[0:61]

        _hdul.append(_hdul[60])
        _hdul.append(hdul[61])

        outname_tmp = fname + '.tmp'
        print(len(_hdul))
        _hdul.writeto(outname_tmp)
        os.rename(outname_tmp, fname) # overwrite

@lru_cache(maxsize=1)
def get_decam_mapper():
    """
    Retrieve DECamMapper object, with caching.

    Notes
    -----
        mapper = DecamMapper() is apparently slow (0.5-1 seconds),
        which I assume is due to some I/O that happens behind the scenes,
        hence the caching.
        Does DECamMapper include non-science (guide, focus) detector mapping?

        N30 is entirely absent??? Probably will need a workaround for this...
        N30 <-> CCDNUM = 61

    """

    from lsst.obs.decam import DecamMapper

    mapper = DecamMapper()

    return mapper


def ccdname_to_ccdnum(ccdname):
    """
    Convert from DECam CCDNAME to CCDNUM.

    Parameters
    ----------
        ccdname : str
            CCDNAME like 'S30'. Needs to be a valid DECam CCDNAME,
            excluding N30.

    Returns
    -------
        ccdnum : int
            CCDNUM, should be an integer in the range [1, 62], excluding 61.

    Notes
    -----
        Not currently vectorized.

    """

    if ccdname == 'N30':
        return 61

    mapper = get_decam_mapper()

    _mapping = mapper.detectorNames

    # swap keys/values
    mapping = dict((v,k) for k,v in _mapping.items())

    if ccdname not in mapping.keys():
        print('INVALID CCDNAME VALUE : ', ccdname)
        assert(False)

    ccdnum = mapping[ccdname]

    return ccdnum

def ccdnum_to_ccdname(ccdnum):
    """
    Convert from DECam CCDNUM to CCDNAME.

    Parameters
    ----------
        ccdnum : int
            CCDNUM, should be an integer in the range [1, 62], excluding 61.

    Returns
    -------
        ccdname : str
            CCDNAME like 'S30'.

    Notes
    -----
        Not currently vectorized.

    """

    if ccdnum == 61:
        return 'N30'

    mapper = get_decam_mapper()

    mapping = mapper.detectorNames

    if ccdnum not in mapping.keys():
        print('INVALID CCDNUM VALUE : ', ccdnum)
        assert(False)

    ccdname = mapping[ccdnum]

    return ccdname
