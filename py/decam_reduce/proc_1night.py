#!/usr/bin/env python

"""
decam_reduce.proc_1night
========================

Pipeline driver for preparing reductions of one observing night of DECam data.

"""

import argparse
import requests
import decam_reduce.util as util
import decam_reduce.common as common
import json
import os
import pandas as pd
from astropy.table import vstack
import numpy as np
import stat
from multiprocessing import Pool

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

    url_short = 'https://astroarchive.noirlab.edu/api/short/ct4m/decam/' + \
                night + '/'
    nightsum = pd.DataFrame(requests.get(url_short).json()[1:])

    return nightsum

def select_raw_science(nightsum, min_exptime_s=None):
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

    # what to do for edge case in which nothing is retained?
    result = nightsum[keep]

    result = result.sort_values('original_filename')

    return result

def select_mastercal(nightsum):
    """
    Select the master calibration records for an observing night.

    Parameters
    ----------
        nightsum : pandas.core.frame.DataFrame
            Observing night summary DataFrame from /short API query.

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
    margin = par['shard_cone_margin_deg'] # deg, not sure...
    tables = []
    for i in range(len(ras)):
        shards = util.getShards_decam_pointing(ras[i], decs[i], depth=7,
                                               margin=margin)
        tables.append(shards)

    table = vstack(tables)

    shards = np.unique(table['shard'])

    # this should be extracted to common.py
    base_url = 'http://tigress-web.princeton.edu/~pprice/ps1_pv3_3pi_20170110/'

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

def write_staging_script(outname, do_ps1_download=False):
    """
    Create a staging script that prepares for reductions with the LSST pipeline.

    Parameters
    ----------
        outname : str
            Name of the shell script to be written.

    Notes
    -----
        Doesn't return anything, but does attempt to write a file.
        Probably needs work in order to propagate/allow various options for the
        underlying commands.

    """

    cmds = []

    cmds.append('mkdir DATA')
    cmds.append('mkdir DATA/CALIB')

    cmds.append('echo lsst.obs.decam.DecamMapper > DATA/_mapper')

    cmds.append('ingestImagesDecam.py DATA --filetype raw raw/*.fz --mode=link')

    cmds.append('ingestCalibs.py DATA --calib DATA/CALIB flats_biases/*.fits.fz --validity 999 --mode=link')

    par = common.decam_params()
    cmds.append('ingestDefects.py DATA ' + par['defect_basedir'] + \
        ' --calib DATA/CALIB')

    ref_cat_dir = 'ps1_pv3_3pi_20170110' if do_ps1_download \
        else par['ps1_fullsky_dir']

    cmds.append('ln -s ' + ref_cat_dir + ' DATA/ref_cats/ps1_pv3_3pi_20170110')

    _cmds = ''
    for cmd in cmds:
        _cmds += cmd + '\n'

    with open(outname, 'wb') as f:
        f.write(_cmds.encode('ascii'))

    add_exec_permission(outname)

def write_launch_script(outname):
    """
    Create a launch script for performing reduction with the LSST pipeline.

    Parameters
    ----------
        outname : str
            Name of the shell script to be written.

    Notes
    -----
        Doesn't return anything, but does attempt to write a file.
        Needs work in order to propagate/allow various options for the
        eventual processCcd.py command.

    """

    cmd = 'processCcd.py DATA --calib DATA/CALIB --rerun processCcdOutputs --id --longlog -j 20'

    with open(outname, 'wb') as f:
        f.write(cmd.encode('ascii'))

    add_exec_permission(outname)

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

def _proc(caldat, limit=None, staging_script_name='stage.sh',
          launch_script_name='launch.sh', do_ps1_download=False):
    """
    Prepare processing for a night of raw DECam data.

    Parameters
    ----------
        caldat : str
            Observing night in format YYYY-MM-DD.
        limit : int, optional
            Only prepare to process first limit raw science exposures (mainly
            for testing purposes).
        staging_script_name : str, optional
            Name of staging shell script to create/write.
        launch_script_name : str, optional
            Name of launch shell script to create/write.
        do_ps1_download : bool, optional
            Whether or not to download PS1 shards from the internet, or
            find them on local disk.

    Notes
    -----
        Move launch.sh, stage.sh defaults into common.py eventually...

    """

    print('WORKING ON NIGHT ' + caldat)

    nightsum = query_night('2018-09-05')
    
    raw = select_raw_science(nightsum)

    if limit is not None:
        raw = raw[0:limit]

    calib = select_mastercal(nightsum)

    download_raw_science(raw)

    download_calibs(calib)

    write_staging_script(staging_script_name, do_ps1_download=do_ps1_download)
    write_launch_script(launch_script_name)

    if do_ps1_download:
        download_ps1_shards(np.array(raw['ra_min']),
                            np.array(raw['dec_min']))
    else:
        print('ATTEMPTING TO USE PS1 REFERENCE CATALOGS ALREADY ON DISK')
    
if __name__ == "__main__":
    descr = 'prepare processing for a night of raw DECam data'

    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('caldat', type=str, nargs=1,
                        help="observing night in YYYY-MM-DD format")

    parser.add_argument('--repo_name', default='DATA', type=str,
                        help="Butler repository name")

    parser.add_argument('--staging_script_name', default='stage.sh', type=str,
                        help="output name for repo staging script")

    parser.add_argument('--launch_script_name', default='launch.sh', type=str,
                        help="output name for processing launch script")

    parser.add_argument('--multiproc', default=None, type=int,
                        help="number of threads for multiprocessing")

    parser.add_argument('--limit', default=None, type=int,
                        help="process only first limit exposures")
    
    # I guess filter here should be an abbreviated name like 'u', 'g', 'r', ...
    # rather than the full filter name?
    parser.add_argument('--filter', default=None, type=str,
                        help="only process raw science data with this filter")

    # PROPID format?
    parser.add_argument('--propid', default=None, type=str,
                        help="only process raw science data with this propid")

    parser.add_argument('--do_ps1_download', default=False, action='store_true',
                        help="download PS1 shard files from the internet?")

    args = parser.parse_args()

    _proc(args.caldat[0], limit=args.limit,
          staging_script_name=args.staging_script_name,
          launch_script_name=args.launch_script_name,
          do_ps1_download=args.do_ps1_download)
