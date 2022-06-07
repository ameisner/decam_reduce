#!/usr/bin/env python

"""
The idea of the 'mono-DECaLS' data set is to process one random DECaLS
(PROPID = 2014B-0404) science exposure per observing night for during which
any DECaLS data were taken.

"""

from decam_reduce.proc_1night import _proc
from pkg_resources import resource_filename
import astropy.io.fits as fits
import os

def stage_mono_decals():
    """
    Set up for mono-DECaLS data set LSST pipeline processing.

    Notes
    -----
        Should also have this create something like a 'stage_all.sh'
        meta-staging script.

    """

    fname_in = resource_filename('decam_reduce',
                                  os.path.join('data',
                                 'decals_unique_nights.fits'))

    tab = fits.getdata(fname_in)

    for row in tab:
        caldat = row['NIGHT']
        expnum = row['EXPNUM_RANDOM']
        expnum_string = str(expnum).zfill(8)
        staging_script_name = 'stage-' + expnum_string + '.sh'
        launch_script_name = 'launch-' + expnum_string + '.sh'

        _proc(caldat, limit=None, staging_script_name=staging_script_name, 
              repo_name='DATA', launch_script_name=launch_script_name, 
              do_ps1_download=False, nmp=None, _filter=None, propid=None, 
              bgal_min=None, expnum=expnum, skip_fringe=True,
              maxRefObjects=3000)

if __name__=="__main__":
    stage_mono_decals()
