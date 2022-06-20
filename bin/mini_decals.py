#!/usr/bin/env python

"""
The idea of the 'mini-DECaLS' data set is to process 12 entire (somewhat)
random DECaLS (PROPID = 2014B-0404) nights of science exposures.

"""

from decam_reduce.proc_1night import _proc
from pkg_resources import resource_filename
import astropy.io.fits as fits
import os

def stage_mini_decals():
    """
    Set up for mini-DECaLS data set LSST pipeline processing.

    Notes
    -----
        Should also have this create something like a 'stage_all.sh'
        meta-staging script.

    """

    mini_decals_nights = ['2014-08-14', '2014-12-30', '2015-03-27', '2015-04-07',
                          '2016-12-20', '2016-06-07', '2017-03-15', '2017-05-21',
                          '2018-02-17', '2018-10-30', '2019-02-06', '2019-01-01']

    for caldat in mini_decals_nights:
        staging_script_name = 'stage-' + caldat + '.sh'
        launch_script_name = 'launch-' + caldat + '.sh'

        _proc(caldat, limit=None, staging_script_name=staging_script_name, 
              repo_name='DATA', launch_script_name=launch_script_name, 
              do_ps1_download=False, nmp=None, _filter=None, propid=None, 
              bgal_min=None, skip_fringe=True, maxRefObjects=3000)

if __name__=="__main__":
    stage_mini_decals()
