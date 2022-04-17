#!/usr/bin/env python

"""
decam_reduce.proc_1night
========================

Pipeline driver for preparing reductions of one observing night of DECam data.

"""

import argparse
import decam_reduce.util as util
import decam_reduce.common as common
import numpy as np

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

    # this would happen if PS1_FULLSKY_DIR environment variable is not set
    assert(ref_cat_dir is not None)

    cmds.append('ln -s ' + ref_cat_dir + ' DATA/ref_cats/ps1_pv3_3pi_20170110')

    _cmds = ''
    for cmd in cmds:
        _cmds += cmd + '\n'

    with open(outname, 'wb') as f:
        f.write(_cmds.encode('ascii'))

    util.add_exec_permission(outname)

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

    util.add_exec_permission(outname)

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

    nightsum = util.query_night('2018-09-05')
    
    raw = util.select_raw_science(nightsum)

    if limit is not None:
        raw = raw[0:limit]

    calib = util.select_mastercal(nightsum)

    util.download_raw_science(raw)

    util.download_calibs(calib)

    write_staging_script(staging_script_name, do_ps1_download=do_ps1_download)
    write_launch_script(launch_script_name)

    if do_ps1_download:
        util.download_ps1_shards(np.array(raw['ra_min']),
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
