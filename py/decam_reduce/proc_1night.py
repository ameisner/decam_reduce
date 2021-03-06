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
import multiprocessing

def write_staging_script(outname, do_ps1_download=False, repo_name='DATA'):
    """
    Create a staging script that prepares for reductions with the LSST pipeline.

    Parameters
    ----------
        outname : str
            Name of the shell script to be written.
        repo_name : str, optional
            Butler repository name. Default is DATA. Perhaps this special 
            default value should be extracted to common.py.

    Notes
    -----
        Doesn't return anything, but does attempt to write a file.
        Probably needs work in order to propagate/allow various options for the
        underlying commands.

    """

    cmds = []

    cmds.append('mkdir -p ' + repo_name)
    cmds.append('mkdir -p ' + repo_name + '/CALIB')

    cmds.append('echo lsst.obs.decam.DecamMapper > ' + repo_name + '/_mapper')

    cmds.append('ingestImagesDecam.py ' + repo_name + ' --filetype raw raw/*.fz --mode=link')

    cmds.append('ingestCalibs.py ' + repo_name + ' --calib ' + repo_name + '/CALIB flats_biases/*.fits.fz --validity 999 --mode=link')

    par = common.decam_params()
    cmds.append('ingestDefects.py ' + repo_name + ' ' + \
        par['defect_basedir'] + ' --calib ' + repo_name + '/CALIB')

    ref_cat_dir = 'ps1_pv3_3pi_20170110' if do_ps1_download \
        else par['ps1_fullsky_dir']

    # this would happen if PS1_FULLSKY_DIR environment variable is not set
    assert(ref_cat_dir is not None)

    cmds.append('mkdir -p ' + repo_name + '/ref_cats')
    cmds.append('ln -s ' + ref_cat_dir + ' ' + repo_name + \
        '/ref_cats/ps1_pv3_3pi_20170110')

    _cmds = ''
    for cmd in cmds:
        _cmds += cmd + '\n'

    with open(outname, 'wb') as f:
        f.write(_cmds.encode('ascii'))

    util.add_exec_permission(outname)

def write_launch_script(outname, caldat, nmp=None, repo_name='DATA',
                        skip_fringe=False, maxRefObjects=None):
    """
    Create a launch script for performing reduction with the LSST pipeline.

    Parameters
    ----------
        outname : str
            Name of the shell script to be written.
        caldat : str
            Observing night in format YYYY-MM-DD.
        nmp : int, optional
            Number of CPUs for processCcd.py to use for reductions. If not 
            specified, then the number of CPUs will be chosen as half of the
            total number of CPUs available on the machine.
        repo_name : str, optional
            Butler repository name. Default is DATA. Perhaps this special 
            default value should be extracted to common.py.
        skip_fringe : bool, optional
            If set True, skip fringe correction step during detrending. Default
            value is False.
        maxRefObjects : int, optional
            Maximum number of astrometric reference objects to retain. Gets
            propagated into calibrate.astrometry.matcher.maxRefObjects config
            parameter in the launch script. LSST pipeline default is 65536,
            and maxRefObjects should not be set to a value larger than that.

    Notes
    -----
        Doesn't return anything, but does attempt to write a file.
        Needs work in order to propagate/allow various options for the
        eventual processCcd.py command.

    """

    if nmp is None:
        nmp = multiprocessing.cpu_count() // 2

    # -j 1 is fine for processCcd.py even though I believe this is
    # the same as just running the CCDs in serial...
    assert((nmp >= 1) and (nmp <= multiprocessing.cpu_count()))

    do_fringe_str = ' --config isr.doFringe=False' if skip_fringe else ''

    maxRefObjects_str = ''
    if maxRefObjects is not None:
        # extract this 65536 special number eventually
        if (maxRefObjects > 65536) or (maxRefObjects < 1):
            print('INVALID maxRefObjects INPUT')
        else:
            maxRefObjects_str = ' --config calibrate.astrometry.matcher.maxRefObjects=' + str(maxRefObjects)

    cmd = 'processCcd.py ' + repo_name + ' --calib ' + repo_name + \
          '/CALIB --rerun processCcdOutputs --id --longlog -j ' + \
          str(nmp) + do_fringe_str + maxRefObjects_str + ' &> processCcd_' + \
          caldat + '.log &\n'

    with open(outname, 'wb') as f:
        f.write(cmd.encode('ascii'))

    util.add_exec_permission(outname)

def _proc(caldat, limit=None, staging_script_name='stage.sh', repo_name='DATA',
          launch_script_name='launch.sh', do_ps1_download=False, nmp=None,
          _filter=None, propid=None, bgal_min=None, expnum=None,
          skip_raw_download=False, skip_fringe=False,
          maxRefObjects=None):
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
        nmp : int, optional 
            Number of CPUs for processCcd.py to use for reductions. If not 
            specified, then the number of CPUs will be chosen as half of the
            total number of CPUs available on the machine.
        repo_name : str, optional
            Butler repository name. Default is DATA. Perhaps this special 
            default value should be extracted to common.py.
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
        skip_raw_download : bool, optional
            If set True, skip actually downloading the raw DECam image files.
            Meant to be used for debugging, to speed things up when perhaps
            only the downloaded master calibration files are of interest, but
            not the raw science images themselves.
        skip_fringe : bool, optional
            If set True, skip fringe correction step during detrending. Default
            value is False.
        maxRefObjects : int, optional
            Maximum number of astrometric reference objects to retain. Gets
            propagated into calibrate.astrometry.matcher.maxRefObjects config
            parameter in the launch script. LSST pipeline default is 65536,
            and maxRefObjects should not be set to a value larger than that.

    Notes
    -----
        Move launch.sh, stage.sh defaults into common.py eventually...

    """

    util.print_hostname()
    util.print_process_id()

    print('WORKING ON NIGHT ' + caldat)

    nightsum = util.query_night(caldat)

    raw = util.select_raw_science(nightsum, _filter=_filter, propid=propid,
                                  bgal_min=bgal_min, expnum=expnum)

    if limit is not None:
        raw = raw[0:limit]

    calib = util.select_mastercal(nightsum, raw=raw)

    if not skip_raw_download:
        util.download_raw_science(raw)

    util.download_calibs(calib)

    write_staging_script(staging_script_name, do_ps1_download=do_ps1_download,
                         repo_name=repo_name)
    write_launch_script(launch_script_name, caldat, nmp=nmp, 
                        repo_name=repo_name, skip_fringe=skip_fringe,
                        maxRefObjects=maxRefObjects)

    if do_ps1_download:
        util.download_ps1_shards(np.array(raw['ra_min']),
                                 np.array(raw['dec_min']))
    else:
        print('ATTEMPTING TO USE PS1 REFERENCE CATALOGS ALREADY ON DISK')

    print('DONE')
    
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
    # ability to have a list of filters e.g., g and r only ?
    parser.add_argument('--filter', default=None, type=str,
                        help="only process raw science data with this filter")

    parser.add_argument('--propid', default=None, type=str,
                        help="only process raw science data with this propid")

    parser.add_argument('--do_ps1_download', default=False, action='store_true',
                        help="download PS1 shard files from the internet?")

    parser.add_argument('--bgal_min', default=None, type=float,
                        help="minimum absolute Galactic latitude in degrees")

    # DECam apparently uses EXPNUM rather than EXPID
    parser.add_argument('--expnum', default=None, type=int,
                        help="restrict to a specific exposure number")

    parser.add_argument('--skip_raw_download', default=False,
                        action='store_true',
                        help="skip downloading raw files; meant for debugging")

    parser.add_argument('--skip_fringe', default=False, action='store_true',
                        help="skip fringe correction")

    parser.add_argument('--maxRefObjects', default=None, type=int,
                        help="maximum number of astrometric reference objects")

    args = parser.parse_args()

    _proc(args.caldat[0], limit=args.limit, repo_name=args.repo_name,
          staging_script_name=args.staging_script_name,
          launch_script_name=args.launch_script_name,
          do_ps1_download=args.do_ps1_download, nmp=args.multiproc,
          _filter=args.filter, propid=args.propid, bgal_min=args.bgal_min,
          expnum=args.expnum, skip_raw_download=args.skip_raw_download,
          skip_fringe=args.skip_fringe, maxRefObjects=args.maxRefObjects)
