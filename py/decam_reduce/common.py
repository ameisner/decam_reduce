"""
decam_reduce.common
===================

Central location for parameters expected to be used in various other places
throughout the codebase. Goal is to factor out special numbers.

"""

from pkg_resources import resource_filename
import os

def decam_params():
    """
    Set of configuration parameters relevant to DECam reductions.

    Returns
    -------
        par : dict
            Dictionary of parameters.

    Notes
    -----
        Eventually generalize somehow for other cameras?

    """

    par = {'max_pixel_radius_deg': 1.0923879,
           'shard_cone_margin_deg': 0.1,
           'min_exptime_s': 1.0,
           'defect_basedir' : resource_filename('decam_reduce',
                                  os.path.join('data', 'defects')),
           'ps1_fullsky_dir' : os.environ['PS1_FULLSKY_DIR']}

    return par
