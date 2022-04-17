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

    base_url_short = 'https://astroarchive.noirlab.edu/api/short/ct4m/decam/'
    base_url_pton = \
        'http://tigress-web.princeton.edu/~pprice/ps1_pv3_3pi_20170110/'

    par = {'max_pixel_radius_deg': 1.0923879,
           'shard_cone_margin_deg': 0.1,
           'min_exptime_s': 1.0,
           'defect_basedir' : resource_filename('decam_reduce',
                                  os.path.join('data', 'defects')),
           'ps1_fullsky_dir' : os.getenv('PS1_FULLSKY_DIR'),
           'base_url_short' : base_url_short,
           'base_url_pton' : base_url_pton}

    return par
