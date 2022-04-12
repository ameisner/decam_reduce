"""
decam_reduce.common
===================

Central location for parameters expected to be used in various other places
throughout the codebase. Goal is to factor out special numbers.

"""

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
           'shard_cone_margin_deg': 0.1}

    return par
