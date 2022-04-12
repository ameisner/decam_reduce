"""
decam_reduce.util
=================

A collection of DECam reduction related utility functions.

"""

from astropy.table import Table

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
    radius = 1.0923879

    radius += margin

    result = getShards(ra, dec, radius, depth=7)

    return result
