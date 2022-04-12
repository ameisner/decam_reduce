from lsst.meas.algorithms.htmIndexer import HtmIndexer
from lsst.geom import degrees
from lsst.geom import SpherePoint
from astropy.table import Table

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
