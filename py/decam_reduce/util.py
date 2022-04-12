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
