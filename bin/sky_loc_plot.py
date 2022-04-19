#!/usr/bin/env python

import decam_reduce.util as util
import decam_reduce.plotting as plotting
import sys

nightsum = util.query_night('2018-09-05')

# restrict to actual science frames

science = util.select_raw_science(nightsum)

if len(sys.argv) == 1:
    coordsys = 'equ'
else:
    coordsys = sys.argv[1]

assert(coordsys in ['equ', 'gal'])

plotting.exp_sky_locations(science, coordsys=coordsys)
