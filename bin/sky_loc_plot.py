import decam_reduce.util as util
import decam_reduce.plotting as plotting

nightsum = util.query_night('2018-09-05')

# restrict to actual science frames

science = util.select_raw_science(nightsum)

plotting.exp_sky_locations(science, coordsys='equ')
