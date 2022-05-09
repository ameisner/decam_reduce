import decam_reduce.plotting as plotting

fname_raw = '/data0/ameisner/test_decam_reduce-maxref/raw/c4d_180905_231618_ori.fits.fz'

rerun_dir = '/data0/ameisner/test_decam_reduce-maxref/DATA/rerun/processCcdOutputs_j6'
rerun_dir_maxref = '/data0/ameisner/test_decam_reduce-maxref/DATA/rerun/processCcdOutputs-calibrate_only'

plotting.outputs_fp_map(fname_raw, rerun_dir, save=True, outname_extra='',
                        title_extra='maxRefObjects=65536')
plotting.outputs_fp_map(fname_raw, rerun_dir_maxref, save=True,
                        outname_extra='-maxref_3000', title_extra='maxRefObjects=3000')
