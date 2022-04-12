# reduce DECam data using the LSST science pipeline

# environment configuration

Add the `decam_reduce/py` directory to your PYTHONPATH. Then you should be able to do stuff like the following:

```
In [1]: import decam_reduce

In [2]: ?decam_reduce
Type:        module
String form: <module 'decam_reduce' from '/home/ameisner/decam_reduce/py/decam_reduce/__init__.py'>
File:        ~/decam_reduce/py/decam_reduce/__init__.py
Docstring:  
decam_reduce
============
A setup for reducing raw DECam images with the LSST pipeline.
In [3]: import decam_reduce.util as util
```
# full pipeline help

    python proc_1night.py --help
    usage: proc_1night.py [-h] [--repo_name REPO_NAME]
                          [--staging_script_name STAGING_SCRIPT_NAME]
                          [--launch_script_name LAUNCH_SCRIPT_NAME]
                          [--multiproc MULTIPROC] [--limit LIMIT]
                          [--filter FILTER] [--propid PROPID]
                          caldat

    process a night of raw DECam data

    positional arguments:
      caldat                observing night in YYYY-MM-DD format

    optional arguments:
      -h, --help            show this help message and exit
      --repo_name REPO_NAME
                            Butler repository name
      --staging_script_name STAGING_SCRIPT_NAME
                            output name for repo staging script
      --launch_script_name LAUNCH_SCRIPT_NAME
                            output name for processing launch script
      --multiproc MULTIPROC
                            number of threads for multiprocessing
      --limit LIMIT         process only first limit exposures
      --filter FILTER       only process raw science data with this filter
      --propid PROPID       only process raw science data with this propid
