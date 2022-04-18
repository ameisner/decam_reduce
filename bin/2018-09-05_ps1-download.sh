PROC_1NIGHT_PATH=$(./proc_1night_path.sh)

python -u $PROC_1NIGHT_PATH 2018-09-05 --limit 1 --do_ps1_download &> 2018-09-05_ps1-download.log &
