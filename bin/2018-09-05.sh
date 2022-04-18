PROC_1NIGHT_PATH=$(echo -e import decam_reduce.proc_1night "\n"print\(decam_reduce.proc_1night.__file__\) |python)

python -u $PROC_1NIGHT_PATH 2018-09-05 --limit 1 &> 2018-09-05.log &
