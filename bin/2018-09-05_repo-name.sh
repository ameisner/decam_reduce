PROC_1NIGHT_PATH=$(./proc_1night_path.sh)

python -u $PROC_1NIGHT_PATH 2018-09-05 --limit 1 --repo_name _DATA_ &> 2018-09-05_repo-name.log &
