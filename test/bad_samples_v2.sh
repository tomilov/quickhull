#/usr/bin/env bash

PROJECT_ROOT=$(dirname $0)/..

function on_exit () {
	for job in `jobs -p`;
	do
		echo "kill job = $job at exit"
		kill -KILL $job
	done
	exit "$?"
}
trap on_exit ERR EXIT SIGINT SIGTERM

echo > log

let delta=10000
for i in $(seq 1 `nproc`)
do
    let 'lower = delta * i'
    let 'upper = lower + delta'
    (
	echo "started job $$"
        for j in $(seq $lower $upper)
        do
		rbox n D3 10000 t$j > /tmp/qh$j
		LHS=`cat /tmp/qh$j | $PROJECT_ROOT/bin/quickhull | head -7 | grep "#number of facets:" | cut -d' ' -f 4`
		RHS=`qconvex Qt Tv s TI /tmp/qh$j 2>&1 | grep "Number of facets:" | cut -d' ' -f6`
		rm /tmp/qh$j
		if [ "$LHS" != "$RHS" ];
		then
			echo "Error: j = $j LHS=$LHS RHS=$RHS"
		fi
        done
	echo "finished job $$"
    ) &> test_$i.txt &
    pids+=("$!")
    echo "started job # $i pid = $!" | tee -a log
done
echo "pids = $pids" | tee log
echo "waiting for jobs = "`jobs -p` | tee -a log
for job in `jobs -p`
do
	echo "waiting for job # $job"
	wait $job || echo "fail waiting for job # $job"
done
