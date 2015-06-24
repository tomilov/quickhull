#!/usr/bin/env bash -vex

#rbox D3 10 s t428
#rbox D3 5 s t30165
#rbox D3 5 s t40466

pushd bin
pids=()
let delta=10000
for i in $(seq 1 $NUMBER_OF_PROCESSORS)
do
    let 'lower = delta * i'
    let 'upper = lower + delta'
    (
        for j in $(seq $lower $upper)
        do
            rbox D3 10 s t$j | /tmp/quickhull
        done
    ) &> test_$i.txt &
    pids+=("$!")
done
echo $pids
popd
