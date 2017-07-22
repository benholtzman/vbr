#!/bin/bash

cwdir=$(pwd)
cd 
hmdir=$(pwd)
cd $cwdir

writename=Job_started_${1}_${2}
echo "Matlab job started! log: $2 PID: $1" > "$hmdir/Dropbox/logs/$writename"
while [ -e /proc/$1 ]
do
    sleep 60s
done
writename=Job_completed_${1}_${2}
echo "Matlab job complete! log: $2 PID: $1" > "$hmdir/Dropbox/logs/$writename"
