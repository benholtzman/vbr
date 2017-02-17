#!/bin/bash
# Usage: ./RunMatBG 'outputfile.txt' 
# outputfile.txt will capture all the matlab screen output, can rename
# to whatever you like.

/usr/local/MATLAB/R2014a/bin/matlab -nodisplay -nodesktop -r "cd ..; DRIVE_TwoPhase_a" > $1 &
#/usr/local/MATLAB/R2014a/bin/matlab -nodisplay -nodesktop -r "DRIVE_VBR" > $1 &
Job_PID=$! # process id
notify-send "Started matlab job, screen output to $1, PID $Job_PID"
stty echo
./sh_check_job.sh $Job_PID $1 &

