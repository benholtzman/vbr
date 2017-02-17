#!/bin/bash
# Usage: ./RunMatBG 'outputfile.txt' 
# outputfile.txt will capture all the matlab screen output, can rename
# to whatever you like.

/usr/local/MATLAB/R2014a/bin/matlab -nodisplay -nodesktop -r "cd ..; DRIVE_TwoPhase_par" > $1 &
#/usr/local/MATLAB/R2014a/bin/matlab -nodisplay -nodesktop -r "DRIVE_VBR" > $1 &
Job_PID=$! # process id
stty echo
./sh_check_jobquiet.sh $Job_PID $1 &

