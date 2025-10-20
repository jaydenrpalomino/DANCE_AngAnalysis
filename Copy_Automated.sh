#!/bin/bash
#set -x

rsync -av /home/jr2514/DANCE/DANCE_Analysis/TimeDeviations_Run_10{4,5,6,8}*.txt ../TimeDeviations/

rsync -av /home/jr2514/DANCE/DANCE_Alpha_Calibrator/ParamOutput/ ~/DANCE/Calibrations/
