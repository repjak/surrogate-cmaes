#!/bin/bash
#
# where
#   EXPID       -- string with experiment's unique ID
#   META_QUEUE  -- string with walltime for Metacentrum (2h/4h/1d/2d/1w)
#   ID1,ID2,... -- integers defining numeric IDs of concrete experiments to run
#               all IDs are taken from file 'allids.txt' if no IDs supplied on
#               command-line (allids.txt is expected in experiment's directory
#               and should be produced by the experiment generator script)
#
# settings within this file:
#   EXPPATH_SHORT  -- $CWD/experiments


# QUEUE = Metacentrum walltime (2h/4h/1d/2d/1w) -- queue will be decided accordingly
QUEUE=1d
export EXPID='modelTest_04'

# CWD = Directory of this particular file
CWD=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
. $CWD/../bash_settings.sh
export EXPPATH_SHORT="$CWD"
export ID
export DIM
export FUNC
export MATLAB_FCN


subtask() {
  export EXPPATH_SHORT
  export ID
  export DIM
  export FUNC=`echo $FUNC | tr '\n' ' '`
  echo ID=$ID : DIM=$DIM : MATLAB_FCN=$MATLAB_FCN : FUNC=$FUNC
  qsub -N "${EXPID}__$1" -l "walltime=$QUEUE" -v FUNC,DIM,MATLAB_FCN,EXPPATH_SHORT $EXPPATH_SHORT/modelTesting_metajob.sh && echo "submitted ok."
  ID=$((ID+1))
}


ID=1;

DIM=2;MATLAB_FCN="modelTesting_ordgp_03";FUNC=`seq 1 12`;
subtask $ID
DIM=2;MATLAB_FCN="modelTesting_ordgp_03";FUNC=`seq 13 24`;
subtask $ID
DIM=5;MATLAB_FCN="modelTesting_ordgp_03";FUNC=`seq 1 8`;
subtask $ID
DIM=5;MATLAB_FCN="modelTesting_ordgp_03";FUNC=`seq 9 16`;
subtask $ID
DIM=5;MATLAB_FCN="modelTesting_ordgp_03";FUNC=`seq 17 24`;
subtask $ID
DIM=10;MATLAB_FCN="modelTesting_ordgp_03";FUNC=`seq 1 6`;
subtask $ID
DIM=10;MATLAB_FCN="modelTesting_ordgp_03";FUNC=`seq 7 12`;
subtask $ID
DIM=10;MATLAB_FCN="modelTesting_ordgp_03";FUNC=`seq 13 18`;
subtask $ID
DIM=10;MATLAB_FCN="modelTesting_ordgp_03";FUNC=`seq 19 24`;
subtask $ID

