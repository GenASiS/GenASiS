#!/bin/bash

#-- Script used by Jenkins to submit test job, wait on it, and extract runtime
#   data

#-- Order of Arguments are: JenkinsBuildID nProcs nThreads Exec ExecArgs
JENKINS_BID=$1
NPROCS=$2
NTHREADS=$3
EXEC=$4
shift 4
EXEC_ARGS=$@

case $GENASIS_MACHINE in

#-- FIXME: Currently only handle 1 node jobs. Need to have smarter logic to 
#          figure out the number of nodes needed based on NPROCS and NTHREADS


Mars*)

  JOB_SIZE=$(( ( $NTHREADS * $NPROCS / 32 + 1 ) * 32 )) 
  
  cat > ${EXEC}.pbs << EOF
#PBS -l size=$JOB_SIZE,walltime=1:00:00
#PBS -A UT-SUPPORT
#PBS -N $EXEC
#PBS -j eo

set -o verbose
cd \$PBS_O_WORKDIR

export OMP_NUM_THREADS=$NTHREADS
aprun -n $NPROCS -d $NTHREADS ./${EXEC}_${GENASIS_MACHINE} $EXEC_ARGS

rm .running
EOF
;;
Beacon*)
  
  #-- Beacon PBS needs work
  N_NODES=$(( $NTHREADS * $NPROCS / 16 + 1 )) 

  cat > ${EXEC}.pbs << EOF
#PBS -l nodes=$N_NODES,walltime=1:00:00
#PBS -A UT-SUPPORT
#PBS -N $EXEC
#PBS -j eo

set -o verbose
cd \$PBS_O_WORKDIR

export OMP_NUM_THREADS=$NTHREADS
mpiexec.hydra -np $NPROCS ./${EXEC}_${GENASIS_MACHINE} $EXEC_ARGS

rm .running
EOF
esac

JOBID=`qsub ${EXEC}.pbs | cut -d "." -f1`
echo $JOBID > .running

while true; do
  if [ -e .running ]
  then
    echo "Waiting for submitted job to finish ... "
    sleep 30
  else
    echo "Submitted job finished."
    echo "Displaying output ... " 
    break
  fi
done

cat ${EXEC}.e${JOBID}
ExtractRuntimeStatistics.sh $EXEC.e${JOBID} $JENKINS_BID
