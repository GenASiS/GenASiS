#!/bin/bash

#-- Script used by Jenkins to submit test job, wait on it, and extract runtime
#   data

#-- Order of Arguments are: JenkinsBuildID nProcs nThreads Exec ExecArgs
NPROCS=$1
NTHREADS=$2
shift 2
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
Linux*)
  
  JOB_NAME=$(basename $PWD) 
  
  #-- Beacon PBS needs work
  N_NODES=$(( $NTHREADS * $NPROCS / 16 + 1 )) 
  cat > ${JOB_NAME}.pbs << EOF
#PBS -l nodes=$N_NODES:ppn=16,walltime=1:00:00
#PBS -A UT-TENN0245
#PBS -N $JOB_NAME
#PBS -j eo

set -o verbose
cd \$PBS_O_WORKDIR
export RUNNING=\${PBS_O_WORKDIR}/.running

export OMP_NUM_THREADS=$NTHREADS

rm $JOB_NAME.out 

EXES=\$(find -type f -executable)
for I_EXE in \$EXES; do
  EXE_DIR=\$(dirname \$I_EXE)
  EXE=\$(basename \$I_EXE)
  
  cd \$EXE_DIR
  
  echo "Executing \${EXE} with $NPROCS MPI tasks and $NTHREADS threads ... " \
    | tee -a \$PBS_O_WORKDIR/$JOB_NAME.out
  
  mpirun -n $NPROCS ./\${EXE} $EXEC_ARGS 2>&1 \
    | tee -a \$RUNNING | tee -a \$PBS_O_WORKDIR/$JOB_NAME.out
  if [[ "\${PIPESTATUS[0]}" -eq 0 ]]; then  
    echo "\${EXE} Verified" >> \$PBS_O_WORKDIR/$JOB_NAME.out
  else
    echo "" >> \$PBS_O_WORKDIR/$JOB_NAME.out
    echo "\${EXE} Failed" >> \$PBS_O_WORKDIR/$JOB_NAME.out
    exit 1
  fi
  echo "" | tee -a \$PBS_O_WORKDIR/$JOB_NAME.out
  echo "--------------------------------------------------------------------" \
    | tee -a \$PBS_O_WORKDIR/$JOB_NAME.out
  
  cd \$PBS_O_WORKDIR
done

rm \$RUNNING
EOF
;;
Darter*)
  
  JOB_NAME=$(basename `dirname $PWD`)
  
  N_NODES=$(( $NTHREADS * $NPROCS / 16 + 1 )) 
  cat > ${JOB_NAME}.pbs << EOF
#PBS -l nodes=$N_NODES,walltime=10:00
#PBS -N $JOB_NAME
#PBS -j eo

set -o verbose
cd \$PBS_O_WORKDIR

export OMP_NUM_THREADS=$NTHREADS

rm $JOB_NAME.out 

for EXE in *_${GENASIS_MACHINE}; do
  echo "Executing \${EXE} with $NPROCS MPI tasks and $NTHREADS threads ... " | tee -a $JOB_NAME.out
  aprun -n $NPROCS ./\${EXE} $EXEC_ARGS 2>&1 \
    | tee -a .running | tee -a $JOB_NAME.out
  if [[ "\${PIPESTATUS[0]}" -eq 0 ]]; then  
    echo Verified >> $JOB_NAME.out
  else
    echo "" >> $JOB_NAME.out
    echo "\${EXE} Failed" >> $JOB_NAME.out
    exit 1
  fi
  echo "" | tee -a $JOB_NAME.out
  echo "--------------------------------------------------------------------" | tee -a $JOB_NAME.out
done

rm .running
EOF
;;
CCondo*)
  
  JOB_NAME=$(basename `dirname $PWD`)
  
  N_NODES=$(( $NTHREADS * $NPROCS / 32 + 1 )) 
  cat > ${JOB_NAME}.pbs << EOF
#PBS -l nodes=$N_NODES:ppn=32,walltime=10:00
#PBS -N $JOB_NAME
#PBS -j eo

set -o verbose
cd \$PBS_O_WORKDIR

export OMP_NUM_THREADS=$NTHREADS
module load PE-gnu

rm $JOB_NAME.out 

for EXE in *_${GENASIS_MACHINE}; do
  echo "Executing \${EXE} with $NPROCS MPI tasks and $NTHREADS threads ... " | tee -a $JOB_NAME.out
  mpirun -n $NPROCS ./\${EXE} $EXEC_ARGS 2>&1 \
    | tee -a .running | tee -a $JOB_NAME.out
  if [[ "\${PIPESTATUS[0]}" -eq 0 ]]; then  
    echo Verified >> $JOB_NAME.out
  else
    echo "" >> $JOB_NAME.out
    echo "\${EXE} Failed" >> $JOB_NAME.out
    exit 1
  fi
  echo "" | tee -a $JOB_NAME.out
  echo "--------------------------------------------------------------------" | tee -a $JOB_NAME.out
done

rm .running
EOF
esac

JOBID=`qsub -V ${JOB_NAME}.pbs | cut -d "." -f1`
echo $JOBID > .running

#-- grep line whose first word is number and assume that's the PBS jobid
if [ ! -e .running ]; then
  echo ".running file does not exist. Failing the test."
  exit -1
fi
jobid=$(grep "^\w[0-9]" .running)

#-- Infinite loop to wait until either $jobid runs or ceases to exist

while true; do

  qstat $jobid &> /dev/null

  DATE=`date +%Y-%m-%d:%H:%M:%S`

  #-- If job is missing for some reason, fail the test
  if [ $? -ne 0 ]; then
    echo "Job $jobid is missing or was deleted"
    echo "Failing this test "
    rm .running
    exit -2
  fi

  job_state=$(qstat -f $jobid | grep job_state | awk '{print $3}')
  if [ "$job_state" == "R" ]; then
    echo "Submitted job $jobid is running. Following .running ... "
    break
  else
    echo "$DATE: Waiting for job $jobid to start ... "
    sleep 10
  fi

done

#-- Start a subshell that check on the existence of .running file
#   && $jobid and terminate when its done. Then we use that subshell PID 
#   with "tail" command to let tail follow the output of the test
(
while true; do
  qstat $jobid &> /dev/null
  if [ -e .running ] && [ $? -eq 0 ]; then
    sleep 30
  else
    break
  fi
done
) &

trackPID=$!
tail -f .running --pid=$trackPID

#-- Display the actual output *.out
echo ""
echo "=================================================================="
echo "Displaying All Output:  "
echo ""
cat $JOB_NAME.out

echo "=================================================================="
echo "  Tests Summmary: "
echo ""
grep Executing $JOB_NAME.out
grep Failed $JOB_NAME.out
grep Verified $JOB_NAME.out

#ExtractRuntimeStatistics.sh $EXEC.e${JOBID} $JENKINS_BID
