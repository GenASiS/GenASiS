JobNameBase = 'SineWaveAdvection_3D'
ProgramName = 'SineWaveAdvection_Beacon_Intel'
Project     = 'UT-TENN0245'
Resolution  = 160
nNodes      = [ 1, 2, 4, 8 ]
nThreads    = 1
WallTime    = [ '0:20:00', '0:10:00', '0:10:00', '0:10:00' ]

nJobs = len ( nNodes )

nCoresNode = 16
nCores = [ nN * nCoresNode for nN in nNodes ]

nProcesses = [ nC / nThreads for nC in nCores ]

from collections import deque
from numpy import prod

JobNumber = range ( nJobs )
for iJob in JobNumber:

  JobName  = JobNameBase
  JobName += '_R{0:04d}'  . format ( Resolution )
  JobName += '_NN{0:03d}' . format ( nNodes [ iJob ] )
  JobName += '_T{0:03d}'  . format ( nThreads )

  Lines  = '#!/bin/bash\n'
  Lines += '#PBS -l nodes={0},walltime={1}\n' \
           . format ( nNodes [ iJob ], WallTime [ iJob ] )
  Lines += '#PBS -A {0}\n' . format ( Project ) 
  Lines += '#PBS -j oe\n'
  Lines += '#PBS -N {0}\n\n' . format ( JobName ) 

  Lines += 'set -o verbose\n'
  Lines += 'date\n\n'

  Lines += 'export OMP_NUM_THREADS={0}\n\n' . format ( nThreads )

  Lines += 'cd $PBS_O_WORKDIR\n'
  Lines += 'svn info\n\n'

  Lines += 'echo $PBS_JOBNAME\n'
  Lines += 'cat $PBS_NODEFILE\n\n'

  nBricks = [ 1, 1, 1 ]
  iaDimension = deque([ 2, 1, 0 ])
  while ( prod ( nBricks ) < nProcesses [ iJob ] ):
    nBricks [ iaDimension [ 0 ] ] *= 2
    iaDimension.rotate ( -1 )  

  print 'nProcesses = ', nProcesses [ iJob ]
  print 'nBricks = ', nBricks

  Lines += 'mpirun -n {0} ./{1} NoWrite=T \\\n' \
           . format ( nProcesses [ iJob ], ProgramName )
  Lines += 'Dimensionality=3D \\\n'
  Lines += 'nCells={0},{1},{2} \\\n' \
           . format ( Resolution, Resolution, Resolution )
  Lines += 'nBricks={0},{1},{2} \\\n' \
           . format ( nBricks [ 0 ], nBricks [ 1 ], nBricks [ 2 ] )
  Lines += 'OutputDirectory=../Output_{0}/$PBS_JOBNAME/ \\\n' \
           . format ( ProgramName )
  Lines += '> $PBS_JOBNAME.$PBS_JOBID.stdout\n'

  file = open ( JobName + '.pbs', 'w' )
  file.write ( Lines )
  file.close ( ) 
