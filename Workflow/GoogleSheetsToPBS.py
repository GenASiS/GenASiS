#!/usr/bin/env python

import os
import sys
import json
import math
import gspread
from oauth2client.service_account import ServiceAccountCredentials

scope = ['https://spreadsheets.google.com/feeds']

ThisDir = os.path.dirname(os.path.realpath(__file__))
credentials = ServiceAccountCredentials.from_json_keyfile_name \
                (os.path.join(ThisDir, 'GenASiS-GApps.json'), scope)
auth = gspread.authorize(credentials)

SheetTitle = sys.argv[1]
fp = auth.open(SheetTitle)

allWorksheet = fp.worksheets()
 
# If given, select only worksheet that match the argument 
worksheets = []
if len(sys.argv) > 2:
  for iArg in sys.argv[2:]:
    for iWS in allWorksheet:
      if iArg in iWS.title:
        worksheets.append(iWS)
else:
  worksheets = allWorksheet
  
Program        = SheetTitle.split("_")[0]
Dimensionality = SheetTitle.split("_")[1]

for iWS in worksheets:

  XRes    = iWS.title.split(" ")[0][1:]
  Machine = iWS.title.split(" ")[1]
  
  ss = iWS.get_all_values()

  ChartType             = None
  MaxLevels             = None
  RefinementThreshold   = None
  CoarseningThreshold   = None
  AdaptionSolveInterval = None
  nCellsArg             = None

  ML = 1

  #-- Assume 1st row is header
  for iR in range(1, len(ss)):

    if not any(ss[iR]):
      continue
    
    if "END" in ss[iR][0]:
      break
    
    if ss[iR][0]:
      nProcs = ss[iR][0]
    
    if ss[iR][1]:
      MaxLevels = ss[iR][1]
      ChartType = 'MULTI_LEVEL'
      ML = int(MaxLevels)
    
    if ss[iR][2]:
      RefinementThreshold = ss[iR][2]
    
    if ss[iR][3]:
      CoarseningThreshold = ss[iR][3]

    if ss[iR][4]:
      AdaptionSolveInterval = ss[iR][4]
    
    BaseRes = int(XRes)/int(math.pow(2,(int(ML)-1)))  
    if ChartType:
      if Dimensionality == '2D':
        nCellsArg = "nCellsBaseLevel=%s,%s" % (BaseRes, BaseRes)
      elif Dimensionality == '3D':
        nCellsArg = "nCellsBaseLevel=%s,%s,%s" % (BaseRes, BaseRes, BaseRes)
    else: 
      if Dimensionality == '2D':
        nCellsArg = "nCells=%s,%s" % (BaseRes, BaseRes)
      elif Dimensionality == '3D':
        nCellsArg = "nCells=%s,%s" % (BaseRes, BaseRes)
    
    #-- Build arguments
    cmd = "aprun -n %s ./%s_%s \\\nVerbosity=INFO_2 \\\nDimensionality=%s \\\n" \
           % (nProcs, Program, Machine, Dimensionality)
    
    if nCellsArg             : cmd += "%s \\\n" % nCellsArg
    if ChartType             : cmd += "ChartType=%s \\\n" % ChartType
    if MaxLevels             : cmd += "MaxLevels=%s \\\n" % MaxLevels
    if RefinementThreshold   : cmd += "RefinementThreshold=%s \\\n" % RefinementThreshold
    if CoarseningThreshold   : cmd += "CoarseningThreshold=%s \\\n" % CoarseningThreshold
    if AdaptionSolveInterval : cmd += "AdaptionSolveInterval=%s \\\n" % AdaptionSolveInterval
    
    cmd += "OutputDirectory=../Output_${PBS_JOBNAME}/ "
    cmd += " 2>&1 | tee ${PBS_JOBNAME}.out\n\n"
    
    if 'Mars' in Machine:
      jobsize = ( ( (int(nProcs)-1) / 32 ) + 1 ) * 32
      PBS_cmd = "#PBS -l size=%s,walltime=1:00:00\n" % jobsize
    elif 'Darter' in Machine:
      jobsize = ( ( (int(nProcs)-1) / 16 ) + 1 ) * 16
      PBS_cmd = "#PBS -l size=%s,walltime=1:00:00\n" % jobsize
    elif 'Beacon' in Machine:
      jobsize = ( (int(nProcs)-1) / 16 ) + 1
      PBS_cmd = "#PBS -l node=%s,walltime=1:00:00\n" % jobsize
    
    PBS_cmd = PBS_cmd + "#PBS -j eo\n"
    if not ChartType:
      PBS_Name = "%s_%s_%s_R%s_L0_RT00_CT00_A0_n%s" \
                  % ( Program, Machine, Dimensionality, XRes, nProcs)
    else:
      if not AdaptionSolveInterval:
        PBS_Name = "%s_%s_%s_R%s_L%s_RT00_CT00_A%s_n%s" \
                    % ( Program, Machine, Dimensionality, XRes, \
                        MaxLevels, '0', nProcs)
      else:
        RT = "%02d" % int(round(float(RefinementThreshold) * 100))
        CT = "%02d" % int(round(float(CoarseningThreshold) * 100))
        PBS_Name = "%s_%s_%s_R%s_L%s_RT%s_CT%s_A%s_n%s" \
                    % ( Program, Machine, Dimensionality, XRes, \
                        MaxLevels, RT, CT, AdaptionSolveInterval, nProcs)
      
    PBS_cmd = PBS_cmd + "#PBS -N %s\n" % PBS_Name
    PBS_cmd = PBS_cmd + "\nset -o verbose\n\ncd $PBS_O_WORKDIR\n\n"
    
#    Epilog  = "module load python/.2.7.11\n"
#    Epilog += "export PYTHONWARNINGS=\"ignore:Unverified HTTPS request\"\n"
#    Epilog += "%s/ProgramOutputToGoogleSheets.py ${PBS_JOBNAME}.out\n\n" % ThisDir
    Epilog = "cp ${PBS_JOBNAME}.out ../Output_${PBS_JOBNAME}/\n"
    
    fp = open("%s.pbs" % PBS_Name, 'w')
    fp.write(PBS_cmd)
    fp.write(cmd)
    fp.write(Epilog)
    fp.close()
    
    #print PBS_cmd
    #print cmd
    #print ''
