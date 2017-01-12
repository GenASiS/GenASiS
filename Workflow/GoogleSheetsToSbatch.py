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

ProgramName = SheetTitle

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
  
for iWS in worksheets:
  ss = iWS.get_all_values()
  
  #-- Assume 1st row is header
  for iR in range(1, len(ss)):
    
    #-- skip empty row
    if not any(ss[iR]):
      continue
    
    nCells         = ss[iR][0]
    nNodes         = ss[iR][1]
    nThreads       = ss[iR][2]
    nBricks        = ss[iR][3]
    WallTime       = ss[iR][4]
    nCores         = ss[iR][5]
    Dimensionality = len(nCells.split(","))
    
    if not WallTime:
      continue
    
    
    JobName  = "%s_%sD_R%04d_NN%03d_T%03d" \
                  % ( ProgramName, Dimensionality, \
                      int(nCells.split(",")[0]), int(nNodes), \
                      int(nThreads) )
    print JobName
    
    S_Header  = "#!/bin/bash\n"
    S_Header += "#SBATCH --job-name=%s\n" % JobName
    S_Header += "#SBATCH --output=\"%s.out.JID%%j\"\n" % JobName
    S_Header += "#SBATCH --partition=compute\n"
    S_Header += "#SBATCH --nodes=%d\n" % int(nNodes)
    S_Header += "#SBATCH --ntasks-per-node=24\n"
    S_Header += "#SBATCH --export=ALL\n"
    S_Header += "#SBATCH -t %s\n\n" % WallTime
    
    S_Cmd  = "export OMP_NUM_THREADS=%d\n" % int(nThreads)
    S_Cmd += "ibrun -v ./%s NoWrite=T \\\n" % ProgramName
    S_Cmd += "Dimensionality=%dD nCells=%s nBricks=%s \\\n| tee %s.out\n" \
              % (Dimensionality, nCells, nBricks, JobName)
              
    fp = open("%s.sbatch" % JobName, 'w')
    fp.write(S_Header)
    fp.write(S_Cmd)
    fp.close()
