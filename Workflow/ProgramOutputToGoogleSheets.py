#!/usr/bin/env python

import os, sys, re, subprocess
import gspread
from oauth2client.service_account import ServiceAccountCredentials

def capture(cmd):
  """
  Capture standard out for a command.
  @param cmd: Command string or array.
  """
  p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

  return p.communicate()[0]
  
def parseValue(str):
  if str != "":
    return str.split("=")[1].strip().split(" ")[0]
  else:
    return ""

Fn         = sys.argv[1]
SheetTitle = sys.argv[2]

(FnBase, FnExt) = os.path.splitext(os.path.basename(Fn))
FnParts = FnBase.split("_")

ProgramName           = FnParts[0]
Machine               = FnParts[1] + "_" + FnParts[2]
Dimensionality        = FnParts[3] 
Resolution            = FnParts[4]
Nodes                 = FnParts[5]
Threads               = FnParts[6]

scope = ['https://spreadsheets.google.com/feeds']

ThisDir = os.path.dirname(os.path.realpath(__file__))
credentials = ServiceAccountCredentials.from_json_keyfile_name \
                (os.path.join(ThisDir, 'GenASiS-GApps.json'), scope)
auth = gspread.authorize(credentials)

#WorksheetTitle  = "%s_%s_%s" % (ProgramName, Machine)
WorksheetTitle  = "Sheet1"

print "Processing %s" % (Fn)

fp = auth.open(SheetTitle)
WS = fp.worksheet(WorksheetTitle)

#-- Collect MaxTimers: 

cmd = "tail -n 300 %s \
         | tac | grep -m1 -B15 \"Max timers\" | tac | grep \"%s\""
cmd2 = "grep \"L1 error\" %s"

Execution = parseValue(capture(cmd % (Fn, "Execution")))
L1Error   = parseValue(capture(cmd2 % (Fn)))

#-- Collect memory usage
#cmd = "tail -n 200 %s \
#         | tac | grep -m1 -B9 \"memory usage\" | tac | grep \"%s\""
#maxHWM = parseValue(capture(cmd % (Fn, "max HWM")))

#-- Collect proper cell ratio
#cmd = "tail -n 200 %s \
#         | tac | grep -m1 -B%d \"Normalized to nLevels\" | tac"
#nProperRatio = capture(cmd % (Fn, MaxLevels)).split("\n")[1:]

#-- Find the row to update

ss = WS.get_all_values()
nCells   = "%d,%d,%d" % (int(Resolution[1:]), int(Resolution[1:]), int(Resolution[1:]))
nNodes   = "%d" % int(Nodes[2:])
nThreads = "%d" % int(Threads[1:])

#-- Assume 1st row is header
for iR in range(1, len(ss)):

  #-- skip empty row
  if not any(ss[iR]):
    continue
  
  nCellsRow   = ss[iR][0]
  nNodesRow   = ss[iR][1]
  nThreadsRow = ss[iR][2]
  
  if (nCellsRow == nCells) and (nNodesRow == nNodes) \
     and (nThreadsRow == nThreads):
    RowNum = iR + 1 #-- spreadsheet cell counts reference starts from 1
    break
  

WS.update_cell(RowNum, 7, Execution)
WS.update_cell(RowNum, 8, L1Error)

