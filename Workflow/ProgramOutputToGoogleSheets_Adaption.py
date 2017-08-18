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
MaxLevels             = int(FnParts[5][1:])
RefinementThreshold   = FnParts[6][2:]
CoarseningThreshold   = FnParts[7][2:]
AdaptionSolveInterval = int(FnParts[8][1:])
nProcsStr             = FnParts[9]

scope = ['https://spreadsheets.google.com/feeds']

ThisDir = os.path.dirname(os.path.realpath(__file__))
credentials = ServiceAccountCredentials.from_json_keyfile_name \
                (os.path.join(ThisDir, 'GenASiS-GApps.json'), scope)
auth = gspread.authorize(credentials)

WorksheetTitle  = "%s %s %s" % (Resolution, Machine, nProcsStr)

print "Processing %s" % (Fn)

fp = auth.open(SheetTitle)
WS = fp.worksheet(WorksheetTitle)

#-- Collect MaxTimers: 

cmd = "tail -n 300 %s \
         | tac | grep -m1 -B25 \"Max timers\" | tac | grep \"%s\""

Execution = parseValue(capture(cmd % (Fn, "Execution")))
IO        = parseValue(capture(cmd % (Fn, "I/O")))
Evolution = parseValue(capture(cmd % (Fn, "Evolution")))
Adaption  = parseValue(capture(cmd % (Fn, "Adaption")))
Solve_1   = parseValue(capture(cmd % (Fn, "Solve_Level_01")))
Solve_2   = parseValue(capture(cmd % (Fn, "Solve_Level_02")))
Solve_3   = parseValue(capture(cmd % (Fn, "Solve_Level_03")))
Solve_4   = parseValue(capture(cmd % (Fn, "Solve_Level_04")))
Adapt_1   = parseValue(capture(cmd % (Fn, "Adapt_Level_01")))
Adapt_2   = parseValue(capture(cmd % (Fn, "Adapt_Level_02")))
Adapt_3   = parseValue(capture(cmd % (Fn, "Adapt_Level_03")))
Adapt_4   = parseValue(capture(cmd % (Fn, "Adapt_Level_04")))

#-- Collect memory usage
cmd = "tail -n 200 %s \
         | tac | grep -m1 -B9 \"memory usage\" | tac | grep \"%s\""
maxHWM = parseValue(capture(cmd % (Fn, "max HWM")))

#-- Collect proper cell ratio
cmd = "tail -n 200 %s \
         | tac | grep -m1 -B%d \"Normalized to nLevels\" | tac"
nProperRatio = capture(cmd % (Fn, MaxLevels)).split("\n")[1:]

#-- Find the row to update
ss = WS.get_all_values()

MaxLevelsRow             = 0
RefinementThresholdRow   = "00"
CoarseningThresholdRow   = "00"
AdaptionSolveIntervalRow = 0

for iR in range(1, len(ss)):

  if not any(ss[iR]):
    continue
    
  if ss[iR][1]:
    MaxLevelsRow = int(ss[iR][1])
  
  if ss[iR][2]:
    RefinementThresholdRow = "%02d" % int(round(float(ss[iR][2]) * 100))
  
  if ss[iR][3]:
    CoarseningThresholdRow = "%02d" % int(round(float(ss[iR][3]) * 100))
  
  if ss[iR][4]:
    AdaptionSolveIntervalRow = int(ss[iR][4])
  
  if (MaxLevelsRow == MaxLevels) \
      and (RefinementThresholdRow == RefinementThreshold) \
      and (CoarseningThresholdRow == CoarseningThreshold) \
      and (AdaptionSolveIntervalRow == AdaptionSolveInterval):
    RowNum = iR + 1 #-- spreadsheet cell counts reference starts from 1
    break

WS.update_cell(RowNum,  6,  Execution)
WS.update_cell(RowNum,  7,  maxHWM)
WS.update_cell(RowNum,  8,  IO)
WS.update_cell(RowNum,  9,  Evolution)
WS.update_cell(RowNum,  10, Adaption)

iL = 0
for iCol in range(12, 12+len(nProperRatio)):
  WS.update_cell(RowNum, iCol, parseValue(nProperRatio[iL]))
  iL = iL + 1

WS.update_cell(RowNum, 16, Solve_1)
WS.update_cell(RowNum, 17, Solve_2)
WS.update_cell(RowNum, 18, Solve_3)
WS.update_cell(RowNum, 19, Solve_4)
WS.update_cell(RowNum, 20, Adapt_1)
WS.update_cell(RowNum, 21, Adapt_2)
WS.update_cell(RowNum, 22, Adapt_3)
WS.update_cell(RowNum, 23, Adapt_4)
