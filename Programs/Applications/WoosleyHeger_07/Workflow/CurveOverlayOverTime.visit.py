
#-- Overlay nPlots number of curves and modify its Hue-ness for each new plot
#-- nPlots

import os, glob, colorsys
from operator import div

WI = GetWindowInformation()

try:
  nPlots
except NameError:
  nPlots         = 10
  
OutputDir      = os.path.dirname((WI.activeSource.split(':')[1]))
OutputDirDB    = os.path.dirname(WI.activeSource)
OutputFileBase = os.path.basename(WI.activeSource)[:-13]

os.chdir(OutputDir)
nFiles = len(glob.glob1(".", "*.silo"))

#Step    = nFiles / nPlots

Checkpoint_1 = 1200
Checkpoint_2 = 8503
Step = ( Checkpoint_2 - Checkpoint_1 ) / nPlots

Hue_1 = 200
Hue_2 = 300        
HueStep = (Hue_2 - Hue_1) / nPlots

Hue = Hue_1
#for iPlot in range(0 + Step, nFiles, Step):
for iPlot in range(Checkpoint_1 + Step, Checkpoint_2, Step):

  Hue += HueStep
  print "Plotting file : %07d" % iPlot
  print "Setting Hue to: %d" % Hue
  
  File = "%s/%s_%07d.silo" % (OutputDirDB, OutputFileBase, iPlot)
  OverlayDatabase(File)
  CA = GetPlotOptions()
  RGB_Value = tuple(map(div, CA.curveColor[:3], (255., 255., 255.)))
  HSV_Value = colorsys.rgb_to_hsv(RGB_Value[0], RGB_Value[1], RGB_Value[2])
  RGB_Value = colorsys.hsv_to_rgb(HSV_Value[0] + (HueStep / 360.), HSV_Value[1], HSV_Value[2])
  
  print RGB_Value
  CA.curveColor = (int(RGB_Value[0]*255), int(RGB_Value[1]*255), int(RGB_Value[2]*255), 255)
  print CA.curveColor
  CA.showLegend = 0
  SetPlotOptions(CA)
  
  
  
  
    
