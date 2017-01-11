import sys
import os

ApplicationName = sys.argv[-1]

Database = 'localhost:' + ApplicationName + '_MultiMesh_*.silo database'
OpenDatabase(Database)

FluidType = 'PressurelessFluid'
if ApplicationName[:-3] == 'RiemannProblem' :
  FluidType = 'PolytropicFluid'


if ApplicationName[-2:] == '3D':
  AddPlot("Pseudocolor", "Grid/%s/ComovingDensity"%FluidType, 1, 1)
  AddOperator("Box", 1)
  BoxAtts = BoxAttributes()
  BoxAtts.amount = BoxAtts.Some  # Some, All
  BoxAtts.minx = 0.01
  BoxAtts.maxx = 0.99
  BoxAtts.miny = 0.01
  BoxAtts.maxy = 0.99
  BoxAtts.minz = 0.01
  BoxAtts.maxz = 0.99
  BoxAtts.inverse = 0
  SetOperatorOptions(BoxAtts, 1)
  DrawPlots()

  View3DAtts = View3DAttributes()
  View3DAtts.viewNormal = (0.665257, 0.59783, 0.44725)
  View3DAtts.focus = (0.5, 0.5, 0.5)
  View3DAtts.viewUp = (-0.339132, -0.291716, 0.894367)
  View3DAtts.viewAngle = 30
  View3DAtts.parallelScale = 0.893089
  View3DAtts.nearPlane = -1.78618
  View3DAtts.farPlane = 1.78618
  View3DAtts.imagePan = (0, 0)
  View3DAtts.imageZoom = 1
  View3DAtts.perspective = 1
  View3DAtts.eyeAngle = 2
  View3DAtts.centerOfRotationSet = 0
  View3DAtts.centerOfRotation = (0.5, 0.5, 0.5)
  View3DAtts.axis3DScaleFlag = 0
  View3DAtts.axis3DScales = (1, 1, 1)
  View3DAtts.shear = (0, 0, 1)
  SetView3D(View3DAtts)

elif ApplicationName[-2:] == '2D':
  AddPlot("Pseudocolor", "Grid/%s/ComovingDensity"%FluidType, 1, 1)
  AddOperator("Elevate", 1)
  AddOperator("Box", 1)
  
  ElevateAtts = ElevateAttributes()
  ElevateAtts.useXYLimits = 1
  ElevateAtts.limitsMode = ElevateAtts.OriginalData  # OriginalData, CurrentPlot
  ElevateAtts.scaling = ElevateAtts.Linear  # Linear, Log, Skew
  ElevateAtts.skewFactor = 1
  ElevateAtts.minFlag = 0
  ElevateAtts.min = 0
  ElevateAtts.maxFlag = 0
  ElevateAtts.max = 1
  ElevateAtts.zeroFlag = 0
  ElevateAtts.variable = "default"
  SetOperatorOptions(ElevateAtts, 1)
  
  BoxAtts = BoxAttributes()
  BoxAtts.amount = BoxAtts.All  # Some, All
  BoxAtts.minx = 0.01
  BoxAtts.maxx = 0.99
  BoxAtts.miny = 0.01
  BoxAtts.maxy = 0.99
  BoxAtts.minz = 0.0
  BoxAtts.maxz = 1.0
  BoxAtts.inverse = 0
  SetOperatorOptions(BoxAtts, 1)
  
  DrawPlots()
  
  View3DAtts = View3DAttributes()
  View3DAtts.viewNormal = (0.372646, -0.805585, 0.460617)
  View3DAtts.focus = (0.5, 0.5, 0.515625)
  View3DAtts.viewUp = (-0.196626, 0.416553, 0.887593)
  View3DAtts.viewAngle = 30
  View3DAtts.parallelScale = 0.892909
  View3DAtts.nearPlane = -1.78582
  View3DAtts.farPlane = 1.78582
  View3DAtts.imagePan = (0, 0)
  View3DAtts.imageZoom = 1
  View3DAtts.perspective = 1
  View3DAtts.eyeAngle = 2
  View3DAtts.centerOfRotationSet = 0
  View3DAtts.centerOfRotation = (0.5, 0.5, 0.515625)
  View3DAtts.axis3DScaleFlag = 0
  View3DAtts.axis3DScales = (1, 1, 1)
  View3DAtts.shear = (0, 0, 1)
  SetView3D(View3DAtts)

elif ApplicationName[-2:] == '1D':
  AddPlot("Curve", "Curves/%s/ComovingDensity"%FluidType, 1, 1)
  DrawPlots()
  
CreateAnnotationObject("TimeSlider")

for state in range ( 0, TimeSliderGetNStates() ):
  SetTimeSliderState(state)

