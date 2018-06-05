import sys
import os

ApplicationName = sys.argv[-1]

Database = 'localhost:' + ApplicationName + '_MultiMesh_*.silo database'
OpenDatabase(Database)

if ApplicationName[-2:] == '3D':
  AddPlot("Pseudocolor", "Chart_01/Fluid/ComovingDensity", 1, 1)
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
  AddPlot("Pseudocolor", "Chart_01/Fluid/ComovingDensity", 1, 1)
  
  if not ("RayleighTaylor" in ApplicationName) \
     and not ("SedovTaylor" in ApplicationName) \
     and not ("FishboneMoncrief" in ApplicationName): 
    AddOperator("Elevate", 1)
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

  if ("FishboneMoncrief" in ApplicationName):
    AddOperator("Transform", 1)
    TA = TransformAttributes ( )
    TA.transformType = TA.Coordinate
    TA.inputCoordSys = TA.Spherical
    TA.outputCoordSys = TA.Rectangular
    SetOperatorOptions(TA, 1)
    
    AddOperator("Project", 1)
    PA = ProjectAttributes ( )
    PA.projectionType = PA.XZRectangular
    SetOperatorOptions(PA, 1)

    p = PseudocolorAttributes ( )
    p.scaling = p.Log
    p.min, p.minFlag = 10.0, 1
    p.max, p.maxFlag = 1.e12, 1
    p.colorTableName = "bluehot"
    SetPlotOptions(p)

    DrawPlots()
           
    View2DAtts = GetView2D()
    View2DAtts.fullFrameActivationMode = View2DAtts.Off
    SetView2D(View2DAtts)
    
  DrawPlots()
  
elif ApplicationName[-2:] == '1D':
  AddPlot("Curve", "Chart_01/Fluid/ComovingDensity", 1, 1)
  DrawPlots()
  
CreateAnnotationObject("TimeSlider")

for state in range ( 0, TimeSliderGetNStates() ):
  SetTimeSliderState(state)
  SaveWindow()
