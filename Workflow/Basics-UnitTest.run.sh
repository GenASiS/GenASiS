#-- Expect: 
#     - GENASIS_WORKSPACE set to location of GenASiS
#     - GENASIS_MACHINE set to get the right Makefile
#     - GENASIS_LAUNCHER set as the job launcher with arguments
#       e.g. 'export GENASIS_LAUNCHER="aprun -n1 "

set -o errexit  #-- exit on any error

function RunExecutables 
  {
  for EXEC in $(find -maxdepth 1 -name "*_${GENASIS_MACHINE}"); do
    echo "============= Running $EXEC =============="
    $GENASIS_LAUNCHER $EXEC
    echo "=============== Done ====================="
    echo ""
  done
  }


echo "======================================================"
echo "               SPECIFIERS UNIT TESTS                  "
echo "======================================================"
cd $GENASIS_WORKSPACE/Programs/UnitTests/Basics/Specifiers/Executables
RunExecutables

echo "======================================================"
echo "                  DEVICES UNIT TESTS                  "
echo "======================================================"
cd $GENASIS_WORKSPACE/Programs/UnitTests/Basics/Devices/Executables
RunExecutables

echo "======================================================"
echo "         DATA MANAGEMENT UNIT TESTS : $PURPOSE        "
echo "======================================================"

cd $GENASIS_WORKSPACE/Programs/UnitTests/Basics/DataManagement/ArrayOperations/Executables
RunExecutables

cd ../../ArrayArrays/Executables
RunExecutables

cd ../../Storages/Executables
RunExecutables

echo "======================================================"
echo "                  DISPLAY UNIT TESTS                  "
echo "======================================================"

cd $GENASIS_WORKSPACE/Programs/UnitTests/Basics/Display/Executables
RunExecutables

echo "======================================================"
echo "             MESSAGE_PASSING UNIT TESTS               "
echo "======================================================"

cd $GENASIS_WORKSPACE/Programs/UnitTests/Basics/MessagePassing/MessagePassingBasics/Executables
RunExecutables

cd ../../PointToPoint/Executables
RunExecutables

cd ../../Collective/Executables
RunExecutables

echo "======================================================"
echo "                FILE SYSTEM UNIT TESTS                "
echo "======================================================"

cd $GENASIS_WORKSPACE/Programs/UnitTests/Basics/FileSystem/FileSystemBasics/Executables
RunExecutables

cd ../../GridImageBasics/Executables
RunExecutables

cd ../../CurveImages/Executables
RunExecutables

cd ../../PointGridImages/Executables
RunExecutables

cd ../../StructuredGridImages/Executables
RunExecutables

cd ../../UnstructuredGridImages/Executables
RunExecutables

echo "======================================================"
echo "                  RUNTIME UNIT TESTS                  "
echo "======================================================"

cd $GENASIS_WORKSPACE/Programs/UnitTests/Basics/Runtime/Executables
RunExecutables
