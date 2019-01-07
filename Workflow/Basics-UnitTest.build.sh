#-- Expect: 
#     - GENASIS_WORKSPACE set to location of GenASiS
#     - GENASIS_MACHINE set to get the right Makefile
#     - DEBUG or OPTIMIZE passed as an argument

set -o errexit  #-- exit on any error

PURPOSE=$1

echo "======================================================"
echo "       VARIABLE MANAGEMENT UNIT TESTS : $PURPOSE      "
echo "======================================================"

cd $GENASIS_WORKSPACE/Programs/UnitTests/Basics/VariableManagement/Specifiers/Executables
make PURPOSE=${PURPOSE} clobber all

cd ../../Devices/Executables
make PURPOSE=${PURPOSE} clobber all

cd ../../ArrayOperations/Executables
make PURPOSE=${PURPOSE} clobber all

cd ../../ArrayArrays/Executables
make PURPOSE=${PURPOSE} clobber all

cd ../../Storages/Executables
make PURPOSE=${PURPOSE} clobber all

echo "======================================================"
echo "            DISPLAY UNIT TESTS : $PURPOSE             "
echo "======================================================"

cd $GENASIS_WORKSPACE/Programs/UnitTests/Basics/Display/Executables
make PURPOSE=${PURPOSE} clobber all

echo "======================================================"
echo "        MESSAGE_PASSING UNIT TESTS : $PURPOSE         "
echo "======================================================"

cd $GENASIS_WORKSPACE/Programs/UnitTests/Basics/MessagePassing/MessagePassingBasics/Executables
make PURPOSE=${PURPOSE} clobber all

cd ../../PointToPoint/Executables
make PURPOSE=${PURPOSE} clobber all

cd ../../Collective/Executables
make PURPOSE=${PURPOSE} clobber all

echo "======================================================"
echo "          FILE SYSTEM UNIT TESTS : $PURPOSE           "
echo "======================================================"

cd $GENASIS_WORKSPACE/Programs/UnitTests/Basics/FileSystem/FileSystemBasics/Executables
make PURPOSE=${PURPOSE} clobber all

cd ../../GridImageBasics/Executables
make PURPOSE=${PURPOSE} clobber all

cd ../../CurveImages/Executables
make PURPOSE=${PURPOSE} clobber all

cd ../../PointGridImages/Executables
make PURPOSE=${PURPOSE} clobber all

cd ../../StructuredGridImages/Executables
make PURPOSE=${PURPOSE} clobber all

cd ../../UnstructuredGridImages/Executables
make PURPOSE=${PURPOSE} clobber all

echo "======================================================"
echo "           RUNTIME UNIT TESTS : $PURPOSE              "
echo "======================================================"

cd $GENASIS_WORKSPACE/Programs/UnitTests/Basics/Runtime/Executables
make PURPOSE=${PURPOSE} clobber all

