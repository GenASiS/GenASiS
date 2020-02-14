#include "Preprocessor"

submodule ( ConservationLawEvolution_Template ) ConservationLawEvolution_Kernel

  use Basics
  implicit none
  
contains


  module procedure ComputeTimeStepKernel
  
    select case ( nDimensions ) 
    case ( 1 ) 
      TimeStepLocal &
        = CellWidth ( 1 ) &
          / maxval ( max ( FEP_1, -FEM_1 ) )
    case ( 2 ) 
      TimeStepLocal &
        = minval ( CellWidth ( 1 : 2 ) ) &
          / maxval ( max ( FEP_1, FEP_2, -FEM_1, -FEM_2 ) )
    case ( 3 ) 
      TimeStepLocal &
        = minval ( CellWidth ) &
          / maxval ( max ( FEP_1, FEP_2, FEP_3, -FEM_1, -FEM_2, -FEM_3 ) )
    end select

  end procedure ComputeTimeStepKernel


end submodule ConservationLawEvolution_Kernel
