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


  module procedure ComputeTimeStepKernel_OMP
  
    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      lV, uV
    real ( KDR ) :: &
      MaxSpeed
      
    lV = 1
    where ( shape ( FEP_1 ) > 1 )
      lV = oV + 1
    end where
    
    uV = 1
    where ( shape ( FEP_1 ) > 1 )
      uV = shape ( FEP_1 ) - oV
    end where
  
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
    MaxSpeed = - huge ( 1.0_KDR )
    !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV, jV, kV) &
    !$OMP& reduction ( max : MaxSpeed )
    do kV = lV ( 3 ) , uV ( 3 )
      do jV = lV ( 2 ), uV ( 2 )    
        do iV = lV ( 1 ), uV ( 1 )
          MaxSpeed &
            = max (  FEP_1 ( iV, jV, kV ),  FEP_2 ( iV, jV, kV ), &
                     FEP_3 ( iV, jV, kV ), -FEM_1 ( iV, jV, kV ), & 
                    -FEM_2 ( iV, jV, kV ), -FEM_3 ( iV, jV, kV ), MaxSpeed  )
        end do
      end do
    end do
    !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    !call Show ( MaxSpeed, 'MaxSpeed OMP', CONSOLE % INFO_2 )
    !call Show &
    !       ( maxval ( max ( FEP_1, FEP_2, FEP_3, -FEM_1, -FEM_2, -FEM_3 ) ), &
    !         'MaxSpeed Fortran', CONSOLE % INFO_2 )
    
    TimeStepLocal = minval ( CellWidth ) / MaxSpeed
    end select

  end procedure ComputeTimeStepKernel_OMP


end submodule ConservationLawEvolution_Kernel
