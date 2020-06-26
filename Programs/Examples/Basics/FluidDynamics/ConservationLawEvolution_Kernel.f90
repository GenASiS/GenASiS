#include "Preprocessor"

submodule ( ConservationLawEvolution_Template ) ConservationLawEvolution_Kernel

  use Basics
  implicit none
  
contains


  module procedure ComputeTimeStepKernel
  
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
    
    MaxSpeed = - huge ( 1.0_KDR )

    !$OMP  OMP_TARGET_DIRECTIVE parallel do simd collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV, jV, kV ) &
    !$OMP& reduction ( max : MaxSpeed ) map ( tofrom: MaxSpeed )
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
    !$OMP  end OMP_TARGET_DIRECTIVE parallel do simd
    
    TimeStepLocal = minval ( CellWidth ( 1 : nDimensions ) ) / MaxSpeed

  end procedure ComputeTimeStepKernel


end submodule ConservationLawEvolution_Kernel
