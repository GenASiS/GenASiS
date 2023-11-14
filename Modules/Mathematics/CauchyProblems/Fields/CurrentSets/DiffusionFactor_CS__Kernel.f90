#include "Preprocessor"

submodule ( DiffusionFactor_CS__Form ) DiffusionFactor_CS__Kernel
  
  use Basics
  
  implicit none
  
contains


  module procedure ComputeKernel

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( DF )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV
        DF ( iV )  =  1.0_KDR
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV  =  1,  nV
        DF ( iV )  =  1.0_KDR
      end do
      !$OMP end parallel do
    end if

  end procedure ComputeKernel


  module procedure Compute_I_CGS_Kernel

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVM, &
      lV, uV
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    lV  =  1
    where ( shape ( DF )  >  1 )
      lV  =  oV
    end where
    
    uV  =  1
    where ( shape ( DF )  >  1 )
      uV  =  shape ( DF )  -  oV
    end where
    uV ( iD )  =  size ( DF, dim = iD )  -  oV  +  1 
      
    iaS  =  0
    iaS ( iD )  =  1
    
    if ( UseDevice ) then
    
      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP private ( iaVM )
      do kV  =  lV ( 3 ),  uV ( 3 ) 
        do jV  =  lV ( 2 ),  uV ( 2 )
          do iV  =  lV ( 1 ),  uV ( 1 )

              iaVM  =  [ iV, jV, kV ]  -  iaS

              DF_I ( iV, jV, kV )  &
                =  min ( 1.0_KDR, &
                         max ( DF ( iV,         jV,         kV         ), &
                               DF ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) ) ) )

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

    else

      !$OMP parallel do collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP private ( iaVM )
      do kV  =  lV ( 3 ),  uV ( 3 ) 
        do jV  =  lV ( 2 ),  uV ( 2 )
          do iV  =  lV ( 1 ),  uV ( 1 )

              iaVM  =  [ iV, jV, kV ]  -  iaS

              DF_I ( iV, jV, kV )  &
                =  min ( 1.0_KDR, &
                         max ( DF ( iV,         jV,         kV         ), &
                               DF ( iaVM ( 1 ), iaVM ( 2 ), iaVM ( 3 ) ) ) )

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end parallel do

    end if

  end procedure Compute_I_CGS_Kernel


end submodule DiffusionFactor_CS__Kernel
