#include "Preprocessor"

submodule ( ConservationLawStep_Form ) ConservationLawStep_Kernel

  use iso_c_binding
  use Basics
  implicit none
  
  include 'ConservationLawStep_Interface.f90'
    
contains

  
  module procedure ComputeDifferencesKernel

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &   
      lV, uV, &
      nSizes
      
    lV = 1
    where ( shape ( V ) > 1 )
      lV = oV + 1
    end where
    lV ( iD ) = oV

    uV = 1
    where ( shape ( V ) > 1 )
      uV = shape ( V ) - oV
    end where
    uV ( iD ) = size ( V, dim = iD ) - 1
    
    if ( UseDirectDevice ) then
      iaS = 0
      nSizes = shape ( V )
      !$OMP target data use_device_ptr ( V, dV_Left, dV_Right )
      call ComputeDifferences_C & 
             ( c_loc ( V ), lV, uV, iaS, iD, &
               c_loc ( dV_Left ), c_loc ( dV_Right ), nSizes )
      !$OMP end target data
      return 
    end if
    
!    dV_Left  = V - cshift ( V, shift = -1, dim = iD )    

    iaS = 0
    iaS ( iD ) = -1
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do simd collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iaVS )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )
        
          !call Show ( [ iV, jV, kV ], '<<< iV, jV, kV: 1' )

          iaVS = [ iV, jV, kV ] + iaS

          dV_Left ( iV, jV, kV )  &
            =  V ( iV, jV, kV )  -  V ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )

        end do !-- iV
      end do !-- jV
    end do !-- kV
    !$OMP end OMP_TARGET_DIRECTIVE parallel do simd
    
!    dV_Right = cshift ( V, shift = 1, dim = iD ) - V

    iaS = 0
    iaS ( iD ) = +1
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do simd collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iaVS )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )
          
          !call Show ( [ iV, jV, kV ], '<<< iV, jV, kV: 2' )

          iaVS = [ iV, jV, kV ] + iaS

          dV_Right ( iV, jV, kV )  &
            =  V ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) - V ( iV, jV, kV )

        end do !-- iV
      end do !-- jV
    end do !-- kV
    !$OMP end OMP_TARGET_DIRECTIVE parallel do simd
    
  end procedure ComputeDifferencesKernel


  module procedure ComputeReconstructionKernel

    !real ( KDR ), dimension ( size ( V ) ) :: &
    !  dV
    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      dV
      
    !where ( dV_Left > 0.0_KDR .and. dV_Right > 0.0_KDR )
    !  dV = min ( Theta * dV_Left, Theta * dV_Right, &
    !               0.5_KDR * ( dV_Left + dV_Right ) )
    !elsewhere ( dV_Left < 0.0_KDR .and. dV_Right < 0.0_KDR )
    !  dV = max ( Theta * dV_Left, Theta * dV_Right, &
    !               0.5_KDR * ( dV_Left + dV_Right ) )      
    !elsewhere
    !  dV = 0.0_KDR
    !endwhere

    !V_Inner = V - 0.5_KDR * dV
    !V_Outer = V + 0.5_KDR * dV
    
    if ( UseDirectDevice ) then
      !$OMP target data use_device_ptr &
      !$OMP   ( V, dV_Left, dV_Right, V_Inner, V_Outer )
      call ComputeReconstruction_C &
             ( c_loc ( V ), c_loc ( dV_Left ), c_loc ( dv_Right ), Theta, &
               c_loc ( V_Inner ), c_loc ( V_Outer ), size ( V ) )
      !$OMP end target data
      return
    end if
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do simd &
    !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( dV )
    do iV = 1, size ( V )
      dV = ( sign ( 0.5_KDR, dV_Left ( iV ) ) &
             + sign ( 0.5_KDR, dV_Right ( iV ) ) ) &
             * min ( abs ( Theta * dV_Left ( iV ) ), &
                     abs ( Theta * dV_Right ( iV ) ), &
                     abs ( 0.5_KDR * ( dV_Left ( iV ) + dV_Right ( iV ) ) ) )
      V_Inner ( iV ) = V ( iV ) - 0.5_KDR * dV
      V_Outer ( iV ) = V ( iV ) + 0.5_KDR * dV
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do simd
    
  end procedure ComputeReconstructionKernel


  module procedure ComputeFluxesKernel

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &   
      lV, uV, &
      nSizes
      
    lV = 1
    where ( shape ( F_I ) > 1 )
      lV = oV + 1
    end where
    lV ( iD ) = oV

    uV = 1
    where ( shape ( F_I ) > 1 )
      uV = shape ( F_I ) - oV
    end where
    uV ( iD ) = size ( F_I, dim = iD ) - 1
    
    !where ( AP_I + AM_I > 0.0_KDR )
    !  F_I = ( AP_I * cshift ( RF_O, shift = -1, dim = iD )  +  AM_I * RF_I &
    !          - AP_I * AM_I * ( U_I - cshift ( U_O, shift = -1, dim = iD ) ) ) &
    !        / ( AP_I + AM_I )
    !elsewhere
    !  F_I = 0.0_KDR
    !end where
    
    if ( UseDirectDevice ) then
      iaS = 0
      nSizes = shape ( F_I )
      
      !$OMP target data use_device_ptr &
      !$OMP   ( AP_I, AP_O, AM_I, AM_O, RF_I, RF_O, U_I, U_O, F_I, F_O )
      call ComputeFluxes_C &
             ( c_loc ( AP_I ), c_loc ( AP_O ), &
               c_loc ( AM_I ), c_loc ( AM_O ), &
               c_loc ( RF_I ), c_loc ( RF_O ), &
               c_loc ( U_I ), c_loc ( U_O ), lV, uV, iaS, iD, &
               c_loc ( F_I ), c_loc ( F_O ), nSizes )
      !$OMP end target data

      return
    end if
    
    
    iaS = 0
    iaS ( iD ) = -1
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do simd collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iaVS )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )
            
          iaVS = [ iV, jV, kV ] + iaS
          
          if ( AP_I ( iV, jV, kV ) + AM_I ( iV, jV, kV ) > 0.0_KDR ) then
            F_I ( iV, jV, kV ) &
              = ( AP_I ( iV, jV, kV ) &
                    * RF_O ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                  + AM_I ( iV, jV, kV ) * RF_I ( iV, jV, kV ) &
                  - AP_I ( iV, jV, kV ) * AM_I ( iV, jV, kV ) &
                    * ( U_I ( iV, jV, kV ) &
                        -  U_O ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) ) ) &
                / ( AP_I ( iV, jV, kV ) + AM_I ( iV, jV, kV ) )
          else
            F_I ( iV, jV, kV ) = 0.0_KDR
          end if

        end do
      end do
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do simd
    
    !where ( AP_O + AM_O > 0.0_KDR )
    !  F_O = ( AP_O * RF_O  +  AM_O * cshift ( RF_I, shift = +1, dim = iD ) &
    !          - AP_O * AM_O * ( cshift ( U_I, shift = +1, dim = iD ) - U_O ) ) &
    !        / ( AP_O + AM_O )
    !elsewhere
    !  F_O = 0.0_KDR
    !end where
    
    iaS = 0
    iaS ( iD ) = +1
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do simd collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iaVS )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )
            
          iaVS = [ iV, jV, kV ] + iaS
          
          if ( AP_O ( iV, jV, kV ) + AM_O ( iV, jV, kV ) > 0.0_KDR ) then
            F_O ( iV, jV, kV ) &
              = ( AP_O ( iV, jV, kV ) * RF_O ( iV, jV, kV ) &
                  +  AM_O ( iV, jV, kV ) &
                     * RF_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                  - AP_O ( iV, jV, kV ) * AM_O ( iV, jV, kV ) &
                    * ( U_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                        - U_O ( iV, jV, kV ) ) ) &
                / ( AP_O ( iV, jV, kV ) + AM_O ( iV, jV, kV ) )
          else
            F_O ( iV, jV, kV ) = 0.0_KDR
          end if
          
        end do
      end do
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do simd
    
  end procedure ComputeFluxesKernel


  module procedure ComputeUpdateKernel
    
    integer ( KDI ) :: &
      iV
      
    if ( UseDirectDevice ) then
      !$OMP target data use_device_ptr ( dU, F_I, F_O )
      call ComputeUpdate_C &
             ( c_loc ( dU ), c_loc ( F_I ), c_loc ( F_O ), &
               V, A, dT, size ( dU ) )
      !$OMP end target data
      return
    end if
      
    !$OMP  OMP_TARGET_DIRECTIVE parallel do simd &
    !$OMP& schedule ( OMP_SCHEDULE_TARGET )
    do iV = 1, size ( dU )
      dU ( iV ) = dU ( iV ) - dT * ( F_O ( iV ) - F_I ( iV ) ) * A / V
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do simd
    
  end procedure ComputeUpdateKernel
  
  
  module procedure AddUpdateKernel
    
    integer ( KDI ) :: &
      iV, &
      nV
    
    nV = size ( O )
    
    if ( UseDirectDevice ) then
      !$OMP target data use_device_ptr ( O, U, C )
      call AddUpdate_C ( c_loc ( O ), c_loc ( U ), c_loc ( C ), nV )
      !$OMP end target data
      return
    end if
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do simd &
    !$OMP& schedule ( OMP_SCHEDULE_TARGET )
    do iV = 1, nV
      C ( iV ) = O ( iV ) + U ( iV )
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do simd
    
  end procedure AddUpdateKernel
  
  
  module procedure CombineUpdatesKernel
    
    integer ( KDI ) :: &
      iV, &
      nV
    
    nV = size ( O )
    
    if ( UseDirectDevice ) then
      !$OMP target data use_device_ptr ( C, O, U )
      call CombineUpdates_C ( c_loc ( C ), c_loc ( O ), c_loc ( U ), nV )
      !$OMP end target data
      return
    end if
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do simd &
    !$OMP& schedule ( OMP_SCHEDULE_TARGET )
    do iV = 1, nV
      C ( iV ) = 0.5_KDR * ( O ( iV ) + ( C ( iV ) + U ( iV ) ) )
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do simd
    
  end procedure CombineUpdatesKernel
  

end submodule ConservationLawStep_Kernel
