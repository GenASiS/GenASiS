#include "Preprocessor"

submodule ( ConservationLawStep_Form ) ConservationLawStep_Kernel

  use Basics
  implicit none
    
contains

  
  module procedure ComputeDifferencesKernel

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &   
      lV, uV
      
    call AssociateHost ( D_V, V )
    call AssociateHost ( D_dV_Left, dV_Left )
    call AssociateHost ( D_dV_Right, dV_Right )

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
    
!    dV_Left  = V - cshift ( V, shift = -1, dim = iD )    

    iaS = 0
    iaS ( iD ) = -1
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iaVS )
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
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
!    dV_Right = cshift ( V, shift = 1, dim = iD ) - V

    iaS = 0
    iaS ( iD ) = +1
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iaVS )
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
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( dV_Right )
    call DisassociateHost ( dV_Left )
    call DisassociateHost ( V )
    
  end procedure ComputeDifferencesKernel


  module procedure ComputeReconstructionKernel

    !real ( KDR ), dimension ( size ( V ) ) :: &
    !  dV
    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      dV
      
    call AssociateHost ( D_V, V )
    call AssociateHost ( D_dV_Left, dV_Left )
    call AssociateHost ( D_dV_Right, dV_Right )
    call AssociateHost ( D_V_Inner, V_Inner )
    call AssociateHost ( D_V_Outer, V_Outer )

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
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( dV )
    do iV = 1, size ( V )
      dV = ( sign ( 0.5_KDR, dV_Left ( iV ) ) &
             + sign ( 0.5_KDR, dV_Right ( iV ) ) ) &
             * min ( abs ( Theta * dV_Left ( iV ) ), &
                     abs ( Theta * dV_Right ( iV ) ), &
                     abs ( 0.5_KDR * ( dV_Left ( iV ) + dV_Right ( iV ) ) ) )
      V_Inner ( iV ) = V ( iV ) - 0.5_KDR * dV
      V_Outer ( iV ) = V ( iV ) + 0.5_KDR * dV
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( V_Outer )
    call DisassociateHost ( V_Inner )
    call DisassociateHost ( dV_Right )
    call DisassociateHost ( dV_Left )
    call DisassociateHost ( V )
    
  end procedure ComputeReconstructionKernel


  module procedure ComputeFluxesKernel

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &   
      lV, uV
      
    call AssociateHost ( D_AP_I, AP_I )
    call AssociateHost ( D_AP_O, AP_O )
    call AssociateHost ( D_AM_I, AM_I )
    call AssociateHost ( D_AM_O, AM_O )
    call AssociateHost ( D_RF_I, RF_I )
    call AssociateHost ( D_RF_O, RF_O )
    call AssociateHost ( D_U_I, U_I )
    call AssociateHost ( D_U_O, U_O )
    call AssociateHost ( D_F_I, F_I )
    call AssociateHost ( D_F_O, F_O )
            
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
    
    iaS = 0
    iaS ( iD ) = -1
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iaVS )
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
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    !where ( AP_O + AM_O > 0.0_KDR )
    !  F_O = ( AP_O * RF_O  +  AM_O * cshift ( RF_I, shift = +1, dim = iD ) &
    !          - AP_O * AM_O * ( cshift ( U_I, shift = +1, dim = iD ) - U_O ) ) &
    !        / ( AP_O + AM_O )
    !elsewhere
    !  F_O = 0.0_KDR
    !end where
    
    iaS = 0
    iaS ( iD ) = +1
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iaVS )
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
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( F_O )
    call DisassociateHost ( F_I )
    call DisassociateHost ( U_O )
    call DisassociateHost ( U_I )
    call DisassociateHost ( RF_O )
    call DisassociateHost ( RF_I )
    call DisassociateHost ( AM_O )
    call DisassociateHost ( AM_I )
    call DisassociateHost ( AP_O )
    call DisassociateHost ( AP_I )
    
  end procedure ComputeFluxesKernel


  module procedure ComputeUpdateKernel
    
    integer ( KDI ) :: &
      iV
      
    call AssociateHost ( D_dU, dU )
    call AssociateHost ( D_F_I, F_I )
    call AssociateHost ( D_F_O, F_O )

    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, size ( dU )
      dU ( iV ) = dU ( iV ) - dT * ( F_O ( iV ) - F_I ( iV ) ) * A / V
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( F_O )
    call DisassociateHost ( F_I )
    call DisassociateHost ( dU )

  end procedure ComputeUpdateKernel
  
  
  module procedure AddUpdateKernel
    
    integer ( KDI ) :: &
      iV, &
      nV
    
    call AssociateHost ( D_O, O )
    call AssociateHost ( D_U, U )
    call AssociateHost ( D_C, C )
    
    nV = size ( O )
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, nV
      C ( iV ) = O ( iV ) + U ( iV )
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( C )
    call DisassociateHost ( U )
    call DisassociateHost ( O )

  end procedure AddUpdateKernel
  
  
  module procedure CombineUpdatesKernel
    
    integer ( KDI ) :: &
      iV, &
      nV
    
    call AssociateHost ( D_C, C )    
    call AssociateHost ( D_O, O )
    call AssociateHost ( D_U, U )
    
    nV = size ( O )
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, nV
      C ( iV ) = 0.5_KDR * ( O ( iV ) + ( C ( iV ) + U ( iV ) ) )
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( U )
    call DisassociateHost ( O )
    call DisassociateHost ( C )

  end procedure CombineUpdatesKernel
  

end submodule ConservationLawStep_Kernel
