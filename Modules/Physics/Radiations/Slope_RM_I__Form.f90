module Slope_RM_I__Form

  !-- Slope_RadiationMoments_Interactions__Form

  use Basics
  use Mathematics
  use Gravitations
  use RadiationMoments_BM__Form
  use Interactions_BM__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_RM_I_Form
    integer ( KDI ) :: &
      iEnergy_B
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum_B
    type ( CommunicatorForm ), pointer :: &
      Communicator_X_1D  =>  null ( )
    type ( CollectiveOperation_R_Form ), dimension ( : ), allocatable :: &
      CO_SplitSource
    class ( Slope_DFV_F_DT_Form ), pointer :: &
      Slope_DFV => null ( )
    class ( RadiationMoments_BM_Form ), pointer :: &
      Radiation => null ( )
    class ( Interactions_BM_Form ), pointer :: &
      Interactions => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_RM_I
    generic, public :: &
      Initialize => InitializeAllocate_RM_I
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_RM_I_Form

    private :: &
      ComputeSource_F

    private :: &
      ComputeKernel, &
      ComputeSimpleKernel

    interface

      module subroutine ComputeKernel &
               ( ProperCell, Xi_J, Xi_H, Chi_J, Chi_H, J, H_1, H_2, H_3, &
                 S_1, S_2, S_3, SF_RD, M_DD_11, M_DD_22, M_DD_33, &
                 S_S_1_D, S_S_2_D, S_S_3_D, dT, S_E, S_S_1, S_S_2, S_S_3, &
                 UseDeviceOption )
        use Basics
        implicit none
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          ProperCell
        real ( KDR ), dimension ( : ), intent ( in ) :: &
           Xi_J,  Xi_H, &
          Chi_J, Chi_H, &
          J, &
          H_1, H_2, H_3, &
          S_1, S_2, S_3, &
          SF_RD, &
          M_DD_11, M_DD_22, M_DD_33, &
          S_S_1_D, S_S_2_D, S_S_3_D  !-- Divergence
        real ( KDR ), intent ( in ) :: &
          dT
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          S_E, &
          S_S_1, S_S_2, S_S_3
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeKernel

      module subroutine ComputeSimpleKernel &
               ( ProperCell, Xi_J, Xi_H, Chi_J, Chi_H, J, H_1, H_2, H_3, &
                 M_DD_11, M_DD_22, M_DD_33, dT, S_E, S_S_1, S_S_2, S_S_3, &
                 UseDeviceOption )
        use Basics
        implicit none
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          ProperCell
        real ( KDR ), dimension ( : ), intent ( in ) :: &
           Xi_J,  Xi_H, &
          Chi_J, Chi_H, &
          J, &
          H_1, H_2, H_3, &
          M_DD_11, M_DD_22, M_DD_33
        real ( KDR ), intent ( in ) :: &
          dT
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          S_E, &
          S_S_1, S_S_2, S_S_3
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeSimpleKernel

    end interface


contains


  subroutine InitializeAllocate_RM_I ( S, R )

    class ( Slope_RM_I_Form ), intent ( inout ) :: &
      S
    class ( RadiationMoments_BM_Form ), intent ( in ), target :: &
      R

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_RM_I' 
    
    Name  =  trim ( R % Name ) // '_Slp_RM_I'

    S % Radiation  =>  R

!    select type ( I  =>  R % Interactions )
!    class is ( Interactions_BM_Form )
!      S % Interactions  =>  I
!    end select

    call Search ( R % iaBalanced, R % ENERGY_DENSITY_B, &
                  S % iEnergy_B )
    call Search ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_1, &
                  S % iMomentum_B ( 1 ) )
    call Search ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_2, &
                  S % iMomentum_B ( 2 ) )
    call Search ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_3, &
                  S % iMomentum_B ( 3 ) )

    call S % Slope_H_Form % Initialize &
           ( R % Atlas, &
             FieldOption = R % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = R % DeviceMemory, &
             PinnedMemoryOption = R % PinnedMemory, &
             DevicesCommunicateOption = R % DevicesCommunicate, &
             nFieldsOption = R % nBalanced, &
             IgnorabilityOption = R % IGNORABILITY + 1 )

  end subroutine InitializeAllocate_RM_I


  subroutine Compute ( S, dT, T_Option )

    class ( Slope_RM_I_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      dT
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum

    call Show ( 'Computing ' // trim ( S % Type ), S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )

    if ( .not. associated ( S % Interactions ) ) then
      call Show ( 'Please set Interactions', &
                  CONSOLE % ERROR )
      call Show ( 'Slope_RM_I__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    associate &
      (  I  =>  S % Interactions, &
         R  =>  S % Radiation )

    call I % Compute ( )

    do iC  =  1,  S % Atlas % nCharts
      select type ( C  =>  S % Atlas % Chart ( iC ) % Element )
        class is ( Chart_GS_Form )

      associate &
        ( IV  =>  I % Storage ( iC ) % Value, &
          RV  =>  R % Storage ( iC ) % Value, &
          SV  =>  S % Storage ( iC ) % Value )
      associate &
        (  Xi_J  =>  IV ( :, I % EMISSIVITY_J ), &
           Xi_H  =>  IV ( :, I % EMISSIVITY_H ), &
          Chi_J  =>  IV ( :, I % OPACITY_J ), &
          Chi_H  =>  IV ( :, I % OPACITY_H ), &
            J    =>  RV ( :, R % ENERGY_DENSITY_C ), &
            H_1  =>  RV ( :, R % MOMENTUM_DENSITY_C_U_1 ), &
            H_2  =>  RV ( :, R % MOMENTUM_DENSITY_C_U_2 ), &
            H_3  =>  RV ( :, R % MOMENTUM_DENSITY_C_U_3 ), &
            S_1  =>  RV ( :, R % MOMENTUM_DENSITY_B_D_1 ), &
            S_2  =>  RV ( :, R % MOMENTUM_DENSITY_B_D_2 ), &
            S_3  =>  RV ( :, R % MOMENTUM_DENSITY_B_D_3 ), &
          S_E    =>  SV ( :, S % iEnergy_B ), &
          S_S_1  =>  SV ( :, S % iMomentum_B ( 1 ) ), &
          S_S_2  =>  SV ( :, S % iMomentum_B ( 2 ) ), &
          S_S_3  =>  SV ( :, S % iMomentum_B ( 3 ) ) )

      select type ( G  =>  R % Geometry )
      class is ( Gravitation_G_Form )

        associate &
          ( GSV  =>  G % Storage ( iC ) % Value )
        associate &
          ( M_DD_11  =>  GSV ( :, G % METRIC_F_DD_11 ), &
            M_DD_22  =>  GSV ( :, G % METRIC_F_DD_22 ), &
            M_DD_33  =>  GSV ( :, G % METRIC_F_DD_33 ) )

        if ( associated ( S % Slope_DFV ) ) then

          call Search &
            ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_1, iMomentum ( 1 ) )
          call Search &
            ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_2, iMomentum ( 2 ) )
          call Search &
            ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_3, iMomentum ( 3 ) )
          
          associate &
            ( SDV  =>  S % Slope_DFV % Storage ( iC ) % Value )
          associate &
            ( SF_RD    =>   RV ( :, R % STRESS_FACTOR_RD ), &
              S_S_1_D  =>  SDV ( :, iMomentum ( 1 ) ), &
              S_S_2_D  =>  SDV ( :, iMomentum ( 2 ) ), &
              S_S_3_D  =>  SDV ( :, iMomentum ( 3 ) ) )

          call ComputeKernel &
                 ( C % ProperCell, Xi_J, Xi_H, Chi_J, Chi_H, J, H_1, H_2, H_3, &
                   S_1, S_2, S_3, SF_RD, M_DD_11, M_DD_22, M_DD_33, &
                   S_S_1_D, S_S_2_D, S_S_3_D, dT, S_E, S_S_1, S_S_2, S_S_3, &
                   UseDeviceOption = S % DeviceMemory )

          end associate !-- S_S_1_D
          end associate !-- SDV

        else

          call ComputeSimpleKernel &
                 ( C % ProperCell, Xi_J, Xi_H, Chi_J, Chi_H, J, H_1, H_2, H_3, &
                   M_DD_11, M_DD_22, M_DD_33, dT, S_E, S_S_1, S_S_2, S_S_3, &
                   UseDeviceOption = S % DeviceMemory )

        end if

        end associate !-- M_DD_11, etc.
        end associate !-- GSV

      class default
        call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
        call Show ( 'Slope_RM_I__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- G

      end associate !-- Xi_J, etc.
      end associate !-- IV, etc.

      class default
        call Show ( 'Chart type not recognized', CONSOLE % ERROR )
        call Show ( 'Slope_RM_I__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C
    end do !-- iC

    end associate !-- I, etc.

    if ( associated ( S % Communicator_X_1D ) ) then
      call ComputeSource_F ( S )
    end if

  end subroutine Compute


  impure elemental subroutine Finalize ( S )

    type ( Slope_RM_I_Form ), intent ( inout ) :: &
      S

    nullify ( S % Interactions )
    nullify ( S % Radiation )
    
    if ( allocated ( S % CO_SplitSource ) ) &
      deallocate ( S % CO_SplitSource )

    nullify ( S % Communicator_X_1D )

  end subroutine Finalize


  subroutine ComputeSource_F ( S )

    class ( Slope_RM_I_Form ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iC, &  !-- iChart
      iEnergy_R, iEnergy_F, &
      nSources, &
      nValues
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum_R, iMomentum_F
    real ( KDR ), dimension ( :, : ), pointer :: &
      RSB, &  !-- 2D alias for outgoing buffer
      FSB     !-- 2D alias for incoming buffer

    associate &
      ( R  =>  S % Radiation, &
        F  =>  S % Radiation % Fluid )

    call Search &
           ( R % iaBalanced, R % ENERGY_DENSITY_B, iEnergy_R )
    call Search &
           ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_1, iMomentum_R ( 1 ) )
    call Search &
           ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_2, iMomentum_R ( 2 ) )
    call Search &
           ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_3, iMomentum_R ( 3 ) )

    call Search &
           ( F % iaBalanced, F % ENERGY_DENSITY_B, iEnergy_F )
    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_1, iMomentum_F ( 1 ) )
    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_2, iMomentum_F ( 2 ) )
    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_3, iMomentum_F ( 3 ) )

    nSources  =  4

    if ( .not. allocated ( S % CO_SplitSource ) ) &
      allocate ( S % CO_SplitSource ( F % Atlas % nCharts ) )

    do iC  =  1,  F % Atlas % nCharts
      associate &
        ( CO  =>  S % CO_SplitSource ( iC ) )
      associate &
        ( RSV  =>  S % Storage ( iC ) % Value, &
          FSV  =>  F % SplitSource % Storage ( iC ) % Value )
      associate &
        ( RS_E    =>  RSV ( :, iEnergy_R ), &
          RS_S_1  =>  RSV ( :, iMomentum_R ( 1 ) ), &
          RS_S_2  =>  RSV ( :, iMomentum_R ( 2 ) ), &
          RS_S_3  =>  RSV ( :, iMomentum_R ( 3 ) ), &
          FS_G    =>  FSV ( :, iEnergy_F ), &
          FS_S_1  =>  FSV ( :, iMomentum_F ( 1 ) ), &
          FS_S_2  =>  FSV ( :, iMomentum_F ( 2 ) ), &
          FS_S_3  =>  FSV ( :, iMomentum_F ( 3 ) ) )

      nValues  =  size ( FSV, dim = 1 )

      if ( .not. allocated ( CO % Outgoing ) ) then
        call CO % Initialize &
               ( S % Communicator_X_1D, &
                 nOutgoing  =  [ nValues * nSources ], &
                 nIncoming  =  [ nValues * nSources ] )
        if ( S % DeviceMemory .and. S % DevicesCommunicate ) then
          call CO % AllocateDevice ( )
        end if
      end if
      
      if ( .not. CO % AllocatedDevice ) &
        call S % UpdateHost ( )
      
      RSB ( 1 : nValues,  1 : nSources )  =>  CO % Outgoing % Value
      FSB ( 1 : nValues,  1 : nSources )  =>  CO % Incoming % Value

      call Copy ( RS_E,   RSB ( :, 1 ), &
                  UseDeviceOption = CO % AllocatedDevice )
      call Copy ( RS_S_1, RSB ( :, 2 ), &
                  UseDeviceOption = CO % AllocatedDevice )
      call Copy ( RS_S_2, RSB ( :, 3 ), &
                  UseDeviceOption = CO % AllocatedDevice )
      call Copy ( RS_S_3, RSB ( :, 4 ), &
                  UseDeviceOption = CO % AllocatedDevice )
      call Multiply ( RSB, -1.0_KDR, &
                      UseDeviceOption = CO % AllocatedDevice )
      call CO % Reduce ( REDUCTION % SUM )
      call Copy ( FSB ( :, 1 ), FS_G, &
                  UseDeviceOption = CO % AllocatedDevice )
      call Copy ( FSB ( :, 2 ), FS_S_1, &
                  UseDeviceOption = CO % AllocatedDevice )
      call Copy ( FSB ( :, 3 ), FS_S_2, &
                  UseDeviceOption = CO % AllocatedDevice )
      call Copy ( FSB ( :, 4 ), FS_S_3, &
                  UseDeviceOption = CO % AllocatedDevice )
      
      if ( .not. CO % AllocatedDevice ) &
        call F % SplitSource % UpdateDevice ( )

      end associate !-- RS_E, etc.
      end associate !-- RSV, etc.
      end associate !-- CO
    end do !-- iC

    end associate !-- R, F

  end subroutine ComputeSource_F


end module Slope_RM_I__Form
