module NeutrinoMoments_G__Form

  !-- NeutrinoMoments_Grey__Form

  use Basics
  use Mathematics
  use Fluids
  use Gravitations
  use Units_R__Form
  use Interactions_BM__Form
  use PhotonMoments_G__Form

  implicit none
  private

  integer ( KDI ), private, parameter :: &
    N_FIELDS_NM    = 7, &
    N_PRIMITIVE_NM = 1, &
    N_BALANCED_NM  = 1

  type, public, extends ( PhotonMoments_G_Form ) :: NeutrinoMoments_G_Form
    integer ( KDI ) :: &
      N_FIELDS_NM    = N_FIELDS_NM, &
      N_PRIMITIVE_NM = N_PRIMITIVE_NM, &
      N_BALANCED_NM  = N_BALANCED_NM
    integer ( KDI ) :: &
      NUMBER_DENSITY_C    = 0, &  !-- Comoving
      NUMBER_DENSITY_C_EQ = 0, &
      NUMBER_DENSITY_C_RD = 0, &  !-- Relative difference from equilibrium
      NUMBER_DENSITY_B    = 0     !-- Balanced
    integer ( KDI ) :: &
      DEGENERACY_GREY, &
      ENERGY_AVERAGE, &
      OCCUPANCY_AVERAGE
  contains
    procedure, private, pass :: &
      InitializeAllocate_RM
    final :: &
      Finalize
    procedure, public, pass ( CS ) :: &
      SetStream
    procedure, public, pass ( CS ) :: &
      ComputeFromPrimitive
    procedure, public, pass :: &
      ComputeFromBalanced
    procedure, private, pass :: &
      ComputeSpectralParametersAll
    procedure, public, pass :: &
      ComputeEquilibriumAll
  end type NeutrinoMoments_G_Form

    private :: &
      Compute_E_S_G_G_Kernel, &
      Compute_J_H_N_G_Kernel, &
      Compute_SP_Kernel, &
      Compute_Eq_Kernel
    
    interface

      module subroutine Compute_E_S_G_G_Kernel &
                 ( E, S_1, S_2, S_3, G, J, H_1, H_2, H_3, N, FF, SF, SF_RD, &
                   M_DD_11, M_DD_22, M_DD_33, V_1, V_2, V_3, UseDeviceOption )
        !-- Compute_BalancedEnergy_Momentum_Galileo_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          E, &
          S_1, S_2, S_3, &
          G
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          J, &
          H_1, H_2, H_3, &
          N, &
          FF, SF, SF_RD
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          M_DD_11, M_DD_22, M_DD_33, &
          V_1, V_2, V_3
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_E_S_G_G_Kernel

      module subroutine Compute_J_H_N_G_Kernel &
               ( J, H_1, H_2, H_3, N, E, S_1, S_2, S_3, G, FF, SF, SF_RD, &!RM,&
                 M_DD_11, M_DD_22, M_DD_33, M_UU_11, M_UU_22, M_UU_33, &
                 V_1, V_2, V_3, UseDeviceOption )
        !-- Compute_ComovingEnergy_Momentum_Galileo_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          J, &
          H_1, H_2, H_3, &
          N
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          E, &
          S_1, S_2, S_3, &
          G, &
          FF, SF, SF_RD
  !      class ( RadiationMomentsForm ), intent ( in ) :: &
  !        RM
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          M_DD_11, M_DD_22, M_DD_33, &
          M_UU_11, M_UU_22, M_UU_33, &
          V_1, V_2, V_3
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_J_H_N_G_Kernel

      module subroutine Compute_SP_Kernel &
               ( T_R, Eta_R, E_Ave, F_Ave, J, N, UseDeviceOption )
        !-- Compute_SpectralParameters_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          T_R, Eta_R, &
          E_Ave, F_Ave
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          J, N
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_SP_Kernel
 
      module subroutine Compute_Eq_Kernel &
               ( J_Eq, N_Eq, J_RD, N_RD, J, N, T, Mu_E, Mu_NP, Sign, &
                 UseDeviceOption )
        !-- Compute_Equilibrium_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          J_Eq, N_Eq, &
          J_RD, N_RD
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          J, N, &
          T, &
          Mu_E, Mu_NP
        real ( KDR ), intent ( in ) :: &
          Sign
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_Eq_Kernel

    end interface


contains


  subroutine InitializeAllocate_RM &
               ( RM, F, Units_R, RadiationType, FieldOption, VectorOption, &
                 NameOption, UnitOption, VectorIndicesOption, &
                 iaPrimitiveOption, iaBalancedOption, nFieldsOption, &
                 IgnorabilityOption )

    class ( NeutrinoMoments_G_Form ), intent ( inout ) :: &
      RM
    class ( Fluid_P_Form ), intent ( in ), target :: &
      F
    class ( Units_R_Form ), dimension ( : ), intent ( in ) :: &
      Units_R
    character ( * ), intent ( in ) :: &
      RadiationType
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( QuantityForm ), dimension ( :, : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaPrimitiveOption, &
      iaBalancedOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    integer ( KDI ) :: &
      iC, &  !-- iChart
      oF, &  !-- oField
      oP, &  !-- oPrimitive
      oB, &  !-- oBalanced
      nFields, &
      nPrimitive, &
      nBalanced
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaPrimitive, &
      iaBalanced
    type ( QuantityForm ), dimension ( :, : ), allocatable :: &
      FieldUnit
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( RM % Type  ==  '' ) &
      RM % Type  =  'a NeutrinoMoments_G' 
    
    Name  =  'NeutrinoMoments'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    !-- Field indices

    oF  =  RM % N_FIELDS_CS  +  RM % N_FIELDS_RM  +  RM % N_FIELDS_PM

    RM % NUMBER_DENSITY_C     =  oF + 1
    RM % NUMBER_DENSITY_C_EQ  =  oF + 2
    RM % NUMBER_DENSITY_C_RD  =  oF + 3
    RM % NUMBER_DENSITY_B     =  oF + 4
    RM % DEGENERACY_GREY      =  oF + 5
    RM % ENERGY_AVERAGE       =  oF + 6
    RM % OCCUPANCY_AVERAGE    =  oF + 7

    nFields  =  oF  +  RM % N_FIELDS_NM
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( oF + 1 : oF + RM % N_FIELDS_NM ) &
      = [ 'NumberDensity_C   ', &
          'NumberDensity_C_Eq', &
          'NumberDensity_C_RD', &
          'NumberDensity_B   ', &
          'DegeneracyGrey    ', &
          'EnergyAverage     ', &
          'OccupancyAverage  ' ]

    !-- Units

    associate ( nC  =>  F % Atlas % nCharts )

    if ( present ( UnitOption ) ) then
      allocate ( FieldUnit, source = UnitOption )
    else
      allocate ( FieldUnit ( nFields, nC ) )
    end if !-- FieldOption

    do iC  =  1, nC
      FieldUnit ( RM % NUMBER_DENSITY_C, iC ) &
        =  Units_R ( iC ) % NumberDensity
      FieldUnit ( RM % NUMBER_DENSITY_C_EQ, iC ) &
        =  Units_R ( iC ) % NumberDensity
      FieldUnit ( RM % NUMBER_DENSITY_B, iC ) &
        =  Units_R ( iC ) % NumberDensity
      FieldUnit ( RM % ENERGY_AVERAGE, iC ) &
        =  Units_R ( iC ) % Coordinate_MS ( 1 )
    end do !-- iC

    end associate !-- nC

    !-- Primitive fields

    oP  =  RM % N_PRIMITIVE_CS  +  RM % N_PRIMITIVE_RM

    if ( present ( iaPrimitiveOption ) ) then
      nPrimitive  =  size ( iaPrimitiveOption )
      allocate ( iaPrimitive, source = iaPrimitiveOption )
    else
      nPrimitive  =  oP  +  RM % N_PRIMITIVE_NM
      allocate ( iaPrimitive ( nPrimitive ) )
    end if !-- iaPrimitiveOption

    iaPrimitive ( oP  +  1 : oP  +  RM % N_PRIMITIVE_NM )  &
      =  [ RM % NUMBER_DENSITY_C ]

    !-- Balanced fields

    oB  =  RM % N_BALANCED_CS  +  RM % N_BALANCED_RM

    if ( present ( iaBalancedOption ) ) then
      nBalanced  =  size ( iaBalancedOption )
      allocate ( iaBalanced, source = iaBalancedOption )
    else
      nBalanced  =  oB  +  RM % N_BALANCED_NM
      allocate ( iaBalanced ( nBalanced ) )
    end if !-- iaBalancedOption

    iaBalanced ( oB  +  1 : oB  +  RM % N_BALANCED_NM )  &
      =  [ RM % NUMBER_DENSITY_B ]

    !-- PhotonMoments_G

    call RM % PhotonMoments_G_Form % Initialize &
           ( F, Units_R, RadiationType, &
             FieldOption = Field, &
             VectorOption = VectorOption, &
             NameOption = Name, &
             UnitOption = FieldUnit, &
             VectorIndicesOption = VectorIndicesOption, &
             iaPrimitiveOption = iaPrimitive, &
             iaBalancedOption = iaBalanced, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine InitializeAllocate_RM


  impure elemental subroutine Finalize ( PM )

    type ( NeutrinoMoments_G_Form ), intent ( inout ) :: &
      PM

  end subroutine Finalize


  subroutine SetStream ( S, CS )

    class ( Stream_BM_Form ), intent ( inout ) :: &
      S
    class ( NeutrinoMoments_G_Form ), intent ( in ) :: &
      CS

    call S % AddFieldSet &
           ( CS, &
             iaSelectedOption &
               =  [ CS % ENERGY_DENSITY_C, &
                    CS % ENERGY_DENSITY_C_EQ, &
                    CS % ENERGY_DENSITY_C_RD, &
                    CS % MOMENTUM_DENSITY_C_U, &
                    CS % FLUX_FACTOR, &
                    CS % STRESS_FACTOR, &
                    CS % STRESS_FACTOR_RD, &
                    CS % DIFFUSION_INDICATOR, &
                    CS % TEMPERATURE_GREY, &
                    CS % NUMBER_DENSITY_C, &
                    CS % NUMBER_DENSITY_C_EQ, &
                    CS % NUMBER_DENSITY_C_RD, &
                    CS % DEGENERACY_GREY, &
                    CS % ENERGY_AVERAGE, &
                    CS % OCCUPANCY_AVERAGE ] )

  end subroutine SetStream


  subroutine ComputeFromPrimitive ( FS_CS, CS )

    class ( FieldSet_BM_Form ), intent ( inout ) :: &
      FS_CS
    class ( NeutrinoMoments_G_Form ), intent ( in ) :: &
      CS

    integer ( KDI ) :: &
      iC

    call Show ( 'ComputeFromPrimitive', CONSOLE % INFO_6 )
    call Show ( CS % Name, 'RadiationMoments', CONSOLE % INFO_6 )

!    call CS % SetFluidVelocity ( FS_CS )

    do iC  =  1, CS % Atlas % nCharts

      associate &
        ( CSV  =>  FS_CS % Storage ( iC ) % Value )
      associate &
        ( J      =>  CSV ( :, CS % ENERGY_DENSITY_C ), &
          H_1    =>  CSV ( :, CS % MOMENTUM_DENSITY_C_U_1 ), &
          H_2    =>  CSV ( :, CS % MOMENTUM_DENSITY_C_U_2 ), &
          H_3    =>  CSV ( :, CS % MOMENTUM_DENSITY_C_U_3 ), &
          N      =>  CSV ( :, CS % NUMBER_DENSITY_C ), &
          E      =>  CSV ( :, CS % ENERGY_DENSITY_B ), &
          S_1    =>  CSV ( :, CS % MOMENTUM_DENSITY_B_D_1 ), &
          S_2    =>  CSV ( :, CS % MOMENTUM_DENSITY_B_D_2 ), &
          S_3    =>  CSV ( :, CS % MOMENTUM_DENSITY_B_D_3 ), &
          G      =>  CSV ( :, CS % NUMBER_DENSITY_B ), &
          FF     =>  CSV ( :, CS % FLUX_FACTOR ), &
          SF     =>  CSV ( :, CS % STRESS_FACTOR ), &
          SF_RD  =>  CSV ( :, CS % STRESS_FACTOR_RD ), &
          V_1    =>  CSV ( :, CS % FLUID_VELOCITY_U_1 ), &
          V_2    =>  CSV ( :, CS % FLUID_VELOCITY_U_2 ), &
          V_3    =>  CSV ( :, CS % FLUID_VELOCITY_U_3 ) )

      select type ( Grvttn  =>  CS % Geometry )
      class is ( Gravitation_G_Form )

        associate &
          ( GSV  =>  Grvttn % Storage ( iC ) % Value )
        associate &
          ( M_DD_11  =>  GSV ( :, Grvttn % METRIC_F_DD_11 ), &
            M_DD_22  =>  GSV ( :, Grvttn % METRIC_F_DD_22 ), &
            M_DD_33  =>  GSV ( :, Grvttn % METRIC_F_DD_33 ) )

        call Compute_E_S_G_G_Kernel &
               ( E, S_1, S_2, S_3, G, J, H_1, H_2, H_3, N, FF, SF, SF_RD, &
                 M_DD_11, M_DD_22, M_DD_33, V_1, V_2, V_3, &
                 UseDeviceOption = CS % DeviceMemory )

        end associate !-- M_DD_11, etc.
        end associate !-- GSV

      class default
        call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
        call Show ( 'RadiationMoments_BM__Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeFromPrimitive', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- Grvttn

      end associate !-- J, etc.
      end associate !-- CSV, etc.

    end do !-- iC

  end subroutine ComputeFromPrimitive


  subroutine ComputeFromBalanced ( CS, T_Option )

    class ( NeutrinoMoments_G_Form ), intent ( inout ) :: &
      CS
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC
    type ( TimerForm ), pointer :: &
      T_K

    call Show ( 'ComputeFromBalanced', CONSOLE % INFO_6 )
    call Show ( CS % Name, 'RadiationMoments', CONSOLE % INFO_6 )

!    call CS % SetFluidVelocity ( CS )

    if ( present ( T_Option ) ) then
      T_K  =>  PROGRAM_HEADER % Timer &
                 ( Handle = CS % iTimer_CFB, &
                   Name = trim ( T_Option % Name ) // '_Krnl', &
                   Level = T_Option % Level + 1 )
    else
      T_K  =>  null ( )
    end if

    if ( associated ( T_K ) ) call T_K % Start ( )
    do iC  =  1, CS % Atlas % nCharts

      associate &
        ( CSV  =>  CS % Storage ( iC ) % Value )
      associate &
        ( J      =>  CSV ( :, CS % ENERGY_DENSITY_C ), &
          H_1    =>  CSV ( :, CS % MOMENTUM_DENSITY_C_U_1 ), &
          H_2    =>  CSV ( :, CS % MOMENTUM_DENSITY_C_U_2 ), &
          H_3    =>  CSV ( :, CS % MOMENTUM_DENSITY_C_U_3 ), &
          N      =>  CSV ( :, CS % NUMBER_DENSITY_C ), &
          E      =>  CSV ( :, CS % ENERGY_DENSITY_B ), &
          S_1    =>  CSV ( :, CS % MOMENTUM_DENSITY_B_D_1 ), &
          S_2    =>  CSV ( :, CS % MOMENTUM_DENSITY_B_D_2 ), &
          S_3    =>  CSV ( :, CS % MOMENTUM_DENSITY_B_D_3 ), &
          G      =>  CSV ( :, CS % NUMBER_DENSITY_B ), &
          FF     =>  CSV ( :, CS % FLUX_FACTOR ), &
          SF     =>  CSV ( :, CS % STRESS_FACTOR ), &
          SF_RD  =>  CSV ( :, CS % STRESS_FACTOR_RD ), &
          V_1    =>  CSV ( :, CS % FLUID_VELOCITY_U_1 ), &
          V_2    =>  CSV ( :, CS % FLUID_VELOCITY_U_2 ), &
          V_3    =>  CSV ( :, CS % FLUID_VELOCITY_U_3 ) )

      select type ( Grvttn  =>  CS % Geometry )
      class is ( Gravitation_G_Form )

        associate &
          ( GSV  =>  Grvttn % Storage ( iC ) % Value )
        associate &
          ( M_DD_11  =>  GSV ( :, Grvttn % METRIC_F_DD_11 ), &
            M_DD_22  =>  GSV ( :, Grvttn % METRIC_F_DD_22 ), &
            M_DD_33  =>  GSV ( :, Grvttn % METRIC_F_DD_33 ), &
            M_UU_11  =>  GSV ( :, Grvttn % METRIC_F_UU_11 ), &
            M_UU_22  =>  GSV ( :, Grvttn % METRIC_F_UU_22 ), &
            M_UU_33  =>  GSV ( :, Grvttn % METRIC_F_UU_33 ) )

        call Compute_J_H_N_G_Kernel &
               ( J, H_1, H_2, H_3, N, E, S_1, S_2, S_3, G, FF, SF, SF_RD, &!RM,&
                 M_DD_11, M_DD_22, M_DD_33, M_UU_11, M_UU_22, M_UU_33, &
                 V_1, V_2, V_3, UseDeviceOption = CS % DeviceMemory )

        end associate !-- M_DD_11, etc.
        end associate !-- GSV

      class default
        call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
        call Show ( 'RadiationMoments_BM__Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeFromBalanced', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- G

      end associate !-- CSV, etc.
      end associate !-- J, etc.

    end do !-- iC
    if ( associated ( T_K ) ) call T_K % Stop ( )

  end subroutine ComputeFromBalanced


  subroutine ComputeSpectralParametersAll ( RM )

    class ( NeutrinoMoments_G_Form ), intent ( inout ) :: &
      RM

    integer ( KDI ) :: &
      iC

    call Show ( 'ComputeSpectralParametersAll', CONSOLE % INFO_6 )
    call Show ( RM % Name, 'NeutrinoMoments', CONSOLE % INFO_6 )

    do iC  =  1, RM % Atlas % nCharts
      associate &
        ( RMV  =>  RM % Storage ( iC ) % Value )
      associate &
        (   T_R    =>  RMV ( :, RM % TEMPERATURE_GREY ), &
          Eta_R    =>  RMV ( :, RM % DEGENERACY_GREY ), &
            E_Ave  =>  RMV ( :, RM % ENERGY_AVERAGE ), &
            F_Ave  =>  RMV ( :, RM % OCCUPANCY_AVERAGE ), &
            J      =>  RMV ( :, RM % ENERGY_DENSITY_C ), &
            N      =>  RMV ( :, RM % NUMBER_DENSITY_C ) )
               
      call Compute_SP_Kernel &
             ( T_R, Eta_R, E_Ave, F_Ave, J, N, &
               UseDeviceOption  =  RM % DeviceMemory )

      end associate !-- T_R, etc.
      end associate !-- RV, etc.
    end do !-- iC

  end subroutine ComputeSpectralParametersAll


  subroutine ComputeEquilibriumAll ( RM )

    class ( NeutrinoMoments_G_Form ), intent ( inout ) :: &
      RM

    integer ( KDI ) :: &
      iC

    call Show ( 'ComputeEquilibriumAll', CONSOLE % INFO_6 )
    call Show ( RM % Name, 'NeutrinoMoments', CONSOLE % INFO_6 )

    select type ( I  =>  RM % Interactions )
      class is ( Interactions_BM_Form )
    select type ( F  =>  I % Fluid )
      class is ( Fluid_P_HN_Form )

    do iC  =  1, RM % Atlas % nCharts
      associate &
        ( RMV  =>  RM % Storage ( iC ) % Value, &
           FV  =>   F % Storage ( iC ) % Value )
      associate &
        (  J_Eq  =>  RMV ( :, RM % ENERGY_DENSITY_C_EQ ), &
           N_Eq  =>  RMV ( :, RM % NUMBER_DENSITY_C_EQ ), &
           J_RD  =>  RMV ( :, RM % ENERGY_DENSITY_C_RD ), &
           N_RD  =>  RMV ( :, RM % NUMBER_DENSITY_C_RD ), &
           J     =>  RMV ( :, RM % ENERGY_DENSITY_C ), &
           N     =>  RMV ( :, RM % NUMBER_DENSITY_C ), &
           T     =>   FV ( :,  F % TEMPERATURE ), &
          Mu_E   =>   FV ( :,  F % CHEMICAL_POTENTIAL_E ), &
          Mu_NP  =>   FV ( :,  F % CHEMICAL_POTENTIAL_N_P ) )

      select case ( trim ( RM % RadiationType ) )
      case ( 'NEUTRINOS_E' )
        call Compute_Eq_Kernel &
               ( J_Eq, N_Eq, J_RD, N_RD, J, N, T, Mu_E, Mu_NP, &
                 Sign = +1.0_KDR, UseDeviceOption = RM % DeviceMemory )
      case ( 'NEUTRINOS_E_BAR' )
        call Compute_Eq_Kernel &
               ( J_Eq, N_Eq, J_RD, N_RD, J, N, T, Mu_E, Mu_NP, &
                 Sign = -1.0_KDR, UseDeviceOption = RM % DeviceMemory )
      case default
        call Show ( 'RadiationType not recognized', CONSOLE % ERROR )
        call Show ( RM % RadiationType, 'RadiationType', CONSOLE % ERROR )
        call Show ( 'NeutrinoMoments_G__Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeEquilibriumAll', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- Name

      end associate !-- T_R, etc.
      end associate !-- RV, etc.
    end do !-- iC

    end select !-- F
    end select !-- I

  end subroutine ComputeEquilibriumAll


end module NeutrinoMoments_G__Form
