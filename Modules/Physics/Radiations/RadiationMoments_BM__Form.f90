module RadiationMoments_BM__Form

  !-- RadiationMoments_BaseManifold__Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluids
  use Units_R__Form

  implicit none
  private

  integer ( KDI ), private, parameter :: &
      N_FIELDS_RM    = 17, &
      N_VECTORS_RM   =  3, &
      N_PRIMITIVE_RM =  7, &
      N_BALANCED_RM  =  4

  type, public, extends ( CurrentSetForm ) :: RadiationMoments_BM_Form
    integer ( KDI ) :: &
      N_FIELDS_RM    = N_FIELDS_RM, &
      N_VECTORS_RM   = N_VECTORS_RM, &
      N_PRIMITIVE_RM = N_PRIMITIVE_RM, &
      N_BALANCED_RM  = N_BALANCED_RM
    integer ( KDI ) :: &
      ENERGY_DENSITY_C    = 0, &  !-- Comoving
      ENERGY_DENSITY_C_EQ = 0, &
      ENERGY_DENSITY_C_RD = 0, &  !-- Relative difference from equilibrium
      ENERGY_DENSITY_B    = 0     !-- Balanced
    integer ( KDI ) :: &
      MOMENTUM_DENSITY_C_U_1 = 0, &    !-- Comoving
      MOMENTUM_DENSITY_C_U_2 = 0, &
      MOMENTUM_DENSITY_C_U_3 = 0, &
      MOMENTUM_DENSITY_B_D_1 = 0, &    !-- Balanced
      MOMENTUM_DENSITY_B_D_2 = 0, &
      MOMENTUM_DENSITY_B_D_3 = 0
    integer ( KDI ), dimension ( 3 ) :: &
      MOMENTUM_DENSITY_C_U = 0, &
      MOMENTUM_DENSITY_B_D = 0
    integer ( KDI ) :: &
      FLUX_FACTOR         = 0, &
      STRESS_FACTOR       = 0, &
      STRESS_FACTOR_RD    = 0, &  !-- Relative difference from 1/3 (diffusion)
      DIFFUSION_INDICATOR = 0
    integer ( KDI ) :: &
      FLUID_VELOCITY_U_1 = 0, &
      FLUID_VELOCITY_U_2 = 0, &
      FLUID_VELOCITY_U_3 = 0
    integer ( KDI ), dimension ( 3 ) :: &
      FLUID_VELOCITY_U = 0
    integer ( KDI ) :: &
      iTimer_CFB = 0
    character ( LDL ) :: &
      RadiationType = ''!, &
  !     MomentsType = ''
    class ( Fluid_P_Form ), pointer :: &
      Fluid => null ( )
    class ( * ), pointer :: &
      Interactions => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_RM
    generic, public :: &
      Initialize => InitializeAllocate_RM
    procedure, public, pass :: &
      SetInteractions
    procedure, public, pass ( CS ) :: &
      SetStream
    procedure, public, pass ( CS ) :: &
      ComputeFromPrimitive
    procedure, public, pass :: &
      ComputeFromBalanced
    procedure, public, pass ( CS ) :: &
      ComputeEigenspeeds
    procedure, public, pass :: &
      ComputeEquilibrium
    procedure, public, pass :: &
      SetFluidVelocity
    final :: &
      Finalize
  end type RadiationMoments_BM_Form

    private :: &
      Compute_E_S_G_Kernel, &
      Compute_J_H_G_Kernel, &
      Compute_ES_G_Kernel

    interface

      module subroutine Compute_E_S_G_Kernel &
                 ( E, S_1, S_2, S_3, J, H_1, H_2, H_3, FF, SF, SF_RD, &
                   M_DD_11, M_DD_22, M_DD_33, V_1, V_2, V_3, UseDeviceOption )
        !-- Compute_BalancedEnergy_Momentum_Galileo_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          E, &
          S_1, S_2, S_3
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          J, &
          H_1, H_2, H_3, &
          FF, SF, SF_RD
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          M_DD_11, M_DD_22, M_DD_33, &
          V_1, V_2, V_3
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_E_S_G_Kernel

      module subroutine Compute_J_H_G_Kernel &
                 ( J, H_1, H_2, H_3, E, S_1, S_2, S_3, FF, SF, SF_RD, & !RM, &
                   M_DD_11, M_DD_22, M_DD_33, M_UU_11, M_UU_22, M_UU_33, &
                   V_1, V_2, V_3, UseDeviceOption )
        !-- Compute_ComovingEnergy_Momentum_Galileo_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          J, &
          H_1, H_2, H_3
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          E, &
          S_1, S_2, S_3, &
          FF, SF, SF_RD
  !      class ( RadiationMomentsForm ), intent ( in ) :: &
  !        RM
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          M_DD_11, M_DD_22, M_DD_33, &
          M_UU_11, M_UU_22, M_UU_33, &
          V_1, V_2, V_3
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_J_H_G_Kernel

      module subroutine Compute_ES_G_Kernel &
               ( H_1, H_2, H_3, H_Dim, FF, M_DD_11, M_DD_22, M_DD_33, &
                 M_UU_Dim, c, EF_P, EF_M, UseDeviceOption )
        !-- Compute_EigenspeedSet_Galileo_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          H_1, H_2, H_3, H_Dim, FF, &
          M_DD_11, M_DD_22, M_DD_33, M_UU_Dim
        real ( KDR ), intent ( in ) :: &
          c
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          EF_P, EF_M
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_ES_G_Kernel
      
    end interface


contains


  subroutine InitializeAllocate_RM &
               ( RM, F, Units_R, RadiationType, FieldOption, VectorOption, &
                 NameOption, UnitOption, VectorIndicesOption, &
                 iaPrimitiveOption, iaBalancedOption, nFieldsOption, &
                 IgnorabilityOption )

    class ( RadiationMoments_BM_Form ), intent ( inout ) :: &
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
      iV, &  !-- iVector
      iP, &  !-- iPrimitive
      iB, &  !-- iBalanced
      iC, &  !-- iChart
      oF, &  !-- oField
      oV, &  !-- oVector
      oP, &  !-- oPrimitive
      oB, &  !-- oBalanced
      nFields, &
      nVectors, &
      nPrimitive, &
      nBalanced
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaPrimitive, &
      iaBalanced
    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
    type ( QuantityForm ), dimension ( :, : ), allocatable :: &
      FieldUnit
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field, &
      Vector

    if ( RM % Type  ==  '' ) &
      RM % Type  =  'a RadiationMoments' 
    
    RM % RadiationType  =  RadiationType

    Name  =  'RadiationMoments'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    RM % Fluid  =>  F

    !-- Field indices

    oF  =  RM % N_FIELDS_CS

    RM % ENERGY_DENSITY_C        =  oF +  1
    RM % ENERGY_DENSITY_C_EQ     =  oF +  2
    RM % ENERGY_DENSITY_C_RD     =  oF +  3
    RM % ENERGY_DENSITY_B        =  oF +  4
    RM % MOMENTUM_DENSITY_C_U_1  =  oF +  5
    RM % MOMENTUM_DENSITY_C_U_2  =  oF +  6
    RM % MOMENTUM_DENSITY_C_U_3  =  oF +  7
    RM % MOMENTUM_DENSITY_B_D_1  =  oF +  8
    RM % MOMENTUM_DENSITY_B_D_2  =  oF +  9
    RM % MOMENTUM_DENSITY_B_D_3  =  oF + 10
    RM % FLUX_FACTOR             =  oF + 11
    RM % STRESS_FACTOR           =  oF + 12
    RM % STRESS_FACTOR_RD        =  oF + 13
    RM % DIFFUSION_INDICATOR     =  oF + 14
    RM % FLUID_VELOCITY_U_1      =  oF + 15
    RM % FLUID_VELOCITY_U_2      =  oF + 16
    RM % FLUID_VELOCITY_U_3      =  oF + 17

    nFields  =  oF  +  RM % N_FIELDS_RM
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    RM % MOMENTUM_DENSITY_C_U  =  [ RM % MOMENTUM_DENSITY_C_U_1, &
                                    RM % MOMENTUM_DENSITY_C_U_2, &
                                    RM % MOMENTUM_DENSITY_C_U_3 ]
    RM % MOMENTUM_DENSITY_B_D  =  [ RM % MOMENTUM_DENSITY_B_D_1, &
                                    RM % MOMENTUM_DENSITY_B_D_2, &
                                    RM % MOMENTUM_DENSITY_B_D_3 ]
    RM % FLUID_VELOCITY_U      =  [ RM % FLUID_VELOCITY_U_1, &
                                    RM % FLUID_VELOCITY_U_2, &
                                    RM % FLUID_VELOCITY_U_3 ]

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( oF + 1 : oF + RM % N_FIELDS_RM ) &
      = [ 'EnergyDensity_C      ', &
          'EnergyDensity_C_Eq   ', &
          'EnergyDensity_C_RD   ', &
          'EnergyDensity_B      ', &
          'MomentumDensity_C_U_1', &
          'MomentumDensity_C_U_2', &
          'MomentumDensity_C_U_3', &
          'MomentumDensity_B_D_1', &
          'MomentumDensity_B_D_2', &
          'MomentumDensity_B_D_3', &
          'FluxFactor           ', &
          'StressFactor         ', &
          'StressFactor_RD      ', &
          'DiffusionIndicator   ', &
          'FluidVelocity_U_1    ', &
          'FluidVelocity_U_2    ', &
          'FluidVelocity_U_3    ' ]
          
    !-- Units

    associate ( nC  =>  F % Atlas % nCharts )

    if ( present ( UnitOption ) ) then
      allocate ( FieldUnit, source = UnitOption )
    else
      allocate ( FieldUnit ( nFields, nC ) )
    end if !-- FieldOption

    do iC  =  1, nC
      FieldUnit ( RM % ENERGY_DENSITY_C, iC ) &
        =  Units_R ( iC ) % EnergyDensity
      FieldUnit ( RM % ENERGY_DENSITY_C_EQ, iC ) &
        =  Units_R ( iC ) % EnergyDensity
      FieldUnit ( RM % ENERGY_DENSITY_B, iC ) &
        =  Units_R ( iC ) % EnergyDensity
      FieldUnit ( RM % MOMENTUM_DENSITY_C_U_1, iC ) &
        =  Units_R ( iC ) % MomentumDensity_U ( 1 )
      FieldUnit ( RM % MOMENTUM_DENSITY_C_U_2, iC ) &
        =  Units_R ( iC ) % MomentumDensity_U ( 2 )
      FieldUnit ( RM % MOMENTUM_DENSITY_C_U_3, iC ) &
        =  Units_R ( iC ) % MomentumDensity_U ( 3 )
      FieldUnit ( RM % MOMENTUM_DENSITY_B_D_1, iC ) &
        =  Units_R ( iC ) % MomentumDensity_D ( 1 )
      FieldUnit ( RM % MOMENTUM_DENSITY_B_D_2, iC ) &
        =  Units_R ( iC ) % MomentumDensity_D ( 2 )
      FieldUnit ( RM % MOMENTUM_DENSITY_B_D_3, iC ) &
        =  Units_R ( iC ) % MomentumDensity_D ( 3 )
      FieldUnit ( RM % FLUID_VELOCITY_U_1, iC ) &
        =  Units_R ( iC ) % Velocity_U ( 1 )
      FieldUnit ( RM % FLUID_VELOCITY_U_2, iC ) &
        =  Units_R ( iC ) % Velocity_U ( 2 )
      FieldUnit ( RM % FLUID_VELOCITY_U_3, iC ) &
        =  Units_R ( iC ) % Velocity_U ( 3 )
    end do !-- iC

    end associate !-- nC

    !-- Vector indices

    oV  =  RM % N_VECTORS_CS

    if ( present ( VectorIndicesOption ) ) then
      nVectors  =  size ( VectorIndicesOption )
      allocate ( VectorIndices ( nVectors ) )
      do iV  =  oV  +  RM % N_VECTORS_RM  +  1,  nVectors 
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do !-- iV
    else
      nVectors  =  oV  +  RM % N_VECTORS_RM
      allocate ( VectorIndices ( nVectors ) )
    end if

    call VectorIndices ( oV + 1 ) % Initialize ( RM % MOMENTUM_DENSITY_C_U )
    call VectorIndices ( oV + 2 ) % Initialize ( RM % MOMENTUM_DENSITY_B_D )
    call VectorIndices ( oV + 3 ) % Initialize ( RM % FLUID_VELOCITY_U )

    !-- Vector names

    if ( present ( VectorOption ) ) then
      allocate ( Vector, source = VectorOption )
    else
      allocate ( Vector ( nVectors ) )
    end if !-- FieldOption

    Vector ( oV  +  1 : oV  +  RM % N_VECTORS_RM ) &
      = [ 'MomentumDensity_C_U', &
          'MomentumDensity_B_D', &
          'FluidVelocity_U    ' ]

    !-- Primitive fields

    oP  =  RM % N_PRIMITIVE_CS

    if ( present ( iaPrimitiveOption ) ) then
      nPrimitive  =  size ( iaPrimitiveOption )
      allocate ( iaPrimitive, source = iaPrimitiveOption )
    else
      nPrimitive  =  oP  +  RM % N_PRIMITIVE_RM
      allocate ( iaPrimitive ( nPrimitive ) )
    end if !-- iaPrimitiveOption

    iaPrimitive ( oP  +  1 : oP  +  RM % N_PRIMITIVE_RM )  &
      =  [ RM % ENERGY_DENSITY_C, RM % MOMENTUM_DENSITY_C_U, &
           RM % FLUID_VELOCITY_U ]

    !-- Balanced fields

    oB  =  RM % N_BALANCED_CS

    if ( present ( iaBalancedOption ) ) then
      nBalanced  =  size ( iaBalancedOption )
      allocate ( iaBalanced, source = iaBalancedOption )
    else
      nBalanced  =  oB  +  RM % N_BALANCED_RM
      allocate ( iaBalanced ( nBalanced ) )
    end if !-- iaBalancedOption

    iaBalanced ( oB  +  1 : oB  +  RM % N_BALANCED_RM )  &
      =  [ RM % ENERGY_DENSITY_B, RM % MOMENTUM_DENSITY_B_D ]

    !-- CurrentSet

    call RM % CurrentSetForm % Initialize &
           ( F % Geometry, &
             FieldOption = Field, &
             VectorOption = Vector, &
             NameOption = Name, &
             UnitOption = FieldUnit, &
             VectorIndicesOption = VectorIndices, &
             iaPrimitiveOption = iaPrimitive, &
             iaBalancedOption = iaBalanced, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

    ! !-- Parameters

    ! if ( Units_F ( 1 ) % BaryonMass  ==  UNIT % IDENTITY ) then
    !   F % BaryonMass  =  1.0_KDR
    ! else
    !   F % BaryonMass  =  CONSTANT % ATOMIC_MASS_UNIT
    ! end if

    ! F % BaryonDensityMin  =  1.0e2_KDR * sqrt ( tiny ( 0.0_KDR ) )
    ! call PROGRAM_HEADER % GetParameter &
    !        ( F % BaryonDensityMin, 'BaryonDensityMin' )

  end subroutine InitializeAllocate_RM


  subroutine SetInteractions ( RM, I )

    class ( RadiationMoments_BM_Form ), intent ( inout ) :: &
      RM
    class ( FieldSet_BM_Form ), intent ( in ), target :: &
      I

    RM % Interactions  =>  I

  end subroutine SetInteractions


  subroutine SetStream ( S, CS )

    class ( Stream_BM_Form ), intent ( inout ) :: &
      S
    class ( RadiationMoments_BM_Form ), intent ( in ) :: &
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
                    CS % DIFFUSION_INDICATOR ] )

  end subroutine SetStream


  subroutine ComputeFromPrimitive ( FS_CS, CS )

    class ( FieldSet_BM_Form ), intent ( inout ) :: &
      FS_CS
    class ( RadiationMoments_BM_Form ), intent ( in ) :: &
      CS

    integer ( KDI ) :: &
      iC

    call Show ( 'ComputeFromPrimitive', CONSOLE % INFO_6 )
    call Show ( CS % Name, 'RadiationMoments', CONSOLE % INFO_6 )

    do iC  =  1, CS % Atlas % nCharts

      associate &
        ( CSV  =>  FS_CS % Storage ( iC ) % Value )
      associate &
        ( J      =>  CSV ( :, CS % ENERGY_DENSITY_C ), &
          H_1    =>  CSV ( :, CS % MOMENTUM_DENSITY_C_U_1 ), &
          H_2    =>  CSV ( :, CS % MOMENTUM_DENSITY_C_U_2 ), &
          H_3    =>  CSV ( :, CS % MOMENTUM_DENSITY_C_U_3 ), &
          E      =>  CSV ( :, CS % ENERGY_DENSITY_B ), &
          S_1    =>  CSV ( :, CS % MOMENTUM_DENSITY_B_D_1 ), &
          S_2    =>  CSV ( :, CS % MOMENTUM_DENSITY_B_D_2 ), &
          S_3    =>  CSV ( :, CS % MOMENTUM_DENSITY_B_D_3 ), &
          FF     =>  CSV ( :, CS % FLUX_FACTOR ), &
          SF     =>  CSV ( :, CS % STRESS_FACTOR ), &
          SF_RD  =>  CSV ( :, CS % STRESS_FACTOR_RD ), &
          V_1    =>  CSV ( :, CS % FLUID_VELOCITY_U_1 ), &
          V_2    =>  CSV ( :, CS % FLUID_VELOCITY_U_2 ), &
          V_3    =>  CSV ( :, CS % FLUID_VELOCITY_U_3 ) )

      select type ( G  =>  CS % Geometry )
      class is ( Gravitation_G_Form )

        associate &
          ( GSV  =>  G % Storage ( iC ) % Value )
        associate &
          ( M_DD_11  =>  GSV ( :, G % METRIC_F_DD_11 ), &
            M_DD_22  =>  GSV ( :, G % METRIC_F_DD_22 ), &
            M_DD_33  =>  GSV ( :, G % METRIC_F_DD_33 ) )

        call Compute_E_S_G_Kernel &
               ( E, S_1, S_2, S_3, J, H_1, H_2, H_3, FF, SF, SF_RD, &
                 M_DD_11, M_DD_22, M_DD_33, V_1, V_2, V_3, &
                 UseDeviceOption = CS % DeviceMemory )

        end associate !-- M_DD_11, etc.
        end associate !-- GSV

      class default
        call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
        call Show ( 'RadiationMoments_BM__Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeFromPrimitive', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- G

      end associate !-- J, etc.
      end associate !-- CSV, etc.

    end do !-- iC

  end subroutine ComputeFromPrimitive


  subroutine ComputeFromBalanced ( CS, T_Option )

    class ( RadiationMoments_BM_Form ), intent ( inout ) :: &
      CS
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC
    type ( TimerForm ), pointer :: &
      T_K

    call Show ( 'ComputeFromBalanced', CONSOLE % INFO_6 )
    call Show ( CS % Name, 'RadiationMoments', CONSOLE % INFO_6 )

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
          E      =>  CSV ( :, CS % ENERGY_DENSITY_B ), &
          S_1    =>  CSV ( :, CS % MOMENTUM_DENSITY_B_D_1 ), &
          S_2    =>  CSV ( :, CS % MOMENTUM_DENSITY_B_D_2 ), &
          S_3    =>  CSV ( :, CS % MOMENTUM_DENSITY_B_D_3 ), &
          FF     =>  CSV ( :, CS % FLUX_FACTOR ), &
          SF     =>  CSV ( :, CS % STRESS_FACTOR ), &
          SF_RD  =>  CSV ( :, CS % STRESS_FACTOR_RD ), &
          V_1    =>  CSV ( :, CS % FLUID_VELOCITY_U_1 ), &
          V_2    =>  CSV ( :, CS % FLUID_VELOCITY_U_2 ), &
          V_3    =>  CSV ( :, CS % FLUID_VELOCITY_U_3 ) )

      select type ( G  =>  CS % Geometry )
      class is ( Gravitation_G_Form )

        associate &
          ( GSV  =>  G % Storage ( iC ) % Value )
        associate &
          ( M_DD_11  =>  GSV ( :, G % METRIC_F_DD_11 ), &
            M_DD_22  =>  GSV ( :, G % METRIC_F_DD_22 ), &
            M_DD_33  =>  GSV ( :, G % METRIC_F_DD_33 ), &
            M_UU_11  =>  GSV ( :, G % METRIC_F_UU_11 ), &
            M_UU_22  =>  GSV ( :, G % METRIC_F_UU_22 ), &
            M_UU_33  =>  GSV ( :, G % METRIC_F_UU_33 ) )

        call Compute_J_H_G_Kernel &
               ( J, H_1, H_2, H_3, E, S_1, S_2, S_3, FF, SF, SF_RD, & !RM, &
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


  subroutine ComputeEigenspeeds ( ES, CS, FS_CS, iaEigenspeeds, iC, iD )

    class ( FieldSet_BM_Form ), intent ( inout ) :: &
      ES
    class ( RadiationMoments_BM_Form ), intent ( in ) :: &
      CS
    class ( FieldSet_BM_Form ), intent ( in ) :: &
      FS_CS
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaEigenspeeds
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension

    associate ( G  =>  CS % Geometry )
    associate &
      ( ESV  =>  ES    % Storage ( iC ) % Value, &
        CSV  =>  FS_CS % Storage ( iC ) % Value, &
        GSV  =>  G     % Storage ( iC ) % Value )
    associate &
      (   EF_P    =>  ESV ( :, iaEigenspeeds ( 1 ) ), &
          EF_M    =>  ESV ( :, iaEigenspeeds ( 2 ) ), &
           H_1    =>  CSV ( :, CS % MOMENTUM_DENSITY_C_U_1 ), &
           H_2    =>  CSV ( :, CS % MOMENTUM_DENSITY_C_U_2 ), &
           H_3    =>  CSV ( :, CS % MOMENTUM_DENSITY_C_U_3 ), &
           H_Dim  =>  CSV ( :, CS % MOMENTUM_DENSITY_C_U ( iD ) ), &
          FF      =>  CSV ( :, CS % FLUX_FACTOR ), &
        M_DD_11   =>  GSV ( :, G % METRIC_F_DD ( 1 ) ), &
        M_DD_22   =>  GSV ( :, G % METRIC_F_DD ( 2 ) ), &
        M_DD_33   =>  GSV ( :, G % METRIC_F_DD ( 3 ) ), &
        M_UU_Dim  =>  GSV ( :, G % METRIC_F_UU ( iD ) ) )
 
    call Compute_ES_G_Kernel &
           ( H_1, H_2, H_3, H_Dim, FF, M_DD_11, M_DD_22, M_DD_33, M_UU_Dim, &
             CONSTANT % SPEED_OF_LIGHT, EF_P, EF_M, &
             UseDeviceOption = CS % DeviceMemory )
  
    end associate !-- EF_P, etc.
    end associate !-- ESV, etc.
    end associate !-- G


  end subroutine ComputeEigenspeeds


  subroutine ComputeEquilibrium ( RM )

    class ( RadiationMoments_BM_Form ), intent ( inout ) :: &
      RM

    call Show ( 'Should be replaced by extension', CONSOLE % ERROR )
    call Show ( 'RadiationMoments_BM__Form', 'module', CONSOLE % ERROR )
    call Show ( 'ComputeEquilibrium', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine ComputeEquilibrium


  subroutine SetFluidVelocity ( RM )

    class ( RadiationMoments_BM_Form ), intent ( inout ) :: &
      RM

    integer ( KDI ) :: &
      iC

    associate ( F  =>  RM % Fluid )
    do iC  =  1, RM % Atlas % nCharts
      associate &
        ( RSV  =>  RM % Storage ( iC ) % Value, &
          FSV  =>  F  % Storage ( iC ) % Value )

      call Copy ( FSV ( :,   F  % VELOCITY_U_1  &
                           : F  % VELOCITY_U_3 ), &
                  RSV ( :,   RM % FLUID_VELOCITY_U_1 &
                           : RM % FLUID_VELOCITY_U_3 ), &
                  UseDeviceOption = RM % DeviceMemory )

      end associate !-- RV, etc.
    end do !-- iC
    end associate !-- F

  end subroutine SetFluidVelocity


  impure elemental subroutine Finalize ( RM )

    type ( RadiationMoments_BM_Form ), intent ( inout ) :: &
      RM

    nullify ( RM % Interactions )
    nullify ( RM % Fluid )

  end subroutine Finalize


end module RadiationMoments_BM__Form
