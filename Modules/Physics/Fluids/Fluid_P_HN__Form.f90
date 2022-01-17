module Fluid_P_HN__Form

  !-- Fluid_Perfect_HeavyNucleus__Form
  
  use Basics
  use Mathematics
  use Gravitations
  use Units_F__Form
  use Fluid_P__Form
  use EOS_P_HN_OConnorOtt__Form
  
  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_HN = 1, &
      N_BALANCED_HN = 1, &
      N_FIELDS_HN    = 11, &
      N_VECTORS_HN   = 0
    
  type, public, extends ( Fluid_P_Form ) :: Fluid_P_HN_Form
    integer ( KDI ) :: &
      N_PRIMITIVE_HN         = N_PRIMITIVE_HN, &
      N_BALANCED_HN          = N_BALANCED_HN, &
      N_FIELDS_HN            = N_FIELDS_HN, &
      N_VECTORS_HN           = N_VECTORS_HN, &
      ELECTRON_FRACTION      = 0, &
      ELECTRON_DENSITY_B     = 0, &
      MASS_FRACTION_PROTON   = 0, &
      MASS_FRACTION_NEUTRON  = 0, &
      MASS_FRACTION_ALPHA    = 0, &
      MASS_FRACTION_HEAVY    = 0, &
      ATOMIC_NUMBER_HEAVY    = 0, &
      MASS_NUMBER_HEAVY      = 0, &
      CHEMICAL_POTENTIAL_N_P = 0, &  
        !-- a.k.a. mu_hat. Includes m_n - m_p. (mu_n and mu_p both
        !   measured with respect to m_n.)
      CHEMICAL_POTENTIAL_E   = 0, &
      UNUSED_VARIABLE        = 0
        !-- Includes m_e.
    real ( KDR ) :: &
      ElectronFractionMin
    logical ( KDL ), private :: &
      Allocated_EOS = .false.
    type ( EOS_P_HN_OConnorOtt_Form ), public, pointer :: &
      EOS => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_F
    procedure, public, pass :: &
      SetElectronFractionMin
    procedure, public, pass ( CS ) :: &
      SetStream
    procedure, public, pass :: &
      Show => Show_FS
    procedure, public, pass :: &
      ComputeFromTemperature
    procedure, public, pass ( CS ) :: &
      ComputeFromPrimitive
    final :: &
      Finalize
!     procedure, public, pass ( C ) :: &
!       ComputeFromConservedCommon
!     procedure, public, pass ( C ) :: &
!       ComputeRawFluxes
!     procedure, public, pass ( C ) :: &
!       ComputeCenterStates
!     procedure, public, nopass :: &
!       Compute_DE_G_Kernel
!     procedure, public, nopass :: &
!       Compute_YE_G_Kernel
!     procedure, public, nopass :: &
!       Apply_EOS_HN_T_Kernel
!     procedure, public, nopass :: &
!       Apply_EOS_HN_SB_E_Kernel
!     procedure, public, nopass :: &
!       Apply_EOS_HN_E_Kernel
!     procedure, private, nopass :: &
!       Apply_EOS_EpilogueKernel
! !    procedure, public, nopass :: &
! !      Apply_EOS_HN_SB_Kernel
  end type Fluid_P_HN_Form

    private :: &
      Apply_EOS_PrologueKernel, &
      Compute_D_S_G_DE_G_Kernel, &
      Apply_EOS_EpilogueKernel

    interface 
  
      module subroutine Apply_EOS_PrologueKernel &
               ( M, N, P, T, E, Y, M_Ref, N_Min, E_Min, T_Min, Y_Min, &
                 UseDeviceOption )
        use Basics
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          M, &
          N, &
          P, &
          T, &
          E, &
          Y
        real ( KDR ), intent ( in ) :: &
          M_Ref, &
          N_Min, &
          E_Min, &
          T_Min, &
          Y_Min
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Apply_EOS_PrologueKernel
    
      module subroutine Compute_D_S_G_DE_G_Kernel & 	 	 
               ( N, V_1, V_2, V_3, E, M, SS, Y, M_DD_11, M_DD_22, M_DD_33, &
                 N_Min, E_Min, D, S_1, S_2, S_3, G, DE, UseDeviceOption )
        !-- Compute_DensityB_Momentum_EnergyB_Galileo_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: & 	 	 
          N, & 	 	 
          V_1, V_2, V_3, &
          E
        real ( KDR ), dimension ( : ), intent ( in ) :: & 	 	 
          M,  &
          SS, &
          Y,  &
          M_DD_11, M_DD_22, M_DD_33
        real ( KDR ), intent ( in ) :: &
          N_Min, &
          E_Min
        real ( KDR ), dimension ( : ), intent ( out ) :: & 	 	 
          D, & 	 	 
          S_1, S_2, S_3, &
          G, &
          DE
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_D_S_G_DE_G_Kernel 	 	 

      module subroutine Apply_EOS_EpilogueKernel &
               ( N, P, T, SS, E, Mu_NP, Mu_E, M, UseDeviceOption )
        use Basics
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          N, &
          P, &
          T, &
          SS, &
          E, &
          Mu_NP, &
          Mu_E
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          M
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Apply_EOS_EpilogueKernel

    end interface

    real ( KDR ), private, protected :: &
      OR_Shift, &
      MassDensity_CGS, &
      SpecificEnergy_CGS, &
      Pressure_CGS, &
      Speed_CGS, &
      MeV
    
    !-- OConnorOtt NucEOS-specific variables
    real ( KDR ), public, protected :: &
      EOS_RF_Accuracy     !-- EOS_RootFinding_Accuracy
    integer ( KDI ), public, parameter :: &
      EOS_Apply_EOS_HN_T = 1_KDI, &   !-- T input
      EOS_Apply_EOS_HN_E = 0_KDI, &   !-- E input, solve for T
      EOS_Apply_EOS_HN_S = 2_KDI      !-- S input, solve for T
    logical ( KDL ), private, protected :: &
      EOS_Initialized = .false.
    type ( EOS_P_HN_OConnorOtt_Form ), pointer, private, protected :: &
      EOS_Pointer
      

contains


  subroutine InitializeAllocate_F &
               ( F, G, Units_F, FieldOption, VectorOption, NameOption, &
                 UnitOption, VectorIndicesOption, iaPrimitiveOption, &
                 iaBalancedOption, nFieldsOption, IgnorabilityOption )

    class ( Fluid_P_HN_Form ), intent ( inout ) :: &
      F
    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    class ( Units_F_Form ), dimension ( : ), intent ( in ) :: &
      Units_F
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( MeasuredValueForm ), dimension ( :, : ), intent ( in ), optional :: &
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
      iP, &  !-- iPrimitive
      iB, &  !-- iBalanced
      iC, &  !-- iChart
      oF, &  !-- oField
      oP, &  !-- oPrimitive
      oB, &  !-- oBalanced
      nFields, &
      nPrimitive, &
      nBalanced
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaPrimitive, &
      iaBalanced, &
      iaFluidOutput, &
      iaSelected_EOS
    type ( MeasuredValueForm ), dimension ( :, : ), allocatable :: &
      FieldUnit
    character ( LDF ) :: &
      EOS_Filename
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( F % Type == '' ) &
      F % Type = 'a Fluid_P_HN'

    !-- Field indices

    oF  =  F % N_FIELDS_CS  +  F % N_FIELDS_D  +  F % N_FIELDS_P

    F % ELECTRON_FRACTION       =  oF +  1
    F % ELECTRON_DENSITY_B      =  oF +  2
    F % MASS_FRACTION_PROTON    =  oF +  3
    F % MASS_FRACTION_NEUTRON   =  oF +  4
    F % MASS_FRACTION_ALPHA     =  oF +  5
    F % MASS_FRACTION_HEAVY     =  oF +  6
    F % ATOMIC_NUMBER_HEAVY     =  oF +  7
    F % MASS_NUMBER_HEAVY       =  oF +  8
    F % CHEMICAL_POTENTIAL_N_P  =  oF +  9
    F % CHEMICAL_POTENTIAL_E    =  oF + 10
    F % UNUSED_VARIABLE         =  of + 11

    nFields  =  oF  +  F % N_FIELDS_P  +  F % N_FIELDS_HN
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( oF + 1 : oF + F % N_FIELDS_HN ) &
      = [ 'ElectronFraction     ', &
          'ElectronDensity_B    ', &
          'MassFractionProton   ', &
          'MassFractionNeutron  ', &
          'MassFractionAlpha    ', &
          'MassFractionHeavy    ', &
          'AtomicNumberHeavy    ', &
          'MassNumberHeavy      ', &
          'ChemicalPotential_N_P', &
          'ChemicalPotential_E  ', &
          'UnusedVariable       ' ]

    !-- Units

    associate ( nC  =>  G % Atlas % nCharts )

    if ( present ( UnitOption ) ) then
      allocate ( FieldUnit, source = UnitOption )
    else
      allocate ( FieldUnit ( nFields, nC ) )
    end if !-- FieldOption

    do iC  =  1, nC
      FieldUnit ( F % ELECTRON_DENSITY_B, iC ) &
        =  Units_F ( iC ) % SqrtDet_M  *  Units_F ( iC ) % NumberDensity
      FieldUnit ( F % CHEMICAL_POTENTIAL_N_P, iC ) &
        =  Units_F ( iC ) % Temperature
      FieldUnit ( F % CHEMICAL_POTENTIAL_E, iC ) &
        =  Units_F ( iC ) % Temperature
    end do !-- iC

    end associate !-- nC

    !-- Vector indices: no additional vectors

    !-- Vector names: no additional vectors

    !-- Primitive fields

    oP  =  F % N_PRIMITIVE_CS  +  F % N_PRIMITIVE_D  +  F % N_PRIMITIVE_P

    if ( present ( iaPrimitiveOption ) ) then
      nPrimitive  =  size ( iaPrimitiveOption )
      allocate ( iaPrimitive, source = iaPrimitiveOption )
    else
      nPrimitive  =  oP  +  F % N_PRIMITIVE_HN
      allocate ( iaPrimitive ( nPrimitive ) )
    end if !-- iaPrimitiveOption

    iaPrimitive ( oP  +  1 : oP  +  F % N_PRIMITIVE_HN )  &
      =  [ F % ELECTRON_FRACTION ]

    !-- Balanced fields

    oP  =  F % N_BALANCED_CS  +  F % N_BALANCED_D  +  F % N_BALANCED_P

    if ( present ( iaBalancedOption ) ) then
      nBalanced  =  size ( iaBalancedOption )
      allocate ( iaBalanced, source = iaBalancedOption )
    else
      nBalanced  =  oP  +  F % N_BALANCED_HN
      allocate ( iaBalanced ( nBalanced ) )
    end if !-- iaBalancedOption

    iaBalanced ( oP  +  1 : oP  +  F % N_BALANCED_HN )  &
      =  [ F % ELECTRON_DENSITY_B ]

    !-- Fluid_P

    call F % Fluid_P_Form % Initialize &
           ( G, Units_F, &
             FieldOption = Field, &
             VectorOption = VectorOption, &
             NameOption = NameOption, &
             UnitOption = FieldUnit, &
             VectorIndicesOption = VectorIndicesOption, &
             iaPrimitiveOption = iaPrimitive, &
             iaBalancedOption = iaBalanced, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

    call F % SetUseInitialTemperature ( .true. )

    !-- Equation of state

    if ( EOS_Initialized ) then
      F % EOS  =>  EOS_Pointer
      return
    end if

    allocate ( F % EOS )
    
    call DelayFileAccess ( PROGRAM_HEADER % Communicator % Rank )
    
    EOS_Filename = '../Parameters/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5'
    call PROGRAM_HEADER % GetParameter ( EOS_Filename, 'EOS_Filename' )
    call F % EOS % Initialize ( EOS_Filename )
    
    allocate ( iaFluidOutput ( 12 ) )
    allocate ( iaSelected_EOS ( 12 ) )
    
    iaFluidOutput &
      = [ F % ENERGY_DENSITY_C, &
          F % PRESSURE, &
          F % ENTROPY_PER_BARYON, &
          F % SOUND_SPEED, &
          F % MASS_FRACTION_ALPHA, &
          F % MASS_FRACTION_HEAVY, &
          F % MASS_FRACTION_NEUTRON, &
          F % MASS_FRACTION_PROTON, &
          F % MASS_NUMBER_HEAVY, &
          F % ATOMIC_NUMBER_HEAVY, &
          F % CHEMICAL_POTENTIAL_E, &
          F % CHEMICAL_POTENTIAL_N_P ]

    iaSelected_EOS &
      = [ F % EOS % LOG_ENERGY, &
          F % EOS % LOG_PRESSURE, &
          F % EOS % ENTROPY, &
          F % EOS % SOUND_SPEED_SQUARE, &
          F % EOS % MASS_FRACTION_A, &
          F % EOS % MASS_FRACTION_H, &
          F % EOS % MASS_FRACTION_N, &
          F % EOS % MASS_FRACTION_P, &
          F % EOS % MASS_NUMBER_BAR, &
          F % EOS % ATOMIC_NUMBER_BAR, &
          F % EOS % CHEMICAL_POTENTIAL_E, &
          F % EOS % CHEMICAL_POTENTIAL_HAT ]
    
    call F % EOS % SelectVariables ( iaFluidOutput, iaSelected_EOS )

    F % Allocated_EOS = .true.
    EOS_Pointer => F % EOS

    !-- Historical Oak Ridge Shift, accounting for nuclear binding energy
    OR_Shift = 8.9_KDR * UNIT % MEGA_ELECTRON_VOLT &
               / CONSTANT % ATOMIC_MASS_UNIT
    
    MassDensity_CGS     =  UNIT % MASS_DENSITY_CGS
    SpecificEnergy_CGS  =  UNIT % ERG  /  UNIT % GRAM
    Pressure_CGS        =  UNIT % BARYE
    Speed_CGS           =  UNIT % CENTIMETER  /  UNIT % SECOND
    MeV                 =  UNIT % MEGA_ELECTRON_VOLT
    EOS_RF_Accuracy     =  1.0e-9_KDR

    EOS_Initialized  =  .true.

    call F % SetBaryonDensityMin  &
           ( ( 10.0_KDR ** F % EOS % MinLogDensity )  *  MassDensity_CGS  &
             /  F % BaryonMass )
    call F % SetEnergyDensityMin  &
           ( ( 10.0_KDR ** F % EOS % MinLogDensity )  *  MassDensity_CGS  &
             *  1.0e-10_KDR )
    call F % SetTemperatureMin  &
           ( ( 10.0_KDR ** F % EOS % MinLogTemperature )  *  MeV )
    call F % SetElectronFractionMin  &
           ( F % EOS % MinElectronFraction )

    if ( F % DeviceMemory ) &
      call F % EOS % AllocateDevice ( )
    
  end subroutine InitializeAllocate_F


  subroutine SetElectronFractionMin ( F, ElectronFractionMin )

    class ( Fluid_P_HN_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      ElectronFractionMin

    F % ElectronFractionMin  =  ElectronFractionMin

    call Show ( 'Setting ElectronFractionMin of a Fluid', F % IGNORABILITY + 1 )
    call Show ( F % Name, 'Name', F % IGNORABILITY + 1 )
    call Show ( F % ElectronFractionMin, &
                F % Unit ( F % ELECTRON_FRACTION, 1 ), 'ElectronFractionMin', &
                F % IGNORABILITY + 1 )

  end subroutine SetElectronFractionMin


  subroutine SetStream ( S, CS )

    class ( StreamForm ), intent ( inout ) :: &
      S
    class ( Fluid_P_HN_Form ), intent ( in ) :: &
      CS

    call S % AddFieldSet &
           ( CS, &
             iaSelectedOption &
               =  [ CS % BARYON_DENSITY_C, CS % VELOCITY_U, &
                    CS % ENERGY_DENSITY_C, CS % PRESSURE, CS % TEMPERATURE, &
                    CS % ENTROPY_PER_BARYON,  CS % ELECTRON_FRACTION, &
                    CS % MASS_FRACTION_PROTON, CS % MASS_FRACTION_NEUTRON, &
                    CS % MASS_FRACTION_ALPHA, CS % MASS_FRACTION_HEAVY, &
                    CS % ATOMIC_NUMBER_HEAVY, CS % MASS_NUMBER_HEAVY, &
                    CS % CHEMICAL_POTENTIAL_N_P, &
                    CS % CHEMICAL_POTENTIAL_E ] )

  end subroutine SetStream


  subroutine Show_FS ( FS )

    class ( Fluid_P_HN_Form ), intent ( in ) :: &
      FS

    call FS % Fluid_P_Form % Show ( )

    call Show ( FS % ElectronFractionMin, 'ElectronFractionMin', &
                FS % IGNORABILITY )

  end subroutine Show_FS


  subroutine ComputeFromTemperature ( F )

    class ( Fluid_P_HN_Form ), intent ( inout ) :: &
      F

    integer ( KDI ) :: &
      iC

    call Show ( 'ComputeFromTemperature', CONSOLE % INFO_6 )
    call Show ( F % Name, 'Fluid', CONSOLE % INFO_6 )

    do iC  =  1, F % Atlas % nCharts

      associate &
        (    FV  =>  F % Storage ( iC ) % Value, &
          M_Ref  =>  F % BaryonMass, &
          N_Min  =>  F % BaryonDensityMin, &
          E_Min  =>  F % EnergyDensityMin, &
          T_Min  =>  F % TemperatureMin, &
          Y_Min  =>  F % ElectronFractionMin )
      associate &
        ( M     =>  FV ( :, F % BARYON_MASS ), &
          N     =>  FV ( :, F % BARYON_DENSITY_C ), &
          V_1   =>  FV ( :, F % VELOCITY_U_1 ), &
          V_2   =>  FV ( :, F % VELOCITY_U_2 ), &
          V_3   =>  FV ( :, F % VELOCITY_U_3 ), &
          D     =>  FV ( :, F % BARYON_DENSITY_B ), &
          S_1   =>  FV ( :, F % MOMENTUM_DENSITY_D_1 ), &
          S_2   =>  FV ( :, F % MOMENTUM_DENSITY_D_2 ), &
          S_3   =>  FV ( :, F % MOMENTUM_DENSITY_D_3 ), &
          E     =>  FV ( :, F % ENERGY_DENSITY_C ), &
          G     =>  FV ( :, F % ENERGY_DENSITY_B ), &
          P     =>  FV ( :, F % PRESSURE ), &
          T     =>  FV ( :, F % TEMPERATURE ), &
          SB    =>  FV ( :, F % ENTROPY_PER_BARYON ), &
          SS    =>  FV ( :, F % SOUND_SPEED ), &
          Y     =>  FV ( :, F % ELECTRON_FRACTION ), &
          DE    =>  FV ( :, F % ELECTRON_DENSITY_B ), &
          Mu_NP =>  FV ( :, F % CHEMICAL_POTENTIAL_N_P ), &
          Mu_E  =>  FV ( :, F % CHEMICAL_POTENTIAL_E ) )

      call Apply_EOS_PrologueKernel &
             ( M, N, P, T, E, Y, M_Ref, N_Min, E_Min, T_Min, Y_Min, &
               UseDeviceOption = F % DeviceMemory )

      associate ( FS  =>  F % Storage ( iC ) )
      call FS % ReassociateHost ( AssociateVariablesOption = .false. )
      call F % EOS % ComputeFromTemperature &
             ( FS, &
               iaFluidInput = [ F % BARYON_DENSITY_C, &
                                F % TEMPERATURE, F % ELECTRON_FRACTION ] )
      call FS % ReassociateHost ( AssociateVariablesOption = .true. )
      end associate !-- FS

      call Apply_EOS_EpilogueKernel &
             ( N, P, T, SS, E, Mu_NP, Mu_E, M, &
               UseDeviceOption = F % DeviceMemory )

      select type ( Gn  =>  F % Geometry )
      class is ( Gravitation_G_Form )

        associate &
          ( GSV  =>  Gn % Storage ( iC ) % Value )
        associate &
          ( M_DD_11  =>  GSV ( :, Gn % METRIC_F_DD_11 ), &
            M_DD_22  =>  GSV ( :, Gn % METRIC_F_DD_22 ), &
            M_DD_33  =>  GSV ( :, Gn % METRIC_F_DD_33 ) )

        call Compute_D_S_G_DE_G_Kernel & 	 	 
               ( N, V_1, V_2, V_3, E, M, SS, Y, M_DD_11, M_DD_22, M_DD_33, &
                 N_Min, E_Min, D, S_1, S_2, S_3, G, DE, &
                 UseDeviceOption = F % DeviceMemory )

        end associate !-- M_DD_11, etc.
        end associate !-- GSV

      class default
        call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
        call Show ( 'Fluid_P_HN__Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeFromTemperature', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- Gn

      end associate !-- M, etc.
      end associate !-- FV, etc.

    end do !-- iC

  end subroutine ComputeFromTemperature


  subroutine ComputeFromPrimitive ( FS_CS, CS )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS_CS
    class ( Fluid_P_HN_Form ), intent ( in ) :: &
      CS

    integer ( KDI ) :: &
      iC

    call Show ( 'ComputeFromPrimitive', CONSOLE % INFO_6 )
    call Show ( CS % Name, 'Fluid', CONSOLE % INFO_6 )

    do iC  =  1, CS % Atlas % nCharts

      associate &
        (    FV  =>  FS_CS % Storage ( iC ) % Value, &
          M_Ref  =>  CS % BaryonMass, &
          N_Min  =>  CS % BaryonDensityMin, &
          E_Min  =>  CS % EnergyDensityMin, &
          T_Min  =>  CS % TemperatureMin, &
          Y_Min  =>  CS % ElectronFractionMin )
      associate &
        ( M     =>  FV ( :, CS % BARYON_MASS ), &
          N     =>  FV ( :, CS % BARYON_DENSITY_C ), &
          V_1   =>  FV ( :, CS % VELOCITY_U_1 ), &
          V_2   =>  FV ( :, CS % VELOCITY_U_2 ), &
          V_3   =>  FV ( :, CS % VELOCITY_U_3 ), &
          D     =>  FV ( :, CS % BARYON_DENSITY_B ), &
          S_1   =>  FV ( :, CS % MOMENTUM_DENSITY_D_1 ), &
          S_2   =>  FV ( :, CS % MOMENTUM_DENSITY_D_2 ), &
          S_3   =>  FV ( :, CS % MOMENTUM_DENSITY_D_3 ), &
          E     =>  FV ( :, CS % ENERGY_DENSITY_C ), &
          G     =>  FV ( :, CS % ENERGY_DENSITY_B ), &
          P     =>  FV ( :, CS % PRESSURE ), &
          T     =>  FV ( :, CS % TEMPERATURE ), &
          SB    =>  FV ( :, CS % ENTROPY_PER_BARYON ), &
          SS    =>  FV ( :, CS % SOUND_SPEED ), &
          Y     =>  FV ( :, CS % ELECTRON_FRACTION ), &
          DE    =>  FV ( :, CS % ELECTRON_DENSITY_B ), &
          Mu_NP =>  FV ( :, CS % CHEMICAL_POTENTIAL_N_P ), &
          Mu_E  =>  FV ( :, CS % CHEMICAL_POTENTIAL_E ) )

      associate ( CSV  =>  CS % Storage ( iC ) % Value )
      call Copy ( CSV ( :, CS % PRESSURE ), P, &
                  UseDeviceOption = CS % DeviceMemory )
      call Copy ( CSV ( :, CS % TEMPERATURE ), T, &
                  UseDeviceOption = CS % DeviceMemory )
      end associate !-- CSV

      call Apply_EOS_PrologueKernel &
             ( M, N, P, T, E, Y, M_Ref, N_Min, E_Min, T_Min, Y_Min, &
               UseDeviceOption = CS % DeviceMemory )

      associate ( FS  =>  FS_CS % Storage ( iC ) )
      call FS % ReassociateHost ( AssociateVariablesOption = .false. )
      call CS % EOS % ComputeFromEnergy &
             ( FS, &
               iaFluidInput = [ CS % BARYON_DENSITY_C, &
                                CS % TEMPERATURE, CS % ELECTRON_FRACTION ], &
               iSolve = CS % ENERGY_DENSITY_C )
      call FS % ReassociateHost ( AssociateVariablesOption = .true. )
      end associate !-- FS

      call Apply_EOS_EpilogueKernel &
             ( N, P, T, SS, E, Mu_NP, Mu_E, M, &
               UseDeviceOption = CS % DeviceMemory )

      select type ( Gn  =>  CS % Geometry )
      class is ( Gravitation_G_Form )

        associate &
          ( GSV  =>  Gn % Storage ( iC ) % Value )
        associate &
          ( M_DD_11  =>  GSV ( :, Gn % METRIC_F_DD_11 ), &
            M_DD_22  =>  GSV ( :, Gn % METRIC_F_DD_22 ), &
            M_DD_33  =>  GSV ( :, Gn % METRIC_F_DD_33 ) )

        call Compute_D_S_G_DE_G_Kernel & 	 	 
               ( N, V_1, V_2, V_3, E, M, SS, Y, M_DD_11, M_DD_22, M_DD_33, &
                 N_Min, E_Min, D, S_1, S_2, S_3, G, DE, &
                 UseDeviceOption = CS % DeviceMemory )

        end associate !-- M_DD_11, etc.
        end associate !-- GSV

      class default
        call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
        call Show ( 'Fluid_P_HN__Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeFromPrimitive', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- Gn

      end associate !-- M, etc.
      end associate !-- FV, etc.

    end do !-- iC

  end subroutine ComputeFromPrimitive


  impure elemental subroutine Finalize ( F )
    
    type ( Fluid_P_HN_Form ), intent ( inout ) :: &
      F
    
    if ( F % Allocated_EOS ) then
      if ( associated ( EOS_Pointer, F % EOS ) ) &
        nullify ( EOS_Pointer )
      deallocate ( F % EOS )
      !-- FIXME: Need to deallocate device memory for EOS
    else
      nullify ( F % EOS )
    end if
    
    nullify ( F % EOS )
    
  end subroutine Finalize 
  
    
end module Fluid_P_HN__Form
