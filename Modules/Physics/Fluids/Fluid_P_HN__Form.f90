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
      N_BALANCED_HN  = 1, &
      N_FIELDS_HN    = 13, &
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
      CHEMICAL_POTENTIAL_N   = 0, &
        !-- Measured relative to m_n.
      CHEMICAL_POTENTIAL_P   = 0, &
        !-- Measured relative to m_n.
      CHEMICAL_POTENTIAL_N_P = 0, &  
        !-- a.k.a. mu_hat. Includes m_n - m_p. (mu_n and mu_p both
        !   measured relative to m_n.)
      CHEMICAL_POTENTIAL_E   = 0, &
      ADIABATIC_INDEX        = 0
        !-- Includes m_e.
    real ( KDR ) :: &
      ElectronFractionMin, &
      ElectronFractionSafe
    logical ( KDL ), private :: &
      Allocated_EOS = .false.
    type ( EOS_P_HN_OConnorOtt_Form ), public, pointer :: &
      EOS => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_F
    procedure, public, pass :: &
      SetElectronFractionMin
    procedure, public, pass :: &
      SetElectronFractionSafe
    procedure, public, pass ( CS ) :: &
      SetStream
    procedure, public, pass :: &
      Show => Show_FS
    procedure, public, pass :: &
      ComputeFromTemperature
    procedure, public, pass ( CS ) :: &
      ComputeFromPrimitive
    procedure, private, pass :: &
      ComputeFromBalancedAll
    procedure, private, pass :: &
      ComputeFromBalancedSingle
    final :: &
      Finalize
!     procedure, public, pass ( C ) :: &
!       ComputeFromConservedCommon
!     procedure, public, pass ( C ) :: &
!       ComputeRawFluxes
!     procedure, public, pass ( C ) :: &
!       ComputeCenterStates
  end type Fluid_P_HN_Form
  
  public :: &
    Compute_N_V_E_YE_G_S_Kernel, &
    Apply_EOS_Prologue_S_Kernel, &
    Apply_EOS_Epilogue_S_Kernel, &
    ComputeFromBalanced_S_Kernel, &
    ComputeFromBalanced_S_V_Kernel

    private :: &
      InitializeModuleVariablesKernel, &
      Apply_EOS_Prologue_A_Kernel, &
      Compute_D_S_G_DE_G_Kernel, &
      Compute_N_V_E_YE_G_A_Kernel, &
      Apply_EOS_Epilogue_A_Kernel


    interface 
    
      module subroutine InitializeModuleVariablesKernel
        use Basics
        implicit none
      end subroutine InitializeModuleVariablesKernel
  
      module subroutine Apply_EOS_Prologue_A_Kernel &
               ( M, N, P, T, E, YE, M_Ref, N_Min, E_Min, T_Min, Y_Min, Y_Safe, &
                 UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          M, &
          N, &
          P, &
          T, &
          E, &
          YE
        real ( KDR ), intent ( in ) :: &
          M_Ref, &
          N_Min, &
          E_Min, &
          T_Min, &
          Y_Min, &
          Y_Safe
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Apply_EOS_Prologue_A_Kernel
    
      module subroutine Apply_EOS_Prologue_S_Kernel &
               ( M, N, P, T, E, YE, M_Ref, N_Min, E_Min, T_Min, Y_Min, Y_Safe, &
                 iV )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          M, &
          N, &
          P, &
          T, &
          E, &
          YE
        real ( KDR ), intent ( in ) :: &
          M_Ref, &
          N_Min, &
          E_Min, &
          T_Min, &
          Y_Min, &
          Y_Safe
        integer ( KDI ), intent ( in ) :: &
          iV
      end subroutine Apply_EOS_Prologue_S_Kernel
    
      module subroutine Compute_D_S_G_DE_G_Kernel & 	 	 
               ( N, V_1, V_2, V_3, E, YE, M, SS, M_DD_11, M_DD_22, M_DD_33, &
                 N_Min, E_Min, Y_Min, Y_Safe, D, S_1, S_2, S_3, G, DE, &
                 UseDeviceOption )
        !-- Compute_DensityB_Momentum_EnergyB_Galileo_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: & 	 	 
          N, & 	 	 
          V_1, V_2, V_3, &
          E, &
          YE
        real ( KDR ), dimension ( : ), intent ( in ) :: & 	 	 
          M,  &
          SS, &
          M_DD_11, M_DD_22, M_DD_33
        real ( KDR ), intent ( in ) :: &
          N_Min, &
          E_Min, &
          Y_Min, &
          Y_Safe
        real ( KDR ), dimension ( : ), intent ( out ) :: & 	 	 
          D, & 	 	 
          S_1, S_2, S_3, &
          G, &
          DE
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_D_S_G_DE_G_Kernel 	 	 

      module subroutine Compute_N_V_E_YE_G_A_Kernel &
               ( D, S_1, S_2, S_3, G, DE, M, M_UU_11, M_UU_22, M_UU_33, &
                 N_Min, E_Min, Y_Min, Y_Safe, N, V_1, V_2, V_3, E, YE, &
                 UseDeviceOption )
        !-- Compute_DensityC_Velocity_EnergyC_ElectronFraction_Galileo
        !    _All_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          D, &
          S_1, S_2, S_3, &
          G, &
          DE
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          M, &
          M_UU_11, M_UU_22, M_UU_33
        real ( KDR ), intent ( in ) :: &
          N_Min, &
          E_Min, &
          Y_Min, &
          Y_Safe
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          N, &
          V_1, V_2, V_3, &
          E, &
          YE
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_N_V_E_YE_G_A_Kernel

      module subroutine Compute_N_V_E_YE_G_S_Kernel &
               ( D, S_1, S_2, S_3, G, DE, M, M_UU_11, M_UU_22, M_UU_33, &
                 N_Min, E_Min, Y_Min, Y_Safe, iV, N, V_1, V_2, V_3, E, YE )
        !-- Compute_DensityC_Velocity_EnergyC_ElectronFraction_Galileo
        !    _Single_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          D, &
          S_1, S_2, S_3, &
          G, &
          DE
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          M, &
          M_UU_11, M_UU_22, M_UU_33
        real ( KDR ), intent ( in ) :: &
          N_Min, &
          E_Min, &
          Y_Min, &
          Y_Safe
        integer ( KDI ), intent ( in ) :: &
          iV
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          N, &
          V_1, V_2, V_3, &
          E, &
          YE
      end subroutine Compute_N_V_E_YE_G_S_Kernel

      module subroutine Apply_EOS_Epilogue_A_Kernel &
               ( N, P, T, SS, E, Mu_N, Mu_P, Mu_NP, Mu_E, M, Gamma, &
                 UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          N, &
          P, &
          T, &
          SS, &
          E, &
          Mu_N, &
          Mu_P, &
          Mu_NP, &
          Mu_E
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          M, &
          Gamma
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Apply_EOS_Epilogue_A_Kernel

      module subroutine Apply_EOS_Epilogue_S_Kernel &
               ( N, P, T, SS, E, Mu_N, Mu_P, Mu_NP, Mu_E, M, Gamma, iV )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          N, &
          P, &
          T, &
          SS, &
          E, &
          Mu_N, &
          Mu_P, &
          Mu_NP, &
          Mu_E
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          M, &
          Gamma
        integer ( KDI ), intent ( in ) :: &
          iV
      end subroutine Apply_EOS_Epilogue_S_Kernel
      
      module subroutine ComputeFromBalanced_S_Kernel &
               ( FV, M, N, V_1, V_2, V_3, D, G, S_1, S_2, S_3, P, T, E, YE, &
                 SS, DE, Mu_N, Mu_P, Mu_NP, Mu_E, Gamma, EOS, &
                 M_UU_11, M_UU_22, M_UU_33, T_L_N, T_L_T, T_Ye, M_Ref, N_Min, &
                 E_Min, T_Min, Y_Min, Y_Safe, E_Shift, ia_F_I, ia_F_O, ia_E, &
                 iSolve, iV )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
          FV
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          M, &
          N, &
          V_1, V_2, V_3, &
          D, &
          G, &
          S_1, S_2, S_3, &
          P, &
          T, &
          E, &
          YE, &
          SS, &
          DE, &
          Mu_N, &
          Mu_P, &
          Mu_NP, &
          Mu_E, &
          Gamma
        real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
          EOS
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          M_UU_11, M_UU_22, M_UU_33, &
          T_L_N, &      !-- TableLogDensity
          T_L_T, &      !-- TableLogTemperature
          T_Ye          !-- TableElectronFraction
        real ( KDR ), intent ( in ) :: &
          M_Ref, &
          N_Min, &
          E_Min, &
          T_Min, &
          Y_Min, &
          Y_Safe, &
          E_Shift
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          ia_F_I, &  !-- iaFluidInput
          ia_F_O, &  !-- iaFluidOutput
          ia_E       !-- iaEOS
        integer ( KDI ), intent ( in ) :: &
          iSolve, &
          iV
      end subroutine ComputeFromBalanced_S_Kernel
      
      module subroutine ComputeFromBalanced_S_V_Kernel &
               ( M, N, V_1, V_2, V_3, D, G, S_1, S_2, S_3, P, T, E, YE, &
                 SB, SS, DE, X_AA, X_A, X_N, X_P, Z, A, &
                 Mu_N, Mu_P, Mu_NP, Mu_E, Gamma, &
                 EOS, &
                 M_UU_11, M_UU_22, M_UU_33, T_L_N, T_L_T, T_Ye, M_Ref, N_Min, &
                 E_Min, T_Min, Y_Min, Y_Safe, E_Shift, ia_F_I, ia_F_O, ia_E, &
                 iSolve, iV )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          M, &
          N, &
          V_1, V_2, V_3, &
          D, &
          G, &
          S_1, S_2, S_3, &
          P, &
          T, &
          E, &
          YE, &
          SB, &
          SS, &
          DE, &
          X_AA, &
          X_A, &
          X_N, &
          X_P, &
          Z, &
          A, &
          Mu_N, &
          Mu_P, &
          Mu_NP, &
          Mu_E, &
          Gamma
        real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
          EOS
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          M_UU_11, M_UU_22, M_UU_33, &
          T_L_N, &      !-- TableLogDensity
          T_L_T, &      !-- TableLogTemperature
          T_Ye          !-- TableElectronFraction
        real ( KDR ), intent ( in ) :: &
          M_Ref, &
          N_Min, &
          E_Min, &
          T_Min, &
          Y_Min, &
          Y_Safe, &
          E_Shift
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          ia_F_I, &  !-- iaFluidInput
          ia_F_O, &  !-- iaFluidOutput
          ia_E       !-- iaEOS
        integer ( KDI ), intent ( in ) :: &
          iSolve, &
          iV
      end subroutine ComputeFromBalanced_S_V_Kernel

    end interface

    
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
    type ( QuantityForm ), dimension ( :, : ), allocatable :: &
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
    F % CHEMICAL_POTENTIAL_N    =  oF +  9
    F % CHEMICAL_POTENTIAL_P    =  oF + 10
    F % CHEMICAL_POTENTIAL_N_P  =  oF + 11
    F % CHEMICAL_POTENTIAL_E    =  oF + 12
    F % ADIABATIC_INDEX         =  of + 13

    nFields  =  oF  +  F % N_FIELDS_HN
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
          'ChemicalPotential_N  ', &
          'ChemicalPotential_P  ', &
          'ChemicalPotential_N_P', &
          'ChemicalPotential_E  ', &
          'AdiabaticIndex       ' ]

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
      FieldUnit ( F % CHEMICAL_POTENTIAL_N, iC ) &
        =  Units_F ( iC ) % Temperature
      FieldUnit ( F % CHEMICAL_POTENTIAL_P, iC ) &
        =  Units_F ( iC ) % Temperature
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

    if ( .not. EOS_Initialized ) then

      allocate ( F % EOS )
      
      call DelayFileAccess ( PROGRAM_HEADER % Communicator % Rank )
      
      EOS_Filename = '../Parameters/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5'
      call PROGRAM_HEADER % GetParameter ( EOS_Filename, 'EOS_Filename' )
      call F % EOS % Initialize ( EOS_Filename )
      
      allocate ( iaFluidOutput ( 15 ) )
      allocate ( iaSelected_EOS ( 15 ) )
      
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
            F % CHEMICAL_POTENTIAL_N_P, &
            F % CHEMICAL_POTENTIAL_N, &
            F % CHEMICAL_POTENTIAL_P, &
            F % ADIABATIC_INDEX ]

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
            F % EOS % CHEMICAL_POTENTIAL_HAT, &
            F % EOS % CHEMICAL_POTENTIAL_N, &
            F % EOS % CHEMICAL_POTENTIAL_P, &
            F % EOS % GAMMA ]
      
      call F % EOS % SelectVariables ( iaFluidOutput, iaSelected_EOS )

      F % Allocated_EOS = .true.
      EOS_Pointer => F % EOS

      call InitializeModuleVariablesKernel ( )

      EOS_RF_Accuracy     =  1.0e-9_KDR

      EOS_Initialized  =  .true.

      if ( F % DeviceMemory ) &
        call F % EOS % AllocateDevice ( )
      
    else
      F % EOS  =>  EOS_Pointer
    end if  !-- EOS_Initialized

    !-- Parameters

    associate &
      ( MassDensity_CGS => UNIT % MASS_DENSITY_CGS % Number, &
                    MeV => UNIT % MEGA_ELECTRON_VOLT % Number )
    
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
    call F % SetElectronFractionSafe  &
           ( 0.45_KDR )
    
    end associate !-- MassDensity_CGS, MeV

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


  subroutine SetElectronFractionSafe ( F, ElectronFractionSafe )

    class ( Fluid_P_HN_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      ElectronFractionSafe

    F % ElectronFractionSafe  =  ElectronFractionSafe

    call Show ( 'Setting ElectronFractionSafe of a Fluid', &
                F % IGNORABILITY + 1 )
    call Show ( F % Name, 'Name', F % IGNORABILITY + 1 )
    call Show ( F % ElectronFractionSafe, &
                F % Unit ( F % ELECTRON_FRACTION, 1 ), 'ElectronFractionSafe', &
                F % IGNORABILITY + 1 )

  end subroutine SetElectronFractionSafe


  subroutine SetStream ( S, CS )

    class ( Stream_BM_Form ), intent ( inout ) :: &
      S
    class ( Fluid_P_HN_Form ), intent ( in ) :: &
      CS

    call S % AddFieldSet &
           ( CS, &
             iaSelectedOption &
               =  [ CS % BARYON_DENSITY_C, &
                    CS % VELOCITY_U, &
                    CS % ENERGY_DENSITY_C, &
                    CS % PRESSURE, &
                    CS % TEMPERATURE, &
                    CS % ENTROPY_PER_BARYON, &
                    CS % SOUND_SPEED, &
                    CS % ELECTRON_FRACTION, &
                    CS % MASS_FRACTION_PROTON, &
                    CS % MASS_FRACTION_NEUTRON, &
                    CS % MASS_FRACTION_ALPHA, &
                    CS % MASS_FRACTION_HEAVY, &
                    CS % ATOMIC_NUMBER_HEAVY, &
                    CS % MASS_NUMBER_HEAVY, &
                    CS % CHEMICAL_POTENTIAL_N, &
                    CS % CHEMICAL_POTENTIAL_P, &
                    CS % CHEMICAL_POTENTIAL_N_P, &
                    CS % CHEMICAL_POTENTIAL_E, &
                    CS % ADIABATIC_INDEX ] )

  end subroutine SetStream


  subroutine Show_FS ( FS )

    class ( Fluid_P_HN_Form ), intent ( in ) :: &
      FS

    call FS % Fluid_P_Form % Show ( )

    call Show ( FS % ElectronFractionMin, 'ElectronFractionMin', &
                FS % IGNORABILITY )
    call Show ( FS % ElectronFractionSafe, 'ElectronFractionSafe', &
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
        (    FV   =>  F % Storage ( iC ) % Value, &
          M_Ref   =>  F % BaryonMass, &
          N_Min   =>  F % BaryonDensityMin, &
          E_Min   =>  F % EnergyDensityMin, &
          T_Min   =>  F % TemperatureMin, &
          Y_Min   =>  F % ElectronFractionMin, &
          Y_Safe  =>  F % ElectronFractionSafe )
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
          YE    =>  FV ( :, F % ELECTRON_FRACTION ), &
          DE    =>  FV ( :, F % ELECTRON_DENSITY_B ), &
          Mu_N  =>  FV ( :, F % CHEMICAL_POTENTIAL_N ), &
          Mu_P  =>  FV ( :, F % CHEMICAL_POTENTIAL_P ), &
          Mu_NP =>  FV ( :, F % CHEMICAL_POTENTIAL_N_P ), &
          Mu_E  =>  FV ( :, F % CHEMICAL_POTENTIAL_E ), &
          Gamma =>  FV ( :, F % ADIABATIC_INDEX ) )

      call Apply_EOS_Prologue_A_Kernel &
             ( M, N, P, T, E, YE, M_Ref, N_Min, E_Min, T_Min, Y_Min, Y_Safe, &
               UseDeviceOption = F % DeviceMemory )

      associate ( FS  =>  F % Storage ( iC ) )
      call FS % ReassociateHost ( AssociateVariablesOption = .false. )
      call F % EOS % ComputeFromTemperature &
             ( FS, &
               iaFluidInput = [ F % BARYON_DENSITY_C, &
                                F % TEMPERATURE, F % ELECTRON_FRACTION ] )
      call FS % ReassociateHost ( AssociateVariablesOption = .true. )
      end associate !-- FS

      call Apply_EOS_Epilogue_A_Kernel &
             ( N, P, T, SS, E, Mu_N, Mu_P, Mu_NP, Mu_E, M, Gamma, &
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
               ( N, V_1, V_2, V_3, E, YE, M, SS, M_DD_11, M_DD_22, M_DD_33, &
                 N_Min, E_Min, Y_Min, Y_Safe, D, S_1, S_2, S_3, G, DE, &
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

    class ( FieldSet_BM_Form ), intent ( inout ) :: &
      FS_CS
    class ( Fluid_P_HN_Form ), intent ( in ) :: &
      CS

    integer ( KDI ) :: &
      iC

    call Show ( 'ComputeFromPrimitive', CONSOLE % INFO_6 )
    call Show ( CS % Name, 'Fluid', CONSOLE % INFO_6 )

    do iC  =  1, CS % Atlas % nCharts

      associate &
        (    FV   =>  FS_CS % Storage ( iC ) % Value, &
          M_Ref   =>  CS % BaryonMass, &
          N_Min   =>  CS % BaryonDensityMin, &
          E_Min   =>  CS % EnergyDensityMin, &
          T_Min   =>  CS % TemperatureMin, &
          Y_Min   =>  CS % ElectronFractionMin, &
          Y_Safe  =>  CS % ElectronFractionSafe )
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
          YE    =>  FV ( :, CS % ELECTRON_FRACTION ), &
          DE    =>  FV ( :, CS % ELECTRON_DENSITY_B ), &
          Mu_N  =>  FV ( :, CS % CHEMICAL_POTENTIAL_N ), &
          Mu_P  =>  FV ( :, CS % CHEMICAL_POTENTIAL_P ), &
          Mu_NP =>  FV ( :, CS % CHEMICAL_POTENTIAL_N_P ), &
          Mu_E  =>  FV ( :, CS % CHEMICAL_POTENTIAL_E ), &
          Gamma =>  FV ( :,  CS % ADIABATIC_INDEX ) )

      associate ( CSV  =>  CS % Storage ( iC ) % Value )
      call Copy ( CSV ( :, CS % PRESSURE ), P, &
                  UseDeviceOption = CS % DeviceMemory )
      call Copy ( CSV ( :, CS % TEMPERATURE ), T, &
                  UseDeviceOption = CS % DeviceMemory )
      end associate !-- CSV

      call Apply_EOS_Prologue_A_Kernel &
             ( M, N, P, T, E, YE, M_Ref, N_Min, E_Min, T_Min, Y_Min, Y_Safe, &
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

      call Apply_EOS_Epilogue_A_Kernel &
             ( N, P, T, SS, E, Mu_N, Mu_P, Mu_NP, Mu_E, M, Gamma, &
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
               ( N, V_1, V_2, V_3, E, YE, M, SS, M_DD_11, M_DD_22, M_DD_33, &
                 N_Min, E_Min, Y_Min, Y_Safe, D, S_1, S_2, S_3, G, DE, &
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


  subroutine ComputeFromBalancedAll ( CS, T_Option )

    class ( Fluid_P_HN_Form ), intent ( inout ) :: &
      CS
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC
    type ( TimerForm ), pointer :: &
      T_G, &
      T_K

    call Show ( 'ComputeFromBalancedAll', CONSOLE % INFO_6 )
    call Show ( CS % Name, 'Fluid', CONSOLE % INFO_6 )

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
        (    FV   =>  CS % Storage ( iC ) % Value, &
          M_Ref   =>  CS % BaryonMass, &
          N_Min   =>  CS % BaryonDensityMin, &
          E_Min   =>  CS % EnergyDensityMin, &
          T_Min   =>  CS % TemperatureMin, &
          Y_Min   =>  CS % ElectronFractionMin, &
          Y_Safe  =>  CS % ElectronFractionSafe )
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
          YE    =>  FV ( :, CS % ELECTRON_FRACTION ), &
          DE    =>  FV ( :, CS % ELECTRON_DENSITY_B ), &
          Mu_N  =>  FV ( :, CS % CHEMICAL_POTENTIAL_N ), &
          Mu_P  =>  FV ( :, CS % CHEMICAL_POTENTIAL_P ), &
          Mu_NP =>  FV ( :, CS % CHEMICAL_POTENTIAL_N_P ), &
          Mu_E  =>  FV ( :, CS % CHEMICAL_POTENTIAL_E ), &
          Gamma =>  FV ( :, CS % ADIABATIC_INDEX ) )

      select type ( Gn  =>  CS % Geometry )
      class is ( Gravitation_G_Form )

        associate &
          ( GSV  =>  Gn % Storage ( iC ) % Value )
        associate &
          ( M_UU_11  =>  GSV ( :, Gn % METRIC_F_UU_11 ), &
            M_UU_22  =>  GSV ( :, Gn % METRIC_F_UU_22 ), &
            M_UU_33  =>  GSV ( :, Gn % METRIC_F_UU_33 ) )

        call Compute_N_V_E_YE_G_A_Kernel &
               ( D, S_1, S_2, S_3, G, DE, M, M_UU_11, M_UU_22, M_UU_33, &
                 N_Min, E_Min, Y_Min, Y_Safe, N, V_1, V_2, V_3, E, YE, &
                 UseDeviceOption = CS % DeviceMemory )

        end associate !-- M_UU_11, etc.
        end associate !-- GSV

      class default
        call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
        call Show ( 'Fluid_P_HN__Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeFromBalancedAll', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- Gn

      call Apply_EOS_Prologue_A_Kernel &
             ( M, N, P, T, E, YE, M_Ref, N_Min, E_Min, T_Min, Y_Min, Y_Safe, &
               UseDeviceOption = CS % DeviceMemory )

      associate ( FS  =>  CS % Storage ( iC ) )
      call FS % ReassociateHost ( AssociateVariablesOption = .false. )
      call CS % EOS % ComputeFromEnergy &
             ( FS, &
               iaFluidInput = [ CS % BARYON_DENSITY_C, &
                                CS % TEMPERATURE, CS % ELECTRON_FRACTION ], &
               iSolve = CS % ENERGY_DENSITY_C )
      call FS % ReassociateHost ( AssociateVariablesOption = .true. )
      end associate !-- FS

      call Apply_EOS_Epilogue_A_Kernel &
             ( N, P, T, SS, E, Mu_N, Mu_P, Mu_NP, Mu_E, M, Gamma, &
               UseDeviceOption = CS % DeviceMemory )

      end associate !-- M, etc.
      end associate !-- FV, etc.

    end do !-- iC
    if ( associated ( T_K ) ) call T_K % Stop ( )

  end subroutine ComputeFromBalancedAll


  subroutine ComputeFromBalancedSingle ( CS, iC, iV, T_Option )

    class ( Fluid_P_HN_Form ), intent ( inout ) :: &
      CS
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option
    integer ( KDI ), intent ( in ) :: &
      iC, &
      iV

    integer ( KDI ), dimension ( 3 ) :: &
      iaFluidInput
    type ( TimerForm ), pointer :: &
      T_G, &
      T_K

    ! call Show ( 'ComputeFromBalancedSingle', CONSOLE % INFO_6 )
    ! call Show ( CS % Name, 'Fluid', CONSOLE % INFO_6 )

    ! if ( present ( T_Option ) ) then
    !   T_K  =>  PROGRAM_HEADER % Timer &
    !              ( Handle = CS % iTimer_CFB, &
    !                Name = trim ( T_Option % Name ) // '_Krnl', &
    !                Level = T_Option % Level + 1 )
    ! else
    !   T_K  =>  null ( )
    ! end if

    ! if ( associated ( T_K ) ) call T_K % Start ( )

    associate &
      (    FV   =>  CS % Storage ( iC ) % Value, &
        M_Ref   =>  CS % BaryonMass, &
        N_Min   =>  CS % BaryonDensityMin, &
        E_Min   =>  CS % EnergyDensityMin, &
        T_Min   =>  CS % TemperatureMin, &
        Y_Min   =>  CS % ElectronFractionMin, &
        Y_Safe  =>  CS % ElectronFractionSafe )
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
        YE    =>  FV ( :, CS % ELECTRON_FRACTION ), &
        DE    =>  FV ( :, CS % ELECTRON_DENSITY_B ), &
        Mu_N  =>  FV ( :, CS % CHEMICAL_POTENTIAL_N ), &
        Mu_P  =>  FV ( :, CS % CHEMICAL_POTENTIAL_P ), &
        Mu_NP =>  FV ( :, CS % CHEMICAL_POTENTIAL_N_P ), &
        Mu_E  =>  FV ( :, CS % CHEMICAL_POTENTIAL_E ), &
        Gamma =>  FV ( :, CS % ADIABATIC_INDEX ) )

    select type ( Gn  =>  CS % Geometry )
    class is ( Gravitation_G_Form )

      associate &
        ( GSV  =>  Gn % Storage ( iC ) % Value )
      associate &
        ( M_UU_11  =>  GSV ( :, Gn % METRIC_F_UU_11 ), &
          M_UU_22  =>  GSV ( :, Gn % METRIC_F_UU_22 ), &
          M_UU_33  =>  GSV ( :, Gn % METRIC_F_UU_33 ) )
      associate &
        ( EOS     => CS % EOS % Table, &
          T_L_N   => CS % EOS % LogDensity, &
          T_L_T   => CS % EOS % LogTemperature, &
          T_Ye    => CS % EOS % ElectronFraction, &
          E_Shift => CS % EOS % EnergyShift, &
          ia_F_I  => [ CS % BARYON_DENSITY_C, &
                       CS % TEMPERATURE, CS % ELECTRON_FRACTION ], &
          ia_F_O  => CS % EOS % iaFluidOutput, &
          ia_E    => CS % EOS % iaSelected, &
          iSolve  => CS % ENERGY_DENSITY_C )
        

!--      call Compute_N_V_E_YE_G_S_Kernel &
!--             ( D, S_1, S_2, S_3, G, DE, M, M_UU_11, M_UU_22, M_UU_33, &
!--               N_Min, E_Min, Y_Min, Y_Safe, iV, N, V_1, V_2, V_3, E, YE )
      
      call ComputeFromBalanced_S_Kernel &
             ( FV, M, N, V_1, V_2, V_3, D, G, S_1, S_2, S_3, P, T, E, YE, &
               SS, DE, Mu_N, Mu_P, Mu_NP, Mu_E, Gamma, EOS, M_UU_11, M_UU_22, &
               M_UU_33, T_L_N, T_L_T, T_Ye, M_Ref, N_Min, E_Min, T_Min, &
               Y_Min, Y_Safe, E_Shift, ia_F_I, ia_F_O, ia_E, iSolve, iV )
      
      end associate !-- EOS, etc.
      end associate !-- M_UU_11, etc.
      end associate !-- GSV

    class default
      call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
      call Show ( 'Fluid_P_HN__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeFromBalancedAll', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- Gn

!--    call Apply_EOS_Prologue_S_Kernel &
!--           ( M, N, P, T, E, YE, M_Ref, N_Min, E_Min, T_Min, Y_Min, Y_Safe, iV )
!--
!--!call Show ( '>>> ComputeFromEnergy' )
!--!call Show ( iV, '>>> iV' )
!--    associate ( FS  =>  CS % Storage ( iC ) )
!--    call FS % ReassociateHost ( AssociateVariablesOption = .false. )
!--    call CS % EOS % ComputeFromEnergy &
!--           ( FS, &
!--             iaFluidInput = [ CS % BARYON_DENSITY_C, &
!--                              CS % TEMPERATURE, CS % ELECTRON_FRACTION ], &
!--             iSolve = CS % ENERGY_DENSITY_C, iV = iV )
!--    
!--    !iaFluidInput = [ CS % BARYON_DENSITY_C, &
!--    !                 CS % TEMPERATURE, CS % ELECTRON_FRACTION ]
!--    !call ComputeFromEnergy_S_Kernel &
!--    !       ( FV, CS % EOS % Table, CS % EOS % LogDensity, &
!--    !         CS % EOS % LogTemperature, CS % EOS % ElectronFraction, &
!--    !         CS % EOS % EnergyShift, iaFluidInput, CS % EOS % iaFluidOutput, &
!--    !         CS % EOS % iaSelected, iSolve = CS % ENERGY_DENSITY_C, iV = iV )
!--    call FS % ReassociateHost ( AssociateVariablesOption = .true. )
!--    end associate !-- FS
!--
!--    call Apply_EOS_Epilogue_S_Kernel &
!--           ( N, P, T, SS, E, Mu_N, Mu_P, Mu_NP, Mu_E, M, iV )

    end associate !-- M, etc.
    end associate !-- FV, etc.

!    if ( associated ( T_K ) ) call T_K % Stop ( )

  end subroutine ComputeFromBalancedSingle


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
