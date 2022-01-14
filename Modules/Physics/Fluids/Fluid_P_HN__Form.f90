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
    logical ( KDL ), private :: &
      Allocated_EOS = .false.
    type ( EOS_P_HN_OConnorOtt_Form ), public, pointer :: &
      EOS => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_F
    procedure, public, pass ( CS ) :: &
      SetStream
    final :: &
      Finalize
    procedure, public, pass :: &
      ComputeFromTemperature
!     procedure, public, pass ( C ) :: &
!       ComputeFromPrimitiveCommon
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
!       Apply_EOS_PrologueKernel
!     procedure, private, nopass :: &
!       Apply_EOS_EpilogueKernel
! !    procedure, public, nopass :: &
! !      Apply_EOS_HN_SB_Kernel
  end type Fluid_P_HN_Form

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

    if ( F % DeviceMemory ) &
      call F % EOS % AllocateDevice ( )
    
  end subroutine InitializeAllocate_F


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
  
    
  subroutine ComputeFromTemperature ( F )

    class ( Fluid_P_HN_Form ), intent ( inout ) :: &
      F

    integer ( KDI ) :: &
      iC

    call Show ( 'ComputeFromTemperature', CONSOLE % INFO_6 )
    call Show ( F % Name, 'Fluid', CONSOLE % INFO_6 )

    do iC  =  1, F % Atlas % nCharts

    end do !-- iC

  end subroutine ComputeFromTemperature


end module Fluid_P_HN__Form
