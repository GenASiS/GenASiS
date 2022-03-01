module Slope_DFV_N__Form

  !-- Slope_DivergenceFiniteVolume_Newton__Form

  use Basics
  use Mathematics
  use Slope_DFV_C_N__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_DFV_N_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_N_DT
    procedure, private, pass :: &
      InitializeAllocate_N_DP
    generic, public :: &
      Initialize => InitializeAllocate_N_DT, InitializeAllocate_N_DP
    final :: &
      Finalize
  end type Slope_DFV_N_Form


contains


  subroutine InitializeAllocate_N_DT &
               ( S, RS, DT, Weight_RK, iVelocity_F, iMomentum_B, &
                 iBaryonMass_F, iBaryonDensity_F, iEnergy_B, SuffixOption )

    class ( Slope_DFV_N_Form ), intent ( inout ) :: &
      S
    class ( RiemannSolver_HLL_Form ), intent ( in ), target :: &
      RS
    class ( DivergencePart_CS_Form ), intent ( in ) :: &
      DT
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Weight_RK
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iVelocity_F, &
      iMomentum_B
    integer ( KDI ), intent ( in ) :: &
      iBaryonMass_F, &
      iBaryonDensity_F, &
      iEnergy_B
    character ( * ), intent ( in ), optional :: &
      SuffixOption

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_DFV_N'

    if ( S % TimerName  ==  '' ) &
      S % TimerName  =  'S_DFV_N_' // trim ( RS % CurrentSet % Name )

    Name  =  'S_DFV_N_' // trim ( RS % CurrentSet % Name )
    if ( present ( SuffixOption ) ) &
      Name  =  trim ( Name ) // '_' // trim ( SuffixOption )

    associate ( F  =>  RS % CurrentSet )

    call S % Slope_H_Form % Initialize &
           ( F % Atlas, &
             FieldOption = F % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = F % DeviceMemory, &
             PinnedMemoryOption = F % PinnedMemory, &
             DevicesCommunicateOption = F % DevicesCommunicate, &
             nFieldsOption = F % nBalanced, &
             IgnorabilityOption = F % IGNORABILITY + 1 )

    associate ( nSC  =>  S % nComponents )

    !-- Slope component: Flat

    nSC  =  nSC + 1
    allocate ( Slope_DFV_F_DT_Form :: S % Component ( nSC ) % Element )
    select type ( SF  =>  S % Component ( nSC ) % Element )
      class is ( Slope_DFV_F_DT_Form )

    call SF % Initialize ( RS, DT, Weight_RK, SuffixOption )

    end select !-- SF

    !-- Slope component: Connection Newton

    nSC  =  nSC + 1
    allocate ( Slope_DFV_C_N_Form :: S % Component ( nSC ) % Element )
    select type ( SCN  =>  S % Component ( nSC ) % Element )
      class is ( Slope_DFV_C_N_Form )

    call SCN % Initialize &
           ( F, iVelocity_F, iMomentum_B, iBaryonMass_F, iBaryonDensity_F, &
             iEnergy_B, SuffixOption )

    end select !-- SCN

    !-- Cleanup

    end associate !-- nSC
    end associate !-- F

  end subroutine InitializeAllocate_N_DT


  subroutine InitializeAllocate_N_DP &
               ( S, RS, DP_1D, Weight_RK, iVelocity_F, iMomentum_B, &
                 iBaryonMass_F, iBaryonDensity_F, iEnergy_B, SuffixOption )

    class ( Slope_DFV_N_Form ), intent ( inout ) :: &
      S
    class ( RiemannSolver_HLL_Form ), intent ( in ), target :: &
      RS
    type ( DivergencePartElement ), dimension ( : ), intent ( in ) :: &
      DP_1D
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Weight_RK
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iVelocity_F, &
      iMomentum_B
    integer ( KDI ), intent ( in ) :: &
      iBaryonMass_F, &
      iBaryonDensity_F, &
      iEnergy_B
    character ( * ), intent ( in ), optional :: &
      SuffixOption

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_DFV_N'

    if ( S % TimerName  ==  '' ) &
      S % TimerName  =  'S_DFV_N_' // trim ( RS % CurrentSet % Name )

    Name  =  'S_DFV_N_' // trim ( RS % CurrentSet % Name )
    if ( present ( SuffixOption ) ) &
      Name  =  trim ( Name ) // '_' // trim ( SuffixOption )

    associate ( F  =>  RS % CurrentSet )

    call S % Slope_H_Form % Initialize &
           ( F % Atlas, &
             FieldOption = F % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = F % DeviceMemory, &
             PinnedMemoryOption = F % PinnedMemory, &
             DevicesCommunicateOption = F % DevicesCommunicate, &
             nFieldsOption = F % nBalanced, &
             IgnorabilityOption = F % IGNORABILITY + 1 )

    associate ( nSC  =>  S % nComponents )

    !-- Slope component: Flat

    nSC  =  nSC + 1
    allocate ( Slope_DFV_F_DP_Form :: S % Component ( nSC ) % Element )
    select type ( SF  =>  S % Component ( nSC ) % Element )
      class is ( Slope_DFV_F_DP_Form )

    call SF % Initialize ( RS, DP_1D, Weight_RK, SuffixOption )

    end select !-- SF

    !-- Slope component: Connection Newton

    nSC  =  nSC + 1
    allocate ( Slope_DFV_C_N_Form :: S % Component ( nSC ) % Element )
    select type ( SCN  =>  S % Component ( nSC ) % Element )
      class is ( Slope_DFV_C_N_Form )

    call SCN % Initialize &
           ( F, iVelocity_F, iMomentum_B, iBaryonMass_F, iBaryonDensity_F, &
             iEnergy_B, SuffixOption )

    end select !-- SCN

    !-- Cleanup

    end associate !-- nSC
    end associate !-- F

  end subroutine InitializeAllocate_N_DP


  impure elemental subroutine Finalize ( S )

    type ( Slope_DFV_N_Form ), intent ( inout ) :: &
      S

  end subroutine Finalize


end module Slope_DFV_N__Form
